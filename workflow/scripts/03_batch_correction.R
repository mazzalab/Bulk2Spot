# Stage 3 -- QC diagnostics and batch correction (standR)
#   sample info + RLE + PCA (before) -> TMM normalization -> optional RUV4 (+ limma) ->
#   final PCA/UMAP (exported for the report) -> diagnostics after correction
#
# Batch correction runs only if `batch_correction.perform` is true AND the batch
# column exists in the data. Otherwise the post-correction PDFs are placeholders.

snakemake@source("common.R")
start_logging(snakemake@log[[1]])

suppressPackageStartupMessages({
  library(standR)
  library(scater)
  library(SpatialExperiment)
  library(ggplot2)
  library(ggsci)
  library(cowplot)
  library(limma)
  library(ggalluvial)
})

out <- snakemake@output
cfg <- snakemake@config
batch_cfg <- cfg$batch_correction
seed <- cfg$dimensionality_reduction$seed
group_col <- unlist(cfg$comparisons$column)[1]
batch_col <- batch_cfg$batch_column

spe <- readRDS(snakemake@input$spe_rds)

has_batch <- batch_col %in% colnames(colData(spe))
if (has_batch) {
  spe[[batch_col]] <- as.factor(spe[[batch_col]])
} else {
  message("Batch column '", batch_col, "' not found: batch-related plots and correction are skipped.")
}
run_correction <- isTRUE(batch_cfg$perform) && has_batch
skip_reason <- if (!has_batch) {
  sprintf("no batch column ('%s') in this dataset", batch_col)
} else {
  "batch_correction.perform = false in config"
}

# Grid of PCA plots coloured by group / Scan / Segment (+ batch, if present)
pca_grid <- function(obj) {
  pca <- reducedDim(obj, "PCA")
  colour_by <- function(col) {
    drawPCA(obj, precomputed = pca, col = "black") +
      geom_point(aes(colour = .data[[col]]), size = 3, alpha = 0.9) + theme_classic()
  }
  plots <- list(colour_by(group_col), colour_by("Scan"), colour_by("Segment"))
  if (has_batch) {
    plots[[4]] <- drawPCA(obj, precomputed = pca, col = "black") +
      geom_point(aes(fill = .data[[batch_col]]), shape = 21, size = 3, alpha = 0.9) + theme_classic()
  }
  plot_grid(plotlist = plots, labels = letters[2 + seq_along(plots)], ncol = 2, align = "hv", axis = "r")
}

scree_plot <- function(obj) {
  p <- plotScreePCA(obj, precomputed = reducedDim(obj, "PCA")) +
    theme_classic() + theme(axis.text.x = element_text(angle = 90))
  plot_grid(p, NULL, rel_widths = c(1, 0.1), ncol = 2)
}

## 1. Sample overview (alluvial plot: first column -> second column -> group) ---------

first_col <- batch_cfg$sample_info_first_column %||% "Patient"
second_col <- batch_cfg$sample_info_column %||% "Segment"
# Show the first column as anonymized labels (<prefix>_01, <prefix>_02 ...)
spe_plot <- spe
first_values <- as.character(spe_plot[[first_col]])
spe_plot[[first_col]] <- paste0(cfg$spe$patient_label_prefix %||% "Pt", "_",
                                sprintf("%02d", match(first_values, unique(first_values))))

pdf(out$sample_info_plot, width = 10, height = 11)
print(plotSampleInfo(spe_plot, column2plot = c(first_col, second_col, group_col)) +
        theme_void() + theme(legend.position = "none"))
dev.off()

## 2. RLE + PCA before correction ------------------------------------------------------

rle_raw <- plotRLExpr(spe, assay = 1, ordannots = "Patient") + xlab("") + ylab("RLE") + theme_classic()
rle_panel <- if (has_batch) {
  rle_log <- plotRLExpr(spe, assay = 2, ordannots = "Patient", fill = .data[[batch_col]]) +
    xlab("") + ylab("RLE") + theme_classic()
  plot_grid(rle_raw, rle_log, labels = "auto", ncol = 2, rel_widths = c(0.4, 0.6))
} else {
  rle_raw
}

set.seed(seed)
spe <- runPCA(spe)

pdf(out$pca_scree_prebatch, height = 13, width = 12)
print(plot_grid(rle_panel, pca_grid(spe), scree_plot(spe), labels = c("", "", "g"),
                ncol = 1, rel_heights = c(1.2, 3, 1.5)))
dev.off()

## 3. TMM normalization + optional RUV4 / limma correction --------------------------------

spe <- geomxNorm(spe, method = cfg$normalization_standR$method)

if (run_correction) {
  factor_of_interest <- batch_cfg$factor_of_interest

  # Negative control genes: least variable across batches
  spe <- findNCGs(spe, batch_name = batch_col, top_n = batch_cfg$ncg_top_n)

  # Choose k (number of unwanted factors) by silhouette score, unless fixed in config
  best_k_plot <- findBestK(spe, factor_of_int = factor_of_interest, factor_batch = batch_col,
                           NCGs = S4Vectors::metadata(spe)$NCGs) +
    ylab("Silhouette Score") + theme_classic()
  best_k <- best_k_plot$data$k[which.max(best_k_plot$data$silhouette)]
  k <- if (identical(batch_cfg$k, "auto")) best_k else as.numeric(batch_cfg$k)
  message("RUV4 with k = ", k, " (best k by silhouette: ", best_k, ")")

  spe_ruv <- geomxBatchCorrection(spe, factors = factor_of_interest,
                                  NCGs = S4Vectors::metadata(spe)$NCGs, k = k)

  # Optional second pass on the log-counts; done before any plotting so all
  # plots and exports show the final corrected data
  if (isTRUE(batch_cfg$apply_limma_residual_correction)) {
    assay(spe_ruv, "logcounts") <- removeBatchEffect(assay(spe_ruv, "logcounts"), batch = spe_ruv[[batch_col]])
    message("Applied limma::removeBatchEffect() on top of RUV4.")
  }
} else {
  message("Batch correction skipped (", skip_reason, "): using the TMM-normalized data.")
  spe_ruv <- spe
}

## 4. Final PCA + UMAP (exported for the HTML report) ----------------------------------------

set.seed(seed)
spe_ruv <- runPCA(spe_ruv, exprs_values = "logcounts")
spe_ruv <- runUMAP(spe_ruv, dimred = "PCA")

meta_cols <- intersect(unique(c(group_col, "Segment", "Scan", "Patient", if (has_batch) batch_col)),
                       colnames(colData(spe_ruv)))

export_coordinates <- function(mat, prefix, path) {
  colnames(mat) <- paste0(prefix, seq_len(ncol(mat)))
  df <- data.frame(SampleID = rownames(mat), mat, check.names = FALSE)
  for (col in meta_cols) df[[col]] <- spe_ruv[[col]]
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
}
pca_mat <- reducedDim(spe_ruv, "PCA")
percent_var <- attr(pca_mat, "percentVar")
export_coordinates(pca_mat, "PC", out$pca_final_txt)
export_coordinates(reducedDim(spe_ruv, "UMAP"), "UMAP", out$umap_final_txt)
write.table(data.frame(PC = paste0("PC", seq_along(percent_var)), percentVar = percent_var),
            out$pca_variance_txt, sep = "\t", row.names = FALSE, quote = FALSE)

## 5. Diagnostics after correction ---------------------------------------------------------

if (run_correction) {
  rle_group <- plotRLExpr(spe_ruv, ordannots = group_col, assay = 2, color = "black") +
    ggtitle("") + ylab("RLE") + xlab("") + theme_classic() + theme(legend.position = "none")
  rle_batch <- plotRLExpr(spe_ruv, ordannots = batch_col, assay = 2, color = "black", fill = .data[[batch_col]]) +
    ggtitle("") + scale_fill_nejm() + ylab("RLE") + xlab("") + theme_classic()

  pdf(out$pca_scree_postbatch, height = 18, width = 20)
  print(plot_grid(plot_grid(rle_group, rle_batch, labels = "auto", ncol = 2, rel_widths = c(0.4, 0.6)),
                  pca_grid(spe_ruv), scree_plot(spe_ruv), labels = c("", "", "g"),
                  ncol = 1, rel_heights = c(1.2, 3, 1.5)))
  dev.off()

  # Pairwise PC1-PC4 scatter plots, filled by one column
  pc_labels <- paste0("PC", 1:4, " (", round(percent_var[1:4], 1), "%)")
  pca_df <- as.data.frame(pca_mat[, 1:4])
  colnames(pca_df) <- paste0("PC", 1:4)
  for (col in unique(c(group_col, "Segment", batch_col))) pca_df[[col]] <- spe_ruv[[col]]

  pairwise_grid <- function(fill_col) {
    scatter <- function(i, j, legend = "none") {
      ggplot(pca_df, aes(x = .data[[paste0("PC", i)]], y = .data[[paste0("PC", j)]], fill = .data[[fill_col]])) +
        geom_point(shape = 21, size = 3, alpha = 0.9, color = "black") +
        labs(x = pc_labels[i], y = pc_labels[j]) +
        theme_classic() + theme(legend.position = legend)
    }
    legend <- get_legend(scatter(1, 2, legend = "right"))
    plot_grid(plot_grid(scatter(1, 2), scatter(1, 3), scatter(1, 4), ncol = 3),
              plot_grid(legend, scatter(2, 3), scatter(2, 4), ncol = 3),
              plot_grid(NULL, NULL, scatter(3, 4), ncol = 3),
              ncol = 1, align = "hv")
  }

  pdf(out$ruv_plot, height = 18, width = 20)
  print(plot_grid(plot_grid(best_k_plot, rle_group, ncol = 1, labels = c("a", "b")),
                  pairwise_grid(group_col), pairwise_grid("Segment"), pairwise_grid(batch_col),
                  ncol = 2, labels = c("", "c", "d", "e")))
  dev.off()
} else {
  placeholder_pdf(out$pca_scree_postbatch, paste0("Post-batch PCA skipped:\n", skip_reason))
  placeholder_pdf(out$ruv_plot, paste0("Batch correction skipped:\n", skip_reason))
}

saveRDS(spe_ruv, out$spe_ruv_rds)
message("Saved SpatialExperiment to ", out$spe_ruv_rds)
