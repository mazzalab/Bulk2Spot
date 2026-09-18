# Stage 5c -- Custom gene-signature scoring (GSVA)
#   one score per signature per ROI -> heatmap, optional test across segments,
#   per grouping column tests (Wilcoxon / Kruskal-Wallis) + boxplots of significant hits

snakemake@source("common.R")
start_logging(snakemake@log[[1]])

suppressPackageStartupMessages({
  library(openxlsx)
  library(GSVA)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(rstatix)
  library(patchwork)
})

inp <- snakemake@input
out <- snakemake@output
cfg <- snakemake@config
gsva_cfg <- cfg$deconvolution_gsva
stats_cfg <- cfg$deconvolution_stats
group_col <- unlist(cfg$comparisons$column)[1]
dir.create(out$gsva_stats_dir, showWarnings = FALSE, recursive = TRUE)

expr_mat <- as.matrix(readRDS(inp$normCPM_rds))
fake_obj <- readRDS(inp$fake_obj_rds)

## 1. Signatures: Excel sheet, one column per signature, gene symbols in rows ----------------

sig_table <- read.xlsx(inp$signature_file, startRow = gsva_cfg$signature_start_row)
signatures <- lapply(sig_table, function(genes) genes[!is.na(genes)])

## 2. GSVA --------------------------------------------------------------------------------------
# GSVA >= 1.50 uses parameter objects; older versions take arguments directly.

method <- gsva_cfg$method %||% "gsva"
kcdf <- gsva_cfg$kcdf %||% "Gaussian"
if ("gsvaParam" %in% getNamespaceExports("GSVA")) {
  param <- switch(method,
    gsva   = gsvaParam(expr_mat, signatures, minSize = gsva_cfg$min_sz, maxSize = gsva_cfg$max_sz, kcdf = kcdf),
    ssgsea = ssgseaParam(expr_mat, signatures, minSize = gsva_cfg$min_sz, maxSize = gsva_cfg$max_sz),
    zscore = zscoreParam(expr_mat, signatures, minSize = gsva_cfg$min_sz, maxSize = gsva_cfg$max_sz),
    plage  = plageParam(expr_mat, signatures, minSize = gsva_cfg$min_sz, maxSize = gsva_cfg$max_sz),
    stop("Unknown deconvolution_gsva.method: ", method))
  gsva_scores <- gsva(param, verbose = TRUE)
} else {
  gsva_scores <- gsva(expr_mat, signatures, method = method, kcdf = kcdf,
                      min.sz = gsva_cfg$min_sz, max.sz = gsva_cfg$max_sz, verbose = TRUE)
}
message(nrow(gsva_scores), " of ", length(signatures), " signatures scored.")

# Long table: Signature, ROI, Score + ROI annotations
meta <- flatten_df(Biobase::pData(fake_obj))
meta$ROI <- rownames(meta)
gsva_long <- as.data.frame(gsva_scores) %>%
  rownames_to_column("Signature") %>%
  pivot_longer(-Signature, names_to = "ROI", values_to = "Score") %>%
  left_join(meta, by = "ROI")

## 3. Scores across segments (pairwise Wilcoxon, BH-adjusted) ------------------------------------

if (isTRUE(gsva_cfg$segment_test_enabled)) {
  seg_stats <- gsva_long %>%
    group_by(Signature) %>%
    wilcox_test(Score ~ Segment) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance() %>%
    add_xy_position(x = "Segment", dodge = 0.8)

  p <- ggplot(gsva_long, aes(x = Segment, y = Score, fill = Segment)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(aes(color = .data[[group_col]]), width = 0.1, size = 2, alpha = 0.7) +
    facet_wrap(~Signature, scales = "free_y") +
    scale_fill_brewer(palette = "Set2") +
    geom_text(data = seg_stats, aes(x = (xmin + xmax) / 2, y = y.position, label = p.adj.signif),
              inherit.aes = FALSE, size = 4) +
    labs(x = "Segment", y = "GSVA score", title = "GSVA enrichment per signature per segment") +
    theme_classic()
  ggsave(out$segment_boxplot, p, width = 30, height = 15)
} else {
  placeholder_pdf(out$segment_boxplot, "Segment-wise GSVA test skipped:\ndeconvolution_gsva.segment_test_enabled = false")
}

## 4. Heatmap -------------------------------------------------------------------------------------

anno_cols <- intersect(unlist(gsva_cfg$heatmap_annotation_columns), colnames(meta))
heatmap_anno <- meta[, anno_cols, drop = FALSE]
pdf(out$heatmap_file, width = 15, height = 10)
pheatmap(gsva_scores, annotation_col = heatmap_anno, scale = gsva_cfg$heatmap_scale, show_colnames = FALSE)
dev.off()

## 5. Tests per grouping column (same logic as stage 5b, on GSVA scores) --------------------------

for (gcol in intersect(unlist(stats_cfg$grouping_columns), colnames(gsva_long))) {
  tests <- test_by_group(gsva_long, gcol, "Segment", "Signature", "Score")
  write.table(tests, file.path(out$gsva_stats_dir, paste0("Stats_", gcol, ".txt")), sep = "\t", row.names = FALSE, quote = FALSE)

  hits <- tests[!is.na(tests$pvalue) & tests$pvalue < stats_cfg$alpha, ]
  message(gcol, ": ", nrow(hits), " significant comparison(s) (p < ", stats_cfg$alpha, ")")
  if (nrow(hits) > 0) {
    paged_boxplots(gsva_long, hits, gcol, "Segment", "Signature", "Score", ylab = "GSVA Score", ylim = c(-1, 1),
                   out_file = file.path(out$gsva_stats_dir, paste0("Signature_significant_boxplots_", gcol, ".pdf")),
                   per_page = stats_cfg$plots_per_page)
  }
}

## 6. Save -------------------------------------------------------------------------------------------

saveRDS(gsva_scores, out$gsva_scores_rds)
saveRDS(gsva_long, out$gsva_long_rds)
write.table(gsva_long, out$gsva_long_txt, sep = "\t", row.names = FALSE, quote = FALSE)
