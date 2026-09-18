# Stage 5a -- Cell-type deconvolution (SpatialDecon)
#   Two independent variants (config `deconvolution.run`):
#     batch_corrected      -> RUV-corrected log-CPM, per-scan barplots
#     no_batch_correction  -> plain CPM
#   plus an optional CIBERSORT-ready export. Skipped variants write placeholders.

snakemake@source("common.R")
start_logging(snakemake@log[[1]])

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(SpatialExperiment)
  library(GeomxTools)
  library(SpatialDecon)
  library(standR)
})

out <- snakemake@output
cfg <- snakemake@config
dec_cfg <- cfg$deconvolution
run_cfg <- dec_cfg$run
neg_target <- cfg$spe$negative_probe_target_name
neg_source <- cfg$spe$negative_probe_source_name

if (!isTRUE(run_cfg$batch_corrected) && !isTRUE(run_cfg$no_batch_correction)) {
  stop("deconvolution.run: set at least one of batch_corrected / no_batch_correction to true.")
}

## 1. Load data --------------------------------------------------------------------------

spe_ruv <- readRDS(snakemake@input$spe_ruv_rds)
raw_counts <- assay(spe_ruv, "counts")
norm_target_data <- readRDS(snakemake@input$norm_target_rds)

# standR drops the negative-control probe, but SpatialDecon needs it for the
# background model: take it back from the stage-1 object. Its name there
# depends on the panel (e.g. WTA: "NegProbe-WTX", CTA: "Negative Probe").
exprs_df <- as.data.frame(exprs(norm_target_data))
neg_name <- intersect(c(neg_target, neg_source), rownames(exprs_df))[1]
if (is.na(neg_name)) {
  stop("Negative-control probe not found as '", neg_target, "' or '", neg_source, "'. Candidates: ",
       toString(grep("neg", rownames(exprs_df), value = TRUE, ignore.case = TRUE)))
}
neg_probe <- exprs_df[neg_name, , drop = FALSE]
rownames(neg_probe) <- neg_target

## 2. Normalized matrices ------------------------------------------------------------------

# log-CPM (TMM factors from stage 3), with RUV factors regressed out if present
dge <- spe2dge(spe_ruv)
dge$samples <- dge$samples[, !grepl("\\.1$", colnames(dge$samples))]
logcpm <- edgeR::cpm(spe_ruv, log = TRUE, prior.count = 1)
ruv_cols <- grep(paste0("^", dec_cfg$ruv_prefix), colnames(dge$samples), value = TRUE)
norm_counts <- if (length(ruv_cols) > 0) {
  limma::removeBatchEffect(logcpm, covariates = dge$samples[, ruv_cols])
} else {
  message("No RUV factor columns found: using log-CPM without further correction.")
  logcpm
}
norm_counts_matrix <- as.matrix(rbind(norm_counts, neg_probe))

# Plain CPM on raw counts (+ negative probe)
normCPM <- edgeR::cpm(as.matrix(rbind(raw_counts, neg_probe)), log = FALSE)

# GeoMxSet carrying the corrected matrix as "q_norm" (reused by stage 5c)
fake_obj <- norm_target_data
assayData(fake_obj) <- list(exprs = exprs(norm_target_data), q_norm = norm_counts_matrix)

## 3. Signature matrix (e.g. safeTME), optionally restricted to some cell types ----------------

data(list = dec_cfg$signature$data_name)
signature_full <- get(dec_cfg$signature$data_name)
cell_types <- unlist(dec_cfg$signature$cell_types_subset)
signature <- if (length(cell_types) > 0) {
  signature_full[intersect(rownames(norm_counts_matrix), rownames(signature_full)), cell_types, drop = FALSE]
} else {
  signature_full
}

celltype_colors <- function(cell_names) {
  if (isTRUE(dec_cfg$use_custom_palette)) {
    cols <- unlist(dec_cfg$custom_palette)[cell_names]
    return(unname(ifelse(is.na(cols), "grey70", cols)))
  }
  n <- length(cell_names)
  if (n <= 12) RColorBrewer::brewer.pal(max(n, 3), "Set3")[seq_len(n)] else
    colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n)
}

deconvolve <- function(norm) {
  bg <- derive_GeoMx_background(norm = norm, probepool = fData(fake_obj)$Module, negnames = neg_target)
  spatialdecon(norm = norm, bg = bg, X = signature)
}

# Stacked barplot of proportions, ROIs ordered by hierarchical clustering
ordered_barplot <- function(prop, path) {
  prop[is.na(prop)] <- 0
  pdf(path, paper = "a4r")
  barplot(prop[, hclust(dist(t(prop)))$order], cex.names = 0.35, col = celltype_colors(rownames(prop)), las = 2)
  dev.off()
}

## 4. Deconvolution -----------------------------------------------------------------------------

res_bc <- NULL
if (isTRUE(run_cfg$batch_corrected)) {
  res_bc <- deconvolve(norm_counts_matrix)
  ordered_barplot(res_bc$prop_of_nontumor, out$barplot_ordered)

  # One panel per scan (column names compared after make.names(): "slide name" == "slide.name")
  pheno <- pData(fake_obj)
  names(pheno) <- make.names(names(pheno))
  scan_ids <- pheno[[make.names(cfg$spe$scan_id_column)]]
  if (is.null(scan_ids)) stop("spe.scan_id_column ('", cfg$spe$scan_id_column, "') not found in the annotation.")
  prop_all <- res_bc$prop_of_all
  prop_all[is.na(prop_all)] <- 0
  colors <- celltype_colors(rownames(prop_all))
  pdf(out$barplot_by_scan, paper = "a4r")
  par(mfrow = c(2, 4), mar = c(5, 4, 3, 1), oma = c(0, 0, 2, 0))
  for (scan in unique(scan_ids)) {
    barplot(prop_all[, scan_ids == scan, drop = FALSE], col = colors, las = 2, cex.names = 0.6,
            main = paste("Scan:", scan), border = NA)
  }
  par(fig = c(0, 1, 0, 1), new = TRUE)
  plot.new()
  legend("bottom", legend = rownames(prop_all), fill = colors, horiz = TRUE, cex = 0.7, bty = "n")
  dev.off()
  saveRDS(res_bc, out$res_dec_batch_corrected_rds)
} else {
  placeholder_pdf(out$barplot_ordered, "Batch-corrected deconvolution skipped:\ndeconvolution.run.batch_corrected = false")
  placeholder_pdf(out$barplot_by_scan, "Batch-corrected deconvolution skipped:\ndeconvolution.run.batch_corrected = false")
  saveRDS(NULL, out$res_dec_batch_corrected_rds)
}

res_nb <- NULL
if (isTRUE(run_cfg$no_batch_correction)) {
  res_nb <- deconvolve(normCPM)
  ordered_barplot(res_nb$prop_of_nontumor, out$barplot_ordered_nb)
  saveRDS(res_nb, out$res_dec_no_batch_correction_rds)
} else {
  placeholder_pdf(out$barplot_ordered_nb, "Deconvolution without batch correction skipped:\ndeconvolution.run.no_batch_correction = false")
  saveRDS(NULL, out$res_dec_no_batch_correction_rds)
}

## 5. Proportions table (batch-corrected variant if it ran) + annotation columns ---------------

res <- if (!is.null(res_bc)) res_bc else res_nb
prop <- as.data.frame(t(res$prop_of_nontumor))
prop$ROI <- rownames(prop)
for (name in names(dec_cfg$export_columns)) {
  prop[[name]] <- pData(fake_obj)[[dec_cfg$export_columns[[name]]]]   # same ROI order
}
write.table(prop, out$proportions_file, sep = "\t", row.names = FALSE, quote = FALSE)

## 6. CIBERSORT export (signature + mixture on shared genes) ------------------------------------

if (isTRUE(run_cfg$export_cibersort)) {
  genes <- intersect(rownames(norm_counts), rownames(signature_full))
  write.table(signature_full[genes, ], out$cibersort_signature_file, sep = "\t", quote = FALSE, col.names = NA)
  write.table(as.matrix(norm_counts)[genes, ], out$cibersort_mixture_file, sep = "\t", quote = FALSE, col.names = NA)
} else {
  writeLines("# deconvolution.run.export_cibersort = false", out$cibersort_signature_file)
  writeLines("# deconvolution.run.export_cibersort = false", out$cibersort_mixture_file)
}

saveRDS(fake_obj, out$fake_obj_rds)
saveRDS(norm_counts_matrix, out$norm_counts_matrix_rds)
saveRDS(normCPM, out$normCPM_rds)
