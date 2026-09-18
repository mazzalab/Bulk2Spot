# Stage 1 -- Preprocessing (GeomxTools)
#   load DCC/PKC/annotation -> segment QC -> probe QC -> gene-level counts ->
#   LOQ-based segment/gene filtering -> Q3 normalization

snakemake@source("common.R")
start_logging(snakemake@log[[1]])

suppressPackageStartupMessages({
  library(NanoStringNCTools)
  library(GeomxTools)
  library(ggplot2)
  library(gridExtra)
  library(cowplot)
  library(openxlsx)
  library(readxl)
})

## 1. Inputs, outputs, config --------------------------------------------------

inp <- snakemake@input
out <- snakemake@output
cfg <- snakemake@config

ann_cfg   <- cfg$annotation
seg_qc    <- cfg$segment_qc
probe_qc  <- cfg$probe_qc
loq_cfg   <- cfg$loq
filt_cfg  <- cfg$gene_filtering
norm_cfg  <- cfg$normalization
plot_cfg  <- cfg$plotting
sheet     <- cfg$annotation_sheet
dcc_col   <- ann_cfg$pheno_data_dcc_colname
group_col <- unlist(cfg$comparisons$column)[1]   # primary biological grouping, used for plots

# Keeps track of features x segments after every filtering step (-> QC_log.txt)
qc_log <- data.frame()
log_step <- function(step, obj) {
  qc_log <<- rbind(qc_log, data.frame(Step = step, Features = nrow(obj), Samples = ncol(obj)))
}

## 2. Load data ----------------------------------------------------------------

dcc_files <- dir(inp$dcc_dir, pattern = "\\.dcc$", full.names = TRUE, recursive = TRUE)
for (pattern in unlist(cfg$exclude_dcc_files)) {
  dcc_files <- dcc_files[!grepl(pattern, dcc_files, fixed = TRUE)]
}

# readNanoStringGeoMxSet() expects the annotation sorted by DCC name; write a
# sorted copy into the output directory (the raw annotation is never modified).
anno <- as.data.frame(read_xlsx(inp$annotation, sheet = sheet))
if (!dcc_col %in% colnames(anno)) {
  stop("annotation.pheno_data_dcc_colname ('", dcc_col, "') is not a column of ", inp$annotation)
}
wb <- createWorkbook()
addWorksheet(wb, sheet)
writeData(wb, sheet, anno[order(anno[[dcc_col]]), ], keepNA = FALSE)
saveWorkbook(wb, out$annotation_reordered, overwrite = TRUE)

raw_data <- readNanoStringGeoMxSet(
  dccFiles = dcc_files,
  pkcFiles = inp$pkc,
  phenoDataFile = out$annotation_reordered,
  phenoDataSheet = sheet,
  phenoDataDccColName = dcc_col,
  protocolDataColNames = unlist(ann_cfg$protocol_data_colnames),
  experimentDataColNames = ann_cfg$experiment_data_colnames
)

# standR (stage 2 onwards) needs a column literally named "Segment".
seg_col <- ann_cfg$segment_column %||% "Segment"
if (!seg_col %in% colnames(pData(raw_data))) {
  stop("annotation.segment_column ('", seg_col, "') not found. Available columns: ",
       paste(colnames(pData(raw_data)), collapse = ", "))
}
pData(raw_data)$Segment <- pData(raw_data)[[seg_col]]

log_step("Raw data (loaded)", raw_data)
modules <- gsub(".pkc", "", annotation(raw_data))

## 3. Segment QC -----------------------------------------------------------------

data <- shiftCountsOne(raw_data, useDALogic = isTRUE(seg_qc$use_da_logic))

qc_cutoffs <- seg_qc[c("minSegmentReads", "percentTrimmed", "percentStitched", "percentAligned",
                       "percentSaturation", "minNegativeCount", "maxNTCCount", "minNuclei", "minArea")]

# Apply minArea only if the annotation has real (non-empty) area values;
# otherwise every segment would be flagged.
has_area <- "area" %in% colnames(pData(raw_data)) &&
  any(!is.na(suppressWarnings(as.numeric(pData(raw_data)$area))))
if (!has_area) {
  message("No usable 'area' column: minArea is not applied.")
  qc_cutoffs$minArea <- NULL
}

data <- setSegmentQCFlags(data, qcCutoffs = qc_cutoffs)

# NA flags mean "check not applicable", not "failed"
qc_flags <- protocolData(data)[["QCFlags"]]
qc_summary <- data.frame(Pass = colSums(!qc_flags, na.rm = TRUE),
                         Warning = colSums(qc_flags, na.rm = TRUE))
qc_status <- ifelse(rowSums(qc_flags, na.rm = TRUE) == 0, "PASS", "WARNING")
qc_summary["TOTAL FLAGS", ] <- c(sum(qc_status == "PASS"), sum(qc_status == "WARNING"))
write.table(qc_summary, out$qc_summary_table, sep = "\t", col.names = NA, quote = FALSE)

# Negative-probe geometric mean per segment and module (used for LOQ below)
neg_geo_means <- esBy(negativeControlSubset(data), GROUP = "Module",
                      FUN = function(x) assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs"))
protocolData(data)[["NegGeoMean"]] <- neg_geo_means

ntc_table <- as.data.frame(table(sData(data)$NTC))
colnames(ntc_table) <- c("NTC Count", "# of Segments")
write.table(ntc_table, out$ntc_summary, sep = "\t", row.names = FALSE, quote = FALSE)

if (isTRUE(seg_qc$remove_flagged_segments)) {
  data <- data[, qc_status == "PASS"]
}
log_step("After Segment QC", data)

# QC histograms (panels a-f of QC_final.pdf)
segment_colors <- unlist(plot_cfg$segment)
if (length(segment_colors) == 0) {
  lv <- sort(unique(as.character(pData(data)$Segment)))
  segment_colors <- setNames(scales::hue_pal()(length(lv)), lv)
}
col_by <- plot_cfg$col_by

qc_histogram <- function(column, threshold, title = column, log_x = FALSE) {
  # sData() stores some QC metrics as nested columns: unlist() gives a plain vector
  df <- data.frame(value = unlist(sData(data)[[column]]), group = sData(data)[[col_by]])
  p <- ggplot(df, aes(x = value, fill = group)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = threshold, lty = "dashed") +
    facet_wrap(~group, nrow = 4) +
    scale_fill_manual(values = segment_colors) +
    guides(fill = "none") +
    labs(x = title, y = "Segments, #", title = title) +
    theme_bw()
  if (log_x) p <- p + scale_x_continuous(trans = "log10")
  p
}

qc_panels <- list(
  qc_histogram("Trimmed (%)", seg_qc$percentTrimmed),
  qc_histogram("Stitched (%)", seg_qc$percentStitched),
  qc_histogram("Aligned (%)", seg_qc$percentAligned),
  qc_histogram("Saturated (%)", seg_qc$percentSaturation, title = "Sequencing Saturation (%)")
)
if (has_area) qc_panels <- c(qc_panels, list(qc_histogram("area", seg_qc$minArea, title = "Area", log_x = TRUE)))
qc_panels <- c(qc_panels, list(tableGrob(qc_summary, theme = ttheme_minimal())))
qc_grid <- plot_grid(plotlist = qc_panels, labels = "auto", ncol = 3)

## 4. Probe QC and gene-level aggregation -----------------------------------------

data <- setBioProbeQCFlags(
  data,
  qcCutoffs = list(minProbeRatio = probe_qc$minProbeRatio, percentFailGrubbs = probe_qc$percentFailGrubbs),
  removeLocalOutliers = isTRUE(probe_qc$removeLocalOutliers)
)
probe_flags <- fData(data)[["QCFlags"]]
data <- subset(data, !probe_flags$LowProbeRatio & !probe_flags$GlobalGrubbsOutlier)
log_step("After Probe QC", data)

target_data <- aggregateCounts(data)
log_step("After aggregateCounts (gene level)", target_data)

## 5. Limit of quantification (LOQ) ------------------------------------------------
# LOQ = max(minLOQ, NegGeoMean * NegGeoSD ^ cutoff), per segment and module

loq <- data.frame(row.names = colnames(target_data))
for (module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"), module)
  if (all(vars %in% colnames(pData(target_data)))) {
    loq[, module] <- pmax(loq_cfg$minLOQ,
                          pData(target_data)[, vars[1]] * pData(target_data)[, vars[2]] ^ loq_cfg$cutoff)
  }
}
pData(target_data)$LOQ <- loq

# Detection matrix: TRUE where a gene is above its segment's LOQ
loq_mat <- c()
for (module in modules) {
  in_module <- fData(target_data)$Module == module
  loq_mat <- rbind(loq_mat, t(esApply(target_data[in_module, ], MARGIN = 1, FUN = function(x) x > loq[, module])))
}
loq_mat <- loq_mat[fData(target_data)$TargetName, ]

pData(target_data)$GenesDetected <- colSums(loq_mat, na.rm = TRUE)
pData(target_data)$GeneDetectionRate <- pData(target_data)$GenesDetected / nrow(target_data)
pData(target_data)$DetectionThreshold <- cut(pData(target_data)$GeneDetectionRate,
                                             breaks = unlist(filt_cfg$detection_breaks),
                                             labels = unlist(filt_cfg$detection_labels))

p_detect_group <- ggplot(pData(target_data), aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = .data[[group_col]])) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate", y = "Samples, #", fill = group_col) +
  theme_bw()

p_detect_segment <- ggplot(pData(target_data), aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = Segment), color = "black") +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = segment_colors) +
  labs(x = "Gene Detection Rate", y = "Samples, #", fill = "Segment") +
  theme_classic()

## 6. Segment filter (minimum fraction of genes detected) ---------------------------

target_data <- target_data[, pData(target_data)$GeneDetectionRate >= filt_cfg$segment_detection_rate_min]
log_step("After segment detection-rate filter", target_data)

loq_mat <- loq_mat[, colnames(target_data)]
fData(target_data)$DetectedSegments <- rowSums(loq_mat, na.rm = TRUE)
fData(target_data)$DetectionRate <- fData(target_data)$DetectedSegments / ncol(target_data)

gene_detection <- data.frame(Freq = unlist(filt_cfg$detection_thresholds_pct))
gene_detection$Number <- sapply(unlist(filt_cfg$detection_thresholds_frac),
                                function(x) sum(fData(target_data)$DetectionRate >= x))
gene_detection$Rate <- gene_detection$Number / nrow(target_data)

p_gene_detection <- ggplot(gene_detection, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = " ")), vjust = -1, size = 4) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1), labels = scales::percent) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.1), expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments", y = "Genes Detected, % of Panel > LOQ") +
  theme_classic()

detection_grid <- plot_grid(p_detect_group, p_detect_segment, p_gene_detection,
                            labels = c("g", "h", "i"), ncol = 3, rel_widths = c(2.5, 2, 2))
pdf(out$qc_final_pdf, height = 12, width = 14)
print(plot_grid(qc_grid, detection_grid, ncol = 1, rel_heights = c(2, 0.9)))
dev.off()

## 7. Gene filter (minimum fraction of segments detected; negative probes always kept) --

neg_probes <- unique(fData(target_data)$TargetName[fData(target_data)$CodeClass == "Negative"])
target_data <- target_data[fData(target_data)$DetectionRate >= filt_cfg$gene_detection_rate_min |
                             fData(target_data)$TargetName %in% neg_probes, ]
log_step("After gene detection-rate filter", target_data)
write.table(exprs(target_data), out$data_raw, quote = FALSE)

## 8. Q3 normalization ----------------------------------------------------------------

norm_target_data <- normalize(target_data, norm_method = norm_cfg$method,
                              desiredQuantile = norm_cfg$desired_quantile, toElt = "q_norm")
log_step("After normalization", norm_target_data)
write.table(assayDataElement(norm_target_data, "q_norm"), out$data_q3, sep = "\t", quote = FALSE)

assayDataElement(norm_target_data, elt = "log_q") <-
  assayDataApply(norm_target_data, 2, FUN = log, base = 2, elt = "q_norm")

## 9. Save -------------------------------------------------------------------------------

write.table(qc_log, out$qc_log, sep = "\t", row.names = FALSE, quote = FALSE)
save.image(file = out$rdata_image)
saveRDS(target_data, out$target_rds)
saveRDS(norm_target_data, out$norm_target_rds)

message(sprintf("Preprocessing complete: %d genes x %d segments.", nrow(target_data), ncol(target_data)))
