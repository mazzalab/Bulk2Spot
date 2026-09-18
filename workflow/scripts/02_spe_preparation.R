# Stage 2 -- SpatialExperiment conversion (standR)
#   normalized GeoMxSet -> SpatialExperiment, plus Patient / Scan / biology annotations

snakemake@source("common.R")
start_logging(snakemake@log[[1]])

suppressPackageStartupMessages({
  library(GeomxTools)
  library(dplyr)
  library(tibble)
  library(standR)
  library(SpatialExperiment)
  library(edgeR)
})

out <- snakemake@output
spe_cfg <- snakemake@config$spe

norm_target_data <- readRDS(snakemake@input$norm_target_rds)

## 1. Flat exports ----------------------------------------------------------------

write.table(exprs(norm_target_data), out$exprs_txt, sep = "\t", quote = FALSE)
write.table(flatten_df(pData(norm_target_data)), out$phenoData_txt, sep = "\t", quote = FALSE)
write.table(fData(norm_target_data), out$featureData_txt, sep = "\t", quote = FALSE, row.names = FALSE)

## 2. Build the three tables standR::readGeoMx() expects ---------------------------

counts <- as.data.frame(exprs(norm_target_data)) %>% rownames_to_column("TargetName")
sample_anno <- as.data.frame(pData(norm_target_data)) %>% rownames_to_column("SegmentDisplayName")
feature_anno <- as.data.frame(fData(norm_target_data))
rownames(feature_anno) <- NULL

# Spatial coordinates: real ones if configured, otherwise sequential placeholders
x_col <- spe_cfg$coordinate_x_column
y_col <- spe_cfg$coordinate_y_column
if (isTRUE(spe_cfg$use_real_coordinates) && all(c(x_col, y_col) %in% colnames(sample_anno))) {
  sample_anno$ROICoordinateX <- sample_anno[[x_col]]
  sample_anno$ROICoordinateY <- sample_anno[[y_col]]
} else {
  sample_anno$ROICoordinateX <- seq_len(nrow(sample_anno))
  sample_anno$ROICoordinateY <- seq_len(nrow(sample_anno))
}

# standR requires a specific name for the negative-control probe
neg_target <- spe_cfg$negative_probe_target_name
counts$TargetName[counts$TargetName == spe_cfg$negative_probe_source_name] <- neg_target
feature_anno$TargetName[feature_anno$Negative == "TRUE"] <- neg_target

## 3. Convert -----------------------------------------------------------------------

spe <- readGeoMx(counts, sample_anno, feature_anno)
message("Negative-control probes: ", sum(rowData(spe)$CodeClass == "Negative"))
spe <- readGeoMxFromDGE(SE2DGEList(spe))
message(sprintf("SpatialExperiment: %d features x %d segments", nrow(spe), ncol(spe)))

## 4. Extra annotations ---------------------------------------------------------------

# Anonymized labels: Scan_01, Scan_02 ... and <prefix>_01, <prefix>_02 ... (in order of appearance)
index_labels <- function(values, prefix) {
  ids <- unique(values)
  paste0(prefix, "_", sprintf("%02d", match(values, ids)))
}
spe$Scan <- index_labels(colData(spe)[[spe_cfg$scan_id_column]], "Scan")
spe$Patient <- index_labels(colData(spe)[[spe_cfg$patient_source_column]], spe_cfg$patient_label_prefix %||% "Pt")

# "biology" = pasted biological columns; preserved during batch correction
biology <- as.data.frame(colData(spe))[, unlist(spe_cfg$biology_columns), drop = FALSE]
spe$biology <- do.call(paste, c(as.list(biology), sep = spe_cfg$biology_separator))

saveRDS(spe, out$spe_rds)
message("Saved SpatialExperiment to ", out$spe_rds)
