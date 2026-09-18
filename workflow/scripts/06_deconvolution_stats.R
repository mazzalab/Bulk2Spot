# Stage 5b -- Statistics on the deconvolution proportions
#   per grouping column: mean-composition barplots, per-ROI barplots,
#   per segment x cell type tests (Wilcoxon / Kruskal-Wallis) + boxplots of significant hits;
#   optional comparison of two segments on summed "meta cell types".

snakemake@source("common.R")
start_logging(snakemake@log[[1]])

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

out <- snakemake@output
cfg <- snakemake@config
stats_cfg <- cfg$deconvolution_stats
dec_cfg <- cfg$deconvolution
dir.create(out$stats_dir, showWarnings = FALSE, recursive = TRUE)

## 1. Proportions in long format (metadata = ROI + exported columns; the rest are cell types) --

prop <- read.delim(snakemake@input$proportions_file, check.names = FALSE)
id_cols <- intersect(c("ROI", names(dec_cfg$export_columns)), colnames(prop))
prop_long <- pivot_longer(prop, cols = -all_of(id_cols), names_to = "celltype", values_to = "proportion")

# Optional custom colours per cell type (missing ones in grey)
fill_scale <- function(cell_types) {
  if (!isTRUE(dec_cfg$use_custom_palette)) return(NULL)
  cols <- unlist(dec_cfg$custom_palette)[cell_types]
  scale_fill_manual(values = setNames(ifelse(is.na(cols), "grey70", cols), cell_types))
}

## 2. Composition barplots per grouping column ------------------------------------------------

group_cols <- intersect(unlist(stats_cfg$grouping_columns), colnames(prop))

for (gcol in group_cols) {
  mean_prop <- prop_long %>%
    group_by(.data[[gcol]], celltype) %>%
    summarise(mean_prop = mean(proportion, na.rm = TRUE), .groups = "drop")

  p_mean <- ggplot(mean_prop, aes(x = .data[[gcol]], y = mean_prop, fill = celltype)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = gcol, y = "Mean proportion (within non-tumor)",
         title = paste("Cell-type composition by", gcol, "(SpatialDecon)")) +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    fill_scale(unique(mean_prop$celltype))
  ggsave(file.path(out$stats_dir, paste0("MeanProportion_by_", gcol, ".pdf")), p_mean, width = 7, height = 7)

  p_roi <- ggplot(prop_long, aes(x = ROI, y = proportion, fill = celltype)) +
    geom_bar(stat = "identity", color = "black", width = 0.8) +
    facet_grid(as.formula(paste("~", gcol)), scales = "free_x") +
    labs(x = "ROI", y = "Proportion") +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.text = element_text(size = 12, face = "bold")) +
    fill_scale(unique(prop_long$celltype))
  ggsave(file.path(out$stats_dir, paste0("Proportion_perROI_by_", gcol, ".pdf")), p_roi, width = 30, height = 8)
}

## 3. Significance tests per segment x cell type -----------------------------------------------
# Needs a `segment` entry in deconvolution.export_columns

if (!"segment" %in% colnames(prop_long)) {
  message("No 'segment' column in proportions.txt: significance tests skipped.")
} else {
  for (gcol in group_cols) {
    tests <- test_by_group(prop_long, gcol, "segment", "celltype", "proportion")
    write.table(tests, file.path(out$stats_dir, paste0("Stats_", gcol, ".txt")), sep = "\t", row.names = FALSE, quote = FALSE)

    hits <- tests[!is.na(tests$pvalue) & tests$pvalue < stats_cfg$alpha, ]
    message(gcol, ": ", nrow(hits), " significant comparison(s) (p < ", stats_cfg$alpha, ")")
    if (nrow(hits) > 0) {
      paged_boxplots(prop_long, hits, gcol, "segment", "celltype", "proportion", ylab = "Proportion", ylim = c(0, 1),
                     out_file = file.path(out$stats_dir, paste0("Celltype_significant_boxplots_", gcol, ".pdf")),
                     per_page = stats_cfg$plots_per_page)
    }
  }
}

## 4. Optional: two segments compared on summed cell-type groups ---------------------------------

pair_cfg <- stats_cfg$segment_pair_comparison
if (isTRUE(pair_cfg$enabled)) {
  segs <- unlist(pair_cfg$segments)
  prop_pair <- prop %>% filter(segment %in% segs) %>% mutate(segment = factor(segment, levels = segs))
  if (length(unique(prop_pair$segment)) != 2) {
    stop("segment_pair_comparison.segments must name exactly 2 segments present in the data. Requested: ",
         toString(segs), "; available: ", toString(unique(prop$segment)))
  }

  plots <- list()
  for (group in names(pair_cfg$cell_type_groups)) {
    cols <- intersect(unlist(pair_cfg$cell_type_groups[[group]]), colnames(prop_pair))
    if (length(cols) == 0) {
      warning("No cell types found for group '", group, "'")
      next
    }
    prop_pair[[group]] <- rowSums(prop_pair[, cols, drop = FALSE], na.rm = TRUE)
    alternative <- pair_cfg$test_alternative[[group]] %||% "two.sided"
    test <- wilcox.test(as.formula(paste0("`", group, "` ~ segment")), data = prop_pair, alternative = alternative)
    y_max <- max(prop_pair[[group]], na.rm = TRUE) * 1.1
    if (!is.finite(y_max) || y_max == 0) y_max <- 1

    plots[[group]] <- ggplot(prop_pair, aes(x = segment, y = .data[[group]], fill = segment)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.5) +
      geom_jitter(width = 0.1, size = 2, alpha = 0.7, color = "darkgrey") +
      ylim(0, y_max) +
      labs(x = "", y = paste("Prop. of", group),
           title = paste0(group, ": ", segs[1], " vs ", segs[2], "\np = ", signif(test$p.value, 3))) +
      theme_classic() + theme(legend.position = "none", plot.title = element_text(size = 12, face = "bold"))
  }
  if (length(plots) > 0) {
    ggsave(out$segment_pair_plot, wrap_plots(plots, ncol = 2), width = 7, height = 7)
  } else {
    placeholder_pdf(out$segment_pair_plot, "Segment-pair comparison: no valid cell_type_groups.")
  }
} else {
  placeholder_pdf(out$segment_pair_plot, "Segment-pair comparison skipped:\ndeconvolution_stats.segment_pair_comparison.enabled = false")
}
