# Small helpers shared by the Bulk2Spot R scripts.
# Load from a Snakemake script with: snakemake@source("common.R")

# Default value for optional config keys: x %||% y returns y when x is NULL.
`%||%` <- function(x, y) if (is.null(x)) y else x

# Send everything a script prints (output, messages, warnings) to the rule's log file.
start_logging <- function(log_path) {
  con <- file(log_path, open = "wt")
  sink(con)
  sink(con, type = "message")
}

# One-page PDF explaining why a plot was skipped. Snakemake requires every
# declared output to exist, even when a step is switched off.
placeholder_pdf <- function(path, text) {
  pdf(path, width = 8, height = 4)
  plot.new()
  text(0.5, 0.5, text, cex = 1.1)
  dev.off()
}

# Drop non-atomic columns (e.g. GeomxTools' per-module LOQ data.frame), which
# would corrupt a flat TSV export.
flatten_df <- function(df) {
  df[, vapply(df, is.atomic, logical(1)), drop = FALSE]
}

# For each segment x feature, compare `value_col` across the levels of `group_col`:
# Wilcoxon rank-sum test for 2 levels, Kruskal-Wallis for more than 2.
# Returns a data.frame with columns <segment_col>, <feature_col>, pvalue.
test_by_group <- function(df, group_col, segment_col, feature_col, value_col) {
  rows <- list()
  for (seg in na.omit(unique(df[[segment_col]]))) {
    for (feat in na.omit(unique(df[[feature_col]]))) {
      tmp <- df[df[[segment_col]] %in% seg & df[[feature_col]] %in% feat, ]
      lv <- unique(na.omit(tmp[[group_col]]))
      if (length(lv) < 2) next
      pval <- tryCatch({
        if (length(lv) == 2) {
          wilcox.test(tmp[[value_col]][tmp[[group_col]] %in% lv[1]],
                      tmp[[value_col]][tmp[[group_col]] %in% lv[2]])$p.value
        } else {
          kruskal.test(tmp[[value_col]], tmp[[group_col]])$p.value
        }
      }, error = function(e) NA)
      rows[[length(rows) + 1]] <- data.frame(seg, feat, pval)
    }
  }
  res <- if (length(rows) > 0) do.call(rbind, rows) else data.frame(character(), character(), numeric())
  setNames(res, c(segment_col, feature_col, "pvalue"))
}

# Boxplots of the significant segment x feature hits, `per_page` plots per PDF page.
paged_boxplots <- function(df, hits, group_col, segment_col, feature_col, value_col,
                           ylab, ylim, out_file, per_page) {
  plots <- lapply(seq_len(nrow(hits)), function(i) {
    tmp <- df[df[[segment_col]] %in% hits[[segment_col]][i] & df[[feature_col]] %in% hits[[feature_col]][i], ]
    ggplot2::ggplot(tmp, ggplot2::aes(x = .data[[group_col]], y = .data[[value_col]], fill = .data[[group_col]])) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5) +
      ggplot2::geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
      ggplot2::labs(x = group_col, y = ylab,
                    title = paste0("Segment: ", hits[[segment_col]][i], " | ", feature_col, ": ",
                                   hits[[feature_col]][i], "\np = ", signif(hits$pvalue[i], 3))) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::coord_cartesian(ylim = ylim)
  })
  # Fixed grid size, so pages with fewer plots keep the same panel size
  n_pages <- ceiling(length(plots) / per_page)
  pdf(out_file, width = 14, height = 10)
  for (page in seq_len(n_pages)) {
    idx <- ((page - 1) * per_page + 1):min(page * per_page, length(plots))
    print(patchwork::wrap_plots(plots[idx], ncol = 2, nrow = ceiling(per_page / 2)))
  }
  dev.off()
  n_pages
}
