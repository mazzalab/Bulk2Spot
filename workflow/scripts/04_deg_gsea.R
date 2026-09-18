# Stage 4 -- Differential expression (limma-voom) and GSEA (clusterProfiler)
#   for every comparison column x segment x pair of groups:
#   DEG table + volcano plot + GO (per ontology) and KEGG GSEA tables/dotplots
#
# Outputs: DEG/<segment>/<comparison>/... and one DEG_GSEA_summary.txt

snakemake@source("common.R")
start_logging(snakemake@log[[1]])

suppressPackageStartupMessages({
  library(standR)
  library(SummarizedExperiment)
  library(limma)
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(openxlsx)
})

out_summary <- snakemake@output$summary_file
out_dir <- snakemake@output$deg_dir
cfg <- snakemake@config
cmp_cfg <- cfg$comparisons
thr <- cfg$thresholds
enr_cfg <- cfg$enrichment
temp_col <- cfg$design$temp_column_name

library(enr_cfg$organism_db, character.only = TRUE)
org_db <- get(enr_cfg$organism_db)

# GSEA p-values are permutation-based: fix the seed for reproducible results
set.seed(cfg$dimensionality_reduction$seed)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

spe_ruv <- readRDS(snakemake@input$spe_ruv_rds)
dge_all <- spe2dge(spe_ruv)   # keeps the TMM normalization factors
dge_all$samples <- dge_all$samples[, !grepl("\\.1$", colnames(dge_all$samples))]

## Helpers ------------------------------------------------------------------------------

# makeContrasts() parses level names as R code: replace spaces and dashes
sanitize <- function(x) gsub("[ -]", "_", x)

# Pairs of groups to compare for one column: all pairs, or the manual list
make_comparisons <- function(column) {
  if (identical(cmp_cfg$mode, "manual")) {
    return(lapply(cmp_cfg$manual_pairs, unlist))
  }
  lv <- sort(unique(as.character(colData(spe_ruv)[[column]])))
  if (length(unlist(cmp_cfg$include_levels)) > 0) lv <- intersect(lv, unlist(cmp_cfg$include_levels))
  lv <- setdiff(lv, unlist(cmp_cfg$exclude_levels))
  if (length(lv) < 2) stop("Column '", column, "' needs at least 2 levels to compare (found: ", toString(lv), ")")
  pairs <- combn(lv, 2, simplify = FALSE)
  setNames(pairs, sapply(pairs, paste, collapse = "_vs_"))
}

# Top-N up and top-N down pathways passing the adjusted p-value cutoff
top_pathways <- function(res) {
  sig <- res[res$p.adjust < thr$pathway_padj_cutoff, ]
  up <- sig[sig$NES > 0, ]
  down <- sig[sig$NES < 0, ]
  tab <- rbind(head(up[order(-up$NES), ], thr$top_n_pathways),
               head(down[order(down$NES), ], thr$top_n_pathways))
  tab[complete.cases(tab), ]
}

pathway_dotplot <- function(tab, title, ylab, file) {
  if (nrow(tab) == 0) return(invisible())
  width <- 8 + max(0, (max(nchar(tab$Description)) - 40) / 20)   # room for long names
  p <- ggplot(tab, aes(x = NES, y = reorder(Description, NES), size = NES, color = p.adjust)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_x_continuous(limits = range(tab$NES) + c(-0.5, 0.5)) +
    scale_color_continuous(low = "red", high = "blue") +
    labs(title = title, x = "NES", y = ylab) +
    theme_bw()
  ggsave(file, p, width = width, height = 8)
}

# GSEA table with core-enrichment Entrez IDs translated to gene symbols, one per column
expand_core_enrichment <- function(res) {
  ids <- unique(unlist(strsplit(res$core_enrichment, "/")))
  id2symbol <- with(bitr(ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org_db), setNames(SYMBOL, ENTREZID))
  genes <- lapply(strsplit(res$core_enrichment, "/"), function(x) unname(na.omit(id2symbol[x])))
  n <- max(lengths(genes))
  gene_cols <- do.call(rbind, lapply(genes, function(g) c(g, rep(NA, n - length(g)))))
  colnames(gene_cols) <- paste0("Gene_", seq_len(n))
  cbind(res[, setdiff(colnames(res), "core_enrichment")],
        core_enrichment_symbol = sapply(genes, paste, collapse = "/"), gene_cols)
}

## Main loop: comparison column x segment x comparison ------------------------------------

summary_rows <- list()

for (column in unlist(cmp_cfg$column)) {
  comparisons <- make_comparisons(column)
  message("Column '", column, "': ", length(comparisons), " comparison(s): ", toString(names(comparisons)))
  spe_ruv[[temp_col]] <- spe_ruv[[column]]

  for (segment in unlist(cfg$segments)) {
    spe_seg <- spe_ruv[, spe_ruv$Segment %in% segment]
    dge <- dge_all[, dge_all$samples$Segment %in% segment]

    for (cmp in names(comparisons)) {
      group_test <- comparisons[[cmp]][1]
      group_ref <- comparisons[[cmp]][2]
      n_test <- sum(spe_seg[[column]] %in% group_test)
      n_ref <- sum(spe_seg[[column]] %in% group_ref)
      message("Segment '", segment, "' | ", cmp, " (n = ", n_test, " vs ", n_ref, ")")

      # Both groups need enough samples in this segment
      if (n_test < cmp_cfg$min_samples_per_group || n_ref < cmp_cfg$min_samples_per_group) {
        message("  -> fewer than ", cmp_cfg$min_samples_per_group, " samples in a group, skipped.")
        next
      }
      cmp_dir <- file.path(out_dir, segment, cmp)
      dir.create(cmp_dir, showWarnings = FALSE, recursive = TRUE)
      tag <- paste0(cmp, "_", segment)

      ## limma-voom -------------------------------------------------------------------
      covariates <- ""
      if (isTRUE(cfg$design$include_ruv_factors)) {
        ruv_cols <- intersect(c("ruv_W1", "ruv_W2", "ruv_W3"), colnames(colData(spe_seg)))
        if (length(ruv_cols) > 0) covariates <- paste0(" + ", paste(ruv_cols, collapse = " + "))
      }
      design <- model.matrix(as.formula(paste0("~0 + ", temp_col, covariates)), data = colData(spe_seg))
      colnames(design) <- sanitize(gsub(temp_col, "", colnames(design)))

      contrast <- paste0(sanitize(group_test), " - ", sanitize(group_ref))
      contr_matrix <- makeContrasts(contrasts = contrast, levels = colnames(design))

      dge_cmp <- estimateDisp(dge, design = design, robust = TRUE)
      fit <- lmFit(voom(dge_cmp, design), design)
      efit <- eBayes(contrasts.fit(fit, contrasts = contr_matrix), robust = TRUE)

      deg <- topTable(efit, coef = contrast, sort.by = "P", n = Inf)
      write.xlsx(deg, file.path(cmp_dir, paste0("TableDEG_", tag, ".xlsx")), rowNames = TRUE)

      is_sig <- deg$adj.P.Val < thr$p_cutoff
      up <- deg[is_sig & deg$logFC > thr$fc_soft, ]
      down <- deg[is_sig & deg$logFC < -thr$fc_soft, ]

      ## Volcano plot (top N up/down genes labelled) -------------------------------------
      labelled <- c(head(rownames(up)[order(-up$logFC)], thr$top_genes_volcano),
                    head(rownames(down)[order(down$logFC)], thr$top_genes_volcano))
      volcano_df <- deg %>%
        dplyr::mutate(gene = rownames(deg), neg_log10_fdr = -log10(adj.P.Val),
                      label = ifelse(gene %in% labelled, gene, ""))
      volcano <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_fdr)) +
        geom_point(aes(color = logFC, alpha = neg_log10_fdr)) +
        scale_color_gradient2(low = "#20854EFF", mid = "white", high = "#EFC000FF", midpoint = 0) +
        scale_alpha_continuous(range = c(0.75, 1), guide = "none") +
        geom_text_repel(aes(label = label), size = 3, max.overlaps = 100) +
        geom_vline(xintercept = c(-thr$fc_soft, thr$fc_soft), linetype = "dashed") +
        geom_vline(xintercept = c(-thr$fc_hard, thr$fc_hard), linetype = "dotted") +
        geom_hline(yintercept = -log10(thr$p_cutoff), linetype = "dashed") +
        xlim(unlist(cfg$plot$volcano_xlim)) + ylim(unlist(cfg$plot$volcano_ylim)) +
        labs(x = expression(log[2] * " FC"), y = expression(-log[10] * " FDR")) +
        theme_classic(base_size = 13) + theme(legend.position = "none")
      ggsave(file.path(cmp_dir, paste0("Volcano_", tag, ".pdf")), volcano, width = 7, height = 5)

      row <- data.frame(Column = column, Segment = segment, Comparison = cmp,
                        Group_Test = group_test, Group_Ref = group_ref, N_Test = n_test, N_Ref = n_ref,
                        N_tested = nrow(deg), N_up = nrow(up), N_down = nrow(down))

      ## GSEA: genes with P < filter and |logFC| > fc_soft, ranked by logFC ------------------
      # One entry per Entrez ID (the probe with the largest |logFC|)
      gene_list <- deg %>%
        dplyr::filter(P.Value < thr$gsea_gene_pvalue_filter, abs(logFC) > thr$fc_soft, !is.na(GeneID)) %>%
        dplyr::arrange(desc(abs(logFC))) %>%
        dplyr::distinct(GeneID, .keep_all = TRUE) %>%
        dplyr::arrange(desc(logFC))
      gene_list <- setNames(gene_list$logFC, gene_list$GeneID)

      if (length(gene_list) < 2) {
        message("  -> not enough genes for GSEA, skipped.")
      } else {
        for (ont in unlist(enr_cfg$ontology)) {
          go <- gseGO(geneList = gene_list, OrgDb = org_db, keyType = "ENTREZID", ont = ont,
                      pAdjustMethod = "BH", pvalueCutoff = thr$gsea_pvalue_cutoff,
                      verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
          res <- if (is.null(go)) data.frame() else go@result
          row[[paste0("N_GO_", ont, "_terms")]] <- nrow(res)
          if (nrow(res) > 0) {
            res$Description <- str_to_title(res$Description)
            pathway_dotplot(top_pathways(res), paste("GO", ont), "GO Term",
                            file.path(cmp_dir, paste0("Dotplot_", ont, "_", tag, ".pdf")))
            write.xlsx(expand_core_enrichment(res), file.path(cmp_dir, paste0("GSEA_", ont, "_results_", tag, ".xlsx")),
                       sheetName = paste0("GSEA_", ont), overwrite = TRUE)
          }
        }

        # KEGG needs internet access (downloads pathway data from KEGG)
        kegg <- tryCatch(
          gseKEGG(geneList = gene_list, organism = enr_cfg$kegg_organism, pAdjustMethod = "BH",
                  pvalueCutoff = thr$gsea_pvalue_cutoff, verbose = FALSE, BPPARAM = BiocParallel::SerialParam()),
          error = function(e) { message("  -> KEGG unavailable: ", conditionMessage(e)); NULL })
        res <- if (is.null(kegg)) data.frame() else kegg@result
        row$N_KEGG_terms <- if (is.null(kegg)) NA else nrow(res)
        if (nrow(res) > 0) {
          res$Description <- str_to_title(res$Description)
          pathway_dotplot(top_pathways(res), "KEGG", "KEGG Pathway",
                          file.path(cmp_dir, paste0("DotplotKEGG_", tag, ".pdf")))
          write.xlsx(expand_core_enrichment(res), file.path(cmp_dir, paste0("GSEA_KEGG_results_", tag, ".xlsx")),
                     sheetName = "GSEA_KEGG", overwrite = TRUE)
        }
      }
      summary_rows[[length(summary_rows) + 1]] <- row
    }
  }
}

## Summary (one row per comparison actually run) ---------------------------------------------

deg_summary <- if (length(summary_rows) > 0) bind_rows(summary_rows) else
  data.frame(Column = character(), Segment = character(), Comparison = character())
write.table(deg_summary, out_summary, sep = "\t", row.names = FALSE, quote = FALSE)
message("Saved DEG/GSEA summary to ", out_summary)
