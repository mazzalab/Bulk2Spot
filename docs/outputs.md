# Outputs

Everything is written under `outputdir`:

```
<outputdir>/
├── results/
│   ├── preprocessing/      stage 1
│   ├── spe/                stage 2
│   ├── batch_correction/   stage 3
│   ├── deg_gsea/           stage 4
│   ├── deconvolution/      stage 5
│   └── report/             stage 6: bulk2spot_report.html
└── logs/
    ├── <rule>.log          analysis log of each rule
    └── cluster/            scheduler stdout/stderr (cluster mode)
```

**Start with `results/report/bulk2spot_report.html`.** It is self-contained: open it in any browser,
or send it to collaborators.

To see real examples of these files, look in [`examples/tutorial_output/`](../examples/README.md),
which contains the tutorial run's outputs.

File types: `.rds` files are R objects (`readRDS()`), `.txt` files are tab-separated tables, and
`.xlsx` files are Excel tables.

## 1. `preprocessing/`

| File | Content |
|---|---|
| `QC_final.pdf` | **Main QC figure**: histograms of trimmed/stitched/aligned/saturation % (and area) per segment type with thresholds, the QC flag summary, and gene detection rates. |
| `qc_summary_table.txt` | Segments passing or flagged for each QC check, plus totals. |
| `ntc_count_summary.txt` | No-template-control read counts. |
| `QC_log.txt` | Genes x segments remaining after each step. |
| `annotation_reordered.xlsx` | The annotation sorted by DCC name (the version that was loaded). |
| `data_raw.txt` | Gene-level counts after filtering (genes x segments). |
| `data_Q3.txt` | Q3-normalized counts. |
| `target_data.rds` | Filtered gene-level `NanoStringGeoMxSet`. |
| `norm_target_data.rds` | Same object with `q_norm` and `log_q` assays. **Input to stages 2 and 5.** |
| `Preprocessing.RData` | The complete R session at the end of the stage, for interactive exploration. |

## 2. `spe/`

| File | Content |
|---|---|
| `spe.rds` | `SpatialExperiment` (standR) with added `Patient`, `Scan`, `biology` columns. |
| `exprs.txt` / `phenoData.txt` / `featureData.txt` | Expression matrix, segment annotation and gene annotation as flat tables. |

## 3. `batch_correction/`

| File | Content |
|---|---|
| `SampleInfo_Plot.pdf` | Alluvial overview: patient → segment → group. |
| `PCAScree_prebatch.pdf` | RLE plots, PCA coloured by group / scan / segment (/ batch), and scree plot, before correction. |
| `PCAScree_postbatch.pdf` | The same after correction (a placeholder if correction was skipped). |
| `RUV.pdf` | RUV4 diagnostics: silhouette score per k, RLE, pairwise PC1-PC4 plots (a placeholder if skipped). |
| `PCA_final.txt`, `PCA_final_variance.txt`, `UMAP_final.txt` | Final PCA and UMAP coordinates with metadata (used by the report). |
| `spe_ruv.rds` | TMM-normalized, and batch-corrected if correction ran, `SpatialExperiment`. **Input to stages 4 and 5.** |

## 4. `deg_gsea/`

| File | Content |
|---|---|
| `DEG_GSEA_summary.txt` | One row per comparison that was run: column, segment, groups, sample sizes, genes tested, up/down genes, number of enriched GO terms per ontology and KEGG pathways. |
| `DEG/<segment>/<comparison>/` | One folder per segment x comparison (below). |

Inside each `DEG/<segment>/<test>_vs_<reference>/` (`<tag>` = `<comparison>_<segment>`):

| File | Content |
|---|---|
| `TableDEG_<tag>.xlsx` | limma-voom results for all genes: `logFC` (test vs reference), `AveExpr`, `t`, `P.Value`, `adj.P.Val` (BH), `B`, and gene annotation. |
| `Volcano_<tag>.pdf` | Volcano plot with the top genes labelled. |
| `GSEA_<ontology>_results_<tag>.xlsx` | GO GSEA results (`NES`, `pvalue`, `p.adjust`, core-enrichment genes as symbols, one per column). |
| `Dotplot_<ontology>_<tag>.pdf` | Top up- and down-regulated GO terms. |
| `GSEA_KEGG_results_<tag>.xlsx` / `DotplotKEGG_<tag>.pdf` | The same for KEGG pathways. |

GSEA files are missing when no term passes the cutoffs, or when too few genes pass the GSEA gene
filters. The summary table shows the counts and the log explains why.

## 5. `deconvolution/`

| File | Content |
|---|---|
| `proportions.txt` | **Cell-type proportions per ROI** (fractions of the non-tumour signal, summing to 1) plus the `export_columns`. |
| `Barplot_ordered_no_batch.pdf` | Stacked proportions per ROI, ROIs ordered by clustering (plain CPM variant). |
| `Barplot_ordered.pdf`, `Barplot_by_Scan.pdf` | The same for the batch-corrected variant, and per scan (placeholders if that variant didn't run). |
| `res_dec_no_batch_correction.rds`, `res_dec_batch_corrected.rds` | Full `spatialdecon()` results: betas, proportions, p-values, residuals (`NULL` if not run). |
| `cibersort_signature.txt`, `cibersort_mixture.txt` | Inputs for the CIBERSORT(x) web tool. |
| `normCPM.rds`, `norm_counts_matrix.rds`, `fake_obj.rds` | Intermediate matrices and the GeoMxSet carrying them (used by the GSVA step). |
| `stats/MeanProportion_by_<column>.pdf` | Mean composition per group. |
| `stats/Proportion_perROI_by_<column>.pdf` | Composition of every ROI, faceted by group. |
| `stats/Stats_<column>.txt` | Wilcoxon / Kruskal-Wallis p-value per segment x cell type. |
| `stats/Celltype_significant_boxplots_<column>.pdf` | Boxplots of the hits with p < `alpha` (only if there are any). |
| `SegmentPair_comparison.pdf` | Two-segment comparison of summed cell-type groups (a placeholder if disabled). |
| `gsva_scores.rds` | GSVA score matrix (signatures x ROIs). |
| `gsva_long.txt` / `.rds` | Scores in long format, with ROI annotation. |
| `GSVA_Heatmap.pdf` | Heatmap of the scores, annotated by segment/group. |
| `GSVA_Segment_boxplot.pdf` | Scores per segment with BH-adjusted pairwise tests. |
| `GSVA_stats/Stats_<column>.txt`, `GSVA_stats/Signature_significant_boxplots_<column>.pdf` | Per-group tests on the GSVA scores. |

## 6. `report/`

`bulk2spot_report.html` has these sections: overview numbers; preprocessing QC (segment/gene funnels,
QC flags, NTC table); sample composition; batch correction (overview PCA before correction, final
PCA/UMAP); DEG & GSEA per segment x comparison (interactive volcano plots, top-gene tables, GO BP and
KEGG plots and tables); deconvolution (composition, significance tables, GSVA heatmap); methods; and a
glossary.
