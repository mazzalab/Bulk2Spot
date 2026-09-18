# Configuration

Every project uses **one YAML file**. Start from [`config/config.yaml`](../config/config.yaml), which
documents each option inline, or from the tutorial config for a filled-in example.

- **Paths.** Absolute paths are used as they are. A relative path is looked up from the directory
  you launch `run.py` from, and then from the Bulk2Spot repository (this is how bundled resources such
  as `resources/signatures/...` are found). `outputdir` is always relative to the launch directory.
- **Column names** must match your annotation workbook exactly, including case.
- Sections are listed below in the order the pipeline uses them.

## General

| Key | Type | Description |
|---|---|---|
| `outputdir` | path | Where `results/` and `logs/` are written. |
| `project_name` | string | Shown in the report and in notification e-mails. |
| `notify_email` | string, optional | E-mail sent on success or failure (needs system `mail`). |
| `conda_envs.analysis` | path/name or `""` | Pre-built analysis environment. `""` = build from `workflow/envs/bulk2spot.yaml`. |
| `conda_envs.report` | path/name or `""` | Pre-built report environment. `""` = build from `workflow/envs/report.yaml`. |

## 1. Raw data

| Key | Description |
|---|---|
| `dcc_dir` | Folder with the `.dcc` files (searched recursively). |
| `pkc_file` | The `.pkc` probe kit file. If you received a zip, unzip it first. |
| `annotation_file` / `annotation_sheet` | Excel workbook with one row per segment, and the sheet to read. |
| `exclude_dcc_files` | DCC file names (without extension) to leave out, e.g. known failed segments. |
| `annotation.pheno_data_dcc_colname` | Column holding each segment's DCC file name (usually `Sample_ID`). |
| `annotation.protocol_data_colnames` | Technical columns (ROI/AOI identifiers, region ...) stored as protocol data. |
| `annotation.experiment_data_colnames` | Experiment-wide columns, or `null`. |
| `annotation.segment_column` | Column with the segment/AOI type (e.g. `segment`). It is copied to a column named `Segment`, which later stages rely on. |

The annotation workbook must also contain the columns that later sections refer to:
the comparison column(s), `spe.patient_source_column`, `spe.scan_id_column`, `spe.biology_columns`,
and `batch_correction.batch_column` (if you use batch correction). If you add an `area` column with
values, `segment_qc.minArea` is applied.

## 2. Segment QC

Passed to `GeomxTools::setSegmentQCFlags()`. A segment is flagged when it fails any check.

| Key | Default | Flag raised when |
|---|---|---|
| `use_da_logic` | `true` | (not a check) `shiftCountsOne(useDALogic=)`: how zero counts are shifted to 1. |
| `minSegmentReads` | 1000 | raw reads < value |
| `percentTrimmed` / `percentStitched` / `percentAligned` | 80 | % of reads trimmed / stitched / aligned < value |
| `percentSaturation` | 50 | sequencing saturation (%) < value |
| `minNegativeCount` | 10 | geometric mean of negative probes < value |
| `maxNTCCount` | 1000 | no-template-control reads > value |
| `minNuclei` | 100 | nuclei count < value |
| `minArea` | 5000 | AOI area < value (only if the annotation has area values) |
| `remove_flagged_segments` | `false` | If `true`, flagged segments are removed. If `false`, they are only reported. |

> These defaults are NanoString's suggested starting points. Look at the histograms in
> `QC_final.pdf` and `qc_summary_table.txt` before accepting them. Low-expressing tissues and WTA
> panels often need lower `minNegativeCount`; the tutorial dataset needs `2` / `maxNTCCount: 9000`.

## 3. Probe QC

| Key | Default | Description |
|---|---|---|
| `minProbeRatio` | 0.1 | Probes whose geometric mean count, divided by that of all probes of the same target, is below this are removed. |
| `percentFailGrubbs` | 20 | Probes that are Grubbs outliers in more than this % of segments are removed. |
| `removeLocalOutliers` | `false` | Also remove per-segment outlier probe counts. |

## 4. Limit of quantification

`LOQ = max(minLOQ, NegGeoMean x NegGeoSD ^ cutoff)`, computed per segment and panel module.
A gene counts as detected in a segment when its count is above that segment's LOQ.

| Key | Default | Description |
|---|---|---|
| `cutoff` | 2 | Number of geometric standard deviations above background. |
| `minLOQ` | 2 | Lower bound for the LOQ. |

## 5. Segment and gene filtering

| Key | Default | Description |
|---|---|---|
| `segment_detection_rate_min` | 0 | Keep segments where at least this fraction of genes is detected (e.g. `0.05`). |
| `gene_detection_rate_min` | 0.1 | Keep genes detected in at least this fraction of segments. Negative probes are always kept. |
| `detection_thresholds_pct` / `_frac` | | Thresholds drawn in the "genes detected" plot (same values in % and as fractions). |
| `detection_breaks` / `detection_labels` | | Bins for the per-segment detection-rate plot (one more break than labels). |

## 6. Q3 normalization

| Key | Default | Description |
|---|---|---|
| `method` | `"quant"` | GeomxTools normalization method (`quant` = quantile/upper-quartile). |
| `desired_quantile` | 0.75 | Quantile used for scaling (0.75 = Q3). |

Q3-normalized values are exported (`data_Q3.txt`) and used for the report's overview PCA. **DEG and
deconvolution use the TMM normalization from section 8 instead.**

## 7. SpatialExperiment conversion

| Key | Description |
|---|---|
| `negative_probe_source_name` / `negative_probe_target_name` | The negative-control probe's name in GeomxTools (`"Negative Probe"`) and the name standR requires (`"NegProbe-WTX"`). |
| `patient_source_column` | Column that identifies the patient/subject. It becomes anonymized labels `Pt_01, Pt_02, ...` in `Patient`. |
| `scan_id_column` | Column that identifies the slide scan. It becomes `Scan_01, ...` in `Scan`. It can be the same column as the patient. |
| `patient_label_prefix` | Prefix of the patient labels (default `Pt`). |
| `biology_columns` / `biology_separator` | Columns pasted together into `biology` (e.g. `Group_Segment`). |
| `use_real_coordinates`, `coordinate_x_column`, `coordinate_y_column` | Use real ROI coordinates from the annotation. Otherwise placeholders 1, 2, 3 ... are used (GeoMx analyses don't need real coordinates). |

Column names containing spaces become dotted after conversion (`slide name` → `slide.name`). Use the
dotted form in this section.

## 8. Diagnostics and batch correction

| Key | Default | Description |
|---|---|---|
| `normalization_standR.method` | `"TMM"` | standR normalization (`TMM`, `upperquartile`, ...). This is the normalization DEG and deconvolution use. |
| `dimensionality_reduction.seed` | 100 | Random seed for PCA, UMAP and the GSEA permutations. |
| `batch_correction.perform` | `true` | Run RUV4 batch correction. It only runs if `batch_column` also exists in the data. |
| `batch_correction.batch_column` | `"Batch"` | Column that defines the batches. |
| `batch_correction.factor_of_interest` | `"Group"` | Biological variable that RUV4 must preserve. |
| `batch_correction.sample_info_first_column` / `sample_info_column` | `Patient` / `Segment` | First two strata of the alluvial sample overview (the third is the comparison column). |
| `batch_correction.ncg_top_n` | 100 | Number of negative control genes (least variable across batches). |
| `batch_correction.k` | `"auto"` | Number of unwanted factors: `"auto"` picks the best silhouette score, or give an integer. |
| `batch_correction.apply_limma_residual_correction` | `true` | Also run `limma::removeBatchEffect()` on the corrected log-counts. |

> **Only correct for a batch that is independent of your biology.** If every batch contains a single
> condition, RUV4 would remove the biological signal. The tutorial turns correction off for this
> reason.

## 9. Differential expression and GSEA

| Key | Default | Description |
|---|---|---|
| `comparisons.column` | `["Group"]` | Column(s) defining the groups. A string or a list. The first one also colours plots. |
| `comparisons.mode` | `"all_pairs"` | `all_pairs`: every pair of levels. `manual`: only `manual_pairs`. |
| `comparisons.include_levels` / `exclude_levels` | `[]` | Restrict the levels used to build pairs. |
| `comparisons.min_samples_per_group` | 2 | Comparisons with fewer segments than this in either group (within a segment type) are skipped. |
| `comparisons.manual_pairs` | `{}` | `name: [test_group, reference_group]`. A positive logFC means higher in the test group. |
| `segments` | | `Segment` values to analyse. Each is analysed separately. |
| `design.temp_column_name` | `"coltest"` | Internal column name used in the model formula (change it only if it clashes). |
| `design.include_ruv_factors` | `false` | Add RUV factors `ruv_W1..3` as covariates in the model. |
| `thresholds.p_cutoff` | 0.1 | Adjusted p-value cutoff for significant genes. |
| `thresholds.fc_soft` | 0.57 | \|log2FC\| cutoff for significant genes and for genes entering GSEA (0.57 ≈ 1.5-fold). |
| `thresholds.fc_hard` | 1.5 | Extra dotted guide line on volcano plots. |
| `thresholds.top_genes_volcano` | 15 | Up and down genes labelled on volcano plots. |
| `thresholds.gsea_gene_pvalue_filter` | 0.05 | Raw p-value filter for genes entering GSEA. |
| `thresholds.gsea_pvalue_cutoff` | 0.05 | `pvalueCutoff` of `gseGO`/`gseKEGG`. |
| `thresholds.pathway_padj_cutoff` | 0.1 | Adjusted p-value cutoff for pathways shown in dotplots. |
| `thresholds.top_n_pathways` | 7 | Up and down pathways shown in dotplots. |
| `enrichment.organism_db` | `"org.Hs.eg.db"` | OrgDb package (`org.Mm.eg.db` for mouse; both are in the environment). |
| `enrichment.ontology` | `["BP","MF","CC"]` | GO ontologies to test. The report shows BP. |
| `enrichment.kegg_organism` | `"hsa"` | KEGG organism code (`mmu` for mouse). Needs internet access. |

> GSEA here ranks only genes that pass the p-value and fold-change filters. This is a focused,
> "pre-filtered" GSEA. To rank all genes, set `gsea_gene_pvalue_filter: 1` and `fc_soft: 0` (note
> that `fc_soft` also defines significant DEGs).

## 10. Plotting

| Key | Description |
|---|---|
| `plotting.col_by` | Column used to colour and facet the QC histograms (normally `Segment`). |
| `plotting.segment` | Optional fixed colours per segment, e.g. `{"PanCK+": "#E41A1C", "PanCK-": "#377EB8"}`. |
| `plot.volcano_xlim` / `volcano_ylim` | Axis limits of volcano plots (genes outside are not drawn). |

## 11. Deconvolution

Remove the three `deconvolution*` sections entirely to leave deconvolution out of the report.

### `deconvolution`

| Key | Default | Description |
|---|---|---|
| `run.batch_corrected` | `false` | Deconvolve RUV-corrected log-CPM (needs RUV factors from section 8). Adds per-scan barplots. |
| `run.no_batch_correction` | `true` | Deconvolve plain CPM. At least one of the two must be `true`. |
| `run.export_cibersort` | `true` | Write signature and mixture matrices for the CIBERSORT(x) web tool. |
| `signature.data_name` | `"safeTME"` | SpatialDecon built-in reference matrix. |
| `signature.cell_types_subset` | `[]` | Restrict to some cell types (`[]` = all 18). |
| `ruv_prefix` | `"ruv_"` | Prefix of the RUV factor columns. |
| `export_columns` | | `column_in_proportions.txt: annotation_column`. Include `segment: "Segment"`, which the statistics need. |
| `use_custom_palette` / `custom_palette` | `false` | Fixed colours per cell type (missing ones are grey). |

The proportions table uses the batch-corrected variant if it ran, and the plain CPM variant otherwise.

### `deconvolution_stats`

| Key | Default | Description |
|---|---|---|
| `grouping_columns` | | Columns of `proportions.txt` to test (must be `export_columns` keys). 2 levels: Wilcoxon; more: Kruskal-Wallis. |
| `alpha` | 0.05 | P-value threshold for plotting significant hits (also used for GSVA). The p-values are not adjusted for multiple testing. |
| `plots_per_page` | 6 | Boxplots per PDF page. |
| `segment_pair_comparison.enabled` | `false` | Compare two segments on summed cell-type groups. |
| `segment_pair_comparison.segments` | | Exactly two `Segment` values. |
| `segment_pair_comparison.cell_type_groups` | | `group_name: [cell types to sum]`. |
| `segment_pair_comparison.test_alternative` | | `group_name: two.sided / greater / less` (tested as first vs second segment). |

### `deconvolution_gsva`

| Key | Default | Description |
|---|---|---|
| `signature_file` | bundled Jerby-Arnon set | Excel file: one column per signature, gene symbols below the header. |
| `signature_start_row` | 2 | Row of the header (signature names). |
| `method` | `"gsva"` | `gsva`, `ssgsea`, `zscore` or `plage`. |
| `kcdf` | `"Poisson"` | `Poisson` for counts/CPM (used here), `Gaussian` for log-scale data (GSVA method only). |
| `min_sz` / `max_sz` | 3 / 500 | Gene-set size limits, after matching to the genes in the data. |
| `segment_test_enabled` | `true` | Pairwise Wilcoxon tests of scores between segments (BH-adjusted), drawn on a boxplot. |
| `heatmap_annotation_columns` | | Annotation columns shown above the heatmap. |
| `heatmap_scale` | `"row"` | `row`, `column` or `none`. |

## Adapting to a new project: checklist

- [ ] `outputdir`, `project_name`
- [ ] Raw-data paths and `annotation_sheet`
- [ ] `annotation.*` column names match the workbook
- [ ] QC thresholds reviewed on `QC_final.pdf`
- [ ] `spe.patient_source_column`, `scan_id_column`, `biology_columns` exist
- [ ] Batch correction only if a non-confounded batch exists
- [ ] `comparisons.column`, `segments` and `factor_of_interest` match your design
- [ ] `enrichment.*` match the organism
- [ ] `deconvolution.export_columns` and `deconvolution_stats.grouping_columns` use existing columns
