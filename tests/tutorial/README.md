# Bulk2Spot tutorial / smoke-test dataset

A small, real dataset for testing Bulk2Spot end to end -- from raw DCC/PKC/annotation files through
QC/normalization, `SpatialExperiment` conversion, batch-correction diagnostics, limma-voom DEG +
GO/KEGG GSEA, deconvolution and the HTML report -- without waiting on a full production run.

## Dataset

NanoString's public **WTA kidney demo dataset** -- the exact dataset used in Bioconductor's
["Analyzing GeoMx-NGS RNA Expression Data with
GeomxTools"](https://bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html)
vignette: 239 segments across 8 tissue slides (Diabetic Kidney Disease [DKD] vs normal), 3 segment
types (`Geometric Segment`, `PanCK+`, `PanCK-`), Whole Transcriptome Atlas (WTA) panel.

The files are fetched from
[Nanostring-Biostats/GeoMxWorkflows](https://github.com/Nanostring-Biostats/GeoMxWorkflows)
(`inst/extdata/WTA_NGS_Example/`, the raw files bundled inside that Bioconductor package) via
`raw.githubusercontent.com` -- no need to install the R package just to unpack them. ~75 MB in total
(239 DCC files, one ~1.8 MB PKC zip, one annotation workbook).

## Usage

All commands run from the **repository root** (the config's relative paths assume it).

```bash
# 1. Download the tutorial data (one-time, well under a minute on a good connection)
tests/tutorial/download_tutorial_data.sh          # default: 8 parallel downloads

# 2. Activate the launcher environment (see docs/installation.md)
conda activate bulk2spot

# 3. Run the whole pipeline, including the HTML report
./run.py -w report -c tests/tutorial/config_tutorial.yaml -q 8
# or on a PBS cluster, one job per rule:
./run.py -w report -c tests/tutorial/config_tutorial.yaml -cl -qu workq -j 8
```

The first run also builds the two conda environments (analysis + report), which takes a while;
later runs reuse them.

`-w report` builds everything upstream of it (preprocessing, SpatialExperiment conversion, batch
diagnostics, DEG/GSEA, deconvolution), so this one command exercises the whole pipeline. Key outputs:

```
tests/tutorial/output/results/preprocessing/norm_target_data.rds
tests/tutorial/output/results/spe/spe.rds
tests/tutorial/output/results/batch_correction/spe_ruv.rds
tests/tutorial/output/results/deg_gsea/DEG_GSEA_summary.txt
tests/tutorial/output/results/deg_gsea/DEG/<segment>/<comparison>/
tests/tutorial/output/results/deconvolution/proportions.txt
tests/tutorial/output/results/report/bulk2spot_report.html
```

The report (`bulk2spot_report.html`) is the easiest way to look through everything at once: QC funnel,
sample composition, PCA/UMAP, interactive volcano/GSEA plots and tables for every segment x comparison,
and the deconvolution results -- self-contained, open it directly in a browser, no server needed.

## This config's dataset-specific choices (vs `config/config.yaml`)

This dataset's annotation columns (`Sample_ID`, `panel`, `slide name`, `class`, `roi`, `segment`,
`aoi`, `area`, `region`, `pathology`, `nuclei`) don't match the generic placeholders in
`config/config.yaml` -- see `config_tutorial.yaml`'s inline comments for exactly what was changed and
why, in particular:

- `comparisons.column: "class"` (DKD vs normal) -- this dataset has no "Group" column; the pipeline
  doesn't hardcode one.
- `segments: ["Geometric Segment", "PanCK+", "PanCK-"]` -- the segment types present.
- `batch_correction.perform: false` -- this dataset has no genuinely independent batch variable (each
  tissue slide is essentially 1:1 with disease class, so using it as "batch" would confound batch
  with biology). Left off rather than fabricated -- this also exercises the pipeline's
  "batch correction skipped" code path.
- `segment_qc.remove_flagged_segments: true` -- this dataset includes 4 "No Template Control"
  segments (a sequencing negative control, no tissue) that are expected to fail segment QC and carry
  `NA` for class/segment/region; they are removed rather than carried into later stages.

## Expected results

**Segment QC**: 212/235 segments pass (90%; 23 flagged, mostly `LowNegatives`) -- the 4 "No Template
Control" segments are removed by `remove_flagged_segments: true`. Final normalized dataset:
**9,922 genes x 212 segments**.

These numbers require recalibrating `segment_qc`'s `minNegativeCount`/`maxNTCCount` away from
`config/config.yaml`'s defaults (10 / 1000) -- checked empirically first (negative-control means per
segment: 1-50, median 6.4; NTC counts: up to 8704), since the defaults would fail >85% of segments
here. This is a normal part of GeoMx analysis, not a pipeline bug: QC thresholds are
dataset/panel-specific.

**Batch correction**: skipped (`batch_correction.perform: false`) -- `RUV.pdf` is a one-page
placeholder (a few KB, vs. hundreds of KB for the real diagnostic plots).

**DEG + GSEA** (limma-voom, DKD vs normal, all 3 segment types):

| Segment | N tested | Up | Down | GO BP terms | KEGG terms |
|---|---|---|---|---|---|
| Geometric Segment | 9,921 | 279 | 264 | 120 | 0 |
| PanCK+ | 9,921 | 398 | 393 | 229 | 31 |
| PanCK- | 9,921 | 86 | 136 | 69 | 1 |

The DEG counts are deterministic. GSEA term counts are not fully reproducible across installations:
`gseGO()`/`gseKEGG()` estimate p-values by permutation (the pipeline fixes the seed via
`dimensionality_reduction.seed`), and KEGG data is downloaded at run time, so results change with
KEGG releases. For example, an earlier run with the same packages found 128/246/69 GO BP terms and
6/32/1 KEGG pathways.

The top PanCK+ DEGs by P-value are real human gene symbols (`MYH10`, `PAPLN`, `ATP6V1C2`, `PLS3`,
`MAN1A1`), confirming the probe annotation -> gene-level counts -> limma-voom -> Entrez ID mapping for
GSEA works end to end.

All 3 segment types produce the full output set (`TableDEG_*.xlsx`, `Volcano_*.pdf`,
`Dotplot_BP_*.pdf`/`DotplotKEGG_*.pdf` and their result tables) under
`output/results/deg_gsea/DEG/<segment>/DKD_vs_normal/`.

**Deconvolution**: **212 ROIs deconvolved against all 18 safeTME cell types**; **24 of the 31
signature columns** of the bundled Jerby-Arnon workbook pass `min_sz` against this WTA gene panel;
**22 significant cell-type comparisons** and **45 significant signature comparisons** (p < 0.05,
grouping column `class`, DKD vs normal), shown in the report's significance tables.

The negative-control probe name in `norm_target_data.rds` is **panel-dependent**: this WTA panel
already uses standR's name (`"NegProbe-WTX"`), while e.g. CTA panels still carry GeomxTools' raw
`"Negative Probe"`. `05_deconvolution.R` looks for both (`spe.negative_probe_target_name` first).

**HTML report**: `bulk2spot_report.html` (~7 MB, self-contained) with 26 interactive Plotly figures
(QC funnels, QC flags, sample composition, PCA/UMAP, volcano and GSEA plots per segment,
deconvolution composition and GSVA heatmap), and no "figure not available" placeholders.

**Run time**: about 20 minutes wall time on a PBS cluster (one job per rule, including queueing),
once the environments exist.

A copy of this run's report, figures and summary tables is in
[`examples/tutorial_output/`](../../examples/README.md), so you can compare your results with it.

## Cleaning up

The downloaded raw data (`data/`) and the pipeline output (`output/`) are git-ignored -- delete them
freely with `rm -rf tests/tutorial/data tests/tutorial/output` and re-run
`download_tutorial_data.sh` any time.
