# Methods

This page describes what each stage computes, so that you can interpret the results and write a
methods section. Config keys are given in `code`.

## Background

The NanoString GeoMx Digital Spatial Profiler (DSP) [1] measures RNA expression in user-selected
regions of a tissue section. Regions of interest (ROIs) are drawn on the slide, and can be split into
areas of illumination (AOIs, "segments") by marker staining, e.g. PanCK+ epithelium vs PanCK− stroma.
Oligonucleotide-tagged probes are released from each segment by UV light and counted by sequencing.
Every segment therefore yields a bulk-like expression profile of a few hundred to thousands of cells:
the "bulk" to "spot" idea behind Bulk2Spot.

## Stage 1: Preprocessing (`01_preprocessing.R`)

1. **Loading.** The DCC (counts), PKC (probe design) and annotation files are read with
   `GeomxTools::readNanoStringGeoMxSet()` [2]. Zero counts are shifted to one (`shiftCountsOne`).
2. **Segment QC.** `setSegmentQCFlags()` flags segments with low reads, poor trimming, stitching or
   alignment, low sequencing saturation, low negative-probe signal, high no-template-control counts,
   few nuclei or a small area (`segment_qc`). Flagged segments are removed if
   `remove_flagged_segments` is set.
3. **Probe QC.** `setBioProbeQCFlags()` removes probes with a low ratio to their target's other
   probes or that are global Grubbs outliers (`probe_qc`). Probe counts are then collapsed to gene
   level (`aggregateCounts`, geometric mean of probes).
4. **Limit of quantification.** For each segment and panel module,
   LOQ = max(`minLOQ`, NegGeoMean × NegGeoSD^`cutoff`), from the negative-control probes. A gene is
   *detected* in a segment when its count exceeds the LOQ.
5. **Filtering.** Segments with too few detected genes (`segment_detection_rate_min`) and genes
   detected in too few segments (`gene_detection_rate_min`) are removed. Negative probes are kept for
   later background modelling.
6. **Q3 normalization.** Counts are scaled to the 75th percentile of each segment
   (`normalization`). This is used for QC exports; downstream analyses use TMM (stage 3).

## Stage 2: SpatialExperiment (`02_spe_preparation.R`)

The filtered data are converted to a `SpatialExperiment` [3] with standR [4] (`readGeoMx`,
`readGeoMxFromDGE`). The negative probe is renamed to standR's convention. Anonymized `Patient` and
`Scan` labels and a `biology` variable (pasted `biology_columns`) are added.

## Stage 3: Diagnostics and batch correction (`03_batch_correction.R`)

1. **Diagnostics before correction.** Relative log expression (RLE) plots and a PCA (scater [5])
   coloured by group, scan, segment and batch.
2. **Normalization.** standR `geomxNorm()` with TMM [6] (`normalization_standR.method`).
3. **Batch correction** (optional; requires `perform: true` and an existing `batch_column`).
   - Negative control genes (NCGs): the `ncg_top_n` genes least variable across batches
     (`findNCGs`).
   - The number of unwanted factors *k* is chosen by silhouette score (`findBestK`) or fixed (`k`).
   - RUV4 [7] (`geomxBatchCorrection`) removes unwanted variation while preserving
     `factor_of_interest`. It adds factors `ruv_W1..k`.
   - Optionally, `limma::removeBatchEffect()` is applied to the corrected log-counts
     (`apply_limma_residual_correction`).
4. **Final embedding.** PCA and UMAP [8] of the final log-counts. Their coordinates are exported for
   the report, and post-correction diagnostics are plotted.

Batch correction is valid only when batches are not confounded with the biology of interest.
RUV4 assumes that the factor of interest varies within batches.

## Stage 4: Differential expression and enrichment (`04_deg_gsea.R`)

For each comparison column, each segment in `segments`, and each pair of groups:

1. **Model.** A no-intercept design `~ 0 + group` is fitted on the segments of that type (optionally
   `+ ruv_W1..3`, `design.include_ruv_factors`), using the TMM-normalized counts. The data are
   transformed with voom [9], fitted with limma [10] (`lmFit`) and tested for the contrast
   *test − reference* with robust empirical Bayes moderation (`eBayes(robust = TRUE)`). P-values are
   adjusted with Benjamini–Hochberg.
2. **Significant genes.** adj.P.Val < `p_cutoff` and |log2FC| > `fc_soft`.
3. **GSEA.** Genes with P.Value < `gsea_gene_pvalue_filter` and |log2FC| > `fc_soft` are ranked by
   log2FC (one probe per Entrez ID, the one with the largest |log2FC|). This ranking is tested with
   clusterProfiler [11] `gseGO` for each GO ontology and `gseKEGG`, based on the GSEA method [12]
   (BH adjustment, `gsea_pvalue_cutoff`, fixed random seed). Core-enrichment genes are reported as
   gene symbols.

Comparisons with fewer than `min_samples_per_group` segments in either group are skipped.

## Stage 5: Deconvolution (`05`–`07`)

### Cell-type deconvolution (`05_deconvolution.R`)

SpatialDecon [13] models each segment's expression as a non-negative combination of reference
cell-type profiles, plus a background term derived from the negative probes
(`derive_GeoMx_background`). The default reference is **safeTME** (18 immune and stromal cell types,
designed for tumour microenvironments and robust to cancer-cell expression). Two inputs can be used:

- **no_batch_correction**: counts per million (CPM) of the raw counts;
- **batch_corrected**: log-CPM with the RUV factors regressed out (`limma::removeBatchEffect`).

The reported proportions are `prop_of_nontumor`: the fractions of the estimated immune/stromal
content in each segment. They are relative abundances, not cell counts. Optionally, a signature and a
mixture matrix are exported for CIBERSORT(x) [14].

### Proportion statistics (`06_deconvolution_stats.R`)

For each grouping column, each segment type and each cell type: Wilcoxon rank-sum test (2 groups) or
Kruskal–Wallis test (more than 2). Hits with p < `alpha` are plotted. These p-values are
**not adjusted for multiple testing** and are meant for exploration. An optional test compares summed
cell-type groups (e.g. all T cells) between two segment types (Wilcoxon).

### Signature scoring (`07_deconvolution_gsva.R`)

GSVA [15] computes one enrichment score per gene signature per segment from CPM values
(`kcdf: Poisson`). By default it uses the tumour-microenvironment programs of Jerby-Arnon *et al.*
[16] (`resources/signatures/`); any Excel file with one gene list per column can be used. Scores are
compared between segment types (pairwise Wilcoxon, BH-adjusted) and between groups within each
segment type (as for the proportions, unadjusted), and shown as a heatmap.

## Stage 6: Report (`generate_report.py`)

A Python script collects the tables written by the previous stages into one HTML page with
interactive Plotly figures. The only computation it performs is a quick overview PCA of the
log2(Q3 + 1) data, shown before correction.

## Reproducibility

- Software versions are pinned (`workflow/envs/*.pin.txt`); tested versions are listed in
  [installation](installation.md#4-verify-the-installation).
- Random seeds (PCA, UMAP, GSEA) come from `dimensionality_reduction.seed`.
- KEGG pathways are downloaded at run time, so KEGG results can change with KEGG releases.
- Snakemake re-runs jobs when their code, parameters or inputs change.

## Example methods paragraph

> GeoMx DSP data were analysed with Bulk2Spot (version X), a Snakemake [17] workflow. Segments and
> probes were quality-controlled with GeomxTools (v3.5.0), using thresholds of …; genes detected
> above the limit of quantification (2 geometric SD above the negative-probe geometric mean) in at
> least 10% of segments were retained. Data were converted with standR (v1.6.0) and TMM-normalized.
> Differential expression between … was assessed per segment type with limma-voom (v3.58.1), with
> genes at BH-adjusted p < 0.1 and |log2FC| > 0.57 considered significant. Gene set enrichment was
> performed with clusterProfiler (v4.10.0) on GO Biological Process and KEGG. Cell-type proportions
> were estimated with SpatialDecon (v1.12.0) using the safeTME matrix, and signature scores were
> computed with GSVA (v1.50.0).

## References

1. Merritt CR, *et al.* Multiplex digital spatial profiling of proteins and RNA in fixed tissue. *Nat Biotechnol* 38, 586–599 (2020).
2. Ortogero N, *et al.* GeomxTools: NanoString GeoMx Tools. Bioconductor R package. doi:10.18129/B9.bioc.GeomxTools
3. Righelli D, *et al.* SpatialExperiment: infrastructure for spatially-resolved transcriptomics data in R using Bioconductor. *Bioinformatics* 38, 3128–3131 (2022).
4. Liu N, *et al.* standR: spatial transcriptomic analysis for GeoMx DSP data. *Nucleic Acids Res* 52, e2 (2024).
5. McCarthy DJ, *et al.* Scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R. *Bioinformatics* 33, 1179–1186 (2017).
6. Robinson MD, Oshlack A. A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biol* 11, R25 (2010).
7. Gagnon-Bartsch JA, Speed TP. Using control genes to correct for unwanted variation in microarray data. *Biostatistics* 13, 539–552 (2012).
8. McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for dimension reduction. arXiv:1802.03426 (2018).
9. Law CW, *et al.* voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biol* 15, R29 (2014).
10. Ritchie ME, *et al.* limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res* 43, e47 (2015).
11. Wu T, *et al.* clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. *The Innovation* 2, 100141 (2021).
12. Subramanian A, *et al.* Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 102, 15545–15550 (2005).
13. Danaher P, *et al.* Advances in mixed cell deconvolution enable quantification of cell types in spatial transcriptomic data. *Nat Commun* 13, 385 (2022).
14. Newman AM, *et al.* Robust enumeration of cell subsets from tissue expression profiles. *Nat Methods* 12, 453–457 (2015).
15. Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. *BMC Bioinformatics* 14, 7 (2013).
16. Jerby-Arnon L, *et al.* A cancer cell program promotes T cell exclusion and resistance to checkpoint blockade. *Cell* 175, 984–997 (2018).
17. Mölder F, *et al.* Sustainable data analysis with Snakemake. *F1000Research* 10, 33 (2021).
