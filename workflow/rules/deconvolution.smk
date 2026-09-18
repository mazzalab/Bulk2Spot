# Stage 5 -- SpatialDecon cell-type deconvolution, proportion statistics and custom-signature GSVA.

# Gene signatures scored by GSVA (default: the bundled Jerby-Arnon set)
GSVA_SIGNATURES = config.get("deconvolution_gsva", {}).get(
    "signature_file", "resources/signatures/Jerby_Arnon_signature.xlsx")

rule deconvolution:
    input:
        norm_target_rds=rules.preprocessing.output.norm_target_rds,
        spe_ruv_rds=rules.batch_correction.output.spe_ruv_rds,
    output:
        proportions_file=DECONV_DIR + "/proportions.txt",
        fake_obj_rds=DECONV_DIR + "/fake_obj.rds",
        norm_counts_matrix_rds=DECONV_DIR + "/norm_counts_matrix.rds",
        normCPM_rds=DECONV_DIR + "/normCPM.rds",
        # Written as real results or placeholders depending on `deconvolution.run` toggles
        res_dec_batch_corrected_rds=DECONV_DIR + "/res_dec_batch_corrected.rds",
        res_dec_no_batch_correction_rds=DECONV_DIR + "/res_dec_no_batch_correction.rds",
        barplot_ordered=DECONV_DIR + "/Barplot_ordered.pdf",
        barplot_by_scan=DECONV_DIR + "/Barplot_by_Scan.pdf",
        barplot_ordered_nb=DECONV_DIR + "/Barplot_ordered_no_batch.pdf",
        cibersort_signature_file=DECONV_DIR + "/cibersort_signature.txt",
        cibersort_mixture_file=DECONV_DIR + "/cibersort_mixture.txt",
    log:
        LOGS + "/deconvolution.log",
    resources:
        mem_mb=scaled(16000),
        runtime=scaled(45),
    conda:
        ANALYSIS_ENV
    script:
        "../scripts/05_deconvolution.R"


rule deconvolution_stats:
    input:
        proportions_file=rules.deconvolution.output.proportions_file,
    output:
        # One set of plots/tables per grouping column -> directory() output
        stats_dir=directory(DECONV_DIR + "/stats"),
        segment_pair_plot=DECONV_DIR + "/SegmentPair_comparison.pdf",
    log:
        LOGS + "/deconvolution_stats.log",
    resources:
        mem_mb=scaled(8000),
        runtime=scaled(30),
    conda:
        ANALYSIS_ENV
    script:
        "../scripts/06_deconvolution_stats.R"


rule deconvolution_gsva:
    input:
        normCPM_rds=rules.deconvolution.output.normCPM_rds,
        fake_obj_rds=rules.deconvolution.output.fake_obj_rds,
        signature_file=resolve(GSVA_SIGNATURES),
    output:
        gsva_scores_rds=DECONV_DIR + "/gsva_scores.rds",
        gsva_long_rds=DECONV_DIR + "/gsva_long.rds",
        gsva_long_txt=DECONV_DIR + "/gsva_long.txt",
        segment_boxplot=DECONV_DIR + "/GSVA_Segment_boxplot.pdf",
        heatmap_file=DECONV_DIR + "/GSVA_Heatmap.pdf",
        gsva_stats_dir=directory(DECONV_DIR + "/GSVA_stats"),
    log:
        LOGS + "/deconvolution_gsva.log",
    resources:
        mem_mb=scaled(16000),
        runtime=scaled(45),
    conda:
        ANALYSIS_ENV
    script:
        "../scripts/07_deconvolution_gsva.R"
