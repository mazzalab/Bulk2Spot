# Stage 4 -- limma-voom DEG + GO/KEGG GSEA for every segment x comparison.
# The set of comparisons depends on the data, so per-comparison files go into one directory() output.

rule deg_gsea:
    input:
        spe_ruv_rds=rules.batch_correction.output.spe_ruv_rds,
    output:
        summary_file=DEG_DIR + "/DEG_GSEA_summary.txt",
        deg_dir=directory(DEG_DIR + "/DEG"),
    log:
        LOGS + "/deg_gsea.log",
    resources:
        mem_mb=scaled(16000),
        runtime=scaled(90),
    conda:
        ANALYSIS_ENV
    script:
        "../scripts/04_deg_gsea.R"
