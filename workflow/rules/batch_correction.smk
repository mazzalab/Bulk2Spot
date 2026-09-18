# Stage 3 -- RLE/PCA diagnostics, TMM normalization, optional RUV4 (+ limma) batch correction.
# Plots that do not apply (e.g. correction switched off) are written as one-page placeholders,
# so every declared output always exists.

rule batch_correction:
    input:
        spe_rds=rules.spe_preparation.output.spe_rds,
    output:
        spe_ruv_rds=BATCH_DIR + "/spe_ruv.rds",
        sample_info_plot=BATCH_DIR + "/SampleInfo_Plot.pdf",
        pca_scree_prebatch=BATCH_DIR + "/PCAScree_prebatch.pdf",
        pca_scree_postbatch=BATCH_DIR + "/PCAScree_postbatch.pdf",
        ruv_plot=BATCH_DIR + "/RUV.pdf",
        # Final PCA/UMAP coordinates, read by the HTML report
        pca_final_txt=BATCH_DIR + "/PCA_final.txt",
        pca_variance_txt=BATCH_DIR + "/PCA_final_variance.txt",
        umap_final_txt=BATCH_DIR + "/UMAP_final.txt",
    log:
        LOGS + "/batch_correction.log",
    resources:
        mem_mb=scaled(16000),
        runtime=scaled(45),
    conda:
        ANALYSIS_ENV
    script:
        "../scripts/03_batch_correction.R"
