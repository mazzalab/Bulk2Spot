# Stage 1 -- load DCC/PKC/annotation, segment + probe QC, LOQ filtering, Q3 normalization.

rule preprocessing:
    input:
        dcc_dir=resolve(config["dcc_dir"]),
        pkc=resolve(config["pkc_file"]),
        annotation=resolve(config["annotation_file"]),
    output:
        target_rds=PREP_DIR + "/target_data.rds",
        norm_target_rds=PREP_DIR + "/norm_target_data.rds",
        rdata_image=PREP_DIR + "/Preprocessing.RData",
        qc_final_pdf=PREP_DIR + "/QC_final.pdf",
        qc_summary_table=PREP_DIR + "/qc_summary_table.txt",
        ntc_summary=PREP_DIR + "/ntc_count_summary.txt",
        qc_log=PREP_DIR + "/QC_log.txt",
        data_raw=PREP_DIR + "/data_raw.txt",
        data_q3=PREP_DIR + "/data_Q3.txt",
        annotation_reordered=PREP_DIR + "/annotation_reordered.xlsx",
    log:
        LOGS + "/preprocessing.log",
    resources:
        mem_mb=scaled(16000),
        runtime=scaled(60),
    conda:
        ANALYSIS_ENV
    script:
        "../scripts/01_preprocessing.R"
