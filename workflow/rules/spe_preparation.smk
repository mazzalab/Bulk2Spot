# Stage 2 -- convert the normalized GeoMxSet into a standR SpatialExperiment.

rule spe_preparation:
    input:
        norm_target_rds=rules.preprocessing.output.norm_target_rds,
    output:
        spe_rds=SPE_DIR + "/spe.rds",
        exprs_txt=SPE_DIR + "/exprs.txt",
        phenoData_txt=SPE_DIR + "/phenoData.txt",
        featureData_txt=SPE_DIR + "/featureData.txt",
    log:
        LOGS + "/spe_preparation.log",
    resources:
        mem_mb=scaled(16000),
        runtime=scaled(30),
    conda:
        ANALYSIS_ENV
    script:
        "../scripts/02_spe_preparation.R"
