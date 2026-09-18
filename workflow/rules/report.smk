# Stage 6 -- self-contained interactive HTML report.
# Runs in its own light Python environment: it only re-renders results computed upstream.

# The deconvolution section is included only when the config has a `deconvolution:` block.
WITH_DECONV = "deconvolution" in config
# The first comparison column is the "primary" grouping used for colouring plots.
_group = config["comparisons"]["column"]
PRIMARY_GROUP = _group if isinstance(_group, str) else _group[0]


rule html_report:
    input:
        qc_log=rules.preprocessing.output.qc_log,
        pheno=rules.spe_preparation.output.phenoData_txt,
        deg_summary=rules.deg_gsea.output.summary_file,
        deg_dir=rules.deg_gsea.output.deg_dir,
        pca_final_txt=rules.batch_correction.output.pca_final_txt,
        umap_final_txt=rules.batch_correction.output.umap_final_txt,
        # `[]` means "no input" when deconvolution is not configured
        deconv_proportions=rules.deconvolution.output.proportions_file if WITH_DECONV else [],
        deconv_gsva_long=rules.deconvolution_gsva.output.gsva_long_txt if WITH_DECONV else [],
        deconv_stats_dir=rules.deconvolution_stats.output.stats_dir if WITH_DECONV else [],
        deconv_gsva_stats_dir=rules.deconvolution_gsva.output.gsva_stats_dir if WITH_DECONV else [],
        # Anchored to the repo, not the launch directory
        script=os.path.join(workflow.basedir, "workflow/scripts/generate_report.py"),
        logo=os.path.join(workflow.basedir, "docs/images/logo-mark.svg"),
    output:
        html=REPORT_DIR + "/bulk2spot_report.html",
    params:
        project_name=config.get("project_name", os.path.basename(OUTDIR)),
        group_column=PRIMARY_GROUP,
        p_cutoff=config["thresholds"]["p_cutoff"],
        fc_soft=config["thresholds"]["fc_soft"],
        batch_flag="--batch-performed" if config["batch_correction"].get("perform") else "",
        deconv_flag=f"--deconv-dir {DECONV_DIR}" if WITH_DECONV else "",
        # Metadata (non cell-type) columns of proportions.txt
        deconv_id_cols=",".join(["ROI"] + list(config.get("deconvolution", {}).get("export_columns", {}))),
    log:
        LOGS + "/html_report.log",
    resources:
        mem_mb=scaled(4000),
        runtime=scaled(30),
    conda:
        REPORT_ENV
    shell:
        """
        python {input.script} \
            --project-name "{params.project_name}" \
            --preprocessing-dir {PREP_DIR} \
            --spe-dir {SPE_DIR} \
            --batch-dir {BATCH_DIR} \
            --deg-dir {DEG_DIR} \
            --group-column "{params.group_column}" \
            --p-cutoff {params.p_cutoff} \
            --fc-soft {params.fc_soft} \
            {params.batch_flag} \
            {params.deconv_flag} \
            --deconv-id-cols "{params.deconv_id_cols}" \
            --logo {input.logo} \
            --output {output.html} > {log} 2>&1
        """
