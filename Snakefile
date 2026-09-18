# =============================================================================
# Bulk2Spot -- Snakemake workflow for NanoString GeoMx DSP spatial transcriptomics
#
#   preprocess -> spe -> batch -> deg
#                             \-> deconv -> report
#
# Launch with run.py (recommended), or with plain Snakemake 7:
#   snakemake --configfile config/config.yaml --use-conda -c 8 report
# =============================================================================

import os


# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

OUTDIR = config["outputdir"].rstrip("/")
RESULTS = OUTDIR + "/results"
LOGS = OUTDIR + "/logs"

PREP_DIR = RESULTS + "/preprocessing"
SPE_DIR = RESULTS + "/spe"
BATCH_DIR = RESULTS + "/batch_correction"
DEG_DIR = RESULTS + "/deg_gsea"
DECONV_DIR = RESULTS + "/deconvolution"
REPORT_DIR = RESULTS + "/report"


def resolve(path):
    """Relative input paths are looked up from the launch directory first, then
    from the Bulk2Spot repository (so bundled resources and the tutorial work
    no matter where run.py is launched from)."""
    if os.path.isabs(path) or os.path.exists(path):
        return path
    return os.path.join(workflow.basedir, path)


# -----------------------------------------------------------------------------
# Software environments
# By default Snakemake builds them from workflow/envs/*.yaml (--use-conda).
# To reuse environments you already built, set their path/name in
# config `conda_envs.analysis` / `conda_envs.report`.
# -----------------------------------------------------------------------------

_envs = config.get("conda_envs") or {}
ANALYSIS_ENV = _envs.get("analysis") or os.path.join(workflow.basedir, "workflow/envs/bulk2spot.yaml")
REPORT_ENV = _envs.get("report") or os.path.join(workflow.basedir, "workflow/envs/report.yaml")

# Fix for Snakemake 7.32: it passes pre-built environment *paths* to conda as
# `--name`, which conda >= 23 rejects. Use `--prefix` for anything that is a path.
from snakemake.deployment.conda import Env
Env.address_argument = property(
    lambda env: f"--prefix '{env.address}'" if os.sep in env.address else f"--name '{env.address}'"
)


# -----------------------------------------------------------------------------
# Retry-scaled resources
# Each rule requests base * attempt, so a job killed for exceeding memory or
# walltime is resubmitted with 2x, 3x ... (see run.py --restart-times).
# -----------------------------------------------------------------------------

def scaled(base):
    return lambda wildcards, attempt: base * attempt


# -----------------------------------------------------------------------------
# Top-level targets (what `run.py -w <target>` builds)
# -----------------------------------------------------------------------------

rule all:
    input:
        REPORT_DIR + "/bulk2spot_report.html",

rule preprocess:
    input:
        PREP_DIR + "/norm_target_data.rds",
        PREP_DIR + "/QC_final.pdf",

rule spe:
    input:
        SPE_DIR + "/spe.rds",

rule batch:
    input:
        BATCH_DIR + "/spe_ruv.rds",
        BATCH_DIR + "/RUV.pdf",

rule deg:
    input:
        DEG_DIR + "/DEG_GSEA_summary.txt",
        DEG_DIR + "/DEG",

rule deconv:
    input:
        DECONV_DIR + "/proportions.txt",
        DECONV_DIR + "/stats",
        DECONV_DIR + "/GSVA_Heatmap.pdf",
        DECONV_DIR + "/GSVA_stats",

rule report:
    input:
        REPORT_DIR + "/bulk2spot_report.html",


# -----------------------------------------------------------------------------
# Modules (order matters: later modules reference earlier rules' outputs)
# -----------------------------------------------------------------------------

include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/spe_preparation.smk"
include: "workflow/rules/batch_correction.smk"
include: "workflow/rules/deg_gsea.smk"
include: "workflow/rules/deconvolution.smk"
include: "workflow/rules/report.smk"


# -----------------------------------------------------------------------------
# Optional e-mail notification (config `notify_email`, uses the system `mail`)
# -----------------------------------------------------------------------------

def notify(subject, snakemake_log):
    email = config.get("notify_email")
    if not email:
        return
    import shutil
    import subprocess
    mail = shutil.which("mail") or shutil.which("mailx")
    if mail:
        body = f"Output directory: {OUTDIR}\nSnakemake log: {snakemake_log}\n"
        subprocess.run([mail, "-s", subject, email], input=body.encode(), check=False, timeout=15)


# `log` (the Snakemake run log path) is provided by Snakemake inside these handlers.
onsuccess:
    notify(f"[Bulk2Spot] SUCCESS: {config.get('project_name', OUTDIR)}", log)

onerror:
    notify(f"[Bulk2Spot] FAILED: {config.get('project_name', OUTDIR)}", log)
