#!/usr/bin/env bash
# Snakemake runs this once, inside the freshly created analysis environment.
# You can also run it by hand after `conda activate <your-env>`.
set -euo pipefail

# 1. Keep R from loading packages out of your personal library (~/R/...),
#    which would shadow this environment's versions.
mkdir -p "${CONDA_PREFIX}/etc/conda/activate.d"
echo "export R_LIBS_USER=\"${CONDA_PREFIX}/lib/R/library\"" \
  > "${CONDA_PREFIX}/etc/conda/activate.d/bulk2spot_r_libs.sh"
export R_LIBS_USER="${CONDA_PREFIX}/lib/R/library"

# 2. standR is not on conda-forge/bioconda: install it from Bioconductor
#    (the release matching this R version), falling back to GitHub.
Rscript -e '
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("standR", quietly = TRUE)) {
  tryCatch(
    BiocManager::install("standR", update = FALSE, ask = FALSE),
    error = function(e) remotes::install_github("DavisLaboratory/standR", upgrade = "never")
  )
}
stopifnot(requireNamespace("standR", quietly = TRUE))
message("standR ", packageVersion("standR"), " is installed.")
'
