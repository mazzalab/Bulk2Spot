# Installation

Bulk2Spot needs only **conda** and **git** on your machine. Everything else (Snakemake, R,
Bioconductor, standR, Python report libraries) is installed into conda environments.

- [Requirements](#requirements)
- [1. Get the code](#1-get-the-code)
- [2. Launcher environment](#2-launcher-environment)
- [3. Analysis environments](#3-analysis-environments)
- [4. Verify the installation](#4-verify-the-installation)
- [HPC notes](#hpc-notes)
- [Troubleshooting](#troubleshooting)

## Requirements

| Requirement | Notes |
|---|---|
| Linux x86-64 | Tested on a RHEL 8 HPC cluster with OpenPBS. The pin files are Linux-only; on other platforms Snakemake solves the `.yaml` files instead (untested). |
| [conda](https://docs.conda.io) >= 23.10 | Miniforge, Miniconda or Anaconda. Recent conda uses the fast libmamba solver; `mamba` is optional. |
| git | To clone the repository. |
| Internet access | To build environments (once), and at run time for KEGG enrichment. |
| Disk | ~5 GB for the environments, plus your data and results. |
| Memory | 16 GB per analysis job is enough for typical GeoMx studies (a few hundred segments). |

Configure the conda channels once if you have not already done so:

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

## 1. Get the code

```bash
git clone https://github.com/<OWNER>/Bulk2Spot.git
cd Bulk2Spot
```

## 2. Launcher environment

The launcher environment contains Snakemake 7.32.4 and what `run.py` needs. You activate it
whenever you run the pipeline.

```bash
conda env create -f environment.yml
conda activate bulk2spot
snakemake --version        # 7.32.4
```

To install it somewhere other than your default conda directory (e.g. a shared project space):

```bash
conda env create -p /shared/condaEnvs/bulk2spot -f environment.yml
conda activate /shared/condaEnvs/bulk2spot
```

> **Why Snakemake 7?** `run.py` uses the Snakemake 7 Python API and cluster interface. The pins on
> `pulp<2.8` and `setuptools<81` in `environment.yml` are needed for Snakemake 7.32 to start at all.

## 3. Analysis environments

The R stages and the report run in two separate environments, defined in `workflow/envs/`:

| File | Environment | Used by |
|---|---|---|
| `bulk2spot.yaml` | R 4.3, GeomxTools, standR, limma, clusterProfiler, SpatialDecon, GSVA, ... | stages 1-5 |
| `report.yaml` | Python 3.10, pandas, plotly | HTML report |

Choose **one** of the following options.

### Option A: automatic (recommended)

Do nothing. Snakemake builds the environments the first time a rule needs them (`run.py` always
passes `--use-conda`):

1. If a pin file (`<env>.linux-64.pin.txt`) is present, the exact package versions it lists are
   installed. These are the versions Bulk2Spot was tested with, and no dependency solving is needed.
   If that fails, Snakemake falls back to solving the `.yaml` file.
2. `bulk2spot.post-deploy.sh` then runs inside the new environment. It installs **standR** from
   Bioconductor (standR is not distributed through conda) and sets `R_LIBS_USER`, so that packages in
   your personal R library can't shadow the environment's versions.

Environments are stored in `.snakemake/conda/` inside the repository. To store them elsewhere, for
example to share them between several clones or users, use `--conda-prefix`:

```bash
./run.py -w report -c my_project.yaml -q 8 --conda-prefix /shared/snakemake_envs
```

To build the environments up front, without running any analysis, use `--create-envs-only` (it needs a
valid config, e.g. the tutorial one):

```bash
./run.py -w report -c tests/tutorial/config_tutorial.yaml --create-envs-only
```

On a shared HPC login node, run this inside a batch job ([HPC notes](#hpc-notes)).

### Option B: pre-built environments

This option suits a facility that maintains central environments, or machines without internet
access on the compute nodes.

```bash
# analysis environment
conda env create -p /shared/condaEnvs/bulk2spot-analysis -f workflow/envs/bulk2spot.yaml
conda activate /shared/condaEnvs/bulk2spot-analysis
bash workflow/envs/bulk2spot.post-deploy.sh           # installs standR
conda deactivate

# report environment
conda env create -p /shared/condaEnvs/bulk2spot-report -f workflow/envs/report.yaml
```

For exactly the tested versions, create them from the pin files instead:
`conda create -p <path> --file workflow/envs/bulk2spot.linux-64.pin.txt` (then run the post-deploy
script as above).

Then point your project config at them:

```yaml
conda_envs:
  analysis: "/shared/condaEnvs/bulk2spot-analysis"
  report: "/shared/condaEnvs/bulk2spot-report"
```

### Option C: no conda at run time

If R, all packages and Python are already available in your session (e.g. inside a container or a
module), run with `--no-conda`. Every rule then runs in the active environment.

## 4. Verify the installation

The tutorial doubles as an installation test:

```bash
conda activate bulk2spot
tests/tutorial/download_tutorial_data.sh
./run.py -w report -c tests/tutorial/config_tutorial.yaml -n        # dry-run: should list 9 jobs
./run.py -w report -c tests/tutorial/config_tutorial.yaml -q 8      # or -cl on a cluster
```

A successful run ends with `Workflow finished successfully` / `(100%) done` and creates
`tests/tutorial/output/results/report/bulk2spot_report.html`. Compare the numbers with the expected
results in [the tutorial README](../tests/tutorial/README.md#expected-results).

To check the analysis environment by hand:

```bash
conda activate .snakemake/conda/<hash>    # or your pre-built environment
Rscript -e 'for (p in c("GeomxTools","standR","SpatialDecon","GSVA","limma","clusterProfiler")) cat(p, format(packageVersion(p)), "\n")'
```

Tested versions: R 4.3.3, Bioconductor 3.18, GeomxTools 3.5.0, standR 1.6.0, SpatialDecon 1.12.0,
GSVA 1.50.0, limma 3.58.1, clusterProfiler 4.10.0.

## HPC notes

- **Don't run heavy work on login nodes.** Use `-cl` (one batch job per rule) for analyses. Build
  environments inside a batch job too, for example:

  ```bash
  cat > build_envs.pbs <<'EOF'
  #PBS -N bulk2spot_envs
  #PBS -l select=1:ncpus=4:mem=16gb
  #PBS -l walltime=04:00:00
  cd $PBS_O_WORKDIR
  source ~/miniconda3/etc/profile.d/conda.sh     # adapt to your conda installation
  conda activate bulk2spot
  ./run.py -w report -c tests/tutorial/config_tutorial.yaml --create-envs-only
  EOF
  qsub build_envs.pbs
  ```

- **Shared filesystem.** The repository, the conda prefix, the config file and `outputdir` must all
  be on a filesystem visible from the compute nodes (not `/tmp`).
- **No internet on compute nodes?** Build the environments first (above, or on a node with internet
  access). KEGG enrichment needs internet at run time: without it, KEGG is skipped with a message and
  GO enrichment still runs.
- **Other schedulers.** See [Usage → cluster execution](usage.md#cluster-execution).

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `ModuleNotFoundError: No module named 'pkg_resources'` or a `pulp` error when starting | The launcher environment was not created from `environment.yml`. Recreate it (it pins `setuptools<81` and `pulp<2.8`). |
| `package X was found, but >= Y is required` in an R log | A package in your personal R library (`~/R/...`) shadows the environment's version. Environments built by Bulk2Spot set `R_LIBS_USER` automatically; for your own environments, run `workflow/envs/bulk2spot.post-deploy.sh` inside them, or `export R_LIBS_USER=$CONDA_PREFIX/lib/R/library`. |
| `there is no package called 'standR'` | The post-deploy step did not run or failed (usually a network problem). Activate the environment and run `bash workflow/envs/bulk2spot.post-deploy.sh`. |
| Environment creation from the pin file fails | Snakemake falls back to the `.yaml` automatically. Pin files are Linux x86-64 only. A few pinned packages come from Anaconda's `defaults` channel: if your institution can't use it, delete the `*.pin.txt` files to solve from conda-forge/bioconda only. |
| "another run.py is already running with <config>" | Only one run per config file can run at a time, and one is still active (check with `ps -ef \| grep run.py`). The lock is released automatically when that process ends, even after a crash. |
| `IncompleteFilesException` after an interrupted run | Add `-ri` (`--rerun-incomplete`). |
| Environment build fails with `post-link script failed for package bioconda::bioconductor-go.db` (or `org.hs.eg.db`) and `curl: (18) ... bytes missing` | Bioconductor *data* packages download their data from bioconductor.org at install time, and the download was cut off. This happens on slow or proxied institutional networks. Try again (ideally at a quieter time), build on a machine with a better connection and point `--conda-prefix` at shared storage, or build the environment once by hand and use it through `conda_envs` (Option B). |
| A Bioconductor data package (e.g. `org.Hs.eg.db`, `GO.db`) fails to load at run time | Its data download failed silently while the environment was built. Remove the environment (`rm -rf .snakemake/conda/<hash>*`) and rebuild it. |
| `SafetyError: The package for ... appears to be corrupted` while creating an environment | A file in your conda package cache (`~/.../pkgs/`) no longer matches the package metadata, e.g. because it was edited by hand. Delete that package's folder and tarball from the cache (the path is printed in the message), so that a clean copy is downloaded. |
