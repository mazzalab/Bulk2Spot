# Usage

- [Workflow targets](#workflow-targets)
- [The launcher: run.py](#the-launcher-runpy)
- [Local execution](#local-execution)
- [Cluster execution](#cluster-execution)
- [Re-running and resuming](#re-running-and-resuming)
- [Running without run.py](#running-without-runpy)
- [Logs and notifications](#logs-and-notifications)
- [Recommended project workflow](#recommended-project-workflow)

Always activate the launcher environment first (`conda activate bulk2spot`, see
[installation](installation.md)).

## Workflow targets

A target builds its own outputs and everything upstream of it.

| Target | Builds | Depends on |
|---|---|---|
| `preprocess` | QC, LOQ filtering, Q3 normalization (`results/preprocessing/`) | raw data |
| `spe` | SpatialExperiment (`results/spe/`) | `preprocess` |
| `batch` | diagnostics, TMM, optional RUV4/limma, PCA/UMAP (`results/batch_correction/`) | `spe` |
| `deg` | DEG + GSEA (`results/deg_gsea/`) | `batch` |
| `deconv` | deconvolution, statistics, GSVA (`results/deconvolution/`) | `batch`, `preprocess` |
| `report` | HTML report (`results/report/`) | all of the above |

```bash
./run.py --list-workflows
```

## The launcher: run.py

```
./run.py -w TARGET [TARGET ...] -c CONFIG [options]
```

| Option | Default | Description |
|---|---|---|
| `-w, --workflow` | (required) | One or more targets (table above). |
| `-c, --configfile` | (required) | Project YAML file. |
| `-q, --cores` | 1 | CPU cores for local execution. |
| `-n, --dry-run` | off | Show the jobs that would run, without running them. |
| `-f, --forceall` | off | Re-run every job of the target(s), even if outputs exist. |
| `-ri, --rerun-incomplete` | off | Re-run jobs whose outputs were left incomplete by a crash. |
| `-k, --keep-going` | off | Keep running independent jobs when one fails. |
| `-p, --printshellcmds` | off | Print shell commands (report rule). |
| `--dag` | off | Print the job graph in DOT format. |
| `--latency-wait` | 60 | Seconds to wait for output files on slow network filesystems. |
| `--rerun-triggers` | all | What makes a finished job re-run: `mtime params input software-env code`. |
| `--no-conda` | off | Run every rule in the currently active environment. |
| `--conda-prefix` | `.snakemake/conda` | Where environments are built and stored. |
| `--conda-frontend` | `conda` | `conda` or `mamba`. |
| `--create-envs-only` | off | Build the conda environments, then exit. |
| `-cl, --cluster` | off | Submit each job to PBS with `qsub`. |
| `-qu, --queue` | `workq` | PBS queue. |
| `--cluster-cmd` | – | Custom submit command for other schedulers (implies cluster mode). |
| `-j, --jobs` | 20 | Maximum number of jobs queued/running at once (cluster mode). |
| `-rt, --restart-times` | 2 | Retries per failed job. Each retry requests more memory and walltime. |

`run.py` allows only one run per config file at a time. Different projects (different config files)
can run at the same time from the same repository.

## Local execution

```bash
./run.py -w report -c my_project.yaml -n       # always dry-run first
./run.py -w report -c my_project.yaml -q 8
```

Jobs run on the current machine, using up to `-q` cores in total. Each rule's `mem_mb` is not
enforced locally, so make sure the machine has enough memory (about 16 GB for typical studies).

## Cluster execution

In cluster mode, `run.py` stays in the foreground and submits every job as a separate batch job. It
needs little CPU and memory itself, but it has to keep running until the workflow ends, so start it
inside `screen`/`tmux` or with `nohup`:

```bash
nohup ./run.py -w report -c my_project.yaml -cl -qu workq -j 10 > run.log 2>&1 &
```

### PBS Pro / OpenPBS (built in)

`-cl` submits each job with:

```
qsub -q <queue> -N b2s.<rule> -l select=1:ncpus=<threads>:mem=<mem_mb>mb -l walltime=<runtime>:00 \
     -o <outputdir>/logs/cluster/<rule>.<jobid>.out -e <outputdir>/logs/cluster/<rule>.<jobid>.err
```

### SLURM, LSF, SGE and other schedulers

Pass your own submit command with `--cluster-cmd`. Snakemake fills in the placeholders `{threads}`,
`{resources.mem_mb}` (MB), `{resources.runtime}` (minutes), `{rule}` and `{jobid}`. Create the log
directory first.

```bash
mkdir -p logs/slurm
# SLURM
./run.py -w report -c my_project.yaml -j 10 \
  --cluster-cmd "sbatch -p short -c {threads} --mem={resources.mem_mb} -t {resources.runtime} -J b2s.{rule} -o logs/slurm/{rule}.{jobid}.out"

# LSF
./run.py -w report -c my_project.yaml -j 10 \
  --cluster-cmd "bsub -n {threads} -M {resources.mem_mb} -W {resources.runtime} -J b2s.{rule} -o logs/lsf/{rule}.{jobid}.out"
```

### Resources and automatic retries

Each rule declares a base memory and walltime (`workflow/rules/*.smk`):

| Rule | Memory | Walltime |
|---|---|---|
| preprocessing | 16 GB | 60 min |
| spe_preparation | 16 GB | 30 min |
| batch_correction | 16 GB | 45 min |
| deg_gsea | 16 GB | 90 min |
| deconvolution | 16 GB | 45 min |
| deconvolution_stats | 8 GB | 30 min |
| deconvolution_gsva | 16 GB | 45 min |
| html_report | 4 GB | 30 min |

If a job fails (for example, killed for exceeding its memory or walltime), it is resubmitted up to
`-rt` times, with 2x and then 3x the base request. For very large studies, raise the base values in
the rule files.

## Re-running and resuming

Snakemake only runs jobs whose outputs are missing or out of date.

| Situation | Command |
|---|---|
| Changed a QC threshold in the config | Re-run the same command; jobs whose parameters changed re-run, along with everything downstream of them. |
| Run was interrupted (crash, `Ctrl-C`, walltime) | Re-run with `-ri`. |
| Re-run one stage and everything after it | Delete that stage's `results/<stage>/` folder and re-run, or use `-f` on the target. |
| See why a job would re-run | `-n` (the dry-run prints a `reason:` for every job). |
| Only timestamps matter (e.g. after copying a project) | `--rerun-triggers mtime` |

## Running without run.py

Bulk2Spot is a standard Snakemake 7 workflow, so you can also call Snakemake directly:

```bash
snakemake -s /path/to/Bulk2Spot/Snakefile --configfile my_project.yaml --use-conda -c 8 report
```

`run.py` adds the per-project lock, the PBS template, the cluster thread handling and sensible defaults
on top of this.

## Logs and notifications

| Log | Location |
|---|---|
| Per-rule analysis log (R output, messages, errors) | `<outputdir>/logs/<rule>.log` |
| Scheduler stdout/stderr (cluster mode) | `<outputdir>/logs/cluster/` |
| Snakemake run log | `.snakemake/log/` in the directory you launched from |

When a job fails, Snakemake prints the path of its log. Look at the end of that file first.

To get an e-mail when a run finishes or fails, set `notify_email` in the config. This needs a working
system `mail`/`mailx` command.

## Recommended project workflow

1. **Copy** `config/config.yaml` and fill in the paths and annotation columns.
2. **Preprocess** (`-w preprocess`) and review `results/preprocessing/QC_final.pdf`,
   `qc_summary_table.txt` and `QC_log.txt`. QC thresholds are dataset- and panel-specific: adjust
   `segment_qc` / `loq` / `gene_filtering` until the retained segments and genes make sense.
3. **Check batch effects** (`-w batch`): `PCAScree_prebatch.pdf` shows whether samples group by
   batch. Turn on `batch_correction.perform` only if a real batch variable, **not confounded with
   your biology**, exists.
4. **Set the design** (`comparisons`, `segments`) and run `-w report`.
5. **Commit your config** next to the results (or in a project repository) for reproducibility.
