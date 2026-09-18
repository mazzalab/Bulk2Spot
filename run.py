#!/usr/bin/env python3
"""Bulk2Spot launcher -- a thin wrapper around the Snakemake 7 Python API.

Examples
  ./run.py --list-workflows
  ./run.py -w report -c config/config.yaml -n                     # dry-run
  ./run.py -w report -c config/config.yaml -q 8                   # run locally on 8 cores
  ./run.py -w report -c config/config.yaml -cl -qu workq -j 10    # one PBS job per rule
  ./run.py -w report -c config/config.yaml --create-envs-only     # only build conda envs
"""

import argparse
import fcntl
import hashlib
import sys
from pathlib import Path

import snakemake
import yaml
from snakemake.resources import DefaultResources

REPO = Path(__file__).resolve().parent

TARGETS = {
    "preprocess": "Load DCC/PKC/annotation, segment + probe QC, LOQ filtering, Q3 normalization",
    "spe": "Convert to a standR SpatialExperiment (Patient / Scan / biology annotations)",
    "batch": "RLE/PCA diagnostics, TMM normalization, optional RUV4 (+ limma) batch correction",
    "deg": "limma-voom DEG + GO/KEGG GSEA for every segment x comparison",
    "deconv": "SpatialDecon deconvolution, proportion statistics, custom-signature GSVA",
    "report": "Self-contained interactive HTML report (runs everything above)",
}

# PBS Pro / OpenPBS submission template. {{...}} are Snakemake per-job placeholders;
# runtime is in minutes, so "{runtime}:00" is a valid MM:SS walltime.
PBS_TEMPLATE = (
    "qsub -q {queue} -N b2s.{{rule}} "
    "-l select=1:ncpus={{threads}}:mem={{resources.mem_mb}}mb "
    "-l walltime={{resources.runtime}}:00 "
    "-o {logdir}/{{rule}}.{{jobid}}.out -e {logdir}/{{rule}}.{{jobid}}.err"
)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-w", "--workflow", nargs="+", metavar="TARGET", help="target(s): " + ", ".join(TARGETS))
    p.add_argument("-c", "--configfile", help="project YAML config (see config/config.yaml)")
    p.add_argument("-q", "--cores", type=int, default=1, help="CPU cores for local execution (default: 1)")
    p.add_argument("--list-workflows", action="store_true", help="list targets and exit")

    run = p.add_argument_group("execution")
    run.add_argument("-n", "--dry-run", action="store_true", help="show what would run, run nothing")
    run.add_argument("-f", "--forceall", action="store_true", help="re-run every job of the target(s)")
    run.add_argument("-ri", "--rerun-incomplete", action="store_true", help="re-run jobs left incomplete by a crash")
    run.add_argument("-k", "--keep-going", action="store_true", help="keep running independent jobs after a failure")
    run.add_argument("-p", "--printshellcmds", action="store_true", help="print shell commands")
    run.add_argument("--dag", action="store_true", help="print the job DAG (pipe into `dot -Tsvg`)")
    run.add_argument("--latency-wait", type=int, default=60, help="seconds to wait for outputs on slow filesystems")
    run.add_argument("--rerun-triggers", nargs="+", choices=["mtime", "params", "input", "software-env", "code"],
                     help="what makes a finished job re-run (Snakemake default: all of them)")

    envs = p.add_argument_group("software environments")
    envs.add_argument("--no-conda", action="store_true", help="do not use conda; run in the active environment")
    envs.add_argument("--conda-prefix", help="where conda envs are built/stored (default: .snakemake/conda)")
    envs.add_argument("--conda-frontend", default="conda", choices=["conda", "mamba"], help="default: conda")
    envs.add_argument("--create-envs-only", action="store_true", help="build the conda envs and exit")

    cl = p.add_argument_group("cluster")
    cl.add_argument("-cl", "--cluster", action="store_true", help="submit every job to PBS with qsub")
    cl.add_argument("-qu", "--queue", default="workq", help="PBS queue (default: workq)")
    cl.add_argument("--cluster-cmd", help="custom submit command for other schedulers, e.g. "
                    "\"sbatch -p short -c {threads} --mem={resources.mem_mb} -t {resources.runtime}\"")
    cl.add_argument("-j", "--jobs", type=int, default=20, help="max jobs queued at once (default: 20)")
    cl.add_argument("-rt", "--restart-times", type=int, default=2,
                    help="retries per failed job, each with more memory/time (default: 2)")
    return p.parse_args()


def lock_config(configfile):
    """Refuse to start twice on the same config. This per-project lock replaces
    Snakemake's directory-wide lock, so different projects can run concurrently
    from the same repository. Released automatically when this process exits."""
    lock_dir = configfile.parent / ".bulk2spot_locks"
    lock_dir.mkdir(exist_ok=True)
    key = hashlib.sha1(str(configfile).encode()).hexdigest()[:16]
    handle = open(lock_dir / f"{key}.lock", "w")
    try:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        sys.exit(f"ERROR: another run.py is already running with {configfile}")
    return handle


def main():
    args = parse_args()

    if args.list_workflows:
        for name, description in TARGETS.items():
            print(f"  {name:<12} {description}")
        return 0
    if not args.configfile or not args.workflow:
        sys.exit("ERROR: -c/--configfile and -w/--workflow are required (see ./run.py -h)")
    unknown = [t for t in args.workflow if t not in TARGETS]
    if unknown:
        sys.exit(f"ERROR: unknown target(s) {unknown}; choose from {list(TARGETS)}")

    configfile = Path(args.configfile).resolve()
    if not configfile.exists():
        sys.exit(f"ERROR: config file not found: {configfile}")
    config = yaml.safe_load(configfile.open())
    _lock = lock_config(configfile)  # held until exit

    kwargs = dict(
        snakefile=str(REPO / "Snakefile"),
        configfiles=[str(configfile)],
        targets=args.workflow,
        cores=args.cores,
        dryrun=args.dry_run,
        forceall=args.forceall,
        force_incomplete=args.rerun_incomplete,
        keepgoing=args.keep_going,
        printshellcmds=args.printshellcmds,
        printdag=args.dag,
        latency_wait=args.latency_wait,
        restart_times=args.restart_times,
        use_conda=not args.no_conda,
        conda_prefix=args.conda_prefix,
        conda_frontend=args.conda_frontend,
        conda_create_envs_only=args.create_envs_only,
        lock=False,  # replaced by lock_config()
    )
    if args.rerun_triggers:
        kwargs["rerun_triggers"] = args.rerun_triggers

    if args.cluster or args.cluster_cmd:
        logdir = Path(config["outputdir"]).resolve() / "logs" / "cluster"
        logdir.mkdir(parents=True, exist_ok=True)
        kwargs["cluster"] = args.cluster_cmd or PBS_TEMPLATE.format(queue=args.queue, logdir=logdir)
        kwargs["nodes"] = args.jobs
        # Fallback resources for jobs that do not declare their own
        kwargs["default_resources"] = DefaultResources(["mem_mb=8000", "runtime=120"])
        # Snakemake caps each job's {threads} at `cores` even on a cluster:
        # use a large value so the submitted jobs get the threads they ask for.
        kwargs["cores"] = 999

    print(f"Bulk2Spot | targets: {' '.join(args.workflow)} | config: {configfile} | "
          f"output: {config.get('outputdir')} | {'cluster' if 'cluster' in kwargs else f'local, {args.cores} core(s)'}"
          f"{' | DRY-RUN' if args.dry_run else ''}", file=sys.stderr)

    ok = snakemake.snakemake(**kwargs)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
