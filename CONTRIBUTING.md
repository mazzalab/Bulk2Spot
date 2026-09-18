# Contributing to Bulk2Spot

Thanks for helping improve Bulk2Spot!

## Reporting problems

Open an issue with:
- what you ran (the `run.py` command and the relevant config sections),
- the error, and the end of the failing rule's log (`<outputdir>/logs/<rule>.log`),
- your platform, and `snakemake --version`.

## Making changes

1. Fork the repository and create a branch: `git checkout -b fix/short-description`.
2. Keep the code simple and in the existing style:
   - one rule per stage in `workflow/rules/`, one R script per rule in `workflow/scripts/`;
   - read parameters from `snakemake@config`, never hard-code paths or column names;
   - put helpers shared by several scripts in `workflow/scripts/common.R`;
   - comment *why*, not *what*.
3. If you add a config option, add it to `config/config.yaml`, `tests/tutorial/config_tutorial.yaml`
   and `docs/configuration.md`.
4. If you change R/Python dependencies, update `workflow/envs/*.yaml` and regenerate the pin file
   from a working environment:
   `conda list -p <env> --explicit > workflow/envs/<name>.linux-64.pin.txt`.
5. Test your change:
   - `./run.py -w report -c tests/tutorial/config_tutorial.yaml -n` (dry-run; CI runs this);
   - a full tutorial run, and compare with the expected results in `tests/tutorial/README.md`.
6. Add a line to `CHANGELOG.md` and open a pull request describing what changed and why.
