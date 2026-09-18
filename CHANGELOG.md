# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [0.1.0] - 2026-09-14

First public version, generalized from an internal GeoMx pipeline.

### Added
- Six-stage workflow: `preprocess`, `spe`, `batch`, `deg`, `deconv`, `report`.
- `run.py` launcher with local, PBS (`-cl`) and custom-scheduler (`--cluster-cmd`) execution,
  retries that request more resources each time, and a per-project run lock.
- Automatic conda environments with tested pin files and a standR post-deploy step.
- Tutorial on NanoString's public WTA kidney dataset, with a download script.
- Documentation: installation, usage, configuration, outputs, methods.
- Example outputs from the tutorial run in `examples/` (report, figures, summary tables).
- Bulk2Spot logo (`docs/images/`), shown in the README and in the HTML report's header and favicon.

### Changed (compared with the internal pipeline)
- All paths are configurable or relative to the repository; there are no site-specific paths.
- Each rule writes a log file to `<outputdir>/logs/`.
- `DEG_GSEA_summary.txt` now includes every comparison column (it used to be overwritten for each
  column), and has a new `Column` field.
- The report's volcano thresholds come from the config; the report no longer claims that batch
  correction was applied when it was not.
- GSVA honours `method` (`ssgsea`, `zscore`, `plage`) with GSVA >= 1.50.
- Per-scan deconvolution barplots use `spe.scan_id_column` instead of a hard-coded `Scan_ID`.
- Output renames: `target_demoData.rds` → `target_data.rds`, `norm_target_demoData.rds` →
  `norm_target_data.rds`, `geomx_report.html` → `bulk2spot_report.html`.
- Removed config keys that were never used: `segment_qc.negGeoMean_threshold`,
  `dimensionality_reduction.n_pca_dimensions`, `deconvolution.negative_flag_column`, `plotting.group`.
