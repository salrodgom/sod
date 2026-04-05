# Modern SOD Additions Relative to the Public SOD 0.52 Release

This repository preserves the original SOD tools while adding a modern ensemble workflow and several extensions for configurational thermodynamics, GULP integration, and post-processing.

## Highlights

- Added a unified executable, `sod_ensemble`, with dedicated modes for `mc`, `exact`, `setup`, `entropy`, `compare`, and `eqmatrix`.
- Refactored the ensemble implementation into reusable Fortran modules instead of keeping all logic in monolithic drivers.
- Extended Monte Carlo support with `uniform`, `metropolis`, and `parallel tempering` samplers.
- Added explicit level selection by lists and ranges for both Monte Carlo and exact workflows.
- Introduced per-level output folders (`xNN`, `mcNN`) for cleaner workflow organization.

## Effective-Hamiltonian and GULP Workflow Improvements

- Added explicit handling of low-Ge and high-Ge effective-Hamiltonian branches.
- Added GULP-based recalibration of the low-side and high-side expansions.
- Replaced the calibration solver based on normal equations with QR least squares for better numerical robustness.
- Added fitted low/high blend weights (`lambda_high`) derived from GULP reference points.
- Calibration coefficients are now printed during execution and stored in per-level reports.
- Added cumulative CSV summaries for both calibration coefficients and fitted blend overrides.

## New Thermodynamic Tools

- Added a dedicated `entropy` mode to generate `sod_entropy_summary.csv` from `OUTSOD`.
- Added a dedicated `compare` mode to compare configurational free-energy curves between two systems, such as a zeolite and quartz.
- The comparison workflow automatically generates `gnuplot` scripts and supports total, low-side, and high-side `DeltaG` estimates.

## Monte Carlo and Sampling Extensions

- Added force-MC execution even when exact enumeration would still be possible.
- Improved OpenMP handling and level scheduling.
- Added persistent MC monitoring outputs such as `MC_TRACE.csv`, per-temperature traces for tempering, and `monitored_OUTSOD`.
- Added helper tooling to tune Metropolis sampling parameters and collect calibration datasets efficiently.

## Symmetry and EQMATRIX

- Refactored EQMATRIX generation, validation, and export into the symmetry layer.
- Added a dedicated `eqmatrix` mode.
- Added validation for problematic high-side mappings and inconsistent complementary-expansion data.

## Input/Output and Helper Scripts

- Added a `setup` mode to prepare `n0X` folders, `OUTSOD`, POSCAR files, and GULP workflow inputs.
- Added support for both the legacy `vasp2gin.sh` workflow and the newer protocol-based GULP pipeline.
- Improved OSDA handling through `--osda-gin` and `--no-osda-gin`.
- Added helper scripts for preparing `compare` folders, checking high-side `OUTSOD` consistency, analyzing pair/clustering trends, tuning Metropolis runs, and collecting calibration datasets.

## Diagnostics and Usability

- Improved runtime logging and output flushing for long calculations.
- Extended GULP extraction so the final optimization `Gnorm` is retained.
- Improved robustness when partial GULP results are available.
- Updated README documentation in English for the unified ensemble workflow.
