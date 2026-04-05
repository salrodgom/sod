# SOD 0.52 - Notes for users

SOD (standing for Site-Occupancy Disorder) is a package of tools for the computer modelling of periodic systems with site disorder, using the supercell ensemble method. 

The package is distributed under the [GPL licence](https://github.com/ypriverol/sod/blob/master/LICENSE.md). 

You can find below the essential info needed to use SOD. Please note that SOD authors can give only limited support to users.


## Functionalities

- Identification of all inequivalent configurations of site substitutions in an arbitrary supercell of an  initial target structure with any space group symmetry.
- Calculation of the degeneracies of configurations.
- Generation of input files for codes like GULP and VASP.
- Simple extrapolation of energies from low to high concentrations within a supercell.
- Statistical mechanics processing of output using either canonical or grand-canonical ensembles.


## Content of the folders

- sod(version)/src contains the source files.
- sod(version)/sgo is a library of space group operators (e.g. 131.sgo contains the operators of the space group 131).
- sod(version)/bin contains the executables. Linux executables are provided here.
- sod(version)/examples contains three examples, based on the cubic perovskite, rutile and rocksalt structures.

## Compiling & installing SOD

- Download the file sod(version).tar.gz (e.g. sod0.52.tar.gz) and copy to a directory, say ROOTSOD:
 
```bash
tar xzvf sod(version).tar.gz
```

- Make compile all the executables into the **bin** folder:
 
```bash 
make all
```

- Add ROOTSOD/sod(version)/bin to your executables path 

```bash 
# add the bin folder to the executables path in your .bashrc file
export PATH=$PATH:ROOTSOD/sod(version)/bin
```

## Unified ensemble workflows

The modern ensemble drivers are available through a unified executable:

```bash
bin/sod <mc|exact|setup|comb|spbe|entropy|compare|eqmatrix> [mode arguments]
bin/sod --mode <mc|exact|setup|comb|spbe|entropy|compare|eqmatrix> [mode arguments]
bin/sod --mode=<mc|exact|setup|comb|spbe|entropy|compare|eqmatrix> [mode arguments]
```

Available modes:

- `mc`: Monte Carlo sampling workflow with optional GULP recalculation and calibration.
- `exact`: exhaustive enumeration of all configurations at the requested substitution levels.
- `setup`: preparation of `n0X` folders, `OUTSOD`, POSCAR files, and external-engine inputs or runs.
- `comb`: pure combinatorial enumeration without external calculations.
- `spbe`: pair-based extrapolation workflow on top of the modern Hamiltonian core.
- `entropy`: configurational entropy from `OUTSOD_Nxxxx`.
- `compare`: comparison of `DeltaG_config` curves between two folders, including a generated `gnuplot` fitting script.
- `eqmatrix`: direct generation of `EQMATRIX` from `INSOD` and `SGO`.

Recommended discovery commands:

```bash
bin/sod --help
bin/sod mc --help
bin/sod exact --help
bin/sod setup --help
bin/sod entropy --help
bin/sod compare --help
```

### `mc` mode

Usage:

```bash
bin/sod mc [-T <K>] [-M <Nmax>] [-C <Nsamples>] [-s <seed>] [-a <sampler>] [--omp|--no-omp] [--force-mc] [-N <spec>]
```

Options:

- `-T, --temperature <K>`: temperature in Kelvin for Boltzmann weights. Default: `1000`.
- `-M, --max-substitutions <N>`: maximum substitution level evaluated when `-N` is absent. Default: `-1` (all levels).
- `-C, --samples <N>`: Monte Carlo samples per level when the combination count exceeds the exact threshold. Default: `5000`.
- `-s, --seed <value>`: random seed. Default: `-1` (derived from `system_clock`).
- `-a, --sampler <uniform|metropolis>`: sampling scheme. Default: `uniform`.
- `--omp`, `--no-omp`: explicitly enable or disable OpenMP.
- `-N <spec>`: level selection. Examples: `-N 5`, `-N 3:8`, `-N 12,30,45`.
- `--parallel-lists`: keep OpenMP enabled even when `-N` is an explicit list.
- `--force-mc`: force Monte Carlo sampling even when exact enumeration would be possible.
- `--template-gin <file>`: append the requested template fragment when generating `.gin` files.
- `--no-template-gin`: skip template fragments when creating `.gin` files. This is the default behavior.

Notes:

- By default every level from `N=0` to `Nmax` is evaluated unless `-N` restricts the set.
- If `C(N,npos) <= 200000`, the level is enumerated exactly; otherwise Monte Carlo sampling is used.
- Results are written to `sod_ensemble_summary.csv` and `sod_ensemble_summary.txt`.

Examples:

```bash
bin/sod mc
bin/sod mc -T 800 -M 6 -C 2000 -s 1234 -a metropolis --omp
bin/sod mc -N 12,30,45
```

### `exact` mode

Usage:

```bash
bin/sod exact [-N <spec>] [-t <tol_eV>] [--just-outsod] [--template-gin <file>]
```

Options:

- `-N <spec>`: requested levels. Examples: `-N -1`, `-N 12`, `-N 1:12`, `-N 3,6,9`.
- `-t, --tolerance <tol>`: energy tolerance in eV used to group minima. Default: `1e-6`.
- `--just-outsod`: write only `xNN/OUTSOD`, skipping `xNN/ENERGIES` and POSCAR files.
- `--template-gin <file>`: append the requested template fragment when generating `.gin` files.
- `--no-template-gin`: skip template fragments when creating `.gin` files. This is the default behavior.

Examples:

```bash
bin/sod exact -N 5:10 -t 1e-5
bin/sod exact -N 8
bin/sod exact -N 3,6,9 --just-outsod
```

### `setup` mode

Usage:

```bash
bin/sod setup -N <spec> [--template-gin <file>]
```

Options:

- `-N, -n <spec>`: levels to prepare. Supports a single value, a range such as `2:6`, or a list such as `0,3,5,7`.
- `--template-gin <file>`: append the requested template fragment when generating `.gin` files.
- `--no-template-gin`: skip template fragments when creating `.gin` files. This is the default behavior.

What it generates:

- `n0X` folders for each selected level.
- `OUTSOD` and POSCAR files inside each folder.
- For `FILER=1`, the GULP execution scripts (`run_jobs.sh`, `extract.sh`, `vasp2gin.sh`, and protocol resources when available).
- For `FILER=2`, the LAMMPS execution scripts plus `template_in.lammps`.
- For `FILER=14`, the ASE execution scripts plus `template_ase.py`.
- An `ENERGIES` file per prepared level after `run_jobs.sh` and `extract.sh` whenever the selected FILER launches an engine.

Example:

```bash
bin/sod setup -N 3,6,9 --template-gin default
```

### `entropy` mode

Usage:

```bash
bin/sod entropy [-N <spec>]
```

Options:

- `-N <spec>`: requested levels. Examples: `-N -1`, `-N 12`, `-N 3:8`, `-N 5,9,11`.

Output:

- Reads `OUTSOD_Nxxxx` from the current folder.
- Writes the aggregated entropy table to `sod_entropy_summary.csv`.

Example:

```bash
bin/sod entropy -N -1
```

### `compare` mode

Usage:

```bash
bin/sod compare --system <folder> --reference <folder> --temperature <K> [options]
```

Required arguments:

- `--system <folder>`: folder containing `sod_ensemble_summary.csv` and `sod_entropy_summary.csv` for the system of interest.
- `--reference <folder>`: reference folder containing `sod_entropy_summary.csv`.
- `-t, --temperature <K>`: temperature used to convert `S_conf` into `TDeltaS`.

Optional arguments:

- `-o, --output-prefix <prefix>`: prefix for the generated `.dat` and `.gnuplot` files.
- `--system-label <text>`: label used in the plots for the system folder.
- `--reference-label <text>`: label used in the plots for the reference folder.

Generated files:

- `<prefix>_system.dat`
- `<prefix>_reference.dat`
- `<prefix>.gnuplot`

The generated `gnuplot` script fits the reference `TDeltaS` curve with order-2/3/4 `x*(1-x)` polynomials and evaluates three `DeltaG` estimates for the study system:

- the total `DeltaG` built from `Delta_exp_total`
- a low-`xY` projection built from `Delta_exp_X_side`
- a high-`xY` projection built from `Delta_exp_Y_side`

This is useful for spotting the intermediate `xY` region where low-branch and high-branch effective-Hamiltonian projections disagree.

Example:

```bash
bin/sod compare --system phase_A --reference phase_B -T 1500 -o compare_1500K
```

### Helper script for `compare`

The repository also provides:

```bash
scripts/prepare_compare_folder.sh [options]
```

This helper prepares the current folder for a later `compare` run by:

1. generating `xNN/OUTSOD` with `sod exact --just-outsod`
2. generating `sod_entropy_summary.csv` with `sod entropy`
3. generating `sod_ensemble_summary.csv` with `sod mc` when needed
4. writing `compare_folder_status.txt`

Options:

- `-N, --levels <spec>`: level specification passed to `exact` and `entropy`. Default: `-1`.
- `-T, --temperature <K>`: temperature for `mc` if `sod_ensemble_summary.csv` must be generated.
- `--sod-bin <path>`: path to the `sod` executable.
- `--label <text>`: optional label stored in `compare_folder_status.txt`.
- `-a, --sampler <name>`: `mc` sampler. Default: `metropolis`.
- `-s, --seed <value>`: `mc` seed. Default: `-1`.
- `--template-gin <file>`: passed through to `sod mc`.
- `--protocol <ver>`: passed through to `sod mc`. Use `1.0` or `2.0` and default to `2.0`.
- `--force-outsod`: regenerate `OUTSOD_Nxxxx` even if some already exist.
- `--force-entropy`: regenerate `sod_entropy_summary.csv`.
- `--force-mc`: regenerate `sod_ensemble_summary.csv` even if it already exists.
- `--skip-mc`: do not run `sod mc`.
- `-h, --help`: show the help message.

Examples:

```bash
scripts/prepare_compare_folder.sh -T 1500
scripts/prepare_compare_folder.sh -T 1500 -N -1 --template-gin template_payload.gin --protocol 2.0
scripts/prepare_compare_folder.sh --sod-bin /path/to/bin/sod -T 1500 --force-outsod --force-mc
```

### Environment variables used by external workflows

- `SOD_GULP_PROTOCOL_VERSION=1.0|2.0`: selects the classic converter-only path or the staged protocol path used by `run_jobs.sh`.
- `SOD_FORCE_RESTART_ACCEPT=<0|1|true|false>`: controls the Metropolis restart-acceptance behaviour in `mc`.
- `SOD_GULP_CPUS` and `SOD_GULP_GLOBAL_LIMIT`: control CPU pinning and the global concurrency limit in `run_jobs.sh`.
- `SOD_LAMMPS_EXECUTABLE=<path>`: selects the LAMMPS executable used when `FILER=2`.
- `SOD_ASE_PYTHON=<path>`: selects the Python interpreter used when `FILER=14`. It must provide ASE and the requested calculator dependencies.

When running long jobs remotely, it can also be useful to enable unbuffered Fortran output:

```bash
export GFORTRAN_UNBUFFERED_ALL=1
```

## Running SOD

- We recommend to create a new folder (say FOLDER_NAME) for each sod project. This will be referred to as the working directory.

- In FOLDER_NAME, you must create a file named *INSOD* which contains all the information for running the combinatorics part of the program. Use the *INSOD* file given in one of the examples as a template. The file is self-explanatory. The format of this file is rigid, so keep the same number of blank lines.

- In FOLDER_NAME, you must also include a file named SGO with the matrix-vector representations of the symmetry operators. First check if your space group is included in the ROOTSOD/sod(version)/sgo library; if this is the case, just copy the file into your working directory, under the name SGO:

```bash
cp ROOTSOD/sod(version)/sgo ./SGO
```

otherwise you have to create the file using the International Tables of Crystallography, or from the website of the Bilbao Crystallographic Server <www.cryst.ehy.es>. The first three numbers in each line are one row of the operator matrix and the fourth number is the component of the operator translation vector.

- If you want to generate Gulp input files for all the independent configurations found by SOD, in addition to setting FILER=1 in the INPUT file, you must provide two files in the working directory:

top.gulp contains the heading of the gulp input file (until the keyword cell).

bottom.gulp contains the tail of the gulp input file (everything after the list of coordinates, including species, potentials, etc).

## Unified workflow

This repository now keeps only the unified `sod` workflow as the supported interface.

- Use `sod comb` for pure combinatorics and generation of inequivalent configurations.
- Use `sod setup` to prepare `nNN` folders, write structures, and run the external workflow that produces `ENERGIES`.
- Use `sod exact` for exhaustive enumeration with the modern effective Hamiltonian.
- Use `sod mc` for stochastic sampling of the physical ensemble.
- Use `sod entropy` to compute configurational entropies from `OUTSOD`.
- Use `sod compare` to compare two systems.
- Use `sod spbe` for the pair-based extrapolation workflow on top of the modern Hamiltonian core.
- Use `sod eqmatrix` to generate `EQMATRIX` directly from `INSOD` and `SGO`.

The legacy standalone frontends and helper scripts (`sod_comb.sh`, `sod_stat.sh`, `sod_gcstat.sh`, `spbesod`, `sod_spbe0.sh`, `sod_spbe1.sh`, and related wrappers) are no longer part of this tree.


## Citing SOD

If you use SOD in your research work, please include a citation to this article:

*Grau-Crespo, R., Hamad, S., Catlow, C. R. A., & De Leeuw, N. H. (2007). Symmetry-adapted configurational modelling of fractional site occupancy in solids. Journal of Physics: Condensed Matter, 19(25), 256201.*
[Original Paper](http://iopscience.iop.org/article/10.1088/0953-8984/19/25/256201/meta) 


Happy SODing!!!

Ricardo Grau-Crespo (r.grau-crespo@reading.ac.uk) and Said Hamad (said@upo.es)
