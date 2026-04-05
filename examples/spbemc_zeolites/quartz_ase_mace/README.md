Quartz example for the ASE workflow with MACE.

This example uses `FILER=14`, so `sod setup` prepares and launches ASE directly.
The provided template relaxes both atomic positions and the unit cell, with `fmax = 0.001 eV/Ang`, on CPU.

Recommended usage:

```bash
conda activate mace
export MPLCONFIGDIR=/tmp/mpl-quartz-mace
../../../bin/sod setup -N 1
```

If you prefer to bypass activation explicitly:

```bash
export SOD_ASE_PYTHON=/home/salvador/miniforge3/envs/mace/bin/python
export MPLCONFIGDIR=/tmp/mpl-quartz-mace
../../../bin/sod setup -N 1
```

Notes:

- `template_ase.py` uses the MACE model stored in `tools/MACE-matpes-r2scan-omat-ft.model`.
- The `tools` entry in this folder is a symlink to the repository-level `tools/` directory, so the model is not duplicated.
- After `setup`, the level directory contains `ENERGIES` and `RELAXED_STRUCTURES`.
- On CPU this relaxation is not instantaneous. The optimizer progress is written to `n01/c00001.ase.log`, while the engine markers and warnings go to `n01/c00001.out.ase`.
