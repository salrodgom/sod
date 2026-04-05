Quartz example for the LAMMPS workflow.

This example uses `FILER=2`, so `sod setup` prepares and launches LAMMPS directly.

Recommended usage:

```bash
../../../bin/sod setup -N 0
../../../bin/sod setup -N 1
../../../bin/sod setup -N 5
```

The same example also supports direct generation for arbitrary substitution levels with:

```bash
../../../bin/sod comb -N 1
../../../bin/sod comb -N 5
```

By default `scripts/run_jobs.sh` looks for LAMMPS in:

- `$SOD_LAMMPS_EXECUTABLE`
- `/home/salvador/bin/lmp_fftw`
- `lmp_fftw` in `PATH`
- `lmp` in `PATH`

Notes:

- The LAMMPS translation uses the same functional form as the current quartz GULP setup where it is directly representable in the direct `FILER=2` path: Buckingham+Coulomb, shell springs, and an angular term.
- The current direct LAMMPS writer still uses a single oxygen label `O`, so this example is a first transferable quartz setup rather than a full one-to-one reproduction of the finer `O2/O3` typing present in the GULP workflow.
- After `setup`, the level directory contains both `ENERGIES` and `RELAXED_STRUCTURES`.
