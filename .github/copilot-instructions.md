# SOD (Site-Occupancy Disorder) Copilot Instructions

This document provides essential context for AI agents working with the SOD codebase - a Fortran package for modeling periodic systems with site disorder using the supercell ensemble method.

## Project Architecture

### Core Components
- **Combinatorics Engine** (`src/combsod.f90`): Main program identifying inequivalent configurations of site substitutions
- **Statistical Processing** (`src/statsod.f90`, `src/gcstatsod.f90`): Handles canonical/grand-canonical ensemble calculations
- **Space Group Operations** (`sgo/*.sgo`): Library of crystallographic symmetry operators 
- **Helper Programs**:
  - `genersod`: Configuration generator
  - `spbesod`: Energy extrapolation
  - `peaks2spec`: Spectrum analysis

### Key Data Files
- `INSOD`: Input configuration (rigid format, template-driven)
- `OUTSOD`: Generated configuration data
- `EQMATRIX`: Symmetry operation transformation matrices
- `SGO`: Space group operator definitions

## Development Workflow

### Building
```bash
make all  # Compiles all executables to bin/
```
Core dependencies: gfortran compiler

### Project Structure Conventions
- Create new analysis in dedicated directories
- Required files per analysis:
  - `INSOD`: Configuration parameters
  - `SGO`: Symmetry operators (copy from `sgo/` or create new)
  - For GULP integration: `top.gulp` and `bottom.gulp`

### Output Organization
- `CALCS/`: Generated input files
- Best practice: Rename calculation directories as `n01`, `n02` etc. based on substitution count
- Configuration indices in `OUTSOD` map to generated structures

## Integration Points

### External Calculator Integration
- GULP: Requires `top.gulp`/`bottom.gulp` templates
- VASP: Direct generation of input files
- Custom calculators: Extend `src/genersod.f90` for new formats

### File Format Conventions
- Strict format requirements for `INSOD` - maintain exact blank line spacing
- Space group operators: 3x3 matrix + translation vector format
- Configuration files: Index-based referencing between `OUTSOD` and generated structures

## Common Workflows

1. Setup new analysis:
   ```bash
   mkdir analysis_dir && cd analysis_dir
   cp /path/to/template/INSOD .
   cp ../sod/sgo/YOUR_SPACE_GROUP.sgo SGO
   ```

2. Run combinatorial analysis:
   ```bash
   sod_comb.sh  # Generates OUTSOD, EQMATRIX, CALCS/
   ```

3. Process statistics:
   ```bash
   sod_stat.sh  # Requires OUTSOD + temperature data
   ```

## Key Files for Understanding
- `examples/example1/`: Basic perovskite structure example
- `src/combsod.f90`: Core algorithm implementation
- `README.md`: User documentation and basic setup