"""ASE template for MACE using the local MatPES r2SCAN OMAT fine-tuned model.

Usage with SOD:

1. Activate the MACE environment:
   conda activate mace

2. In the calculation directory, provide this template as ``template_ase.py``.
   A symlink works well:
   ln -s /home/salvador/coding/sod/sod/tools/template_ase_mace.py template_ase.py

3. Run ``sod setup`` with ``FILER=14`` in INSOD.

Optional override:
   export SOD_MACE_MODEL=/full/path/to/another.model
"""

from __future__ import annotations

import os
from pathlib import Path

from mace.calculators import MACECalculator


MODEL_FILENAME = "MACE-matpes-r2scan-omat-ft.model"


def resolve_model_path() -> Path:
    """Return the MACE model path, favoring an explicit environment override."""

    env_model = os.environ.get("SOD_MACE_MODEL", "").strip()
    if env_model:
        model_path = Path(env_model).expanduser()
        if model_path.is_file():
            return model_path.resolve()
        raise FileNotFoundError(f"SOD_MACE_MODEL does not point to a valid file: {model_path}")

    template_dir = Path(__file__).resolve().parent
    candidates = [
        template_dir / MODEL_FILENAME,
        template_dir / "tools" / MODEL_FILENAME,
        template_dir.parent / "tools" / MODEL_FILENAME,
        Path.cwd() / MODEL_FILENAME,
        Path.cwd() / "tools" / MODEL_FILENAME,
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()

    searched = "\n".join(str(path) for path in candidates)
    raise FileNotFoundError(
        "Could not locate the MACE model file. Set SOD_MACE_MODEL or place "
        f"{MODEL_FILENAME} in one of these locations:\n{searched}"
    )


def build_calculator(_atoms=None) -> MACECalculator:
    """Create the CPU MACE calculator used by the generic ASE driver."""

    model_path = resolve_model_path()
    return MACECalculator(
        model_paths=str(model_path),
        device="cpu",
        default_dtype="float64",
    )


# Default relaxation controls for the generic ASE driver in scripts/vasp2ase.py.
relax_positions = True
relax_cell = False
fmax = 0.01
steps = 500
trajectory = None
