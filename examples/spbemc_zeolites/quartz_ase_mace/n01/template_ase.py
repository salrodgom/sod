from pathlib import Path
from mace.calculators import MACECalculator
MODEL_PATH = (Path(__file__).resolve().parent.parent / "tools" / "MACE-matpes-r2scan-omat-ft.model").resolve()

def build_calculator(_atoms=None):
    return MACECalculator(
        model_paths=str(MODEL_PATH),
        device="cpu",
        default_dtype="float64",
    )

relax_positions = True
relax_cell = True
fmax = 0.001
steps = 500
trajectory = None
