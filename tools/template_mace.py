from pathlib import Path
import torch
from mace.calculators import MACECalculator
MODEL_PATH = (Path(__file__).resolve().parent / "MACE-matpes-r2scan-omat-ft.model").resolve()
def _detect_device():
    if torch.cuda.is_available():
        return "cuda"
    return "cpu"

def build_calculator(_atoms=None):
    device = _detect_device()
    return MACECalculator(
        model_paths=str(MODEL_PATH),
        device=device,
        default_dtype="float64",
    )

relax_positions = True
relax_cell = True

fmax = 0.01
steps = 500
trajectory = None
