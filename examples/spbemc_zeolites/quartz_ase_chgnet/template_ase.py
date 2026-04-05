import os
from chgnet.model.model import CHGNet
from chgnet.model.dynamics import CHGNetCalculator
MODEL_NAME = os.environ.get("CHGNET_MODEL", "0.3.0")   # "0.3.0" o "r2scan"
def _detect_device():
    if torch.cuda.is_available():
        return "cuda"
    return "cpu"

def build_calculator(_atoms=None):
    device = _detect_device()
    if MODEL_NAME == "0.3.0":
        model = CHGNet.load()
    else:
        model = CHGNet.load(model_name=MODEL_NAME)

    return CHGNetCalculator(
        model=model,
        use_device=device,
        on_isolated_atoms="warn",
    )

relax_positions = True
relax_cell = True
fmax = 0.01
steps = 500
trajectory = None
