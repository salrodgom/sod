#!/usr/bin/env python3
"""Run an ASE relaxation from a VASP/POSCAR file using a user template.

The template is a regular Python file. It must provide either:

* ``calculator``: a ready-made ASE calculator object, or
* ``build_calculator([atoms])``: a factory that returns a calculator.

Optional template hooks:

* ``prepare_atoms(atoms)``: mutate or replace the ``Atoms`` object before the run.
* ``relax_positions``: bool, defaults to ``True``.
* ``relax_cell``: bool, defaults to ``False``.
* ``optimizer_class``: ASE optimizer class, defaults to ``ase.optimize.BFGS``.
* ``optimizer_kwargs``: dict with extra keyword arguments for the optimizer.
* ``fmax``: convergence threshold in eV/Ang.
* ``steps``: maximum number of optimization steps.
* ``trajectory``: trajectory filename. If omitted, no trajectory is written.
* ``logfile``: optimizer log filename. Defaults to ``"<prefix>.ase.log"``.
"""

from __future__ import annotations

import argparse
import inspect
import runpy
import sys
from pathlib import Path
from typing import Any, Callable

from ase.filters import UnitCellFilter
from ase.io import read, write
from ase.optimize import BFGS


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Read a VASP/POSCAR structure, attach an ASE calculator from a user "
            "template, relax it, and print a parseable final energy marker."
        )
    )
    parser.add_argument("structure", help="Input structure readable by ASE, typically cXXXXX.vasp.")
    parser.add_argument(
        "template",
        nargs="?",
        default="template_ase.py",
        help="Python template that constructs the ASE calculator [template_ase.py].",
    )
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Prefix for relaxed outputs. Defaults to the input filename stem.",
    )
    parser.add_argument(
        "--fmax",
        type=float,
        default=0.01,
        help="Default force threshold in eV/Ang when the template does not override it [0.01].",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=500,
        help="Default maximum number of optimizer steps when the template does not override it [500].",
    )
    return parser


def call_maybe_with_atoms(func: Callable[..., Any], atoms: Any) -> Any:
    """Call a user hook either with zero arguments or with the Atoms object."""

    signature = inspect.signature(func)
    if len(signature.parameters) == 0:
        return func()
    return func(atoms)


def load_template(template_path: Path) -> dict[str, Any]:
    if not template_path.is_file():
        raise FileNotFoundError(f"Template file not found: {template_path}")
    return runpy.run_path(str(template_path))


def resolve_calculator(namespace: dict[str, Any], atoms: Any) -> Any:
    if "calculator" in namespace:
        return namespace["calculator"]
    if "build_calculator" in namespace:
        return call_maybe_with_atoms(namespace["build_calculator"], atoms)
    raise RuntimeError(
        "The ASE template must define either a 'calculator' object or a "
        "'build_calculator' factory."
    )


def maybe_prepare_atoms(namespace: dict[str, Any], atoms: Any) -> Any:
    if "prepare_atoms" not in namespace:
        return atoms
    prepared = call_maybe_with_atoms(namespace["prepare_atoms"], atoms)
    return atoms if prepared is None else prepared


def main() -> int:
    args = build_parser().parse_args()
    structure_path = Path(args.structure)
    template_path = Path(args.template)

    atoms = read(structure_path)
    namespace = load_template(template_path)
    atoms = maybe_prepare_atoms(namespace, atoms)
    atoms.calc = resolve_calculator(namespace, atoms)

    output_prefix = args.output_prefix or structure_path.stem
    relax_positions = bool(namespace.get("relax_positions", True))
    relax_cell = bool(namespace.get("relax_cell", False))
    fmax = float(namespace.get("fmax", args.fmax))
    steps = int(namespace.get("steps", args.steps))
    optimizer_class = namespace.get("optimizer_class", BFGS)
    optimizer_kwargs = dict(namespace.get("optimizer_kwargs", {}))
    trajectory = namespace.get("trajectory", None)
    logfile = namespace.get("logfile", f"{output_prefix}.ase.log")

    if trajectory is not None:
        optimizer_kwargs.setdefault("trajectory", trajectory)
    optimizer_kwargs.setdefault("logfile", logfile)

    print(f"SOD_ASE_STRUCTURE {structure_path}", flush=True)
    print(f"SOD_ASE_TEMPLATE {template_path}", flush=True)
    print(f"SOD_ASE_RELAX_POSITIONS {int(relax_positions)}", flush=True)
    print(f"SOD_ASE_RELAX_CELL {int(relax_cell)}", flush=True)
    print(f"SOD_ASE_FMAX {fmax:.8f}", flush=True)
    print(f"SOD_ASE_STEPS {steps}", flush=True)

    if relax_positions or relax_cell:
        optimizer_target = atoms
        if relax_cell:
            optimizer_target = UnitCellFilter(atoms)
        optimizer = optimizer_class(optimizer_target, **optimizer_kwargs)
        optimizer.run(fmax=fmax, steps=steps)

    final_energy = atoms.get_potential_energy()
    relaxed_cif = f"{output_prefix}.relaxed.cif"
    write(relaxed_cif, atoms)

    print(f"SOD_FINAL_ENERGY {final_energy:.16f}", flush=True)
    print(f"SOD_FINAL_STRUCTURE {relaxed_cif}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - failure path for external engines
        print(f"SOD_ASE_ERROR {exc}", file=sys.stderr, flush=True)
        raise
