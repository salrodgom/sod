#!/usr/bin/env python3
"""Tune Metropolis Monte Carlo sample counts for a SOD calculation directory.

The script runs repeated `sod mc -a metropolis` jobs for one level N
and one temperature, increasing `-C` over a user-provided grid. Each run is
executed in an isolated sandbox directory that symlinks the required inputs
from the original calculation folder, so the main folder is not overwritten.

For every run, the script:
  - parses `MC_TRACE.csv`
  - discards the same 20% burn-in used internally by the code
  - estimates the integrated autocorrelation time of the energy trace
  - estimates the effective sample size (ESS)
  - computes block means and a block-to-block energy spread
  - records the runtime and the summary values from `sod_ensemble_summary.csv`

Across replicates for each `-C`, the script reports:
  - mean runtime
  - mean and spread of `<E>`
  - mean and minimum ESS
  - mean acceptance ratio
  - maximum block spread

It then recommends the smallest `-C` that satisfies simple convergence
criteria based on replicate-to-replicate agreement, ESS, and block stability.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import shlex
import shutil
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[1]
DEFAULT_SOD = REPO_ROOT / "bin" / "sod"
DEFAULT_SCRIPTS_DIR = REPO_ROOT / "scripts"
DEFAULT_SAMPLE_GRID = (250, 500, 1000, 2000, 5000, 10000)


@dataclass
class RunMetrics:
    sample_count: int
    replicate: int
    seed: int
    runtime_s: float
    accepted_samples: int
    post_burn_samples: int
    tau_int: float
    ess: float
    block_spread_ev: float
    block_std_ev: float
    e_exp_total_ev: float
    e_min_total_ev: float
    var_total_ev2: float
    delta_exp_total_kjmol: float | None
    acceptance_ratio: float
    run_dir: Path
    trace_file: Path
    summary_file: Path
    visited_energies_file: Path | None


def looks_like_scripts_dir(path: Path) -> bool:
    required = ("run_jobs.sh", "extract.sh", "vasp2gin.sh")
    return path.is_dir() and all((path / item).exists() for item in required)


def infer_scripts_dir(sod_executable: Path) -> Path:
    exe = sod_executable.resolve()
    repo_root = exe.parent.parent
    candidate = repo_root / "scripts"
    if looks_like_scripts_dir(candidate):
        return candidate
    if looks_like_scripts_dir(DEFAULT_SCRIPTS_DIR):
        return DEFAULT_SCRIPTS_DIR
    return candidate


def parse_int_grid(text: str) -> list[int]:
    values: list[int] = []
    for part in text.split(","):
        token = part.strip()
        if not token:
            continue
        try:
            value = int(token)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"Invalid integer in --samples-grid: {token!r}") from exc
        if value <= 0:
            raise argparse.ArgumentTypeError("--samples-grid values must be positive.")
        values.append(value)
    if not values:
        raise argparse.ArgumentTypeError("--samples-grid cannot be empty.")
    return sorted(set(values))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run and analyze repeated Metropolis MC jobs to estimate the minimum "
            "sample count (-C) needed for stable ensemble averages."
        )
    )
    parser.add_argument("calc_dir", type=Path, help="Calculation directory containing INSOD/SGO and nXX folders.")
    parser.add_argument("-N", "--level", type=int, required=True, help="Substitution level to tune.")
    parser.add_argument("-T", "--temperature", type=float, required=True, help="Temperature in Kelvin.")
    parser.add_argument(
        "--samples-grid",
        type=parse_int_grid,
        default=list(DEFAULT_SAMPLE_GRID),
        help="Comma-separated list of -C values to test [250,500,1000,2000,5000,10000].",
    )
    parser.add_argument("--replicates", type=int, default=3, help="Independent runs per -C [3].")
    parser.add_argument("--seed-base", type=int, default=1001, help="Base seed for replicate generation [1001].")
    parser.add_argument(
        "--burn-fraction",
        type=float,
        default=0.20,
        help="Burn-in fraction removed from MC_TRACE analysis [0.20].",
    )
    parser.add_argument(
        "--target-ess",
        type=float,
        default=200.0,
        help="Minimum effective sample size (ESS) required for a run to be considered reliable [200].",
    )
    parser.add_argument(
        "--energy-tolerance-mev",
        type=float,
        default=5.0,
        help="Maximum 95%% CI half-width across replicates for <E>, in meV [5.0].",
    )
    parser.add_argument(
        "--block-tolerance-mev",
        type=float,
        default=10.0,
        help="Maximum within-run block spread accepted, in meV [10.0].",
    )
    parser.add_argument(
        "--blocks",
        type=int,
        default=5,
        help="Number of post-burn blocks used for block-mean diagnostics [5].",
    )
    parser.add_argument(
        "--sod-executable",
        type=Path,
        default=DEFAULT_SOD,
        help=f"Path to sod executable [{DEFAULT_SOD}].",
    )
    parser.add_argument(
        "--scripts-dir",
        type=Path,
        default=None,
        help="Path to the SOD scripts directory. By default it is inferred from --sod-executable.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory where tuning runs and reports are written [calc_dir/metropolis_tuning_*].",
    )
    parser.add_argument(
        "--restart-accept",
        choices=("default", "on", "off"),
        default="default",
        help="Set SOD_FORCE_RESTART_ACCEPT during tuning [default].",
    )
    parser.add_argument(
        "--no-force-mc",
        action="store_true",
        help="Do not append --force-mc automatically.",
    )
    parser.add_argument(
        "--rerun",
        action="store_true",
        help="Delete and rerun existing sample/replicate sandboxes.",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="If gnuplot is available, also render the generated gnuplot script to PNG.",
    )
    args, extra_args = parser.parse_known_args()
    if args.replicates <= 0:
        parser.error("--replicates must be positive.")
    if not (0.0 <= args.burn_fraction < 1.0):
        parser.error("--burn-fraction must lie in [0, 1).")
    if args.blocks <= 0:
        parser.error("--blocks must be positive.")
    if extra_args and extra_args[0] == "--":
        extra_args = extra_args[1:]
    args.extra_args = extra_args
    return args


def ensure_path_exists(path: Path, label: str) -> None:
    if not path.exists():
        raise SystemExit(f"Error: {label} was not found: {path}")


def default_output_dir(calc_dir: Path, level: int, temperature: float) -> Path:
    temp_label = f"{temperature:.2f}".rstrip("0").rstrip(".").replace(".", "p")
    return calc_dir / f"metropolis_tuning_N{level:02d}_T{temp_label}K"


def level_directory_name(level: int) -> str:
    return f"mc0{level}"


def create_sandbox(calc_dir: Path, run_dir: Path, tuning_root: Path) -> None:
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    for entry in calc_dir.iterdir():
        if entry.resolve() == tuning_root.resolve():
            continue
        name = entry.name
        if entry.is_dir():
            if name.startswith("n"):
                (run_dir / name).symlink_to(entry.resolve(), target_is_directory=True)
        else:
            if name in {"INSOD", "SGO", "EQMATRIX"} or entry.suffix in {".gin", ".include", ".lib"}:
                (run_dir / name).symlink_to(entry.resolve())


def read_trace(trace_path: Path) -> list[float]:
    energies: list[float] = []
    with trace_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split(";")
            if len(parts) < 3:
                continue
            energies.append(float(parts[2]))
    return energies


def post_burn_trace(values: list[float], burn_fraction: float) -> list[float]:
    if not values:
        return []
    keep = max(1, math.ceil((1.0 - burn_fraction) * len(values)))
    start = max(0, len(values) - keep)
    return values[start:]


def estimate_tau_int(values: list[float]) -> float:
    n = len(values)
    if n < 4:
        return 1.0
    mean_val = statistics.fmean(values)
    centered = [value - mean_val for value in values]
    var0 = sum(value * value for value in centered) / float(n)
    if var0 <= 0.0:
        return 1.0

    tau = 1.0
    max_lag = min(n // 2, 1000)
    for lag in range(1, max_lag + 1):
        cov = 0.0
        limit = n - lag
        for idx in range(limit):
            cov += centered[idx] * centered[idx + lag]
        cov /= float(limit)
        rho = cov / var0
        if not math.isfinite(rho) or rho <= 0.0:
            break
        tau += 2.0 * rho
    return max(1.0, tau)


def block_diagnostics(values: list[float], block_count: int) -> tuple[float, float]:
    n = len(values)
    if n < 2:
        return 0.0, 0.0
    blocks = min(block_count, n)
    if blocks <= 1:
        return 0.0, 0.0

    block_means: list[float] = []
    start = 0
    for block_idx in range(blocks):
        end = round((block_idx + 1) * n / blocks)
        chunk = values[start:end]
        start = end
        if not chunk:
            continue
        block_means.append(statistics.fmean(chunk))

    if len(block_means) <= 1:
        return 0.0, 0.0
    spread = max(block_means) - min(block_means)
    stddev = statistics.pstdev(block_means)
    return spread, stddev


def parse_summary_row(summary_path: Path, level: int) -> dict[str, str]:
    with summary_path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        header = next(reader)
        header[0] = header[0].lstrip("#")
        for row in reader:
            if not row:
                continue
            if int(row[0]) == level:
                return dict(zip(header, row))
    raise ValueError(f"Level {level} was not found in {summary_path}")


def to_optional_float(text: str | None) -> float | None:
    if text is None:
        return None
    token = text.strip()
    if not token or token == "--":
        return None
    return float(token)


def analyze_run(
    run_dir: Path,
    level: int,
    sample_count: int,
    replicate: int,
    seed: int,
    runtime_s: float,
    burn_fraction: float,
    block_count: int,
) -> RunMetrics:
    level_dir = run_dir / level_directory_name(level)
    trace_path = level_dir / "MC_TRACE.csv"
    summary_path = run_dir / "sod_ensemble_summary.csv"
    visited_energies = level_dir / "ENERGIES"
    ensure_path_exists(trace_path, "MC trace")
    ensure_path_exists(summary_path, "summary CSV")

    trace = read_trace(trace_path)
    if not trace:
        raise ValueError(f"No trace samples found in {trace_path}")
    trace_post = post_burn_trace(trace, burn_fraction)
    tau_int = estimate_tau_int(trace_post)
    ess = float(len(trace_post)) / tau_int if tau_int > 0.0 else float(len(trace_post))
    block_spread, block_std = block_diagnostics(trace_post, block_count)

    if visited_energies.exists():
        convenience = run_dir / "visited_ENERGIES"
        if convenience.exists() or convenience.is_symlink():
            convenience.unlink()
        convenience.symlink_to(visited_energies)
    else:
        visited_energies = None

    row = parse_summary_row(summary_path, level)
    return RunMetrics(
        sample_count=sample_count,
        replicate=replicate,
        seed=seed,
        runtime_s=runtime_s,
        accepted_samples=len(trace),
        post_burn_samples=len(trace_post),
        tau_int=tau_int,
        ess=ess,
        block_spread_ev=block_spread,
        block_std_ev=block_std,
        e_exp_total_ev=float(row["E_exp_total"]),
        e_min_total_ev=float(row["E_min_total"]),
        var_total_ev2=float(row["Var_total"]),
        delta_exp_total_kjmol=to_optional_float(row.get("Delta_exp_total")),
        acceptance_ratio=float(row["Acceptance_ratio"]),
        run_dir=run_dir,
        trace_file=trace_path,
        summary_file=summary_path,
        visited_energies_file=visited_energies,
    )


def run_metropolis_case(
    calc_dir: Path,
    tuning_root: Path,
    sod_executable: Path,
    level: int,
    temperature: float,
    sample_count: int,
    replicate: int,
    seed: int,
    extra_args: list[str],
    force_mc: bool,
    restart_accept: str,
    scripts_dir: Path,
    rerun: bool,
) -> RunMetrics:
    run_dir = tuning_root / f"C{sample_count:06d}" / f"rep{replicate:02d}"
    if rerun and run_dir.exists():
        shutil.rmtree(run_dir)
    if not run_dir.exists():
        create_sandbox(calc_dir, run_dir, tuning_root)

    stdout_path = run_dir / "run.stdout"
    stderr_path = run_dir / "run.stderr"

    command = [
        str(sod_executable.resolve()),
        "mc",
        "-N",
        str(level),
        "-T",
        str(temperature),
        "-C",
        str(sample_count),
        "-s",
        str(seed),
        "-a",
        "metropolis",
    ]
    if force_mc:
        command.append("--force-mc")
    command.extend(extra_args)

    env = os.environ.copy()
    env["SOD_SCRIPTS"] = str(scripts_dir.resolve())
    if restart_accept == "on":
        env["SOD_FORCE_RESTART_ACCEPT"] = "1"
    elif restart_accept == "off":
        env["SOD_FORCE_RESTART_ACCEPT"] = "0"

    started = time.monotonic()
    with stdout_path.open("w", encoding="utf-8") as stdout_handle, stderr_path.open("w", encoding="utf-8") as stderr_handle:
        completed = subprocess.run(
            command,
            cwd=run_dir,
            env=env,
            stdout=stdout_handle,
            stderr=stderr_handle,
            text=True,
            check=False,
        )
    runtime_s = time.monotonic() - started
    if completed.returncode != 0:
        raise RuntimeError(
            f"Metropolis run failed for C={sample_count}, replicate={replicate}. "
            f"See {stdout_path} and {stderr_path}."
        )

    return analyze_run(
        run_dir=run_dir,
        level=level,
        sample_count=sample_count,
        replicate=replicate,
        seed=seed,
        runtime_s=runtime_s,
        burn_fraction=args.burn_fraction,
        block_count=args.blocks,
    )


def mean(values: Iterable[float]) -> float:
    values = list(values)
    if not values:
        return float("nan")
    return statistics.fmean(values)


def stdev(values: Iterable[float]) -> float:
    values = list(values)
    if len(values) <= 1:
        return 0.0
    return statistics.stdev(values)


def write_run_csv(metrics: list[RunMetrics], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "C",
                "replicate",
                "seed",
                "runtime_s",
                "accepted_samples",
                "post_burn_samples",
                "tau_int",
                "ess",
                "block_spread_ev",
                "block_std_ev",
                "E_exp_total_eV",
                "E_min_total_eV",
                "Var_total_eV2",
                "Delta_exp_total_kJmol",
                "Acceptance_ratio",
                "run_dir",
                "trace_file",
                "summary_file",
                "visited_energies_file",
            ]
        )
        for item in metrics:
            writer.writerow(
                [
                    item.sample_count,
                    item.replicate,
                    item.seed,
                    f"{item.runtime_s:.6f}",
                    item.accepted_samples,
                    item.post_burn_samples,
                    f"{item.tau_int:.6f}",
                    f"{item.ess:.6f}",
                    f"{item.block_spread_ev:.8f}",
                    f"{item.block_std_ev:.8f}",
                    f"{item.e_exp_total_ev:.10f}",
                    f"{item.e_min_total_ev:.10f}",
                    f"{item.var_total_ev2:.10f}",
                    "" if item.delta_exp_total_kjmol is None else f"{item.delta_exp_total_kjmol:.10f}",
                    f"{item.acceptance_ratio:.8f}",
                    str(item.run_dir),
                    str(item.trace_file),
                    str(item.summary_file),
                    "" if item.visited_energies_file is None else str(item.visited_energies_file),
                ]
            )


def summarize_by_c(metrics: list[RunMetrics], args: argparse.Namespace) -> tuple[list[dict[str, float | int | bool]], int | None]:
    energy_tol_ev = args.energy_tolerance_mev * 1.0e-3
    block_tol_ev = args.block_tolerance_mev * 1.0e-3
    summaries: list[dict[str, float | int | bool]] = []
    recommended_c: int | None = None

    for sample_count in sorted({item.sample_count for item in metrics}):
        rows = [item for item in metrics if item.sample_count == sample_count]
        e_values = [item.e_exp_total_ev for item in rows]
        ess_values = [item.ess for item in rows]
        runtime_values = [item.runtime_s for item in rows]
        accept_values = [item.acceptance_ratio for item in rows]
        tau_values = [item.tau_int for item in rows]
        block_spreads = [item.block_spread_ev for item in rows]
        missing_visited = any(item.visited_energies_file is None for item in rows)

        e_mean = mean(e_values)
        e_sd = stdev(e_values)
        ci95 = 1.96 * e_sd / math.sqrt(len(rows)) if rows else float("nan")
        ess_mean = mean(ess_values)
        ess_min = min(ess_values) if ess_values else float("nan")
        runtime_mean = mean(runtime_values)
        block_spread_max = max(block_spreads) if block_spreads else float("nan")
        tau_mean = mean(tau_values)
        accept_mean = mean(accept_values)

        stable = (
            ess_min >= args.target_ess
            and ci95 <= energy_tol_ev
            and block_spread_max <= block_tol_ev
        )
        if stable and recommended_c is None:
            recommended_c = sample_count

        summaries.append(
            {
                "C": sample_count,
                "runs": len(rows),
                "runtime_mean_s": runtime_mean,
                "E_exp_mean_eV": e_mean,
                "E_exp_sd_eV": e_sd,
                "E_exp_ci95_eV": ci95,
                "tau_int_mean": tau_mean,
                "ESS_mean": ess_mean,
                "ESS_min": ess_min,
                "acceptance_mean": accept_mean,
                "block_spread_max_eV": block_spread_max,
                "missing_visited_energies": missing_visited,
                "stable": stable,
            }
        )

    return summaries, recommended_c


def write_summary_csv(summaries: list[dict[str, float | int | bool]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "C",
                "runs",
                "runtime_mean_s",
                "E_exp_mean_eV",
                "E_exp_sd_eV",
                "E_exp_ci95_eV",
                "tau_int_mean",
                "ESS_mean",
                "ESS_min",
                "acceptance_mean",
                "block_spread_max_eV",
                "missing_visited_energies",
                "stable",
            ]
        )
        for row in summaries:
            writer.writerow(
                [
                    row["C"],
                    row["runs"],
                    f"{float(row['runtime_mean_s']):.6f}",
                    f"{float(row['E_exp_mean_eV']):.10f}",
                    f"{float(row['E_exp_sd_eV']):.10f}",
                    f"{float(row['E_exp_ci95_eV']):.10f}",
                    f"{float(row['tau_int_mean']):.6f}",
                    f"{float(row['ESS_mean']):.6f}",
                    f"{float(row['ESS_min']):.6f}",
                    f"{float(row['acceptance_mean']):.8f}",
                    f"{float(row['block_spread_max_eV']):.10f}",
                    "yes" if bool(row["missing_visited_energies"]) else "no",
                    "yes" if bool(row["stable"]) else "no",
                ]
            )


def write_recommendation(
    path: Path,
    args: argparse.Namespace,
    summaries: list[dict[str, float | int | bool]],
    recommended_c: int | None,
) -> None:
    energy_tol_ev = args.energy_tolerance_mev * 1.0e-3
    block_tol_ev = args.block_tolerance_mev * 1.0e-3
    with path.open("w", encoding="utf-8") as handle:
        handle.write("Metropolis tuning report\n")
        handle.write("=======================\n\n")
        handle.write(f"Calculation directory : {args.calc_dir.resolve()}\n")
        handle.write(f"sod executable: {args.sod_executable.resolve()}\n")
        handle.write(f"SOD scripts directory : {args.scripts_dir.resolve()}\n")
        handle.write(f"Level N               : {args.level}\n")
        handle.write(f"Temperature (K)       : {args.temperature}\n")
        handle.write(f"Replicates per C      : {args.replicates}\n")
        handle.write(f"Sample grid           : {', '.join(str(x) for x in args.samples_grid)}\n")
        handle.write(f"Burn fraction         : {args.burn_fraction:.3f}\n")
        handle.write(f"Target ESS            : {args.target_ess:.1f}\n")
        handle.write(f"Energy tolerance      : {energy_tol_ev:.6f} eV ({args.energy_tolerance_mev:.3f} meV)\n")
        handle.write(f"Block tolerance       : {block_tol_ev:.6f} eV ({args.block_tolerance_mev:.3f} meV)\n\n")

        if recommended_c is None:
            handle.write("Recommendation: no tested C satisfies all convergence criteria.\n")
            handle.write("Use the largest tested C or extend --samples-grid.\n\n")
        else:
            handle.write(f"Recommendation: use C = {recommended_c}\n")
            handle.write("This is the smallest tested value that satisfied:\n")
            handle.write("  - minimum ESS across replicates >= target ESS\n")
            handle.write("  - 95% CI of <E> across replicates <= energy tolerance\n")
            handle.write("  - maximum within-run block spread <= block tolerance\n\n")

        handle.write("Summary by C\n")
        handle.write("------------\n")
        for row in summaries:
            handle.write(
                "C={C:6d}  runtime={runtime:8.2f}s  <E>={emean: .8f} eV  CI95={ci: .6f} eV  "
                "tau={tau:6.2f}  ESS(min/mean)=({ess_min:6.1f}/{ess_mean:6.1f})  "
                "acc={acc: .4f}  block_spread_max={spread: .6f} eV  stable={stable}\n".format(
                    C=int(row["C"]),
                    runtime=float(row["runtime_mean_s"]),
                    emean=float(row["E_exp_mean_eV"]),
                    ci=float(row["E_exp_ci95_eV"]),
                    tau=float(row["tau_int_mean"]),
                    ess_min=float(row["ESS_min"]),
                    ess_mean=float(row["ESS_mean"]),
                    acc=float(row["acceptance_mean"]),
                    spread=float(row["block_spread_max_eV"]),
                    stable="yes" if bool(row["stable"]) else "no",
                )
            )
        missing_visited = [row for row in summaries if int(row["runs"]) > 0 and bool(row.get("missing_visited_energies", False))]
        if missing_visited:
            handle.write("\nWarnings\n")
            handle.write("--------\n")
            handle.write("Some runs ended without mcNN/ENERGIES. That usually means either:\n")
            handle.write("  - the executable does not include the MC snapshot logic yet, or\n")
            handle.write("  - the GULP/extract stage failed before producing ENERGIES.\n")


def write_gnuplot_script(summary_csv: Path, output_script: Path, output_png: Path) -> None:
    script = f"""set datafile separator ","
set terminal pngcairo size 1400,1000 enhanced font "Helvetica,10"
set output "{output_png.name}"
set multiplot layout 2,2 rowsfirst title "Metropolis tuning summary"

set key left top
set grid
set xlabel "C"
set ylabel "Runtime (s)"
plot "{summary_csv.name}" using 1:3 with linespoints lw 2 pt 7 title "mean runtime"

set xlabel "C"
set ylabel "<E> (eV)"
plot "{summary_csv.name}" using 1:4:6 with yerrorlines lw 2 pt 7 title "<E> +/- CI95"

set xlabel "C"
set ylabel "ESS"
plot "{summary_csv.name}" using 1:8 with linespoints lw 2 pt 7 title "ESS mean", \\
     "{summary_csv.name}" using 1:9 with linespoints lw 2 pt 5 title "ESS min"

set xlabel "C"
set ylabel "Acceptance ratio"
plot "{summary_csv.name}" using 1:10 with linespoints lw 2 pt 7 title "acceptance mean"

unset multiplot
"""
    output_script.write_text(script, encoding="utf-8")


def maybe_render_plot(output_dir: Path, gnuplot_script: Path) -> bool:
    gnuplot = shutil.which("gnuplot")
    if gnuplot is None:
        return False
    completed = subprocess.run(
        [gnuplot, gnuplot_script.name],
        cwd=output_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    return completed.returncode == 0


def print_console_summary(
    args: argparse.Namespace,
    summaries: list[dict[str, float | int | bool]],
    recommended_c: int | None,
    output_dir: Path,
) -> None:
    print(f"Calculation directory : {args.calc_dir.resolve()}")
    print(f"sod exec              : {args.sod_executable.resolve()}")
    print(f"SOD scripts dir       : {args.scripts_dir.resolve()}")
    print(f"Level N               : {args.level}")
    print(f"Temperature (K)       : {args.temperature}")
    print(f"Replicates per C      : {args.replicates}")
    print(f"Sample grid           : {', '.join(str(x) for x in args.samples_grid)}")
    print(f"Reports written in    : {output_dir}")
    print()
    for row in summaries:
        print(
            "C={C:6d}  runtime={runtime:8.2f}s  CI95={ci: .6f} eV  "
            "tau={tau:6.2f}  ESS(min/mean)=({ess_min:6.1f}/{ess_mean:6.1f})  "
            "acc={acc: .4f}  stable={stable}{warn}".format(
                C=int(row["C"]),
                runtime=float(row["runtime_mean_s"]),
                ci=float(row["E_exp_ci95_eV"]),
                tau=float(row["tau_int_mean"]),
                ess_min=float(row["ESS_min"]),
                ess_mean=float(row["ESS_mean"]),
                acc=float(row["acceptance_mean"]),
                stable="yes" if bool(row["stable"]) else "no",
                warn="  [missing visited_ENERGIES]" if bool(row.get("missing_visited_energies", False)) else "",
            )
        )
    print()
    if recommended_c is None:
        print("Recommendation: no tested C satisfied all convergence criteria.")
        print("Try a larger --samples-grid or relax tolerances if needed.")
    else:
        print(f"Recommendation: use C = {recommended_c}")


args = parse_args()


def main() -> int:
    calc_dir = args.calc_dir.resolve()
    ensure_path_exists(calc_dir, "calculation directory")
    ensure_path_exists(args.sod_executable.resolve(), "sod executable")
    ensure_path_exists(calc_dir / "INSOD", "INSOD")
    ensure_path_exists(calc_dir / "SGO", "SGO")
    if args.scripts_dir is None:
        args.scripts_dir = infer_scripts_dir(args.sod_executable)
    args.scripts_dir = args.scripts_dir.resolve()
    if not looks_like_scripts_dir(args.scripts_dir):
        raise SystemExit(
            "Error: SOD scripts directory not found or incomplete: "
            f"{args.scripts_dir}\n"
            "Use --scripts-dir /path/to/sod/scripts."
        )

    tuning_root = args.output_dir.resolve() if args.output_dir is not None else default_output_dir(calc_dir, args.level, args.temperature)
    tuning_root.mkdir(parents=True, exist_ok=True)

    all_metrics: list[RunMetrics] = []
    force_mc = not args.no_force_mc

    for sample_count in args.samples_grid:
        for replicate in range(1, args.replicates + 1):
            seed = args.seed_base + 1000 * replicate + sample_count
            print(
                f"Running metropolis tuning case: N={args.level}, T={args.temperature}, "
                f"C={sample_count}, replicate={replicate}, seed={seed}"
            )
            metrics = run_metropolis_case(
                calc_dir=calc_dir,
                tuning_root=tuning_root,
                sod_executable=args.sod_executable,
                level=args.level,
                temperature=args.temperature,
                sample_count=sample_count,
                replicate=replicate,
                seed=seed,
                extra_args=args.extra_args,
                force_mc=force_mc,
                restart_accept=args.restart_accept,
                scripts_dir=args.scripts_dir,
                rerun=args.rerun,
            )
            all_metrics.append(metrics)

    run_csv = tuning_root / "metropolis_runs.csv"
    summary_csv = tuning_root / "metropolis_tuning_summary.csv"
    recommendation_txt = tuning_root / "metropolis_recommendation.txt"
    gnuplot_script = tuning_root / "metropolis_tuning.gnuplot"
    png_path = tuning_root / "metropolis_tuning.png"

    write_run_csv(all_metrics, run_csv)
    summaries, recommended_c = summarize_by_c(all_metrics, args)
    write_summary_csv(summaries, summary_csv)
    write_recommendation(recommendation_txt, args, summaries, recommended_c)
    write_gnuplot_script(summary_csv, gnuplot_script, png_path)
    if args.plot:
        maybe_render_plot(tuning_root, gnuplot_script)

    print_console_summary(args, summaries, recommended_c, tuning_root)
    print(f"Detailed runs CSV     : {run_csv}")
    print(f"Summary CSV           : {summary_csv}")
    print(f"Recommendation text   : {recommendation_txt}")
    print(f"Gnuplot script        : {gnuplot_script}")
    if args.plot and png_path.exists():
        print(f"PNG figure            : {png_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
