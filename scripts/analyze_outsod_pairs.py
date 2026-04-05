#!/usr/bin/env python3
"""Analyze pair-energy correlations from a SOD OUTSOD/ENERGIES level directory.

The main outputs are expressed in terms of substituted-species pair distances
in Angstrom and distance-cutoff cluster descriptors. Legacy shell-based CSV
output is still written when geometric data are available, but the terminal
summary and plots focus on distance bins and cutoffs because they are easier
to interpret.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Iterable


ASE_FALLBACK_PYTHON = Path("/home/salvador/miniforge3/bin/python")
COORD_TOL = 1.0e-3


def ensure_ase_runtime() -> None:
    try:
        import ase  # noqa: F401
    except ModuleNotFoundError as exc:
        fallback = ASE_FALLBACK_PYTHON
        current = Path(sys.executable).resolve()
        if fallback.exists() and current != fallback:
            os.execv(str(fallback), [str(fallback), str(Path(__file__).resolve()), *sys.argv[1:]])
        raise SystemExit(
            "Error: ASE is not available. Run this script with "
            f"{fallback} or another Python environment that provides ASE."
        ) from exc


ensure_ase_runtime()

import numpy as np
from ase import Atoms
from ase.cell import Cell
from ase.io import read


def parse_outsod(path: Path) -> tuple[int, int, list[int], list[tuple[int, ...]]]:
    with path.open("r", encoding="utf-8") as handle:
        header1 = handle.readline().split()
        if len(header1) < 5:
            raise ValueError(f"Invalid OUTSOD header in {path}")
        substitutions = int(header1[0])
        total_sites = int(header1[3])

        header2 = handle.readline().split()
        if len(header2) < 2:
            raise ValueError(f"Invalid OUTSOD configuration count in {path}")
        expected_configs = int(header2[0])

        degeneracies: list[int] = []
        subsets: list[tuple[int, ...]] = []
        for line in handle:
            parts = line.split()
            if not parts:
                continue
            degeneracies.append(int(parts[1]))
            subsets.append(tuple(int(value) for value in parts[2:]))

    return substitutions, total_sites, degeneracies, subsets


def parse_energies(path: Path) -> list[float]:
    energies: list[float] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if not parts:
                continue
            energies.append(float(parts[0]))
    return energies


def weighted_mean(values: Iterable[float], weights: Iterable[int]) -> float:
    num = 0.0
    den = 0.0
    for value, weight in zip(values, weights):
        num += float(value) * float(weight)
        den += float(weight)
    if den == 0.0:
        raise ValueError("Weighted mean with zero total weight")
    return num / den


def weighted_linear_stats(xvals: list[int], yvals: list[float], weights: list[int]) -> tuple[float, float]:
    if len(xvals) != len(yvals) or len(xvals) != len(weights):
        raise ValueError("Incompatible lengths in weighted_linear_stats")

    x_mean = weighted_mean(xvals, weights)
    y_mean = weighted_mean(yvals, weights)
    sum_w = float(sum(weights))
    var_x = sum(float(w) * (float(x) - x_mean) ** 2 for x, w in zip(xvals, weights)) / sum_w
    var_y = sum(float(w) * (float(y) - y_mean) ** 2 for y, w in zip(yvals, weights)) / sum_w
    if var_x <= 0.0 or var_y <= 0.0:
        return float("nan"), float("nan")

    cov_xy = (
        sum(float(w) * (float(x) - x_mean) * (float(y) - y_mean) for x, y, w in zip(xvals, yvals, weights))
        / sum_w
    )
    slope = cov_xy / var_x
    corr = cov_xy / (var_x * var_y) ** 0.5
    return slope, corr


def next_data_line(lines: list[str], start: int) -> tuple[int, str]:
    index = start
    while index < len(lines):
        stripped = lines[index].strip()
        if stripped and not stripped.startswith("#"):
            return index, stripped
        index += 1
    raise ValueError("Unexpected end of INSOD while reading data.")


def parse_insod(path: Path) -> dict[str, object]:
    lines = path.read_text(encoding="utf-8").splitlines()
    idx = 0

    idx, _ = next_data_line(lines, idx)
    idx += 1  # title

    idx, line = next_data_line(lines, idx)
    cell_params = [float(value) for value in line.split()[:6]]
    idx += 1

    idx, line = next_data_line(lines, idx)
    nsp = int(line.split()[0])
    idx += 1

    idx, line = next_data_line(lines, idx)
    atom_types = line.split()
    idx += 1

    idx, line = next_data_line(lines, idx)
    natsp0 = [int(value) for value in line.split()[:nsp]]
    nat0 = sum(natsp0)
    idx += 1

    coords0: list[list[float]] = []
    for _ in range(nat0):
        idx, line = next_data_line(lines, idx)
        coords0.append([float(value) for value in line.split()[:3]])
        idx += 1

    idx, line = next_data_line(lines, idx)
    na, nb, nc = (int(value) for value in line.split()[:3])
    idx += 1

    idx, line = next_data_line(lines, idx)
    sptarget = int(line.split()[0])

    spat0: list[int] = []
    for species_idx, count in enumerate(natsp0, start=1):
        spat0.extend([species_idx] * count)

    return {
        "cell_params": cell_params,
        "atom_types": atom_types,
        "natsp0": natsp0,
        "coords0": np.array(coords0, dtype=float),
        "spat0": np.array(spat0, dtype=int),
        "na": na,
        "nb": nb,
        "nc": nc,
        "sptarget": sptarget,
    }


def parse_sgo(path: Path) -> list[tuple[np.ndarray, np.ndarray]]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    ops: list[tuple[np.ndarray, np.ndarray]] = []
    idx = 1
    while idx < len(lines):
        if lines[idx] == "0":
            break
        int(lines[idx])  # operator label, used only for validation
        idx += 1
        rotation = np.zeros((3, 3), dtype=float)
        translation = np.zeros(3, dtype=float)
        for row in range(3):
            parts = lines[idx].split()
            idx += 1
            rotation[row, :] = [float(parts[0]), float(parts[1]), float(parts[2])]
            translation[row] = float(parts[3])
        ops.append((rotation, translation))
    if not ops:
        raise ValueError(f"No symmetry operators found in {path}")
    return ops


def wrap_fractional_vector(vector: np.ndarray) -> np.ndarray:
    wrapped = np.array(vector, dtype=float, copy=True)
    wrapped -= np.floor(wrapped)
    wrapped[np.isclose(wrapped, 1.0, atol=1.0e-12)] = 0.0
    return wrapped


def build_substitutable_site_geometry(structure_dir: Path) -> tuple[list[float], np.ndarray]:
    insod = parse_insod(structure_dir / "INSOD")
    ops = parse_sgo(structure_dir / "SGO")

    coords0 = insod["coords0"]
    spat0 = insod["spat0"]
    na = int(insod["na"])
    nb = int(insod["nb"])
    nc = int(insod["nc"])
    sptarget = int(insod["sptarget"])

    coords1r: list[np.ndarray] = []
    spat1r: list[int] = []
    for atom_index in range(len(coords0)):
        for rotation, translation in ops:
            coords1r.append(wrap_fractional_vector(rotation @ coords0[atom_index] + translation))
            spat1r.append(int(spat0[atom_index]))

    coords1: list[np.ndarray] = []
    spat_unit: list[int] = []
    for coord, species in zip(coords1r, spat1r):
        found = False
        for coord_ref in coords1:
            diff = coord - coord_ref
            if float(np.dot(diff, diff)) <= COORD_TOL:
                found = True
                break
        if not found:
            coords1.append(coord)
            spat_unit.append(species)

    vt: list[np.ndarray] = []
    for ia in range(na):
        for ib in range(nb):
            for ic in range(nc):
                vt.append(np.array([ia / float(na), ib / float(nb), ic / float(nc)], dtype=float))

    coords_super: list[np.ndarray] = []
    spat_super: list[int] = []
    for coord, species in zip(coords1, spat_unit):
        for shift in vt:
            coords_super.append(
                wrap_fractional_vector(
                    shift + np.array([coord[0] / na, coord[1] / nb, coord[2] / nc], dtype=float)
                )
            )
            spat_super.append(species)

    target_coords = np.array(
        [coord for coord, species in zip(coords_super, spat_super) if species == sptarget],
        dtype=float,
    )

    cell_params = list(insod["cell_params"])
    cell_params[0] *= float(na)
    cell_params[1] *= float(nb)
    cell_params[2] *= float(nc)
    return cell_params, target_coords


def assign_shells_from_pair_distances(
    pair_distances: dict[tuple[int, int], float], shell_tol: float
) -> tuple[dict[tuple[int, int], int], dict[int, float]]:
    sorted_pairs = sorted(pair_distances.items(), key=lambda item: item[1])
    pair_shells: dict[tuple[int, int], int] = {}
    shell_distances: dict[int, float] = {}
    shell_members: list[list[float]] = []
    shell_index = 0
    last_distance = None
    for pair, distance in sorted_pairs:
        if last_distance is None or abs(distance - last_distance) > shell_tol:
            shell_index += 1
            shell_members.append([])
            last_distance = distance
        pair_shells[pair] = shell_index
        shell_members[shell_index - 1].append(distance)

    for idx, distances in enumerate(shell_members, start=1):
        shell_distances[idx] = sum(distances) / float(len(distances))

    return pair_shells, shell_distances


def build_pair_geometry(site_coords: np.ndarray, cell_params: list[float], shell_tol: float) -> tuple[dict[tuple[int, int], float], dict[tuple[int, int], int], dict[int, float]]:
    atoms = Atoms(
        symbols=["H"] * len(site_coords),
        cell=Cell.fromcellpar(cell_params),
        scaled_positions=site_coords,
        pbc=True,
    )

    pair_distances: dict[tuple[int, int], float] = {}
    for idx_i, idx_j in combinations(range(len(site_coords)), 2):
        pair = (idx_i + 1, idx_j + 1)
        distance = float(atoms.get_distance(idx_i, idx_j, mic=True))
        pair_distances[pair] = distance
    pair_shells, shell_distances = assign_shells_from_pair_distances(pair_distances, shell_tol)
    return pair_distances, pair_shells, shell_distances


def load_optimized_pair_distances(
    level_dir: Path, subsets: list[tuple[int, ...]], degeneracies: list[int], y_symbol: str
) -> tuple[dict[tuple[int, int], float], int, int]:
    weighted_sums: dict[tuple[int, int], float] = {}
    weighted_counts: dict[tuple[int, int], int] = {}
    found_cifs = 0
    used_cifs = 0

    for cfg_index, (subset, degeneracy) in enumerate(zip(subsets, degeneracies), start=1):
        cif_path = level_dir / f"c{cfg_index:05d}.vasp.cif"
        if not cif_path.is_file():
            continue

        found_cifs += 1
        atoms = read(cif_path)
        symbols = atoms.get_chemical_symbols()
        effective_y_symbol = y_symbol
        if not effective_y_symbol:
            symbol_counts = Counter(symbols)
            matching_symbols = [
                symbol
                for symbol, count in symbol_counts.items()
                if count == len(subset) and symbol.upper() not in {"O", "H"}
            ]
            if len(matching_symbols) == 1:
                effective_y_symbol = matching_symbols[0]
        if not effective_y_symbol:
            continue

        y_indices = [idx for idx, symbol in enumerate(symbols) if symbol == effective_y_symbol]
        if len(y_indices) != len(subset):
            continue

        used_cifs += 1
        sorted_subset = tuple(sorted(subset))
        for local_i, local_j in combinations(range(len(sorted_subset)), 2):
            pair = (sorted_subset[local_i], sorted_subset[local_j])
            distance = float(atoms.get_distance(y_indices[local_i], y_indices[local_j], mic=True))
            weighted_sums[pair] = weighted_sums.get(pair, 0.0) + float(degeneracy) * distance
            weighted_counts[pair] = weighted_counts.get(pair, 0) + int(degeneracy)

    pair_distances = {
        pair: weighted_sums[pair] / float(weighted_counts[pair]) for pair in weighted_sums if weighted_counts[pair] > 0
    }
    return pair_distances, found_cifs, used_cifs


def summarize_shells(
    subsets: list[tuple[int, ...]],
    energies: list[float],
    degeneracies: list[int],
    pair_shells: dict[tuple[int, int], int],
    shell_distances: dict[int, float],
    weighted_global_mean: float,
) -> list[dict[str, float]]:
    shell_count = len(shell_distances)
    config_shell_counts = np.zeros((len(subsets), shell_count), dtype=int)

    for cfg_idx, subset in enumerate(subsets):
        for pair in combinations(sorted(subset), 2):
            shell_idx = pair_shells.get(pair)
            if shell_idx is not None:
                config_shell_counts[cfg_idx, shell_idx - 1] += 1

    shell_rows: list[dict[str, float]] = []
    for shell_idx in range(1, shell_count + 1):
        counts = config_shell_counts[:, shell_idx - 1].tolist()
        present_indices = [idx for idx, count in enumerate(counts) if count > 0]
        if present_indices:
            present_energies = [energies[idx] for idx in present_indices]
            present_weights = [degeneracies[idx] for idx in present_indices]
            presence_mean = weighted_mean(present_energies, present_weights)
            weighted_configs = sum(present_weights)
        else:
            presence_mean = float("nan")
            weighted_configs = 0

        slope, corr = weighted_linear_stats(counts, energies, degeneracies)
        shell_rows.append(
            {
                "shell_index": shell_idx,
                "shell_distance_ang": shell_distances[shell_idx],
                "pairs_in_shell": sum(1 for value in pair_shells.values() if value == shell_idx),
                "weighted_configs_with_shell": weighted_configs,
                "presence_delta_eV": presence_mean - weighted_global_mean if present_indices else float("nan"),
                "slope_eV_per_pair": slope,
                "weighted_corr": corr,
            }
        )

    return shell_rows


def summarize_radial_bins(
    subsets: list[tuple[int, ...]],
    energies: list[float],
    degeneracies: list[int],
    pair_distances: dict[tuple[int, int], float],
    weighted_global_mean: float,
    bin_width: float,
) -> list[dict[str, float]]:
    if not pair_distances:
        return []

    min_dist = min(pair_distances.values())
    max_dist = max(pair_distances.values())
    start = math.floor(min_dist / bin_width) * bin_width
    end = math.ceil(max_dist / bin_width) * bin_width
    if end <= start:
        end = start + bin_width

    nbins = max(1, int(math.ceil((end - start) / bin_width)))
    bin_pairs: list[list[tuple[int, int]]] = [[] for _ in range(nbins)]
    pair_to_bin: dict[tuple[int, int], int] = {}

    for pair, distance in pair_distances.items():
        bin_idx = int((distance - start) / bin_width)
        if bin_idx >= nbins:
            bin_idx = nbins - 1
        if bin_idx < 0:
            bin_idx = 0
        pair_to_bin[pair] = bin_idx
        bin_pairs[bin_idx].append(pair)

    config_bin_counts = np.zeros((len(subsets), nbins), dtype=int)
    for cfg_idx, subset in enumerate(subsets):
        for pair in combinations(sorted(subset), 2):
            bin_idx = pair_to_bin.get(pair)
            if bin_idx is not None:
                config_bin_counts[cfg_idx, bin_idx] += 1

    total_weighted_pair_occurrences = 0.0
    for bin_idx in range(nbins):
        total_weighted_pair_occurrences += sum(
            float(degeneracies[cfg_idx]) * float(config_bin_counts[cfg_idx, bin_idx]) for cfg_idx in range(len(subsets))
        )

    radial_rows: list[dict[str, float]] = []
    for bin_idx in range(nbins):
        lower = start + float(bin_idx) * bin_width
        upper = lower + bin_width
        center = 0.5 * (lower + upper)
        counts = config_bin_counts[:, bin_idx].tolist()
        present_indices = [idx for idx, count in enumerate(counts) if count > 0]
        if present_indices:
            present_energies = [energies[idx] for idx in present_indices]
            present_weights = [degeneracies[idx] for idx in present_indices]
            presence_mean = weighted_mean(present_energies, present_weights)
            weighted_configs = sum(present_weights)
        else:
            presence_mean = float("nan")
            weighted_configs = 0

        weighted_occurrences = sum(
            float(degeneracies[cfg_idx]) * float(config_bin_counts[cfg_idx, bin_idx]) for cfg_idx in range(len(subsets))
        )
        slope, corr = weighted_linear_stats(counts, energies, degeneracies)
        radial_rows.append(
            {
                "r_lower_ang": lower,
                "r_upper_ang": upper,
                "r_center_ang": center,
                "pairs_in_bin": len(bin_pairs[bin_idx]),
                "weighted_pair_occurrences": weighted_occurrences,
                "weighted_pair_fraction": (
                    weighted_occurrences / total_weighted_pair_occurrences if total_weighted_pair_occurrences > 0.0 else 0.0
                ),
                "weighted_configs_with_bin": weighted_configs,
                "presence_delta_eV": presence_mean - weighted_global_mean if present_indices else float("nan"),
                "slope_eV_per_pair": slope,
                "weighted_corr": corr,
            }
        )

    return radial_rows


def summarize_radial_curves(
    subsets: list[tuple[int, ...]],
    energies: list[float],
    degeneracies: list[int],
    pair_distances: dict[tuple[int, int], float],
    bin_width: float,
    low_energy_fraction: float,
) -> tuple[list[dict[str, float]], float]:
    if not pair_distances:
        return [], float("nan")

    min_dist = min(pair_distances.values())
    max_dist = max(pair_distances.values())
    start = math.floor(min_dist / bin_width) * bin_width
    end = math.ceil(max_dist / bin_width) * bin_width
    if end <= start:
        end = start + bin_width

    nbins = max(1, int(math.ceil((end - start) / bin_width)))
    bin_pairs: list[list[tuple[int, int]]] = [[] for _ in range(nbins)]
    pair_to_bin: dict[tuple[int, int], int] = {}
    for pair, distance in pair_distances.items():
        bin_idx = int((distance - start) / bin_width)
        if bin_idx >= nbins:
            bin_idx = nbins - 1
        if bin_idx < 0:
            bin_idx = 0
        pair_to_bin[pair] = bin_idx
        bin_pairs[bin_idx].append(pair)

    config_bin_counts = np.zeros((len(subsets), nbins), dtype=int)
    for cfg_idx, subset in enumerate(subsets):
        for pair in combinations(sorted(subset), 2):
            bin_idx = pair_to_bin.get(pair)
            if bin_idx is not None:
                config_bin_counts[cfg_idx, bin_idx] += 1

    total_pair_types = float(len(pair_distances))
    total_all_occurrences = sum(
        float(degeneracies[cfg_idx]) * float(config_bin_counts[cfg_idx, bin_idx])
        for cfg_idx in range(len(subsets))
        for bin_idx in range(nbins)
    )

    cfg_order = sorted(range(len(subsets)), key=lambda idx: energies[idx])
    total_deg = float(sum(degeneracies))
    target_low_deg = max(1.0, min(1.0, low_energy_fraction) * total_deg)
    low_indices: list[int] = []
    cumulative_deg = 0.0
    low_energy_threshold = float("nan")
    for idx in cfg_order:
        low_indices.append(idx)
        cumulative_deg += float(degeneracies[idx])
        low_energy_threshold = energies[idx]
        if cumulative_deg >= target_low_deg:
            break
    low_index_set = set(low_indices)

    total_low_occurrences = sum(
        float(degeneracies[cfg_idx]) * float(config_bin_counts[cfg_idx, bin_idx])
        for cfg_idx in low_indices
        for bin_idx in range(nbins)
    )

    rows: list[dict[str, float]] = []
    for bin_idx in range(nbins):
        lower = start + float(bin_idx) * bin_width
        upper = lower + bin_width
        center = 0.5 * (lower + upper)
        n_pair_types = len(bin_pairs[bin_idx])
        random_pair_fraction = n_pair_types / total_pair_types if total_pair_types > 0.0 else float("nan")
        all_occurrences = sum(
            float(degeneracies[cfg_idx]) * float(config_bin_counts[cfg_idx, bin_idx]) for cfg_idx in range(len(subsets))
        )
        low_occurrences = sum(
            float(degeneracies[cfg_idx]) * float(config_bin_counts[cfg_idx, bin_idx]) for cfg_idx in low_index_set
        )
        all_pair_fraction = all_occurrences / total_all_occurrences if total_all_occurrences > 0.0 else float("nan")
        low_pair_fraction = low_occurrences / total_low_occurrences if total_low_occurrences > 0.0 else float("nan")

        rows.append(
            {
                "r_lower_ang": lower,
                "r_upper_ang": upper,
                "r_center_ang": center,
                "n_lattice_pairs": n_pair_types,
                "random_pair_fraction": random_pair_fraction,
                "all_weighted_pair_occurrences": all_occurrences,
                "all_pair_fraction": all_pair_fraction,
                "low_weighted_pair_occurrences": low_occurrences,
                "low_pair_fraction": low_pair_fraction,
                "all_over_random": (
                    all_pair_fraction / random_pair_fraction
                    if random_pair_fraction and not math.isnan(random_pair_fraction)
                    else float("nan")
                ),
                "low_over_random": (
                    low_pair_fraction / random_pair_fraction
                    if random_pair_fraction and not math.isnan(random_pair_fraction)
                    else float("nan")
                ),
            }
        )

    return rows, low_energy_threshold


def count_cluster_sizes(
    subset: tuple[int, ...],
    pair_shells: dict[tuple[int, int], int],
    shell_cutoff: int,
    cluster_max_size: int,
) -> dict[int, int]:
    if len(subset) < 2:
        return {size: 0 for size in range(2, cluster_max_size + 1)}

    neighbors = {site: set() for site in subset}
    for pair in combinations(sorted(subset), 2):
        shell_idx = pair_shells.get(pair)
        if shell_idx is not None and shell_idx <= shell_cutoff:
            neighbors[pair[0]].add(pair[1])
            neighbors[pair[1]].add(pair[0])

    visited: set[int] = set()
    counts = {size: 0 for size in range(2, cluster_max_size + 1)}
    for site in subset:
        if site in visited:
            continue
        stack = [site]
        visited.add(site)
        component_size = 0
        while stack:
            current = stack.pop()
            component_size += 1
            for neighbor in neighbors[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        if 2 <= component_size <= cluster_max_size:
            counts[component_size] += 1

    return counts


def count_cluster_sizes_distance(
    subset: tuple[int, ...],
    pair_distances: dict[tuple[int, int], float],
    cutoff_distance: float,
    cluster_max_size: int,
) -> dict[int, int]:
    if len(subset) < 2:
        return {size: 0 for size in range(2, cluster_max_size + 1)}

    neighbors = {site: set() for site in subset}
    for pair in combinations(sorted(subset), 2):
        distance = pair_distances.get(pair)
        if distance is not None and distance <= cutoff_distance:
            neighbors[pair[0]].add(pair[1])
            neighbors[pair[1]].add(pair[0])

    visited: set[int] = set()
    counts = {size: 0 for size in range(2, cluster_max_size + 1)}
    for site in subset:
        if site in visited:
            continue
        stack = [site]
        visited.add(site)
        component_size = 0
        while stack:
            current = stack.pop()
            component_size += 1
            for neighbor in neighbors[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        if 2 <= component_size <= cluster_max_size:
            counts[component_size] += 1

    return counts


def summarize_clusters(
    subsets: list[tuple[int, ...]],
    energies: list[float],
    degeneracies: list[int],
    pair_shells: dict[tuple[int, int], int],
    shell_distances: dict[int, float],
    weighted_global_mean: float,
    cluster_shell_max: int,
    cluster_max_size: int,
) -> list[dict[str, float]]:
    if not pair_shells or not shell_distances:
        return []

    shell_cap = min(cluster_shell_max, len(shell_distances))
    cluster_rows: list[dict[str, float]] = []
    for shell_cutoff in range(1, shell_cap + 1):
        motif_series = {size: [] for size in range(2, cluster_max_size + 1)}
        for subset in subsets:
            counts = count_cluster_sizes(subset, pair_shells, shell_cutoff, cluster_max_size)
            for size in range(2, cluster_max_size + 1):
                motif_series[size].append(counts[size])

        for size in range(2, cluster_max_size + 1):
            counts = motif_series[size]
            present_indices = [idx for idx, count in enumerate(counts) if count > 0]
            if present_indices:
                present_energies = [energies[idx] for idx in present_indices]
                present_weights = [degeneracies[idx] for idx in present_indices]
                presence_mean = weighted_mean(present_energies, present_weights)
                weighted_configs = sum(present_weights)
            else:
                presence_mean = float("nan")
                weighted_configs = 0

            mean_count = weighted_mean(counts, degeneracies)
            slope, corr = weighted_linear_stats(counts, energies, degeneracies)
            cluster_rows.append(
                {
                    "shell_cutoff": shell_cutoff,
                    "shell_cutoff_distance_ang": shell_distances[shell_cutoff],
                    "cluster_size": size,
                    "weighted_configs_with_cluster": weighted_configs,
                    "weighted_mean_cluster_count": mean_count,
                    "presence_delta_eV": presence_mean - weighted_global_mean if present_indices else float("nan"),
                    "slope_eV_per_cluster": slope,
                    "weighted_corr": corr,
                }
            )

    return cluster_rows


def summarize_clusters_distance(
    subsets: list[tuple[int, ...]],
    energies: list[float],
    degeneracies: list[int],
    pair_distances: dict[tuple[int, int], float],
    cutoff_distances: list[float],
    weighted_global_mean: float,
    cluster_max_size: int,
) -> list[dict[str, float]]:
    if not pair_distances or not cutoff_distances:
        return []

    cluster_rows: list[dict[str, float]] = []
    for cutoff_distance in cutoff_distances:
        motif_series = {size: [] for size in range(2, cluster_max_size + 1)}
        for subset in subsets:
            counts = count_cluster_sizes_distance(subset, pair_distances, cutoff_distance, cluster_max_size)
            for size in range(2, cluster_max_size + 1):
                motif_series[size].append(counts[size])

        for size in range(2, cluster_max_size + 1):
            counts = motif_series[size]
            present_indices = [idx for idx, count in enumerate(counts) if count > 0]
            if present_indices:
                present_energies = [energies[idx] for idx in present_indices]
                present_weights = [degeneracies[idx] for idx in present_indices]
                presence_mean = weighted_mean(present_energies, present_weights)
                weighted_configs = sum(present_weights)
            else:
                presence_mean = float("nan")
                weighted_configs = 0

            mean_count = weighted_mean(counts, degeneracies)
            slope, corr = weighted_linear_stats(counts, energies, degeneracies)
            cluster_rows.append(
                {
                    "cutoff_distance_ang": cutoff_distance,
                    "cluster_size": size,
                    "weighted_configs_with_cluster": weighted_configs,
                    "weighted_mean_cluster_count": mean_count,
                    "presence_delta_eV": presence_mean - weighted_global_mean if present_indices else float("nan"),
                    "slope_eV_per_cluster": slope,
                    "weighted_corr": corr,
                }
            )

    return cluster_rows


def write_gnuplot_script(
    script_path: Path,
    png_path: Path,
    pair_csv: Path,
    radial_curve_csv: Path | None,
    radial_csv: Path | None,
    cluster_csv: Path | None,
    level_label: str,
    substitutions: int,
) -> None:
    radial_curve_plot = ""
    enrichment_plot = ""
    radial_plot = ""
    cluster_slope_plot = ""
    cluster_presence_plot = ""

    if radial_curve_csv is not None:
        radial_curve_plot = (
            "set ytics nomirror\n"
            "unset y2tics\n"
            f"plot '{radial_curve_csv}' using 3:5 with lines lw 2 title 'Random', \\\n"
            f"     '{radial_curve_csv}' using 3:7 with lines lw 2 title 'All configs', \\\n"
            f"     '{radial_curve_csv}' using 3:9 with lines lw 2 title 'Low-energy'\n"
        )
        enrichment_plot = (
            f"plot '{radial_curve_csv}' using 3:10 with lines lw 2 title 'All / random', \\\n"
            f"     '{radial_curve_csv}' using 3:11 with lines lw 2 title 'Low-energy / random', \\\n"
            "     1.0 with lines dt 2 lw 1 title 'Ratio = 1'\n"
        )
    else:
        radial_curve_plot = "plot 1/0 notitle\n"
        enrichment_plot = "plot 1/0 notitle\n"

    if radial_csv is not None:
        radial_plot = (
            "set y2tics\n"
            "set ytics nomirror\n"
            f"plot '{radial_csv}' using 3:9 with linespoints lw 2 pt 7 title 'Slope per pair', \\\n"
            f"     '{radial_csv}' using 3:8 with linespoints lw 2 pt 5 title 'Presence delta', \\\n"
            f"     '{radial_csv}' using 3:6 axes x1y2 with impulses lw 2 title 'Pair fraction'\n"
        )
    else:
        radial_plot = "plot 1/0 notitle\n"

    if cluster_csv is not None:
        slope_lines = []
        presence_lines = []
        for size in range(2, max(2, substitutions) + 1):
            slope_lines.append(
                f"'{cluster_csv}' using 1:($2=={size} ? $6 : 1/0) with linespoints lw 2 pt 7 title 'size {size}'"
            )
            presence_lines.append(
                f"'{cluster_csv}' using 1:($2=={size} ? $5 : 1/0) with linespoints lw 2 pt 7 title 'size {size}'"
            )
        cluster_slope_plot = "plot " + ", \\\n     ".join(slope_lines) + "\n"
        cluster_presence_plot = "plot " + ", \\\n     ".join(presence_lines) + "\n"
    else:
        cluster_slope_plot = "plot 1/0 notitle\n"
        cluster_presence_plot = "plot 1/0 notitle\n"

    script_text = f"""set datafile separator comma
set terminal pngcairo size 1800,1200 enhanced font ',10'
set output '{png_path}'
set key outside
set grid
set tics out
set border lw 1.2
set multiplot layout 3,2 title 'OUTSOD / ENERGIES analysis: {level_label}'

set title 'Pair effect vs distance'
set xlabel 'Y-Y distance (A)'
set ylabel 'Delta E vs weighted mean (eV)'
set cblabel 'Weighted count'
set palette rgb 33,13,10
plot '{pair_csv}' using 9:8:4 with points pt 7 ps 1.3 lc palette title 'pairs'

set title 'Radial pair curves'
set xlabel 'Y-Y distance (A)'
set ylabel 'Pair fraction'
unset cblabel
{radial_curve_plot.rstrip()}

set title 'Radial enrichment'
set xlabel 'Y-Y distance (A)'
set ylabel 'Observed / random'
{enrichment_plot.rstrip()}

set title 'Radial-bin trends'
set xlabel 'Pair distance (A)'
set ylabel 'Energy trend (eV)'
unset cblabel
set y2label 'Weighted pair fraction'
{radial_plot.rstrip()}
unset y2tics
set ytics auto
unset y2label

set title 'Cluster slopes'
set xlabel 'Cluster cutoff distance (A)'
set ylabel 'Slope per cluster (eV)'
{cluster_slope_plot.rstrip()}

set title 'Cluster presence delta'
set xlabel 'Cluster cutoff distance (A)'
set ylabel 'Presence delta (eV)'
{cluster_presence_plot.rstrip()}

unset multiplot
"""
    script_path.write_text(script_text, encoding="utf-8")


def maybe_run_gnuplot(script_path: Path) -> tuple[bool, str]:
    gnuplot = shutil.which("gnuplot")
    if gnuplot is None:
        return False, "gnuplot not found in PATH"

    completed = subprocess.run(
        [gnuplot, str(script_path)],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip() or "gnuplot failed"
        return False, message
    return True, "ok"


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze how site pairs in OUTSOD correlate with ENERGIES, using "
            "Y-Y distances in Angstrom and distance-cutoff cluster summaries "
            "when INSOD/SGO or optimized GULP CIFs are available."
        )
    )
    parser.add_argument(
        "level_dir",
        help="Directory containing OUTSOD and ENERGIES, e.g. examples/case/n04",
    )
    parser.add_argument(
        "--structure-dir",
        help="Directory containing INSOD and SGO. Defaults to the parent of level_dir.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=12,
        help="Number of top stabilizing/destabilizing pairs to print (default: 12).",
    )
    parser.add_argument(
        "--shell-tol",
        type=float,
        default=1.0e-3,
        help="Tolerance in Angstrom used to group pair distances into shells (default: 1e-3).",
    )
    parser.add_argument(
        "--output",
        default="pair_energy_report.csv",
        help="Pair-level CSV filename to write inside the level directory (default: pair_energy_report.csv).",
    )
    parser.add_argument(
        "--shell-output",
        default="shell_energy_report.csv",
        help="Shell-level CSV filename to write inside the level directory (default: shell_energy_report.csv).",
    )
    parser.add_argument(
        "--radial-output",
        default="radial_energy_report.csv",
        help="Distance-bin CSV filename written inside the level directory (default: radial_energy_report.csv).",
    )
    parser.add_argument(
        "--radial-curve-output",
        default="radial_curve_report.csv",
        help="Radial pair-curve CSV filename written inside the level directory (default: radial_curve_report.csv).",
    )
    parser.add_argument(
        "--bin-width",
        type=float,
        default=0.1,
        help="Distance bin width in Angstrom for radial summaries (default: 0.1).",
    )
    parser.add_argument(
        "--cluster-output",
        default="cluster_energy_report.csv",
        help="Cluster-level CSV filename to write inside the level directory (default: cluster_energy_report.csv).",
    )
    parser.add_argument(
        "--cluster-shell-max",
        type=int,
        default=3,
        help="Maximum cumulative shell index used to build Y clusters (default: 3).",
    )
    parser.add_argument(
        "--cluster-max-size",
        type=int,
        default=4,
        help="Largest connected-cluster size to summarize (default: 4).",
    )
    parser.add_argument(
        "--cluster-cutoff-max-distance",
        type=float,
        default=4.0,
        help="Largest Y-Y cutoff distance in Angstrom used for cluster summaries (default: 4.0).",
    )
    parser.add_argument(
        "--y-symbol",
        default="",
        help="Optional substituted-species symbol used in optimized CIFs. If omitted, the script tries to infer it.",
    )
    parser.add_argument(
        "--low-energy-fraction",
        type=float,
        default=0.10,
        help=(
            "Weighted fraction of lowest-energy microstates used to build the low-energy radial curve "
            "(default: 0.10)."
        ),
    )
    parser.add_argument(
        "--plot-script",
        default="analysis_plots.gnuplot",
        help="gnuplot script filename written inside the level directory (default: analysis_plots.gnuplot).",
    )
    parser.add_argument(
        "--plot-png",
        default="analysis_plots.png",
        help="PNG filename written inside the level directory when gnuplot succeeds (default: analysis_plots.png).",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Do not generate or run the gnuplot summary script.",
    )
    parser.add_argument(
        "--no-geometry",
        action="store_true",
        help="Skip geometric analysis even if INSOD/SGO are available.",
    )
    args = parser.parse_args()

    level_dir = Path(args.level_dir).resolve()
    structure_dir = Path(args.structure_dir).resolve() if args.structure_dir else level_dir.parent.resolve()
    outsod_path = level_dir / "OUTSOD"
    energies_path = level_dir / "ENERGIES"

    if not outsod_path.is_file():
        raise SystemExit(f"Error: {outsod_path} not found")
    if not energies_path.is_file():
        raise SystemExit(f"Error: {energies_path} not found")

    substitutions, total_sites, degeneracies, subsets = parse_outsod(outsod_path)
    energies = parse_energies(energies_path)

    if len(subsets) != len(energies):
        raise SystemExit(
            f"Error: OUTSOD has {len(subsets)} configurations but ENERGIES has {len(energies)} entries"
        )

    global_mean = sum(energies) / float(len(energies))
    global_weighted_mean = weighted_mean(energies, degeneracies)

    pair_distances: dict[tuple[int, int], float] = {}
    pair_shells: dict[tuple[int, int], int] = {}
    shell_distances: dict[int, float] = {}
    shell_rows: list[dict[str, float]] = []
    radial_curve_rows: list[dict[str, float]] = []
    radial_rows: list[dict[str, float]] = []
    cluster_rows: list[dict[str, float]] = []
    geometry_warning = ""
    geometry_status = "skipped by user"
    optimized_found = 0
    optimized_used = 0
    low_energy_threshold = float("nan")

    if not args.no_geometry:
        insod_path = structure_dir / "INSOD"
        sgo_path = structure_dir / "SGO"
        if insod_path.is_file() and sgo_path.is_file():
            try:
                cell_params, site_coords = build_substitutable_site_geometry(structure_dir)
                if len(site_coords) != total_sites:
                    geometry_warning = (
                        f"geometry skipped: OUTSOD expects {total_sites} sites but INSOD/SGO reconstruct {len(site_coords)}"
                    )
                else:
                    pair_distances, pair_shells, shell_distances = build_pair_geometry(
                        site_coords, cell_params, args.shell_tol
                    )
                    optimized_pair_distances, optimized_found, optimized_used = load_optimized_pair_distances(
                        level_dir, subsets, degeneracies, args.y_symbol
                    )
                    if optimized_pair_distances:
                        pair_distances.update(optimized_pair_distances)
                        pair_shells, shell_distances = assign_shells_from_pair_distances(pair_distances, args.shell_tol)
                    if optimized_found > 0:
                        geometry_status = (
                            f"used optimized GULP CIFs for {optimized_used}/{optimized_found} available configurations"
                        )
                    else:
                        geometry_status = "used INSOD/SGO site geometry only"
                    shell_rows = summarize_shells(
                        subsets,
                        energies,
                        degeneracies,
                        pair_shells,
                        shell_distances,
                        global_weighted_mean,
                    )
                    cluster_rows = summarize_clusters(
                        subsets,
                        energies,
                        degeneracies,
                        pair_shells,
                        shell_distances,
                        global_weighted_mean,
                        max(1, args.cluster_shell_max),
                        max(2, args.cluster_max_size),
                    )
            except Exception as exc:  # noqa: BLE001
                geometry_warning = f"geometry skipped: {type(exc).__name__}: {exc}"
        else:
            geometry_warning = f"geometry skipped: INSOD/SGO not found in {structure_dir}"

    pair_rows = []
    for pair in combinations(range(1, total_sites + 1), 2):
        idxs = [idx for idx, subset in enumerate(subsets) if pair[0] in subset and pair[1] in subset]
        if not idxs:
            continue

        pair_energies = [energies[idx] for idx in idxs]
        pair_weights = [degeneracies[idx] for idx in idxs]
        pair_mean = sum(pair_energies) / float(len(pair_energies))
        pair_weighted_mean = weighted_mean(pair_energies, pair_weights)
        pair_row = {
            "site_i": pair[0],
            "site_j": pair[1],
            "count_unique": len(idxs),
            "count_weighted": sum(pair_weights),
            "mean_energy_eV": pair_mean,
            "weighted_mean_energy_eV": pair_weighted_mean,
            "delta_from_global_mean_eV": pair_mean - global_mean,
            "delta_from_weighted_global_mean_eV": pair_weighted_mean - global_weighted_mean,
            "distance_ang": pair_distances.get(pair, float("nan")),
            "shell_index": pair_shells.get(pair, -1),
            "shell_distance_ang": shell_distances.get(pair_shells.get(pair, -1), float("nan")),
        }
        pair_rows.append(pair_row)

    if pair_distances:
        radial_curve_rows, low_energy_threshold = summarize_radial_curves(
            subsets,
            energies,
            degeneracies,
            pair_distances,
            max(args.bin_width, 1.0e-6),
            args.low_energy_fraction,
        )
        radial_rows = summarize_radial_bins(
            subsets,
            energies,
            degeneracies,
            pair_distances,
            global_weighted_mean,
            max(args.bin_width, 1.0e-6),
        )
        cutoff_distances = [
            row["r_upper_ang"]
            for row in radial_rows
            if row["r_upper_ang"] <= max(args.cluster_cutoff_max_distance, 0.0) + 1.0e-12
        ]
        cluster_rows = summarize_clusters_distance(
            subsets,
            energies,
            degeneracies,
            pair_distances,
            cutoff_distances,
            global_weighted_mean,
            max(2, args.cluster_max_size),
        )

    pair_rows.sort(key=lambda row: row["delta_from_weighted_global_mean_eV"])

    out_path = level_dir / args.output
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "site_i",
                "site_j",
                "count_unique",
                "count_weighted",
                "mean_energy_eV",
                "weighted_mean_energy_eV",
                "delta_from_global_mean_eV",
                "delta_from_weighted_global_mean_eV",
                "distance_ang",
                "shell_index",
                "shell_distance_ang",
            ],
        )
        writer.writeheader()
        writer.writerows(pair_rows)

    shell_out_path = level_dir / args.shell_output
    if shell_rows:
        with shell_out_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "shell_index",
                    "shell_distance_ang",
                    "pairs_in_shell",
                    "weighted_configs_with_shell",
                    "presence_delta_eV",
                    "slope_eV_per_pair",
                    "weighted_corr",
                ],
            )
            writer.writeheader()
            writer.writerows(shell_rows)

    radial_out_path = level_dir / args.radial_output
    if radial_rows:
        with radial_out_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "r_lower_ang",
                    "r_upper_ang",
                    "r_center_ang",
                    "pairs_in_bin",
                    "weighted_pair_occurrences",
                    "weighted_pair_fraction",
                    "weighted_configs_with_bin",
                    "presence_delta_eV",
                    "slope_eV_per_pair",
                    "weighted_corr",
                ],
            )
            writer.writeheader()
            writer.writerows(radial_rows)

    radial_curve_out_path = level_dir / args.radial_curve_output
    if radial_curve_rows:
        with radial_curve_out_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "r_lower_ang",
                    "r_upper_ang",
                    "r_center_ang",
                    "n_lattice_pairs",
                    "random_pair_fraction",
                    "all_weighted_pair_occurrences",
                    "all_pair_fraction",
                    "low_weighted_pair_occurrences",
                    "low_pair_fraction",
                    "all_over_random",
                    "low_over_random",
                ],
            )
            writer.writeheader()
            writer.writerows(radial_curve_rows)

    cluster_out_path = level_dir / args.cluster_output
    if cluster_rows:
        with cluster_out_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "cutoff_distance_ang",
                    "cluster_size",
                    "weighted_configs_with_cluster",
                    "weighted_mean_cluster_count",
                    "presence_delta_eV",
                    "slope_eV_per_cluster",
                    "weighted_corr",
                ],
            )
            writer.writeheader()
            writer.writerows(cluster_rows)

    plot_script_path = level_dir / args.plot_script
    plot_png_path = level_dir / args.plot_png
    plot_status = "skipped by user"
    if not args.no_plots and radial_rows:
        write_gnuplot_script(
            plot_script_path,
            plot_png_path,
            out_path,
            radial_curve_out_path if radial_curve_rows else None,
            radial_out_path if radial_rows else None,
            cluster_out_path if cluster_rows else None,
            str(level_dir),
            substitutions,
        )
        ran_ok, message = maybe_run_gnuplot(plot_script_path)
        if ran_ok:
            plot_status = f"generated {plot_png_path.name}"
        else:
            plot_status = f"script only ({message})"

    print(f"Level directory           : {level_dir}")
    print(f"Structure directory       : {structure_dir}")
    print(f"Configurations            : {len(energies)}")
    print(f"Substitutions             : {substitutions}")
    print(f"Total substitutable sites : {total_sites}")
    print(f"Global mean energy (eV)   : {global_mean:.8f}")
    print(f"Weighted global mean (eV) : {global_weighted_mean:.8f}")
    print(f"Pair CSV report           : {out_path}")
    if radial_rows:
        if radial_curve_rows:
            print(f"Radial curve CSV report   : {radial_curve_out_path}")
        print(f"Radial CSV report         : {radial_out_path}")
        if cluster_rows:
            print(f"Cluster CSV report        : {cluster_out_path}")
        print(f"Plot script               : {plot_script_path}")
        print(f"Plot status               : {plot_status}")
        print(f"Geometry basis            : {geometry_status}")
        if not math.isnan(low_energy_threshold):
            print(
                f"Low-energy radial window  : lowest {100.0 * max(0.0, min(1.0, args.low_energy_fraction)):.1f}% "
                f"(threshold <= {low_energy_threshold:.8f} eV)"
            )
    elif geometry_warning:
        print(f"Geometry status           : {geometry_warning}")
    else:
        print("Geometry status           : skipped by user")

    print()
    print("Most stabilizing pairs (weighted delta, eV):")
    for row in pair_rows[: args.top]:
        extra = ""
        if not math.isnan(row["distance_ang"]):
            extra = f"  d={row['distance_ang']:.4f} A"
        print(
            f"  ({row['site_i']:>2d},{row['site_j']:>2d})  "
            f"delta={row['delta_from_weighted_global_mean_eV']:+.6f}  "
            f"weighted_count={row['count_weighted']}{extra}"
        )

    print()
    print("Most destabilizing pairs (weighted delta, eV):")
    for row in pair_rows[-args.top :]:
        extra = ""
        if not math.isnan(row["distance_ang"]):
            extra = f"  d={row['distance_ang']:.4f} A"
        print(
            f"  ({row['site_i']:>2d},{row['site_j']:>2d})  "
            f"delta={row['delta_from_weighted_global_mean_eV']:+.6f}  "
            f"weighted_count={row['count_weighted']}{extra}"
        )

    if radial_rows:
        print()
        print("Radial-bin trends (weighted slope in eV per extra Y-Y pair):")
        for row in radial_rows[: min(8, len(radial_rows))]:
            print(
                f"  {row['r_lower_ang']:.2f}-{row['r_upper_ang']:.2f} A  "
                f"slope={row['slope_eV_per_pair']:+.6f}  corr={row['weighted_corr']:+.3f}  "
                f"presence_delta={row['presence_delta_eV']:+.6f}  "
                f"pair_fraction={row['weighted_pair_fraction']:.4f}"
            )

    if radial_curve_rows:
        print()
        print("Radial enrichment (ratio vs random):")
        shown = 0
        for row in radial_curve_rows:
            if row["n_lattice_pairs"] <= 0:
                continue
            print(
                f"  {row['r_lower_ang']:.2f}-{row['r_upper_ang']:.2f} A  "
                f"all/random={row['all_over_random']:.3f}  "
                f"low/random={row['low_over_random']:.3f}"
            )
            shown += 1
            if shown >= 8:
                break

    if cluster_rows:
        print()
        print("Cluster trends (connected Y clusters below a distance cutoff):")
        for row in cluster_rows:
            if row["cluster_size"] > substitutions:
                continue
            if row["weighted_configs_with_cluster"] <= 0:
                continue
            print(
                f"  r_cut<={row['cutoff_distance_ang']:.2f} A  "
                f"size={int(row['cluster_size'])}  "
                f"slope={row['slope_eV_per_cluster']:+.6f}  "
                f"corr={row['weighted_corr']:+.3f}  "
                f"presence_delta={row['presence_delta_eV']:+.6f}"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
