#!/usr/bin/env python3
"""
Comprehensive Comparison of Two MD Trajectories
===============================================

This script performs an extensive structural, dynamic, and energetic comparison between
two molecular dynamics (MD) trajectories — typically a "Complete Structure" vs a
"Published Structure" — to evaluate relative stability, convergence, and predictive utility.

Analyses Performed
------------------
1. Full-protein RMSF (all residues)                       → *_rmsf_full.pdf
2. Aligned-residue RMSF (based on sequence alignment)     → *_rmsf_aligned.pdf
3. Radius of Gyration (Rg) over time                      → *_rg.pdf
4. Backbone RMSD vs. time                                 → *_rmsd_time.pdf
5. Cα-based k-means cluster population bar chart          → *_cluster_populations.pdf
6. PCA scatter plot (PC1 vs PC2, aligned Cα)              → *_pca_scatter.pdf
7. PCA statistics: variance explained & standard deviation → printed to screen
8. Persistent H-bonds (≥30% occupancy)                    → *_hbonds_persistent.csv
   - Summary of shared and exclusive H-bonds printed to screen
   - Top 20 persistent H-bonds printed to screen
9. DSSP secondary structure analysis                      → *_dssp.pdf
10. Per-residue DSSP difference heatmap                   → *_dssp_difference.pdf
11. DCCM heatmaps (Cα-only) for each structure            → *_dccm_ca_LABEL.pdf
12. Differential DCCM (Δ-correlation matrix)              → *_dccm_delta.pdf
13. Ligand–loop DCCM correlation summary                  → *_dccm_ligand_loop_residue.pdf
14. MD output/screenout validation (if provided)            → printed / PDFs / JSON
    - Auto-detects Python/OpenMM, DivCon/C++/OpenMM, or legacy three-column output
    - Compares system construction, energy decomposition, protocol, and thermodynamics
    - If OpenMM XmlSerializer files are beside the MD logs, compares the serialized
      System, Integrator, and initial State semantically
15. Ligand Analyses (if `--ligand` is specified):
    a. Ligand RMSD vs. time                               → *_ligand_rmsd.pdf
    b. Ligand SASA vs. time                               → *_ligand_sasa.pdf
    c. Ligand per-atom RMSF                               → *_ligand_rmsf.pdf
    d. Ligand-pocket H-bond persistence                   → printed to screen
    e. Ligand-pocket contact fingerprint matrix           → printed to screen
    f. Ligand–loop DCCM statistics                        → printed to screen
16. Summary JSON file containing all key metrics          → *_summary.json
17. Trajectory-validation additions (non-destructive):
    a. Aligned RMSF profile similarity + ΔRMSF              → *_rmsf_delta.pdf
    b. Common-reference aligned-Cα RMSD                     → *_rmsd_common_reference.pdf
    c. Rg/RMSD first-vs-last-window drift statistics        → printed to screen / JSON
    d. Mean aligned-Cα internal-distance comparison         → *_ca_distance_delta.pdf

JSON Output: *_summary.json
---------------------------
The JSON file contains all major results in structured format, supporting downstream parsing by
automated scoring tools, AI assistants, or cross-comparative scripts.

Top-level fields include:

  - label1 / label2: human-readable names for the input trajectories
  - out_prefix: basename used for all generated files
  - generated_on: ISO timestamp of run
  - units: dictionary of measurement units for each metric
  - n_frames / n_atoms / n_residues: structure metadata
  - rmsf: RMSF statistics (full and aligned)
  - rg: radius of gyration statistics + KS p-value
  - rmsd_backbone: RMSD (mean, std, max) vs time
  - trajectory_validation: paired trajectory-comparison metrics (RMSF profile similarity,
    common-reference RMSD, drift/stability, and aligned-Cα internal distances)
  - clustering: number of clusters + frame counts
  - pca: PC1/PC2 explained variance, stddev, centroid distance
  - hbonds: shared and exclusive persistent hydrogen bonds
  - top_hbonds: top persistent H-bond interactions
  - dssp_diff: residue-wise secondary structure disagreements (>50% occupancy)
  - dccm_overall: DCCM correlations and Δ-correlation highlights
  - ligand: (if applicable) ligand RMSD, RMSF, SASA, H-bonds, contact fingerprint, and DCCM
  - out1_summary / out2_summary: backwards-compatible thermodynamic summaries
  - openmm_validation: normalized MD-log data plus system, energy, protocol, production, and
    optional serialized OpenMM System/Integrator/State comparisons

Usage Examples
--------------
Basic:
    python analysisMD.py traj1.dcd.gz top1.prmtop.gz traj2.dcd.gz top2.prmtop.gz

With MD screenout/log files (plain text or gzip-compressed):
    python analysisMD.py traj1.dcd.gz top1.prmtop.gz traj2.dcd.gz top2.prmtop.gz \
        --out1 run1/md.screenout.gz --out2 run2/md.screenout

With ligand analysis (resname e.g. CBN):
    python analysisMD.py traj1.dcd.gz top1.prmtop.gz traj2.dcd.gz top2.prmtop.gz \
        --ligand CBN

Custom labeling and output prefix:
    python analysisMD.py traj1.dcd.gz top1.prmtop.gz traj2.dcd.gz top2.prmtop.gz \
        --label1 "Complete" --label2 "Published" -o 1kzn

Dependencies
------------
  - mdtraj
  - numpy
  - matplotlib
  - seaborn
  - Biopython (Bio.Align)
  - scikit-learn (PCA and clustering)
  - gzip
  - DSSP (external binary; install via conda: `conda install -c salilab dssp`)

Developer Notes
---------------
- Temporary files (from gzipped input) are auto-cleaned on exit.
- Sequence alignment uses Biopython’s global aligner.
- All plots are DPI-configurable via CONFIG.
- JSON summary is machine- and AI-readable for downstream scoring models.
- Ligand analysis requires ligand to be present and properly identified in both topologies.

Interpretation Notes
---------------------
- High RMSF/RMSD and poor DCCM agreement suggest unstable or misprepared structures.
- Tight PCA clustering and low radius of gyration imply convergence and compactness.
- Shared hydrogen bonds and consistent DSSP assignments are positive indicators.
- Ligand metrics (RMSD, SASA, DCCM, H-bonding) help assess binding pose stability and consistency.
- Always examine ΔDCCM and fingerprint matrices to pinpoint structure-specific behaviors.

This script is suitable for pre-screening MD simulations before ensemble averaging, AI-based prediction, or free energy calculation workflows.
"""

from __future__ import annotations

# ─── Standard Library ────────────────────────────────────────────────────────
import argparse
import gzip
import heapq
import hashlib
import shutil
import csv
import logging
import re
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple
import json
from datetime import datetime, timezone

# ─── Third-Party Libraries ───────────────────────────────────────────────────
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ks_2samp, pearsonr, spearmanr
from Bio.Align import PairwiseAligner
from Bio.Data import IUPACData
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

try:
    import openmm as mm
    from openmm import unit as openmm_unit
except ImportError:
    mm = None
    openmm_unit = None

# ─────────────────────────────  CONFIG  ──────────────────────────────────────
residue_mapping = {             # standard 3‑letter + Amber variants
    **{aa: aa for aa in
       ("ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS",
        "MET","PHE","PRO","SER","THR","TRP","TYR","VAL")},
    "HIE":"HIS","HID":"HIS","HIP":"HIS","CYX":"CYS","GLH":"GLU","ASH":"ASP",
}
aa3to1 = {k.upper(): v for k, v in IUPACData.protein_letters_3to1_extended.items()}
_trim3 = re.compile(r"([A-Z]{3})").match

CONFIG = {
    "rmsf_units": "angstrom",        # units to report RMSF (typically Å)
    "rmsf_scale": 10.0,              # scale factor from nm to Å

    "min_frames_for_clustering": 10,
    "max_reference_frames": 50,
    "max_clusters": 10,
    "default_k_clusters": 5,

    "top_n_persistent_hbonds": 20,
    "hbond_threshold": 0.30,         # % persistence threshold

    "dccm": {
        "top_n": 5,                  # number of top DCCM pairs (positive/negative)
        "ligand_loop_top_n": 5       # number of top ligand-loop residues to print
    },

    "contact_fingerprint": {
        "distance_cutoff_nm": 0.35,  # 3.5 Å cutoff for defining contact
        "min_occupancy": 0.30        # minimum contact occupancy to report
    },

    # Additional trajectory-only validation.  These metrics compare structural
    # behavior derived from the trajectories themselves; they do not depend on
    # MD output/screen output and therefore remain useful when the two simulations
    # were run with different hardware, random seeds, or solvent realizations.
    "trajectory_validation": {
        "stability_fraction": 0.25,        # compare first/last 25% of each trajectory
        "top_n_rmsf_differences": 10,
        "top_n_distance_differences": 10,
        "distance_pair_chunk_size": 5000   # bound temporary memory for large proteins
    },

    # Semantic comparison of OpenMM XmlSerializer output.  These tolerances are
    # intentionally tight: the XML files describe the actual OpenMM objects, not
    # noisy trajectory observables.  Larger differences should be surfaced for
    # inspection rather than silently accepted.
    "openmm_xml_validation": {
        "absolute_tolerance": 1.0e-10,
        "relative_tolerance": 1.0e-10,
        "position_tolerance_nm": 1.0e-7,
        "energy_tolerance_kj_per_mol": 1.0e-6,
        "constraint_geometry_thresholds_nm": (1.0e-6, 1.0e-5, 1.0e-4),
        "top_n_differences": 10
    },

    "plot": {
        "dpi": 300,
        "figsize": {
            "dccm": (8, 6),
            "ligand_loop": (5, 4),
            "dssp": (10, 6),
            "dssp_diff": (12, 6)
        }
    }
}

# ──────────────────────────  HELPERS  ────────────────────────────────────────
def decompress_if_gz(path: Path, tmp_files: list[Path]) -> Path:
    if path.suffix != ".gz":
        return path

    # Extract extension from e.g. output.dcd.gz → .dcd
    inner_ext = "".join(Path(path.stem).suffixes) or ".dat"
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=inner_ext, prefix="mdtraj_")
    with gzip.open(path, "rb") as gz_in, open(tmp.name, "wb") as out:
        shutil.copyfileobj(gz_in, out)

    tmp_path = Path(tmp.name)
    logging.info("Decompressing %s → %s", path, tmp_path)
    tmp_files.append(tmp_path)
    return tmp_path

def trim3(name:str)->str:
    m=_trim3(name.upper()); return m.group(1) if m else name[:3].upper()

def standardise_names(traj: md.Trajectory)->None:
    for res in traj.topology.residues:
        core = residue_mapping.get(trim3(res.name), trim3(res.name))
        res.name = core

def load_traj(traj: Path, top: Path) -> md.Trajectory:
    t = md.load(traj, top=top)
    n_residues = t.topology.n_residues
    logging.info("Loaded %-20s  frames:%5d  atoms:%6d  residues:%4d",
                 traj.name, t.n_frames, t.n_atoms, n_residues)
    return t

def seq_and_ca(traj: md.Trajectory)->Tuple[str,List[int]]:
    seq, ca = [], []
    for res in traj.topology.residues:
        if not res.is_protein: continue
        letter = aa3to1.get(res.name);  # skip unknown AA
        if not letter: continue
        try:
            ca_atom = next(a for a in res.atoms if a.name=="CA")
        except StopIteration:
            continue
        seq.append(letter); ca.append(ca_atom.index)
    return "".join(seq), ca

def compute_rmsf(traj: md.Trajectory, idx: np.ndarray)->np.ndarray:
    # md.rmsf() centers the target trajectory in place.  Always operate on a
    # private copy so one analysis cannot silently alter coordinates used by
    # later analyses (PCA, DCCM, clustering, etc.).
    work = traj[:]  # Trajectory slicing copies coordinate arrays by default
    return md.rmsf(work, reference=work, atom_indices=idx)*CONFIG["rmsf_scale"]  # nm→Å

def aligned_trajectory_copy(target: md.Trajectory, reference: md.Trajectory,
                            target_indices: np.ndarray, ref_indices: np.ndarray,
                            frame: int = 0) -> md.Trajectory:
    """Return a copy of *target* rigid-body aligned to one common reference frame."""
    aligned = target[:]  # never mutate the caller's trajectory
    aligned.superpose(
        reference,
        frame=frame,
        atom_indices=np.asarray(target_indices, dtype=int),
        ref_atom_indices=np.asarray(ref_indices, dtype=int)
    )
    return aligned

def paired_profile_stats(values1: np.ndarray, values2: np.ndarray) -> dict:
    """Return paired similarity metrics for two equally sized structural profiles."""
    a = np.asarray(values1, dtype=float)
    b = np.asarray(values2, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"Paired profiles have different shapes: {a.shape} vs {b.shape}")

    finite = np.isfinite(a) & np.isfinite(b)
    a = a[finite]
    b = b[finite]
    if a.size == 0:
        return {
            "n": 0, "pearson_r": None, "spearman_rho": None,
            "rmse": None, "mae": None, "max_abs_difference": None
        }

    diff = a - b

    # Pearson/Spearman are undefined for constant profiles.  Preserve the other
    # distance metrics and represent undefined correlations as None in JSON.
    pearson = None
    spearman = None
    if a.size > 1 and np.std(a) > 0.0 and np.std(b) > 0.0:
        pearson = float(pearsonr(a, b)[0])
        spearman = float(spearmanr(a, b)[0])

    return {
        "n": int(a.size),
        "pearson_r": pearson,
        "spearman_rho": spearman,
        "rmse": float(np.sqrt(np.mean(diff * diff))),
        "mae": float(np.mean(np.abs(diff))),
        "max_abs_difference": float(np.max(np.abs(diff)))
    }

def trajectory_stability_stats(values: np.ndarray, fraction: float = 0.25) -> dict:
    """Summarize slow drift without imposing a pass/fail convergence criterion."""
    y = np.asarray(values, dtype=float)
    y = y[np.isfinite(y)]
    if y.size == 0:
        return {
            "n": 0, "window_fraction": float(fraction),
            "first_window_mean": None, "last_window_mean": None,
            "last_minus_first": None, "slope_per_frame": None,
            "projected_change_over_trajectory": None
        }

    fraction = min(max(float(fraction), 0.0), 0.5)
    n_window = max(1, int(np.ceil(y.size * fraction)))
    x = np.arange(y.size, dtype=float)
    slope = float(np.polyfit(x, y, 1)[0]) if y.size > 1 else 0.0
    first_mean = float(np.mean(y[:n_window]))
    last_mean = float(np.mean(y[-n_window:]))

    return {
        "n": int(y.size),
        "window_fraction": float(fraction),
        "first_window_mean": first_mean,
        "last_window_mean": last_mean,
        "last_minus_first": float(last_mean - first_mean),
        "slope_per_frame": slope,
        "projected_change_over_trajectory": float(slope * max(y.size - 1, 0))
    }

def mean_pair_distances(traj: md.Trajectory, atom_pairs: np.ndarray, chunk_size: int = 5000) -> np.ndarray:
    """Compute mean pair distances in Å while bounding temporary memory usage."""
    pairs = np.asarray(atom_pairs, dtype=int)
    means = np.empty(len(pairs), dtype=float)
    chunk_size = max(1, int(chunk_size))
    for start in range(0, len(pairs), chunk_size):
        stop = min(start + chunk_size, len(pairs))
        distances_nm = md.compute_distances(traj, pairs[start:stop], periodic=False)
        means[start:stop] = distances_nm.mean(axis=0) * 10.0
    return means

def safe(label:str)->str:  # filename‑safe
    return re.sub(r"[^A-Za-z0-9]+","_",label.strip()).lower()

def res_label(res_or_atom):
    """Return residue label 'A:ARG42' for either an Atom or Residue."""
    res = res_or_atom.residue if hasattr(res_or_atom, "residue") else res_or_atom
    chain = chr(65 + res.chain.index)
    return f"{chain}:{res.name}{res.resSeq}"
    
_FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


def open_text_auto(path: Path):
    """Open a plain-text or gzip-compressed log by inspecting its file signature."""
    path = Path(path)
    with open(path, "rb") as probe:
        is_gzip = probe.read(2) == b"\x1f\x8b"
    if is_gzip:
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")


def _finite_float(value):
    """Convert a numeric token to float, representing NaN/Inf as None for valid JSON."""
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if np.isfinite(parsed) else None


def _parse_bool(value: str):
    value = value.strip().lower()
    if value == "true":
        return True
    if value == "false":
        return False
    return None


def _first_match(text: str, pattern: str, cast=float, flags=0):
    match = re.search(pattern, text, flags)
    if not match:
        return None
    try:
        return cast(match.group(1))
    except (TypeError, ValueError):
        return None


def _first_box(text: str, pattern: str):
    match = re.search(pattern, text)
    if not match:
        return None
    values = re.findall(_FLOAT_RE, match.group(1))
    if len(values) < 3:
        return None
    return [float(v) for v in values[:3]]


def _parse_energy_decompositions(lines: list[str]) -> dict:
    sections = {}
    key_map = {
        "Bond": "bond",
        "Angle": "angle",
        "Torsion": "torsion",
        "Nonbonded": "nonbonded",
        "4-part sum": "four_part_sum",
        "Total": "total",
        "Residual": "residual"
    }

    for i, line in enumerate(lines):
        header = re.search(
            r"OpenMM energy decomposition - (initial|post-minimization) \(kJ/mol\):",
            line,
            re.IGNORECASE
        )
        if not header:
            continue

        section_name = header.group(1).lower().replace("-", "_")
        values = {}
        for candidate in lines[i + 1:i + 10]:
            cleaned = re.sub(r"^\s*\[DEBUG\]\s*", "", candidate).strip()
            match = re.match(
                rf"(Bond|Angle|Torsion|Nonbonded|4-part sum|Total|Residual)\s*:\s*({_FLOAT_RE})",
                cleaned
            )
            if match:
                values[key_map[match.group(1)]] = float(match.group(2))
        if values:
            sections[section_name] = values

    return sections


def _parse_representative_particles(text: str) -> dict:
    particles = {}
    pattern = re.compile(
        rf"index=\s*(\d+)\s+(\S+)\s+DivConType=(\S+)\s+"
        rf"DivConQ=({_FLOAT_RE})\s+OpenMMQ=({_FLOAT_RE})\s+"
        rf"sigma=({_FLOAT_RE})\s+nm\s+epsilon=({_FLOAT_RE})\s+kJ/mol"
    )

    for match in pattern.finditer(text):
        index, identity, divcon_type, divcon_q, openmm_q, sigma, epsilon = match.groups()
        identity_upper = identity.upper()
        type_upper = divcon_type.upper()

        if identity_upper in {"HOH:O", "WAT:O"}:
            role = "water_o"
        elif identity_upper in {"HOH:H1", "WAT:H1"}:
            role = "water_h1"
        elif identity_upper in {"HOH:H2", "WAT:H2"}:
            role = "water_h2"
        elif identity_upper.startswith("NA:") or type_upper == "NA+":
            role = "sodium"
        elif identity_upper.startswith("CL:") or type_upper == "CL-":
            role = "chloride"
        else:
            continue

        particles[role] = {
            "index": int(index),
            "identity": identity,
            "divcon_type": divcon_type,
            "divcon_charge_e": float(divcon_q),
            "openmm_charge_e": float(openmm_q),
            "sigma_nm": float(sigma),
            "epsilon_kj_per_mol": float(epsilon)
        }

    return particles


def _parse_water_constraints(lines: list[str]) -> list[dict]:
    constraints = []
    in_section = False
    for line in lines:
        if "Representative water constraints:" in line:
            in_section = True
            continue
        if not in_section:
            continue

        cleaned = re.sub(r"^\s*\[DEBUG\]\s*", "", line).strip()
        match = re.match(rf"(\d+)-(\d+)\s*:\s*({_FLOAT_RE})\s+nm$", cleaned)
        if match:
            constraints.append({
                "particle1": int(match.group(1)),
                "particle2": int(match.group(2)),
                "distance_nm": float(match.group(3))
            })
            continue

        if constraints:
            break

    return constraints


def _parse_python_production(lines: list[str]) -> list[dict]:
    rows = []
    in_table = False

    for line in lines:
        if line.startswith("Production metrics (step,"):
            in_table = True
            continue
        if not in_table:
            continue
        if rows and (line.startswith("Elapsed time:") or line.startswith("Simulation summary:")):
            break

        parts = [part.strip() for part in line.strip().split(",")]
        if len(parts) != 11 or not parts[0].isdigit():
            continue

        rows.append({
            "step": int(parts[0]),
            "time_ps": None,
            "potential_energy_kj_per_mol": _finite_float(parts[1]),
            "kinetic_energy_kj_per_mol": _finite_float(parts[2]),
            "temperature_k": _finite_float(parts[3]),
            "volume_nm3": _finite_float(parts[4]),
            "density_g_per_ml": _finite_float(parts[5]),
            "target_pressure_bar": _finite_float(parts[6]),
            "instantaneous_pressure_bar": _finite_float(parts[7]),
            "wall_time_s": _finite_float(parts[10])
        })

    return rows


def _parse_cpp_production(lines: list[str]) -> list[dict]:
    rows = []
    in_table = False

    for line in lines:
        if "Steps" in line and "PE(kcal/mol)" in line and "Vol(A^3)" in line:
            in_table = True
            continue
        if not in_table:
            continue
        if "OpenMM Molecular Dynamics Summary" in line:
            break

        parts = [part.strip() for part in line.split("|")]
        if len(parts) < 9 or not parts[0].isdigit():
            continue

        pe_kcal = _finite_float(parts[2])
        ke_kcal = _finite_float(parts[3])
        volume_a3 = _finite_float(parts[5])
        rows.append({
            "step": int(parts[0]),
            "time_ps": _finite_float(parts[1]),
            "potential_energy_kj_per_mol": pe_kcal * 4.184 if pe_kcal is not None else None,
            "kinetic_energy_kj_per_mol": ke_kcal * 4.184 if ke_kcal is not None else None,
            "temperature_k": _finite_float(parts[4]),
            "volume_nm3": volume_a3 / 1000.0 if volume_a3 is not None else None,
            "density_g_per_ml": _finite_float(parts[6]),
            "target_pressure_bar": None,
            "instantaneous_pressure_bar": None,
            "wall_time_s": None
        })

    return rows


def _parse_legacy_production(lines: list[str]) -> list[dict]:
    """Retain support for the older three-column step,energy,temp OUT format."""
    rows = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = [part.strip() for part in stripped.split(",")]
        if len(parts) != 3 or not parts[0].isdigit():
            continue

        energy = _finite_float(parts[1])
        temperature = _finite_float(parts[2])
        if energy is None or temperature is None:
            continue
        rows.append({
            "step": int(parts[0]),
            "time_ps": None,
            "potential_energy_kj_per_mol": energy,
            "kinetic_energy_kj_per_mol": None,
            "temperature_k": temperature,
            "volume_nm3": None,
            "density_g_per_ml": None,
            "target_pressure_bar": None,
            "instantaneous_pressure_bar": None,
            "wall_time_s": None
        })
    return rows


def _parse_wall_time_hms(value: str):
    """Convert HH:MM:SS wall-time text to seconds."""
    match = re.fullmatch(r"\s*(\d+):(\d+):(\d+(?:\.\d+)?)\s*", value)
    if not match:
        return None
    hours, minutes, seconds = match.groups()
    return 3600.0 * int(hours) + 60.0 * int(minutes) + float(seconds)


def _finalize_timing(timing: dict) -> None:
    """Add derived timing fields without inventing missing stage timings."""
    stage_fields = (
        "minimization_wall_time_s",
        "nvt_equilibration_wall_time_s",
        "npt_equilibration_wall_time_s",
        "production_wall_time_s",
    )
    stage_values = [timing.get(field) for field in stage_fields]
    if all(value is not None for value in stage_values):
        timing["timed_md_stages_wall_time_s"] = float(sum(stage_values))

    production = timing.get("production_wall_time_s")
    total = timing.get("total_execution_wall_time_s")
    if production is not None and total is not None and total > 0.0:
        timing["non_production_execution_wall_time_s"] = float(total - production)
        timing["production_fraction_of_total_percent"] = float(100.0 * production / total)

    timed_md = timing.get("timed_md_stages_wall_time_s")
    if timed_md is not None and total is not None and total > 0.0:
        timing["non_timed_md_execution_wall_time_s"] = float(total - timed_md)
        timing["timed_md_fraction_of_total_percent"] = float(100.0 * timed_md / total)


def _parse_python_timing(text: str) -> dict:
    """Parse wall-clock timings emitted by runMD-fromDivCon.py."""
    timing = {}
    stage_patterns = {
        "minimization_wall_time_s": (
            rf"Minimizing with OpenMM:.*?Post-minimization Potential Energy:.*?"
            rf"Elapsed time:\s*({_FLOAT_RE})\s+seconds"
        ),
        "nvt_equilibration_wall_time_s": (
            rf"Equilibrating \(NVT\).*?Elapsed time:\s*({_FLOAT_RE})\s+seconds"
        ),
        "npt_equilibration_wall_time_s": (
            rf"Equilibrating \(NPT\).*?Elapsed time:\s*({_FLOAT_RE})\s+seconds"
        ),
        "production_wall_time_s": (
            rf"Running Production NPT Simulation.*?Elapsed time:\s*({_FLOAT_RE})\s+seconds"
        ),
    }
    for field, pattern in stage_patterns.items():
        value = _first_match(text, pattern, float, re.DOTALL)
        if value is not None:
            timing[field] = value

    # Prefer the explicit production summary if present.  It is the same timer
    # used for the production table and is not the total Python process runtime.
    summary_wall_time = _first_match(text, rf"^\s*wall_time_s:\s*({_FLOAT_RE})", float, re.MULTILINE)
    if summary_wall_time is not None:
        timing["production_wall_time_s"] = summary_wall_time

    # Current Python logs do not emit a whole-process timer, but accept one if
    # future versions add it.  Do not substitute the sum of stage timers: that
    # would omit setup, force-field construction, and output work.
    total = _first_match(
        text,
        rf"(?:Total Computation Time \(Seconds\):|Total execution time(?: \(seconds\))?:)\s*({_FLOAT_RE})",
        float,
        re.IGNORECASE,
    )
    if total is not None:
        timing["total_execution_wall_time_s"] = total

    _finalize_timing(timing)
    return timing


def _parse_cpp_timing(text: str) -> dict:
    """Parse wall-clock timings emitted by the DivCon/OpenMM MD path."""
    timing = {}

    minimization = _first_match(
        text, rf"\[MD\]\s+Minimization time:\s*({_FLOAT_RE})\s+seconds", float
    )
    if minimization is not None:
        timing["minimization_wall_time_s"] = minimization

    nvt_time = _first_match(
        text, rf"\[MD\]\s+NVT equilibration time:\s*({_FLOAT_RE})\s+seconds", float
    )
    if nvt_time is not None:
        timing["nvt_equilibration_wall_time_s"] = nvt_time

    npt_time = _first_match(
        text, rf"\[MD\]\s+NPT equilibration time:\s*({_FLOAT_RE})\s+seconds", float
    )
    if npt_time is None:
        npt_time = _first_match(text, rf"NPT time:\s*({_FLOAT_RE})", float)
    if npt_time is not None:
        timing["npt_equilibration_wall_time_s"] = npt_time

    production = _first_match(
        text, rf"Production total time:\s*({_FLOAT_RE})\s+seconds", float
    )
    if production is None:
        formatted = _first_match(
            text, r"Elapsed wall time\s*:\s*(\d+:\d+:\d+(?:\.\d+)?)", str
        )
        if formatted is not None:
            production = _parse_wall_time_hms(formatted)
    if production is not None:
        timing["production_wall_time_s"] = production

    timed_stages = _first_match(
        text, rf"\[MD\]\s+Timed MD stages total:\s*({_FLOAT_RE})\s+seconds", float
    )
    if timed_stages is not None:
        timing["timed_md_stages_wall_time_s"] = timed_stages

    total = _first_match(
        text, rf"Total Computation Time \(Seconds\):\s*({_FLOAT_RE})", float
    )
    if total is not None:
        timing["total_execution_wall_time_s"] = total

    _finalize_timing(timing)
    return timing


def _infer_production_timing(protocol: dict, production: list[dict]) -> None:
    steps = sorted({row["step"] for row in production})
    positive_diffs = [b - a for a, b in zip(steps, steps[1:]) if b > a]
    if positive_diffs and len(set(positive_diffs)) == 1 and "production_report_interval_steps" not in protocol:
        protocol["production_report_interval_steps"] = positive_diffs[0]
        protocol["production_report_interval_source"] = "inferred_from_production_rows"

    duration = protocol.get("production_duration_ps")
    nsteps = protocol.get("production_steps")
    if duration is not None and nsteps:
        protocol["timestep_ps"] = float(duration) / float(nsteps)
        protocol["timestep_source"] = "inferred_from_reported_production_duration"
        for row in production:
            if row.get("time_ps") is None:
                row["time_ps"] = float(row["step"]) * protocol["timestep_ps"]


def parse_out_file(out_path: Path):
    """
    Parse a Python/OpenMM, DivCon/C++/OpenMM, or legacy MD output file.

    Input may be plain text or gzip-compressed. Returned values use a normalized
    representation so output 1 and output 2 can be compared independent of which
    implementation produced them.
    """
    try:
        with open_text_auto(out_path) as fh:
            text = fh.read()
    except Exception as exc:
        logging.warning("Failed to read MD output %s: %s", out_path, exc)
        return None

    lines = text.splitlines()
    if "Production metrics (step," in text or "OpenMM execution platform:" in text:
        implementation = "python_openmm"
    elif "[DEBUG] OpenMM force/nonbonded configuration:" in text or "OpenMM Molecular Dynamics Summary" in text:
        implementation = "cpp_openmm"
    else:
        implementation = "legacy_or_unknown"

    data = {
        "source": str(out_path),
        "implementation": implementation,
        "system": {},
        "protocol": {},
        "representative_particles": _parse_representative_particles(text),
        "representative_water_constraints": _parse_water_constraints(lines),
        "energy_decomposition_kj_per_mol": _parse_energy_decompositions(lines),
        "equilibration": {},
        "production": [],
        "timing": {}
    }
    system = data["system"]
    protocol = data["protocol"]

    if implementation == "python_openmm":
        system["platform"] = _first_match(text, r"OpenMM execution platform:\s*(\S+)", str)
        system["particles"] = _first_match(text, r"OpenMM particles:\s*(\d+)", int)
        system["constraints"] = _first_match(text, r"OpenMM constraints:\s*(\d+)", int)
        system["nonbonded_particles"] = _first_match(text, r"nonbonded particles\s*:\s*(\d+)", int)
        system["nonbonded_exceptions"] = _first_match(text, r"nonbonded exceptions\s*:\s*(\d+)", int)
        system["water_particles"] = _first_match(text, r"DivCon water particles aligned:\s*(\d+)", int)
        system["total_particle_charge_e"] = _first_match(
            text, rf"total particle charge \(e\):\s*({_FLOAT_RE})"
        )
        system["nonbonded_method"] = _first_match(text, r"^\s*method\s*:\s*(\S+)", str, re.MULTILINE)
        system["cutoff_nm"] = _first_match(text, rf"cutoff \(nm\)\s*:\s*({_FLOAT_RE})")
        system["ewald_error_tolerance"] = _first_match(
            text, rf"Ewald error tolerance\s*:\s*({_FLOAT_RE})"
        )
        dispersion = _first_match(text, r"dispersion correction\s*:\s*(True|False)", str, re.IGNORECASE)
        switching = _first_match(text, r"switching function\s*:\s*(True|False)", str, re.IGNORECASE)
        system["dispersion_correction"] = _parse_bool(dispersion) if dispersion is not None else None
        system["switching_function"] = _parse_bool(switching) if switching is not None else None
        system["initial_box_nm"] = _first_box(text, r"OpenMM initial box \(nm\):\s*\[([^\]]+)\]")

        protocol["minimization_limit"] = _first_match(
            text, r"Minimization step override:\s*(\d+)\s+maximum iterations", int
        )
        if protocol["minimization_limit"] is not None:
            protocol["minimization_limit_kind"] = "max_iterations"
        equil_steps = _first_match(
            text, r"Equilibration step override:\s*(\d+)\s+steps for both NVT and NPT", int
        )
        if equil_steps is not None:
            protocol["nvt_steps"] = equil_steps
            protocol["npt_steps"] = equil_steps
        protocol["production_steps"] = _first_match(
            text, r"Production step override:\s*(\d+)\s+steps", int
        )
        protocol["production_report_interval_steps"] = _first_match(
            text, r"Production report interval override:\s*(\d+)\s+steps", int
        )
        if protocol["production_report_interval_steps"] is not None:
            protocol["production_report_interval_source"] = "reported"
        protocol["production_duration_ps"] = _first_match(
            text, rf"Running Production NPT Simulation -\s*({_FLOAT_RE})\s+ps"
        )

        data["production"] = _parse_python_production(lines)
        pressures = {
            row["target_pressure_bar"] for row in data["production"]
            if row["target_pressure_bar"] is not None
        }
        if len(pressures) == 1:
            protocol["target_pressure_bar"] = next(iter(pressures))
        data["timing"] = _parse_python_timing(text)

    elif implementation == "cpp_openmm":
        system["particles"] = _first_match(
            text, r"^\[DEBUG\]\s+particles\s*:\s*(\d+)", int, re.MULTILINE
        )
        system["constraints"] = _first_match(
            text, r"^\[DEBUG\]\s+constraints\s*:\s*(\d+)", int, re.MULTILINE
        )
        system["harmonic_bonds"] = _first_match(
            text, r"^\[DEBUG\]\s+harmonic bonds\s*:\s*(\d+)", int, re.MULTILINE
        )
        system["harmonic_angles"] = _first_match(
            text, r"^\[DEBUG\]\s+harmonic angles\s*:\s*(\d+)", int, re.MULTILINE
        )
        system["periodic_torsions"] = _first_match(
            text, r"^\[DEBUG\]\s+periodic torsions\s*:\s*(\d+)", int, re.MULTILINE
        )
        system["nonbonded_particles"] = _first_match(
            text, r"^\[DEBUG\]\s+nonbonded particles\s*:\s*(\d+)", int, re.MULTILINE
        )
        system["nonbonded_exceptions"] = _first_match(
            text, r"^\[DEBUG\]\s+nonbonded exceptions\s*:\s*(\d+)", int, re.MULTILINE
        )
        system["zero_energy_exceptions"] = _first_match(
            text, r"^\[DEBUG\]\s+zero-energy\s*:\s*(\d+)", int, re.MULTILINE
        )
        system["nonzero_energy_exceptions"] = _first_match(
            text, r"^\[DEBUG\]\s+nonzero-energy\s*:\s*(\d+)", int, re.MULTILINE
        )
        system["water_molecules"] = _first_match(
            text, r"^\[DEBUG\]\s+water molecules\s*:\s*(\d+)", int, re.MULTILINE
        )
        system["water_particles"] = _first_match(
            text, r"^\[DEBUG\]\s+water particles\s*:\s*(\d+)", int, re.MULTILINE
        )
        system["total_particle_charge_e"] = _first_match(
            text, rf"^\[DEBUG\]\s+total particle charge\s*:\s*({_FLOAT_RE})\s+e", float, re.MULTILINE
        )
        system["nonbonded_method"] = _first_match(
            text, r"^\[DEBUG\]\s+nonbonded method\s*:\s*([^\s(]+)", str, re.MULTILINE
        )
        system["cutoff_nm"] = _first_match(
            text, rf"^\[DEBUG\]\s+cutoff \(nm\)\s*:\s*({_FLOAT_RE})", float, re.MULTILINE
        )
        system["ewald_error_tolerance"] = _first_match(
            text, rf"^\[DEBUG\]\s+Ewald error tolerance\s*:\s*({_FLOAT_RE})", float, re.MULTILINE
        )
        dispersion = _first_match(
            text, r"^\[DEBUG\]\s+dispersion correction\s*:\s*(true|false)", str,
            re.MULTILINE | re.IGNORECASE
        )
        switching = _first_match(
            text, r"^\[DEBUG\]\s+switching function\s*:\s*(true|false)", str,
            re.MULTILINE | re.IGNORECASE
        )
        system["dispersion_correction"] = _parse_bool(dispersion) if dispersion is not None else None
        system["switching_function"] = _parse_bool(switching) if switching is not None else None
        system["initial_box_nm"] = _first_box(text, r"mdSwitches\.periodicBox:\s*([^\n]+)")

        protocol["minimization_limit"] = _first_match(
            text, r"\[MD\] Minimizing before equilibration for\s*(\d+)\s+steps", int
        )
        if protocol["minimization_limit"] is not None:
            protocol["minimization_limit_kind"] = "reported_steps"
        nvt_match = re.search(
            rf"\[MD\] Starting NVT equilibration for\s*(\d+)\s+steps at\s*({_FLOAT_RE})\s+K",
            text
        )
        if nvt_match:
            protocol["nvt_steps"] = int(nvt_match.group(1))
            protocol["target_temperature_k"] = float(nvt_match.group(2))
        npt_match = re.search(
            rf"\[MD\] Starting NPT equilibration for\s*(\d+)\s+steps at\s*({_FLOAT_RE})\s+K "
            rf"and\s*({_FLOAT_RE})\s+bar",
            text
        )
        if npt_match:
            protocol["npt_steps"] = int(npt_match.group(1))
            protocol["target_temperature_k"] = float(npt_match.group(2))
            protocol["target_pressure_bar"] = float(npt_match.group(3))
        production_match = re.search(
            rf"\[MD\] Starting production MD for\s*(\d+)\s+steps\s*\(({_FLOAT_RE})\s+ps\)",
            text
        )
        if production_match:
            protocol["production_steps"] = int(production_match.group(1))
            protocol["production_duration_ps"] = float(production_match.group(2))

        nvt_energy = _first_match(
            text, rf"\[MD\] NVT equilibration done\. Potential Energy:\s*({_FLOAT_RE})"
        )
        if nvt_energy is not None:
            data["equilibration"]["nvt"] = {"potential_energy_kj_per_mol": nvt_energy}
        npt_done = re.search(
            rf"\[MD\] NPT equilibration done\. Potential Energy:\s*({_FLOAT_RE}).*?"
            rf"\[MD\] Final box after NPT equilibration \(nm\):\s*"
            rf"({_FLOAT_RE})\s+({_FLOAT_RE})\s+({_FLOAT_RE})\s+"
            rf"Volume\(nm\^3\):\s*({_FLOAT_RE})\s+density \(g/ml\):\s*({_FLOAT_RE})\s+T:\s*({_FLOAT_RE})",
            text,
            re.DOTALL
        )
        if npt_done:
            data["equilibration"]["npt"] = {
                "potential_energy_kj_per_mol": float(npt_done.group(1)),
                "box_nm": [float(npt_done.group(i)) for i in (2, 3, 4)],
                "volume_nm3": float(npt_done.group(5)),
                "density_g_per_ml": float(npt_done.group(6)),
                "temperature_k": float(npt_done.group(7))
            }

        data["production"] = _parse_cpp_production(lines)
        data["timing"] = _parse_cpp_timing(text)

    else:
        data["production"] = _parse_legacy_production(lines)
        if data["production"]:
            data["implementation"] = "legacy_three_column"
        else:
            logging.warning("No recognized MD data found in %s", out_path)
            return None

    initial_pe = _first_match(
        text, rf"(?:OpenMM initial potential energy \(kJ/mol\):|\[DEBUG\] Initial Potential Energy:)\s*({_FLOAT_RE})",
        float,
        re.IGNORECASE
    )
    if initial_pe is not None:
        data["initial_potential_energy_kj_per_mol"] = initial_pe

    post_min_pe = _first_match(
        text, rf"(?:Post-minimization Potential Energy:|\[MD\] Post-minimization Potential Energy:)\s*({_FLOAT_RE})",
        float
    )
    if post_min_pe is not None:
        data["post_minimization_potential_energy_kj_per_mol"] = post_min_pe

    _infer_production_timing(protocol, data["production"])
    return data


def _symmetric_relative_difference_percent(value1, value2, signed=False):
    """Symmetric percentage difference; unlike percent-vs-input2 it is invariant to input order."""
    if value1 is None or value2 is None:
        return None
    denominator = abs(value1) + abs(value2)
    if denominator <= 1.0e-15:
        return 0.0 if value1 == value2 else None
    numerator = 2.0 * (value1 - value2)
    result = 100.0 * numerator / denominator
    return result if signed else abs(result)


def _numeric_comparison(value1, value2) -> dict:
    delta = float(value1 - value2)
    return {
        "input1": float(value1),
        "input2": float(value2),
        "delta_input1_minus_input2": delta,
        "absolute_difference": abs(delta),
        "symmetric_relative_difference_percent": _symmetric_relative_difference_percent(value1, value2)
    }


def _compare_mapping(mapping1: dict, mapping2: dict, fields: list[str]) -> dict:
    comparison = {}
    for field in fields:
        value1 = mapping1.get(field)
        value2 = mapping2.get(field)
        if value1 is None or value2 is None:
            continue

        if isinstance(value1, bool) or isinstance(value1, str) or isinstance(value1, int):
            comparison[field] = {
                "input1": value1,
                "input2": value2,
                "match": value1 == value2
            }
        elif isinstance(value1, (float, np.floating)) and isinstance(value2, (float, np.floating)):
            comparison[field] = _numeric_comparison(value1, value2)
        elif isinstance(value1, list) and isinstance(value2, list) and len(value1) == len(value2):
            array1 = np.asarray(value1, dtype=float)
            array2 = np.asarray(value2, dtype=float)
            delta = array1 - array2
            comparison[field] = {
                "input1": [float(v) for v in array1],
                "input2": [float(v) for v in array2],
                "delta_input1_minus_input2": [float(v) for v in delta],
                "max_absolute_difference": float(np.max(np.abs(delta))),
                "rmse_difference": float(np.sqrt(np.mean(delta * delta)))
            }

    return comparison


def _compare_energy_sections(data1: dict, data2: dict) -> dict:
    result = {}
    sections1 = data1.get("energy_decomposition_kj_per_mol", {})
    sections2 = data2.get("energy_decomposition_kj_per_mol", {})

    for section_name in ("initial", "post_minimization"):
        values1 = sections1.get(section_name)
        values2 = sections2.get(section_name)
        if not values1 or not values2:
            continue

        section = {}
        for component in ("bond", "angle", "torsion", "nonbonded", "four_part_sum", "total", "residual"):
            if component in values1 and component in values2:
                section[component] = _numeric_comparison(values1[component], values2[component])
        result[section_name] = section

    return result


def _compare_representative_particles(data1: dict, data2: dict) -> dict:
    result = {}
    particles1 = data1.get("representative_particles", {})
    particles2 = data2.get("representative_particles", {})

    for role in ("water_o", "water_h1", "water_h2", "sodium", "chloride"):
        p1 = particles1.get(role)
        p2 = particles2.get(role)
        if not p1 or not p2:
            continue

        result[role] = {
            "index": {"input1": p1["index"], "input2": p2["index"], "match": p1["index"] == p2["index"]},
            "divcon_type": {
                "input1": p1["divcon_type"],
                "input2": p2["divcon_type"],
                "match": p1["divcon_type"] == p2["divcon_type"]
            },
            "openmm_charge_e": _numeric_comparison(p1["openmm_charge_e"], p2["openmm_charge_e"]),
            "sigma_nm": _numeric_comparison(p1["sigma_nm"], p2["sigma_nm"]),
            "epsilon_kj_per_mol": _numeric_comparison(
                p1["epsilon_kj_per_mol"], p2["epsilon_kj_per_mol"]
            )
        }

    return result


def _compare_water_constraints(data1: dict, data2: dict) -> dict:
    distances1 = sorted(
        item["distance_nm"] for item in data1.get("representative_water_constraints", [])
    )
    distances2 = sorted(
        item["distance_nm"] for item in data2.get("representative_water_constraints", [])
    )
    if not distances1 or len(distances1) != len(distances2):
        return {}

    return {
        "n_constraints": len(distances1),
        "sorted_distance_comparisons_nm": [
            _numeric_comparison(value1, value2)
            for value1, value2 in zip(distances1, distances2)
        ],
        "max_absolute_distance_difference_nm": float(
            np.max(np.abs(np.asarray(distances1) - np.asarray(distances2)))
        )
    }


def _production_rows_by_step(data: dict) -> dict[int, dict]:
    return {int(row["step"]): row for row in data.get("production", [])}


def _production_observable_stats(values1: np.ndarray, values2: np.ndarray) -> dict:
    delta = values1 - values2
    relative = np.array([
        _symmetric_relative_difference_percent(a, b) for a, b in zip(values1, values2)
    ], dtype=float)
    return {
        "n": int(len(delta)),
        "mean_input1": float(np.mean(values1)),
        "mean_input2": float(np.mean(values2)),
        "mean_signed_difference": float(np.mean(delta)),
        "mean_absolute_difference": float(np.mean(np.abs(delta))),
        "rmse_difference": float(np.sqrt(np.mean(delta * delta))),
        "max_absolute_difference": float(np.max(np.abs(delta))),
        "mean_absolute_symmetric_relative_difference_percent": float(np.mean(relative)),
        "max_absolute_symmetric_relative_difference_percent": float(np.max(relative))
    }


def compare_production(data1: dict, data2: dict) -> dict:
    rows1 = _production_rows_by_step(data1)
    rows2 = _production_rows_by_step(data2)
    common_steps = sorted(set(rows1) & set(rows2))
    observables = (
        "potential_energy_kj_per_mol",
        "kinetic_energy_kj_per_mol",
        "temperature_k",
        "volume_nm3",
        "density_g_per_ml"
    )

    result = {
        "matched_steps": common_steps,
        "n_matched_steps": len(common_steps),
        "observables": {},
        "matched_rows": []
    }

    for observable in observables:
        paired = [
            (rows1[step].get(observable), rows2[step].get(observable))
            for step in common_steps
        ]
        paired = [(a, b) for a, b in paired if a is not None and b is not None]
        if not paired:
            continue
        values1 = np.asarray([a for a, _ in paired], dtype=float)
        values2 = np.asarray([b for _, b in paired], dtype=float)
        result["observables"][observable] = _production_observable_stats(values1, values2)

    for step in common_steps:
        row = {"step": step}
        for observable in observables:
            value1 = rows1[step].get(observable)
            value2 = rows2[step].get(observable)
            if value1 is not None and value2 is not None:
                row[observable] = _numeric_comparison(value1, value2)
        result["matched_rows"].append(row)

    return result


def compare_timing(data1: dict, data2: dict) -> dict:
    """Compare available wall-clock timings while preserving missing values."""
    fields = (
        "minimization_wall_time_s",
        "nvt_equilibration_wall_time_s",
        "npt_equilibration_wall_time_s",
        "production_wall_time_s",
        "timed_md_stages_wall_time_s",
        "total_execution_wall_time_s",
        "non_production_execution_wall_time_s",
        "production_fraction_of_total_percent",
    )
    timing1 = data1.get("timing", {})
    timing2 = data2.get("timing", {})
    result = {}

    for field in fields:
        value1 = timing1.get(field)
        value2 = timing2.get(field)
        if value1 is None and value2 is None:
            continue

        values = {
            "input1": float(value1) if value1 is not None else None,
            "input2": float(value2) if value2 is not None else None,
        }
        if value1 is not None and value2 is not None:
            values.update(_numeric_comparison(value1, value2))
        result[field] = values

    return result


def compare_openmm_outputs(data1: dict, data2: dict) -> dict:
    system_fields = [
        "particles", "constraints", "nonbonded_particles", "nonbonded_exceptions",
        "water_particles", "total_particle_charge_e", "nonbonded_method", "cutoff_nm",
        "ewald_error_tolerance", "dispersion_correction", "switching_function", "initial_box_nm"
    ]
    protocol_fields = [
        "minimization_limit", "nvt_steps", "npt_steps", "production_steps",
        "production_report_interval_steps", "target_temperature_k", "target_pressure_bar",
        "production_duration_ps", "timestep_ps"
    ]

    return {
        "system": _compare_mapping(data1.get("system", {}), data2.get("system", {}), system_fields),
        "protocol": _compare_mapping(data1.get("protocol", {}), data2.get("protocol", {}), protocol_fields),
        "representative_particles": _compare_representative_particles(data1, data2),
        "representative_water_constraints": _compare_water_constraints(data1, data2),
        "energy_decomposition_kj_per_mol": _compare_energy_sections(data1, data2),
        "timing": compare_timing(data1, data2),
        "production": compare_production(data1, data2)
    }


def compatibility_out_summary(data: dict):
    """Return the historical OUT-summary fields for existing JSON consumers."""
    if not data or not data.get("production"):
        return None

    rows = data["production"]
    energies = [
        row["potential_energy_kj_per_mol"]
        for row in rows
        if row["potential_energy_kj_per_mol"] is not None
    ]
    temps = [row["temperature_k"] for row in rows if row["temperature_k"] is not None]
    if not energies or not temps:
        return None

    times = [row["time_ps"] for row in rows if row.get("time_ps") is not None]
    total_ps = max(times) if times else data.get("protocol", {}).get("production_duration_ps")
    if total_ps is None and len(rows) > 1:
        steps = [row["step"] for row in rows]
        interval = steps[1] - steps[0]
        total_ps = steps[-1] * (interval / 1000.0)

    return {
        "n_steps": int(len(rows)),
        "total_ps": float(total_ps) if total_ps is not None else 0.0,
        "avg_energy": float(np.mean(energies)),
        "std_energy": float(np.std(energies)),
        "avg_temp": float(np.mean(temps)),
        "std_temp": float(np.std(temps))
    }


def _format_implementation(value: str) -> str:
    return {
        "python_openmm": "Python/OpenMM",
        "cpp_openmm": "DivCon/C++/OpenMM",
        "legacy_three_column": "Legacy three-column OUT"
    }.get(value, value)


def print_single_output_summary(label: str, data: dict) -> None:
    print(f"\nMD output summary: {label}")
    if not data:
        print("  No recognized MD output data found.")
        return

    print(f"  Detected format            : {_format_implementation(data['implementation'])}")
    system = data.get("system", {})
    protocol = data.get("protocol", {})
    for field, title in (
        ("platform", "OpenMM platform"),
        ("particles", "Particles"),
        ("constraints", "Constraints"),
        ("nonbonded_exceptions", "Nonbonded exceptions"),
        ("water_particles", "Water particles")
    ):
        if system.get(field) is not None:
            print(f"  {title:<27}: {system[field]}")

    if protocol.get("production_steps") is not None:
        print(f"  {'Production steps':<27}: {protocol['production_steps']}")
    if protocol.get("production_duration_ps") is not None:
        print(f"  {'Production duration':<27}: {protocol['production_duration_ps']:.3f} ps")
    print(f"  {'Production rows':<27}: {len(data.get('production', []))}")

    timing = data.get("timing", {})
    if timing.get("production_wall_time_s") is not None:
        print(f"  {'Production wall time':<27}: {timing['production_wall_time_s']:.6f} s")
    if timing.get("total_execution_wall_time_s") is not None:
        print(f"  {'Total execution wall time':<27}: {timing['total_execution_wall_time_s']:.6f} s")
    if timing.get("production_fraction_of_total_percent") is not None:
        print(
            f"  {'Production / total':<27}: "
            f"{timing['production_fraction_of_total_percent']:.3f}%"
        )


def print_system_comparison(label1: str, label2: str, comparison: dict) -> None:
    print("\nOpenMM System / Hamiltonian Validation")
    print("=" * 72)
    system = comparison.get("system", {})
    if not system:
        print("  No common system-configuration fields were available for comparison.")
        return

    for field, values in system.items():
        name = field.replace("_", " ")
        if "match" in values:
            status = "MATCH" if values["match"] else "DIFFER"
            print(f"  {name:<30}: {status:<6}  {values['input1']}  vs  {values['input2']}")
        elif "max_absolute_difference" in values:
            print(f"  {name:<30}: max |Δ|={values['max_absolute_difference']:.6g}")
        else:
            print(
                f"  {name:<30}: Δ({label1}-{label2})="
                f"{values['delta_input1_minus_input2']:+.6g}"
            )


def print_representative_parameter_comparison(comparison: dict) -> None:
    particles = comparison.get("representative_particles", {})
    constraints = comparison.get("representative_water_constraints", {})
    if not particles and not constraints:
        return

    print("\nRepresentative Effective Parameter Validation")
    print("=" * 72)
    for role, values in particles.items():
        charge = values["openmm_charge_e"]
        sigma = values["sigma_nm"]
        epsilon = values["epsilon_kj_per_mol"]
        index_status = "MATCH" if values["index"]["match"] else "DIFFER"
        type_status = "MATCH" if values["divcon_type"]["match"] else "DIFFER"
        print(
            f"  {role:<10} index={index_status:<6} type={type_status:<6} "
            f"Δq={charge['delta_input1_minus_input2']:+.3e} e  "
            f"Δsigma={sigma['delta_input1_minus_input2']:+.3e} nm  "
            f"Δepsilon={epsilon['delta_input1_minus_input2']:+.3e} kJ/mol"
        )

    if constraints:
        print(
            f"  water constraints: n={constraints['n_constraints']}  "
            f"max |Δdistance|={constraints['max_absolute_distance_difference_nm']:.3e} nm"
        )


def print_energy_comparison(label1: str, label2: str, comparison: dict) -> None:
    energy = comparison.get("energy_decomposition_kj_per_mol", {})
    component_labels = {
        "bond": "Bond",
        "angle": "Angle",
        "torsion": "Torsion",
        "nonbonded": "Nonbonded",
        "four_part_sum": "4-part sum",
        "total": "Total",
        "residual": "Residual"
    }

    for section_name, title in (
        ("initial", "OpenMM Energy Validation — Initial (diagnostic)"),
        ("post_minimization", "OpenMM Energy Validation — Post-minimization")
    ):
        section = energy.get(section_name)
        if not section:
            continue

        print(f"\n{title}")
        print("=" * len(title))
        print(
            f"  {'Component':<14} {label1[:18]:>18} {label2[:18]:>18} "
            f"{'Δ (1-2)':>14} {'sym |Δ| %':>11}"
        )
        print("  " + "-" * 79)
        for component, display_name in component_labels.items():
            if component not in section:
                continue
            values = section[component]
            relative = values["symmetric_relative_difference_percent"]
            relative_text = f"{relative:.4f}" if relative is not None else "n/a"
            print(
                f"  {display_name:<14} {values['input1']:>18.4f} {values['input2']:>18.4f} "
                f"{values['delta_input1_minus_input2']:>+14.4f} {relative_text:>11}"
            )


def print_protocol_comparison(label1: str, label2: str, comparison: dict) -> None:
    protocol = comparison.get("protocol", {})
    if not protocol:
        return

    print("\nMD Protocol Validation")
    print("=" * 72)
    for field, values in protocol.items():
        name = field.replace("_", " ")
        if "match" in values:
            status = "MATCH" if values["match"] else "DIFFER"
            print(f"  {name:<34}: {status:<6}  {values['input1']}  vs  {values['input2']}")
        else:
            print(
                f"  {name:<34}: {values['input1']:.6g} vs {values['input2']:.6g} "
                f"(Δ {values['delta_input1_minus_input2']:+.6g})"
            )


def print_timing_comparison(label1: str, label2: str, comparison: dict) -> None:
    timing = comparison.get("timing", {})
    if not timing:
        return

    print("\nTiming Analysis")
    print("=" * 72)
    print(
        f"  {'Timing metric':<31} {label1[:16]:>16} {label2[:16]:>16} "
        f"{'Delta (1-2)':>14} {'sym |Delta| %':>13}"
    )
    print("  " + "-" * 94)

    labels = (
        ("minimization_wall_time_s", "Minimization wall time", "s"),
        ("nvt_equilibration_wall_time_s", "NVT equilibration wall time", "s"),
        ("npt_equilibration_wall_time_s", "NPT equilibration wall time", "s"),
        ("production_wall_time_s", "Production MD wall time", "s"),
        ("timed_md_stages_wall_time_s", "Timed MD stages (complete sum)", "s"),
        ("total_execution_wall_time_s", "Total execution wall time", "s"),
        ("non_timed_md_execution_wall_time_s", "Total minus timed MD stages", "s"),
        ("timed_md_fraction_of_total_percent", "Timed MD stages / total execution", "%"),
        ("non_production_execution_wall_time_s", "Total minus production MD", "s"),
        ("production_fraction_of_total_percent", "Production MD / total execution", "%"),
    )

    def fmt(value, unit):
        if value is None:
            return "n/a"
        if unit == "%":
            return f"{value:.3f}%"
        return f"{value:.6f}"

    for field, title, unit in labels:
        values = timing.get(field)
        if not values:
            continue
        value1 = values.get("input1")
        value2 = values.get("input2")
        delta = values.get("delta_input1_minus_input2")
        relative = values.get("symmetric_relative_difference_percent")
        delta_text = "n/a" if delta is None else f"{delta:+.6f}"
        relative_text = "n/a" if relative is None else f"{relative:.3f}"
        print(
            f"  {title:<31} {fmt(value1, unit):>16} {fmt(value2, unit):>16} "
            f"{delta_text:>14} {relative_text:>13}"
        )

    print(f"  Differences are defined as {label1} - {label2}.")
    print(
        "  Timed MD stages / total execution sums minimization, NVT, NPT, and production timers "
        "when all four are available."
    )
    print(
        "  Production MD / total execution uses the production timer alone divided by the whole-process timer."
    )
    if any(
        timing.get(field, {}).get("input1") is None or timing.get(field, {}).get("input2") is None
        for field in ("total_execution_wall_time_s", "production_fraction_of_total_percent")
        if field in timing
    ):
        print(
            "  n/a means that the corresponding screenout does not report the required whole-process timing; "
            "stage sums are not substituted for total execution time."
        )


def print_production_comparison(label1: str, label2: str, comparison: dict) -> None:
    production = comparison.get("production", {})
    print("\nProduction Thermodynamic Validation")
    print("=" * 72)
    if not production or not production.get("matched_steps"):
        print("  No common production steps with thermodynamic data were found.")
        return

    steps = production["matched_steps"]
    print(f"  Matched production checkpoints: {len(steps)} ({steps[0]} through {steps[-1]})")
    print(f"  Differences are defined as {label1} - {label2}.")
    print("  Relative values use symmetric percent difference: 200*|x1-x2|/(|x1|+|x2|).")
    print()
    print(f"  {'Observable':<22} {'mean |Δ|':>14} {'RMSE Δ':>14} {'mean sym %':>12} {'max sym %':>12}")
    print("  " + "-" * 78)

    names = {
        "potential_energy_kj_per_mol": "Potential energy",
        "kinetic_energy_kj_per_mol": "Kinetic energy",
        "temperature_k": "Temperature",
        "volume_nm3": "Volume",
        "density_g_per_ml": "Density"
    }
    for key, title in names.items():
        stats = production.get("observables", {}).get(key)
        if not stats:
            continue
        print(
            f"  {title:<22} {stats['mean_absolute_difference']:>14.6g} "
            f"{stats['rmse_difference']:>14.6g} "
            f"{stats['mean_absolute_symmetric_relative_difference_percent']:>12.4f} "
            f"{stats['max_absolute_symmetric_relative_difference_percent']:>12.4f}"
        )


def plot_production_validation(
    data1: dict, data2: dict, comparison: dict, label1: str, label2: str, out_prefix: str
) -> None:
    production = comparison.get("production", {})
    common_steps = production.get("matched_steps", [])
    if not common_steps:
        return

    rows1 = _production_rows_by_step(data1)
    rows2 = _production_rows_by_step(data2)
    plot_specs = (
        ("potential_energy_kj_per_mol", "Potential energy (kJ/mol)", "potential_energy_validation"),
        ("temperature_k", "Temperature (K)", "temperature_validation"),
        ("volume_nm3", "Volume (nm³)", "volume_validation"),
        ("density_g_per_ml", "Density (g/mL)", "density_validation")
    )

    for key, ylabel, suffix in plot_specs:
        points = [
            (step, rows1[step].get(key), rows2[step].get(key))
            for step in common_steps
            if rows1[step].get(key) is not None and rows2[step].get(key) is not None
        ]
        if not points:
            continue

        x = np.asarray([point[0] for point in points], dtype=int)
        values1 = np.asarray([point[1] for point in points], dtype=float)
        values2 = np.asarray([point[2] for point in points], dtype=float)
        plt.figure()
        plt.plot(x, values1, marker="o", label=label1)
        plt.plot(x, values2, marker="o", label=label2)
        plt.xlabel("Production step")
        plt.ylabel(ylabel)
        plt.title(ylabel.split(" (")[0] + " validation")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{out_prefix}_{suffix}.pdf", dpi=CONFIG["plot"]["dpi"])
        plt.close()

    pe_points = [
        (
            step,
            rows1[step].get("potential_energy_kj_per_mol"),
            rows2[step].get("potential_energy_kj_per_mol")
        )
        for step in common_steps
        if rows1[step].get("potential_energy_kj_per_mol") is not None
        and rows2[step].get("potential_energy_kj_per_mol") is not None
    ]
    if pe_points:
        x = np.asarray([point[0] for point in pe_points], dtype=int)
        signed_relative = np.asarray([
            _symmetric_relative_difference_percent(point[1], point[2], signed=True)
            for point in pe_points
        ])
        plt.figure()
        plt.axhline(0.0, linewidth=0.8)
        plt.plot(x, signed_relative, marker="o")
        plt.xlabel("Production step")
        plt.ylabel(f"Signed symmetric ΔPE (%) [{label1} - {label2}]")
        plt.title("Potential-energy difference")
        plt.tight_layout()
        plt.savefig(f"{out_prefix}_potential_energy_delta.pdf", dpi=CONFIG["plot"]["dpi"])
        plt.close()


# ───────── SERIALIZED OPENMM XML VALIDATION ───────────────────────────────────
_OPENMM_SERIALIZED_FILENAMES = {
    "system": "openmm_system.xml",
    "integrator": "openmm_integrator.xml",
    "initial_state": "openmm_initial_state.xml",
}


def _openmm_float(value, target_unit=None) -> float:
    """Convert an OpenMM scalar/Quantity to a plain float in a documented unit."""
    if target_unit is not None and hasattr(value, "value_in_unit"):
        return float(value.value_in_unit(target_unit))
    try:
        return float(value)
    except (TypeError, ValueError):
        if hasattr(value, "_value"):
            return float(value._value)
        raise


def _openmm_numbers_match(value1: float, value2: float, absolute_tolerance=None,
                           relative_tolerance=None) -> bool:
    cfg = CONFIG["openmm_xml_validation"]
    atol = cfg["absolute_tolerance"] if absolute_tolerance is None else absolute_tolerance
    rtol = cfg["relative_tolerance"] if relative_tolerance is None else relative_tolerance
    return abs(value1 - value2) <= atol + rtol * max(abs(value1), abs(value2))


def _openmm_scalar_comparison(value1, value2, absolute_tolerance=None,
                              relative_tolerance=None) -> dict:
    if isinstance(value1, (bool, str, int)) or isinstance(value2, (bool, str, int)):
        return {"input1": value1, "input2": value2, "match": value1 == value2}

    value1 = float(value1)
    value2 = float(value2)
    return {
        "input1": value1,
        "input2": value2,
        "delta_input1_minus_input2": value1 - value2,
        "absolute_difference": abs(value1 - value2),
        "match": _openmm_numbers_match(
            value1, value2, absolute_tolerance, relative_tolerance
        )
    }


def _openmm_compare_settings(settings1: dict, settings2: dict) -> dict:
    fields = {}
    missing1 = []
    missing2 = []
    n_different = 0

    for field in sorted(set(settings1) | set(settings2)):
        if field not in settings1:
            missing1.append(field)
            n_different += 1
            continue
        if field not in settings2:
            missing2.append(field)
            n_different += 1
            continue
        comparison = _openmm_scalar_comparison(settings1[field], settings2[field])
        fields[field] = comparison
        if not comparison["match"]:
            n_different += 1

    return {
        "fields": fields,
        "missing_from_input1": missing1,
        "missing_from_input2": missing2,
        "n_different": n_different,
        "match": n_different == 0
    }


def _openmm_compare_keyed_records(records1: list, records2: list, field_names: tuple[str, ...]) -> dict:
    """Compare parameter records as keyed multisets, independent of record order."""
    grouped1 = defaultdict(list)
    grouped2 = defaultdict(list)
    for key, values in records1:
        grouped1[tuple(key)].append(tuple(float(v) for v in values))
    for key, values in records2:
        grouped2[tuple(key)].append(tuple(float(v) for v in values))

    keys1 = set(grouped1)
    keys2 = set(grouped2)
    only1 = sorted(keys1 - keys2)
    only2 = sorted(keys2 - keys1)
    multiplicity_mismatches = []
    parameter_mismatches = []
    max_abs = {name: 0.0 for name in field_names}
    field_accumulators = {
        name: {
            "compared_count": 0,
            "mismatch_count": 0,
            "sum_signed": 0.0,
            "sum_absolute": 0.0,
            "sum_squared": 0.0,
            "max_absolute": 0.0,
        }
        for name in field_names
    }

    cfg = CONFIG["openmm_xml_validation"]
    atol = cfg["absolute_tolerance"]
    rtol = cfg["relative_tolerance"]

    for key in sorted(keys1 & keys2):
        values1 = sorted(grouped1[key])
        values2 = sorted(grouped2[key])
        if len(values1) != len(values2):
            multiplicity_mismatches.append({
                "key": list(key),
                "count_input1": len(values1),
                "count_input2": len(values2)
            })

        for params1, params2 in zip(values1, values2):
            differences = {}
            mismatch = False
            tolerance_score = 0.0
            for field, value1, value2 in zip(field_names, params1, params2):
                delta = value1 - value2
                absolute = abs(delta)
                max_abs[field] = max(max_abs[field], absolute)
                allowed = atol + rtol * max(abs(value1), abs(value2))
                score = absolute / allowed if allowed > 0.0 else (0.0 if absolute == 0.0 else float("inf"))
                tolerance_score = max(tolerance_score, score)
                match = absolute <= allowed
                if not match:
                    mismatch = True

                acc = field_accumulators[field]
                acc["compared_count"] += 1
                acc["sum_signed"] += delta
                acc["sum_absolute"] += absolute
                acc["sum_squared"] += delta * delta
                acc["max_absolute"] = max(acc["max_absolute"], absolute)
                if not match:
                    acc["mismatch_count"] += 1

                differences[field] = {
                    "input1": value1,
                    "input2": value2,
                    "delta_input1_minus_input2": delta,
                    "absolute_difference": absolute,
                    "match": match
                }

            if mismatch:
                parameter_mismatches.append({
                    "key": list(key),
                    "fields": differences,
                    "tolerance_score": tolerance_score
                })

    parameter_mismatches.sort(key=lambda item: item["tolerance_score"], reverse=True)
    top_n = CONFIG["openmm_xml_validation"]["top_n_differences"]

    field_statistics = {}
    top_parameter_mismatches_by_field = {}
    for field, acc in field_accumulators.items():
        count = acc["compared_count"]
        field_statistics[field] = {
            "compared_count": count,
            "mismatch_count": acc["mismatch_count"],
            "mean_signed_difference": acc["sum_signed"] / count if count else 0.0,
            "mean_absolute_difference": acc["sum_absolute"] / count if count else 0.0,
            "rms_difference": float(np.sqrt(acc["sum_squared"] / count)) if count else 0.0,
            "max_absolute_difference": acc["max_absolute"],
        }

        field_differences = []
        for mismatch in parameter_mismatches:
            values = mismatch["fields"][field]
            if values["match"]:
                continue
            field_differences.append({
                "key": mismatch["key"],
                "input1": values["input1"],
                "input2": values["input2"],
                "delta_input1_minus_input2": values["delta_input1_minus_input2"],
                "absolute_difference": values["absolute_difference"],
            })
        field_differences.sort(key=lambda item: item["absolute_difference"], reverse=True)
        top_parameter_mismatches_by_field[field] = field_differences[:top_n]

    return {
        "count_input1": len(records1),
        "count_input2": len(records2),
        "count_match": len(records1) == len(records2),
        "keys_only_input1_count": len(only1),
        "keys_only_input2_count": len(only2),
        "keys_only_input1": [list(key) for key in only1[:top_n]],
        "keys_only_input2": [list(key) for key in only2[:top_n]],
        "multiplicity_mismatch_count": len(multiplicity_mismatches),
        "multiplicity_mismatches": multiplicity_mismatches[:top_n],
        "parameter_mismatch_count": len(parameter_mismatches),
        "top_parameter_mismatches": parameter_mismatches[:top_n],
        "top_parameter_mismatches_by_field": top_parameter_mismatches_by_field,
        "field_statistics": field_statistics,
        "max_absolute_difference_by_field": max_abs,
        "match": (
            len(records1) == len(records2)
            and not only1
            and not only2
            and not multiplicity_mismatches
            and not parameter_mismatches
        )
    }


def _openmm_path_metadata(path: Path | None) -> dict:
    if path is None:
        return {"path": None, "exists": False}
    path = Path(path)
    result = {"path": str(path), "exists": path.is_file()}
    if path.is_file():
        payload = path.read_bytes()
        result["size_bytes"] = len(payload)
        result["sha256"] = hashlib.sha256(payload).hexdigest()
    return result


def _discover_openmm_serialized_paths(out_path: Path | None) -> dict[str, Path | None]:
    if out_path is None:
        return {key: None for key in _OPENMM_SERIALIZED_FILENAMES}
    base = Path(out_path).parent
    return {key: base / filename for key, filename in _OPENMM_SERIALIZED_FILENAMES.items()}


def _load_openmm_serialized_object(path: Path):
    with open(path, "r") as handle:
        return mm.XmlSerializer.deserialize(handle.read())


def _openmm_force_inventory(system) -> dict[str, list]:
    inventory = defaultdict(list)
    for force_index in range(system.getNumForces()):
        force = system.getForce(force_index)
        inventory[type(force).__name__].append((force_index, force))
    return dict(inventory)


def _openmm_force_metadata(force1, force2) -> dict:
    metadata = {
        "force_group": _openmm_scalar_comparison(
            int(force1.getForceGroup()), int(force2.getForceGroup())
        ),
        "name": _openmm_scalar_comparison(str(force1.getName()), str(force2.getName()))
    }
    if hasattr(force1, "usesPeriodicBoundaryConditions") and hasattr(force2, "usesPeriodicBoundaryConditions"):
        metadata["uses_periodic_boundary_conditions"] = _openmm_scalar_comparison(
            bool(force1.usesPeriodicBoundaryConditions()),
            bool(force2.usesPeriodicBoundaryConditions())
        )
    return metadata


def _openmm_force_metadata_difference_count(metadata: dict) -> int:
    # Force group/name differences are useful diagnostics but do not by themselves
    # alter the Hamiltonian.  PBC usage is physics-affecting and is counted by the
    # force-specific comparison where applicable.
    return sum(
        1 for key in ("force_group", "name")
        if key in metadata and not metadata[key]["match"]
    )


def _openmm_harmonic_bond_records(force) -> list:
    records = []
    k_unit = openmm_unit.kilojoule_per_mole / openmm_unit.nanometer**2
    for index in range(force.getNumBonds()):
        p1, p2, length, k = force.getBondParameters(index)
        key = tuple(sorted((int(p1), int(p2))))
        records.append((key, (
            _openmm_float(length, openmm_unit.nanometer),
            _openmm_float(k, k_unit)
        )))
    return records


def _openmm_harmonic_angle_records(force) -> list:
    records = []
    k_unit = openmm_unit.kilojoule_per_mole / openmm_unit.radian**2
    for index in range(force.getNumAngles()):
        p1, p2, p3, theta, k = force.getAngleParameters(index)
        forward = (int(p1), int(p2), int(p3))
        reverse = (int(p3), int(p2), int(p1))
        key = min(forward, reverse)
        records.append((key, (
            _openmm_float(theta, openmm_unit.radian),
            _openmm_float(k, k_unit)
        )))
    return records


def _openmm_periodic_torsion_records(force) -> list:
    records = []
    for index in range(force.getNumTorsions()):
        p1, p2, p3, p4, periodicity, phase, k = force.getTorsionParameters(index)
        # Keep atom order explicit.  Reversing a torsion can change the sign
        # convention for its phase, so treating reversed tuples as identical can
        # hide a real representation difference.
        key = (int(p1), int(p2), int(p3), int(p4), int(periodicity))
        records.append((key, (
            _openmm_float(phase, openmm_unit.radian),
            _openmm_float(k, openmm_unit.kilojoule_per_mole)
        )))
    return records


def _openmm_nonbonded_method_name(value) -> str:
    method_map = {
        int(mm.NonbondedForce.NoCutoff): "NoCutoff",
        int(mm.NonbondedForce.CutoffNonPeriodic): "CutoffNonPeriodic",
        int(mm.NonbondedForce.CutoffPeriodic): "CutoffPeriodic",
        int(mm.NonbondedForce.Ewald): "Ewald",
        int(mm.NonbondedForce.PME): "PME",
        int(mm.NonbondedForce.LJPME): "LJPME",
    }
    return method_map.get(int(value), str(int(value)))


def _openmm_nonbonded_settings(force) -> dict:
    settings = {
        "method": _openmm_nonbonded_method_name(force.getNonbondedMethod()),
        "cutoff_nm": _openmm_float(force.getCutoffDistance(), openmm_unit.nanometer),
        "reaction_field_dielectric": float(force.getReactionFieldDielectric()),
        "ewald_error_tolerance": float(force.getEwaldErrorTolerance()),
        "use_dispersion_correction": bool(force.getUseDispersionCorrection()),
        "use_switching_function": bool(force.getUseSwitchingFunction()),
        "exceptions_use_periodic_boundary_conditions": bool(
            force.getExceptionsUsePeriodicBoundaryConditions()
        )
    }
    if force.getUseSwitchingFunction():
        settings["switching_distance_nm"] = _openmm_float(
            force.getSwitchingDistance(), openmm_unit.nanometer
        )

    alpha, nx, ny, nz = force.getPMEParameters()
    settings.update({
        "pme_alpha_per_nm": _openmm_float(alpha, openmm_unit.nanometer**-1),
        "pme_grid_x": int(nx),
        "pme_grid_y": int(ny),
        "pme_grid_z": int(nz)
    })
    if hasattr(force, "getLJPMEParameters"):
        alpha_lj, nx_lj, ny_lj, nz_lj = force.getLJPMEParameters()
        settings.update({
            "ljpme_alpha_per_nm": _openmm_float(alpha_lj, openmm_unit.nanometer**-1),
            "ljpme_grid_x": int(nx_lj),
            "ljpme_grid_y": int(ny_lj),
            "ljpme_grid_z": int(nz_lj)
        })
    return settings


def _openmm_nonbonded_particle_records(force) -> list:
    records = []
    for index in range(force.getNumParticles()):
        charge, sigma, epsilon = force.getParticleParameters(index)
        records.append(((index,), (
            _openmm_float(charge, openmm_unit.elementary_charge),
            _openmm_float(sigma, openmm_unit.nanometer),
            _openmm_float(epsilon, openmm_unit.kilojoule_per_mole)
        )))
    return records


def _openmm_nonbonded_exception_records(force) -> list:
    records = []
    charge_product_unit = openmm_unit.elementary_charge**2
    for index in range(force.getNumExceptions()):
        p1, p2, charge_product, sigma, epsilon = force.getExceptionParameters(index)
        key = tuple(sorted((int(p1), int(p2))))
        records.append((key, (
            _openmm_float(charge_product, charge_product_unit),
            _openmm_float(sigma, openmm_unit.nanometer),
            _openmm_float(epsilon, openmm_unit.kilojoule_per_mole)
        )))
    return records


def _openmm_nonbonded_global_parameter_records(force) -> list:
    return [
        ((str(force.getGlobalParameterName(index)),),
         (float(force.getGlobalParameterDefaultValue(index)),))
        for index in range(force.getNumGlobalParameters())
    ]


def _openmm_nonbonded_particle_offset_records(force) -> list:
    records = []
    for index in range(force.getNumParticleParameterOffsets()):
        parameter, particle, charge_scale, sigma_scale, epsilon_scale = (
            force.getParticleParameterOffset(index)
        )
        records.append(((str(parameter), int(particle)), (
            _openmm_float(charge_scale, openmm_unit.elementary_charge),
            _openmm_float(sigma_scale, openmm_unit.nanometer),
            _openmm_float(epsilon_scale, openmm_unit.kilojoule_per_mole)
        )))
    return records


def _openmm_nonbonded_exception_offset_records(force) -> list:
    records = []
    charge_product_unit = openmm_unit.elementary_charge**2
    for index in range(force.getNumExceptionParameterOffsets()):
        parameter, exception_index, charge_scale, sigma_scale, epsilon_scale = (
            force.getExceptionParameterOffset(index)
        )
        p1, p2, _, _, _ = force.getExceptionParameters(int(exception_index))
        pair = tuple(sorted((int(p1), int(p2))))
        records.append(((str(parameter), *pair), (
            _openmm_float(charge_scale, charge_product_unit),
            _openmm_float(sigma_scale, openmm_unit.nanometer),
            _openmm_float(epsilon_scale, openmm_unit.kilojoule_per_mole)
        )))
    return records


def _openmm_atom_annotation(topology, particle_index: int) -> dict | None:
    """Return a JSON-safe chemical identity for one OpenMM particle index."""
    if topology is None:
        return None
    try:
        atom = topology.atom(int(particle_index))
    except Exception:
        return None

    residue = atom.residue
    chain_index = int(residue.chain.index)
    chain_label = chr(ord("A") + chain_index) if 0 <= chain_index < 26 else str(chain_index)
    residue_label = f"{chain_label}:{residue.name}{residue.resSeq}"
    return {
        "particle_index": int(particle_index),
        "atom_name": str(atom.name),
        "residue_name": str(residue.name),
        "residue_sequence": int(residue.resSeq),
        "chain_index": chain_index,
        "residue_label": residue_label,
        "label": f"{residue_label}:{atom.name}",
    }


def _openmm_key_annotation(key, topology, atom_key_count: int | None) -> str | None:
    if topology is None or atom_key_count is None:
        return None
    labels = []
    for particle_index in list(key)[:atom_key_count]:
        annotation = _openmm_atom_annotation(topology, particle_index)
        labels.append(annotation["label"] if annotation else f"particle {particle_index}")
    suffix = list(key)[atom_key_count:]
    text = " / ".join(labels)
    if suffix:
        text += f" ; key-extra={suffix}"
    return text


def _openmm_finalize_delta_groups(groups: dict) -> list[dict]:
    result = []
    for name, values in groups.items():
        count = values["count"]
        result.append({
            "name": name,
            "count": count,
            "mean_signed_delta_e": values["sum_signed"] / count if count else 0.0,
            "mean_absolute_delta_e": values["sum_absolute"] / count if count else 0.0,
            "rms_delta_e": float(np.sqrt(values["sum_squared"] / count)) if count else 0.0,
            "max_absolute_delta_e": values["max_absolute"],
        })
    result.sort(key=lambda item: (item["max_absolute_delta_e"], item["count"]), reverse=True)
    return result


def _openmm_nonbonded_charge_residue_summary(records1: list, records2: list,
                                               topology1, topology2) -> dict | None:
    """Localize per-particle charge disagreements by residue type and instance."""
    if topology1 is None or topology2 is None:
        return None

    by_type = defaultdict(lambda: {
        "count": 0, "sum_signed": 0.0, "sum_absolute": 0.0,
        "sum_squared": 0.0, "max_absolute": 0.0,
    })
    by_instance = defaultdict(lambda: {
        "count": 0, "sum_signed": 0.0, "sum_absolute": 0.0,
        "sum_squared": 0.0, "max_absolute": 0.0,
    })
    individual = []

    for (key1, values1), (key2, values2) in zip(records1, records2):
        if tuple(key1) != tuple(key2):
            continue
        q1 = float(values1[0])
        q2 = float(values2[0])
        if _openmm_numbers_match(q1, q2):
            continue

        particle_index = int(key1[0])
        annotation1 = _openmm_atom_annotation(topology1, particle_index)
        annotation2 = _openmm_atom_annotation(topology2, particle_index)
        type1 = annotation1["residue_name"] if annotation1 else "<unknown>"
        type2 = annotation2["residue_name"] if annotation2 else "<unknown>"
        residue_type = type1 if type1 == type2 else f"{type1}/{type2}"
        label1 = annotation1["residue_label"] if annotation1 else f"particle {particle_index}"
        label2 = annotation2["residue_label"] if annotation2 else f"particle {particle_index}"
        residue_instance = label1 if label1 == label2 else f"{label1}/{label2}"

        delta = q1 - q2
        absolute = abs(delta)
        for groups, group_name in ((by_type, residue_type), (by_instance, residue_instance)):
            acc = groups[group_name]
            acc["count"] += 1
            acc["sum_signed"] += delta
            acc["sum_absolute"] += absolute
            acc["sum_squared"] += delta * delta
            acc["max_absolute"] = max(acc["max_absolute"], absolute)

        individual.append({
            "particle_index": particle_index,
            "input1": annotation1,
            "input2": annotation2,
            "charge_input1_e": q1,
            "charge_input2_e": q2,
            "delta_input1_minus_input2_e": delta,
            "absolute_difference_e": absolute,
        })

    individual.sort(key=lambda item: item["absolute_difference_e"], reverse=True)
    top_n = CONFIG["openmm_xml_validation"]["top_n_differences"]
    return {
        "mismatch_count": len(individual),
        "by_residue_type": _openmm_finalize_delta_groups(by_type),
        "by_residue_instance": _openmm_finalize_delta_groups(by_instance),
        "top_particle_mismatches": individual[:top_n],
    }


def _openmm_constraint_geometry(system, state, topology) -> dict:
    """Measure how closely a serialized State satisfies its System constraints."""
    positions = np.asarray(
        state.getPositions(asNumpy=True).value_in_unit(openmm_unit.nanometer), dtype=float
    )
    box = np.asarray(
        state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(openmm_unit.nanometer), dtype=float
    )
    inverse_box = None
    try:
        if abs(float(np.linalg.det(box))) > 1.0e-15:
            inverse_box = np.linalg.inv(box)
    except np.linalg.LinAlgError:
        inverse_box = None

    buckets = {"all": [], "water": [], "nonwater": []}
    top_n = CONFIG["openmm_xml_validation"]["top_n_differences"]
    top_heap = []
    for constraint_index in range(system.getNumConstraints()):
        p1, p2, target = system.getConstraintParameters(constraint_index)
        p1 = int(p1)
        p2 = int(p2)
        target_nm = _openmm_float(target, openmm_unit.nanometer)
        displacement = positions[p2] - positions[p1]
        if inverse_box is not None:
            fractional = displacement @ inverse_box
            fractional -= np.round(fractional)
            displacement = fractional @ box
        actual_nm = float(np.linalg.norm(displacement))
        delta_nm = actual_nm - target_nm
        absolute_nm = abs(delta_nm)

        is_water = False
        annotation1 = _openmm_atom_annotation(topology, p1)
        annotation2 = _openmm_atom_annotation(topology, p2)
        if annotation1 and annotation2:
            is_water = (
                annotation1["residue_label"] == annotation2["residue_label"]
                and annotation1["residue_name"].upper() in {"HOH", "WAT"}
            )

        buckets["all"].append(delta_nm)
        buckets["water" if is_water else "nonwater"].append(delta_nm)
        record = {
            "constraint_index": int(constraint_index),
            "particle1": p1,
            "particle2": p2,
            "input_atom1": annotation1,
            "input_atom2": annotation2,
            "is_water": bool(is_water),
            "target_nm": target_nm,
            "actual_nm": actual_nm,
            "delta_actual_minus_target_nm": delta_nm,
            "absolute_violation_nm": absolute_nm,
        }
        heap_item = (absolute_nm, int(constraint_index), record)
        if len(top_heap) < top_n:
            heapq.heappush(top_heap, heap_item)
        elif absolute_nm > top_heap[0][0]:
            heapq.heapreplace(top_heap, heap_item)

    thresholds = CONFIG["openmm_xml_validation"]["constraint_geometry_thresholds_nm"]
    summaries = {}
    for name, values in buckets.items():
        array = np.asarray(values, dtype=float)
        absolute = np.abs(array)
        summaries[name] = {
            "count": int(array.size),
            "mean_signed_violation_nm": float(array.mean()) if array.size else 0.0,
            "mean_absolute_violation_nm": float(absolute.mean()) if array.size else 0.0,
            "rms_violation_nm": float(np.sqrt(np.mean(array * array))) if array.size else 0.0,
            "max_absolute_violation_nm": float(absolute.max()) if array.size else 0.0,
            "counts_above_threshold_nm": {
                f"{threshold:.1e}": int(np.count_nonzero(absolute > threshold))
                for threshold in thresholds
            },
        }

    top_records = [item[2] for item in sorted(top_heap, reverse=True)]
    return {
        "minimum_image_used": inverse_box is not None,
        "summaries": summaries,
        "top_violations": top_records,
    }


def _compare_openmm_force(force1, force2, constraint_pairs1=None, constraint_pairs2=None,
                          topology1=None, topology2=None) -> dict:
    force_type = type(force1).__name__
    result = {
        "type": force_type,
        "metadata": _openmm_force_metadata(force1, force2),
        "physics_differences": 0,
        "metadata_differences": 0,
        "fully_compared": True
    }
    result["metadata_differences"] = _openmm_force_metadata_difference_count(result["metadata"])

    pbc = result["metadata"].get("uses_periodic_boundary_conditions")
    if pbc is not None and not pbc["match"]:
        result["physics_differences"] += 1

    if isinstance(force1, mm.HarmonicBondForce) and isinstance(force2, mm.HarmonicBondForce):
        bond_records1 = _openmm_harmonic_bond_records(force1)
        bond_records2 = _openmm_harmonic_bond_records(force2)
        terms = _openmm_compare_keyed_records(
            bond_records1,
            bond_records2,
            ("length_nm", "k_kj_per_mol_nm2")
        )
        result["terms"] = terms

        constraint_pairs1 = set() if constraint_pairs1 is None else set(constraint_pairs1)
        constraint_pairs2 = set() if constraint_pairs2 is None else set(constraint_pairs2)
        bond_keys1 = [tuple(key) for key, _ in bond_records1]
        bond_keys2 = [tuple(key) for key, _ in bond_records2]
        unique_keys1 = set(bond_keys1)
        unique_keys2 = set(bond_keys2)
        only1 = unique_keys1 - unique_keys2
        only2 = unique_keys2 - unique_keys1
        result["constraint_overlap"] = {
            "harmonic_terms_on_constrained_pairs_input1": sum(
                1 for key in bond_keys1 if key in constraint_pairs1
            ),
            "harmonic_terms_on_constrained_pairs_input2": sum(
                1 for key in bond_keys2 if key in constraint_pairs2
            ),
            "input1_only_terms_on_constrained_pairs": sum(
                1 for key in only1 if key in constraint_pairs1
            ),
            "input1_only_terms_not_constrained": sum(
                1 for key in only1 if key not in constraint_pairs1
            ),
            "input2_only_terms_on_constrained_pairs": sum(
                1 for key in only2 if key in constraint_pairs2
            ),
            "input2_only_terms_not_constrained": sum(
                1 for key in only2 if key not in constraint_pairs2
            ),
        }
        if not terms["match"]:
            result["physics_differences"] += 1

    elif isinstance(force1, mm.HarmonicAngleForce) and isinstance(force2, mm.HarmonicAngleForce):
        terms = _openmm_compare_keyed_records(
            _openmm_harmonic_angle_records(force1),
            _openmm_harmonic_angle_records(force2),
            ("theta_rad", "k_kj_per_mol_rad2")
        )
        result["terms"] = terms
        if not terms["match"]:
            result["physics_differences"] += 1

    elif isinstance(force1, mm.PeriodicTorsionForce) and isinstance(force2, mm.PeriodicTorsionForce):
        terms = _openmm_compare_keyed_records(
            _openmm_periodic_torsion_records(force1),
            _openmm_periodic_torsion_records(force2),
            ("phase_rad", "k_kj_per_mol")
        )
        result["terms"] = terms
        if not terms["match"]:
            result["physics_differences"] += 1

    elif isinstance(force1, mm.NonbondedForce) and isinstance(force2, mm.NonbondedForce):
        settings = _openmm_compare_settings(
            _openmm_nonbonded_settings(force1), _openmm_nonbonded_settings(force2)
        )
        particle_records1 = _openmm_nonbonded_particle_records(force1)
        particle_records2 = _openmm_nonbonded_particle_records(force2)
        particles = _openmm_compare_keyed_records(
            particle_records1,
            particle_records2,
            ("charge_e", "sigma_nm", "epsilon_kj_per_mol")
        )
        exception_records1 = _openmm_nonbonded_exception_records(force1)
        exception_records2 = _openmm_nonbonded_exception_records(force2)
        exceptions = _openmm_compare_keyed_records(
            exception_records1,
            exception_records2,
            ("charge_product_e2", "sigma_nm", "epsilon_kj_per_mol")
        )
        globals_result = _openmm_compare_keyed_records(
            _openmm_nonbonded_global_parameter_records(force1),
            _openmm_nonbonded_global_parameter_records(force2),
            ("default_value",)
        )
        particle_offsets = _openmm_compare_keyed_records(
            _openmm_nonbonded_particle_offset_records(force1),
            _openmm_nonbonded_particle_offset_records(force2),
            ("charge_scale_e", "sigma_scale_nm", "epsilon_scale_kj_per_mol")
        )
        exception_offsets = _openmm_compare_keyed_records(
            _openmm_nonbonded_exception_offset_records(force1),
            _openmm_nonbonded_exception_offset_records(force2),
            ("charge_product_scale_e2", "sigma_scale_nm", "epsilon_scale_kj_per_mol")
        )
        charge1 = sum(values[0] for _, values in particle_records1)
        charge2 = sum(values[0] for _, values in particle_records2)
        result.update({
            "settings": settings,
            "particles": particles,
            "particle_charge_totals_e": {
                "input1": charge1,
                "input2": charge2,
                "delta_input1_minus_input2": charge1 - charge2,
            },
            "charge_residue_summary": _openmm_nonbonded_charge_residue_summary(
                particle_records1, particle_records2, topology1, topology2
            ),
            "exceptions": exceptions,
            "global_parameters": globals_result,
            "particle_parameter_offsets": particle_offsets,
            "exception_parameter_offsets": exception_offsets
        })
        for part in (settings, particles, exceptions, globals_result, particle_offsets, exception_offsets):
            if not part["match"]:
                result["physics_differences"] += 1

    elif isinstance(force1, mm.CMMotionRemover) and isinstance(force2, mm.CMMotionRemover):
        frequency = _openmm_scalar_comparison(
            int(force1.getFrequency()), int(force2.getFrequency())
        )
        result["frequency"] = frequency
        if not frequency["match"]:
            result["physics_differences"] += 1

    else:
        # Preserve visibility of an unexpected force type.  Raw serialization is
        # still useful as a strict fallback, but do not claim semantic equivalence
        # when we have not explicitly compared that Force API.
        result["fully_compared"] = False
        xml1 = mm.XmlSerializer.serialize(force1)
        xml2 = mm.XmlSerializer.serialize(force2)
        result["exact_force_xml_match"] = xml1 == xml2
        if xml1 != xml2:
            result["physics_differences"] += 1

    return result


def _compare_openmm_systems(system1, system2, topology1=None, topology2=None) -> dict:
    problems = []
    warnings = []
    result = {}

    n_particles1 = int(system1.getNumParticles())
    n_particles2 = int(system2.getNumParticles())
    result["particles"] = _openmm_scalar_comparison(n_particles1, n_particles2)
    if n_particles1 != n_particles2:
        problems.append(
            f"particle count differs: {n_particles1} vs {n_particles2}"
        )

    mass_records1 = [
        ((index,), (_openmm_float(system1.getParticleMass(index), openmm_unit.dalton),))
        for index in range(n_particles1)
    ]
    mass_records2 = [
        ((index,), (_openmm_float(system2.getParticleMass(index), openmm_unit.dalton),))
        for index in range(n_particles2)
    ]
    masses = _openmm_compare_keyed_records(mass_records1, mass_records2, ("mass_dalton",))
    result["particle_masses"] = masses
    mass_total1 = sum(values[0] for _, values in mass_records1)
    mass_total2 = sum(values[0] for _, values in mass_records2)
    result["particle_mass_totals_dalton"] = {
        "input1": mass_total1,
        "input2": mass_total2,
        "delta_input1_minus_input2": mass_total1 - mass_total2,
    }
    if not masses["match"]:
        problems.append(
            f"particle masses differ: {masses['parameter_mismatch_count']} parameter mismatches, "
            f"{masses['keys_only_input1_count']} input1-only, {masses['keys_only_input2_count']} input2-only"
        )

    constraint_records1 = []
    for index in range(system1.getNumConstraints()):
        p1, p2, distance = system1.getConstraintParameters(index)
        constraint_records1.append((
            tuple(sorted((int(p1), int(p2)))),
            (_openmm_float(distance, openmm_unit.nanometer),)
        ))
    constraint_records2 = []
    for index in range(system2.getNumConstraints()):
        p1, p2, distance = system2.getConstraintParameters(index)
        constraint_records2.append((
            tuple(sorted((int(p1), int(p2)))),
            (_openmm_float(distance, openmm_unit.nanometer),)
        ))
    constraints = _openmm_compare_keyed_records(
        constraint_records1, constraint_records2, ("distance_nm",)
    )
    result["constraints"] = constraints
    constraint_pairs1 = {tuple(key) for key, _ in constraint_records1}
    constraint_pairs2 = {tuple(key) for key, _ in constraint_records2}
    if not constraints["match"]:
        problems.append(
            f"constraints differ: {constraints['count_input1']} vs {constraints['count_input2']} terms; "
            f"{constraints['keys_only_input1_count']} input1-only, "
            f"{constraints['keys_only_input2_count']} input2-only, "
            f"{constraints['parameter_mismatch_count']} distance mismatches"
        )

    box1 = system1.getDefaultPeriodicBoxVectors()
    box2 = system2.getDefaultPeriodicBoxVectors()
    box_values1 = np.asarray([
        [_openmm_float(component, openmm_unit.nanometer) for component in vector]
        for vector in box1
    ], dtype=float)
    box_values2 = np.asarray([
        [_openmm_float(component, openmm_unit.nanometer) for component in vector]
        for vector in box2
    ], dtype=float)
    box_delta = box_values1 - box_values2
    box_match = all(
        _openmm_numbers_match(float(a), float(b))
        for a, b in zip(box_values1.ravel(), box_values2.ravel())
    )
    result["default_periodic_box_nm"] = {
        "input1": box_values1.tolist(),
        "input2": box_values2.tolist(),
        "delta_input1_minus_input2": box_delta.tolist(),
        "max_absolute_difference_nm": float(np.max(np.abs(box_delta))),
        "match": box_match
    }
    if not box_match:
        problems.append(
            "default periodic box differs: max |delta|="
            f"{result['default_periodic_box_nm']['max_absolute_difference_nm']:.6g} nm"
        )

    virtual_sites1 = [index for index in range(n_particles1) if system1.isVirtualSite(index)]
    virtual_sites2 = [index for index in range(n_particles2) if system2.isVirtualSite(index)]
    result["virtual_sites"] = {
        "indices_input1": virtual_sites1,
        "indices_input2": virtual_sites2,
        "match": virtual_sites1 == virtual_sites2,
        "fully_compared": not virtual_sites1 and not virtual_sites2
    }
    if virtual_sites1 != virtual_sites2:
        problems.append(
            f"virtual-site particle indices differ: {len(virtual_sites1)} vs {len(virtual_sites2)}"
        )
    elif virtual_sites1:
        warnings.append(
            "virtual sites are present; indices/types should be inspected because their internal parameters "
            "are not yet semantically expanded by this comparator"
        )

    inventory1 = _openmm_force_inventory(system1)
    inventory2 = _openmm_force_inventory(system2)
    inventory_result = {}
    for force_type in sorted(set(inventory1) | set(inventory2)):
        count1 = len(inventory1.get(force_type, []))
        count2 = len(inventory2.get(force_type, []))
        inventory_result[force_type] = {
            "input1": count1,
            "input2": count2,
            "match": count1 == count2
        }
        if count1 != count2:
            problems.append(
                f"force inventory differs for {force_type}: {count1} vs {count2}"
            )
    result["force_inventory"] = inventory_result

    forces_result = {}
    metadata_difference_count = 0
    fully_compared = True
    for force_type in sorted(set(inventory1) | set(inventory2)):
        forces_result[force_type] = []
        list1 = inventory1.get(force_type, [])
        list2 = inventory2.get(force_type, [])
        for occurrence, ((index1, force1), (index2, force2)) in enumerate(zip(list1, list2)):
            comparison = _compare_openmm_force(
                force1, force2,
                constraint_pairs1=constraint_pairs1, constraint_pairs2=constraint_pairs2,
                topology1=topology1, topology2=topology2
            )
            comparison["force_index_input1"] = int(index1)
            comparison["force_index_input2"] = int(index2)
            forces_result[force_type].append(comparison)
            metadata_difference_count += comparison["metadata_differences"]
            fully_compared = fully_compared and comparison["fully_compared"]
            if comparison["physics_differences"]:
                problems.append(
                    f"{force_type}[{occurrence}] differs in "
                    f"{comparison['physics_differences']} physics-affecting section(s)"
                )
            if not comparison["fully_compared"]:
                warnings.append(
                    f"{force_type}[{occurrence}] used exact Force XML as a fallback; "
                    "no force-specific semantic comparator is implemented"
                )

    result["forces"] = forces_result
    result["metadata_difference_count"] = metadata_difference_count
    result["fully_compared"] = fully_compared and result["virtual_sites"]["fully_compared"]
    result["problems"] = problems
    result["warnings"] = warnings
    result["physics_match"] = len(problems) == 0 and result["fully_compared"]
    return result


def _compare_openmm_integrators(integrator1, integrator2) -> dict:
    problems = []
    type1 = type(integrator1).__name__
    type2 = type(integrator2).__name__
    result = {
        "type": {"input1": type1, "input2": type2, "match": type1 == type2},
        "settings": {}
    }
    if type1 != type2:
        problems.append(f"integrator type differs: {type1} vs {type2}")

    settings1 = {
        "step_size_ps": _openmm_float(integrator1.getStepSize(), openmm_unit.picosecond),
        "constraint_tolerance": float(integrator1.getConstraintTolerance()),
        "integration_force_groups": int(integrator1.getIntegrationForceGroups())
    }
    settings2 = {
        "step_size_ps": _openmm_float(integrator2.getStepSize(), openmm_unit.picosecond),
        "constraint_tolerance": float(integrator2.getConstraintTolerance()),
        "integration_force_groups": int(integrator2.getIntegrationForceGroups())
    }

    for name, getter, target_unit in (
        ("temperature_k", "getTemperature", openmm_unit.kelvin),
        ("friction_per_ps", "getFriction", openmm_unit.picosecond**-1),
        ("random_number_seed", "getRandomNumberSeed", None),
    ):
        if hasattr(integrator1, getter) and hasattr(integrator2, getter):
            value1 = getattr(integrator1, getter)()
            value2 = getattr(integrator2, getter)()
            if name == "random_number_seed":
                settings1[name] = int(value1)
                settings2[name] = int(value2)
            else:
                settings1[name] = _openmm_float(value1, target_unit)
                settings2[name] = _openmm_float(value2, target_unit)

    settings = _openmm_compare_settings(settings1, settings2)
    result["settings"] = settings
    if not settings["match"]:
        problems.append(f"integrator settings differ in {settings['n_different']} field(s)")
    result["problems"] = problems
    result["match"] = len(problems) == 0
    return result


def _compare_openmm_states(state1, state2) -> dict:
    problems = []
    cfg = CONFIG["openmm_xml_validation"]
    result = {}

    positions1 = np.asarray(
        state1.getPositions(asNumpy=True).value_in_unit(openmm_unit.nanometer), dtype=float
    )
    positions2 = np.asarray(
        state2.getPositions(asNumpy=True).value_in_unit(openmm_unit.nanometer), dtype=float
    )
    position_summary = {
        "shape_input1": list(positions1.shape),
        "shape_input2": list(positions2.shape),
        "match": False
    }
    if positions1.shape == positions2.shape:
        delta = positions1 - positions2
        per_atom = np.sqrt(np.sum(delta * delta, axis=1))
        bad = per_atom > cfg["position_tolerance_nm"]
        order = np.argsort(per_atom)[::-1]
        top_n = cfg["top_n_differences"]
        centroid1 = positions1.mean(axis=0)
        centroid2 = positions2.mean(axis=0)
        centered1 = positions1 - centroid1
        centered2 = positions2 - centroid2
        centered_delta = centered1 - centered2
        centered_per_atom = np.sqrt(np.sum(centered_delta * centered_delta, axis=1))

        covariance = centered1.T @ centered2
        u_mat, _, v_t = np.linalg.svd(covariance)
        rotation = u_mat @ v_t
        if np.linalg.det(rotation) < 0.0:
            u_mat[:, -1] *= -1.0
            rotation = u_mat @ v_t
        aligned1 = centered1 @ rotation
        aligned_delta = aligned1 - centered2
        aligned_per_atom = np.sqrt(np.sum(aligned_delta * aligned_delta, axis=1))

        position_summary.update({
            "n_atoms": int(positions1.shape[0]),
            "position_tolerance_nm": cfg["position_tolerance_nm"],
            "n_atoms_above_tolerance": int(np.count_nonzero(bad)),
            "max_atom_displacement_nm": float(np.max(per_atom)) if len(per_atom) else 0.0,
            "median_atom_displacement_nm": float(np.median(per_atom)) if len(per_atom) else 0.0,
            "rms_atom_displacement_nm": float(np.sqrt(np.mean(per_atom * per_atom))) if len(per_atom) else 0.0,
            "max_absolute_coordinate_difference_nm": float(np.max(np.abs(delta))) if delta.size else 0.0,
            "centroid_delta_input1_minus_input2_nm": [float(v) for v in (centroid1 - centroid2)],
            "centered_rms_atom_displacement_nm": (
                float(np.sqrt(np.mean(centered_per_atom * centered_per_atom)))
                if len(centered_per_atom) else 0.0
            ),
            "centered_max_atom_displacement_nm": (
                float(np.max(centered_per_atom)) if len(centered_per_atom) else 0.0
            ),
            "kabsch_rms_atom_displacement_nm": (
                float(np.sqrt(np.mean(aligned_per_atom * aligned_per_atom)))
                if len(aligned_per_atom) else 0.0
            ),
            "kabsch_max_atom_displacement_nm": (
                float(np.max(aligned_per_atom)) if len(aligned_per_atom) else 0.0
            ),
            "top_atom_displacements": [
                {
                    "particle_index": int(index),
                    "displacement_nm": float(per_atom[index]),
                    "delta_xyz_nm": [float(v) for v in delta[index]]
                }
                for index in order[:top_n]
                if per_atom[index] > cfg["position_tolerance_nm"]
            ],
            "match": not np.any(bad)
        })
    result["positions"] = position_summary
    if not position_summary["match"]:
        problems.append(
            "initial positions differ: "
            f"{position_summary.get('n_atoms_above_tolerance', 'shape mismatch')} atoms above "
            f"{cfg['position_tolerance_nm']:.1e} nm tolerance"
        )

    box1 = np.asarray(
        state1.getPeriodicBoxVectors(asNumpy=True).value_in_unit(openmm_unit.nanometer), dtype=float
    )
    box2 = np.asarray(
        state2.getPeriodicBoxVectors(asNumpy=True).value_in_unit(openmm_unit.nanometer), dtype=float
    )
    box_delta = box1 - box2
    box_match = all(
        _openmm_numbers_match(float(a), float(b))
        for a, b in zip(box1.ravel(), box2.ravel())
    )
    result["periodic_box_nm"] = {
        "input1": box1.tolist(),
        "input2": box2.tolist(),
        "delta_input1_minus_input2": box_delta.tolist(),
        "max_absolute_difference_nm": float(np.max(np.abs(box_delta))),
        "match": box_match
    }
    if not box_match:
        problems.append(
            "initial State periodic box differs: max |delta|="
            f"{result['periodic_box_nm']['max_absolute_difference_nm']:.6g} nm"
        )

    potential1 = _openmm_float(state1.getPotentialEnergy(), openmm_unit.kilojoule_per_mole)
    potential2 = _openmm_float(state2.getPotentialEnergy(), openmm_unit.kilojoule_per_mole)
    result["potential_energy_kj_per_mol"] = _openmm_scalar_comparison(
        potential1, potential2,
        absolute_tolerance=cfg["energy_tolerance_kj_per_mol"],
        relative_tolerance=0.0
    )
    if not result["potential_energy_kj_per_mol"]["match"]:
        problems.append(
            "initial State potential energy differs by "
            f"{result['potential_energy_kj_per_mol']['absolute_difference']:.6g} kJ/mol"
        )

    kinetic1 = _openmm_float(state1.getKineticEnergy(), openmm_unit.kilojoule_per_mole)
    kinetic2 = _openmm_float(state2.getKineticEnergy(), openmm_unit.kilojoule_per_mole)
    result["kinetic_energy_kj_per_mol"] = _openmm_scalar_comparison(
        kinetic1, kinetic2,
        absolute_tolerance=cfg["energy_tolerance_kj_per_mol"],
        relative_tolerance=0.0
    )
    result["kinetic_energy_kj_per_mol"]["diagnostic_only"] = True

    result["time_ps"] = _openmm_scalar_comparison(
        _openmm_float(state1.getTime(), openmm_unit.picosecond),
        _openmm_float(state2.getTime(), openmm_unit.picosecond)
    )
    if not result["time_ps"]["match"]:
        problems.append("initial State time differs")

    if hasattr(state1, "getStepCount") and hasattr(state2, "getStepCount"):
        result["step_count"] = _openmm_scalar_comparison(
            int(state1.getStepCount()), int(state2.getStepCount())
        )
        if not result["step_count"]["match"]:
            problems.append("initial State step count differs")

    try:
        parameters1 = {str(key): float(value) for key, value in state1.getParameters().items()}
        parameters2 = {str(key): float(value) for key, value in state2.getParameters().items()}
        result["context_parameters"] = _openmm_compare_settings(parameters1, parameters2)
        if not result["context_parameters"]["match"]:
            problems.append(
                f"initial State Context parameters differ in "
                f"{result['context_parameters']['n_different']} field(s)"
            )
    except Exception as exc:
        result["context_parameters"] = {"available": False, "error": str(exc)}

    result["problems"] = problems
    result["match"] = len(problems) == 0
    return result


def compare_serialized_openmm_xml(out1: Path | None, out2: Path | None,
                                  label1: str, label2: str,
                                  topology1=None, topology2=None):
    """Discover and semantically compare OpenMM XmlSerializer files beside MD logs."""
    paths1 = _discover_openmm_serialized_paths(out1)
    paths2 = _discover_openmm_serialized_paths(out2)
    any_xml = any(
        path is not None and Path(path).is_file()
        for path in list(paths1.values()) + list(paths2.values())
    )
    if not any_xml:
        return None

    result = {
        "comparison_version": 3,
        "input1_label": label1,
        "input2_label": label2,
        "tolerances": dict(CONFIG["openmm_xml_validation"]),
        "files": {"input1": {}, "input2": {}},
        "problems": [],
        "warnings": []
    }
    for key in _OPENMM_SERIALIZED_FILENAMES:
        result["files"]["input1"][key] = _openmm_path_metadata(paths1[key])
        result["files"]["input2"][key] = _openmm_path_metadata(paths2[key])

    for key in _OPENMM_SERIALIZED_FILENAMES:
        exists1 = result["files"]["input1"][key]["exists"]
        exists2 = result["files"]["input2"][key]["exists"]
        if exists1 != exists2:
            missing_label = label2 if exists1 else label1
            result["problems"].append(f"{key} XML is missing for {missing_label}")

    if mm is None:
        result["warnings"].append(
            "OpenMM Python module is unavailable; serialized XML files were discovered but can not be deserialized."
        )
        result["fully_compared"] = False
        return result

    comparisons = (
        ("system", lambda obj1, obj2: _compare_openmm_systems(
            obj1, obj2, topology1=topology1, topology2=topology2
        )),
        ("integrator", _compare_openmm_integrators),
        ("initial_state", _compare_openmm_states),
    )
    loaded_objects = {}
    for key, comparator in comparisons:
        path1 = paths1[key]
        path2 = paths2[key]
        exists1 = path1 is not None and Path(path1).is_file()
        exists2 = path2 is not None and Path(path2).is_file()
        if not exists1 or not exists2:
            continue

        result["files"]["input1"][key]["exact_bytes_match"] = (
            Path(path1).read_bytes() == Path(path2).read_bytes()
        )
        result["files"]["input2"][key]["exact_bytes_match"] = (
            result["files"]["input1"][key]["exact_bytes_match"]
        )

        try:
            object1 = _load_openmm_serialized_object(Path(path1))
            object2 = _load_openmm_serialized_object(Path(path2))
            loaded_objects[key] = (object1, object2)
            comparison = comparator(object1, object2)
            result[key] = comparison
            for problem in comparison.get("problems", []):
                result["problems"].append(f"{key}: {problem}")
            for warning in comparison.get("warnings", []):
                result["warnings"].append(f"{key}: {warning}")
        except Exception as exc:
            result["problems"].append(
                f"failed to compare {key} XML: {type(exc).__name__}: {exc}"
            )
            result[key] = {"comparison_error": f"{type(exc).__name__}: {exc}"}

    if "system" in loaded_objects and "initial_state" in loaded_objects:
        system1, system2 = loaded_objects["system"]
        state1, state2 = loaded_objects["initial_state"]
        try:
            result["initial_constraint_geometry"] = {
                "input1": _openmm_constraint_geometry(system1, state1, topology1),
                "input2": _openmm_constraint_geometry(system2, state2, topology2),
            }
        except Exception as exc:
            result["warnings"].append(
                f"initial constraint geometry diagnostic failed: {type(exc).__name__}: {exc}"
            )

    result["fully_compared"] = all(
        key in result
        and "comparison_error" not in result[key]
        and result[key].get("fully_compared", True)
        for key, _ in comparisons
        if result["files"]["input1"][key]["exists"]
        or result["files"]["input2"][key]["exists"]
    )
    result["match"] = len(result["problems"]) == 0 and result["fully_compared"]
    return result


def _print_openmm_record_summary(title: str, result: dict) -> None:
    print(
        f"    {title:<28}: {result['count_input1']} vs {result['count_input2']}  "
        f"only1={result['keys_only_input1_count']} only2={result['keys_only_input2_count']} "
        f"param_diff={result['parameter_mismatch_count']}"
    )


def _print_openmm_field_statistics(title: str, result: dict) -> None:
    for field, stats in result.get("field_statistics", {}).items():
        if stats["mismatch_count"] == 0:
            continue
        extra = ""
        if field == "theta_rad":
            extra = (
                f"; max={np.degrees(stats['max_absolute_difference']):.4g} deg"
                f"; rms={np.degrees(stats['rms_difference']):.4g} deg"
            )
        print(
            f"      {title} {field}: mismatches={stats['mismatch_count']}/"
            f"{stats['compared_count']} mean_abs={stats['mean_absolute_difference']:.6g} "
            f"rms={stats['rms_difference']:.6g} max={stats['max_absolute_difference']:.6g}{extra}"
        )


def _print_openmm_top_record_differences(title: str, result: dict, limit: int = 3,
                                         topology1=None, topology2=None,
                                         atom_key_count: int | None = None) -> None:
    shown = 0
    for side, key_name, topology in (
        ("input1-only", "keys_only_input1", topology1),
        ("input2-only", "keys_only_input2", topology2),
    ):
        for key in result.get(key_name, []):
            if shown >= limit:
                return
            annotation = _openmm_key_annotation(key, topology, atom_key_count)
            suffix = f"  [{annotation}]" if annotation else ""
            print(f"      {title} {side} key={key}{suffix}")
            shown += 1

    for mismatch in result.get("top_parameter_mismatches", []):
        if shown >= limit:
            return
        changed = []
        for field, values in mismatch["fields"].items():
            if values["match"]:
                continue
            changed.append(
                f"{field}: {values['input1']:.8g} vs {values['input2']:.8g} "
                f"(delta={values['delta_input1_minus_input2']:+.3e})"
            )
        annotation1 = _openmm_key_annotation(mismatch["key"], topology1, atom_key_count)
        annotation2 = _openmm_key_annotation(mismatch["key"], topology2, atom_key_count)
        annotations = []
        if annotation1:
            annotations.append(f"input1={annotation1}")
        if annotation2 and annotation2 != annotation1:
            annotations.append(f"input2={annotation2}")
        suffix = f"  [{' ; '.join(annotations)}]" if annotations else ""
        print(f"      {title} key={mismatch['key']}{suffix}  " + "; ".join(changed))
        shown += 1


def print_serialized_openmm_xml_comparison(label1: str, label2: str, comparison: dict,
                                              topology1=None, topology2=None) -> None:
    if not comparison:
        return

    print("\nSerialized OpenMM XML Validation")
    print("=" * 72)
    for key, filename in _OPENMM_SERIALIZED_FILENAMES.items():
        file1 = comparison["files"]["input1"][key]
        file2 = comparison["files"]["input2"][key]
        status1 = file1["path"] if file1["exists"] else "MISSING"
        status2 = file2["path"] if file2["exists"] else "MISSING"
        print(f"  {filename}")
        print(f"    {label1:<20}: {status1}")
        print(f"    {label2:<20}: {status2}")

    system = comparison.get("system")
    if system and "comparison_error" not in system:
        print("\n  System semantic comparison")
        particle_status = "MATCH" if system["particles"]["match"] else "DIFFER"
        print(
            f"    particles                    : {particle_status:<6}  "
            f"{system['particles']['input1']} vs {system['particles']['input2']}"
        )
        _print_openmm_record_summary("particle masses", system["particle_masses"])
        _print_openmm_field_statistics("particle masses", system["particle_masses"])
        print(
            f"      total mass delta: "
            f"{system['particle_mass_totals_dalton']['delta_input1_minus_input2']:+.6g} Da"
        )
        if not system["particle_masses"]["match"]:
            _print_openmm_top_record_differences(
                "particle masses", system["particle_masses"],
                topology1=topology1, topology2=topology2, atom_key_count=1
            )
        _print_openmm_record_summary("constraints", system["constraints"])
        _print_openmm_field_statistics("constraints", system["constraints"])
        if not system["constraints"]["match"]:
            _print_openmm_top_record_differences(
                "constraints", system["constraints"],
                topology1=topology1, topology2=topology2, atom_key_count=2
            )
        print(
            f"    default periodic box         : max |delta|="
            f"{system['default_periodic_box_nm']['max_absolute_difference_nm']:.6g} nm"
        )
        print("    force inventory:")
        for force_type, values in system["force_inventory"].items():
            status = "MATCH" if values["match"] else "DIFFER"
            print(
                f"      {force_type:<24}: {status:<6}  "
                f"{values['input1']} vs {values['input2']}"
            )
        for force_type, force_comparisons in system["forces"].items():
            for occurrence, force_result in enumerate(force_comparisons):
                prefix = f"{force_type}[{occurrence}]"
                if "terms" in force_result:
                    _print_openmm_record_summary(prefix, force_result["terms"])
                    _print_openmm_field_statistics(prefix, force_result["terms"])
                    if force_type == "HarmonicBondForce" and "constraint_overlap" in force_result:
                        overlap = force_result["constraint_overlap"]
                        print(
                            f"      constrained harmonic terms: "
                            f"{overlap['harmonic_terms_on_constrained_pairs_input1']} vs "
                            f"{overlap['harmonic_terms_on_constrained_pairs_input2']}"
                        )
                        print(
                            f"      input1-only bonds: constrained="
                            f"{overlap['input1_only_terms_on_constrained_pairs']} "
                            f"unconstrained={overlap['input1_only_terms_not_constrained']}; "
                            f"input2-only: constrained="
                            f"{overlap['input2_only_terms_on_constrained_pairs']} "
                            f"unconstrained={overlap['input2_only_terms_not_constrained']}"
                        )
                    if not force_result["terms"]["match"]:
                        atom_key_count = {
                            "HarmonicBondForce": 2,
                            "HarmonicAngleForce": 3,
                            "PeriodicTorsionForce": 4,
                        }.get(force_type)
                        _print_openmm_top_record_differences(
                            prefix, force_result["terms"],
                            topology1=topology1, topology2=topology2,
                            atom_key_count=atom_key_count
                        )
                if force_type == "NonbondedForce":
                    print(
                        f"    {prefix + ' settings':<28}: "
                        f"different={force_result['settings']['n_different']}"
                    )
                    _print_openmm_record_summary(prefix + " particles", force_result["particles"])
                    _print_openmm_field_statistics(prefix + " particles", force_result["particles"])
                    charge_totals = force_result["particle_charge_totals_e"]
                    print(
                        f"      total particle charge (e): "
                        f"{charge_totals['input1']:+.9g} vs {charge_totals['input2']:+.9g} "
                        f"(delta={charge_totals['delta_input1_minus_input2']:+.6g})"
                    )
                    if not force_result["particles"]["match"]:
                        _print_openmm_top_record_differences(
                            prefix + " particles", force_result["particles"],
                            topology1=topology1, topology2=topology2, atom_key_count=1
                        )
                    charge_summary = force_result.get("charge_residue_summary")
                    if charge_summary and charge_summary["mismatch_count"]:
                        print("      charge mismatches by residue type (largest |delta q|):")
                        for entry in charge_summary["by_residue_type"][:5]:
                            print(
                                f"        {entry['name']:<20} n={entry['count']:<5d} "
                                f"mean_abs={entry['mean_absolute_delta_e']:.6g} e "
                                f"max={entry['max_absolute_delta_e']:.6g} e"
                            )
                        print("      largest residue-instance charge mismatches:")
                        for entry in charge_summary["by_residue_instance"][:3]:
                            print(
                                f"        {entry['name']:<28} n={entry['count']:<4d} "
                                f"max={entry['max_absolute_delta_e']:.6g} e"
                            )
                    _print_openmm_record_summary(prefix + " exceptions", force_result["exceptions"])
                    _print_openmm_field_statistics(prefix + " exceptions", force_result["exceptions"])
                    if not force_result["exceptions"]["match"]:
                        _print_openmm_top_record_differences(
                            prefix + " exceptions", force_result["exceptions"],
                            topology1=topology1, topology2=topology2, atom_key_count=2
                        )
        if system["metadata_difference_count"]:
            print(
                f"    force metadata differences   : {system['metadata_difference_count']} "
                "(group/name; reported separately from Hamiltonian differences)"
            )

    integrator = comparison.get("integrator")
    if integrator and "comparison_error" not in integrator:
        status = "MATCH" if integrator["match"] else "DIFFER"
        print(f"\n  Integrator                     : {status}")
        print(
            f"    type                         : {integrator['type']['input1']} vs "
            f"{integrator['type']['input2']}"
        )
        for field, values in integrator["settings"]["fields"].items():
            if not values["match"]:
                print(f"    {field:<28}: {values['input1']} vs {values['input2']}")

    state = comparison.get("initial_state")
    if state and "comparison_error" not in state:
        status = "MATCH" if state["match"] else "DIFFER"
        positions = state["positions"]
        print(f"\n  Initial State                  : {status}")
        if positions.get("n_atoms") is not None:
            print(
                f"    positions                    : {positions['n_atoms_above_tolerance']} atoms above "
                f"{positions['position_tolerance_nm']:.1e} nm; max displacement="
                f"{positions['max_atom_displacement_nm']:.6g} nm; rms="
                f"{positions['rms_atom_displacement_nm']:.6g} nm"
            )
            centroid_delta = positions["centroid_delta_input1_minus_input2_nm"]
            print(
                "    centroid delta (1-2) nm      : "
                + " ".join(f"{value:+.6g}" for value in centroid_delta)
            )
            print(
                f"    centered coordinate delta    : rms="
                f"{positions['centered_rms_atom_displacement_nm']:.6g} nm; max="
                f"{positions['centered_max_atom_displacement_nm']:.6g} nm"
            )
            print(
                f"    Kabsch-aligned delta          : rms="
                f"{positions['kabsch_rms_atom_displacement_nm']:.6g} nm; max="
                f"{positions['kabsch_max_atom_displacement_nm']:.6g} nm"
            )
        else:
            print(
                f"    positions                    : shape {positions['shape_input1']} vs "
                f"{positions['shape_input2']}"
            )
        print(
            f"    periodic box                 : max |delta|="
            f"{state['periodic_box_nm']['max_absolute_difference_nm']:.6g} nm"
        )
        print(
            f"    potential energy             : delta="
            f"{state['potential_energy_kj_per_mol']['delta_input1_minus_input2']:+.6g} kJ/mol"
        )
        print(
            f"    kinetic energy (diagnostic)  : delta="
            f"{state['kinetic_energy_kj_per_mol']['delta_input1_minus_input2']:+.6g} kJ/mol"
        )

    constraint_geometry = comparison.get("initial_constraint_geometry")
    if constraint_geometry:
        print("\n  Initial constraint-geometry diagnostic")
        for side_key, label in (("input1", label1), ("input2", label2)):
            diagnostic = constraint_geometry[side_key]
            print(f"    {label}:")
            for bucket_name in ("all", "water", "nonwater"):
                stats = diagnostic["summaries"][bucket_name]
                thresholds = stats["counts_above_threshold_nm"]
                print(
                    f"      {bucket_name:<9} n={stats['count']:<7d} "
                    f"rms={stats['rms_violation_nm']:.6g} nm "
                    f"max={stats['max_absolute_violation_nm']:.6g} nm "
                    f">1e-5nm={thresholds.get('1.0e-05', 0)}"
                )
            if diagnostic["top_violations"]:
                print("      largest violations:")
                for item in diagnostic["top_violations"][:3]:
                    atom1 = item["input_atom1"]
                    atom2 = item["input_atom2"]
                    label_atom1 = atom1["label"] if atom1 else str(item["particle1"])
                    label_atom2 = atom2["label"] if atom2 else str(item["particle2"])
                    print(
                        f"        {label_atom1} -- {label_atom2}: "
                        f"actual-target={item['delta_actual_minus_target_nm']:+.6g} nm "
                        f"(target={item['target_nm']:.8g} nm)"
                    )

    if comparison.get("problems"):
        print("\n  Serialized XML differences requiring review:")
        for problem in comparison["problems"]:
            print(f"    - {problem}")
    elif comparison.get("fully_compared"):
        print("\n  No semantic differences detected in the serialized OpenMM objects.")
    else:
        print("\n  No differences detected in compared fields, but the XML comparison is incomplete.")

    if comparison.get("warnings"):
        print("\n  Serialized XML comparison warnings:")
        for warning in comparison["warnings"]:
            print(f"    - {warning}")


# ──────────────────────────  MAIN  ───────────────────────────────────────────
def main()->None:
    ap=argparse.ArgumentParser()
    ap.add_argument("traj1"); ap.add_argument("top1")
    ap.add_argument("traj2"); ap.add_argument("top2")
    ap.add_argument("-o","--out-prefix", default="comparison", help="Use as a basename for the comparison (e.g. PDBid)")
    ap.add_argument("--label1", default="Complete Structure")
    ap.add_argument("--label2", default="Published Structure")
    ap.add_argument("--out1", type=Path, help="MD output/screenout for trajectory 1; plain text or gzip (optional)")
    ap.add_argument("--out2", type=Path, help="MD output/screenout for trajectory 2; plain text or gzip (optional)")
    ap.add_argument("--ligand", type=str, help="3-letter ligand residue name (e.g., CBN or LIG)")
    ap.add_argument("-q","--quiet",action="store_true")
    args=ap.parse_args()

    logging.basicConfig(level=logging.WARNING if args.quiet else logging.INFO,
                        format="%(levelname)s: %(message)s")

    tmp_files=[]
    try:
        traj1_path = decompress_if_gz(Path(args.traj1), tmp_files)
        top1_path  = decompress_if_gz(Path(args.top1), tmp_files)
        traj2_path = decompress_if_gz(Path(args.traj2), tmp_files)
        top2_path  = decompress_if_gz(Path(args.top2), tmp_files)
        
        t1 = load_traj(traj1_path, top1_path)
        t2 = load_traj(traj2_path, top2_path)

        # Standardise residue names
        standardise_names(t1); standardise_names(t2)

        # Slice protein only
        prot1=t1.atom_slice([a.index for a in t1.topology.atoms if a.residue.is_protein])
        prot2=t2.atom_slice([a.index for a in t2.topology.atoms if a.residue.is_protein])
        
        summary = {
            "label1": args.label1,
            "label2": args.label2,
            "out_prefix": args.out_prefix
        }
        summary["trajectory_validation"] = {}
        summary["generated_on"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
        summary["units"] = {}
        summary["units"]["temp"] = "K"
        summary["units"]["energy"] = "kJ/mol"
        summary["units"]["energy"] = "kJ/mol"
        summary["n_frames"] = {
            args.label1: t1.n_frames,
            args.label2: t2.n_frames
        }
        summary["n_atoms"] = {
            args.label1: t1.n_atoms,
            args.label2: t2.n_atoms
        }
        summary["n_residues"] = {
            args.label1: len([r for r in t1.topology.residues if r.is_protein]),
            args.label2: len([r for r in t2.topology.residues if r.is_protein])
        }
        
        # ───────── Snapshot Output (solute-centered PDBs) ─────────
        print("\nSaving representative snapshots (solute-centered)...")
        
        def save_centered_snapshots(traj: md.Trajectory, solute_sel: np.ndarray, label: str):
            # Compute solute center-of-geometry per frame
            xyz_solute = traj.xyz[:, solute_sel, :]
            com = xyz_solute.mean(axis=1, keepdims=True)  # (n_frames, 1, 3)
        
            # Apply COM-centering to all atoms
            centered_xyz = traj.xyz - com
        
            # Create new trajectory with solute atoms only
            traj_solute = traj.atom_slice(solute_sel)
            traj_solute.xyz = centered_xyz[:, solute_sel, :]
        
            # Select 10 evenly spaced frame indices
            frame_indices = np.linspace(0, traj_solute.n_frames - 1, 10, dtype=int)
            snapshots = traj_solute.slice(frame_indices)
        
            out_path = f"{args.out_prefix}_snapshots_{safe(label)}.pdb"
            snapshots.save_pdb(out_path)
            print(f"  {label:<20} → {out_path}")
        
        # Select solute atoms (e.g., all protein)
        solute1 = t1.topology.select("protein")
        solute2 = t2.topology.select("protein")
        
        save_centered_snapshots(t1, solute1, args.label1)
        save_centered_snapshots(t2, solute2, args.label2)

        
        # ───────── FULL RMSF ─────────
        ca1_all = prot1.topology.select("name CA")
        ca2_all = prot2.topology.select("name CA")
        rmsf1_full = compute_rmsf(prot1, ca1_all)
        rmsf2_full = compute_rmsf(prot2, ca2_all)
        ks_full = ks_2samp(rmsf1_full,rmsf2_full)
        logging.info("FULL RMSF | %s mean %.3f±%.3f | %s mean %.3f±%.3f | KS‑p %.3e",
                     args.label1,rmsf1_full.mean(),rmsf1_full.std(),
                     args.label2,rmsf2_full.mean(),rmsf2_full.std(),ks_full.pvalue)
        plt.figure()
        plt.plot(rmsf1_full,label=args.label1); plt.plot(rmsf2_full,label=args.label2)
        plt.xlabel("Residue index"); plt.ylabel("RMSF (Å)")
        plt.title(f"Full RMSF  KS p={ks_full.pvalue:.3e}"); plt.legend(); plt.tight_layout()
        plt.savefig(f"{args.out_prefix}_rmsf_full.pdf",dpi=CONFIG["plot"]["dpi"]); plt.close()

        # ───────── SEQUENCE ALIGNMENT & ALIGNED RMSF ─────────
        seq1,caL1 = seq_and_ca(prot1); seq2,caL2 = seq_and_ca(prot2)
        aln = PairwiseAligner().align(seq1,seq2)[0]
        m1,m2=[],[]; i1=i2=0
        for a,b in zip(aln.target,aln.query):
            if a!="-" and b!="-": m1.append(caL1[i1]); m2.append(caL2[i2])
            if a!="-": i1+=1
            if b!="-": i2+=1

        m1_arr = np.asarray(m1, dtype=int)
        m2_arr = np.asarray(m2, dtype=int)

        # Alignment-dependent analyses must be symmetric.  Previously prot2 was
        # superposed in place while prot1 was left in its original orientation,
        # which contaminated RMSF, PCA and DCCM comparisons.  Build independent
        # copies and align BOTH trajectories to the same reference: trajectory 1,
        # frame 0, using the sequence-matched Cα atoms.  The raw prot1/prot2
        # trajectories remain untouched for all pre-existing non-aligned analyses.
        aligned_prot1 = aligned_trajectory_copy(prot1, prot1, m1_arr, m1_arr, frame=0)
        aligned_prot2 = aligned_trajectory_copy(prot2, prot1, m2_arr, m1_arr, frame=0)

        rmsf1_aln=compute_rmsf(aligned_prot1,m1_arr)
        rmsf2_aln=compute_rmsf(aligned_prot2,m2_arr)
        ks_aln=ks_2samp(rmsf1_aln,rmsf2_aln)
        logging.info("ALIGNED RMSF | %s %.3f±%.3f | %s %.3f±%.3f | KS‑p %.3e",
                     args.label1,rmsf1_aln.mean(),rmsf1_aln.std(),
                     args.label2,rmsf2_aln.mean(),rmsf2_aln.std(),ks_aln.pvalue)
        plt.figure()
        plt.plot(rmsf1_aln,label=args.label1); plt.plot(rmsf2_aln,label=args.label2)
        plt.xlabel("Aligned residue index"); plt.ylabel("RMSF (Å)")
        plt.title(f"Aligned‑core RMSF  KS p={ks_aln.pvalue:.3e}")
        plt.legend(); plt.tight_layout()
        plt.savefig(f"{args.out_prefix}_rmsf_aligned.pdf",dpi=CONFIG["plot"]["dpi"]); plt.close()

        # Aligned residue lists (built earlier from m1 / m2)
        aligned_res1 = [prot1.topology.atom(ca).residue for ca in m1]
        aligned_res2 = [prot2.topology.atom(ca).residue for ca in m2]
        n_aligned    = len(aligned_res1)

        # ───────── TRAJECTORY VALIDATION: RMSF PROFILE SIMILARITY ─────────
        # KS remains above for backwards compatibility.  For a residue-by-residue
        # comparison, however, paired metrics answer the more useful question:
        # do the same aligned residues show the same relative flexibility?
        rmsf_profile_stats = paired_profile_stats(rmsf1_aln, rmsf2_aln)
        delta_rmsf = rmsf1_aln - rmsf2_aln
        top_n_rmsf = min(CONFIG["trajectory_validation"]["top_n_rmsf_differences"], n_aligned)
        top_rmsf_idx = np.argsort(np.abs(delta_rmsf))[::-1][:top_n_rmsf]

        print("\nTrajectory Validation — Aligned RMSF Profile")
        print("=" * 48)
        if rmsf_profile_stats["pearson_r"] is not None:
            print(f"  Pearson r                    : {rmsf_profile_stats['pearson_r']:.4f}")
            print(f"  Spearman rho                 : {rmsf_profile_stats['spearman_rho']:.4f}")
        else:
            print("  Pearson/Spearman             : undefined (constant profile)")
        print(f"  RMSF RMSE                    : {rmsf_profile_stats['rmse']:.3f} Å")
        print(f"  Mean |ΔRMSF|                 : {rmsf_profile_stats['mae']:.3f} Å")
        print(f"  Maximum |ΔRMSF|              : {rmsf_profile_stats['max_abs_difference']:.3f} Å")
        print(f"  Top {top_n_rmsf} residue differences:")
        for i in top_rmsf_idx:
            print(f"    {res_label(aligned_res1[i]):<12} / {res_label(aligned_res2[i]):<12} "
                  f"ΔRMSF={delta_rmsf[i]:+.3f} Å")

        plt.figure()
        plt.axhline(0.0, linewidth=0.8)
        plt.plot(delta_rmsf)
        plt.xlabel("Aligned residue index"); plt.ylabel(f"ΔRMSF (Å)  [{args.label1} - {args.label2}]")
        plt.title("Aligned-residue RMSF difference")
        plt.tight_layout()
        plt.savefig(f"{args.out_prefix}_rmsf_delta.pdf", dpi=CONFIG["plot"]["dpi"]); plt.close()

        summary["trajectory_validation"]["rmsf_profile"] = {
            **rmsf_profile_stats,
            "delta_definition": f"{args.label1} - {args.label2}",
            "top_absolute_differences": [
                {
                    "aligned_index": int(i),
                    "residue1": res_label(aligned_res1[i]),
                    "residue2": res_label(aligned_res2[i]),
                    "rmsf1": float(rmsf1_aln[i]),
                    "rmsf2": float(rmsf2_aln[i]),
                    "delta": float(delta_rmsf[i])
                }
                for i in top_rmsf_idx
            ]
        }
        
        summary["rmsf"] = {
            "full": {
                "mean1": float(rmsf1_full.mean()),
                "std1": float(rmsf1_full.std()),
                "mean2": float(rmsf2_full.mean()),
                "std2": float(rmsf2_full.std()),
                "ks_p": float(ks_full.pvalue)
            },
            "aligned": {
                "mean1": float(rmsf1_aln.mean()),
                "std1": float(rmsf1_aln.std()),
                "mean2": float(rmsf2_aln.mean()),
                "std2": float(rmsf2_aln.std()),
                "ks_p": float(ks_aln.pvalue)
            }
        }
            
        summary["units"]["rmsf"] = "Å"

        # ───────── Rg ─────────
        rg1=md.compute_rg(prot1)*10.0; rg2=md.compute_rg(prot2)*10.0
        ks_rg=ks_2samp(rg1,rg2)
        logging.info("Rg | %s %.2f±%.2f | %s %.2f±%.2f | KS‑p %.3e",
                     args.label1,rg1.mean(),rg1.std(),
                     args.label2,rg2.mean(),rg2.std(),ks_rg.pvalue)
        plt.figure()
        plt.plot(rg1,label=args.label1); plt.plot(rg2,label=args.label2)
        plt.xlabel("Frame"); plt.ylabel("Rg (Å)")
        plt.title(f"Radius of gyration  KS p={ks_rg.pvalue:.3e}")
        plt.legend(); plt.tight_layout()
        plt.savefig(f"{args.out_prefix}_rg.pdf",dpi=CONFIG["plot"]["dpi"]); plt.close()

        summary["rg"] = {
            "mean1": float(rg1.mean()), "std1": float(rg1.std()),
            "mean2": float(rg2.mean()), "std2": float(rg2.std()),
            "ks_p": float(ks_rg.pvalue)
        }

        summary["units"]["rg"] = "Å"

        # ───────── TRAJECTORY VALIDATION: GLOBAL-SIZE STABILITY ─────────
        stability_fraction = CONFIG["trajectory_validation"]["stability_fraction"]
        rg_stability1 = trajectory_stability_stats(rg1, stability_fraction)
        rg_stability2 = trajectory_stability_stats(rg2, stability_fraction)
        rg_mean_delta = float(rg1.mean() - rg2.mean())
        rg_relative_delta_pct = (
            100.0 * rg_mean_delta / rg2.mean() if abs(rg2.mean()) > 1.0e-12 else None
        )

        print("\nTrajectory Validation — Radius of Gyration Stability")
        print("=" * 52)
        print(f"  Mean ΔRg ({args.label1} - {args.label2}) : {rg_mean_delta:+.3f} Å")
        if rg_relative_delta_pct is not None:
            print(f"  Relative mean difference       : {rg_relative_delta_pct:+.3f}%")
        for label, stats in ((args.label1, rg_stability1), (args.label2, rg_stability2)):
            print(f"  {label:<20} first→last window Δ : {stats['last_minus_first']:+.3f} Å "
                  f"(linear projected Δ {stats['projected_change_over_trajectory']:+.3f} Å)")

        summary["trajectory_validation"]["rg_stability"] = {
            "mean_difference": rg_mean_delta,
            "relative_mean_difference_percent": rg_relative_delta_pct,
            args.label1: rg_stability1,
            args.label2: rg_stability2
        }

        # ───────── BACKBONE RMSD vs TIME ─────────
        bb_idx1 = prot1.topology.select("backbone")
        bb_idx2 = prot2.topology.select("backbone")
        rmsd1 = md.rmsd(prot1, prot1, 0, bb_idx1)*10.0
        rmsd2 = md.rmsd(prot2, prot2, 0, bb_idx2)*10.0
        plt.figure()
        plt.plot(rmsd1,label=args.label1); plt.plot(rmsd2,label=args.label2)
        plt.xlabel("Frame"); plt.ylabel("Backbone RMSD to start (Å)")
        plt.title("Backbone RMSD vs time"); plt.legend(); plt.tight_layout()
        plt.savefig(f"{args.out_prefix}_rmsd_time.pdf",dpi=CONFIG["plot"]["dpi"]); plt.close()
        
        # Record backbone RMSD statistics in summary
        summary["rmsd_backbone"] = {
            args.label1: {
                "mean": float(rmsd1.mean()),
                "std": float(rmsd1.std()),
                "max": float(rmsd1.max())
            },
            args.label2: {
                "mean": float(rmsd2.mean()),
                "std": float(rmsd2.std()),
                "max": float(rmsd2.max())
            }
        }

        summary["units"]["rmsd_backbone"] = "Å"

        # ───────── TRAJECTORY VALIDATION: RMSD STABILITY + COMMON REFERENCE ─────────
        rmsd_stability1 = trajectory_stability_stats(rmsd1, stability_fraction)
        rmsd_stability2 = trajectory_stability_stats(rmsd2, stability_fraction)
        summary["trajectory_validation"]["backbone_rmsd_stability"] = {
            "reference": "each trajectory's own frame 0 (existing RMSD definition)",
            args.label1: rmsd_stability1,
            args.label2: rmsd_stability2
        }

        # The existing RMSD intentionally remains unchanged.  This additional RMSD
        # uses the common aligned-Cα set and frame 0 of trajectory 1 as the same
        # structural reference for BOTH trajectories.
        rmsd_common1 = md.rmsd(
            prot1, prot1, 0, atom_indices=m1_arr, ref_atom_indices=m1_arr
        ) * 10.0
        rmsd_common2 = md.rmsd(
            prot2, prot1, 0, atom_indices=m2_arr, ref_atom_indices=m1_arr
        ) * 10.0

        plt.figure()
        plt.plot(rmsd_common1, label=args.label1)
        plt.plot(rmsd_common2, label=args.label2)
        plt.xlabel("Frame"); plt.ylabel("Aligned Cα RMSD to common reference (Å)")
        plt.title(f"Common-reference RMSD ({args.label1} frame 0)")
        plt.legend(); plt.tight_layout()
        plt.savefig(f"{args.out_prefix}_rmsd_common_reference.pdf", dpi=CONFIG["plot"]["dpi"]); plt.close()

        common_rmsd_stability1 = trajectory_stability_stats(rmsd_common1, stability_fraction)
        common_rmsd_stability2 = trajectory_stability_stats(rmsd_common2, stability_fraction)
        print("\nTrajectory Validation — Common-reference Cα RMSD")
        print("=" * 52)
        print(f"  Common reference              : {args.label1}, frame 0")
        print(f"  {args.label1:<20} mean±std     : {rmsd_common1.mean():.3f} ± {rmsd_common1.std():.3f} Å")
        print(f"  {args.label2:<20} mean±std     : {rmsd_common2.mean():.3f} ± {rmsd_common2.std():.3f} Å")
        print(f"  {args.label1:<20} first→last Δ : {common_rmsd_stability1['last_minus_first']:+.3f} Å")
        print(f"  {args.label2:<20} first→last Δ : {common_rmsd_stability2['last_minus_first']:+.3f} Å")

        summary["trajectory_validation"]["common_reference_rmsd"] = {
            "reference": {"trajectory": args.label1, "frame": 0, "selection": "aligned Cα"},
            args.label1: {
                "mean": float(rmsd_common1.mean()),
                "std": float(rmsd_common1.std()),
                "max": float(rmsd_common1.max()),
                "stability": common_rmsd_stability1
            },
            args.label2: {
                "mean": float(rmsd_common2.mean()),
                "std": float(rmsd_common2.std()),
                "max": float(rmsd_common2.max()),
                "stability": common_rmsd_stability2
            }
        }
        summary["units"]["common_reference_rmsd"] = "Å"

        # ───────── TRAJECTORY VALIDATION: INTERNAL Cα DISTANCE MAP ─────────
        # Internal distances are invariant to overall translation/rotation and are
        # therefore a useful direct comparison of the sampled protein architecture.
        if n_aligned >= 2:
            local_pairs = np.array(
                [(i, j) for i in range(n_aligned) for j in range(i + 1, n_aligned)],
                dtype=int
            )
            atom_pairs1 = np.column_stack((m1_arr[local_pairs[:, 0]], m1_arr[local_pairs[:, 1]]))
            atom_pairs2 = np.column_stack((m2_arr[local_pairs[:, 0]], m2_arr[local_pairs[:, 1]]))

            # periodic=False deliberately compares intramolecular geometry without
            # applying a minimum-image shortcut to long intraprotein distances.
            # Compute in pair chunks so this addition does not allocate a potentially
            # enormous (n_frames × n_pairs) matrix for larger proteins.
            pair_chunk = CONFIG["trajectory_validation"]["distance_pair_chunk_size"]
            mean_dist1 = mean_pair_distances(prot1, atom_pairs1, pair_chunk)
            mean_dist2 = mean_pair_distances(prot2, atom_pairs2, pair_chunk)
            delta_dist = mean_dist1 - mean_dist2
            distance_stats = paired_profile_stats(mean_dist1, mean_dist2)

            delta_matrix = np.zeros((n_aligned, n_aligned), dtype=float)
            delta_matrix[local_pairs[:, 0], local_pairs[:, 1]] = delta_dist
            delta_matrix[local_pairs[:, 1], local_pairs[:, 0]] = delta_dist
            vmax = max(float(np.max(np.abs(delta_matrix))), 1.0e-6)

            plt.figure(figsize=CONFIG["plot"]["figsize"]["dccm"])
            sns.heatmap(delta_matrix, cmap="bwr", center=0.0, vmin=-vmax, vmax=vmax, square=True,
                        cbar_kws={"label": f"Δ mean Cα distance (Å): {args.label1} - {args.label2}"})
            plt.xlabel("Aligned residue index"); plt.ylabel("Aligned residue index")
            plt.title("Difference in mean aligned-Cα internal distances")
            plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_ca_distance_delta.pdf", dpi=CONFIG["plot"]["dpi"]); plt.close()

            top_n_dist = min(
                CONFIG["trajectory_validation"]["top_n_distance_differences"], len(delta_dist)
            )
            top_dist_idx = np.argsort(np.abs(delta_dist))[::-1][:top_n_dist]

            print("\nTrajectory Validation — Mean Internal Cα Distances")
            print("=" * 52)
            if distance_stats["pearson_r"] is not None:
                print(f"  Pearson r                    : {distance_stats['pearson_r']:.5f}")
                print(f"  Spearman rho                 : {distance_stats['spearman_rho']:.5f}")
            print(f"  Distance RMSE                : {distance_stats['rmse']:.3f} Å")
            print(f"  Mean |Δdistance|             : {distance_stats['mae']:.3f} Å")
            print(f"  Maximum |Δdistance|          : {distance_stats['max_abs_difference']:.3f} Å")
            print(f"  Top {top_n_dist} pair differences:")
            for k in top_dist_idx:
                i, j = local_pairs[k]
                print(f"    {res_label(aligned_res1[i])}–{res_label(aligned_res1[j])} / "
                      f"{res_label(aligned_res2[i])}–{res_label(aligned_res2[j])} "
                      f"Δ={delta_dist[k]:+.3f} Å")

            summary["trajectory_validation"]["ca_internal_distances"] = {
                **distance_stats,
                "n_aligned_ca": int(n_aligned),
                "n_pairs": int(len(local_pairs)),
                "delta_definition": f"{args.label1} - {args.label2}",
                "top_absolute_differences": [
                    {
                        "residue1_system1": res_label(aligned_res1[i]),
                        "residue2_system1": res_label(aligned_res1[j]),
                        "residue1_system2": res_label(aligned_res2[i]),
                        "residue2_system2": res_label(aligned_res2[j]),
                        "mean_distance1": float(mean_dist1[k]),
                        "mean_distance2": float(mean_dist2[k]),
                        "delta": float(delta_dist[k])
                    }
                    for k in top_dist_idx
                    for i, j in [local_pairs[k]]
                ]
            }
            summary["units"]["ca_internal_distances"] = "Å"
        else:
            logging.warning("Internal Cα distance comparison skipped: fewer than two aligned Cα atoms.")
            summary["trajectory_validation"]["ca_internal_distances"] = None

        # ───────── CLUSTERING (k‑means, CA‑RMSD feature matrix) ─────────
        min_frames = CONFIG["min_frames_for_clustering"]
        k_opt, pop1, pop2 = None, [], []
        if prot1.n_frames >= min_frames and prot2.n_frames >= min_frames \
           and len(ca1_all) and len(ca2_all):
        
            ref_frames = min(CONFIG["max_reference_frames"], prot1.n_frames, prot2.n_frames) 
            ref_idx    = np.arange(ref_frames)
        
            # Feature matrix: for each frame, RMSD to each reference frame
            feat1 = np.column_stack([
                md.rmsd(prot1, prot1, i, atom_indices=ca1_all) * 10.0
                for i in ref_idx
            ])                                # shape: (n_frames1, ref_frames)
            feat2 = np.column_stack([
                md.rmsd(prot2, prot2, i, atom_indices=ca2_all) * 10.0
                for i in ref_idx
            ])                                # shape: (n_frames2, ref_frames)
        
            X = np.vstack([feat1, feat2])     # samples = frames (rows)
        
            n_samples = X.shape[0]
            max_k     = min(CONFIG["max_clusters"], n_samples)    # never ask for more clusters than samples
        
            # --- choose k by a simple elbow test ---
            inertias = []
            for k in range(2, max_k + 1):
                km = KMeans(n_clusters=k, n_init="auto", random_state=0).fit(X)
                inertias.append(km.inertia_)
                if k > 3 and inertias[-2] - inertias[-1] < 0.05 * inertias[1]:
                    k_opt = k
                    break
            else:
                k_opt = min(5, max_k)
        
            km = KMeans(n_clusters=k_opt, n_init="auto", random_state=1).fit(X)
            labels = km.labels_
        
            pop1 = np.bincount(labels[:prot1.n_frames], minlength=k_opt)
            pop2 = np.bincount(labels[prot1.n_frames:], minlength=k_opt)
        
            logging.info("Cluster populations (k=%d):", k_opt)
            for i, (p1, p2) in enumerate(zip(pop1, pop2)):
                logging.info("  Cluster %-2d | %-18s %4d (%.1f%%) "
                             "| %-18s %4d (%.1f%%)",
                             i,
                             args.label1, p1, 100 * p1 / prot1.n_frames,
                             args.label2, p2, 100 * p2 / prot2.n_frames)
        
            # Bar‑chart
            idx   = np.arange(k_opt)
            width = 0.4
            plt.figure()
            plt.bar(idx - width/2, pop1 / prot1.n_frames * 100, width, label=args.label1)
            plt.bar(idx + width/2, pop2 / prot2.n_frames * 100, width, label=args.label2)
            plt.xlabel("Cluster"); plt.ylabel("Population (%)")
            plt.title("k‑means cluster populations (CA‑RMSD features)")
            plt.legend(); plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_cluster_populations.pdf", dpi=CONFIG["plot"]["dpi"])
            plt.close()
        else:
            k_opt = None
            pop1 = pop2 = []
            logging.warning("Clustering skipped: insufficient frames or atoms.")

        summary["clustering"] = {
            "k": k_opt,
            "populations": {
                args.label1: list(map(int, pop1)),
                args.label2: list(map(int, pop2))
            }
        }

        # ───────── PCA (aligned Cα only) ─────────
        # Use the same symmetrically common-reference-aligned copies as the
        # aligned RMSF analysis.  This removes rigid-body translation/rotation
        # from BOTH trajectories without modifying the original trajectories.
        n_ca_common = len(m1)
        coords1 = aligned_prot1.xyz[:, m1_arr, :].reshape(aligned_prot1.n_frames, n_ca_common * 3)
        coords2 = aligned_prot2.xyz[:, m2_arr, :].reshape(aligned_prot2.n_frames, n_ca_common * 3)
        
        X_pca = np.vstack([coords1, coords2])  # (n_frames1 + n_frames2) × (n_ca_common*3)
        pca = PCA(n_components=2, random_state=0).fit(X_pca)
        
        pc1 = pca.transform(coords1)
        pc2 = pca.transform(coords2)
        
        plt.figure()
        plt.scatter(pc1[:, 0], pc1[:, 1], s=5, alpha=0.5, label=args.label1)
        plt.scatter(pc2[:, 0], pc2[:, 1], s=5, alpha=0.5, label=args.label2)
        plt.xlabel("PC1"); plt.ylabel("PC2")
        plt.title("PCA of aligned common Cα coordinates")
        plt.legend(); plt.tight_layout()
        plt.savefig(f"{args.out_prefix}_pca_scatter.pdf", dpi=CONFIG["plot"]["dpi"])
        plt.close()

        # Variance explained
        var_exp = pca.explained_variance_ratio_
        print("\nPrincipal Component Variance Explained:")
        print(f"  PC1: {var_exp[0]*100:.1f}%")
        print(f"  PC2: {var_exp[1]*100:.1f}%")
        
        # Spread in PCA space
        def pca_stats(pc):
            return np.std(pc[:, 0]), np.std(pc[:, 1]), np.mean(pc[:, 0]), np.mean(pc[:, 1])
        
        sd1_pc1, sd1_pc2, mu1_pc1, mu1_pc2 = pca_stats(pc1)
        sd2_pc1, sd2_pc2, mu2_pc1, mu2_pc2 = pca_stats(pc2)
        
        print("\nConformational spread (std. dev. along PC1/PC2):")
        print(f"  {args.label1:<18} PC1={sd1_pc1:.2f}  PC2={sd1_pc2:.2f}")
        print(f"  {args.label2:<18} PC1={sd2_pc1:.2f}  PC2={sd2_pc2:.2f}")
        
        # Distance between PCA centroids
        centroid_dist = np.linalg.norm([mu1_pc1 - mu2_pc1, mu1_pc2 - mu2_pc2])
        print(f"\nCentroid distance ({args.label1} vs {args.label2}): {centroid_dist:.2f}\n")

        summary["pca"] = {
            "explained_variance": [float(v) for v in var_exp],
            "std_pc1": [float(sd1_pc1), float(sd2_pc1)],
            "std_pc2": [float(sd1_pc2), float(sd2_pc2)],
            "centroid_distance": float(centroid_dist)
        }

        summary["units"]["pca"] = "unitless"

        # ───────── Hydrogen‑bond persistence ─────────
        hb1 = md.wernet_nilsson(prot1, periodic=False, sidechain_only=False)
        hb2 = md.wernet_nilsson(prot2, periodic=False, sidechain_only=False)
        
        def occupancy_by_label(hb_frame_list, topology):
            """
            Compute residue-pair H-bond occupancy as the fraction of frames in
            which at least one atom-level H-bond connects the two residues.

            Multiple atom-level H-bonds between the same residue pair in a single
            frame count only once, so every returned occupancy is in [0, 1].
            Residue labels are sorted deliberately: this is an undirected residue
            interaction rather than an atom-level donor→acceptor identity.
            """
            if not hb_frame_list:
                return {}
        
            counts = {}
            n_frames = len(hb_frame_list)
        
            for frame_array in hb_frame_list:
                frame_pairs = set()
                for bond in frame_array:
                    donor = int(bond[0])
                    acceptor = int(bond[-1])
                    d_lbl = res_label(topology.atom(donor))
                    a_lbl = res_label(topology.atom(acceptor))
                    frame_pairs.add(tuple(sorted((d_lbl, a_lbl))))

                for key in frame_pairs:
                    counts[key] = counts.get(key, 0) + 1
        
            return {k: v / n_frames for k, v in counts.items()}
        
        occ1 = occupancy_by_label(hb1, prot1.topology)
        occ2 = occupancy_by_label(hb2, prot2.topology)
        
        persistent = {
            k: (occ1.get(k, 0.0), occ2.get(k, 0.0))
            for k in set(occ1) | set(occ2)
            if max(occ1.get(k, 0.0), occ2.get(k, 0.0)) >= CONFIG["hbond_threshold"]
        }
        
        # Write CSV output
        if persistent:
            csvfile = f"{args.out_prefix}_hbonds_persistent.csv"
            with open(csvfile, "w", newline="") as fh:
                writer = csv.writer(fh)
                writer.writerow(["Donor", "Acceptor", args.label1, args.label2])
                for (don, acc), (f1, f2) in sorted(persistent.items(), key=lambda x: -max(x[1])):
                    writer.writerow([don, acc, f"{f1:.2f}", f"{f2:.2f}"])
            logging.info("Wrote persistent H‑bonds to %s (%d entries)", csvfile, len(persistent))
        else:
            logging.info("No hydrogen bonds ≥30%% occupancy in either trajectory")
        
        # Print summary to screen
        top_n = CONFIG["top_n_persistent_hbonds"]
        hb_threshold_pct = int(CONFIG["hbond_threshold"] * 100)
        sorted_persistent = sorted(persistent.items(), key=lambda x: -max(x[1]))
        
        print(f"\nTop persistent hydrogen bonds (≥{hb_threshold_pct}% occupancy in either trajectory):")
        print(f"{'Residue 1':<18} ↔ {'Residue 2':<18} | {args.label1:^18} | {args.label2:^18}")
        print("-" * 70)
        
        for (don, acc), (f1, f2) in sorted_persistent[:top_n]:
            print(f"{don:<18} ↔ {acc:<18} | {f1:>8.2f}         | {f2:>8.2f}")
        
        # Summary stats
        n1 = sum(1 for v in persistent.values() if v[0] >= CONFIG["hbond_threshold"])
        n2 = sum(1 for v in persistent.values() if v[1] >= CONFIG["hbond_threshold"])
        shared = sum(1 for v in persistent.values() if v[0] >= CONFIG["hbond_threshold"] and v[1] >= CONFIG["hbond_threshold"])
        only1 = sum(1 for v in persistent.values() if v[0] >= CONFIG["hbond_threshold"] and v[1] < CONFIG["hbond_threshold"])
        only2 = sum(1 for v in persistent.values() if v[1] >= CONFIG["hbond_threshold"] and v[0] < CONFIG["hbond_threshold"])
        
        print(f"\n{args.label1}: {n1} persistent H-bonds (≥{hb_threshold_pct}%)")
        print(f"{args.label2}: {n2} persistent H-bonds (≥{hb_threshold_pct}%)")
        print(f"Shared: {shared}")
        print(f"Exclusive to {args.label1}: {only1}")
        print(f"Exclusive to {args.label2}: {only2}\n")
        
        summary["hbonds"] = {
            "n1": n1,
            "n2": n2,
            "shared": shared,
            "exclusive1": only1,
            "exclusive2": only2,
            "occupancy_definition": "fraction of frames with >=1 atom-level H-bond for the residue pair"
        }
        
        summary["top_hbonds"] = [
            {
                # residue1/residue2 are the scientifically correct labels for this
                # order-independent residue-pair occupancy.  donor/acceptor are
                # retained for backwards compatibility with existing JSON consumers.
                "residue1": don,
                "residue2": acc,
                "donor": don,
                "acceptor": acc,
                args.label1: float(f1),
                args.label2: float(f2)
            }
            for (don, acc), (f1, f2) in sorted_persistent
        ]
        # ───────── DSSP heat‑map ─────────
        n_disagree = 0
        try:
            ss1 = md.compute_dssp(prot1)
            ss2 = md.compute_dssp(prot2)
            struct_map={'H':0,'E':1,'C':2,'G':3,'I':4,'B':5,'T':6,'S':7}
            ss_num1=np.vectorize(struct_map.get)(ss1)
            ss_num2=np.vectorize(struct_map.get)(ss2)
        
            plt.figure(figsize=CONFIG["plot"]["figsize"]["dssp"])
            for i,(ss,label) in enumerate([(ss_num1,args.label1),(ss_num2,args.label2)],1):
                plt.subplot(2,1,i)
                sns.heatmap(ss.T,cbar=i==1,cmap="tab20",
                            cbar_kws={"label":"DSSP"})
                plt.ylabel("Residue"); plt.xlabel("Frame")
                plt.title(f"{label} – DSSP")
            plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_dssp.pdf", dpi=CONFIG["plot"]["dpi"])
            plt.close()
        
            # Difference map (on aligned residues only)
            r1_idx = [prot1.topology.atom(i).residue.index for i in m1]
            r2_idx = [prot2.topology.atom(i).residue.index for i in m2]
            aligned_ss1 = ss1[:, r1_idx]
            aligned_ss2 = ss2[:, r2_idx]
            dssp_diff = (aligned_ss1 != aligned_ss2).astype(np.int_)
        
            plt.figure(figsize=CONFIG["plot"]["figsize"]["dssp_diff"])
            sns.heatmap(dssp_diff.T, cmap="Reds", cbar_kws={"label": "Mismatch (1 = different)"})
            plt.xlabel("Frame")
            plt.ylabel("Aligned residue index")
            plt.title(f"DSSP secondary structure difference ({args.label1} vs {args.label2})")
            plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_dssp_difference.pdf", dpi=CONFIG["plot"]["dpi"])
            plt.close()
        
            diff_summary = dssp_diff.mean(axis=0)
            n_disagree = np.sum(diff_summary > 0.5)
            print(f"\n{n_disagree} aligned residues differed in secondary structure "
                  f"in >50% of frames.\n")

            summary["dssp_diff"] = {
                "residues_differing_gt_50pct": int(n_disagree)
            }
            summary["dssp_diff_pct"] = float(n_disagree) / float(len(aligned_res1))
        
        except Exception as e:
            logging.warning("DSSP calculation failed (%s)", e)
            summary["dssp_diff"] = {"residues_differing_gt_50pct": None}

        # ───────── DCCM (Dynamic Cross-Correlation Matrix) for Aligned Cα Atoms ─────────
        TOP_N_DCCM = CONFIG["dccm"]["top_n"]
        
        def compute_dccm(traj: md.Trajectory, atom_indices: list[int]) -> np.ndarray:
            """Compute normalized DCCM for a given set of atoms (typically Cα)."""
            xyz = traj.xyz[:, atom_indices, :]  # (n_frames, n_atoms, 3)
            disp = xyz - xyz.mean(axis=0, keepdims=True)  # displacement vectors
            n_atoms = xyz.shape[1]
            corr = np.zeros((n_atoms, n_atoms))
            for i in range(n_atoms):
                for j in range(n_atoms):
                    vi = disp[:, i, :].reshape(len(traj), 3)
                    vj = disp[:, j, :].reshape(len(traj), 3)
                    numerator = np.sum(np.sum(vi * vj, axis=1))
                    denom = np.sqrt(np.sum(vi ** 2) * np.sum(vj ** 2))
                    corr[i, j] = numerator / denom if denom > 0 else 0.0
            return np.clip(corr, -1.0, 1.0)
        
        def top_correlated_residues(dccm_matrix, residues, top_n=5, sign="positive"):
            """Return top correlated or anti-correlated residue pairs."""
            N = dccm_matrix.shape[0]
            pairs = []
            for i in range(N):
                for j in range(i + 1, N):
                    corr = dccm_matrix[i, j]
                    if sign == "positive":
                        score = corr
                    elif sign == "negative":
                        score = -corr
                    else:
                        score = abs(corr)
                    pairs.append(((i, j), corr, score))
            sorted_pairs = sorted(pairs, key=lambda x: x[2], reverse=True)[:top_n]
            results = []
            for (i, j), corr_val, _ in sorted_pairs:
                res1 = residues[i]
                res2 = residues[j]
                label1 = f"{chr(65 + res1.chain.index)}:{res1.name}{res1.resSeq}"
                label2 = f"{chr(65 + res2.chain.index)}:{res2.name}{res2.resSeq}"
                results.append((label1, label2, corr_val))
            return results
        
        # Compute DCCMs using the same symmetrically common-reference-aligned
        # trajectory copies used for aligned RMSF and PCA.
        dccm1 = compute_dccm(aligned_prot1, m1_arr)
        dccm2 = compute_dccm(aligned_prot2, m2_arr)
        
        # Plot DCCM matrices
        for mat, label in zip([dccm1, dccm2], [args.label1, args.label2]):
            plt.figure(figsize=CONFIG["plot"]["figsize"]["dccm"])
            sns.heatmap(mat, vmin=-1, vmax=1, cmap="coolwarm", square=True, cbar_kws={"label": "DCCM correlation"})
            plt.title(f"DCCM (aligned Cα) — {label}")
            plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_dccm_ca_{safe(label)}.pdf", dpi=CONFIG["plot"]["dpi"])
            plt.close()

        # ───────── Differential DCCM (ΔDCCM) Heatmap ─────────
        delta_dccm = dccm1 - dccm2  # ΔDCCM = label1 - label2

        plt.figure(figsize=CONFIG["plot"]["figsize"]["dccm"])
        sns.heatmap(delta_dccm, cmap="bwr", center=0, vmin=-1, vmax=1,
                    square=True, cbar_kws={"label": f"ΔDCCM ({args.label1} - {args.label2})"})
        plt.title(f"Differential DCCM ({args.label1} - {args.label2}; Aligned Cα Atoms)")
        plt.tight_layout()
        plt.savefig(f"{args.out_prefix}_dccm_delta.pdf", dpi=CONFIG["plot"]["dpi"])
        plt.close()
        print("  Differential DCCM saved to PDF.")
        
        # DCCM summary stats
        avg_corr_overall = np.mean(dccm1 - dccm2)
        max_diff = np.max(np.abs(dccm1 - dccm2))
        print(f"\nDCCM Summary (Aligned Cα Atoms):")
        print(f"  Mean difference: {avg_corr_overall:.3f}")
        print(f"  Max abs. difference: {max_diff:.3f}")
        
        print(f"\nTop {TOP_N_DCCM} most positively correlated residue pairs ({args.label1}):")
        for a, b, c in top_correlated_residues(dccm1, aligned_res1, top_n=TOP_N_DCCM, sign="positive"):
            print(f"  {a} ↔ {b}  |  corr = {c:.3f}")
        
        print(f"\nTop {TOP_N_DCCM} most negatively correlated residue pairs ({args.label1}):")
        for a, b, c in top_correlated_residues(dccm1, aligned_res1, top_n=TOP_N_DCCM, sign="negative"):
            print(f"  {a} ↔ {b}  |  corr = {c:.3f}")
        
        print(f"\nTop {TOP_N_DCCM} most positively correlated residue pairs ({args.label2}):")
        for a, b, c in top_correlated_residues(dccm2, aligned_res2, top_n=TOP_N_DCCM, sign="positive"):
            print(f"  {a} ↔ {b}  |  corr = {c:.3f}")
        
        print(f"\nTop {TOP_N_DCCM} most negatively correlated residue pairs ({args.label2}):")
        for a, b, c in top_correlated_residues(dccm2, aligned_res2, top_n=TOP_N_DCCM, sign="negative"):
            print(f"  {a} ↔ {b}  |  corr = {c:.3f}")

        # Compute ΔDCCM top changes
        delta_pairs = []
        for i in range(delta_dccm.shape[0]):
            for j in range(i + 1, delta_dccm.shape[1]):
                delta_pairs.append(((i, j), delta_dccm[i, j]))
        
        delta_top_positive = sorted(delta_pairs, key=lambda x: x[1], reverse=True)[:TOP_N_DCCM]
        delta_top_negative = sorted(delta_pairs, key=lambda x: x[1])[:TOP_N_DCCM]
        
        def format_delta_pair(i, j, val):
            res1 = aligned_res1[i]
            res2 = aligned_res1[j]
            label1 = f"{chr(65 + res1.chain.index)}:{res1.name}{res1.resSeq}"
            label2 = f"{chr(65 + res2.chain.index)}:{res2.name}{res2.resSeq}"
            return {"res1": label1, "res2": label2, "delta_corr": float(val)}
        
        print(f"\nTop {TOP_N_DCCM} ΔDCCM (positive):")
        for (i, j), val in delta_top_positive:
            p = format_delta_pair(i, j, val)
            print(f"  {p['res1']} ↔ {p['res2']}  |  Δcorr = {p['delta_corr']:.3f}")
        
        print(f"\nTop {TOP_N_DCCM} ΔDCCM (negative):")
        for (i, j), val in delta_top_negative:
            p = format_delta_pair(i, j, val)
            print(f"  {p['res1']} ↔ {p['res2']}  |  Δcorr = {p['delta_corr']:.3f}")

        dccm_stats = {
            "alignment_length": len(aligned_res1),
            "mean_delta": float(avg_corr_overall),
            "max_abs_delta": float(max_diff),
            "top_positive": {
                args.label1: [
                    {"res1": a, "res2": b, "corr": float(c)}
                    for a, b, c in top_correlated_residues(
                        dccm1, aligned_res1, top_n=TOP_N_DCCM, sign="positive")
                ],
                args.label2: [
                    {"res1": a, "res2": b, "corr": float(c)}
                    for a, b, c in top_correlated_residues(
                        dccm2, aligned_res2, top_n=TOP_N_DCCM, sign="positive")
                ]
            },
            "top_negative": {
                args.label1: [
                    {"res1": a, "res2": b, "corr": float(c)}
                    for a, b, c in top_correlated_residues(
                        dccm1, aligned_res1, top_n=TOP_N_DCCM, sign="negative")
                ],
                args.label2: [
                    {"res1": a, "res2": b, "corr": float(c)}
                    for a, b, c in top_correlated_residues(
                        dccm2, aligned_res2, top_n=TOP_N_DCCM, sign="negative")
                ]
            },
            "delta_top_positive": [
                format_delta_pair(i, j, val) for (i, j), val in delta_top_positive
            ],
            "delta_top_negative": [
                format_delta_pair(i, j, val) for (i, j), val in delta_top_negative
            ]
        }
        summary["dccm_overall"] = dccm_stats
        summary["units"]["correlation"] = "unitless"

        # ───────── Ligand detection or reporting ─────────
        ligand_resname = args.ligand.upper() if args.ligand else None
        loop_residues = None 
        ligand_res = None
        
        if ligand_resname:
            def find_ligand_atoms(top, ligand_name):
                atoms = []
                for res in top.residues:
                    resname_norm = res.name[:3].upper()
                    if resname_norm == ligand_name.upper():
                        atoms.extend([atom.index for atom in res.atoms])
                return atoms
            
            ligand1 = find_ligand_atoms(t1.topology, ligand_resname)
            ligand2 = find_ligand_atoms(t2.topology, ligand_resname)
            
            if len(ligand1) == 0 or len(ligand2) == 0:
                print(f"\nWARNING: Ligand '{ligand_resname}' not found in one or both trajectories.\n")
                print(f"  Ligand candidates in t1: {[res.name for res in t1.topology.residues if res.name[:3].upper() == ligand_resname.upper()]}")
                print(f"  Ligand candidates in t2: {[res.name for res in t2.topology.residues if res.name[:3].upper() == ligand_resname.upper()]}")
            else:
                print(f"\nLigand '{ligand_resname}' found in both structures.")
                print(f"  {args.label1}: {len(ligand1)} atoms")
                print(f"  {args.label2}: {len(ligand2)} atoms")
                # Optional: Add analysis code here (RMSD, SASA, etc.)
        else:
            # No ligand specified — report all non-standard residues
            def get_nonprotein_resnames(top):
                return sorted(set(
                    res.name for res in top.residues
                    if not res.is_protein and res.name not in {"HOH", "NA", "CL"}
                ))
        
            other1 = get_nonprotein_resnames(t1.topology)
            other2 = get_nonprotein_resnames(t2.topology)
            unknowns = sorted(set(other1) | set(other2))
        
            if unknowns:
                print("\nLigand not specified.")
                print("These are the unknown / non-protein residues found:")
                print("  " + " ".join(unknowns) + "\n")
            else:
                print("\nLigand not specified, and no non-protein residues found.\n")

        if ligand_resname and len(ligand1) > 0 and len(ligand2) > 0:
            print("\nLigand Stability Analysis")
            print("=" * 30)

            summary["ligand"] = {
                "resname": ligand_resname
            }
        
            # Slice trajectories
            lig_traj1 = t1.atom_slice(ligand1)
            lig_traj2 = t2.atom_slice(ligand2)
        
            # ───── Ligand RMSD (vs. frame 0) ─────
            lig_rmsd1 = md.rmsd(lig_traj1, lig_traj1, 0) * CONFIG["rmsf_scale"]  # Å
            lig_rmsd2 = md.rmsd(lig_traj2, lig_traj2, 0) * CONFIG["rmsf_scale"]  # Å
        
            plt.figure()
            plt.plot(lig_rmsd1, label=f"{args.label1}")
            plt.plot(lig_rmsd2, label=f"{args.label2}")
            plt.xlabel("Frame"); plt.ylabel("Ligand RMSD (Å)")
            plt.title(f"Ligand RMSD (resname {ligand_resname})")
            plt.legend(); plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_ligand_rmsd.pdf", dpi=CONFIG["plot"]["dpi"])
            plt.close()
        
            print(f"  Ligand RMSD (Å):")
            print(f"    {args.label1:<20} mean={lig_rmsd1.mean():.2f}  std={lig_rmsd1.std():.2f}  max={lig_rmsd1.max():.2f}")
            print(f"    {args.label2:<20} mean={lig_rmsd2.mean():.2f}  std={lig_rmsd2.std():.2f}  max={lig_rmsd2.max():.2f}")

            summary["ligand"]["rmsd"] = {
                args.label1: {
                    "mean": float(lig_rmsd1.mean()),
                    "std": float(lig_rmsd1.std()),
                    "max": float(lig_rmsd1.max())
                },
                args.label2: {
                    "mean": float(lig_rmsd2.mean()),
                    "std": float(lig_rmsd2.std()),
                    "max": float(lig_rmsd2.max())
                }
            }

            # ───── Ligand SASA (Shrake-Rupley) ─────
            sasa1 = md.shrake_rupley(lig_traj1).sum(axis=1)
            sasa2 = md.shrake_rupley(lig_traj2).sum(axis=1)
        
            plt.figure()
            plt.plot(sasa1, label=args.label1)
            plt.plot(sasa2, label=args.label2)
            plt.xlabel("Frame"); plt.ylabel("Ligand SASA (Å²)")
            plt.title(f"Ligand Solvent Accessibility (resname {ligand_resname})")
            plt.legend(); plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_ligand_sasa.pdf", dpi=CONFIG["plot"]["dpi"])
            plt.close()
        
            print(f"  Ligand SASA (Å²):")
            print(f"    {args.label1:<20} mean={sasa1.mean():.1f}  std={sasa1.std():.1f}")
            print(f"    {args.label2:<20} mean={sasa2.mean():.1f}  std={sasa2.std():.1f}")

            summary["ligand"]["sasa"] = {
                args.label1: {
                    "mean": float(sasa1.mean()),
                    "std": float(sasa1.std())
                },
                args.label2: {
                    "mean": float(sasa2.mean()),
                    "std": float(sasa2.std())
                }
            }
            summary["units"]["sasa"] = "Å²"
            
            # Ligand RMSF (per atom)
            rmsf_lig1 = md.rmsf(lig_traj1, reference=lig_traj1) * 10.0
            rmsf_lig2 = md.rmsf(lig_traj2, reference=lig_traj2) * 10.0
            
            # Perform KS test
            ks_lig = ks_2samp(rmsf_lig1, rmsf_lig2)
            
            # Plot
            plt.figure()
            plt.plot(rmsf_lig1, label=args.label1)
            plt.plot(rmsf_lig2, label=args.label2)
            plt.xlabel("Ligand Atom Index"); plt.ylabel("RMSF (Å)")
            plt.title(f"Ligand Atom Flexibility (RMSF)")
            plt.legend(); plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_ligand_rmsf.pdf", dpi=CONFIG["plot"]["dpi"])
            plt.close()
            
            # Print
            print("  Ligand RMSF:")
            print(f"    {args.label1:<20} mean={rmsf_lig1.mean():.2f}  std={rmsf_lig1.std():.2f}  max={rmsf_lig1.max():.2f}")
            print(f"    {args.label2:<20} mean={rmsf_lig2.mean():.2f}  std={rmsf_lig2.std():.2f}  max={rmsf_lig2.max():.2f}")
            print(f"    KS‑p (ligand RMSF): {ks_lig.pvalue:.2e}")

            summary["ligand"]["rmsf"] = {
                args.label1: {
                    "mean": float(rmsf_lig1.mean()),
                    "std": float(rmsf_lig1.std()),
                    "max": float(rmsf_lig1.max())
                },
                args.label2: {
                    "mean": float(rmsf_lig2.mean()),
                    "std": float(rmsf_lig2.std()),
                    "max": float(rmsf_lig2.max())
                },
                "ks_p": float(ks_lig.pvalue)
            }

            # ───── Ligand–Pocket H‑bonds (≥30 % persistence) ─────
            print("\nLigand–Pocket H‑bond Analysis")
            print("=" * 30)
            
            hb_threshold = 0.30          # 30 %
            hb_threshold_pct = int(hb_threshold * 100)
            
            # --- helper -----------------------------------------------------------
            
            def collect_lig_hbonds(traj, ligand_atom_ids_in_slice):
                """
                Return ligand↔protein residue-pair H-bond occupancy fractions.

                A residue pair is counted at most once per frame even when several
                atom-level H-bonds connect the same ligand/protein residue pair.
                This keeps persistence values in the physically meaningful [0, 1]
                range and matches the global residue-pair H-bond analysis.
                """
                lig_set   = set(ligand_atom_ids_in_slice)
                n_frames  = len(traj)
                counts    = {}
            
                for frame_hb in md.wernet_nilsson(traj):
                    frame_pairs = set()
                    for hbond in frame_hb:
                        don, acc = map(int, hbond[[0, -1]])
                        in_lig = (don in lig_set, acc in lig_set)
            
                        if in_lig.count(True) == 1:          # XOR → one ligand, one protein
                            don_lbl = res_label(traj.topology.atom(don))
                            acc_lbl = res_label(traj.topology.atom(acc))
                            frame_pairs.add(tuple(sorted((don_lbl, acc_lbl))))

                    for key in frame_pairs:
                        counts[key] = counts.get(key, 0) + 1
            
                return {k: v / n_frames for k, v in counts.items()}
            # ----------------------------------------------------------------------
            
            # Build “protein + ligand” sliced trajectories for BOTH systems
            lp1_ids = np.unique(np.concatenate([prot1.topology.select("all"), ligand1]))
            lp2_ids = np.unique(np.concatenate([prot2.topology.select("all"), ligand2]))
            
            traj_lp1 = t1.atom_slice(lp1_ids)
            traj_lp2 = t2.atom_slice(lp2_ids)
            
            # Map ligand atom IDs *within each slice*
            lig1_in_slice = np.where(np.isin(lp1_ids, ligand1))[0]
            lig2_in_slice = np.where(np.isin(lp2_ids, ligand2))[0]
            
            hb1 = collect_lig_hbonds(traj_lp1, lig1_in_slice)
            hb2 = collect_lig_hbonds(traj_lp2, lig2_in_slice)
            
            # Union & threshold filter
            pairs = {
                k: (hb1.get(k, 0.0), hb2.get(k, 0.0))
                for k in set(hb1) | set(hb2)
                if max(hb1.get(k, 0.0), hb2.get(k, 0.0)) >= hb_threshold
            }
            
            if pairs:
                print(f"  {len(pairs)} ligand–protein H‑bonds ≥{hb_threshold_pct}% persistence")
                for (resA, resB), (f1, f2) in sorted(pairs.items(),
                                                     key=lambda x: -max(x[1])):
                    print(f"    {resA} ↔ {resB} | {args.label1}: {f1:.2f}  {args.label2}: {f2:.2f}")
            
                # Store in JSON summary -------------------------------------------
                summary["ligand"]["hbond_persistence"] = [
                    {
                        "residue1": resA,
                        "residue2": resB,
                        # Retain the legacy keys so existing JSON consumers do not
                        # break; these residue-level pairs are intentionally undirected.
                        "donor": resA,
                        "acceptor": resB,
                        args.label1: float(f1),
                        args.label2: float(f2)
                    }
                    for (resA, resB), (f1, f2) in pairs.items()
                ]
            else:
                print("  No persistent ligand–pocket H‑bonds found.")

            # ───── Contact Fingerprint (aligned residues, ligand ↔ protein) ─────
            cutoff_nm = CONFIG["contact_fingerprint"]["distance_cutoff_nm"]
                        
            def build_unique_pairs(residues, ligand_atoms):
                """Return (unique_pairs, residue_index_list)."""
                seen = set()
                pairs = []
                res_cols = []
                for i, res in enumerate(residues):
                    for a in res.atoms:
                        if a.element.symbol == 'H':
                            continue                       # heavy atoms only
                        for b in ligand_atoms:
                            p = (a.index, b)
                            if p in seen:                  # drop duplicate pair
                                continue
                            seen.add(p)
                            pairs.append(list(p))
                            res_cols.append(i)            # which residue owns this column
                return pairs, np.array(res_cols, dtype=int)
            
            # ---- Complete Structure ---------------------------------------------------
            pairs1, rescol1 = build_unique_pairs(aligned_res1, ligand1)
            dist1, _ = md.compute_contacts(t1, pairs1, scheme='closest-heavy')
            
            # ---- Published Structure --------------------------------------------------
            pairs2, rescol2 = build_unique_pairs(aligned_res2, ligand2)
            dist2, _ = md.compute_contacts(t2, pairs2, scheme='closest-heavy')
            
            # ---- Occupancy (fraction of frames with any atom < cutoff) ----------------
            fp1 = np.zeros(n_aligned)
            fp2 = np.zeros(n_aligned)
            
            for i in range(n_aligned):
                cols1 = np.where(rescol1 == i)[0]
                cols2 = np.where(rescol2 == i)[0]
                if cols1.size:
                    fp1[i] = (dist1[:, cols1] < cutoff_nm).any(axis=1).mean()
                if cols2.size:
                    fp2[i] = (dist2[:, cols2] < cutoff_nm).any(axis=1).mean()
            
            # ---- Pretty print ---------------------------------------------------------
            print(f"\nLigand–Residue Contact Fingerprint (aligned residues; contact ≥{int(CONFIG['contact_fingerprint']['distance_cutoff_nm']*10):.1f} Å)")
            print(f"{f'Residue ({args.label1})':<22} | {args.label1:^6} | {args.label2:^6}")
            print("-"*44)
            for i, (r1, r2) in enumerate(zip(aligned_res1, aligned_res2)):
                if max(fp1[i], fp2[i]) >= CONFIG["contact_fingerprint"]["min_occupancy"]:           # show only “interesting” residues
                    name1 = f"{chr(65 + r1.chain.index)}:{r1.name}{r1.resSeq}"
                    name2 = f"{chr(65 + r2.chain.index)}:{r2.name}{r2.resSeq}"
                    print(f"{name1:<22} | {fp1[i]:>5.2f} | {fp2[i]:>5.2f}   ({name2})")

            summary["ligand"]["contact_fingerprint"] = [
                {
                    "residue": f"{chr(65 + res.chain.index)}:{res.name}{res.resSeq}",
                    args.label1: float(fp1[i]),
                    args.label2: float(fp2[i])
                }
                for i, res in enumerate(aligned_res1)
                if max(fp1[i], fp2[i]) > CONFIG["contact_fingerprint"]["min_occupancy"]
            ]
            summary["ligand"]["contact_fingerprint_summary"] = {
                "n_contacts_over_50pct": sum(
                    1 for i in range(len(aligned_res1))
                    if max(fp1[i], fp2[i]) > CONFIG["contact_fingerprint"]["min_occupancy"]
                )
            }
            summary["units"]["occupancy"] = "fraction"
            # ───── Ligand–Loop DCCM (Complete Structure, residue‑level) ─────
            print("\nLigand–Loop DCCM Analysis")
            print("=" * 30)

            # Set how many top residues to report
            N_TOP = CONFIG["dccm"]["ligand_loop_top_n"]

            # 1. loop residues present in prot1 but absent from prot2
            loop_residues = [
                res for res in prot1.topology.residues
                if res.is_protein and res.index not in {r.index for r in prot2.topology.residues}
            ]

            # 2. the ligand (as a residue object).  We assume one copy.
            ligand_res = next((res for res in t1.topology.residues
                               if res.name == ligand_resname), None)

            if not loop_residues:
                print(f"  Skipped: no unique loop residues in {args.label1}.")
            elif ligand_res is None:
                print("  Skipped: ligand residue not found.")
            else:
                all_residues = [ligand_res] + loop_residues     # ligand is index 0

                # ---------- helper functions ---------------------------------
                def residue_displacements(traj, residues):
                    """Return (n_frames × 3·n_res) matrix of COM displacements."""
                    n_f, n_r = traj.n_frames, len(residues)
                    com = np.zeros((n_f, n_r, 3), dtype=np.float32)
                    for i, res in enumerate(residues):
                        idx = [a.index for a in res.atoms if a.element.symbol != "H"]
                        com[:, i, :] = traj.atom_slice(idx).xyz.mean(axis=1)
                    disp = com - com.mean(axis=0, keepdims=True)
                    return disp.reshape(n_f, -1)

                def dccm_from_disp(disp):
                    """Return (n_res × n_res) DCCM from displacement matrix."""
                    n_res = disp.shape[1] // 3
                    dccm  = np.corrcoef(disp.T).reshape(n_res, 3, n_res, 3).mean(axis=(1, 3))
                    return dccm

                # 3. build displacement matrix & DCCM
                disp_mat = residue_displacements(t1, all_residues)
                dccm_res = dccm_from_disp(disp_mat)             # (n_res × n_res)

                # 4. extract ligand (row 0) ↔ loop block
                lig_loop_corr = dccm_res[0, 1:]                  # vector length = len(loop_residues)
                avg_corr = lig_loop_corr.mean()
                max_corr = np.abs(lig_loop_corr).max()
                idx_max  = np.argmax(np.abs(lig_loop_corr))

                # 5. save figure
                plt.figure(figsize=CONFIG["plot"]["figsize"]["ligand_loop"])
                plt.imshow(dccm_res, vmin=-1, vmax=1, cmap="bwr")
                plt.colorbar(label="Cross‑correlation")
                plt.title(f"Residue‑level DCCM  ({ligand_resname} ↔ new loops)")
                plt.tight_layout()
                plt.savefig(f"{args.out_prefix}_dccm_ligand_loop_residue.pdf", dpi=CONFIG["plot"]["dpi"])
                plt.close()
                print("  Ligand–Loop DCCM saved to PDF.")

                # 6. summary
                print("\nLigand–Loop DCCM Summary (Residue‑level):")
                print(f"  Average ligand–loop correlation : {avg_corr:.3f}")
                print(f"  Max |correlation|                : {max_corr:.3f}")
                print(f"  Most‑correlated loop residue     : {res_label(loop_residues[idx_max])}")

                # Top N positively correlated
                top_pos_idx = np.argsort(-lig_loop_corr)[:N_TOP]
                # Top N negatively correlated
                top_neg_idx = np.argsort(lig_loop_corr)[:N_TOP]

                print(f"  Top {N_TOP} loop residues by positive correlation:")
                for i in top_pos_idx:
                    print(f"    {res_label(loop_residues[i]):<12} |  corr = {lig_loop_corr[i]:.3f}")

                print(f"  Top {N_TOP} loop residues by negative correlation:")
                for i in top_neg_idx:
                    print(f"    {res_label(loop_residues[i]):<12} |  corr = {lig_loop_corr[i]:.3f}")

        lig_loop_stats = None
        if ligand_resname and loop_residues and ligand_res is not None:
            lig_loop_stats = {
                "avg_corr"         : float(avg_corr),
                "max_abs_corr"     : float(max_corr),
                "most_correlated"  : res_label(loop_residues[idx_max]),
                "top_positive"     : [
                    {"residue": res_label(loop_residues[i]),
                     "corr"   : float(lig_loop_corr[i])}
                    for i in top_pos_idx
                ],
                "top_negative"     : [
                    {"residue": res_label(loop_residues[i]),
                     "corr"   : float(lig_loop_corr[i])}
                    for i in top_neg_idx
                ]
            }
            summary["ligand"]["loop_dccm"] = lig_loop_stats
        
        # ───────── MD output / OpenMM system + energy validation (optional) ─────────
        out1_data = parse_out_file(args.out1) if args.out1 else None
        out2_data = parse_out_file(args.out2) if args.out2 else None

        if args.out1 or args.out2:
            print("\n" + "=" * 72)
            print("MD Output / OpenMM Validation")
            print("=" * 72)
            if args.out1:
                print_single_output_summary(args.label1, out1_data)
            if args.out2:
                print_single_output_summary(args.label2, out2_data)

            openmm_validation = {
                "delta_definition": f"{args.label1} - {args.label2}",
                "relative_difference_definition": (
                    "symmetric percent difference = 200*|x1-x2|/(|x1|+|x2|)"
                ),
                "input1": out1_data,
                "input2": out2_data,
                "comparison": None
            }

            serialized_xml_comparison = compare_serialized_openmm_xml(
                args.out1, args.out2, args.label1, args.label2,
                topology1=t1.topology, topology2=t2.topology
            )
            if serialized_xml_comparison is not None:
                openmm_validation["serialized_xml"] = serialized_xml_comparison
                print_serialized_openmm_xml_comparison(
                    args.label1, args.label2, serialized_xml_comparison,
                    topology1=t1.topology, topology2=t2.topology
                )

            if out1_data:
                compat1 = compatibility_out_summary(out1_data)
                if compat1:
                    summary["out1_summary"] = compat1
            if out2_data:
                compat2 = compatibility_out_summary(out2_data)
                if compat2:
                    summary["out2_summary"] = compat2

            if out1_data and out2_data:
                output_comparison = compare_openmm_outputs(out1_data, out2_data)
                openmm_validation["comparison"] = output_comparison

                print_system_comparison(args.label1, args.label2, output_comparison)
                print_representative_parameter_comparison(output_comparison)
                print_protocol_comparison(args.label1, args.label2, output_comparison)
                print_timing_comparison(args.label1, args.label2, output_comparison)
                print_energy_comparison(args.label1, args.label2, output_comparison)
                print_production_comparison(args.label1, args.label2, output_comparison)
                plot_production_validation(
                    out1_data, out2_data, output_comparison,
                    args.label1, args.label2, args.out_prefix
                )

            summary["openmm_validation"] = openmm_validation

        json_file = f"{args.out_prefix}_summary.json"
        with open(json_file, "w") as jfh:
            json.dump(summary, jfh, indent=2, ensure_ascii=False)
            print(f"INFO: Wrote summary JSON to {json_file}", flush=True)

    finally:
        for f in tmp_files:
            try:
                f.unlink(missing_ok=True)  # Python 3.8+
            except Exception as e:
                logging.warning("Could not delete temp file %s: %s", f, e)

if __name__=="__main__":
    main()
