#!/usr/bin/env python3
"""
Aggregate and triage analysisMD.py *_summary.json files.

This script combines the original cross-system comparison summary with the
newer DivCon/C++ vs Python/OpenMM engineering-validation triage.

It has three input modes:

1. Explicit JSON files / globs (backwards-compatible style):

       python summarize_MD_comparisons_v3.py */*_summary.json \
           --csv all_systems_summary.csv \
           --plots summary_plots.pdf

2. A PDB-ID list file.  Each JSON is expected at
   <root>/<PDBID>/<PDBID>_summary.json:

       python summarize_MD_comparisons_v3.py --root . --list anal.list

3. Automatic discovery of immediate four-character alphanumeric directories:

       python summarize_MD_comparisons_v3.py --root .

Unless overridden, the script writes:

    md_comparison_summary.csv
    md_comparison_report.txt
    md_comparison_review.list
    md_comparison_good.list
    md_comparison_good_limited.list
    md_comparison_pending.list
    md_comparison_errors.list

Use --plots FILE.pdf to additionally write a multi-page PDF with descriptive
cross-system distributions.

Triage categories are intentionally REVIEW / GOOD / GOOD_LIMITED / PENDING / ERROR rather
than PASS / FAIL.  The thresholds below are conservative smoke-test review
thresholds, not claims of ensemble equivalence or statistical convergence.

Dependencies:
    Python >= 3.8, numpy, pandas, matplotlib
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import re
import statistics
import sys
from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd


LOG_TAIL_BYTES = 256 * 1024

# These are deliberately high-confidence terminal failure signatures rather
# than an exhaustive list of Python/OpenMM exception names.  An unfinished log
# with no terminal signature remains PENDING because it may still be running.
FATAL_LOG_PATTERNS = (
    ("Python traceback", re.compile(r"Traceback \(most recent call last\):")),
    ("fatal Python error", re.compile(r"Fatal Python error:", re.IGNORECASE)),
    (
        "production exception",
        re.compile(r"^ERROR during production metrics at step ", re.MULTILINE),
    ),
    ("segmentation fault", re.compile(r"segmentation fault", re.IGNORECASE)),
    ("bus error", re.compile(r"\bbus error\b", re.IGNORECASE)),
    (
        "floating-point exception",
        re.compile(r"\bfloating point exception\b", re.IGNORECASE),
    ),
    ("aborted process", re.compile(r"\baborted(?: \(core dumped\))?\b", re.IGNORECASE)),
    ("core dump", re.compile(r"\bcore dumped\b", re.IGNORECASE)),
    (
        "uncaught C++ exception",
        re.compile(r"terminate called after throwing an instance", re.IGNORECASE),
    ),
    (
        "killed process",
        re.compile(
            r"^(?:\s*Killed\s*|.*:\s*line\s+\d+:\s+\d+\s+Killed\b.*)$",
            re.MULTILINE | re.IGNORECASE,
        ),
    ),
    (
        "out-of-memory termination",
        re.compile(r"(?:out of memory|oom-kill|killed process)", re.IGNORECASE),
    ),
)


THRESHOLDS = {
    # Hamiltonian/system invariants.
    "charge_abs_e": 1.0e-4,
    "charge_delta_e": 1.0e-4,
    "initial_box_max_delta_nm": 1.0e-3,      # 0.01 Angstrom
    "representative_charge_delta_e": 1.0e-5,
    "representative_sigma_delta_nm": 1.0e-6,
    "representative_epsilon_delta_kj_mol": 1.0e-5,
    "water_constraint_delta_nm": 1.0e-4,

    # Post-minimization Hamiltonian comparison.  Initial energies are diagnostic
    # only because the two workflows can enter OpenMM from different geometries.
    "postmin_total_sym_pct": 0.50,
    "postmin_nonbonded_sym_pct": 0.50,

    # Short stochastic production comparison.  These are triage limits, not
    # ensemble-equivalence criteria.
    "production_pe_mean_sym_pct": 0.75,
    "production_pe_max_sym_pct": 1.00,
    "production_temperature_mean_sym_pct": 5.0,
    "production_volume_mean_sym_pct": 5.0,
    "production_density_mean_sym_pct": 5.0,

    # Structural comparison.  Internal C-alpha distances are useful because
    # they are invariant to overall translation/rotation.
    "ca_distance_rmse_angstrom": 1.0,
    "ca_distance_pearson_min": 0.98,
    "rg_relative_difference_pct": 5.0,
    "rmsf_profile_rmse_angstrom": 0.50,

    # Secondary short-trajectory indicators: reported as notes only.
    "dssp_difference_fraction_report": 0.20,
    "hbond_jaccard_report": 0.50,
}


EXACT_SYSTEM_FIELDS = (
    "particles",
    "constraints",
    "nonbonded_particles",
    "nonbonded_exceptions",
    "water_particles",
    "nonbonded_method",
    "dispersion_correction",
    "switching_function",
)

NUMERIC_SYSTEM_FIELDS = (
    "cutoff_nm",
    "ewald_error_tolerance",
)

PROTOCOL_FIELDS = (
    "minimization_limit",
    "nvt_steps",
    "npt_steps",
    "production_steps",
    "production_report_interval_steps",
    "target_temperature_k",
    "target_pressure_bar",
    "production_duration_ps",
    "timestep_ps",
)

CORE_DIAGNOSTICS = (
    "charge_input1_e",
    "charge_input2_e",
    "postmin_total_sym_pct",
    "postmin_nonbonded_sym_pct",
    "production_pe_mean_sym_pct",
)

# Compatibility with the analysisMD.py version currently generating the running
# validation set.  That script computes Rg directly from the raw protein
# coordinates and computes C-alpha internal distances with periodic=False before
# explicitly making molecules whole/reimaging them.  For wrapped trajectories,
# those JSON structural metrics can therefore become tens of Angstroms different
# even when the Hamiltonian/thermodynamic validation is otherwise sound.
#
# Keep all structural values in the CSV/report for later diagnosis, but do NOT
# use them as hard REVIEW triggers until analysisMD.py is fixed and the affected
# summaries are regenerated.
CURRENT_ANALYSISMD_STRUCTURAL_ADVISORY_ONLY = True

# Metrics retained from the original summarizer.  These are descriptive and are
# not primary short-run validation criteria.
LEGACY_DELTA_BASES = (
    "rmsf_full_mean",
    "rmsf_aln_mean",
    "rg_mean",
    "bb_rmsd_mean",
    "lig_rmsd_mean",
    "lig_rmsf_mean",
)

VALIDATION_DISTRIBUTION_FIELDS = (
    "postmin_total_sym_pct",
    "postmin_nonbonded_sym_pct",
    "production_pe_mean_sym_pct",
    "production_pe_max_sym_pct",
    "production_temperature_mean_sym_pct",
    "production_volume_mean_sym_pct",
    "production_density_mean_sym_pct",
    "rmsf_profile_rmse_angstrom",
    "ca_distance_rmse_angstrom",
    "ca_distance_pearson_r",
    "ca_distance_max_abs_angstrom",
    "rg_relative_difference_pct",
    "charge_input1_e",
    "charge_input2_e",
    "charge_delta_e",
)


CSV_FIELDS = [
    # Identity / status.
    "system",
    "status",
    "json_path",
    "label1",
    "label2",
    "generated_on",

    # Basic trajectory/system sizes.
    "frames_input1",
    "frames_input2",
    "atoms_input1",
    "atoms_input2",
    "n_residues_1",
    "n_residues_2",

    # Original structural summary metrics.
    "rmsf_full_mean_1",
    "rmsf_full_mean_2",
    "rmsf_full_ks_p",
    "rmsf_full_mean_delta",
    "rmsf_full_mean_abs_delta",
    "rmsf_aln_mean_1",
    "rmsf_aln_mean_2",
    "rmsf_aln_ks_p",
    "rmsf_aln_mean_delta",
    "rmsf_aln_mean_abs_delta",
    "rg_mean_1",
    "rg_mean_2",
    "rg_ks_p",
    "rg_mean_delta",
    "rg_mean_abs_delta",
    "bb_rmsd_mean_1",
    "bb_rmsd_mean_2",
    "bb_rmsd_mean_delta",
    "bb_rmsd_mean_abs_delta",
    "hbonds_shared",
    "hbonds_excl_1",
    "hbonds_excl_2",
    "hbond_jaccard",
    "dssp_diff_pct",
    "dccm_mean_delta",
    "dccm_max_abs_delta",

    # Optional ligand metrics retained from the original summarizer.
    "lig_present",
    "lig_rmsd_mean_1",
    "lig_rmsd_mean_2",
    "lig_rmsd_mean_delta",
    "lig_rmsd_mean_abs_delta",
    "lig_rmsf_mean_1",
    "lig_rmsf_mean_2",
    "lig_rmsf_mean_delta",
    "lig_rmsf_mean_abs_delta",

    # Newer trajectory-validation metrics.
    "rmsf_profile_rmse_angstrom",
    "rmsf_profile_mae_angstrom",
    "rmsf_profile_max_abs_angstrom",
    "rmsf_profile_pearson_r",
    "rmsf_profile_spearman_rho",
    "rg_relative_difference_pct",
    "common_rmsd_input1_mean_angstrom",
    "common_rmsd_input2_mean_angstrom",
    "common_rmsd_mean_delta_angstrom",
    "common_rmsd_mean_abs_delta_angstrom",
    "ca_distance_rmse_angstrom",
    "ca_distance_mae_angstrom",
    "ca_distance_pearson_r",
    "ca_distance_spearman_rho",
    "ca_distance_max_abs_angstrom",

    # OpenMM system/Hamiltonian construction.
    "charge_input1_e",
    "charge_input2_e",
    "charge_delta_e",
    "charge_abs_delta_e",
    "initial_box_max_delta_nm",
    "system_mismatch_count",
    "protocol_mismatch_count",

    # Post-minimization energy decomposition.
    "postmin_total_sym_pct",
    "postmin_nonbonded_sym_pct",
    "postmin_bond_sym_pct",
    "postmin_angle_sym_pct",
    "postmin_torsion_sym_pct",

    # Production thermodynamics.
    "production_matched_steps",
    "production_pe_mean_sym_pct",
    "production_pe_max_sym_pct",
    "production_temperature_mean_sym_pct",
    "production_volume_mean_sym_pct",
    "production_density_mean_sym_pct",

    # Triage metadata.
    "missing_core_diagnostics",
    "review_categories",
    "review_reasons",
    "structural_advisories",
    "notes",
]


def nested(mapping: Any, *keys: str) -> Any:
    value = mapping
    for key in keys:
        if not isinstance(value, dict):
            return None
        value = value.get(key)
    return value


def finite_float(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def fmt(value: Any, digits: int = 4) -> str:
    value = finite_float(value)
    if value is None:
        return "-"
    return f"{value:.{digits}g}"


def add_reason(reasons: list[str], text: str) -> None:
    if text not in reasons:
        reasons.append(text)


def set_pair_delta(row: dict[str, Any], base: str) -> None:
    value1 = finite_float(row.get(f"{base}_1"))
    value2 = finite_float(row.get(f"{base}_2"))
    if value1 is None or value2 is None:
        row[f"{base}_delta"] = None
        row[f"{base}_abs_delta"] = None
        return
    delta = value1 - value2
    row[f"{base}_delta"] = delta
    row[f"{base}_abs_delta"] = abs(delta)


def compare_exact_system_fields(summary: dict, reasons: list[str]) -> int:
    comparison = nested(summary, "openmm_validation", "comparison", "system") or {}
    mismatches = 0

    for field in EXACT_SYSTEM_FIELDS:
        item = comparison.get(field)
        if isinstance(item, dict) and item.get("match") is False:
            mismatches += 1
            add_reason(
                reasons,
                f"system {field} mismatch ({item.get('input1')} vs {item.get('input2')})",
            )

    for field in NUMERIC_SYSTEM_FIELDS:
        item = comparison.get(field)
        difference = finite_float(nested(item, "absolute_difference"))
        if difference is not None and difference > 1.0e-12:
            mismatches += 1
            add_reason(reasons, f"system {field} differs by {difference:.3g}")

    return mismatches


def compare_protocol(summary: dict, reasons: list[str]) -> int:
    comparison = nested(summary, "openmm_validation", "comparison", "protocol") or {}
    mismatches = 0

    for field in PROTOCOL_FIELDS:
        item = comparison.get(field)
        if not isinstance(item, dict):
            continue
        if item.get("match") is False:
            mismatches += 1
            add_reason(
                reasons,
                f"protocol {field} mismatch ({item.get('input1')} vs {item.get('input2')})",
            )
            continue

        difference = finite_float(item.get("absolute_difference"))
        if difference is not None and difference > 1.0e-12:
            mismatches += 1
            add_reason(reasons, f"protocol {field} differs by {difference:.3g}")

    return mismatches


def representative_parameter_checks(summary: dict, reasons: list[str]) -> None:
    comparison = nested(
        summary, "openmm_validation", "comparison", "representative_particles"
    ) or {}

    for role, particle in comparison.items():
        if not isinstance(particle, dict):
            continue

        charge_delta = finite_float(
            nested(particle, "openmm_charge_e", "absolute_difference")
        )
        sigma_delta = finite_float(nested(particle, "sigma_nm", "absolute_difference"))
        epsilon_delta = finite_float(
            nested(particle, "epsilon_kj_per_mol", "absolute_difference")
        )

        if (
            charge_delta is not None
            and charge_delta > THRESHOLDS["representative_charge_delta_e"]
        ):
            add_reason(reasons, f"{role} representative charge delta={charge_delta:.3g} e")
        if (
            sigma_delta is not None
            and sigma_delta > THRESHOLDS["representative_sigma_delta_nm"]
        ):
            add_reason(reasons, f"{role} representative sigma delta={sigma_delta:.3g} nm")
        if (
            epsilon_delta is not None
            and epsilon_delta > THRESHOLDS["representative_epsilon_delta_kj_mol"]
        ):
            add_reason(
                reasons,
                f"{role} representative epsilon delta={epsilon_delta:.3g} kJ/mol",
            )

    constraint_delta = finite_float(
        nested(
            summary,
            "openmm_validation",
            "comparison",
            "representative_water_constraints",
            "max_absolute_distance_difference_nm",
        )
    )
    if (
        constraint_delta is not None
        and constraint_delta > THRESHOLDS["water_constraint_delta_nm"]
    ):
        add_reason(reasons, f"water constraint delta={constraint_delta:.3g} nm")


def extract_metrics(system: str, json_path: Path, summary: dict) -> dict[str, Any]:
    label1 = str(summary.get("label1") or "input1")
    label2 = str(summary.get("label2") or "input2")
    openmm = summary.get("openmm_validation") or {}
    input1 = openmm.get("input1") or {}
    input2 = openmm.get("input2") or {}
    comparison = openmm.get("comparison") or {}
    trajectory = summary.get("trajectory_validation") or {}

    row: dict[str, Any] = {field: None for field in CSV_FIELDS}
    row.update(
        {
            "system": system,
            "json_path": str(json_path),
            "label1": label1,
            "label2": label2,
            "generated_on": summary.get("generated_on"),
        }
    )

    # Basic trajectory identity and sizes.
    row["frames_input1"] = nested(summary, "n_frames", label1)
    row["frames_input2"] = nested(summary, "n_frames", label2)
    row["atoms_input1"] = nested(summary, "n_atoms", label1)
    row["atoms_input2"] = nested(summary, "n_atoms", label2)
    row["n_residues_1"] = nested(summary, "n_residues", label1)
    row["n_residues_2"] = nested(summary, "n_residues", label2)

    # Original structural summary metrics.
    row["rmsf_full_mean_1"] = finite_float(nested(summary, "rmsf", "full", "mean1"))
    row["rmsf_full_mean_2"] = finite_float(nested(summary, "rmsf", "full", "mean2"))
    row["rmsf_full_ks_p"] = finite_float(nested(summary, "rmsf", "full", "ks_p"))
    row["rmsf_aln_mean_1"] = finite_float(nested(summary, "rmsf", "aligned", "mean1"))
    row["rmsf_aln_mean_2"] = finite_float(nested(summary, "rmsf", "aligned", "mean2"))
    row["rmsf_aln_ks_p"] = finite_float(nested(summary, "rmsf", "aligned", "ks_p"))
    row["rg_mean_1"] = finite_float(nested(summary, "rg", "mean1"))
    row["rg_mean_2"] = finite_float(nested(summary, "rg", "mean2"))
    row["rg_ks_p"] = finite_float(nested(summary, "rg", "ks_p"))
    row["bb_rmsd_mean_1"] = finite_float(nested(summary, "rmsd_backbone", label1, "mean"))
    row["bb_rmsd_mean_2"] = finite_float(nested(summary, "rmsd_backbone", label2, "mean"))

    # Persistent H-bonds and secondary structure.
    hbonds = summary.get("hbonds") or {}
    row["hbonds_shared"] = hbonds.get("shared")
    row["hbonds_excl_1"] = hbonds.get("exclusive1")
    row["hbonds_excl_2"] = hbonds.get("exclusive2")
    try:
        n1 = int(hbonds.get("n1"))
        n2 = int(hbonds.get("n2"))
        shared = int(hbonds.get("shared"))
        union = n1 + n2 - shared
        row["hbond_jaccard"] = float(shared / union) if union > 0 else 1.0
    except (TypeError, ValueError):
        row["hbond_jaccard"] = None

    row["dssp_diff_pct"] = finite_float(summary.get("dssp_diff_pct"))
    row["dccm_mean_delta"] = finite_float(nested(summary, "dccm_overall", "mean_delta"))
    row["dccm_max_abs_delta"] = finite_float(
        nested(summary, "dccm_overall", "max_abs_delta")
    )

    # Optional ligand metrics retained from the original summarizer.
    ligand = summary.get("ligand")
    row["lig_present"] = isinstance(ligand, dict)
    row["lig_rmsd_mean_1"] = finite_float(nested(ligand, "rmsd", label1, "mean"))
    row["lig_rmsd_mean_2"] = finite_float(nested(ligand, "rmsd", label2, "mean"))
    row["lig_rmsf_mean_1"] = finite_float(nested(ligand, "rmsf", label1, "mean"))
    row["lig_rmsf_mean_2"] = finite_float(nested(ligand, "rmsf", label2, "mean"))

    for base in LEGACY_DELTA_BASES:
        set_pair_delta(row, base)

    # Newer trajectory-validation metrics.  Residue-specific top-difference labels
    # are deliberately not used for triage; only numeric metrics are consumed.
    rmsf_profile = trajectory.get("rmsf_profile") or {}
    row["rmsf_profile_rmse_angstrom"] = finite_float(rmsf_profile.get("rmse"))
    row["rmsf_profile_mae_angstrom"] = finite_float(rmsf_profile.get("mae"))
    row["rmsf_profile_max_abs_angstrom"] = finite_float(
        rmsf_profile.get("max_abs_difference")
    )
    row["rmsf_profile_pearson_r"] = finite_float(rmsf_profile.get("pearson_r"))
    row["rmsf_profile_spearman_rho"] = finite_float(rmsf_profile.get("spearman_rho"))

    row["rg_relative_difference_pct"] = finite_float(
        nested(trajectory, "rg_stability", "relative_mean_difference_percent")
    )

    common = trajectory.get("common_reference_rmsd") or {}
    common1 = finite_float(nested(common, label1, "mean"))
    common2 = finite_float(nested(common, label2, "mean"))
    row["common_rmsd_input1_mean_angstrom"] = common1
    row["common_rmsd_input2_mean_angstrom"] = common2
    if common1 is not None and common2 is not None:
        delta = common1 - common2
        row["common_rmsd_mean_delta_angstrom"] = delta
        row["common_rmsd_mean_abs_delta_angstrom"] = abs(delta)

    ca = trajectory.get("ca_internal_distances") or {}
    row["ca_distance_rmse_angstrom"] = finite_float(ca.get("rmse"))
    row["ca_distance_mae_angstrom"] = finite_float(ca.get("mae"))
    row["ca_distance_pearson_r"] = finite_float(ca.get("pearson_r"))
    row["ca_distance_spearman_rho"] = finite_float(ca.get("spearman_rho"))
    row["ca_distance_max_abs_angstrom"] = finite_float(ca.get("max_abs_difference"))

    # System/Hamiltonian construction.
    row["charge_input1_e"] = finite_float(nested(input1, "system", "total_particle_charge_e"))
    row["charge_input2_e"] = finite_float(nested(input2, "system", "total_particle_charge_e"))
    if row["charge_input1_e"] is not None and row["charge_input2_e"] is not None:
        delta = row["charge_input1_e"] - row["charge_input2_e"]
        row["charge_delta_e"] = delta
        row["charge_abs_delta_e"] = abs(delta)

    row["initial_box_max_delta_nm"] = finite_float(
        nested(comparison, "system", "initial_box_nm", "max_absolute_difference")
    )

    # Post-minimization energy decomposition.  Initial energies remain descriptive
    # in analysisMD and are intentionally not triage criteria here.
    postmin = nested(comparison, "energy_decomposition_kj_per_mol", "post_minimization") or {}
    for component, column in (
        ("total", "postmin_total_sym_pct"),
        ("nonbonded", "postmin_nonbonded_sym_pct"),
        ("bond", "postmin_bond_sym_pct"),
        ("angle", "postmin_angle_sym_pct"),
        ("torsion", "postmin_torsion_sym_pct"),
    ):
        row[column] = finite_float(
            nested(postmin, component, "symmetric_relative_difference_percent")
        )

    # Production thermodynamics.
    production = comparison.get("production") or {}
    row["production_matched_steps"] = production.get("n_matched_steps")
    observables = production.get("observables") or {}
    production_columns = {
        "potential_energy_kj_per_mol": (
            "production_pe_mean_sym_pct",
            "production_pe_max_sym_pct",
        ),
        "temperature_k": ("production_temperature_mean_sym_pct", None),
        "volume_nm3": ("production_volume_mean_sym_pct", None),
        "density_g_per_ml": ("production_density_mean_sym_pct", None),
    }
    for observable, (mean_column, max_column) in production_columns.items():
        stats = observables.get(observable) or {}
        row[mean_column] = finite_float(
            stats.get("mean_absolute_symmetric_relative_difference_percent")
        )
        if max_column:
            row[max_column] = finite_float(
                stats.get("max_absolute_symmetric_relative_difference_percent")
            )

    return row


def classify(summary: dict, row: dict[str, Any]) -> None:
    reasons: list[str] = []
    categories: set[str] = set()
    structural_advisories: list[str] = []
    notes: list[str] = []

    def review(category: str, text: str) -> None:
        categories.add(category)
        add_reason(reasons, text)

    frames1 = row.get("frames_input1")
    frames2 = row.get("frames_input2")
    atoms1 = row.get("atoms_input1")
    atoms2 = row.get("atoms_input2")

    if frames1 is not None and frames2 is not None and frames1 != frames2:
        review("SYSTEM", f"trajectory frame mismatch ({frames1} vs {frames2})")
    if atoms1 is not None and atoms2 is not None and atoms1 != atoms2:
        review("SYSTEM", f"trajectory atom mismatch ({atoms1} vs {atoms2})")

    before = len(reasons)
    row["system_mismatch_count"] = compare_exact_system_fields(summary, reasons)
    if len(reasons) > before:
        categories.add("SYSTEM")

    before = len(reasons)
    row["protocol_mismatch_count"] = compare_protocol(summary, reasons)
    if len(reasons) > before:
        categories.add("PROTOCOL")

    before = len(reasons)
    representative_parameter_checks(summary, reasons)
    if len(reasons) > before:
        categories.add("SYSTEM")

    q1 = finite_float(row.get("charge_input1_e"))
    q2 = finite_float(row.get("charge_input2_e"))
    if q1 is not None and abs(q1) > THRESHOLDS["charge_abs_e"]:
        review("CHARGE", f"{row['label1']} charge={q1:+.6g} e")
    if q2 is not None and abs(q2) > THRESHOLDS["charge_abs_e"]:
        review("CHARGE", f"{row['label2']} charge={q2:+.6g} e")
    if q1 is not None and q2 is not None:
        qdelta = abs(q1 - q2)
        if qdelta > THRESHOLDS["charge_delta_e"]:
            review("CHARGE", f"charge disagreement={qdelta:.6g} e")

    box_delta = finite_float(row.get("initial_box_max_delta_nm"))
    if box_delta is not None and box_delta > THRESHOLDS["initial_box_max_delta_nm"]:
        review("SYSTEM", f"initial box max delta={box_delta:.4g} nm")

    post_total = finite_float(row.get("postmin_total_sym_pct"))
    if post_total is not None and post_total > THRESHOLDS["postmin_total_sym_pct"]:
        review("ENERGY", f"post-min total energy={post_total:.3f}%")

    post_nb = finite_float(row.get("postmin_nonbonded_sym_pct"))
    if post_nb is not None and post_nb > THRESHOLDS["postmin_nonbonded_sym_pct"]:
        review("ENERGY", f"post-min nonbonded energy={post_nb:.3f}%")

    pe_mean = finite_float(row.get("production_pe_mean_sym_pct"))
    if pe_mean is not None and pe_mean > THRESHOLDS["production_pe_mean_sym_pct"]:
        review("ENERGY", f"production PE mean={pe_mean:.3f}%")

    pe_max = finite_float(row.get("production_pe_max_sym_pct"))
    if pe_max is not None and pe_max > THRESHOLDS["production_pe_max_sym_pct"]:
        review("ENERGY", f"production PE max={pe_max:.3f}%")

    for column, threshold_key, label in (
        ("production_temperature_mean_sym_pct", "production_temperature_mean_sym_pct", "production temperature mean"),
        ("production_volume_mean_sym_pct", "production_volume_mean_sym_pct", "production volume mean"),
        ("production_density_mean_sym_pct", "production_density_mean_sym_pct", "production density mean"),
    ):
        value = finite_float(row.get(column))
        if value is not None and value > THRESHOLDS[threshold_key]:
            review("THERMODYNAMICS", f"{label}={value:.3f}%")

    # Current analysisMD.py compatibility: preserve structural metrics, but treat
    # threshold excursions as advisory rather than REVIEW.  The current producer
    # can compare wrapped coordinates without first making protein molecules whole.
    ca_rmse = finite_float(row.get("ca_distance_rmse_angstrom"))
    if ca_rmse is not None and ca_rmse > THRESHOLDS["ca_distance_rmse_angstrom"]:
        qualifier = "SEVERE " if ca_rmse > 3.0 else ""
        structural_advisories.append(f"{qualifier}C-alpha distance RMSE={ca_rmse:.3f} A")

    ca_r = finite_float(row.get("ca_distance_pearson_r"))
    if ca_r is not None and ca_r < THRESHOLDS["ca_distance_pearson_min"]:
        structural_advisories.append(f"C-alpha distance Pearson r={ca_r:.4f}")

    rg_relative = finite_float(row.get("rg_relative_difference_pct"))
    if rg_relative is not None and abs(rg_relative) > THRESHOLDS["rg_relative_difference_pct"]:
        qualifier = "SEVERE " if abs(rg_relative) > 20.0 else ""
        structural_advisories.append(f"{qualifier}Rg mean delta={rg_relative:+.3f}%")

    rmsf_rmse = finite_float(row.get("rmsf_profile_rmse_angstrom"))
    if rmsf_rmse is not None and rmsf_rmse > THRESHOLDS["rmsf_profile_rmse_angstrom"]:
        structural_advisories.append(f"RMSF profile RMSE={rmsf_rmse:.3f} A")

    hbond_jaccard = finite_float(row.get("hbond_jaccard"))
    if hbond_jaccard is not None and hbond_jaccard < THRESHOLDS["hbond_jaccard_report"]:
        notes.append(f"low H-bond Jaccard={hbond_jaccard:.3f}")

    dssp = finite_float(row.get("dssp_diff_pct"))
    if dssp is not None and dssp > THRESHOLDS["dssp_difference_fraction_report"]:
        notes.append(f"DSSP disagreement={100.0 * dssp:.1f}%")

    missing = [name for name in CORE_DIAGNOSTICS if row.get(name) is None]
    row["missing_core_diagnostics"] = ";".join(missing)

    if reasons:
        row["status"] = "REVIEW"
    elif missing:
        row["status"] = "GOOD_LIMITED"
    else:
        row["status"] = "GOOD"

    row["review_categories"] = ";".join(sorted(categories))
    row["review_reasons"] = "; ".join(reasons)
    row["structural_advisories"] = "; ".join(structural_advisories)
    row["notes"] = "; ".join(notes)


def pending_row(
    system: str,
    json_path: Path,
    message: str = "summary JSON not available yet",
) -> dict[str, Any]:
    row = {field: None for field in CSV_FIELDS}
    row.update(
        {
            "system": system,
            "status": "PENDING",
            "json_path": str(json_path),
            "review_categories": "",
            "review_reasons": message,
            "structural_advisories": "",
            "notes": "",
            "missing_core_diagnostics": ";".join(CORE_DIAGNOSTICS),
        }
    )
    return row


def error_row(system: str, json_path: Path, message: str) -> dict[str, Any]:
    row = {field: None for field in CSV_FIELDS}
    row.update(
        {
            "system": system,
            "status": "ERROR",
            "json_path": str(json_path),
            "review_categories": "ERROR",
            "review_reasons": message,
            "structural_advisories": "",
            "notes": "",
            "missing_core_diagnostics": ";".join(CORE_DIAGNOSTICS),
        }
    )
    return row


def read_log_tail(path: Path, max_bytes: int = LOG_TAIL_BYTES) -> str:
    """Read only the tail of a potentially large screenout file."""
    with path.open("rb") as handle:
        handle.seek(0, 2)
        size = handle.tell()
        handle.seek(max(0, size - max_bytes))
        return handle.read().decode("utf-8", errors="replace")


def inspect_terminal_log(path: Path, success_markers: tuple[str, ...]) -> tuple[str, str | None]:
    """Best-effort terminal-state detection from a screenout file.

    SUCCESS requires all supplied completion markers.  ERROR requires a strong
    fatal signature.  Anything else is INCOMPLETE because the process may still
    be running and a screenout file alone cannot distinguish that safely.
    """
    if not path.is_file():
        return "MISSING", None

    try:
        text = read_log_tail(path)
    except OSError as exc:
        return "ERROR", f"cannot read {path}: {exc}"

    success_position = -1
    if success_markers and all(marker in text for marker in success_markers):
        success_position = max(text.rfind(marker) for marker in success_markers)

    failure_position = -1
    failure_label = None
    for label, pattern in FATAL_LOG_PATTERNS:
        matches = list(pattern.finditer(text))
        if matches and matches[-1].start() > failure_position:
            failure_position = matches[-1].start()
            failure_label = label

    if failure_position > success_position:
        return "ERROR", failure_label
    if success_position >= 0:
        return "SUCCESS", None
    return "INCOMPLETE", None


def classify_missing_summary(system: str, json_path: Path) -> dict[str, Any]:
    """Classify a missing analysis summary from available run screenouts."""
    run_dir = json_path.parent
    log_specs = (
        (
            "Python MD",
            run_dir / "python" / "md.screenout",
            ("Simulation summary:", "  wall_time_s:"),
        ),
        (
            "C++ MD",
            run_dir / "cpp" / "md.screenout",
            ("Job Complete", "Total Computation Time (Seconds):"),
        ),
        (
            "analysis",
            run_dir / "analysis.screenout",
            ("INFO: Wrote summary JSON to",),
        ),
    )

    states = {}
    for label, path, success_markers in log_specs:
        state, detail = inspect_terminal_log(path, success_markers)
        states[label] = (state, detail, path)
        if state == "ERROR":
            reason = detail or "terminal failure signature"
            return error_row(
                system,
                json_path,
                f"{label} failed: {reason} detected in {path}",
            )

    analysis_state, _, analysis_path = states["analysis"]
    if analysis_state == "SUCCESS":
        return error_row(
            system,
            json_path,
            f"analysis reports successful summary write but JSON is missing: {analysis_path}",
        )

    incomplete = [
        label for label in ("Python MD", "C++ MD", "analysis")
        if states[label][0] == "INCOMPLETE"
    ]
    completed = [
        label for label in ("Python MD", "C++ MD")
        if states[label][0] == "SUCCESS"
    ]

    if incomplete:
        message = (
            "summary JSON not available; no terminal success/failure marker yet in "
            + ", ".join(incomplete)
        )
    elif completed:
        message = (
            "summary JSON not available yet; completed run logs: "
            + ", ".join(completed)
        )
    else:
        message = "summary JSON not available yet; no terminal run evidence found"

    return pending_row(system, json_path, message)


def infer_system_name(json_path: Path, summary: dict | None = None) -> str:
    if isinstance(summary, dict):
        out_prefix = summary.get("out_prefix")
        if out_prefix:
            return Path(str(out_prefix)).stem

    stem = json_path.name
    suffix = "_summary.json"
    if stem.endswith(suffix):
        return stem[:-len(suffix)]
    return json_path.stem


def analyze_one(json_path: Path, system_hint: str | None = None) -> dict[str, Any]:
    system = system_hint or infer_system_name(json_path)
    if not json_path.is_file():
        return classify_missing_summary(system, json_path)

    try:
        with json_path.open() as handle:
            summary = json.load(handle)
    except Exception as exc:
        return error_row(system, json_path, f"cannot parse JSON: {exc}")

    if not isinstance(summary, dict):
        return error_row(system, json_path, "JSON root is not an object")

    system = system_hint or infer_system_name(json_path, summary)
    try:
        row = extract_metrics(system, json_path, summary)
        classify(summary, row)
        return row
    except Exception as exc:
        return error_row(system, json_path, f"analysis error: {exc}")


def read_pdbids(path: Path) -> list[str]:
    pdbids: list[str] = []
    seen: set[str] = set()
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            pdbid = line.split()[0]
            if pdbid not in seen:
                pdbids.append(pdbid)
                seen.add(pdbid)
    return pdbids


def discover_pdbids(root: Path) -> list[str]:
    pattern = re.compile(r"^[A-Za-z0-9]{4}$")
    return sorted(
        path.name
        for path in root.iterdir()
        if path.is_dir() and pattern.fullmatch(path.name)
    )


def expand_json_patterns(patterns: Iterable[str]) -> list[Path]:
    paths: list[Path] = []
    seen: set[Path] = set()

    for pattern in patterns:
        candidate = Path(pattern)
        matches: list[str]
        if candidate.is_file():
            matches = [str(candidate)]
        else:
            matches = glob.glob(pattern)

        if not matches:
            print(f"[WARN] No JSON files matched: {pattern}", file=sys.stderr)
            continue

        for match in matches:
            path = Path(match)
            if path not in seen:
                paths.append(path)
                seen.add(path)

    return paths


def resolve_inputs(args: argparse.Namespace) -> list[tuple[str | None, Path]]:
    root = args.root.resolve()

    if args.json_files and args.list_file:
        raise ValueError("positional JSON inputs and --list are mutually exclusive")

    if args.json_files:
        return [(None, path) for path in expand_json_patterns(args.json_files)]

    if not root.is_dir():
        raise ValueError(f"root directory does not exist: {root}")

    if args.list_file:
        if not args.list_file.is_file():
            raise ValueError(f"list file does not exist: {args.list_file}")
        pdbids = read_pdbids(args.list_file)
    else:
        pdbids = discover_pdbids(root)

    return [
        (pdbid, root / pdbid / f"{pdbid}_summary.json")
        for pdbid in pdbids
    ]


def clean_numeric_values(rows: list[dict[str, Any]], field: str) -> list[float]:
    values = []
    for row in rows:
        if row.get("status") == "ERROR":
            continue
        value = finite_float(row.get(field))
        if value is not None:
            values.append(value)
    return values


def distribution(values: list[float]) -> dict[str, float] | None:
    if not values:
        return None
    arr = np.asarray(values, dtype=float)
    result = {
        "n": float(len(arr)),
        "min": float(np.min(arr)),
        "mean": float(np.mean(arr)),
        "median": float(np.median(arr)),
        "max": float(np.max(arr)),
        "p90": float(np.percentile(arr, 90)),
        "p95": float(np.percentile(arr, 95)),
        "p99": float(np.percentile(arr, 99)),
    }
    if len(arr) >= 2:
        result["std"] = float(np.std(arr, ddof=1))
    else:
        result["std"] = 0.0
    return result


def signed_delta_fields(rows: list[dict[str, Any]]) -> list[str]:
    if not rows:
        return []
    fields = []
    for base in LEGACY_DELTA_BASES:
        field = f"{base}_delta"
        if clean_numeric_values(rows, field):
            fields.append(field)
    if clean_numeric_values(rows, "dccm_mean_delta"):
        fields.append("dccm_mean_delta")
    return fields


def stats_table(rows: list[dict[str, Any]], fields: Iterable[str]) -> pd.DataFrame:
    records = []
    for field in fields:
        stats = distribution(clean_numeric_values(rows, field))
        if stats is None:
            continue
        records.append({"metric": field, **stats})
    if not records:
        return pd.DataFrame(
            columns=["metric", "n", "min", "mean", "median", "std", "p90", "p95", "p99", "max"]
        )
    return pd.DataFrame(records)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    frame = pd.DataFrame(rows)
    frame = frame.reindex(columns=CSV_FIELDS)
    frame.to_csv(path, index=False)


def write_list(path: Path, rows: list[dict[str, Any]], statuses: set[str]) -> None:
    with path.open("w") as handle:
        for row in rows:
            if row.get("status") in statuses:
                handle.write(f"{row['system']}\n")


def write_report(path: Path, rows: list[dict[str, Any]], input_description: str) -> None:
    counts = {status: 0 for status in ("REVIEW", "GOOD", "GOOD_LIMITED", "PENDING", "ERROR")}
    for row in rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1

    legacy_stats = stats_table(rows, signed_delta_fields(rows))
    validation_stats = stats_table(rows, VALIDATION_DISTRIBUTION_FIELDS)

    with path.open("w") as handle:
        handle.write("MD comparison / OpenMM validation summary\n")
        handle.write("=========================================\n\n")
        handle.write(f"Input: {input_description}\n")
        handle.write(f"Structures considered: {len(rows)}\n")
        handle.write(f"REVIEW       : {counts.get('REVIEW', 0)}\n")
        handle.write(f"GOOD         : {counts.get('GOOD', 0)}\n")
        handle.write(f"GOOD_LIMITED : {counts.get('GOOD_LIMITED', 0)}\n")
        handle.write(f"PENDING      : {counts.get('PENDING', 0)}\n")
        handle.write(f"ERROR        : {counts.get('ERROR', 0)}\n\n")

        handle.write("Interpretation\n")
        handle.write("--------------\n")
        handle.write("REVIEW: one or more conservative engineering-review thresholds were exceeded.\n")
        handle.write("GOOD: no review threshold was exceeded and all core diagnostics were present.\n")
        handle.write("GOOD_LIMITED: available metrics did not trigger review, but core diagnostics were absent.\n")
        handle.write(
            "PENDING: summary JSON does not exist and available screenouts do not show a "
            "terminal failure; the run may still be active.\n"
        )
        handle.write(
            "ERROR: a terminal failure was detected in a screenout, or an existing summary "
            "JSON was unreadable or could not be processed.\n"
        )
        handle.write("These labels are triage categories, not proof of ensemble equivalence or convergence.\n\n")
        handle.write("Current analysisMD compatibility\n")
        handle.write("------------------------------\n")
        handle.write(
            "Structural trajectory metrics are retained as advisories only for this running batch. "
            "The current analysisMD.py computes Rg and internal C-alpha distances from raw trajectory "
            "coordinates without an explicit make-whole/reimaging step; wrapped proteins can therefore "
            "create very large apparent structural differences. These metrics do not trigger REVIEW here.\n\n"
        )

        handle.write("Primary triage thresholds\n")
        handle.write("-------------------------\n")
        for key, value in THRESHOLDS.items():
            handle.write(f"{key}: {value}\n")
        handle.write("\n")

        handle.write("Structures requiring review\n")
        handle.write("---------------------------\n")
        review_rows = [row for row in rows if row["status"] == "REVIEW"]
        if not review_rows:
            handle.write("None\n")
        else:
            for row in review_rows:
                handle.write(f"{row['system']}: {row['review_reasons']}\n")
        handle.write("\n")

        handle.write("Structural advisories from current analysisMD output\n")
        handle.write("---------------------------------------------------\n")
        advisory_rows = [row for row in rows if row.get("structural_advisories")]
        if not advisory_rows:
            handle.write("None\n")
        else:
            for row in advisory_rows:
                handle.write(f"{row['system']}: {row['structural_advisories']}\n")
        handle.write("\n")

        handle.write("Good but limited diagnostics\n")
        handle.write("----------------------------\n")
        limited_rows = [row for row in rows if row["status"] == "GOOD_LIMITED"]
        if not limited_rows:
            handle.write("None\n")
        else:
            for row in limited_rows:
                handle.write(
                    f"{row['system']}: missing {row['missing_core_diagnostics']}\n"
                )
        handle.write("\n")

        handle.write("Pending summaries\n")
        handle.write("-----------------\n")
        pending_rows = [row for row in rows if row["status"] == "PENDING"]
        if not pending_rows:
            handle.write("None\n")
        else:
            for row in pending_rows:
                handle.write(f"{row['system']}\n")
        handle.write("\n")

        handle.write("Errors / failed runs\n")
        handle.write("--------------------\n")
        error_rows = [row for row in rows if row["status"] == "ERROR"]
        if not error_rows:
            handle.write("None\n")
        else:
            for row in error_rows:
                handle.write(f"{row['system']}: {row['review_reasons']}\n")
        handle.write("\n")

        handle.write("Cross-system signed structural deltas\n")
        handle.write("-------------------------------------\n")
        handle.write(
            "These retain the original summarizer's label1-label2 signed-delta view. "
            "Positive and negative values can cancel in the mean.\n"
        )
        if legacy_stats.empty:
            handle.write("No data\n")
        else:
            handle.write(legacy_stats.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
            handle.write("\n")
        handle.write("\n")

        handle.write("Validation metric distributions\n")
        handle.write("-------------------------------\n")
        if validation_stats.empty:
            handle.write("No data\n")
        else:
            handle.write(validation_stats.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
            handle.write("\n")
        handle.write("\n")

        handle.write("Compact validation table\n")
        handle.write("------------------------\n")
        header = (
            f"{'PDB':<8} {'STATUS':<13} {'postTot%':>9} {'postNB%':>9} "
            f"{'prodPE%':>9} {'CA_RMSE':>9} {'CA_r':>7} {'Rg%':>8} "
            f"{'q1':>10} {'q2':>10}\n"
        )
        handle.write(header)
        handle.write("-" * len(header.rstrip()) + "\n")
        for row in rows:
            handle.write(
                f"{str(row['system']):<8} {str(row['status']):<13} "
                f"{fmt(row.get('postmin_total_sym_pct')):>9} "
                f"{fmt(row.get('postmin_nonbonded_sym_pct')):>9} "
                f"{fmt(row.get('production_pe_mean_sym_pct')):>9} "
                f"{fmt(row.get('ca_distance_rmse_angstrom')):>9} "
                f"{fmt(row.get('ca_distance_pearson_r')):>7} "
                f"{fmt(row.get('rg_relative_difference_pct')):>8} "
                f"{fmt(row.get('charge_input1_e')):>10} "
                f"{fmt(row.get('charge_input2_e')):>10}\n"
            )


def finite_dataframe_columns(df: pd.DataFrame, columns: Iterable[str]) -> list[str]:
    result = []
    for column in columns:
        if column not in df.columns:
            continue
        numeric = pd.to_numeric(df[column], errors="coerce")
        if np.isfinite(numeric.to_numpy(dtype=float)).any():
            result.append(column)
    return result


def write_plots(path: Path, rows: list[dict[str, Any]]) -> None:
    df = pd.DataFrame(rows)
    if df.empty:
        return

    with PdfPages(path) as pdf:
        legacy_columns = finite_dataframe_columns(
            df,
            [
                "rmsf_full_mean_delta",
                "rmsf_aln_mean_delta",
                "rg_mean_delta",
                "bb_rmsd_mean_delta",
                "dccm_mean_delta",
            ],
        )
        if legacy_columns:
            fig, ax = plt.subplots(figsize=(8, 5))
            df[legacy_columns].apply(pd.to_numeric, errors="coerce").boxplot(ax=ax)
            ax.axhline(0.0, ls="--", lw=1, zorder=0)
            ax.set_title("Signed label1 - label2 structural differences")
            ax.set_ylabel("Difference (mixed native units)")
            ax.tick_params(axis="x", rotation=30)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

        percent_columns = finite_dataframe_columns(
            df,
            [
                "postmin_total_sym_pct",
                "postmin_nonbonded_sym_pct",
                "production_pe_mean_sym_pct",
                "production_temperature_mean_sym_pct",
                "production_volume_mean_sym_pct",
                "production_density_mean_sym_pct",
            ],
        )
        if percent_columns:
            fig, ax = plt.subplots(figsize=(8, 5))
            df[percent_columns].apply(pd.to_numeric, errors="coerce").boxplot(ax=ax)
            ax.set_title("OpenMM validation relative differences")
            ax.set_ylabel("Symmetric absolute difference (%)")
            ax.tick_params(axis="x", rotation=30)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

        structural_columns = finite_dataframe_columns(
            df,
            [
                "rmsf_profile_rmse_angstrom",
                "ca_distance_rmse_angstrom",
                "common_rmsd_mean_abs_delta_angstrom",
            ],
        )
        if structural_columns:
            fig, ax = plt.subplots(figsize=(8, 5))
            df[structural_columns].apply(pd.to_numeric, errors="coerce").boxplot(ax=ax)
            ax.set_title("Trajectory-validation structural differences")
            ax.set_ylabel("Angstrom")
            ax.tick_params(axis="x", rotation=30)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def print_console_summary(rows: list[dict[str, Any]], outputs: dict[str, Path]) -> None:
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1

    print(f"Structures: {len(rows)}")
    print(
        "Status: "
        f"REVIEW={counts.get('REVIEW', 0)}  "
        f"GOOD={counts.get('GOOD', 0)}  "
        f"GOOD_LIMITED={counts.get('GOOD_LIMITED', 0)}  "
        f"PENDING={counts.get('PENDING', 0)}  "
        f"ERROR={counts.get('ERROR', 0)}"
    )
    print()

    review_rows = [row for row in rows if row["status"] == "REVIEW"]
    if review_rows:
        print("Manual review candidates:")
        for row in review_rows:
            print(f"  {row['system']}: {row['review_reasons']}")
        print()

    advisory_rows = [row for row in rows if row.get("structural_advisories")]
    if advisory_rows:
        print(
            f"Structural advisories (not REVIEW triggers with current analysisMD.py): "
            f"{len(advisory_rows)}"
        )
        print()

    delta_stats = stats_table(rows, signed_delta_fields(rows))
    if not delta_stats.empty:
        print("Cross-system signed deltas:")
        print(delta_stats.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
        print()

    print("Wrote:")
    for label, path in outputs.items():
        print(f"  {label:<13}: {path}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate analysisMD *_summary.json files, report cross-system metrics, "
            "and flag conservative manual-review candidates."
        )
    )
    parser.add_argument(
        "json_files",
        nargs="*",
        help=(
            "Optional summary JSON files or glob patterns. If omitted, use --list "
            "or discover four-character directories under --root."
        ),
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("."),
        help="Parent directory containing PDBID directories (default: current directory).",
    )
    parser.add_argument(
        "--list",
        dest="list_file",
        type=Path,
        help="Optional file containing PDB IDs; blank/comment lines are ignored.",
    )
    parser.add_argument(
        "--prefix",
        default="md_comparison",
        help="Prefix for generated report/list/CSV files (default: md_comparison).",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        metavar="OUT.csv",
        help="Override the default <prefix>_summary.csv path.",
    )
    parser.add_argument(
        "--plots",
        type=Path,
        metavar="OUT.pdf",
        help="Optionally write cross-system distribution plots to a multi-page PDF.",
    )
    args = parser.parse_args(argv)

    try:
        input_pairs = resolve_inputs(args)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if not input_pairs:
        print("ERROR: no JSON files or PDB IDs found to process.", file=sys.stderr)
        return 1

    rows = [analyze_one(path, system_hint) for system_hint, path in input_pairs]

    # Review candidates first, then limited/good results, then data-flow errors.
    rank = {"REVIEW": 0, "GOOD_LIMITED": 1, "GOOD": 2, "PENDING": 3, "ERROR": 4}
    rows.sort(key=lambda row: (rank.get(str(row.get("status")), 99), str(row["system"])))

    prefix = Path(args.prefix)
    csv_path = args.csv or prefix.with_name(prefix.name + "_summary.csv")
    report_path = prefix.with_name(prefix.name + "_report.txt")
    review_path = prefix.with_name(prefix.name + "_review.list")
    good_path = prefix.with_name(prefix.name + "_good.list")
    limited_path = prefix.with_name(prefix.name + "_good_limited.list")
    pending_path = prefix.with_name(prefix.name + "_pending.list")
    error_path = prefix.with_name(prefix.name + "_errors.list")

    for output_path in (csv_path, report_path, review_path, good_path, limited_path, pending_path, error_path):
        output_path.parent.mkdir(parents=True, exist_ok=True)
    if args.plots:
        args.plots.parent.mkdir(parents=True, exist_ok=True)

    write_csv(csv_path, rows)

    if args.json_files:
        input_description = "explicit JSON files/globs"
    elif args.list_file:
        input_description = f"PDB list {args.list_file} under root {args.root.resolve()}"
    else:
        input_description = f"auto-discovered PDB directories under {args.root.resolve()}"

    write_report(report_path, rows, input_description)
    write_list(review_path, rows, {"REVIEW"})
    write_list(good_path, rows, {"GOOD"})
    write_list(limited_path, rows, {"GOOD_LIMITED"})
    write_list(pending_path, rows, {"PENDING"})
    write_list(error_path, rows, {"ERROR"})

    outputs = {
        "CSV": csv_path,
        "report": report_path,
        "review list": review_path,
        "good list": good_path,
        "limited list": limited_path,
        "pending list": pending_path,
        "error list": error_path,
    }

    if args.plots:
        write_plots(args.plots, rows)
        outputs["plots"] = args.plots

    print_console_summary(rows, outputs)

    # A missing JSON is PENDING unless available screenouts contain strong terminal
    # failure evidence. Malformed/unreadable existing JSON remains ERROR. REVIEW
    # is also not a process failure.
    return 1 if any(row["status"] == "ERROR" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
