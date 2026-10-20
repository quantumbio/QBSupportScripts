#!/usr/bin/env python3
"""
Triage DivCon/C++ vs Python/OpenMM validation JSON summaries.

This script is intentionally a review triage tool, not a scientific pass/fail
classifier.  It reads the *_summary.json files written by analysisMD.py and
extracts the strongest engineering-validation metrics into one table.

Typical use from the directory containing the four-character PDBID directories:

    python summarize_MD_validation_v1.py --root . --list anal.list

or discover all immediate four-character alphanumeric directories automatically:

    python summarize_MD_validation_v1.py --root .

For each PDBID the expected JSON path is:

    <root>/<PDBID>/<PDBID>_summary.json

Outputs (prefix defaults to "md_validation"):

    md_validation_summary.csv
    md_validation_report.txt
    md_validation_review.list
    md_validation_good.list
    md_validation_good_limited.list
    md_validation_errors.list

The default thresholds are conservative smoke-test triage thresholds chosen for
short C++/Python validation runs.  They are deliberately kept in one dictionary
near the top of the script so they can be adjusted after more empirical data are
available.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
import sys
from pathlib import Path
from typing import Any


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
    # only because the two workflows can enter OpenMM from slightly different
    # starting geometries.
    "postmin_total_sym_pct": 0.50,
    "postmin_nonbonded_sym_pct": 0.50,

    # Short stochastic production comparison.  These are triage limits, not
    # ensemble-equivalence criteria.
    "production_pe_mean_sym_pct": 0.75,
    "production_pe_max_sym_pct": 1.00,
    "production_temperature_mean_sym_pct": 5.0,
    "production_volume_mean_sym_pct": 5.0,
    "production_density_mean_sym_pct": 5.0,

    # Structural comparison.  Internal C-alpha distances are especially useful
    # because they are invariant to overall translation/rotation.
    "ca_distance_rmse_angstrom": 1.0,
    "ca_distance_pearson_min": 0.98,
    "rg_relative_difference_pct": 5.0,
    "rmsf_profile_rmse_angstrom": 0.50,

    # Secondary short-trajectory indicators.  These are reported but intentionally
    # not used as primary review triggers by default.
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
    "ca_distance_rmse_angstrom",
    "rg_relative_difference_pct",
)


CSV_FIELDS = [
    "pdbid",
    "status",
    "json_path",
    "label1",
    "label2",
    "frames_input1",
    "frames_input2",
    "atoms_input1",
    "atoms_input2",
    "charge_input1_e",
    "charge_input2_e",
    "charge_delta_e",
    "initial_box_max_delta_nm",
    "postmin_total_sym_pct",
    "postmin_nonbonded_sym_pct",
    "postmin_bond_sym_pct",
    "postmin_angle_sym_pct",
    "postmin_torsion_sym_pct",
    "production_matched_steps",
    "production_pe_mean_sym_pct",
    "production_pe_max_sym_pct",
    "production_temperature_mean_sym_pct",
    "production_volume_mean_sym_pct",
    "production_density_mean_sym_pct",
    "rmsf_profile_rmse_angstrom",
    "rg_relative_difference_pct",
    "common_rmsd_input1_mean_angstrom",
    "common_rmsd_input2_mean_angstrom",
    "common_rmsd_mean_delta_angstrom",
    "ca_distance_rmse_angstrom",
    "ca_distance_pearson_r",
    "ca_distance_max_abs_angstrom",
    "hbond_input1",
    "hbond_input2",
    "hbond_shared",
    "hbond_jaccard",
    "dssp_diff_fraction",
    "system_mismatch_count",
    "protocol_mismatch_count",
    "missing_core_diagnostics",
    "review_reasons",
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

    # Numeric Hamiltonian settings should be effectively identical when both are
    # available.  Use a tight absolute tolerance rather than a relative percent.
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
        add_reason(
            reasons,
            f"water constraint delta={constraint_delta:.3g} nm",
        )


def extract_metrics(pdbid: str, json_path: Path, summary: dict) -> dict[str, Any]:
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
            "pdbid": pdbid,
            "json_path": str(json_path),
            "label1": label1,
            "label2": label2,
        }
    )

    # Basic trajectory identity.
    row["frames_input1"] = nested(summary, "n_frames", label1)
    row["frames_input2"] = nested(summary, "n_frames", label2)
    row["atoms_input1"] = nested(summary, "n_atoms", label1)
    row["atoms_input2"] = nested(summary, "n_atoms", label2)

    # System/Hamiltonian construction.
    row["charge_input1_e"] = finite_float(
        nested(input1, "system", "total_particle_charge_e")
    )
    row["charge_input2_e"] = finite_float(
        nested(input2, "system", "total_particle_charge_e")
    )
    if row["charge_input1_e"] is not None and row["charge_input2_e"] is not None:
        row["charge_delta_e"] = row["charge_input1_e"] - row["charge_input2_e"]

    row["initial_box_max_delta_nm"] = finite_float(
        nested(comparison, "system", "initial_box_nm", "max_absolute_difference")
    )

    # Post-minimization energy decomposition.  Do not classify on the initial
    # decomposition because starting geometries can differ before minimization.
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

    # Structural metrics.  Deliberately use numeric metrics only; residue-specific
    # top-difference labels are not needed for triage.
    row["rmsf_profile_rmse_angstrom"] = finite_float(
        nested(trajectory, "rmsf_profile", "rmse")
    )
    row["rg_relative_difference_pct"] = finite_float(
        nested(trajectory, "rg_stability", "relative_mean_difference_percent")
    )

    common = trajectory.get("common_reference_rmsd") or {}
    common1 = finite_float(nested(common, label1, "mean"))
    common2 = finite_float(nested(common, label2, "mean"))
    row["common_rmsd_input1_mean_angstrom"] = common1
    row["common_rmsd_input2_mean_angstrom"] = common2
    if common1 is not None and common2 is not None:
        row["common_rmsd_mean_delta_angstrom"] = common1 - common2

    ca = trajectory.get("ca_internal_distances") or {}
    row["ca_distance_rmse_angstrom"] = finite_float(ca.get("rmse"))
    row["ca_distance_pearson_r"] = finite_float(ca.get("pearson_r"))
    row["ca_distance_max_abs_angstrom"] = finite_float(ca.get("max_abs_difference"))

    hbonds = summary.get("hbonds") or {}
    row["hbond_input1"] = hbonds.get("n1")
    row["hbond_input2"] = hbonds.get("n2")
    row["hbond_shared"] = hbonds.get("shared")
    try:
        n1 = int(row["hbond_input1"])
        n2 = int(row["hbond_input2"])
        shared = int(row["hbond_shared"])
        union = n1 + n2 - shared
        row["hbond_jaccard"] = float(shared / union) if union > 0 else 1.0
    except (TypeError, ValueError):
        row["hbond_jaccard"] = None

    row["dssp_diff_fraction"] = finite_float(summary.get("dssp_diff_pct"))

    return row


def classify(summary: dict, row: dict[str, Any]) -> None:
    reasons: list[str] = []
    notes: list[str] = []

    frames1 = row.get("frames_input1")
    frames2 = row.get("frames_input2")
    atoms1 = row.get("atoms_input1")
    atoms2 = row.get("atoms_input2")

    if frames1 is not None and frames2 is not None and frames1 != frames2:
        add_reason(reasons, f"trajectory frame mismatch ({frames1} vs {frames2})")
    if atoms1 is not None and atoms2 is not None and atoms1 != atoms2:
        add_reason(reasons, f"trajectory atom mismatch ({atoms1} vs {atoms2})")

    row["system_mismatch_count"] = compare_exact_system_fields(summary, reasons)
    row["protocol_mismatch_count"] = compare_protocol(summary, reasons)
    representative_parameter_checks(summary, reasons)

    q1 = finite_float(row.get("charge_input1_e"))
    q2 = finite_float(row.get("charge_input2_e"))
    if q1 is not None and abs(q1) > THRESHOLDS["charge_abs_e"]:
        add_reason(reasons, f"{row['label1']} charge={q1:+.6g} e")
    if q2 is not None and abs(q2) > THRESHOLDS["charge_abs_e"]:
        add_reason(reasons, f"{row['label2']} charge={q2:+.6g} e")
    if q1 is not None and q2 is not None:
        qdelta = abs(q1 - q2)
        if qdelta > THRESHOLDS["charge_delta_e"]:
            add_reason(reasons, f"charge disagreement={qdelta:.6g} e")

    box_delta = finite_float(row.get("initial_box_max_delta_nm"))
    if (
        box_delta is not None
        and box_delta > THRESHOLDS["initial_box_max_delta_nm"]
    ):
        add_reason(reasons, f"initial box max delta={box_delta:.4g} nm")

    post_total = finite_float(row.get("postmin_total_sym_pct"))
    if post_total is not None and post_total > THRESHOLDS["postmin_total_sym_pct"]:
        add_reason(reasons, f"post-min total energy={post_total:.3f}%")

    post_nb = finite_float(row.get("postmin_nonbonded_sym_pct"))
    if post_nb is not None and post_nb > THRESHOLDS["postmin_nonbonded_sym_pct"]:
        add_reason(reasons, f"post-min nonbonded energy={post_nb:.3f}%")

    pe_mean = finite_float(row.get("production_pe_mean_sym_pct"))
    if pe_mean is not None and pe_mean > THRESHOLDS["production_pe_mean_sym_pct"]:
        add_reason(reasons, f"production PE mean={pe_mean:.3f}%")

    pe_max = finite_float(row.get("production_pe_max_sym_pct"))
    if pe_max is not None and pe_max > THRESHOLDS["production_pe_max_sym_pct"]:
        add_reason(reasons, f"production PE max={pe_max:.3f}%")

    for column, threshold_key, label in (
        (
            "production_temperature_mean_sym_pct",
            "production_temperature_mean_sym_pct",
            "production temperature mean",
        ),
        (
            "production_volume_mean_sym_pct",
            "production_volume_mean_sym_pct",
            "production volume mean",
        ),
        (
            "production_density_mean_sym_pct",
            "production_density_mean_sym_pct",
            "production density mean",
        ),
    ):
        value = finite_float(row.get(column))
        if value is not None and value > THRESHOLDS[threshold_key]:
            add_reason(reasons, f"{label}={value:.3f}%")

    ca_rmse = finite_float(row.get("ca_distance_rmse_angstrom"))
    if ca_rmse is not None and ca_rmse > THRESHOLDS["ca_distance_rmse_angstrom"]:
        qualifier = "SEVERE " if ca_rmse > 3.0 else ""
        add_reason(reasons, f"{qualifier}C-alpha distance RMSE={ca_rmse:.3f} A")

    ca_r = finite_float(row.get("ca_distance_pearson_r"))
    if ca_r is not None and ca_r < THRESHOLDS["ca_distance_pearson_min"]:
        add_reason(reasons, f"C-alpha distance Pearson r={ca_r:.4f}")

    rg_relative = finite_float(row.get("rg_relative_difference_pct"))
    if (
        rg_relative is not None
        and abs(rg_relative) > THRESHOLDS["rg_relative_difference_pct"]
    ):
        qualifier = "SEVERE " if abs(rg_relative) > 20.0 else ""
        add_reason(reasons, f"{qualifier}Rg mean delta={rg_relative:+.3f}%")

    rmsf_rmse = finite_float(row.get("rmsf_profile_rmse_angstrom"))
    if (
        rmsf_rmse is not None
        and rmsf_rmse > THRESHOLDS["rmsf_profile_rmse_angstrom"]
    ):
        add_reason(reasons, f"RMSF profile RMSE={rmsf_rmse:.3f} A")

    # Secondary metrics are shown as notes only.  Five-frame smoke tests do not
    # justify turning H-bond/DSSP differences into hard validation criteria.
    hbond_jaccard = finite_float(row.get("hbond_jaccard"))
    if (
        hbond_jaccard is not None
        and hbond_jaccard < THRESHOLDS["hbond_jaccard_report"]
    ):
        notes.append(f"low H-bond Jaccard={hbond_jaccard:.3f}")

    dssp = finite_float(row.get("dssp_diff_fraction"))
    if dssp is not None and dssp > THRESHOLDS["dssp_difference_fraction_report"]:
        notes.append(f"DSSP disagreement={100.0 * dssp:.1f}%")

    # Track missing critical diagnostics separately rather than converting every
    # absent field into a false scientific failure.  In particular, some C++ logs
    # do not currently print all system-construction diagnostics.
    missing = [name for name in CORE_DIAGNOSTICS if row.get(name) is None]
    row["missing_core_diagnostics"] = ";".join(missing)

    if reasons:
        row["status"] = "REVIEW"
    elif missing:
        row["status"] = "GOOD_LIMITED"
    else:
        row["status"] = "GOOD"

    row["review_reasons"] = "; ".join(reasons)
    row["notes"] = "; ".join(notes)


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


def error_row(pdbid: str, json_path: Path, message: str) -> dict[str, Any]:
    row = {field: None for field in CSV_FIELDS}
    row.update(
        {
            "pdbid": pdbid,
            "status": "ERROR",
            "json_path": str(json_path),
            "review_reasons": message,
            "notes": "",
            "missing_core_diagnostics": ";".join(CORE_DIAGNOSTICS),
        }
    )
    return row


def analyze_one(pdbid: str, json_path: Path) -> dict[str, Any]:
    if not json_path.is_file():
        return error_row(pdbid, json_path, "summary JSON missing")

    try:
        with json_path.open() as handle:
            summary = json.load(handle)
    except Exception as exc:
        return error_row(pdbid, json_path, f"cannot parse JSON: {exc}")

    if not isinstance(summary, dict):
        return error_row(pdbid, json_path, "JSON root is not an object")

    try:
        row = extract_metrics(pdbid, json_path, summary)
        classify(summary, row)
        return row
    except Exception as exc:
        return error_row(pdbid, json_path, f"analysis error: {exc}")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field) for field in CSV_FIELDS})


def write_list(path: Path, rows: list[dict[str, Any]], statuses: set[str]) -> None:
    with path.open("w") as handle:
        for row in rows:
            if row.get("status") in statuses:
                handle.write(f"{row['pdbid']}\n")


def numeric_distribution(rows: list[dict[str, Any]], field: str) -> dict[str, float] | None:
    values = [finite_float(row.get(field)) for row in rows if row.get("status") != "ERROR"]
    clean = [value for value in values if value is not None]
    if not clean:
        return None

    result = {
        "n": float(len(clean)),
        "min": min(clean),
        "median": statistics.median(clean),
        "max": max(clean),
    }
    if len(clean) >= 2:
        result["mean"] = statistics.fmean(clean)
    return result


def write_report(path: Path, rows: list[dict[str, Any]], root: Path) -> None:
    counts = {status: 0 for status in ("REVIEW", "GOOD", "GOOD_LIMITED", "ERROR")}
    for row in rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1

    with path.open("w") as handle:
        handle.write("DivCon/OpenMM validation triage report\n")
        handle.write("====================================\n\n")
        handle.write(f"Root: {root}\n")
        handle.write(f"Structures considered: {len(rows)}\n")
        handle.write(f"REVIEW       : {counts.get('REVIEW', 0)}\n")
        handle.write(f"GOOD         : {counts.get('GOOD', 0)}\n")
        handle.write(f"GOOD_LIMITED : {counts.get('GOOD_LIMITED', 0)}\n")
        handle.write(f"ERROR        : {counts.get('ERROR', 0)}\n\n")

        handle.write("Interpretation\n")
        handle.write("--------------\n")
        handle.write("REVIEW: one or more conservative triage thresholds were exceeded.\n")
        handle.write("GOOD: no review threshold was exceeded and all core diagnostics were present.\n")
        handle.write("GOOD_LIMITED: no available metric triggered review, but one or more core diagnostics were absent.\n")
        handle.write("ERROR: summary JSON was missing, unreadable, or could not be processed.\n")
        handle.write("These categories are engineering triage labels, not proof of ensemble equivalence.\n\n")

        handle.write("Primary thresholds\n")
        handle.write("------------------\n")
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
                handle.write(f"{row['pdbid']}: {row['review_reasons']}\n")
        handle.write("\n")

        handle.write("Errors / missing summaries\n")
        handle.write("--------------------------\n")
        error_rows = [row for row in rows if row["status"] == "ERROR"]
        if not error_rows:
            handle.write("None\n")
        else:
            for row in error_rows:
                handle.write(f"{row['pdbid']}: {row['review_reasons']}\n")
        handle.write("\n")

        handle.write("Dataset metric ranges\n")
        handle.write("---------------------\n")
        distribution_fields = (
            "postmin_total_sym_pct",
            "postmin_nonbonded_sym_pct",
            "production_pe_mean_sym_pct",
            "production_pe_max_sym_pct",
            "ca_distance_rmse_angstrom",
            "ca_distance_pearson_r",
            "rg_relative_difference_pct",
            "rmsf_profile_rmse_angstrom",
        )
        for field in distribution_fields:
            stats = numeric_distribution(rows, field)
            if not stats:
                handle.write(f"{field}: no data\n")
                continue
            handle.write(
                f"{field}: n={int(stats['n'])} min={stats['min']:.6g} "
                f"median={stats['median']:.6g} max={stats['max']:.6g}"
            )
            if "mean" in stats:
                handle.write(f" mean={stats['mean']:.6g}")
            handle.write("\n")
        handle.write("\n")

        handle.write("Compact table\n")
        handle.write("-------------\n")
        header = (
            f"{'PDB':<6} {'STATUS':<13} {'postTot%':>9} {'postNB%':>9} "
            f"{'prodPE%':>9} {'CA_RMSE':>9} {'CA_r':>7} {'Rg%':>8} {'q1':>10} {'q2':>10}\n"
        )
        handle.write(header)
        handle.write("-" * (len(header.rstrip()) + 2) + "\n")
        for row in rows:
            handle.write(
                f"{row['pdbid']:<6} {row['status']:<13} "
                f"{fmt(row.get('postmin_total_sym_pct')):>9} "
                f"{fmt(row.get('postmin_nonbonded_sym_pct')):>9} "
                f"{fmt(row.get('production_pe_mean_sym_pct')):>9} "
                f"{fmt(row.get('ca_distance_rmse_angstrom')):>9} "
                f"{fmt(row.get('ca_distance_pearson_r')):>7} "
                f"{fmt(row.get('rg_relative_difference_pct')):>8} "
                f"{fmt(row.get('charge_input1_e')):>10} "
                f"{fmt(row.get('charge_input2_e')):>10}\n"
            )


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
        f"ERROR={counts.get('ERROR', 0)}"
    )
    print()

    review_rows = [row for row in rows if row["status"] == "REVIEW"]
    if review_rows:
        print("Manual review candidates:")
        for row in review_rows:
            print(f"  {row['pdbid']}: {row['review_reasons']}")
        print()

    print("Wrote:")
    for label, path in outputs.items():
        print(f"  {label:<13}: {path}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize <PDBID>/<PDBID>_summary.json validation results and "
            "flag conservative manual-review candidates."
        )
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
        help="Optional file containing PDB IDs to process; blank/comment lines are ignored.",
    )
    parser.add_argument(
        "--prefix",
        default="md_validation",
        help="Prefix for generated report/list/CSV files (default: md_validation).",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    if not root.is_dir():
        print(f"ERROR: root directory does not exist: {root}", file=sys.stderr)
        return 2

    if args.list_file:
        if not args.list_file.is_file():
            print(f"ERROR: list file does not exist: {args.list_file}", file=sys.stderr)
            return 2
        pdbids = read_pdbids(args.list_file)
    else:
        pdbids = discover_pdbids(root)

    if not pdbids:
        print("ERROR: no PDB IDs found to process.", file=sys.stderr)
        return 1

    rows = [
        analyze_one(pdbid, root / pdbid / f"{pdbid}_summary.json")
        for pdbid in pdbids
    ]

    # Put manual-review candidates first, then limited/good results, then errors.
    rank = {"REVIEW": 0, "GOOD_LIMITED": 1, "GOOD": 2, "ERROR": 3}
    rows.sort(key=lambda row: (rank.get(str(row.get("status")), 99), row["pdbid"]))

    prefix = Path(args.prefix)
    csv_path = prefix.with_name(prefix.name + "_summary.csv")
    report_path = prefix.with_name(prefix.name + "_report.txt")
    review_path = prefix.with_name(prefix.name + "_review.list")
    good_path = prefix.with_name(prefix.name + "_good.list")
    limited_path = prefix.with_name(prefix.name + "_good_limited.list")
    error_path = prefix.with_name(prefix.name + "_errors.list")

    write_csv(csv_path, rows)
    write_report(report_path, rows, root)
    write_list(review_path, rows, {"REVIEW"})
    write_list(good_path, rows, {"GOOD"})
    write_list(limited_path, rows, {"GOOD_LIMITED"})
    write_list(error_path, rows, {"ERROR"})

    outputs = {
        "CSV": csv_path,
        "report": report_path,
        "review list": review_path,
        "good list": good_path,
        "limited list": limited_path,
        "error list": error_path,
    }
    print_console_summary(rows, outputs)

    # Missing/invalid JSON is a scripting/data-flow failure and should produce a
    # nonzero status.  REVIEW is intentionally not a process failure.
    return 1 if any(row["status"] == "ERROR" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
