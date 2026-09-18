#!/usr/bin/env python3
"""Replay serialized OpenMM systems/states to diagnose minimization divergence.

Expected layout (default root='.'):

    cpp/openmm_system.xml
    cpp/openmm_integrator.xml
    cpp/openmm_initial_state.xml
    cpp/md.screenout                 # optional, used only for validation

    python/openmm_system.xml
    python/openmm_integrator.xml
    python/openmm_initial_state.xml
    python/md.screenout              # optional, used only for validation

The script performs three related experiments for both systems:

1. Reconstruct the exact serialized System/Integrator/initial State and report
   the raw energy decomposition and constraint violations.
2. Apply Context.applyConstraints() only, then report the coordinate/energy
   change caused by projection onto each system's constraint manifold.
3. Reset to the raw initial coordinates and replay LocalEnergyMinimizer with
   the same tolerance/max-iteration semantics used by the MD drivers.

It then cross-evaluates both initial geometries and both replay-minimized
geometries under both Hamiltonians.  For each geometry row, the same box and
positions are used for both systems, so cross-system energy differences isolate
Hamiltonian differences at a fixed physical state.

Optionally, --also-minimize-projected performs a second minimization beginning
from the explicitly constraint-projected coordinates.  This is useful for
checking whether pre-projecting constraints changes the problematic result.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

try:
    import openmm as mm
    from openmm import unit
except ImportError as exc:  # pragma: no cover - runtime environment dependent
    raise SystemExit(
        "ERROR: OpenMM is required. Run this script from the same Python/OpenMM "
        "environment used for the validation calculations."
    ) from exc


ENERGY_GROUP_NAMES = {
    0: "Bond",
    1: "Angle",
    2: "Torsion",
    3: "Nonbonded",
    4: "Other",
}


@dataclass
class ReplayInput:
    label: str
    directory: Path
    system_xml: str
    integrator_xml: str
    initial_state_xml: str
    system: Any
    integrator: Any
    initial_state: Any
    context: Any
    raw_positions: Any
    raw_box: tuple
    constraint_tolerance: float
    max_iterations: int
    original_postmin_energy: float | None
    projected_positions: Any | None = None
    projected_box: tuple | None = None
    postmin_positions: Any | None = None
    postmin_box: tuple | None = None
    projected_postmin_positions: Any | None = None
    projected_postmin_box: tuple | None = None


class Reporter:
    def __init__(self) -> None:
        self.lines: list[str] = []

    def write(self, text: str = "") -> None:
        print(text, flush=True)
        self.lines.append(text)

    def save(self, path: Path) -> None:
        path.write_text("\n".join(self.lines) + "\n")


def read_text(path: Path) -> str:
    if not path.is_file():
        raise FileNotFoundError(path)
    return path.read_text()


def deserialize(xml_text: str):
    return mm.XmlSerializer.deserialize(xml_text)


def quantity_positions_nm(positions) -> np.ndarray:
    if hasattr(positions, "value_in_unit"):
        return np.asarray(positions.value_in_unit(unit.nanometer), dtype=float)
    return np.asarray(positions, dtype=float)


def state_positions(state):
    return state.getPositions(asNumpy=True)


def box_tuple(state) -> tuple:
    box = state.getPeriodicBoxVectors()
    return (box[0], box[1], box[2])


def set_context_geometry(context, positions, box: tuple) -> None:
    context.setPeriodicBoxVectors(box[0], box[1], box[2])
    context.setPositions(positions)


def copy_state_parameters(context, state) -> None:
    try:
        parameters = state.getParameters()
    except Exception:
        return
    for name, value in parameters.items():
        try:
            context.setParameter(name, value)
        except Exception:
            # A serialized initial State may contain a parameter that is not
            # present in the reconstructed System.  That should not prevent
            # replay of position-only/minimization diagnostics.
            pass


def configure_context_from_state(context, state) -> None:
    set_context_geometry(context, state_positions(state), box_tuple(state))
    copy_state_parameters(context, state)
    try:
        context.setTime(state.getTime())
    except Exception:
        pass
    try:
        context.setStepCount(state.getStepCount())
    except Exception:
        pass


def force_group_inventory(system) -> dict[int, list[str]]:
    inventory: dict[int, list[str]] = {}
    for index in range(system.getNumForces()):
        force = system.getForce(index)
        group = int(force.getForceGroup())
        inventory.setdefault(group, []).append(force.__class__.__name__)
    return inventory


def energy_components(context, system) -> dict[str, float]:
    values: dict[str, float] = {}
    groups = sorted(force_group_inventory(system))
    for group in groups:
        state = context.getState(getEnergy=True, groups=(1 << group))
        value = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        values[ENERGY_GROUP_NAMES.get(group, f"Group{group}")] = float(value)
    total = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )
    values["Total"] = float(total)
    return values


def vector_displacement(reference_positions, test_positions) -> dict[str, float]:
    ref = quantity_positions_nm(reference_positions)
    test = quantity_positions_nm(test_positions)
    if ref.shape != test.shape:
        raise ValueError(f"Position-shape mismatch: {ref.shape} vs {test.shape}")
    delta = test - ref
    per_atom = np.sqrt(np.sum(delta * delta, axis=1))
    return {
        "rms_nm": float(np.sqrt(np.mean(per_atom * per_atom))),
        "mean_nm": float(np.mean(per_atom)),
        "max_nm": float(np.max(per_atom)),
    }


def constraint_stats(system, positions) -> dict[str, Any]:
    xyz = quantity_positions_nm(positions)
    n_constraints = int(system.getNumConstraints())
    if n_constraints == 0:
        return {
            "count": 0,
            "rms_violation_nm": 0.0,
            "mean_abs_violation_nm": 0.0,
            "max_abs_violation_nm": 0.0,
            "count_gt_1e-5_nm": 0,
            "largest": [],
        }

    violations = np.empty(n_constraints, dtype=float)
    largest_rows: list[tuple[float, int, int, int, float, float]] = []
    for index in range(n_constraints):
        atom1, atom2, target = system.getConstraintParameters(index)
        atom1 = int(atom1)
        atom2 = int(atom2)
        target_nm = float(target.value_in_unit(unit.nanometer))
        actual_nm = float(np.linalg.norm(xyz[atom1] - xyz[atom2]))
        violation = actual_nm - target_nm
        violations[index] = violation
        largest_rows.append(
            (abs(violation), index, atom1, atom2, actual_nm, target_nm)
        )

    largest_rows.sort(reverse=True)
    largest = [
        {
            "constraint_index": int(index),
            "atom1": int(atom1),
            "atom2": int(atom2),
            "actual_nm": float(actual_nm),
            "target_nm": float(target_nm),
            "violation_nm": float(actual_nm - target_nm),
        }
        for _, index, atom1, atom2, actual_nm, target_nm in largest_rows[:10]
    ]

    return {
        "count": n_constraints,
        "rms_violation_nm": float(np.sqrt(np.mean(violations * violations))),
        "mean_abs_violation_nm": float(np.mean(np.abs(violations))),
        "max_abs_violation_nm": float(np.max(np.abs(violations))),
        "count_gt_1e-5_nm": int(np.count_nonzero(np.abs(violations) > 1.0e-5)),
        "largest": largest,
    }


def detect_max_iterations(screenout: Path, fallback: int) -> int:
    if not screenout.is_file():
        return fallback
    text = screenout.read_text(errors="replace")
    patterns = (
        r"maxIterations\s*=\s*(\d+)",
        r"Minimizing before equilibration for\s+(\d+)\s+steps",
        r"Minimizing(?:\s+with\s+OpenMM)?[^\n]*?(\d+)\s+(?:steps|maximum iterations)",
    )
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return int(match.group(1))
    return fallback


def parse_original_postmin_energy(screenout: Path) -> float | None:
    if not screenout.is_file():
        return None
    text = screenout.read_text(errors="replace")
    match = re.search(
        r"Post-minimization\s+Potential\s+Energy:\s*"
        r"([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)",
        text,
        flags=re.IGNORECASE,
    )
    if match:
        return float(match.group(1))
    return None


def create_context(system, integrator, platform_name: str, threads: int | None):
    if platform_name.lower() == "auto":
        return mm.Context(system, integrator)

    platform = mm.Platform.getPlatformByName(platform_name)
    properties: dict[str, str] = {}
    if threads is not None and platform.getName() == "CPU":
        properties["Threads"] = str(threads)
    if properties:
        return mm.Context(system, integrator, platform, properties)
    return mm.Context(system, integrator, platform)


def load_replay_input(
    root: Path,
    dirname: str,
    label: str,
    platform_name: str,
    threads: int | None,
    max_iterations_override: int | None,
    fallback_max_iterations: int,
) -> ReplayInput:
    directory = root / dirname
    system_path = directory / "openmm_system.xml"
    integrator_path = directory / "openmm_integrator.xml"
    state_path = directory / "openmm_initial_state.xml"
    screenout_path = directory / "md.screenout"

    system_xml = read_text(system_path)
    integrator_xml = read_text(integrator_path)
    state_xml = read_text(state_path)

    system = deserialize(system_xml)
    integrator = deserialize(integrator_xml)
    state = deserialize(state_xml)
    context = create_context(system, integrator, platform_name, threads)
    configure_context_from_state(context, state)

    try:
        constraint_tolerance = float(integrator.getConstraintTolerance())
    except Exception:
        constraint_tolerance = 1.0e-5

    max_iterations = (
        max_iterations_override
        if max_iterations_override is not None
        else detect_max_iterations(screenout_path, fallback_max_iterations)
    )

    return ReplayInput(
        label=label,
        directory=directory,
        system_xml=system_xml,
        integrator_xml=integrator_xml,
        initial_state_xml=state_xml,
        system=system,
        integrator=integrator,
        initial_state=state,
        context=context,
        raw_positions=state_positions(state),
        raw_box=box_tuple(state),
        constraint_tolerance=constraint_tolerance,
        max_iterations=max_iterations,
        original_postmin_energy=parse_original_postmin_energy(screenout_path),
    )


def print_energy_table(reporter: Reporter, title: str, values: dict[str, float]) -> None:
    reporter.write(title)
    for name, value in values.items():
        reporter.write(f"    {name:<12s} {value:18.6f} kJ/mol")


def print_constraint_summary(reporter: Reporter, title: str, stats: dict[str, Any]) -> None:
    reporter.write(title)
    reporter.write(
        "    "
        f"n={stats['count']}  rms={stats['rms_violation_nm']:.9g} nm  "
        f"mean|d|={stats['mean_abs_violation_nm']:.9g} nm  "
        f"max={stats['max_abs_violation_nm']:.9g} nm  "
        f">1e-5nm={stats['count_gt_1e-5_nm']}"
    )
    for row in stats["largest"][:3]:
        reporter.write(
            "      "
            f"constraint {row['constraint_index']}: {row['atom1']}-{row['atom2']}  "
            f"actual={row['actual_nm']:.9g} target={row['target_nm']:.9g}  "
            f"delta={row['violation_nm']:+.9g} nm"
        )


def energy_delta(after: dict[str, float], before: dict[str, float]) -> dict[str, float]:
    keys = sorted(set(after) | set(before))
    return {key: after.get(key, 0.0) - before.get(key, 0.0) for key in keys}


def replay_one(
    run: ReplayInput,
    reporter: Reporter,
    tolerance: float,
    platform_name: str,
    threads: int | None,
    also_minimize_projected: bool,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "label": run.label,
        "platform": run.context.getPlatform().getName(),
        "constraint_tolerance": run.constraint_tolerance,
        "minimization_tolerance_kj_per_mol_nm": tolerance,
        "max_iterations": run.max_iterations,
    }

    reporter.write()
    reporter.write("=" * 80)
    reporter.write(f"{run.label}: serialized replay")
    reporter.write("=" * 80)
    reporter.write(f"Platform             : {run.context.getPlatform().getName()}")
    reporter.write(f"Constraint tolerance : {run.constraint_tolerance:.9g}")
    reporter.write(f"Minimizer tolerance  : {tolerance:.9g} kJ/mol/nm")
    reporter.write(f"Max iterations       : {run.max_iterations}")

    group_inventory = force_group_inventory(run.system)
    reporter.write("Force groups:")
    for group, force_names in sorted(group_inventory.items()):
        reporter.write(f"    {group}: {', '.join(force_names)}")
    result["force_groups"] = {str(k): v for k, v in group_inventory.items()}

    raw_energy = energy_components(run.context, run.system)
    raw_constraints = constraint_stats(run.system, run.raw_positions)
    print_energy_table(reporter, "Raw serialized initial energy:", raw_energy)
    print_constraint_summary(reporter, "Raw initial constraint violations:", raw_constraints)

    result["raw_initial"] = {
        "energy_kj_per_mol": raw_energy,
        "constraints": raw_constraints,
    }

    run.context.applyConstraints(run.constraint_tolerance)
    projected_state = run.context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=False)
    run.projected_positions = state_positions(projected_state)
    run.projected_box = box_tuple(projected_state)
    projected_energy = energy_components(run.context, run.system)
    projected_constraints = constraint_stats(run.system, run.projected_positions)
    projection_displacement = vector_displacement(run.raw_positions, run.projected_positions)

    reporter.write()
    reporter.write("After Context.applyConstraints():")
    reporter.write(
        "    coordinate displacement: "
        f"rms={projection_displacement['rms_nm']:.9g} nm  "
        f"mean={projection_displacement['mean_nm']:.9g} nm  "
        f"max={projection_displacement['max_nm']:.9g} nm"
    )
    print_energy_table(reporter, "Projected energy:", projected_energy)
    print_energy_table(
        reporter,
        "Projection-only energy change (projected - raw):",
        energy_delta(projected_energy, raw_energy),
    )
    print_constraint_summary(reporter, "Projected constraint violations:", projected_constraints)

    result["projected"] = {
        "energy_kj_per_mol": projected_energy,
        "energy_change_kj_per_mol": energy_delta(projected_energy, raw_energy),
        "constraints": projected_constraints,
        "coordinate_displacement_nm": projection_displacement,
    }

    # Replay the actual MD-driver minimizer semantics from the original raw
    # serialized coordinates.  LocalEnergyMinimizer itself handles constraints
    # through temporary harmonic restraints; do not pre-project here.
    set_context_geometry(run.context, run.raw_positions, run.raw_box)
    mm.LocalEnergyMinimizer.minimize(run.context, tolerance, run.max_iterations)
    postmin_state = run.context.getState(
        getPositions=True,
        getEnergy=True,
        getParameters=True,
        enforcePeriodicBox=False,
    )
    run.postmin_positions = state_positions(postmin_state)
    run.postmin_box = box_tuple(postmin_state)
    postmin_energy = energy_components(run.context, run.system)
    postmin_constraints = constraint_stats(run.system, run.postmin_positions)
    postmin_displacement = vector_displacement(run.raw_positions, run.postmin_positions)

    reporter.write()
    reporter.write("Replay LocalEnergyMinimizer from raw initial positions:")
    reporter.write(
        "    coordinate displacement from raw: "
        f"rms={postmin_displacement['rms_nm']:.9g} nm  "
        f"mean={postmin_displacement['mean_nm']:.9g} nm  "
        f"max={postmin_displacement['max_nm']:.9g} nm"
    )
    print_energy_table(reporter, "Replay post-min energy:", postmin_energy)
    print_constraint_summary(reporter, "Replay post-min constraint violations:", postmin_constraints)

    if run.original_postmin_energy is not None:
        replay_total = postmin_energy["Total"]
        delta = replay_total - run.original_postmin_energy
        reporter.write(
            "Original md.screenout post-min PE: "
            f"{run.original_postmin_energy:.6f} kJ/mol"
        )
        reporter.write(
            f"Replay - original post-min PE  : {delta:+.6f} kJ/mol"
        )
        result["original_postmin_energy_kj_per_mol"] = run.original_postmin_energy
        result["replay_minus_original_postmin_kj_per_mol"] = delta

    replay_state_path = run.directory / "openmm_replay_postmin_state.xml"
    replay_state_path.write_text(mm.XmlSerializer.serialize(postmin_state))
    reporter.write(f"Wrote replay post-min State: {replay_state_path}")

    result["replay_postmin"] = {
        "energy_kj_per_mol": postmin_energy,
        "constraints": postmin_constraints,
        "coordinate_displacement_nm": postmin_displacement,
        "state_xml": str(replay_state_path),
    }

    if also_minimize_projected:
        # Reuse the existing Context rather than allocate another large PME
        # Context.  This matters for the 300k+ atom validation systems.
        set_context_geometry(run.context, run.raw_positions, run.raw_box)
        run.context.applyConstraints(run.constraint_tolerance)
        mm.LocalEnergyMinimizer.minimize(run.context, tolerance, run.max_iterations)
        projected_postmin_state = run.context.getState(
            getPositions=True,
            getEnergy=True,
            getParameters=True,
            enforcePeriodicBox=False,
        )
        run.projected_postmin_positions = state_positions(projected_postmin_state)
        run.projected_postmin_box = box_tuple(projected_postmin_state)
        projected_postmin_energy = energy_components(run.context, run.system)
        projected_postmin_constraints = constraint_stats(
            run.system, run.projected_postmin_positions
        )
        displacement_from_regular_postmin = vector_displacement(
            run.postmin_positions, run.projected_postmin_positions
        )

        reporter.write()
        reporter.write("Alternate replay: applyConstraints() BEFORE minimization:")
        print_energy_table(
            reporter,
            "Projected-start post-min energy:",
            projected_postmin_energy,
        )
        reporter.write(
            "    vs regular replay post-min positions: "
            f"rms={displacement_from_regular_postmin['rms_nm']:.9g} nm  "
            f"max={displacement_from_regular_postmin['max_nm']:.9g} nm"
        )
        print_constraint_summary(
            reporter,
            "Projected-start post-min constraint violations:",
            projected_postmin_constraints,
        )

        alt_state_path = run.directory / "openmm_replay_projected_postmin_state.xml"
        alt_state_path.write_text(mm.XmlSerializer.serialize(projected_postmin_state))
        reporter.write(f"Wrote projected-start post-min State: {alt_state_path}")

        result["projected_start_postmin"] = {
            "energy_kj_per_mol": projected_postmin_energy,
            "constraints": projected_postmin_constraints,
            "position_delta_vs_regular_postmin_nm": displacement_from_regular_postmin,
            "state_xml": str(alt_state_path),
        }

    return result


def cross_evaluate(
    geometry_name: str,
    geometries: dict[str, tuple[Any, tuple]],
    runs: dict[str, ReplayInput],
    reporter: Reporter,
) -> dict[str, Any]:
    reporter.write()
    reporter.write("=" * 80)
    reporter.write(f"Cross-Hamiltonian evaluation: {geometry_name}")
    reporter.write("=" * 80)
    reporter.write(
        "Each row uses one geometry source's positions AND box for both Hamiltonians."
    )

    results: dict[str, Any] = {}
    component_names: list[str] | None = None

    for geometry_label, (positions, box) in geometries.items():
        results[geometry_label] = {}
        reporter.write()
        reporter.write(f"Geometry source: {geometry_label}")
        for system_label, run in runs.items():
            set_context_geometry(run.context, positions, box)
            energies = energy_components(run.context, run.system)
            constraints = constraint_stats(run.system, positions)
            results[geometry_label][system_label] = {
                "energy_kj_per_mol": energies,
                "constraints": constraints,
            }
            if component_names is None:
                component_names = list(energies)
            reporter.write(
                f"  under {system_label:<8s}: "
                + "  ".join(
                    f"{name}={value:.6f}"
                    for name, value in energies.items()
                )
            )
            reporter.write(
                f"      constraint rms={constraints['rms_violation_nm']:.9g} nm  "
                f"max={constraints['max_abs_violation_nm']:.9g} nm"
            )

    # Add compact 2x2 matrices by energy component for easier visual diagnosis.
    if component_names is not None and len(geometries) == 2 and len(runs) == 2:
        geometry_labels = list(geometries)
        system_labels = list(runs)
        reporter.write()
        reporter.write("2x2 energy matrices (kJ/mol):")
        for component in component_names:
            reporter.write(f"  {component}:")
            reporter.write(
                f"    {'geometry':<12s} {system_labels[0]:>18s} {system_labels[1]:>18s}"
            )
            for geometry_label in geometry_labels:
                left = results[geometry_label][system_labels[0]]["energy_kj_per_mol"][component]
                right = results[geometry_label][system_labels[1]]["energy_kj_per_mol"][component]
                reporter.write(
                    f"    {geometry_label:<12s} {left:18.6f} {right:18.6f}"
                )

    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Replay serialized C++/Python OpenMM systems and diagnose minimization divergence."
    )
    parser.add_argument(
        "root",
        nargs="?",
        default=".",
        type=Path,
        help="Comparison root containing cpp/ and python/ (default: current directory)",
    )
    parser.add_argument(
        "--platform",
        default="CPU",
        help="OpenMM platform for replay, e.g. CPU, CUDA, OpenCL, or auto (default: CPU)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="CPU platform thread count (default: OpenMM/platform default)",
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=None,
        help=(
            "Override minimizer max iterations for both systems. If omitted, "
            "detect from each md.screenout, falling back to 500."
        ),
    )
    parser.add_argument(
        "--fallback-max-iterations",
        type=int,
        default=500,
        help="Fallback when max iterations cannot be read from md.screenout (default: 500)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=10.0,
        help="LocalEnergyMinimizer RMS-force tolerance in kJ/mol/nm (default: 10)",
    )
    parser.add_argument(
        "--also-minimize-projected",
        action="store_true",
        help=(
            "Also applyConstraints() before a second minimization for each system. "
            "This roughly doubles replay minimization cost."
        ),
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=None,
        help="Human-readable report path (default: ROOT/openmm_replay_diagnostic.txt)",
    )
    parser.add_argument(
        "--json",
        dest="json_path",
        type=Path,
        default=None,
        help="JSON result path (default: ROOT/openmm_replay_diagnostic.json)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    report_path = args.report or (root / "openmm_replay_diagnostic.txt")
    json_path = args.json_path or (root / "openmm_replay_diagnostic.json")

    if args.max_iterations is not None and args.max_iterations <= 0:
        raise SystemExit("ERROR: --max-iterations must be greater than zero")
    if args.fallback_max_iterations <= 0:
        raise SystemExit("ERROR: --fallback-max-iterations must be greater than zero")
    if args.tolerance <= 0:
        raise SystemExit("ERROR: --tolerance must be greater than zero")
    if args.threads is not None and args.threads <= 0:
        raise SystemExit("ERROR: --threads must be greater than zero")

    reporter = Reporter()
    reporter.write("OpenMM serialized replay diagnostic")
    reporter.write("=" * 80)
    reporter.write(f"Root                  : {root}")
    reporter.write(f"Replay platform       : {args.platform}")
    reporter.write(f"Minimizer tolerance   : {args.tolerance} kJ/mol/nm")
    reporter.write(
        "Projected minimization: "
        + ("enabled" if args.also_minimize_projected else "disabled")
    )

    try:
        cpp = load_replay_input(
            root,
            "cpp",
            "C++",
            args.platform,
            args.threads,
            args.max_iterations,
            args.fallback_max_iterations,
        )
        python_run = load_replay_input(
            root,
            "python",
            "Python",
            args.platform,
            args.threads,
            args.max_iterations,
            args.fallback_max_iterations,
        )
    except FileNotFoundError as exc:
        raise SystemExit(f"ERROR: required serialized file not found: {exc}") from exc

    runs = {"C++": cpp, "Python": python_run}
    output: dict[str, Any] = {
        "root": str(root),
        "platform_requested": args.platform,
        "minimization_tolerance_kj_per_mol_nm": args.tolerance,
        "also_minimize_projected": args.also_minimize_projected,
        "runs": {},
    }

    # Cross-evaluate the raw initial states before either context is minimized.
    initial_geometries = {
        "C++": (cpp.raw_positions, cpp.raw_box),
        "Python": (python_run.raw_positions, python_run.raw_box),
    }
    output["cross_initial"] = cross_evaluate(
        "raw initial geometries", initial_geometries, runs, reporter
    )

    # Restore each Context to its own initial state, then replay independently.
    configure_context_from_state(cpp.context, cpp.initial_state)
    configure_context_from_state(python_run.context, python_run.initial_state)
    output["runs"]["C++"] = replay_one(
        cpp,
        reporter,
        args.tolerance,
        args.platform,
        args.threads,
        args.also_minimize_projected,
    )
    output["runs"]["Python"] = replay_one(
        python_run,
        reporter,
        args.tolerance,
        args.platform,
        args.threads,
        args.also_minimize_projected,
    )

    postmin_geometries = {
        "C++": (cpp.postmin_positions, cpp.postmin_box),
        "Python": (python_run.postmin_positions, python_run.postmin_box),
    }
    output["cross_postmin"] = cross_evaluate(
        "replay-minimized geometries", postmin_geometries, runs, reporter
    )

    if args.also_minimize_projected:
        projected_postmin_geometries = {
            "C++": (cpp.projected_postmin_positions, cpp.projected_postmin_box),
            "Python": (
                python_run.projected_postmin_positions,
                python_run.projected_postmin_box,
            ),
        }
        output["cross_projected_start_postmin"] = cross_evaluate(
            "constraint-projected-then-minimized geometries",
            projected_postmin_geometries,
            runs,
            reporter,
        )

    reporter.write()
    reporter.write("=" * 80)
    reporter.write("Interpretation guide")
    reporter.write("=" * 80)
    reporter.write(
        "* Large applyConstraints() energy shifts indicate that projection onto the "
        "constraint manifold materially changes the ordinary Hamiltonian energy."
    )
    reporter.write(
        "* If replay post-min energies reproduce md.screenout, the serialized "
        "System/State are sufficient to reproduce the original minimization behavior."
    )
    reporter.write(
        "* In a cross matrix, a large C++-vs-Python difference within the SAME geometry "
        "row is a Hamiltonian difference."
    )
    reporter.write(
        "* If each geometry is low only under the Hamiltonian that produced it, the two "
        "Systems are driving minimization toward different basins."
    )
    reporter.write(
        "* If --also-minimize-projected substantially changes one result, explicit "
        "constraint projection before minimization is a concrete protocol variable to test."
    )

    reporter.save(report_path)
    json_path.write_text(json.dumps(output, indent=2))
    reporter.write(f"Wrote report: {report_path}")
    reporter.write(f"Wrote JSON  : {json_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
