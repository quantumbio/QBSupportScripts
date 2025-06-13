#!/usr/bin/env python3
"""
Compare two MD trajectories:

* Standardise residue names (handles Amber protonation codes).
* Align sequences (one‑letter codes) and superpose on Cα atoms.
* Compute per‑residue RMSF, radius of gyration, DSSP.
* Perform Kolmogorov–Smirnov tests and save comparison plots.

The script accepts plain or .gz‑compressed files; compressed inputs
are transparently decompressed to temporary files before loading.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path
from typing import List

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
from Bio.Align import PairwiseAligner
from Bio.Data import IUPACData

###############################################################################
# ---------------------------  CONFIGURATION  ---------------------------------
###############################################################################

# Map non‑standard / protonation‑state names to canonical 3‑letter codes.
residue_mapping = {
    # standard AAs
    "ALA": "ALA", "ARG": "ARG", "ASN": "ASN", "ASP": "ASP", "CYS": "CYS",
    "GLU": "GLU", "GLN": "GLN", "GLY": "GLY", "HIS": "HIS", "ILE": "ILE",
    "LEU": "LEU", "LYS": "LYS", "MET": "MET", "PHE": "PHE", "PRO": "PRO",
    "SER": "SER", "THR": "THR", "TRP": "TRP", "TYR": "TYR", "VAL": "VAL",

    # Amber protonation / variants
    "HIE": "HIS", "HID": "HIS", "HIP": "HIS",  # Histidine
    "CYX": "CYS",                              # Disulfide Cys
    "GLH": "GLU",                              # Protonated Glu
    "ASH": "ASP",                              # Protonated Asp
}

aa3to1 = {k.upper(): v for k, v in IUPACData.protein_letters_3to1_extended.items()}

###############################################################################
# ---------------------------  HELPER FUNCTIONS  ------------------------------
###############################################################################


def decompress_if_gz(path: Path) -> Path:
    """
    Return a path to an uncompressed copy of *path*.

    If *path* ends with '.gz', its contents are decompressed into a named
    temporary file (kept until program exit).  Otherwise the original path is
    returned unchanged.
    """
    if not path.suffix == ".gz":
        return path

    # Determine original (undecompressed) extension (e.g. '.xtc', '.pdb')
    orig_suffix = "".join(path.stem.split(".")[1:]) or path.stem.split(".")[0]
    tmp = tempfile.NamedTemporaryFile(
        delete=False, suffix=f".{orig_suffix}", prefix="mdtraj_"
    )
    logging.info("Decompressing %s → %s", path, tmp.name)
    with gzip.open(path, "rb") as gz_in, open(tmp.name, "wb") as out:
        shutil.copyfileobj(gz_in, out)
    return Path(tmp.name)


def trim_to_three_letters(name: str) -> str:
    """
    Keep the *first* three consecutive uppercase letters of *name*.

    Examples
    --------
    'LYS1' → 'LYS'
    'CYX'  → 'CYX'
    'GLH'  → 'GLH'
    'SEP'  → 'SEP'
    """
    m = re.match(r"([A-Z]{3})", name.upper())
    return m.group(1) if m else name[:3].upper()


def standardise_residue_names(traj: md.Trajectory) -> List[str]:
    """Return list of canonical 3‑letter codes corresponding to *traj* residues."""
    std: List[str] = []
    for res in traj.topology.residues:
        core = trim_to_three_letters(res.name)
        mapped = residue_mapping.get(core, core)
        if mapped not in residue_mapping.values() and mapped not in {"HOH", "NA", "CL"}:
            logging.warning("Unmatched residue name '%s'", res.name)
        std.append(mapped)

    return std


def to_one_letter(seq3: List[str]) -> str:
    """Translate 3‑letter codes → one‑letter codes, dropping unknowns."""
    letters: List[str] = []
    for aa in seq3:
        letter = aa3to1.get(aa, None)
        if letter:
            letters.append(letter)
    return "".join(letters)


def load_trajectory(traj_path: Path, top_path: Path) -> md.Trajectory:
    """Load a trajectory with MDTraj, printing basic info."""
    traj = md.load(traj_path, top=top_path)
    logging.info("Loaded %s | frames: %d | atoms: %d",
                 traj_path.name, traj.n_frames, traj.n_atoms)
    return traj


def compute_rmsf(traj: md.Trajectory, ca_idx: np.ndarray) -> np.ndarray:
    """
    Compute RMSF per Cα atom (Å) for *traj*.

    MDTraj returns nm; convert to Å for easier human digestion.
    """
    rmsf = md.rmsf(traj, reference=traj, atom_indices=ca_idx)  # nm
    return rmsf * 10.0  # → Å


###############################################################################
# -------------------------------  MAIN  --------------------------------------
###############################################################################


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare two MD trajectories (RMSF, Rg, DSSP).")
    parser.add_argument("traj1", type=Path, help="Trajectory 1 file  (DCD/XTC/… or .gz)")
    parser.add_argument("top1",  type=Path, help="Topology    1 file  (PDB/PRMTOP/… or .gz)")
    parser.add_argument("traj2", type=Path, help="Trajectory 2 file  (DCD/XTC/… or .gz)")
    parser.add_argument("top2",  type=Path, help="Topology    2 file  (PDB/PRMTOP/… or .gz)")
    parser.add_argument("-o", "--out-prefix", default="comparison",
                        help="Prefix for output PDFs (default: %(default)s)")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Reduce log output.")
    parser.add_argument("--label1", default="Complete Structure",
                        help="Label for trajectory 1 (default: %(default)s)")
    parser.add_argument("--label2", default="Published Structure",
                        help="Label for trajectory 2 (default: %(default)s)")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="%(levelname)s: %(message)s")

    # Track any temporary files we create so we can delete them later.
    temp_paths: List[Path] = []
    try:
        traj1_path = decompress_if_gz(args.traj1)
        top1_path  = decompress_if_gz(args.top1)
        traj2_path = decompress_if_gz(args.traj2)
        top2_path  = decompress_if_gz(args.top2)
        temp_paths.extend(p for p in (traj1_path, top1_path, traj2_path, top2_path)
                          if p != args.traj1 and p != args.top1 and p != args.traj2 and p != args.top2)

        # Load trajectories
        t1_raw = load_trajectory(traj1_path, top1_path)
        t2_raw = load_trajectory(traj2_path, top2_path)

        # ---------------------------------------------------------------------
        # Standardise residue names and update topology (in‑place for MDTraj).
        # ---------------------------------------------------------------------
        std1 = standardise_residue_names(t1_raw)
        std2 = standardise_residue_names(t2_raw)
        for res, new in zip(t1_raw.topology.residues, std1):
            res.name = new
        for res, new in zip(t2_raw.topology.residues, std2):
            res.name = new

        # ---------------------------------------------------------------------
        # Sequence alignment (one‑letter codes)
        # ---------------------------------------------------------------------
        seq1 = to_one_letter(std1)
        seq2 = to_one_letter(std2)
        aligner = PairwiseAligner()
        aln = aligner.align(seq1, seq2)[0]
        logging.info("Global alignment score: %.2f  length: %d",
                     aln.score, aln.shape[1])

        # -----------------------------------------------------------------
        # Keep *protein only* (drops solvent ions AND ligands/ cofactors)
        # -----------------------------------------------------------------
        prot_atoms_1 = [a.index for a in t1_raw.topology.atoms if a.residue.is_protein]
        prot_atoms_2 = [a.index for a in t2_raw.topology.atoms if a.residue.is_protein]
        
        t1 = t1_raw.atom_slice(prot_atoms_1)
        t2 = t2_raw.atom_slice(prot_atoms_2)
#        prot_sel = "not (resname HOH CL NA)"
#        t1 = t1_raw.atom_slice(t1_raw.topology.select(prot_sel))
#        t2 = t2_raw.atom_slice(t2_raw.topology.select(prot_sel))

        # ---------------------------------------------------------------------
        # Superpose on Cα atoms
        # ---------------------------------------------------------------------
        ca1 = t1.topology.select("protein and name CA")
        ca2 = t2.topology.select("protein and name CA")
        common_len = min(len(ca1), len(ca2))
        if common_len == 0:
            logging.error("No shared Cα atoms found; aborting.")
            sys.exit(1)
        t2.superpose(t1, atom_indices=ca2[:common_len], ref_atom_indices=ca1[:common_len])

        # ---------------------------------------------------------------------
        # RMSF
        # ---------------------------------------------------------------------
        rmsf1 = compute_rmsf(t1, ca1[:common_len])
        rmsf2 = compute_rmsf(t2, ca2[:common_len])

        avg1, std1_r = rmsf1.mean(), rmsf1.std()
        avg2, std2_r = rmsf2.mean(), rmsf2.std()
        ks_rmsf = ks_2samp(rmsf1, rmsf2)

        logging.info("RMSF Å | %s mean %.3f ± %.3f, %s mean %.3f ± %.3f, "
                     "KS‑p %.3e", args.label1, avg1, std1_r, args.label2, avg2, std2_r, ks_rmsf.pvalue)

        # Plot RMSF
        plt.figure()
        plt.plot(np.arange(common_len), rmsf1, label=args.label1)
        plt.plot(np.arange(common_len), rmsf2, label=args.label2)
        plt.xlabel("Cα index (aligned)")
        plt.ylabel("RMSF (Å)")
        plt.title(f"RMSF per residue   KS p={ks_rmsf.pvalue:.3e}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{args.out_prefix}_rmsf.pdf", dpi=300)
        plt.close()

        # ---------------------------------------------------------------------
        # Radius of gyration
        # ---------------------------------------------------------------------
        rg1 = md.compute_rg(t1) * 10.0  # nm→Å
        rg2 = md.compute_rg(t2) * 10.0
        ks_rg = ks_2samp(rg1, rg2)

        logging.info("Rg Å | %s mean %.2f ± %.2f, %s mean %.2f ± %.2f, "
                     "KS‑p %.3e", args.label1, rg1.mean(), rg1.std(), args.label2, rg2.mean(), rg2.std(),
                     ks_rg.pvalue)

        plt.figure()
        plt.plot(rg1, label=args.label1)
        plt.plot(rg2, label=args.label2)
        plt.xlabel("Frame")
        plt.ylabel("Radius of gyration (Å)")
        plt.title(f"Radius of gyration   KS p={ks_rg.pvalue:.3e}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{args.out_prefix}_rg.pdf", dpi=300)
        plt.close()

        # ---------------------------------------------------------------------
        # DSSP analysis (optional – requires DSSP binary)
        # ---------------------------------------------------------------------
        try:
            ss1 = md.compute_dssp(t1)
            ss2 = md.compute_dssp(t2)
            struct_map = {'H': 0, 'E': 1, 'C': 2, 'G': 3, 'I': 4, 'B': 5, 'T': 6, 'S': 7}
            ss_num1 = np.vectorize(struct_map.get)(ss1)
            ss_num2 = np.vectorize(struct_map.get)(ss2)

            import seaborn as sns  # optional heavy import

            plt.figure(figsize=(10, 6))
            for i, (ss, label) in enumerate([(ss_num1, args.label1),
                                             (ss_num2, args.label2)], 1):
                plt.subplot(2, 1, i)
                sns.heatmap(ss.T, cbar=i == 1, cmap="tab20",
                            cbar_kws={"label": "DSSP code"})
                plt.ylabel("Residue")
                plt.xlabel("Frame")
                plt.title(f"{label} – DSSP")
            plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_dssp.pdf", dpi=300)
            plt.close()
        except Exception as e:
            logging.warning("DSSP calculation skipped (%s)", e)

    finally:
        # Clean up any temporary decompressed files.
        for p in temp_paths:
            try:
                p.unlink()
            except Exception:
                pass


if __name__ == "__main__":
    main()
