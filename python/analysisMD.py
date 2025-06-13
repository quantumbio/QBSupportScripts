#!/usr/bin/env python3
"""
Compare two MD trajectories (Complete vs Published Structure).

Produces:
* Full‑structure RMSF profile for each trajectory
* Aligned‑residue RMSF comparison (per‑residue benchmark)
* Radius‑of‑gyration comparison
* Optional DSSP heat‑maps
* KS‑tests for all numeric comparisons

Handles .gz‑compressed inputs transparently.
"""

from __future__ import annotations
import argparse, gzip, logging, os, re, shutil, sys, tempfile
from pathlib import Path
from typing import List, Tuple

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
from Bio.Align import PairwiseAligner
from Bio.Data import IUPACData

# ──────────────────────────────────────────────────────────────────────────────
#  Configuration
# ──────────────────────────────────────────────────────────────────────────────
residue_mapping = {
    # standard
    "ALA":"ALA","ARG":"ARG","ASN":"ASN","ASP":"ASP","CYS":"CYS","GLU":"GLU","GLN":"GLN",
    "GLY":"GLY","HIS":"HIS","ILE":"ILE","LEU":"LEU","LYS":"LYS","MET":"MET","PHE":"PHE",
    "PRO":"PRO","SER":"SER","THR":"THR","TRP":"TRP","TYR":"TYR","VAL":"VAL",
    # Amber variants
    "HIE":"HIS","HID":"HIS","HIP":"HIS","CYX":"CYS","GLH":"GLU","ASH":"ASP",
}
aa3to1 = {k.upper(): v for k, v in IUPACData.protein_letters_3to1_extended.items()}

# ──────────────────────────────────────────────────────────────────────────────
#  Helper functions
# ──────────────────────────────────────────────────────────────────────────────
def decompress_if_gz(path: Path) -> Path:
    if path.suffix != ".gz":
        return path
    orig_suffix = "".join(path.stem.split(".")[1:]) or path.stem.split(".")[0]
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=f".{orig_suffix}",
                                      prefix="mdtraj_")
    logging.info("Decompressing %s → %s", path, tmp.name)
    with gzip.open(path, "rb") as gz_in, open(tmp.name, "wb") as out:
        shutil.copyfileobj(gz_in, out)
    return Path(tmp.name)

_trim3 = re.compile(r"([A-Z]{3})").match
def trim_to_three(name: str) -> str:
    m = _trim3(name.upper())
    return m.group(1) if m else name[:3].upper()

def standardise_residue_names(traj: md.Trajectory) -> List[str]:
    std = []
    for res in traj.topology.residues:
        core = residue_mapping.get(trim_to_three(res.name), trim_to_three(res.name))
        if core not in residue_mapping.values() and core not in {"HOH","NA","CL"}:
            logging.warning("Unmatched residue %s", res.name)
        std.append(core)
    return std

def load_traj(traj_path: Path, top_path: Path) -> md.Trajectory:
    t = md.load(traj_path, top=top_path)
    logging.info("Loaded %-20s  frames=%5d  atoms=%d",
                 traj_path.name, t.n_frames, t.n_atoms)
    return t

def compute_rmsf(traj: md.Trajectory, ca_idx: np.ndarray) -> np.ndarray:
    """Return RMSF in Å (MDTraj gives nm)."""
    return md.rmsf(traj, reference=traj, atom_indices=ca_idx) * 10.0

def protein_seq_and_ca(traj: md.Trajectory) -> Tuple[str, List[int]]:
    seq, ca_idx = [], []
    for res in traj.topology.residues:
        if not res.is_protein:
            continue
        letter = aa3to1.get(res.name)
        if letter is None:
            continue
        try:
            ca_atom = next(a for a in res.atoms if a.name == "CA")
        except StopIteration:
            continue
        seq.append(letter)
        ca_idx.append(ca_atom.index)
    return "".join(seq), ca_idx

def safe_name(label: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", label.strip()).lower()

# ──────────────────────────────────────────────────────────────────────────────
#  Main
# ──────────────────────────────────────────────────────────────────────────────
def main() -> None:
    p = argparse.ArgumentParser(description="Full vs aligned RMSF comparison")
    p.add_argument("traj1", type=Path); p.add_argument("top1", type=Path)
    p.add_argument("traj2", type=Path); p.add_argument("top2", type=Path)
    p.add_argument("-o","--out-prefix", default="comparison")
    p.add_argument("--label1", default="Complete Structure")
    p.add_argument("--label2", default="Published Structure")
    p.add_argument("-q","--quiet", action="store_true")
    args = p.parse_args()

    logging.basicConfig(level=logging.WARNING if args.quiet else logging.INFO,
                        format="%(levelname)s: %(message)s")

    tmp_files: List[Path] = []
    try:
        t1 = load_traj(decompress_if_gz(args.traj1), decompress_if_gz(args.top1))
        t2 = load_traj(decompress_if_gz(args.traj2), decompress_if_gz(args.top2))
        tmp_dir = Path(tempfile.gettempdir())
        tmp_files = list(tmp_dir.iterdir())

        # Standardise residue names in‑place
        for res,new in zip(t1.topology.residues, standardise_residue_names(t1)): res.name=new
        for res,new in zip(t2.topology.residues, standardise_residue_names(t2)): res.name=new

        # Slice protein‑only trajectories
        prot1 = t1.atom_slice([a.index for a in t1.topology.atoms if a.residue.is_protein])
        prot2 = t2.atom_slice([a.index for a in t2.topology.atoms if a.residue.is_protein])

        # ------------------------------------------------------------------ #
        #  FULL RMSF (each trajectory independently)
        # ------------------------------------------------------------------ #
        ca1_all = prot1.topology.select("name CA")
        ca2_all = prot2.topology.select("name CA")
        rmsf1_full = compute_rmsf(prot1, ca1_all)
        rmsf2_full = compute_rmsf(prot2, ca2_all)
        ks_full = ks_2samp(rmsf1_full, rmsf2_full)

        logging.info("FULL RMSF Å | %-18s mean %.3f ± %.3f | %-18s mean %.3f ± %.3f | KS‑p %.3e",
                     args.label1, rmsf1_full.mean(), rmsf1_full.std(),
                     args.label2, rmsf2_full.mean(), rmsf2_full.std(),
                     ks_full.pvalue)

        plt.figure()
        plt.plot(np.arange(len(rmsf1_full)), rmsf1_full, label=args.label1)
        plt.plot(np.arange(len(rmsf2_full)), rmsf2_full, label=args.label2)
        plt.xlabel("Residue index (own numbering)")
        plt.ylabel("RMSF (Å)")
        plt.title(f"Full RMSF profile   KS p={ks_full.pvalue:.3e}")
        plt.legend(); plt.tight_layout()
        plt.savefig(f"{args.out_prefix}_rmsf_full.pdf", dpi=300); plt.close()

        # ------------------------------------------------------------------ #
        #  ALIGNED‑RESIDUE RMSF
        # ------------------------------------------------------------------ #
        seq1, ca_list1 = protein_seq_and_ca(prot1)
        seq2, ca_list2 = protein_seq_and_ca(prot2)
        aln = PairwiseAligner().align(seq1, seq2)[0]
        logging.info("Global alignment score %.2f  length %d", aln.score, aln.shape[1])

        # Build matched CA index lists
        matched1, matched2 = [], []
        i1 = i2 = 0
        for a,b in zip(aln.target, aln.query):
            if a != "-" and b != "-":
                matched1.append(ca_list1[i1])
                matched2.append(ca_list2[i2])
            if a != "-": i1 += 1
            if b != "-": i2 += 1
        if not matched1:
            logging.error("No aligned residues with CA atoms found."); sys.exit(1)

        # Superpose Traj‑2 onto Traj‑1 using aligned Cα atoms
        prot2.superpose(prot1,
                        atom_indices=np.array(matched2),
                        ref_atom_indices=np.array(matched1))

        rmsf1_aln = compute_rmsf(prot1, np.array(matched1))
        rmsf2_aln = compute_rmsf(prot2, np.array(matched2))
        ks_aln = ks_2samp(rmsf1_aln, rmsf2_aln)

        logging.info("ALIGNED RMSF Å | %-18s mean %.3f ± %.3f | %-18s mean %.3f ± %.3f | KS‑p %.3e",
                     args.label1, rmsf1_aln.mean(), rmsf1_aln.std(),
                     args.label2, rmsf2_aln.mean(), rmsf2_aln.std(),
                     ks_aln.pvalue)

        plt.figure()
        x = np.arange(len(matched1))
        plt.plot(x, rmsf1_aln, label=args.label1)
        plt.plot(x, rmsf2_aln, label=args.label2)
        plt.xlabel("Aligned residue index")
        plt.ylabel("RMSF (Å)")
        plt.title(f"Aligned‑core RMSF   KS p={ks_aln.pvalue:.3e}")
        plt.legend(); plt.tight_layout()
        plt.savefig(f"{args.out_prefix}_rmsf_aligned.pdf", dpi=300); plt.close()

        # ------------------------------------------------------------------ #
        #  Radius of gyration (global)
        # ------------------------------------------------------------------ #
        rg1 = md.compute_rg(prot1)*10.0
        rg2 = md.compute_rg(prot2)*10.0
        ks_rg = ks_2samp(rg1, rg2)
        logging.info("Rg Å | %-18s %.2f ± %.2f | %-18s %.2f ± %.2f | KS‑p %.3e",
                     args.label1, rg1.mean(), rg1.std(),
                     args.label2, rg2.mean(), rg2.std(), ks_rg.pvalue)
        plt.figure()
        plt.plot(rg1, label=args.label1); plt.plot(rg2, label=args.label2)
        plt.xlabel("Frame"); plt.ylabel("Radius of gyration (Å)")
        plt.title(f"Radius of gyration   KS p={ks_rg.pvalue:.3e}")
        plt.legend(); plt.tight_layout()
        plt.savefig(f"{args.out_prefix}_rg.pdf", dpi=300); plt.close()

        # ------------------------------------------------------------------ #
        #  DSSP (optional)
        # ------------------------------------------------------------------ #
        try:
            ss1 = md.compute_dssp(prot1)
            ss2 = md.compute_dssp(prot2)
            struct_map = {'H':0,'E':1,'C':2,'G':3,'I':4,'B':5,'T':6,'S':7}
            ss_num1 = np.vectorize(struct_map.get)(ss1)
            ss_num2 = np.vectorize(struct_map.get)(ss2)
            import seaborn as sns
            plt.figure(figsize=(10,6))
            for i,(ss,label) in enumerate([(ss_num1,args.label1),(ss_num2,args.label2)],1):
                plt.subplot(2,1,i)
                sns.heatmap(ss.T,cbar=i==1,cmap="tab20",
                            cbar_kws={"label":"DSSP"})
                plt.ylabel("Residue"); plt.xlabel("Frame")
                plt.title(f"{label} – DSSP")
            plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_dssp.pdf", dpi=300); plt.close()
        except Exception as e:
            logging.warning("DSSP calculation skipped (%s)", e)

    finally:
        for f in tmp_files:
            try: f.unlink()
            except Exception: pass

if __name__ == "__main__":
    main()
