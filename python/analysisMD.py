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
14. OUT.gz trajectory summaries (if provided)             → printed to screen
    - Includes: number of steps, simulated ps, average energy & temperature, drift
15. Ligand Analyses (if `--ligand` is specified):
    a. Ligand RMSD vs. time                               → *_ligand_rmsd.pdf
    b. Ligand SASA vs. time                               → *_ligand_sasa.pdf
    c. Ligand per-atom RMSF                               → *_ligand_rmsf.pdf
    d. Ligand-pocket H-bond persistence                   → printed to screen
    e. Ligand-pocket contact fingerprint matrix           → printed to screen
    f. Ligand–loop DCCM statistics                        → printed to screen
16. Summary JSON file containing all key metrics          → *_summary.json

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
  - clustering: number of clusters + frame counts
  - pca: PC1/PC2 explained variance, stddev, centroid distance
  - hbonds: shared and exclusive persistent hydrogen bonds
  - top_hbonds: top persistent H-bond interactions
  - dssp_diff: residue-wise secondary structure disagreements (>50% occupancy)
  - dccm_overall: DCCM correlations and Δ-correlation highlights
  - ligand: (if applicable) ligand RMSD, RMSF, SASA, H-bonds, contact fingerprint, and DCCM
  - out1_summary / out2_summary: thermodynamic summaries parsed from OUT.gz logs

Usage Examples
--------------
Basic:
    python analysisMD.py traj1.dcd.gz top1.prmtop.gz traj2.dcd.gz top2.prmtop.gz

With OUT.gz simulation logs:
    python analysisMD.py traj1.dcd.gz top1.prmtop.gz traj2.dcd.gz top2.prmtop.gz \
        --out1 OUT1.gz --out2 OUT2.gz

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
import shutil
import csv
import logging
import re
import tempfile
from pathlib import Path
from typing import List, Tuple
import json
from datetime import datetime, timezone

# ─── Third-Party Libraries ───────────────────────────────────────────────────
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ks_2samp
from Bio.Align import PairwiseAligner
from Bio.Data import IUPACData
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

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
    return md.rmsf(traj, reference=traj, atom_indices=idx)*CONFIG["rmsf_scale"]  # nm→Å

def safe(label:str)->str:  # filename‑safe
    return re.sub(r"[^A-Za-z0-9]+","_",label.strip()).lower()

def res_label(res_or_atom):
    """Return residue label 'A:ARG42' for either an Atom or Residue."""
    res = res_or_atom.residue if hasattr(res_or_atom, "residue") else res_or_atom
    chain = chr(65 + res.chain.index)
    return f"{chain}:{res.name}{res.resSeq}"
    
def parse_out_file(out_gz_path: Path):
    """
    Parse OUT.gz file for energy and temperature summary.
    """
    steps, energies, temps = [], [], []

    try:
        with gzip.open(out_gz_path, 'rt') as fh:
            for line in fh:
                if line.strip().startswith('#') or not line.strip():
                    continue
                parts = line.strip().split(',')
                if len(parts) == 3 and parts[0].strip().isdigit():
                    try:
                        steps.append(int(parts[0]))
                        energies.append(float(parts[1]))
                        temps.append(float(parts[2]))
                    except ValueError:
                        continue
    except Exception as e:
        logging.warning("Failed to parse OUT.gz: %s", e)
        return None

    if not steps:
        return None

    interval = steps[1] - steps[0] if len(steps) > 1 else 1000
    total_ps = steps[-1] * (interval / 1000)

    return {
        "n_steps": len(steps),
        "total_ps": total_ps,
        "avg_energy": np.mean(energies),
        "std_energy": np.std(energies),
        "avg_temp": np.mean(temps),
        "std_temp": np.std(temps)
    }

def print_sim_summary(label: str, data: dict):
    if not data:
        print(f"\nSimulation Summary: {label}")
        print("  No valid OUT.gz data found.")
        return

    print(f"\nSimulation Summary: {label}")
    print(f"  Total simulated time     : {data['total_ps']:.1f} ps ({data['total_ps'] / 1000:.2f} ns)")
    print(f"  Number of steps reported : {data['n_steps']}")
    print(f"  Avg. potential energy    : {data['avg_energy']:.1f} ± {data['std_energy']:.1f} kJ/mol")
    print(f"  Avg. temperature         : {data['avg_temp']:.1f} ± {data['std_temp']:.1f} K")

    drift_ratio = data['std_energy'] / abs(data['avg_energy'])
    if drift_ratio < 0.01:
        print("  Energy drift             : Stable (std dev < 1%)")
    elif drift_ratio < 0.03:
        print("  Energy drift             : Mild instability (std dev ~2%)")
    else:
        print("  Energy drift             : Significant (std dev > 3%)")


# ──────────────────────────  MAIN  ───────────────────────────────────────────
def main()->None:
    ap=argparse.ArgumentParser()
    ap.add_argument("traj1"); ap.add_argument("top1")
    ap.add_argument("traj2"); ap.add_argument("top2")
    ap.add_argument("-o","--out-prefix", default="comparison", help="Use as a basename for the comparison (e.g. PDBid)")
    ap.add_argument("--label1", default="Complete Structure")
    ap.add_argument("--label2", default="Published Structure")
    ap.add_argument("--out1", type=Path, help="OUT.gz file for trajectory 1 (optional)")
    ap.add_argument("--out2", type=Path, help="OUT.gz file for trajectory 2 (optional)")
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
        prot2.superpose(prot1, atom_indices=np.array(m2), ref_atom_indices=np.array(m1))
        rmsf1_aln=compute_rmsf(prot1,np.array(m1)); rmsf2_aln=compute_rmsf(prot2,np.array(m2))
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
        aligned_res1 = [t1.topology.atom(ca).residue for ca in m1]
        aligned_res2 = [t2.topology.atom(ca).residue for ca in m2]
        n_aligned    = len(aligned_res1)
        
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
        n_ca_common = len(m1)
        coords1 = prot1.atom_slice(np.array(m1)).xyz.reshape(prot1.n_frames, n_ca_common * 3)
        coords2 = prot2.atom_slice(np.array(m2)).xyz.reshape(prot2.n_frames, n_ca_common * 3)
        
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
        print(f"\nCentroid distance (Complete vs Published): {centroid_dist:.2f}\n")

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
            Compute per-bond occupancy using residue labels instead of raw atom indices.
        
            Returns:
                dict: (donor_label, acceptor_label) → occupancy fraction
            """
            if not hb_frame_list:
                return {}
        
            counts = {}
            n_frames = len(hb_frame_list)
        
            for frame_array in hb_frame_list:
                for bond in frame_array:
                    donor = int(bond[0])
                    acceptor = int(bond[-1])
                    d_lbl = res_label(topology.atom(donor))
                    a_lbl = res_label(topology.atom(acceptor))
                    key = tuple(sorted((d_lbl, a_lbl)))
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
        print(f"{'Donor':<18} → {'Acceptor':<18} | {args.label1:^18} | {args.label2:^18}")
        print("-" * 70)
        
        for (don, acc), (f1, f2) in sorted_persistent[:top_n]:
            print(f"{don:<18} → {acc:<18} | {f1:>8.2f}         | {f2:>8.2f}")
        
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
            "exclusive2": only2
        }
        
        summary["top_hbonds"] = [
            {
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
            plt.title("DSSP secondary structure difference (Complete vs Published)")
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
        
        # Compute DCCMs using aligned Cα atoms
        dccm1 = compute_dccm(prot1, np.array(m1))
        dccm2 = compute_dccm(prot2, np.array(m2))
        
        # Plot DCCM matrices
        for mat, label in zip([dccm1, dccm2], [args.label1, args.label2]):
            plt.figure(figsize=CONFIG["plot"]["figsize"]["dccm"])
            sns.heatmap(mat, vmin=-1, vmax=1, cmap="coolwarm", square=True, cbar_kws={"label": "DCCM correlation"})
            plt.title(f"DCCM (aligned Cα) — {label}")
            plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_dccm_ca_{safe(label)}.pdf", dpi=CONFIG["plot"]["dpi"])
            plt.close()

        # ───────── Differential DCCM (ΔDCCM) Heatmap ─────────
        delta_dccm = dccm1 - dccm2  # ΔDCCM = Complete - Published

        plt.figure(figsize=CONFIG["plot"]["figsize"]["dccm"])
        sns.heatmap(delta_dccm, cmap="bwr", center=0, vmin=-1, vmax=1,
                    square=True, cbar_kws={"label": "ΔDCCM (Complete - Published)"})
        plt.title("Differential DCCM (Aligned Cα Atoms)")
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
                print(f"  Complete Structure: {len(ligand1)} atoms")
                print(f"  Published Structure: {len(ligand2)} atoms")
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
                Count H‑bonds where exactly one atom is in `ligand_atom_ids_in_slice`
                (indices refer to *traj*).  Returns { (don_label, acc_label) : fraction }.
                The donor/acceptor labels are sorted so direction does not create duplicates.
                """
                lig_set   = set(ligand_atom_ids_in_slice)
                n_frames  = len(traj)
                counts    = {}
            
                for frame_hb in md.wernet_nilsson(traj):
                    for hbond in frame_hb:
                        don, acc = map(int, hbond[[0, -1]])
                        in_lig = (don in lig_set, acc in lig_set)
            
                        if in_lig.count(True) == 1:          # XOR → one ligand, one protein
                            don_lbl = res_label(traj.topology.atom(don))
                            acc_lbl = res_label(traj.topology.atom(acc))
                            key     = tuple(sorted((don_lbl, acc_lbl)))  # order‑independent
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
            print(f"{'Residue (Complete)':<22} | {args.label1:^6} | {args.label2:^6}")
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
                print("  Skipped: no unique loop residues in Complete Structure.")
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
        
        # ───────── Simulation summary from OUT.gz (optional) ─────────
        if args.out1 or args.out2:
            print("\n" + "=" * 60)
            print("Production Energy & Temperature Summary")
            print("=" * 60)
        
        if args.out1:
            out1_data = parse_out_file(args.out1)
            if out1_data:
                summary["out1_summary"] = {k: float(v) for k, v in out1_data.items()}
            print_sim_summary(args.label1, out1_data)
        
        if args.out2:
            out2_data = parse_out_file(args.out2)
            if out2_data:
                summary["out2_summary"] = {k: float(v) for k, v in out2_data.items()}
            print_sim_summary(args.label2, out2_data)

        json_file = f"{args.out_prefix}_summary.json"
        with open(json_file, "w") as jfh:
            json.dump(summary, jfh, indent=2, ensure_ascii=False)
            logging.info("Wrote summary JSON to %s", json_file)

    finally:
        for f in tmp_files:
            try:
                f.unlink(missing_ok=True)  # Python 3.8+
            except Exception as e:
                logging.warning("Could not delete temp file %s: %s", f, e)

if __name__=="__main__":
    main()
