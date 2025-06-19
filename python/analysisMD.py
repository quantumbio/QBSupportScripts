#!/usr/bin/env python3
"""
Comprehensive Comparison of Two MD Trajectories
===============================================

This script performs an extensive structural and energetic comparison between two
molecular dynamics (MD) trajectories — typically a "Complete Structure" vs a
"Published Structure" — to evaluate relative stability, convergence, and predictive utility.

Outputs
-------
1. Full RMSF (all protein residues)                    → *_rmsf_full.pdf
2. Aligned-core RMSF (based on sequence alignment)     → *_rmsf_aligned.pdf
3. Radius of gyration (Rg) over time                   → *_rg.pdf
4. Backbone RMSD vs time                               → *_rmsd_time.pdf
5. k-means cluster population bar chart (CA-RMSD)      → *_cluster_populations.pdf
6. PCA scatter plot (PC1 vs PC2, aligned Cα)           → *_pca_scatter.pdf
7. PCA variance and spread summary                     → printed to screen
8. Persistent hydrogen bonds (≥30% occupancy)          → *_hbonds_persistent.csv
   - Includes top 20 persistent H-bonds printed to screen
   - Summary of shared vs exclusive bonds also printed
9. DSSP secondary structure heatmaps                   → *_dssp.pdf
10. DSSP difference map (per-residue disagreement)     → *_dssp_difference.pdf
11. OUT.gz production summary (if provided)            → printed to screen
    - Includes simulated time, avg. energy, temperature, and energy drift
12. Ligand Stability Analysis (if --ligand supplied)   → printed to screen and saved as:
    a. Ligand RMSD over time                           → *_ligand_rmsd.pdf
    b. Ligand SASA over time                           → *_ligand_sasa.pdf
    c. Ligand per-atom RMSF                            → *_ligand_rmsf.pdf
    d. Ligand-pocket hydrogen bonds (≥30%)             → printed to screen
    e. Ligand-residue contact fingerprint map          → printed to screen
    f. All ligand data added to *_summary.json
13. JSON summary for all metrics                       → *_summary.json
    - See field descriptions below

JSON Output Description
-----------------------
The *_summary.json file contains a structured record of all computed metrics:

{
  "label1": Name of trajectory 1 (e.g., "Complete Structure"),
  "label2": Name of trajectory 2 (e.g., "Published Structure"),

  "rmsf": {
    "full":     { "mean1", "std1", "mean2", "std2", "ks_p" },
    "aligned":  { "mean1", "std1", "mean2", "std2", "ks_p" }
  },

  "rg": {
    "mean1", "std1", "mean2", "std2", "ks_p"
  },

  "pca": {
    "explained_variance": [ PC1%, PC2% ],
    "std_pc1": [ traj1_std, traj2_std ],
    "std_pc2": [ traj1_std, traj2_std ],
    "centroid_distance": Euclidean distance between PCA centroids
  },

  "clustering": {
    "k": number of clusters selected,
    "populations": {
      "label1": [ cluster sizes in traj1 ],
      "label2": [ cluster sizes in traj2 ]
    }
  },

  "hbonds": {
    "n1": persistent H-bonds in traj1 (≥30%),
    "n2": persistent H-bonds in traj2 (≥30%),
    "shared": number of shared persistent H-bonds,
    "exclusive1": count only in traj1,
    "exclusive2": count only in traj2
  },

  "top_hbonds": [  # sorted by highest max persistence
    {
      "donor":    "A:ARG61",
      "acceptor": "A:GLU35",
      "label1":   1.00,
      "label2":   1.00
    },
    ...
  ],

  "dssp_diff": {
    "residues_differing_gt_50pct": Number of aligned residues whose secondary structure disagrees >50% of the time
  },

  "out1_summary": {
    "total_ps": total simulated time (ps),
    "avg_energy": mean potential energy (kJ/mol),
    "std_energy": std dev of potential energy,
    "avg_temp": mean temperature (K),
    "std_temp": temperature fluctuation
  },

  "out2_summary": {
    (same fields as out1_summary)
  },

  "ligand": {   # Included only if --ligand is provided AND the ligand is present in both structures
    "resname": 3-letter residue name (e.g., "CBN"),

    "rmsd": {
      "label1": { "mean", "std", "max" },
      "label2": { "mean", "std", "max" }
    },

    "sasa": {
      "label1": { "mean", "std" },
      "label2": { "mean", "std" }
    },

    "rmsf": {
      "label1": { "mean", "std", "max" },
      "label2": { "mean", "std", "max" }
    },

    "hbond_persistence": [
      {
        "donor": "A:ARG32",
        "acceptor": "B:GLU45",
        "label1": 0.95,
        "label2": 0.00
      },
      ...
    ],

    "contact_fingerprint": [
      {
        "residue": "A:ASP17",
        "label1": 0.44,
        "label2": 0.12
      },
      ...
    ]
  }
}

Usage
-----
Basic:
    python analysisMD.py traj1.dcd.gz top1.prmtop.gz traj2.dcd.gz top2.prmtop.gz

With optional OUT.gz files (adds simulation summary):
    python analysisMD.py traj1.dcd.gz top1.prmtop.gz traj2.dcd.gz top2.prmtop.gz \
        --out1 OUT1.gz --out2 OUT2.gz

With ligand analysis:
    python analysisMD.py traj1.dcd.gz top1.prmtop.gz traj2.dcd.gz top2.prmtop.gz \
        --ligand LIG

Optional labels for reporting:
    --label1 "Complete Structure" --label2 "Published Structure"

Output prefix (e.g., "1kzn_comparison"):
    -o 1kzn_comparison

Developer Notes
---------------
- All temporary decompressed files are tracked and deleted on exit.
- Sequence alignment is performed using Biopython’s PairwiseAligner.
- PCA and clustering use scikit-learn (requires installation).
- Hydrogen bond persistence is based on MDTraj’s Wernet-Nilsson criteria.
- DSSP analysis requires the DSSP binary (e.g., `conda install -c salilab dssp`).
- OUT.gz files are parsed as CSV with 3 columns: step, potential energy, temperature.
- The ligand section requires --ligand RESNAME to be provided and found in both structures.
- The *_summary.json file is suitable for downstream batch analysis.

Dependencies:
    mdtraj, numpy, matplotlib, seaborn, Bio, scikit-learn, gzip

Recommended:
    Use this script to benchmark MD stability, convergence, and reliability prior
    to ensemble averaging, binding free energy calculations, or scoring pipelines.
"""

from __future__ import annotations
import argparse, gzip, logging, os, re, shutil, sys, tempfile, csv
from pathlib import Path
from typing import List, Tuple

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
from Bio.Align import PairwiseAligner
from Bio.Data import IUPACData
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import json

# ─────────────────────────────  CONFIG  ──────────────────────────────────────
residue_mapping = {             # standard 3‑letter + Amber variants
    **{aa: aa for aa in
       ("ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS",
        "MET","PHE","PRO","SER","THR","TRP","TYR","VAL")},
    "HIE":"HIS","HID":"HIS","HIP":"HIS","CYX":"CYS","GLH":"GLU","ASH":"ASP",
}
aa3to1 = {k.upper(): v for k, v in IUPACData.protein_letters_3to1_extended.items()}
_trim3 = re.compile(r"([A-Z]{3})").match

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
    return md.rmsf(traj, reference=traj, atom_indices=idx)*10.0  # nm→Å

def safe(label:str)->str:  # filename‑safe
    return re.sub(r"[^A-Za-z0-9]+","_",label.strip()).lower()
    
def parse_out_file(out_gz_path: Path):
    """
    Parse OUT.gz file for energy and temperature summary.
    """
    import gzip
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
    ap.add_argument("-o","--out-prefix", default="comparison")
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
        plt.savefig(f"{args.out_prefix}_rmsf_full.pdf",dpi=300); plt.close()

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
        plt.savefig(f"{args.out_prefix}_rmsf_aligned.pdf",dpi=300); plt.close()

        # Aligned residue lists (built earlier from m1 / m2)
        aligned_res1 = [t1.topology.atom(ca).residue for ca in m1]
        aligned_res2 = [t2.topology.atom(ca).residue for ca in m2]
        n_aligned    = len(aligned_res1)

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
        plt.savefig(f"{args.out_prefix}_rg.pdf",dpi=300); plt.close()

        # ───────── BACKBONE RMSD vs TIME ─────────
        bb_idx1 = prot1.topology.select("backbone")
        bb_idx2 = prot2.topology.select("backbone")
        rmsd1 = md.rmsd(prot1, prot1, 0, bb_idx1)*10.0
        rmsd2 = md.rmsd(prot2, prot2, 0, bb_idx2)*10.0
        plt.figure()
        plt.plot(rmsd1,label=args.label1); plt.plot(rmsd2,label=args.label2)
        plt.xlabel("Frame"); plt.ylabel("Backbone RMSD to start (Å)")
        plt.title("Backbone RMSD vs time"); plt.legend(); plt.tight_layout()
        plt.savefig(f"{args.out_prefix}_rmsd_time.pdf",dpi=300); plt.close()

        # ───────── CLUSTERING (k‑means, CA‑RMSD feature matrix) ─────────
        min_frames = 10
        k_opt, pop1, pop2 = None, [], []
        if prot1.n_frames >= min_frames and prot2.n_frames >= min_frames \
           and len(ca1_all) and len(ca2_all):
        
            ref_frames = min(50, prot1.n_frames, prot2.n_frames)      # up to 50 refs
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
            max_k     = min(10, n_samples)    # never ask for more clusters than samples
        
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
            plt.savefig(f"{args.out_prefix}_cluster_populations.pdf", dpi=300)
            plt.close()
        else:
            k_opt = None
            pop1 = pop2 = []
            logging.warning("Clustering skipped: insufficient frames or atoms.")

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
        plt.savefig(f"{args.out_prefix}_pca_scatter.pdf", dpi=300)
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

        # ───────── Hydrogen‑bond persistence ─────────
        hb1 = md.wernet_nilsson(prot1, periodic=False, sidechain_only=False)
        hb2 = md.wernet_nilsson(prot2, periodic=False, sidechain_only=False)
        def occupancy(hb_frame_list):
            """
            Compute per‑bond occupancy fractions.
        
            Returns
            -------
            dict
                key   = (donor_atom_index, acceptor_atom_index)
                value = occupancy fraction (0–1)
            """
            if not hb_frame_list:
                return {}
        
            n_frames = len(hb_frame_list)
            counts = {}
        
            for frame_array in hb_frame_list:            # one array per frame
                for bond in frame_array:                 # bond length = 2 or 3
                    donor = int(bond[0])
                    acceptor = int(bond[-1])             # last column is always acceptor
                    key = (donor, acceptor)
                    counts[key] = counts.get(key, 0) + 1
        
            return {k: v / n_frames for k, v in counts.items()}
        occ1=occupancy(hb1); occ2=occupancy(hb2)
        persistent = {k:(occ1.get(k,0.0), occ2.get(k,0.0))
                      for k in set(occ1)|set(occ2)
                      if max(occ1.get(k,0.0), occ2.get(k,0.0))>=0.30}
        if persistent:
            csvfile=f"{args.out_prefix}_hbonds_persistent.csv"
            with open(csvfile,"w",newline="") as fh:
                w=csv.writer(fh); w.writerow(["Donor","Acceptor",args.label1,args.label2])
                for (don_idx, acc_idx), (f1, f2) in sorted(persistent.items(),
                                                           key=lambda x: -max(x[1])):
                    don_res = prot1.topology.atom(don_idx).residue
                    acc_res = prot1.topology.atom(acc_idx).residue
                    w.writerow([f"{don_res}-{don_idx}", f"{acc_res}-{acc_idx}",
                                f"{f1:.2f}", f"{f2:.2f}"])
            logging.info("Wrote persistent H‑bonds to %s (%d entries)",csvfile,len(persistent))
        else:
            logging.info("No hydrogen bonds ≥30%% occupancy in either trajectory")

        # Print summary table to screen
        top_n = 20
        sorted_persistent = sorted(persistent.items(), key=lambda x: -max(x[1]))
        
        print("\nTop persistent hydrogen bonds (≥30% occupancy in either trajectory):")
        print(f"{'Donor':<18} → {'Acceptor':<18} | {args.label1:^18} | {args.label2:^18}")
        print("-" * 70)
        
        for (don_idx, acc_idx), (f1, f2) in sorted(persistent.items(), key=lambda x: -max(x[1]))[:20]:
            top_d = prot1 if don_idx < prot1.n_atoms else prot2
            top_a = prot1 if acc_idx < prot1.n_atoms else prot2
        
            don_atom = top_d.topology.atom(don_idx)
            acc_atom = top_a.topology.atom(acc_idx)
        
            don_res = don_atom.residue
            acc_res = acc_atom.residue
        
            don_str = f"{chr(65 + don_res.chain.index)}:{don_res.name}{don_res.resSeq}"
            acc_str = f"{chr(65 + acc_res.chain.index)}:{acc_res.name}{acc_res.resSeq}"
        
            print(f"{don_str:<18} → {acc_str:<18} | {f1:>8.2f}         | {f2:>8.2f}")
        
        # Summary stats
        n1 = sum(1 for v in persistent.values() if v[0] >= 0.30)
        n2 = sum(1 for v in persistent.values() if v[1] >= 0.30)
        shared = sum(1 for v in persistent.values() if v[0] >= 0.30 and v[1] >= 0.30)
        only1 = sum(1 for v in persistent.values() if v[0] >= 0.30 and v[1] < 0.30)
        only2 = sum(1 for v in persistent.values() if v[1] >= 0.30 and v[0] < 0.30)
        
        print(f"\n{args.label1}: {n1} persistent H-bonds (≥30%)")
        print(f"{args.label2}: {n2} persistent H-bonds (≥30%)")
        print(f"Shared: {shared}")
        print(f"Exclusive to {args.label1}: {only1}")
        print(f"Exclusive to {args.label2}: {only2}\n")

        # ───────── DSSP heat‑map ─────────
        n_disagree = 0
        try:
            ss1 = md.compute_dssp(prot1)
            ss2 = md.compute_dssp(prot2)
            struct_map={'H':0,'E':1,'C':2,'G':3,'I':4,'B':5,'T':6,'S':7}
            ss_num1=np.vectorize(struct_map.get)(ss1)
            ss_num2=np.vectorize(struct_map.get)(ss2)
        
            import seaborn as sns
            plt.figure(figsize=(10,6))
            for i,(ss,label) in enumerate([(ss_num1,args.label1),(ss_num2,args.label2)],1):
                plt.subplot(2,1,i)
                sns.heatmap(ss.T,cbar=i==1,cmap="tab20",
                            cbar_kws={"label":"DSSP"})
                plt.ylabel("Residue"); plt.xlabel("Frame")
                plt.title(f"{label} – DSSP")
            plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_dssp.pdf", dpi=300)
            plt.close()
        
            # Difference map (on aligned residues only)
            r1_idx = [prot1.topology.atom(i).residue.index for i in m1]
            r2_idx = [prot2.topology.atom(i).residue.index for i in m2]
            aligned_ss1 = ss1[:, r1_idx]
            aligned_ss2 = ss2[:, r2_idx]
            dssp_diff = (aligned_ss1 != aligned_ss2).astype(np.int_)
        
            plt.figure(figsize=(12, 6))
            sns.heatmap(dssp_diff.T, cmap="Reds", cbar_kws={"label": "Mismatch (1 = different)"})
            plt.xlabel("Frame")
            plt.ylabel("Aligned residue index")
            plt.title("DSSP secondary structure difference (Complete vs Published)")
            plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_dssp_difference.pdf", dpi=300)
            plt.close()
        
            diff_summary = dssp_diff.mean(axis=0)
            n_disagree = np.sum(diff_summary > 0.5)
            print(f"\n{n_disagree} aligned residues differed in secondary structure "
                  f"in >50% of frames.\n")
        
        except Exception as e:
            logging.warning("DSSP calculation failed (%s)", e)

        # ───────── DCCM (Dynamic Cross-Correlation Matrix) for Aligned Cα Atoms ─────────
        TOP_N_DCCM = 5  # Number of top correlated and anti-correlated residue pairs to report
        
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
            plt.figure(figsize=(8, 6))
            sns.heatmap(mat, vmin=-1, vmax=1, cmap="coolwarm", square=True, cbar_kws={"label": "DCCM correlation"})
            plt.title(f"DCCM (aligned Cα) — {label}")
            plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_dccm_ca_{safe(label)}.pdf", dpi=300)
            plt.close()

        # ───────── Differential DCCM (ΔDCCM) Heatmap ─────────
        delta_dccm = dccm1 - dccm2  # ΔDCCM = Complete - Published

        plt.figure(figsize=(8, 6))
        sns.heatmap(delta_dccm, cmap="bwr", center=0, vmin=-1, vmax=1,
                    square=True, cbar_kws={"label": "ΔDCCM (Complete - Published)"})
        plt.title("Differential DCCM (Aligned Cα Atoms)")
        plt.tight_layout()
        plt.savefig(f"{args.out_prefix}_dccm_delta.pdf", dpi=300)
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

        dccm_stats = {
            "mean_delta"      : float(avg_corr_overall),
            "max_abs_delta"   : float(max_diff),
            "top_positive"    : {
                args.label1: [ {"res1": a, "res2": b, "corr": float(c)}
                               for a, b, c in top_correlated_residues(
                                   dccm1, aligned_res1, top_n=TOP_N_DCCM, sign="positive") ],
                args.label2: [ {"res1": a, "res2": b, "corr": float(c)}
                               for a, b, c in top_correlated_residues(
                                   dccm2, aligned_res2, top_n=TOP_N_DCCM, sign="positive") ],
            },
            "top_negative"    : {
                args.label1: [ {"res1": a, "res2": b, "corr": float(c)}
                               for a, b, c in top_correlated_residues(
                                   dccm1, aligned_res1, top_n=TOP_N_DCCM, sign="negative") ],
                args.label2: [ {"res1": a, "res2": b, "corr": float(c)}
                               for a, b, c in top_correlated_residues(
                                   dccm2, aligned_res2, top_n=TOP_N_DCCM, sign="negative") ],
            }
        }

        # ───────── Ligand detection or reporting ─────────
        ligand_resname = args.ligand.upper() if args.ligand else None
        
        if ligand_resname:
            ligand1 = t1.topology.select(f"resname {ligand_resname}")
            ligand2 = t2.topology.select(f"resname {ligand_resname}")
        
            if len(ligand1) == 0 or len(ligand2) == 0:
                print(f"\nWARNING: Ligand '{ligand_resname}' not found in one or both trajectories.\n")
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
        
            # Slice trajectories
            lig_traj1 = t1.atom_slice(ligand1)
            lig_traj2 = t2.atom_slice(ligand2)
        
            # ───── Ligand RMSD (vs. frame 0) ─────
            lig_rmsd1 = md.rmsd(lig_traj1, lig_traj1, 0) * 10.0  # Å
            lig_rmsd2 = md.rmsd(lig_traj2, lig_traj2, 0) * 10.0  # Å
        
            plt.figure()
            plt.plot(lig_rmsd1, label=f"{args.label1}")
            plt.plot(lig_rmsd2, label=f"{args.label2}")
            plt.xlabel("Frame"); plt.ylabel("Ligand RMSD (Å)")
            plt.title(f"Ligand RMSD (resname {ligand_resname})")
            plt.legend(); plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_ligand_rmsd.pdf", dpi=300)
            plt.close()
        
            print(f"  Ligand RMSD (Å):")
            print(f"    {args.label1:<20} mean={lig_rmsd1.mean():.2f}  std={lig_rmsd1.std():.2f}  max={lig_rmsd1.max():.2f}")
            print(f"    {args.label2:<20} mean={lig_rmsd2.mean():.2f}  std={lig_rmsd2.std():.2f}  max={lig_rmsd2.max():.2f}")
        
            # ───── Ligand SASA (Shrake-Rupley) ─────
            sasa1 = md.shrake_rupley(lig_traj1).sum(axis=1)
            sasa2 = md.shrake_rupley(lig_traj2).sum(axis=1)
        
            plt.figure()
            plt.plot(sasa1, label=args.label1)
            plt.plot(sasa2, label=args.label2)
            plt.xlabel("Frame"); plt.ylabel("Ligand SASA (Å²)")
            plt.title(f"Ligand Solvent Accessibility (resname {ligand_resname})")
            plt.legend(); plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_ligand_sasa.pdf", dpi=300)
            plt.close()
        
            print(f"  Ligand SASA (Å²):")
            print(f"    {args.label1:<20} mean={sasa1.mean():.1f}  std={sasa1.std():.1f}")
            print(f"    {args.label2:<20} mean={sasa2.mean():.1f}  std={sasa2.std():.1f}")
        
            # ───── Ligand RMSF (per atom) ─────
            rmsf_lig1 = md.rmsf(lig_traj1, reference=lig_traj1) * 10.0
            rmsf_lig2 = md.rmsf(lig_traj2, reference=lig_traj2) * 10.0
        
            plt.figure()
            plt.plot(rmsf_lig1, label=args.label1)
            plt.plot(rmsf_lig2, label=args.label2)
            plt.xlabel("Ligand Atom Index"); plt.ylabel("RMSF (Å)")
            plt.title(f"Ligand Atom Flexibility (RMSF)")
            plt.legend(); plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_ligand_rmsf.pdf", dpi=300)
            plt.close()

            print("  Ligand RMSF:")
            print(f"    {args.label1:<20} mean={rmsf_lig1.mean():.2f}  std={rmsf_lig1.std():.2f}  max={rmsf_lig1.max():.2f}")
            print(f"    {args.label2:<20} mean={rmsf_lig2.mean():.2f}  std={rmsf_lig2.std():.2f}  max={rmsf_lig2.max():.2f}")

            # ───── Ligand–Pocket H-bonds (≥30%) [Aligned Structures] ─────
            print("\nLigand–Pocket H-bond Analysis")
            print("=" * 30)
            
            # Union of protein + ligand (from ALIGNED prot1/prot2)
            ligand_protein1 = np.unique(np.concatenate([prot1.topology.select("all"), ligand1]))
            ligand_protein2 = np.unique(np.concatenate([prot2.topology.select("all"), ligand2]))
            
            traj_lp1 = t1.atom_slice(ligand_protein1)
            traj_lp2 = t2.atom_slice(ligand_protein2)
            
            hb1_all = md.wernet_nilsson(traj_lp1)
            hb2_all = md.wernet_nilsson(traj_lp2)
            
            def extract_ligand_hbonds(hb_frame_list, ligand_slice_indices):
                ligand_set = set(ligand_slice_indices)
                counts = {}
                n_frames = len(hb_frame_list)
                for frame in hb_frame_list:
                    for hbond in frame:
                        don, acc = map(int, hbond[[0, -1]])
                        if don in ligand_set or acc in ligand_set:
                            key = (don, acc)
                            counts[key] = counts.get(key, 0) + 1
                return {k: v / n_frames for k, v in counts.items()}
            
            # Map ligand atom indices to sliced traj index space
            ligand1_map = np.where(np.isin(ligand_protein1, ligand1))[0]
            ligand2_map = np.where(np.isin(ligand_protein2, ligand2))[0]
            
            lig_hb1 = extract_ligand_hbonds(hb1_all, ligand1_map)
            lig_hb2 = extract_ligand_hbonds(hb2_all, ligand2_map)
            
            # Merge and filter shared H-bonds
            lig_hb_shared = {
                k: (lig_hb1.get(k, 0), lig_hb2.get(k, 0))
                for k in set(lig_hb1) | set(lig_hb2)
                if max(lig_hb1.get(k, 0), lig_hb2.get(k, 0)) >= 0.3
            }
            
            if lig_hb_shared:
                print(f"  {len(lig_hb_shared)} ligand–protein H-bonds ≥30% persistence")
                for (don, acc), (f1, f2) in sorted(lig_hb_shared.items(), key=lambda x: -max(x[1])):
                    # Use aligned trajectories for residue labeling
                    don_atom = traj_lp1.topology.atom(don) if don < traj_lp1.n_atoms else traj_lp2.topology.atom(don)
                    acc_atom = traj_lp1.topology.atom(acc) if acc < traj_lp1.n_atoms else traj_lp2.topology.atom(acc)
                    r1 = don_atom.residue
                    r2 = acc_atom.residue
                    name1 = f"{chr(65 + r1.chain.index)}:{r1.name}{r1.resSeq}"
                    name2 = f"{chr(65 + r2.chain.index)}:{r2.name}{r2.resSeq}"
                    print(f"    {name1} ↔ {name2} | {args.label1}: {f1:.2f}  {args.label2}: {f2:.2f}")
            else:
                print("  No persistent ligand–pocket H-bonds found.")

            # ───── Contact Fingerprint (aligned residues, ligand ↔ protein) ─────
            cutoff_nm = 0.35                              # 3.5 Å
                        
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
            print("\nLigand–Residue Contact Fingerprint (aligned residues; contact ≥3.5 Å)")
            print(f"{'Residue (Complete)':<22} | {args.label1:^6} | {args.label2:^6}")
            print("-"*44)
            for i, (r1, r2) in enumerate(zip(aligned_res1, aligned_res2)):
                if max(fp1[i], fp2[i]) >= 0.30:           # show only “interesting” residues
                    name1 = f"{chr(65 + r1.chain.index)}:{r1.name}{r1.resSeq}"
                    name2 = f"{chr(65 + r2.chain.index)}:{r2.name}{r2.resSeq}"
                    print(f"{name1:<22} | {fp1[i]:>5.2f} | {fp2[i]:>5.2f}   ({name2})")

            # ───── Ligand–Loop DCCM (Complete Structure, residue‑level) ─────
            print("\nLigand–Loop DCCM Analysis")
            print("=" * 30)

            # Set how many top residues to report
            N_TOP = 5

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

                def res_label(res):
                    return f"{chr(65 + res.chain.index)}:{res.name}{res.resSeq}"

                # 3. build displacement matrix & DCCM
                disp_mat = residue_displacements(t1, all_residues)
                dccm_res = dccm_from_disp(disp_mat)             # (n_res × n_res)

                # 4. extract ligand (row 0) ↔ loop block
                lig_loop_corr = dccm_res[0, 1:]                  # vector length = len(loop_residues)
                avg_corr = lig_loop_corr.mean()
                max_corr = np.abs(lig_loop_corr).max()
                idx_max  = np.argmax(np.abs(lig_loop_corr))

                # 5. save figure
                plt.figure(figsize=(5, 4))
                plt.imshow(dccm_res, vmin=-1, vmax=1, cmap="bwr")
                plt.colorbar(label="Cross‑correlation")
                plt.title("Residue‑level DCCM  (CBN ↔ new loops)")
                plt.tight_layout()
                plt.savefig(f"{args.out_prefix}_dccm_ligand_loop_residue.pdf", dpi=300)
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

        summary = {
            "label1": args.label1,
            "label2": args.label2,
            "rmsf": {
                "full": {
                    "mean1": float(rmsf1_full.mean()), "std1": float(rmsf1_full.std()),
                    "mean2": float(rmsf2_full.mean()), "std2": float(rmsf2_full.std()),
                    "ks_p": float(ks_full.pvalue)
                },
                "aligned": {
                    "mean1": float(rmsf1_aln.mean()), "std1": float(rmsf1_aln.std()),
                    "mean2": float(rmsf2_aln.mean()), "std2": float(rmsf2_aln.std()),
                    "ks_p": float(ks_aln.pvalue)
                }
            },
            "rg": {
                "mean1": float(rg1.mean()), "std1": float(rg1.std()),
                "mean2": float(rg2.mean()), "std2": float(rg2.std()),
                "ks_p": float(ks_rg.pvalue)
            },
            "pca": {
                "explained_variance": [float(v) for v in var_exp],
                "std_pc1": [float(sd1_pc1), float(sd2_pc1)],
                "std_pc2": [float(sd1_pc2), float(sd2_pc2)],
                "centroid_distance": float(centroid_dist)
            },
            "clustering": {
                "k": k_opt,
                "populations": {
                    args.label1: list(map(int, pop1)),
                    args.label2: list(map(int, pop2))
                }
            },
            "hbonds": {
                "n1": n1, "n2": n2,
                "shared": shared,
                "exclusive1": only1,
                "exclusive2": only2
            },
            "top_hbonds" : [
                {
                    "donor": f"{chr(65 + don_res.chain.index)}:{don_res.name}{don_res.resSeq}",
                    "acceptor": f"{chr(65 + acc_res.chain.index)}:{acc_res.name}{acc_res.resSeq}",
                    args.label1: float(f1),
                    args.label2: float(f2)
                }
                for (don_idx, acc_idx), (f1, f2) in sorted_persistent
                for don_res, acc_res in [
                    (
                        (prot1 if don_idx < prot1.n_atoms else prot2).topology.atom(don_idx).residue,
                        (prot1 if acc_idx < prot1.n_atoms else prot2).topology.atom(acc_idx).residue
                    )
                ]
            ],
            "dssp_diff": {
                "residues_differing_gt_50pct": int(n_disagree)
            },
            "dccm_overall" : dccm_stats
        }

        # ───── Add Ligand Data to JSON Summary (if analyzed) ─────
        if ligand_resname and len(ligand1) > 0 and len(ligand2) > 0:
            summary["ligand"] = {
                "resname": ligand_resname,
                "rmsd": {
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
                },
                "sasa": {
                    args.label1: {
                        "mean": float(sasa1.mean()),
                        "std": float(sasa1.std())
                    },
                    args.label2: {
                        "mean": float(sasa2.mean()),
                        "std": float(sasa2.std())
                    }
                },
                "rmsf": {
                    args.label1: {
                        "mean": float(rmsf_lig1.mean()),
                        "std": float(rmsf_lig1.std()),
                        "max": float(rmsf_lig1.max())
                    },
                    args.label2: {
                        "mean": float(rmsf_lig2.mean()),
                        "std": float(rmsf_lig2.std()),
                        "max": float(rmsf_lig2.max())
                    }
                },
                "hbond_persistence": [
                    {
                        "donor": f"{chr(65 + t1.topology.atom(don).residue.chain.index)}:{t1.topology.atom(don).residue.name}{t1.topology.atom(don).residue.resSeq}" if don < t1.n_atoms else
                                 f"{chr(65 + t2.topology.atom(don).residue.chain.index)}:{t2.topology.atom(don).residue.name}{t2.topology.atom(don).residue.resSeq}",
                        "acceptor": f"{chr(65 + t1.topology.atom(acc).residue.chain.index)}:{t1.topology.atom(acc).residue.name}{t1.topology.atom(acc).residue.resSeq}" if acc < t1.n_atoms else
                                    f"{chr(65 + t2.topology.atom(acc).residue.chain.index)}:{t2.topology.atom(acc).residue.name}{t2.topology.atom(acc).residue.resSeq}",
                        args.label1: float(f1),
                        args.label2: float(f2)
                    }
                    for (don, acc), (f1, f2) in lig_hb_shared.items()
                ],
                "contact_fingerprint": [
                    {
                        "residue": f"{chr(65 + res.chain.index)}:{res.name}{res.resSeq}",
                        args.label1: float(fp1[i]),
                        args.label2: float(fp2[i])
                    }
                    for i, res in enumerate(aligned_res1)  # use aligned_res1 here
                    if max(fp1[i], fp2[i]) > 0.3  # keep only significant contacts
                ]
            }
        if lig_loop_stats is not None:
            # ensure "ligand" key exists (created earlier)
            summary.setdefault("ligand", {})["loop_dccm"] = lig_loop_stats
        
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


        with open(f"{args.out_prefix}_summary.json", "w") as jfh:
            json.dump(summary, jfh, indent=2)
            logging.info("Wrote summary JSON to %s", f"{args.out_prefix}_summary.json")


    finally:
        for f in tmp_files:
            try:
                f.unlink(missing_ok=True)  # Python 3.8+
            except Exception as e:
                logging.warning("Could not delete temp file %s: %s", f, e)

if __name__=="__main__":
    main()
