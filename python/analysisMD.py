#!/usr/bin/env python3
"""
Comprehensive comparison of two MD trajectories
===============================================
Outputs
-------
1. Full RMSF (all protein residues)                → *_rmsf_full.pdf
2. Aligned‑core RMSF                               → *_rmsf_aligned.pdf
3. Radius of gyration                              → *_rg.pdf
4. Backbone RMSD vs time                           → *_rmsd_time.pdf
5. k‑means cluster population bar‑chart            → *_cluster_populations.pdf
6. PCA (first two PCs, aligned Cα)                 → *_pca_scatter.pdf
7. Persistent hydrogen bonds (≥30 % occupancy)     → *_hbonds_persistent.csv
8. DSSP heat‑map (if DSSP binary is present)       → *_dssp.pdf
All numeric comparisons are logged with KS‑tests.
Handles .gz‑compressed inputs transparently.
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

def load_traj(traj:Path, top:Path)->md.Trajectory:
    t = md.load(traj, top=top)
    logging.info("Loaded %-20s  frames:%5d  atoms:%d", traj.name, t.n_frames, t.n_atoms)
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

# ──────────────────────────  MAIN  ───────────────────────────────────────────
def main()->None:
    ap=argparse.ArgumentParser()
    ap.add_argument("traj1"); ap.add_argument("top1")
    ap.add_argument("traj2"); ap.add_argument("top2")
    ap.add_argument("-o","--out-prefix", default="comparison")
    ap.add_argument("--label1", default="Complete Structure")
    ap.add_argument("--label2", default="Published Structure")
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

        # ───────── DSSP heat‑map (unchanged) ─────────
        try:
            ss1 = md.compute_dssp(prot1); ss2 = md.compute_dssp(prot2)
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
            plt.tight_layout(); plt.savefig(f"{args.out_prefix}_dssp.pdf",dpi=300); plt.close()
        except Exception as e:
            logging.warning("DSSP calculation skipped (%s)",e)

    finally:
        for f in tmp_files:
            try:
                f.unlink(missing_ok=True)  # Python 3.8+
            except Exception as e:
                logging.warning("Could not delete temp file %s: %s", f, e)

if __name__=="__main__":
    main()
