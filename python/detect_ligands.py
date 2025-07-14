import mdtraj as md
import gzip
import shutil
import tempfile
import sys
from pathlib import Path

# Known protein residues (includes protonation variants)
KNOWN_PROTEIN_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",
    "HIS", "HIE", "HIP", "HID", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "ASH", "GLH", "CYX", "CYM"
}

def decompress_if_gz(path: Path, tmp_files: list[Path]) -> Path:
    if path.suffix != ".gz":
        return path
    inner_ext = "".join(Path(path.stem).suffixes) or ".dat"
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=inner_ext, prefix="mdtraj_")
    with gzip.open(path, "rb") as gz_in, open(tmp.name, "wb") as out:
        shutil.copyfileobj(gz_in, out)
    tmp_files.append(Path(tmp.name))
    return Path(tmp.name)

def get_ligand_resnames(top):
    return sorted(set(
        res.name[:3] for res in top.residues
        if not res.is_water
        and res.name[:3] not in KNOWN_PROTEIN_RESNAMES
        and not res.name.startswith("NA")  # skip ions like NA+, CL-
        and not res.name.startswith("CL")
    ))

def main():
    if len(sys.argv) != 5:
        print("Usage: python detect_ligands.py traj1.dcd.gz top1.prmtop.gz traj2.dcd.gz top2.prmtop.gz")
        sys.exit(1)

    tmp_files = []
    try:
        traj1_path = decompress_if_gz(Path(sys.argv[1]), tmp_files)
        top1_path = decompress_if_gz(Path(sys.argv[2]), tmp_files)
        traj2_path = decompress_if_gz(Path(sys.argv[3]), tmp_files)
        top2_path = decompress_if_gz(Path(sys.argv[4]), tmp_files)

        t1 = md.load(traj1_path, top=top1_path)
        t2 = md.load(traj2_path, top=top2_path)

        lig1 = get_ligand_resnames(t1.topology)
        lig2 = get_ligand_resnames(t2.topology)
        common = sorted(set(lig1) | set(lig2))

        if common:
            print("\nLigand not specified.")
            print("These are the normalized non-protein 3-letter residue codes found:")
            print("  " + " ".join(common) + "\n")
        else:
            print("\nNo non-protein residues found.\n")

    finally:
        for f in tmp_files:
            try:
                f.unlink()
            except Exception:
                pass

if __name__ == "__main__":
    main()
