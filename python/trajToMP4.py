import sys
import os
import glob
import math
import argparse
import tempfile
import shutil
import gzip
import logging
from pathlib import Path
from pymol import cmd
from pymol.cgo import *
import mdtraj as md

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
    print(f"Decompressed {path} → {tmp_path}")
    tmp_files.append(tmp_path)
    return tmp_path

def create_gap_spheres(frame, object_name="trajectory", distance_threshold=6.0, previous_spheres=[]):
    cmd.frame(frame)
    for name in previous_spheres:
        cmd.delete(name)

    atoms = cmd.get_model(f"{object_name} and name CA", frame).atom
    prev_atom = None
    spherelist = []

    for atom in atoms:
        if prev_atom:
            dist = math.dist(atom.coord, prev_atom.coord)
            if dist > distance_threshold:
                midpoint = [(a + b) / 2 for a, b in zip(atom.coord, prev_atom.coord)]
                radius = dist / 2
                spherelist.append([
                    CYLINDER,
                    atom.coord[0], atom.coord[1], atom.coord[2],
                    prev_atom.coord[0], prev_atom.coord[1], prev_atom.coord[2],
                    radius, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0
                ])
        prev_atom = atom

    current_spheres = []
    for i, sphere in enumerate(spherelist):
        name = f'gap_sphere_{frame}_{i}'
        cmd.load_cgo(sphere, name)
        cmd.set("cgo_transparency", 0.75, name)
        current_spheres.append(name)
    return current_spheres

# ──────────────────────────  CLI ARGUMENTS  ──────────────────────────────────
parser = argparse.ArgumentParser(description="Render PyMOL movie with optional gap highlighting.")
parser.add_argument("--top", required=True, help="Topology file (.prmtop or .pdb)")
parser.add_argument("--traj", required=True, help="Trajectory file (.dcd or multi-model .pdb)")
parser.add_argument("--ligand_code", required=True, help="3-letter ligand code")
parser.add_argument("--frame", type=int, default=1, help="Stride between frames (default=1)")
parser.add_argument("--addgap", action="store_true", help="Highlight gaps with magenta cylinders")
args = parser.parse_args()

# ──────────────────────────  TRAJECTORY SETUP  ───────────────────────────────
tmp_files = []
top_path = decompress_if_gz(Path(args.top), tmp_files)
traj_path = decompress_if_gz(Path(args.traj), tmp_files)

# Convert DCD+PRMTOP to multi-model PDB if needed
if traj_path.suffix == ".dcd":
    print("Converting DCD + PRMTOP to PDB trajectory...")
    traj = md.load(traj_path, top=top_path)
    traj = traj[::args.frame]
    pdb_out = Path(tempfile.mktemp(suffix=".pdb", prefix="pymol_traj_"))
    traj.save(str(pdb_out))
    tmp_files.append(pdb_out)
    pymol_traj_file = str(pdb_out)
else:
    pymol_traj_file = str(traj_path)

# ──────────────────────────  PYMOL SETUP  ────────────────────────────────────
output_folder = "movie_frames"
output_movie = "trajectory_movie.mp4"
viewwidth = 1920
viewheight = 2160

# Clean output folder
os.makedirs(output_folder, exist_ok=True)
for f in glob.glob(f"{output_folder}/*.png"):
    os.remove(f)

# Load trajectory
cmd.load(pymol_traj_file, "trajectory")
cmd.hide("everything", "all")
cmd.show("cartoon", "polymer")
cmd.remove("resn HOH")
cmd.show("spheres", f"resn {args.ligand_code}")

cmd.dss("trajectory")
cmd.color("red", "ss H")
cmd.color("yellow", "ss S")
cmd.color("green", "ss L")

cmd.set("ray_trace_frames", 1)
cmd.set("ray_shadow", "on")
cmd.set("ambient_occlusion_mode", 1)
cmd.set("light_count", 8)
cmd.set("antialias", 2)

cmd.viewport(viewwidth, viewheight)
cmd.clip("near", 10)
cmd.clip("far", -10)

cmd.select("zoom_selection", "polymer or organic and not resn HOH+NA+CL")
cmd.zoom("zoom_selection", buffer=5)
view = cmd.get_view()

# ──────────────────────────  FRAME RENDERING  ────────────────────────────────
frame_count = cmd.count_states("trajectory")
previous_spheres = []
for frame in range(1, frame_count + 1):
    cmd.frame(frame)
    if args.addgap:
        previous_spheres = create_gap_spheres(frame, previous_spheres=previous_spheres)
    cmd.set_view(view)
    cmd.png(f"{output_folder}/frame{frame:04d}.png", width=viewwidth, height=viewheight, ray=1)
    # Uncomment to render only the first N frames
    # if frame >= 10:
    #     break

# ──────────────────────────  FFMPEG ENCODING  ────────────────────────────────
print("Encoding video with FFmpeg (2-pass)...")
os.system(f"ffmpeg -y -framerate 24 -i {output_folder}/frame%04d.png -c:v libx264 "
          f"-pix_fmt yuv420p -crf 18 -preset slow -b:v 20M -pass 1 -an -f mp4 /dev/null")
os.system(f"ffmpeg -framerate 24 -i {output_folder}/frame%04d.png -c:v libx264 "
          f"-pix_fmt yuv420p -crf 18 -preset slow -b:v 20M -pass 2 {output_movie}")

print(f"Movie saved as {output_movie}")

# ──────────────────────────  CLEANUP TEMP FILES  ─────────────────────────────
for tmp in tmp_files:
    try:
        tmp.unlink()
        print(f"Deleted temp file: {tmp}")
    except Exception as e:
        print(f"Warning: Could not delete temp file {tmp}: {e}")


# Usage Note: Upon completion of the LEFT.mp4 and the RIGHT.mp4, use the following to generate the combined:
# % ffmpeg -i LEFT.mp4 -i RIGHT.mp4 -filter_complex "\
# [0:v]scale=1920:2160,drawtext=text='MD With Gap':x=(w-text_w)/2:y=300:fontsize=96:fontcolor=white[v0]; \
# [1:v]scale=1920:2160,drawtext=text='MD Without Gap':x=(w-text_w)/2:y=300:fontsize=96:fontcolor=white[v1]; \
# [v0][v1]hstack=inputs=2[stacked]; \
# [stacked]drawtext=text='Snapshot\\: %{eif\\:2*n+1\\:d} / 999':x=(w-text_w)/2:y=h-100:fontsize=48:fontcolor=white[v]" -map "[v]" output.mp4

