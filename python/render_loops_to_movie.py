import os
import sys
import glob
import time
import random
import argparse
import subprocess
import math
from itertools import cycle, islice
from pymol import cmd

# ──────────────── ARGUMENTS ────────────────
parser = argparse.ArgumentParser(description="Render loop ensembles (left/right + serial) into PyMOL movie.")
parser.add_argument("--main_pdb",    required=True, help="Path to the main protein structure (e.g., main.pdb)")
parser.add_argument("--left_glob",   required=True, help="Glob for left fragments (e.g., 'loops/left/**/*.pdb')")
parser.add_argument("--right_glob",  required=True, help="Glob for right fragments (e.g., 'loops/right/**/*.pdb')")
parser.add_argument("--serial_glob", required=True, help="Glob for serial fragments (e.g., 'loops/serial/**/*.pdb')")
parser.add_argument("--output",      default="loop_ensemble_movie.mp4", help="Output MP4 filename")
parser.add_argument("--output_dir",  default="movie_frames",            help="Directory to save PNG frames")
parser.add_argument("--limit",       type=int, default=None,            help="If set, only use first N fragments for quick testing")
args = parser.parse_args()

# ──────────────── FILE DISCOVERY ────────────────
main_pdb     = args.main_pdb
left_files   = sorted(glob.glob(args.left_glob,   recursive=True))
right_files  = sorted(glob.glob(args.right_glob,  recursive=True))
serial_files = sorted(glob.glob(args.serial_glob, recursive=True))

if not os.path.exists(main_pdb):
    sys.exit("❌ main.pdb not found")
if not left_files or not right_files or not serial_files:
    sys.exit("❌ One or more fragment sets not found")

# optional limit
if args.limit:
    left_files   = left_files[:2*args.limit]
    right_files  = right_files[:2*args.limit]
    serial_files = serial_files[:args.limit]

# initial‐phase steps
INIT_STEPS = 20

# total frames = init(rotate+trans+zoom) + left/right + serial
total_frames  = INIT_STEPS*2 + max(len(left_files), len(right_files)) + len(serial_files)
frame_digits  = len(str(total_frames))
ffmpeg_pattern = f"frame%0{frame_digits}d.png"

# ──────────────── OUTPUT DIRECTORY ────────────────
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)
else:
    for f in glob.glob(os.path.join(args.output_dir, "*.png")):
        os.remove(f)

# ──────────────── PYMOL SETUP ────────────────
cmd.load(main_pdb, "main")
cmd.hide("everything", "all")
cmd.show("cartoon", "main")
cmd.bg_color("black")

viewwidth, viewheight = 1920, 1080
cmd.viewport(viewwidth, viewheight)
cmd.orient("main")
current_view = cmd.get_view()

cmd.set("ray_trace_frames",       1)
cmd.set("ray_shadow",             "on")
cmd.set("ambient_occlusion_mode", 1)
cmd.set("antialias",              2)
cmd.set("light_count",            2)

# ──────────────── FRAME SAVING ────────────────
frame_idx = 0
def save_frame():
    global frame_idx
    path = os.path.join(args.output_dir, f"frame{frame_idx:0{frame_digits}d}.png")
    cmd.png(path, width=viewwidth, height=viewheight, ray=1)
    frame_idx += 1

# ──────────────── 1) 180° ROTATION (INIT_STEPS) ────────────────
cmd.set_view(current_view)
for i in range(INIT_STEPS):
    cmd.turn("x", 180.0 / INIT_STEPS)
    save_frame()
rotated_view = cmd.get_view()

# ──────────────── 2) LOAD FIRST FRAGMENTS & CENTROIDS ────────────────
cmd.load(left_files[0],  "tmp_left")
cmd.load(right_files[0], "tmp_right")
cmd.hide("everything", "tmp_left or tmp_right")
cmd.show("cartoon", "tmp_left and name N+CA+C+O")
cmd.show("cartoon", "tmp_right and name N+CA+C+O")
cmd.select("focus_frag", "tmp_left or tmp_right")

# compute main centroid
main_coords    = [a.coord for a in cmd.get_model("main and polymer").atom]
main_centroid  = [sum(c[i] for c in main_coords)/len(main_coords) for i in range(3)]
# compute fragment centroid
frag_coords    = [a.coord for a in cmd.get_model("focus_frag").atom]
frag_centroid  = [sum(c[i] for c in frag_coords)/len(frag_coords) for i in range(3)]
vec_shift      = [frag_centroid[i] - main_centroid[i] for i in range(3)]

# ─────────────── INIT_STEPS ───────────────
#INIT_STEPS = 10

# 1) start from the 180°-rotated view
start_view = list(rotated_view)

# 2) compute the end view by doing one zoom+center on the fragment
cmd.set_view(start_view)
cmd.zoom("focus_frag", complete=1)
end_view = list(cmd.get_view())

# 3) one loop that both shifts pivot (translation) and adjusts zoom
for i in range(1, INIT_STEPS+1):
    t = i / float(INIT_STEPS)
    # blend every entry of the view matrix
    v = [ start_view[j] + (end_view[j] - start_view[j]) * t
          for j in range(len(start_view)) ]
    cmd.set_view(v)
    save_frame()

# 4) capture final view and clean up
current_view = tuple(end_view)
cmd.delete("tmp_left")
cmd.delete("tmp_right")
cmd.delete("focus_frag")


# ──────────────── RANDOMIZED LOOP MATCHING ────────────────
max_len = max(len(left_files), len(right_files))
if len(left_files) < max_len:
    random.shuffle(left_files)
    left_pool = list(islice(cycle(left_files), max_len))
else:
    left_pool = left_files

if len(right_files) < max_len:
    random.shuffle(right_files)
    right_pool = list(islice(cycle(right_files), max_len))
else:
    right_pool = right_files

# ──────────────── PHASE 1: LEFT + RIGHT ────────────────
for left, right in zip(left_pool, right_pool):
    cmd.delete("left_frag")
    cmd.delete("right_frag")
    cmd.load(left,  "left_frag")
    cmd.load(right, "right_frag")

    cmd.hide("everything", "left_frag or right_frag")
    cmd.show("cartoon", "left_frag and name N+CA+C+O")
    cmd.show("cartoon", "right_frag and name N+CA+C+O")
    cmd.color("cyan",   "left_frag")
    cmd.color("orange", "right_frag")

    cmd.set_view(current_view)
    save_frame()

# ──────────────── PHASE 2: SERIAL ────────────────
cmd.delete("left_frag")
cmd.delete("right_frag")

for serial in serial_files:
    cmd.delete("serial_frag")
    cmd.load(serial, "serial_frag")
    cmd.hide("everything", "serial_frag")
    cmd.show("cartoon", "serial_frag and name N+CA+C+O")
    cmd.color("purple", "serial_frag")

    cmd.set_view(current_view)
    save_frame()

# ──────────────── FFMPEG RENDER ────────────────
print(f"Rendering {frame_idx} frames into {args.output} in background...")
log_file = os.path.join(args.output_dir, "ffmpeg.log")
with open(log_file, "w") as log, open(os.devnull, "r") as devnull:
    subprocess.Popen([
        "ffmpeg", "-y", "-framerate", "24",
        "-i", f"{args.output_dir}/{ffmpeg_pattern}",
        "-c:v", "libx264", "-pix_fmt", "yuv420p",
        "-crf", "18", "-preset", "slow",
        args.output
    ], stdin=devnull, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)

print("✅ Done.")
