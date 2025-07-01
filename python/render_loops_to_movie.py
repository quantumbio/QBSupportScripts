import os
import sys
import glob
import time
import random
import argparse
import subprocess
from itertools import cycle, islice
from pymol import cmd

random.seed()

# ──────────────── ARGUMENTS ────────────────
parser = argparse.ArgumentParser(description="Render loop ensembles (left/right + serial) into PyMOL movie.")

parser.add_argument("--main_pdb", required=True, help="Path to the main protein structure (e.g., main.pdb)")
parser.add_argument("--left_glob", required=True, help="Glob pattern for left fragments (e.g., 'loops/left/**/*.pdb')")
parser.add_argument("--right_glob", required=True, help="Glob pattern for right fragments (e.g., 'loops/right/**/*.pdb')")
parser.add_argument("--serial_glob", required=True, help="Glob pattern for serial fragments (e.g., 'loops/serial/**/*.pdb')")

parser.add_argument("--output", default="loop_ensemble_movie.mp4", help="Output MP4 file name")
parser.add_argument("--output_dir", default="movie_frames", help="Directory to save PNG frames")
parser.add_argument("--limit", type=int, default=None,
                    help="If set, only use the first N fragments from each glob for quick testing")

args = parser.parse_args()

# ──────────────── FILE DISCOVERY ────────────────
main_pdb = args.main_pdb
left_files = sorted(glob.glob(args.left_glob, recursive=True))
right_files = sorted(glob.glob(args.right_glob, recursive=True))
serial_files = sorted(glob.glob(args.serial_glob, recursive=True))

if not os.path.exists(main_pdb):
    sys.exit("❌ main.pdb not found")

if not left_files or not right_files or not serial_files:
    sys.exit("❌ One or more fragment sets not found")

# Optional limit for testing
if args.limit:
    left_files = left_files[:2*args.limit]
    right_files = right_files[:2*args.limit]
    serial_files = serial_files[:args.limit]

total_frames = max(len(left_files), len(right_files)) + len(serial_files)
frame_digits = len(str(total_frames))
ffmpeg_pattern = f"frame%0{frame_digits}d.png"

# ──────────────── OUTPUT DIRECTORY ────────────────
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)
else:
    for f in glob.glob(os.path.join(args.output_dir, "*.png")):
        os.remove(f)

# ──────────────── PyMOL CONFIGURATION ────────────────
cmd.load(main_pdb, "main")
cmd.hide("everything", "all")
cmd.show("cartoon", "main")
cmd.bg_color("black")

viewwidth, viewheight = 1920, 1080
cmd.viewport(viewwidth, viewheight)
cmd.orient("main")
current_view = cmd.get_view()

cmd.set("ray_trace_frames", 1)
cmd.set("ray_shadow", "on")
cmd.set("ambient_occlusion_mode", 1)
cmd.set("antialias", 2)
cmd.set("light_count", 2)

# ──────────────── RANDOMIZED LOOP MATCHING ────────────────
frame_idx = 0
max_len = max(len(left_files), len(right_files))

# Proper random cycling of the smaller pool
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
did_initial_orient = False
for left, right in zip(left_pool, right_pool):
    cmd.delete("left_frag")
    cmd.delete("right_frag")
    cmd.load(left, "left_frag")
    cmd.load(right, "right_frag")

    cmd.hide("everything", "left_frag")
    cmd.hide("everything", "right_frag")
    cmd.show("cartoon", "left_frag and name N+CA+C+O")
    cmd.show("cartoon", "right_frag and name N+CA+C+O")
    cmd.color("cyan", "left_frag")
    cmd.color("orange", "right_frag")

#    cmd.set_view(current_view)

    # Only orient on first pair
    if not did_initial_orient:
        cmd.orient("left_frag or right_frag")
        cmd.turn("x", 180)
        current_view = cmd.get_view()  # Capture this new view to reuse
        did_initial_orient = True
    else:
        cmd.set_view(current_view)

    out_path = os.path.join(args.output_dir, f"frame{frame_idx:0{frame_digits}d}.png")
    cmd.png(out_path, width=viewwidth, height=viewheight, ray=1)
    frame_idx += 1

# ──────────────── PHASE 2: SERIAL ────────────────
# Finished all left/right pairs
cmd.delete("left_frag")
cmd.delete("right_frag")

for serial in serial_files:
    cmd.delete("serial_frag")
    cmd.load(serial, "serial_frag")
    cmd.color("purple", "serial_frag")

    cmd.set_view(current_view)
    out_path = os.path.join(args.output_dir, f"frame{frame_idx:0{frame_digits}d}.png")
    cmd.png(out_path, width=viewwidth, height=viewheight, ray=1)
    frame_idx += 1

# ──────────────── FFMPEG RENDER ────────────────
print(f"Rendering {frame_idx} frames into {args.output} in background...")

log_file = os.path.join(args.output_dir, "ffmpeg.log")
with open(log_file, "w") as log, open(os.devnull, "r") as devnull:
    subprocess.Popen(
        [
            "ffmpeg",
            "-y",
            "-framerate", "24",
            "-i", f"{args.output_dir}/{ffmpeg_pattern}",
            "-c:v", "libx264",
            "-pix_fmt", "yuv420p",
            "-crf", "18",
            "-preset", "slow",
            args.output
        ],
        stdin=devnull,
        stdout=log,
        stderr=subprocess.STDOUT,
        start_new_session=True  # Optional: detach from terminal session
    )

print(f"Background FFmpeg started. Log: {log_file}")
print("✅ Done.")
