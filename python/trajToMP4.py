import sys
from pymol import cmd
import os
import glob
import argparse
import math

from pymol.cgo import *
import math

def create_gap_spheres(frame, object_name="trajectory", distance_threshold=6.0, previous_spheres=[]):
    """
    Detects gaps in a protein by checking CA-CA distances greater than the specified threshold,
    and creates a magenta, transparent sphere centered on the gap with a radius based on the CA-CA distance.
    
    Parameters:
    frame (int): The trajectory frame to analyze.
    object_name (str): The name of the object to use (default is "trajectory").
    distance_threshold (float): The distance (in Å) above which to consider a gap (default is 6.0 Å).
    previous_spheres (list): List of previously created sphere names to delete.
    """
    # Ensure the correct frame is set
    cmd.frame(frame)

    # Delete spheres from the previous frame
    for sphere_name in previous_spheres:
        cmd.delete(sphere_name)
    
    # Get the list of all Cα (CA) atoms in the current frame
    atoms = cmd.get_model(f"{object_name} and name CA", frame).atom

    # Initialize the list of CGO commands for spheres
    spherelist = []

    # Iterate over the CA atoms to calculate distances between consecutive atoms
    prev_atom = None
    for atom in atoms:
        if prev_atom is not None:
            # Calculate the distance between the previous CA atom and the current CA atom
            dist = math.sqrt((atom.coord[0] - prev_atom.coord[0]) ** 2 +
                             (atom.coord[1] - prev_atom.coord[1]) ** 2 +
                             (atom.coord[2] - prev_atom.coord[2]) ** 2)

            # Check if the distance exceeds the gap threshold
            if dist > distance_threshold:
                # Calculate the midpoint of the gap
                midpoint = [(atom.coord[0] + prev_atom.coord[0]) / 2,
                            (atom.coord[1] + prev_atom.coord[1]) / 2,
                            (atom.coord[2] + prev_atom.coord[2]) / 2]

                # Calculate the radius to be 50% of the CA-CA distance
                radius = dist / 2

                # Add the sphere to the spherelist
                spherelist.append(
#                    [COLOR, 1.0, 0.0, 1.0,  # Magenta color
#                     SPHERE, midpoint[0], midpoint[1], midpoint[2], radius]
                     [CYLINDER, atom.coord[0], atom.coord[1], atom.coord[2], prev_atom.coord[0], prev_atom.coord[1], prev_atom.coord[2], radius, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0,]
                )

        # Update prev_atom for next iteration
        prev_atom = atom

    # Load the spheres into PyMOL
    current_spheres = []
    for i, sphere in enumerate(spherelist):
        sphere_name = f'gap_sphere_{frame}_{i}'
        cmd.load_cgo(sphere, sphere_name)

        # Set the transparency of the sphere (75% transparent)
        cmd.set("cgo_transparency", 0.75, sphere_name)

        # Add the sphere name to the current_spheres list
        current_spheres.append(sphere_name)

    # Return the current list of sphere names for the next frame
    return current_spheres

# Parse Command-Line Arguments
parser = argparse.ArgumentParser(description="Process a PDB trajectory and generate a movie with optional gap highlighting.")
parser.add_argument("trajectory_file", help="The trajectory file (PDB format).")
parser.add_argument("ligand_code", help="The 3-letter code of the ligand.")
parser.add_argument("--addgap", action="store_true", help="Optional flag to add gap highlighting.")
args = parser.parse_args()

# Output configuration
output_folder = "movie_frames"
output_movie = "trajectory_movie.mp4"

# Step 1: Clean Output Folder
if not os.path.exists(output_folder):
    os.mkdir(output_folder)
else:
    # Remove existing PNG files in the output folder
    for file in glob.glob(f"{output_folder}/*.png"):
        os.remove(file)

# Step 2: Load the PDB Trajectory
cmd.load(args.trajectory_file, "trajectory")

# Step 3: Set up representations
# Hide all and show cartoons for protein
cmd.hide("everything", "all")
cmd.show("cartoon", "polymer")  # Show protein as cartoon

# Remove water for clarity
cmd.remove("resn HOH")  # Removes water molecules

# Highlight ligand in CPK representation
ligand_selection = f"resn {args.ligand_code}"  # Use the input ligand code
cmd.show("spheres", ligand_selection)  # Show ligand as spheres

# Step 4: Setup DSSP for Secondary Structure Assignment
cmd.dss("trajectory")
cmd.color("red", "ss H")
cmd.color("yellow", "ss S")
cmd.color("green", "ss L")

# Ray tracing and rendering settings
cmd.set("ray_trace_frames", 1)
cmd.set("ray_shadow", "on")  # Enable ray tracing shadows
cmd.set("ambient_occlusion_mode", 1)  # Enable ambient occlusion for softer lighting
cmd.set("light_count", 1)  # Increase the number of light sources for smoother shadows
cmd.set("antialias", 2)

# this modification is for 1/2 a widescreen (for a side-by-side view)
viewwidth=3840/2
viewheight=2160

cmd.viewport(viewwidth, viewheight)  # 16:9 aspect ratio for 4K
cmd.clip("near", 10)     # Move the near clipping plane closer
cmd.clip("far", -10)     # Move the far clipping plane further away

cmd.select("zoom_selection", "polymer or organic and not resn HOH+NA+CL")
#cmd.orient("zoom_selection")  # Orient the selected atoms
cmd.zoom("zoom_selection")  # Zoom with a buffer around the selection

current_view = cmd.get_view()
cmd.png(f"{output_folder}/test.png", width=viewwidth, height=viewheight, ray=1)

# Step 5: Iterate through Frames and Save 4K Image
frame_count = cmd.count_states("trajectory")
previous_spheres = []
for frame in range(1, frame_count + 1):
    cmd.frame(frame)
    if args.addgap:
        previous_spheres = create_gap_spheres(frame, previous_spheres=previous_spheres)
    # Save the PNG frame at 4K resolution
    cmd.set_view(current_view)
    cmd.png(f"{output_folder}/frame{frame:04d}.png", width=viewwidth, height=viewheight, ray=1)
    
    # Break the loop after 10 frames
#    if frame >= 1:
#        break

# Step 6: Convert PNG Frames to 4K MP4 Using Two-Pass FFmpeg
# First pass
os.system(f"ffmpeg -y -framerate 24 -i {output_folder}/frame%04d.png -c:v libx264 -pix_fmt yuv420p -crf 18 -preset slow -b:v 20M -pass 1 -an -f mp4 /dev/null")

# Second pass
os.system(f"ffmpeg -framerate 24 -i {output_folder}/frame%04d.png -c:v libx264 -pix_fmt yuv420p -crf 18 -preset slow -b:v 20M -pass 2 {output_movie}")

# Optional: Print completion message
print(f"Movie saved as {output_movie}")

# Usage Note: Upon completion of the LEFT.mp4 and the RIGHT.mp4, use the following to generate the combined:
# % ffmpeg -i LEFT.mp4 -i RIGHT.mp4 -filter_complex "\
# [0:v]scale=1920:2160,drawtext=text='MD With Gap':x=(w-text_w)/2:y=300:fontsize=96:fontcolor=white[v0]; \
# [1:v]scale=1920:2160,drawtext=text='MD Without Gap':x=(w-text_w)/2:y=300:fontsize=96:fontcolor=white[v1]; \
# [v0][v1]hstack=inputs=2[stacked]; \
# [stacked]drawtext=text='Snapshot\\: %{eif\\:2*n+1\\:d} / 999':x=(w-text_w)/2:y=h-100:fontsize=48:fontcolor=white[v]" -map "[v]" output.mp4

