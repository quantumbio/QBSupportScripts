import argparse
import gzip
import mdtraj as md
import os

# Set up argument parser
parser = argparse.ArgumentParser(description="Process a trajectory file and save as PDB.")
parser.add_argument("topology_file", type=str, help="Path to the topology file (e.g., PDB).")
parser.add_argument("trajectory_file", type=str, help="Path to the trajectory file (e.g., DCD).")
parser.add_argument("output_pdb_file", type=str, help="Path to save the output PDB file.")
parser.add_argument("--frame", type=int, default=1, help="Write every Nth frame (default: 1).")

# Parse arguments
args = parser.parse_args()

# Validate --frame argument
if args.frame < 1:
    raise ValueError("--frame must be an integer greater than or equal to 1.")

# Load the trajectory and topology
try:
    trajectory = md.load(args.trajectory_file, top=args.topology_file)
except Exception as e:
    raise ValueError(f"Error loading trajectory: {e}")

# Print trajectory information
print("Trajectory Report:")
print(f"Number of frames: {trajectory.n_frames}")
print(f"Number of atoms: {trajectory.n_atoms}")
print(f"Number of residues: {trajectory.topology.n_residues}")
print(f"Number of chains: {trajectory.topology.n_chains}")
print(f"Time between frames (ps): {trajectory.timestep:.2f}" if trajectory.timestep else "Timestep unavailable")
print(f"Unit cell dimensions (if present): {trajectory.unitcell_lengths}")

# Slice the trajectory to include every Nth frame
sliced_trajectory = trajectory[::args.frame]

print(f"UPDATED: Number of frames: {sliced_trajectory.n_frames}")
print(f"UPDATED: Time between frames (ps): {sliced_trajectory.timestep:.2f}" if sliced_trajectory.timestep else "Timestep unavailable")

# Image molecules into the primary box and make them whole
sliced_trajectory.image_molecules(inplace=True, make_whole=True)

# Write every Nth frame to a regular PDB
print(f"Writing to {args.output_pdb_file}...")
sliced_trajectory.save_pdb(args.output_pdb_file)
print(f"PDB saved as {args.output_pdb_file}")

# Compress the PDB file manually
with open(args.output_pdb_file, 'rb') as f_in:
    with gzip.open(args.output_pdb_file + '.gz', 'wb') as f_out:
        f_out.writelines(f_in)
    print(f"Compressed PDB saved as {args.output_pdb_file}.gz")

# Delete the original PDB file
os.remove(args.output_pdb_file)
