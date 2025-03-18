import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ks_2samp  # For Kolmogorov-Smirnov test
import sys
from Bio.Align import PairwiseAligner

# With the need to perform sequence analysis, DSSP, etc, we need to rename our residues to match the standards.
#       This renaming should account for not only the A, B, and C (etc) suffixes added to the residue names but equivelentces like HIE, HID, and HIP
residue_mapping = {
    # Standard amino acids
    'ALA': 'ALA', 'ARG': 'ARG', 'ASN': 'ASN', 'ASP': 'ASP', 'CYS': 'CYS',
    'GLU': 'GLU', 'GLN': 'GLN', 'GLY': 'GLY', 'HIS': 'HIS', 'ILE': 'ILE',
    'LEU': 'LEU', 'LYS': 'LYS', 'MET': 'MET', 'PHE': 'PHE', 'PRO': 'PRO',
    'SER': 'SER', 'THR': 'THR', 'TRP': 'TRP', 'TYR': 'TYR', 'VAL': 'VAL',

    # Amber protonation states
    'HIE': 'HIS', 'HID': 'HIS', 'HIP': 'HIS',  # Histidine
    'CYX': 'CYS',  # Cysteine
    'GLH': 'GLU',  # Glutamate
    'ASH': 'ASP',  # Aspartate
}

def standardize_residue_names(trajectory):
    # Create a new list for the standardized residue names
    standardized_residues = []
    
    for residue in trajectory.topology.residues:
        residue_name = residue.name
        # Only remove the last character if the length is greater than 3
        if len(residue_name) > 3:
            trimmed_residue_name = residue_name[:-1]  # Remove the last character
        else:
            trimmed_residue_name = residue_name
        
        # Map the trimmed residue name to the standard name
        standardized_name = residue_mapping.get(trimmed_residue_name, None)

        # Print if there is no match in the mapping (excluding solvent)
        if standardized_name is None and trimmed_residue_name not in ["HOH", "NA", "CL"]:
            print(f"Warning: Unmatched residue name '{trimmed_residue_name}' in trajectory.")
        
        standardized_residues.append(standardized_name if standardized_name is not None else trimmed_residue_name)
    
    return standardized_residues

trajFile1 = sys.argv[1]
topoFile1 = sys.argv[2]
trajFile2 = sys.argv[3]
topoFile2 = sys.argv[4]

trajLabel1 = trajFile1
trajLabel2 = trajFile2

# Load two trajectories
input_trajectory_1 = md.load(trajFile1,top=topoFile1)

# Print some basic information about the loaded trajectory
print("Number of frames (1):", input_trajectory_1.n_frames)
print("Number of atoms (1):", input_trajectory_1.n_atoms)

input_trajectory_2 = md.load(trajFile2,top=topoFile2)
print("Number of frames (2):", input_trajectory_2.n_frames)
print("Number of atoms (2):", input_trajectory_2.n_atoms)

# Standardize residue names
standardized_residues_1 = standardize_residue_names(input_trajectory_1)
standardized_residues_2 = standardize_residue_names(input_trajectory_2)

# Create an aligner
aligner = PairwiseAligner()

# Define a set of valid standard residues from the residue_mapping keys
valid_residues = set(residue_mapping.keys())

# Filter to keep only standardized amino acid residues for alignment
filtered_residues_1 = [res for res in standardized_residues_1 if res in valid_residues]
filtered_residues_2 = [res for res in standardized_residues_2 if res in valid_residues]

# Convert lists of filtered residues to strings for alignment
sequence_1 = "".join(filtered_residues_1)
sequence_2 = "".join(filtered_residues_2)

# Perform the alignment with the filtered sequences
alignment = aligner.align(sequence_1, sequence_2)

# Assuming we want the first alignment result
aligned_seq_1 = alignment[0].target  # The aligned sequence from sequence_1
aligned_seq_2 = alignment[0].query   # The aligned sequence from sequence_2

# Convert the aligned sequences to strings if they are still in array form
aligned_seq_1_str = ''.join([str(res) for res in aligned_seq_1])
aligned_seq_2_str = ''.join([str(res) for res in aligned_seq_2])

# Print the aligned sequences for verification
#print("Aligned Sequence 1:", aligned_seq_1_str)
#print("Aligned Sequence 2:", aligned_seq_2_str)

# Replace the residue names in the topology
for residue, new_name in zip(input_trajectory_1.topology.residues, standardized_residues_1):
    residue.name = new_name

for residue, new_name in zip(input_trajectory_2.topology.residues, standardized_residues_2):
    residue.name = new_name

# Select protein atoms (exclude water)
protein_atoms_1 = input_trajectory_1.topology.select("not (resname HOH CL NA)")
protein_atoms_2 = input_trajectory_2.topology.select("not (resname HOH CL NA)")

# Slice the trajectories to include only protein atoms
trajectory_1 = input_trajectory_1.atom_slice(protein_atoms_1)
trajectory_2 = input_trajectory_2.atom_slice(protein_atoms_2)

# Align based on the CA carbons
ca_indices_1 = trajectory_1.topology.select("protein and name CA")
ca_indices_2 = trajectory_2.topology.select("protein and name CA")
trajectory_2 = trajectory_2.superpose(trajectory_1, atom_indices=ca_indices_1)

# Create masks for non-gap residues
mask_1 = np.isin(trajectory_1.topology.residues, trajectory_1.topology.select("not (resname HOH CL NA)"))
mask_2 = np.isin(trajectory_2.topology.residues, trajectory_2.topology.select("not (resname HOH CL NA)"))

# ---------------------------------------
# Compute RMSF for both trajectories
# ---------------------------------------

# Compute RMSF for alpha carbons only (for both trajectories)
rmsf_1 = np.sqrt(np.mean((trajectory_1.xyz[:, ca_indices_1, :] - np.mean(trajectory_1.xyz[:, ca_indices_1, :], axis=0))**2, axis=0))
rmsf_2 = np.sqrt(np.mean((trajectory_2.xyz[:, ca_indices_2, :] - np.mean(trajectory_2.xyz[:, ca_indices_2, :], axis=0))**2, axis=0))

# Average over dimensions (x, y, z)
rmsf_per_residue_1 = np.sqrt(np.sum(rmsf_1**2, axis=1))
rmsf_per_residue_2 = np.sqrt(np.sum(rmsf_2**2, axis=1))

# Calculate average and standard deviation
avg_rmsf_1 = np.mean(rmsf_per_residue_1)
std_rmsf_1 = np.std(rmsf_per_residue_1)
avg_rmsf_2 = np.mean(rmsf_per_residue_2)
std_rmsf_2 = np.std(rmsf_per_residue_2)

# Kolmogorov-Smirnov test on RMSF values
ks_stat_rmsf, p_value_rmsf = ks_2samp(rmsf_per_residue_1, rmsf_per_residue_2)

# Report results
print(f"Trajectory 1: Average RMSF = {avg_rmsf_1:.3f}, Standard Deviation = {std_rmsf_1:.3f}")
print(f"Trajectory 2: Average RMSF = {avg_rmsf_2:.3f}, Standard Deviation = {std_rmsf_2:.3f}")
print(f"K-S test for RMSF: D = {ks_stat_rmsf}, p-value = {p_value_rmsf}")

# Plot RMSF for both trajectories
plt.figure()
plt.plot(ca_indices_1, rmsf_per_residue_1, label=trajLabel1, color='b')
plt.plot(ca_indices_2, rmsf_per_residue_2, label=trajLabel2, color='r')
plt.xlabel('Residue Index (Alpha Carbon)')
plt.ylabel('RMSF (nm)')
plt.title(f'Root Mean Square Fluctuation per Residue\nK-S test p-value: {p_value_rmsf:.4f}')
plt.legend()
plt.savefig('rmsf_comparison.pdf', dpi=300, bbox_inches='tight')

# ---------------------------------------
# Compute Radius of Gyration (Rg) for both trajectories
# ---------------------------------------
rg_1 = md.compute_rg(trajectory_1)
rg_2 = md.compute_rg(trajectory_2)

# Calculate average and standard deviation of Rg
avg_rg_1 = np.mean(rg_1)
std_rg_1 = np.std(rg_1)
avg_rg_2 = np.mean(rg_2)
std_rg_2 = np.std(rg_2)

# Kolmogorov-Smirnov test on Rg values
ks_stat_rg, p_value_rg = ks_2samp(rg_1, rg_2)

# Report results
print(f"Trajectory 1: Average Rg = {avg_rg_1:.3f}, Standard Deviation = {std_rg_1:.3f}")
print(f"Trajectory 2: Average Rg = {avg_rg_2:.3f}, Standard Deviation = {std_rg_2:.3f}")
print(f"K-S test for Rg: D = {ks_stat_rg:.3f}, p-value = {p_value_rg:.3e}")

plt.figure()
plt.plot(rg_1, label=trajLabel1, color='b')
plt.plot(rg_2, label=trajLabel2, color='r')
plt.xlabel('Frame (time)')
plt.ylabel('Radius of Gyration (nm)')
plt.title(f'Radius of Gyration Over Time\nK-S test p-value: {p_value_rg:.4f}')
plt.legend()
plt.savefig('rg_comparison.pdf', dpi=300, bbox_inches='tight')

sys.exit()

# ---------------------------------------
# Compute Secondary Structure (DSSP) for both trajectories
# ---------------------------------------
ss_1 = md.compute_dssp(trajectory_1)
ss_2 = md.compute_dssp(trajectory_2)

# DSSP structure mapping
structure_mapping = {'H': 0, 'E': 1, 'C': 2, 'G': 3, 'I': 4, 'B': 5, 'T': 6, 'S': 7}
structure_labels = ['Helix (H)', 'Beta sheet (E)', 'Coil (C)', '3-10 Helix (G)', 'Pi Helix (I)', 
                    'Beta Bridge (B)', 'Turn (T)', 'Bend (S)']

# Convert DSSP characters to numerical values
ss_numeric_1 = np.array([[structure_mapping[ss_elem] for ss_elem in frame] for frame in ss_1])
ss_numeric_2 = np.array([[structure_mapping[ss_elem] for ss_elem in frame] for frame in ss_2])

# Create heatmap for Trajectory 1
plt.figure(figsize=(12, 6))
plt.subplot(2, 1, 1)
sns.heatmap(ss_numeric_1.T, cmap='tab20', cbar_kws={'ticks': list(structure_mapping.values()), 
                                                    'label': 'Secondary Structure'},
            yticklabels=False)
plt.xlabel('Frame (time)')
plt.ylabel('Residue')
plt.title(f"{trajLabel1}: Secondary Structure Evolution")

# Create heatmap for Trajectory 2
plt.subplot(2, 1, 2)
sns.heatmap(ss_numeric_2.T, cmap='tab20', cbar_kws={'ticks': list(structure_mapping.values()), 
                                                    'label': 'Secondary Structure'},
            yticklabels=False)
plt.xlabel('Frame (time)')
plt.ylabel('Residue')
plt.title(f"{trajLabel2}: Secondary Structure Evolution")

# Customize colorbar to show structure labels instead of numbers
cbar = plt.gca().collections[0].colorbar
cbar.set_ticks(list(structure_mapping.values()))
cbar.set_ticklabels(structure_labels)

# Save the DSSP heatmap comparison
plt.savefig('dssp_comparison.pdf', dpi=300, bbox_inches='tight')

