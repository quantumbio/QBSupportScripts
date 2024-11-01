import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ks_2samp  # For Kolmogorov-Smirnov test
import sys

trajFile1 = sys.argv[1]
topoFile1 = sys.argv[2]
trajFile2 = sys.argv[3]
topoFile2 = sys.argv[4]

trajLabel1 = trajFile1
trajLabel2 = trajFile2

# Load two trajectories
input_trajectory_1 = md.load(trajFile1,top=topoFile1)
input_trajectory_2 = md.load(trajFile2,top=topoFile2)

# Select protein atoms (exclude water)
protein_atoms_1 = input_trajectory_1.topology.select("protein")
protein_atoms_2 = input_trajectory_2.topology.select("protein")

# Slice the trajectories to include only protein atoms
trajectory_1 = input_trajectory_1.atom_slice(protein_atoms_1)
trajectory_2 = input_trajectory_2.atom_slice(protein_atoms_2)

# ---------------------------------------
# Compute RMSF for both trajectories
# ---------------------------------------
ca_indices_1 = trajectory_1.topology.select("protein and name CA")
ca_indices_2 = trajectory_2.topology.select("protein and name CA")

# Compute RMSF for alpha carbons only (for both trajectories)
rmsf_1 = np.sqrt(np.mean((trajectory_1.xyz[:, ca_indices_1, :] - np.mean(trajectory_1.xyz[:, ca_indices_1, :], axis=0))**2, axis=0))
rmsf_2 = np.sqrt(np.mean((trajectory_2.xyz[:, ca_indices_2, :] - np.mean(trajectory_2.xyz[:, ca_indices_2, :], axis=0))**2, axis=0))

# Average over dimensions (x, y, z)
rmsf_per_residue_1 = np.sqrt(np.sum(rmsf_1**2, axis=1))
rmsf_per_residue_2 = np.sqrt(np.sum(rmsf_2**2, axis=1))

# Kolmogorov-Smirnov test on RMSF values
ks_stat_rmsf, p_value_rmsf = ks_2samp(rmsf_per_residue_1, rmsf_per_residue_2)
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

# Kolmogorov-Smirnov test on Rg values
ks_stat_rg, p_value_rg = ks_2samp(rg_1, rg_2)
print(f"K-S test for Radius of Gyration: D = {ks_stat_rg}, p-value = {p_value_rg}")

plt.figure()
plt.plot(rg_1, label=trajLabel1, color='b')
plt.plot(rg_2, label=trajLabel2, color='r')
plt.xlabel('Frame (time)')
plt.ylabel('Radius of Gyration (nm)')
plt.title(f'Radius of Gyration Over Time\nK-S test p-value: {p_value_rg:.4f}')
plt.legend()
plt.savefig('rg_comparison.pdf', dpi=300, bbox_inches='tight')

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

