#  /***********************************************************************
#     Copyright (c) 2024 QuantumBio Inc. and/or its affiliates.
# 	
#  This source code is the property of QuantumBio Inc. and/or its affiliates
#  and is provided AS IS.
# 
#  This source code may contain proprietary and Confidential Information, 
#  including trade secrets, belonging to QuantumBio Inc. and/or its 
#  affiliates.
# 
#  Please see http://www.quantumbioinc.com/ for more information.
# 
#  ***********************************************************************/
# 
# The script provided below will do the following:
# 
#     * Read in a supplied *.prmtop *.inpcrd file pair.
#     * Generate an XML file including all required OpenMM parameters for the AMBER potential encoded in the *.prmtop file.
#     * Read in this XML file along with the *.inpcrd file and generate a standard waterbox + counterions for the structure.
#     * Perform a minimization on the new system.
#     * Perform a NVT+NPT equilibration.
#     * Perform a 2 ns MD simulation.
# 
# Use of this script requires the installation OpenMM, MDTraj, and ParmED.
# You may (optionally) use conda to install this environment:
# 
#     % conda create --name openmm
#     % conda activate openmm
#     % conda install -c conda-forge mdtraj openmm pdbfixer openmm-setup parmed
# 
# To generate the input used for this script, using qmechanic, perform the following step:
# 
#     % /path/to/DivConSuite/bin/qmechanic 1lri --prepare all -v 2 -h amberff14sb -O -p prmtop inpcrd
#     This will perform the following steps:
#         * Complete the addition of any missing residues / R-groups / etc using the included SEQRES in the 1lri PDB file
#         * Complete the addition of any missing protons
#         * Determine GAFF types and AM1-BCC charges for any unknown residues/ligands
#         * Crystallographically refine any added atoms/residues along with any significant bond/torsion outliers using amberff14sb.
#         * Write a final pair of AMBER-compatible prmtop and inpcrd
# 
# These two files are then provided to this runMD-fromDivCon.py script within an OpenMM environment:
# 
#     % conda activate openmm
#     % python runMD-fromDivCon.py 1lri_out.prmtop 1lri_out.inpcrd
#

import openmm as mm
import openmm.app as app
import simtk.unit as unit
import sys
import os
import numpy as np
import parmed as pmd
import xml.etree.ElementTree as ET
import hashlib
import string
import importlib.resources
import mdtraj as md
import time

minimize_nsteps     = round (5000      /   10)
nvt_equil_nsteps    = round (100000    /   10)
ntp_equil_nsteps    = round (100000    /   10)
production_nsteps   = round (1000000   /   100)

# Load the Amber topology and coordinate files
prmtopFile = sys.argv[1]
inpcrdFile = sys.argv[2]
prmtop = pmd.load_file(prmtopFile, inpcrdFile)

# Identify and modify HOH residues
for residue in prmtop.residues:
    if residue.name == "HOH":
        # Identify H1 and H2
        h1 = residue.atoms[0]  # Adjust index if necessary
        h2 = residue.atoms[1]  # Adjust index if necessary
        
        # Find and remove the bond between H1 and H2
        for bond in h1.bonds:
            if bond.atom2 == h2 or bond.atom1 == h2:
                # Remove the bond from both atoms
                h1.bonds.remove(bond)
                h2.bonds.remove(bond)
                break  # Exit loop after removing the bond

prmtop.save('initial.pdb')

# Initialize for unique residue names
alphabet = string.ascii_uppercase  # 'A', 'B', 'C', ..., 'Z'
unique_bond_fingerprints = {}
new_residue_names = {}

# Function to generate a unique hash of the bonding pattern of a residue
def generate_bond_fingerprint(residue):
    bond_list = []
    for atom in residue.atoms:
        for bond in atom.bonds:
            bond_list.append((bond.atom1.name, bond.atom2.name))
    
    # Sort bonds and create a hash of the bonding pattern
    bond_list.sort()
    bond_fingerprint = hashlib.md5(str(bond_list).encode('utf-8')).hexdigest()
    return bond_fingerprint

# Identify all unique residue types in the prmtop file
unique_residue_types = set(res.name for res in prmtop.residues)

# Initialize counters for unique names per residue type
residue_counters = {res_type: 0 for res_type in unique_residue_types}

# Initialize a mapping of unique names back to original names
original_residue_names = {}

# Initialize counter
counter = 0

# Loop through all residues and analyze their bonding patterns for all residue types
for residue in prmtop.residues:
    residue_type = residue.name

    # Generate the bond fingerprint for this residue
    bond_fingerprint = generate_bond_fingerprint(residue)
    
    # Retrieve chain ID from the first atom in the residue
    chain_id = residue.atoms[0].tree if residue.atoms else 'Unknown'
    if chain_id == "BLA":
        # Convert counter to a letter
        code = chr(ord('A') + (counter % 26))
        counter += 1

    # Check if we have seen this bonding pattern before for the specific residue type
    if (residue_type, bond_fingerprint) not in unique_bond_fingerprints or residue.name == "PO4":
        # Assign a unique identifier (A, B, C, etc.) per residue type
        unique_residue_code = alphabet[residue_counters[residue_type] % len(alphabet)]
        unique_residue_name = residue_type + unique_residue_code

        # Print Chain, Residue Name, and Residue UID information
        print(f"New unique residue {unique_residue_name}: Chain {code}, Name {residue.name}, UID {residue.idx}")
        
        # Save the unique fingerprint and name
        unique_bond_fingerprints[(residue_type, bond_fingerprint)] = unique_residue_name
        new_residue_names[residue] = unique_residue_name
                
        # Increment the counter for the residue type
        residue_counters[residue_type] += 1
    else:
        # This bonding pattern has been seen before for this residue type, assign the same name
        new_residue_names[residue] = unique_bond_fingerprints[(residue_type, bond_fingerprint)]

# Extract defined residues from the forcefield using XML parsing
defined_residues = set()
forcefield_xml_file = importlib.resources.files('openmm.app.data') / 'amber14/tip3pfb.xml'

# Parse the XML to find all residue names
tree = ET.parse(forcefield_xml_file)
root = tree.getroot()

# Look for <Residue> elements in the XML
for residue in root.findall(".//Residue"):
    defined_residues.add(residue.get('name'))

# Initialize skip_residues with the defined residues from the forcefield
skip_residues = {"WAT", "HOH"}.union(defined_residues)

# Rename residues in the topology based on unique bonding patterns
for residue in prmtop.residues:
    residue_name = residue.name
    
    # Skip renaming for residues in the skip_residues set
    if residue_name in skip_residues:
        continue

    if residue in new_residue_names:
        # Store the original name to revert later
        original_residue_names[new_residue_names[residue]] = residue.name
        residue.name = new_residue_names[residue]

# Save the modified topology to check
prmtop.save('modified_with_unique_residues.prmtop')
prmtop.save('modified_with_unique_residues.pdb')
prmtop.save('modified_with_unique_residues.inpcrd')

# Now create the XML force field with the unique residues
param_set = pmd.openmm.OpenMMParameterSet.from_structure(prmtop)
param_set.write('complete_forcefield_with_unique_residues.xml')

# Load the force field XML file
tree = ET.parse('complete_forcefield_with_unique_residues.xml')
root = tree.getroot()

# Create a set to track unique residue names already added to XML
written_residues = set()

# Create a new <Residues> section
residues_element = ET.Element("Residues")

# Loop over the unique residues and their parameters
for residue in prmtop.residues:
    residue_name = residue.name
    
    # Skip writing residue information for residues in the skip_residues set
    if residue_name in skip_residues:
        continue  # Skip this residue if it's in the skip set
    
    # Check if this residue has already been added to the XML
    if residue_name in written_residues:
        continue  # Skip this residue if it has already been processed

    residue_element = ET.Element("Residue", name=residue_name)

    # Add atoms
    for atom in residue.atoms:
        atom_element = ET.SubElement(residue_element, "Atom", name=atom.name, type=atom.type, charge=str(atom.charge))
    
    # Add bonds and external bonds
    bonded_atoms = set()
    for atom in residue.atoms:
        for bond in atom.bonds:
            other_atom = bond.atom2 if bond.atom1 == atom else bond.atom1
            
            if other_atom.residue == atom.residue:
                # Internal bond: atoms belong to the same residue
                bond_tuple = tuple(sorted([bond.atom1.name, bond.atom2.name]))
                if bond_tuple not in bonded_atoms:
                    bond_element = ET.SubElement(residue_element, "Bond", atomName1=bond.atom1.name, atomName2=bond.atom2.name)
                    bonded_atoms.add(bond_tuple)
            else:
                # External bond: atoms belong to different residues
                external_bond_element = ET.SubElement(residue_element, "ExternalBond", atomName=atom.name)

    # Append residue to the <Residues> section
    residues_element.append(residue_element)
    
    # Add the residue name to the written set
    written_residues.add(residue_name)

# Append <Residues> section to the root of the force field XML
root.append(residues_element)

# Function to indent XML for pretty printing
def indent(elem, level=0):
    i = "\n" + level * "  "  # Use two spaces for indentation
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for subelem in elem:
            indent(subelem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

# Indent the root for pretty printing
indent(root)

# Write the modified and indented XML to a file
tree.write('complete_forcefield_with_unique_residues.xml', xml_declaration=True, encoding='utf-8', method="xml")

# Load the modified prmtop and forcefield XML for simulation
prmtop = app.AmberPrmtopFile('modified_with_unique_residues.prmtop')
inpcrd = app.AmberInpcrdFile('modified_with_unique_residues.inpcrd')

# Extract positions from inpcrd
positions = inpcrd.getPositions()

# Convert the positions to a NumPy array for easy manipulation
positions_array = np.array([[pos[0].value_in_unit(unit.nanometers),
                              pos[1].value_in_unit(unit.nanometers),
                              pos[2].value_in_unit(unit.nanometers)] for pos in positions])

# Calculate bounding box dimensions based on XYZ coordinates
x_coords = positions_array[:, 0]
y_coords = positions_array[:, 1]
z_coords = positions_array[:, 2]

min_x, max_x = min(x_coords), max(x_coords)
min_y, max_y = min(y_coords), max(y_coords)
min_z, max_z = min(z_coords), max(z_coords)

# Create box dimensions with padding
padding = 10.0 * unit.nanometers  # Example padding
box_size_x = (max_x - min_x) + 2 * padding.value_in_unit(unit.nanometers)
box_size_y = (max_y - min_y) + 2 * padding.value_in_unit(unit.nanometers)
box_size_z = (max_z - min_z) + 2 * padding.value_in_unit(unit.nanometers)

# Create box vectors
box_vectors = [
    mm.Vec3(box_size_x, 0.0, 0.0),
    mm.Vec3(0.0, box_size_y, 0.0),
    mm.Vec3(0.0, 0.0, box_size_z),
]

# Set periodic box vectors for the topology
prmtop.topology.setPeriodicBoxVectors(box_vectors)

# Get the periodic box vectors
box_vectors = prmtop.topology.getPeriodicBoxVectors()

# Check if the box vectors are defined
if box_vectors is not None:
    # Calculate the box dimensions
    box_dimensions = [
        box_vectors[0][0].value_in_unit(unit.nanometers),
        box_vectors[1][1].value_in_unit(unit.nanometers),
        box_vectors[2][2].value_in_unit(unit.nanometers)
    ]

    print("Box Dimensions (nm):", box_dimensions)
else:
    print("No periodic box dimensions are defined in the topology.")

forcefield = app.ForceField('amber14/tip3pfb.xml','complete_forcefield_with_unique_residues.xml')

# Use Modeller to create a waterbox
modeller = app.Modeller(prmtop.topology, inpcrd.positions)

# Create box dimensions with padding
padding = 1.0 * unit.nanometers  # Example padding
box_size_x = (max_x - min_x) + 2 * padding.value_in_unit(unit.nanometers)
box_size_y = (max_y - min_y) + 2 * padding.value_in_unit(unit.nanometers)
box_size_z = (max_z - min_z) + 2 * padding.value_in_unit(unit.nanometers)

# Create box vector for solvent addition
box_vector = mm.Vec3(box_size_x, box_size_y, box_size_z)

# Add a water box and neutralize the system with counterions
modeller.addSolvent(forcefield, model='tip3p', boxSize=box_vector, ionicStrength=0.15*unit.molar)

# Now create the system again after modifying the topology
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometers,
    constraints=app.HBonds
)

# Set up the integrator
integrator = mm.LangevinIntegrator(
    300*unit.kelvin,
    1.0/unit.picoseconds,
    0.002*unit.picoseconds
)

# Set up the simulation
simulation = app.Simulation(modeller.topology, system, integrator)
#simulation = app.Simulation(prmtop.topology, system, integrator)

# Set initial positions from Amber coordinates
simulation.context.setPositions(modeller.positions)


# Double check system
# Extract the system from the simulation
system_tmp = simulation.system

# Identify the HarmonicBondForce (which stores bonded interactions)
bond_force = None
for i in range(system_tmp.getNumForces()):
    force = system_tmp.getForce(i)
    if isinstance(force, mm.HarmonicBondForce):
        bond_force = force
        break

# Ensure bond force is found
if bond_force is not None:
    total_bonds = bond_force.getNumBonds()
    missing_bond_count = 0

    # Iterate over all bonds in the system_tmp
    for i in range(total_bonds):
        particle1, particle2, bond_length, bond_k = bond_force.getBondParameters(i)

        # Check if the bond force constant (bond_k) is zero or missing
        if bond_k._value == 0.0:
            missing_bond_count += 1

    # Compute missing bond percentage
    missing_percentage = (missing_bond_count / total_bonds) * 100 if total_bonds > 0 else 0

    # Print bond report
    print(f"Total number of bonds in the simulation: {total_bonds}")
    print(f"Number of bonds with missing force constants: {missing_bond_count}")
    print(f"Percentage of missing bond parameters: {missing_percentage:.2f}%")
else:
    print("No HarmonicBondForce found in the simulation system_tmp.")



# Minimize energy
print(f'Minimizing {minimize_nsteps} steps ...', flush=True)
start_time = time.time()

# Minimize until energy convergence
prev_energy = float('inf')
tolerance = 10.0  # kJ/mol threshold for convergence
for i in range(0, minimize_nsteps, 500):
    simulation.minimizeEnergy(maxIterations=500)
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    if abs(prev_energy - energy) < tolerance:
        print(f'Converged (tolerance: {tolerance:.2f} kJ/mol) at step: {i}')
        break
    prev_energy = energy
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.6f} seconds")

def save_imaged_pdb(simulation, filename):
    """
    Extracts topology and positions from an OpenMM simulation, applies imaging, and saves a PDB file.
    
    Parameters:
    - simulation: OpenMM Simulation object
    - filename: str, Output filename for the PDB file
    """
    # Extract topology from OpenMM simulation
    topology = md.Topology.from_openmm(simulation.topology)

    # Get positions and periodic box vectors from OpenMM
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)  # Convert to NumPy (nm)

    # Get periodic box vectors
    box_vectors = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometers)  # Convert to nm
    box_lengths = np.array([box_vectors[0][0], box_vectors[1][1], box_vectors[2][2]])  # Extract box dimensions

    # OpenMM assumes an orthorhombic box, so angles are always 90 degrees
    box_angles = np.array([[90.0, 90.0, 90.0]])  # Degrees

    # Create MDTraj trajectory with periodic box information
    traj = md.Trajectory(positions[np.newaxis, :, :], topology)
    traj.unitcell_lengths = box_lengths[np.newaxis, :]
    traj.unitcell_angles = box_angles  # Explicitly set unit cell angles

    # Ensure molecules are imaged correctly
    traj.image_molecules(inplace=True, make_whole=True)

    # Save to PDB
    traj.save_pdb(filename)

    print(f"Saved imaged PDB: {filename}")

save_imaged_pdb(simulation,"minimized_with_unique_residues.pdb")

# Save initial state to XML
with open("initial_state.xml", "w") as f:
    f.write(mm.XmlSerializer.serialize(simulation.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)))

def monitor_rmsd_equilibration(simulation, reference_xml, threshold=0.05, max_steps=100000, report_interval=1000):
    """Monitor RMSD to determine equilibration convergence."""
    rmsd_values = []
    converged = False

    for step in range(0, max_steps, report_interval):
        simulation.step(report_interval)

        # Extract positions and save temporary XML
        state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
        with open("tmp_state.xml", "w") as f:
            f.write(mm.XmlSerializer.serialize(state))

        # Load reference state
        with open(reference_xml, "r") as f:
            reference_state = mm.XmlSerializer.deserialize(f.read())

        ref_positions = np.array(reference_state.getPositions(asNumpy=True))
        cur_positions = np.array(state.getPositions(asNumpy=True))
        rmsd = np.sqrt(np.mean((cur_positions - ref_positions) ** 2))
        rmsd_values.append(rmsd)

        print(f"Step {step + report_interval}: RMSD = {rmsd:.3f} nm", flush=True)

        if len(rmsd_values) > 5 and all(r < threshold for r in rmsd_values[-5:]):
            print(f'Equilibration converged at step {step + report_interval} with RMSD {rmsd:.3f} nm')
            converged = True
            break

    os.remove("tmp_state.xml")
    if not converged:
        print(f"Maximum equilibration steps {max_steps} reached without full convergence.")

# NVT Equilibration with RMSD monitoring
print('Equilibrating (NVT) with RMSD monitoring...', flush=True)
start_time = time.time()
monitor_rmsd_equilibration(simulation, "initial_state.xml", 0.05, nvt_equil_nsteps)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.6f} seconds")

# Save the state after NVT equilibration
with open("nvt_equilibrated.xml", "w") as f:
    f.write(mm.XmlSerializer.serialize(simulation.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)))
    
# NPT Equilibration with RMSD monitoring
print('Equilibrating (NPT) with RMSD monitoring...', flush=True)
system.addForce(mm.MonteCarloBarostat(1 * unit.atmospheres, 300 * unit.kelvin, 25))
simulation.context.reinitialize(preserveState=True)
start_time = time.time()
monitor_rmsd_equilibration(simulation, "nvt_equilibrated.xml", 0.05, ntp_equil_nsteps)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.6f} seconds")

# Save final equilibrated positions
save_imaged_pdb(simulation,"equilibrated_with_NVT+NPT.pdb")

report_interval = min(production_nsteps,1000)
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
#simulation.reporters.append(app.PDBReporter('output.pdb', report_interval))
simulation.reporters.append(app.StateDataReporter(sys.stdout, report_interval, step=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(app.StateDataReporter('energies.csv', report_interval, step=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(app.DCDReporter('output.dcd', report_interval))
print (f'Running Production NPT Simulation - {production_nsteps * 0.002} ps ....', flush=True)
start_time = time.time()
simulation.step(production_nsteps)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.6f} seconds")

# Save to PDB
filename = f'final-{production_nsteps * 0.002}ps.pdb'  # Format to one decimal place
save_imaged_pdb(simulation,"equilibrated_with_NVT+NPT.pdb")

# Note: to ouput a parmtop file, use ParMed - which requires a conversion to non-rigid Water
#       Reason: PDB files are often used to provide a topology file to MDTraj.
#           However, it would seem that when the list of residues is long (e.g. water box) PDB file limits imprede.
tmpSystem = forcefield.createSystem (simulation.topology, rigidWater=False)
positions = simulation.context.getState(getPositions=True).getPositions()   # update to remove enforcePeriodicBox=True
structure = pmd.openmm.topsystem.load_topology(simulation.topology, system=tmpSystem, xyz=positions)
structure.save("output.prmtop", format="amber") 
structure.save("output.inpcrd", format="rst7") 

sys.exit()

t = md.load('output.pdb')
print (t)
t.image_molecules(inplace=True,make_whole=True)
print (t)
t.save_pdb('trajectory.pdb')
os.remove("output.pdb")

import gzip
import shutil
import threading

# Function to compress a file
def compress_file(input_filename, output_filename):
    print(f"Compressing {input_filename}...")
    with open(input_filename, 'rb') as f_in:
        with gzip.open(output_filename, 'wb', compresslevel=9) as f_out:
            shutil.copyfileobj(f_in, f_out)
    print(f"{input_filename} compressed to {output_filename}.")

# Create threads for each compression task
thread1 = threading.Thread(target=compress_file, args=('output.dcd', 'output.dcd.gz'))
thread2 = threading.Thread(target=compress_file, args=('trajectory.pdb', 'trajectory.pdb.gz'))

# Start both threads
thread1.start()
thread2.start()

# Wait for both threads to complete
thread1.join()
thread2.join()

print("Both files compressed successfully.")
os.remove("output.dcd")
os.remove("trajectory.pdb")

