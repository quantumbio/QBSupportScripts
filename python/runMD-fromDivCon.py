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
import csv

import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run OpenMM MD simulation with optional test mode.")
parser.add_argument("prmtop", help="Path to AMBER prmtop file")
parser.add_argument("inpcrd", help="Path to AMBER inpcrd file")
parser.add_argument("--test-run", action="store_true", help="Run a short test simulation")
parser.add_argument("--medium-run", action="store_true", help="Run a medium simulation")
parser.add_argument(
    "--minimization-steps",
    type=int,
    metavar="N",
    help="Override the maximum number of OpenMM minimization iterations",
)
parser.add_argument(
    "--equilibration-steps",
    type=int,
    metavar="N",
    help="Override the number of steps used for each NVT and NPT equilibration stage",
)
parser.add_argument(
    "--production-steps",
    type=int,
    metavar="N",
    help="Override the number of production MD steps",
)
parser.add_argument(
    "--report-interval",
    type=int,
    metavar="N",
    help="Override the production trajectory/metrics reporting interval in MD steps",
)
parser.add_argument("--skip-waterbox", action="store_true", help="Skip adding the water box")

args = parser.parse_args()
prmtopFile = args.prmtop
inpcrdFile = args.inpcrd
test_run = args.test_run
medium_run = args.medium_run
skip_waterbox = args.skip_waterbox

if test_run:
    minimize_nsteps   = 500  # ~quick minimization
    nvt_equil_nsteps  = 1000
    ntp_equil_nsteps  = 1000
    production_nsteps = 5000  # 10 ps at 0.002 ps timestep
    preport_interval  = 100
    print("Test run mode enabled: using reduced step counts for quick validation.")
elif medium_run:
    minimize_nsteps   = 50000
    nvt_equil_nsteps  = 50000
    ntp_equil_nsteps  = 50000
    production_nsteps = 200000
    preport_interval  = 10000
    print("Medium run mode enabled: using intermediate step counts.")
else:
    minimize_nsteps   = round(5000 / 10)
    nvt_equil_nsteps  = round(100000 / 10)
    ntp_equil_nsteps  = round(100000 / 10)
    production_nsteps = 25000000  # 50 ns for 0.002 ps timestep
    preport_interval  = 25000

# Explicit step counts override the preset selected above.  A single
# equilibration value is intentionally applied to both NVT and NPT so this
# interface mirrors qmechanic's --moleculardynamics EQUIL,PRODUCTION settings.
if args.minimization_steps is not None:
    if args.minimization_steps <= 0:
        parser.error("--minimization-steps must be greater than zero")
    minimize_nsteps = args.minimization_steps
    print(f"Minimization step override: {args.minimization_steps} maximum iterations.")

if args.equilibration_steps is not None:
    if args.equilibration_steps <= 0:
        parser.error("--equilibration-steps must be greater than zero")
    nvt_equil_nsteps = args.equilibration_steps
    ntp_equil_nsteps = args.equilibration_steps
    print(f"Equilibration step override: {args.equilibration_steps} steps for both NVT and NPT.")

if args.production_steps is not None:
    if args.production_steps <= 0:
        parser.error("--production-steps must be greater than zero")
    production_nsteps = args.production_steps
    print(f"Production step override: {args.production_steps} steps.")

if args.report_interval is not None:
    if args.report_interval <= 0:
        parser.error("--report-interval must be greater than zero")
    preport_interval = args.report_interval
    print(f"Production report interval override: {args.report_interval} steps.")

# Never configure a reporter interval beyond the whole production run.
# This guarantees at least one trajectory/report frame for short validation runs.
if preport_interval > production_nsteps:
    print(
        f"Production report interval {preport_interval} exceeds production length "
        f"{production_nsteps}; using {production_nsteps} steps."
    )
    preport_interval = production_nsteps

# Load the Amber topology and coordinate files
prmtop = pmd.load_file(prmtopFile, inpcrdFile)

# Identify and modify HOH residues
for residue in prmtop.residues:
    if residue.name == "HOH":
        # Identify H1 and H2
        h1 = residue.atoms[1]  # Adjust index if necessary
        h2 = residue.atoms[2]  # Adjust index if necessary
        
        # Find and remove the bond between H1 and H2
        for bond in h1.bonds:
            if bond.atom2 == h2 or bond.atom1 == h2:
                # Remove the bond from both atoms
                h1.bonds.remove(bond)
                h2.bonds.remove(bond)
                break  # Exit loop after removing the bond

#prmtop.save('initial.pdb')
from openmm.app import PDBFile
with open('initial.pdb', 'w') as f:
    PDBFile.writeFile(prmtop.topology, prmtop.positions, f)

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

# Use the same TIP3P water family as the DivCon waterbox builder.  Keep the
# filename in one place so residue-template discovery and final ForceField
# construction can not silently drift to different water models.
WATER_MODEL_XML = 'amber14/tip3p.xml'
GENERATED_FORCEFIELD_XML = 'complete_forcefield_with_unique_residues.xml'
WATER_RESIDUE_NAMES = {"HOH", "WAT"}

# Preserve the atom-indexed DivCon/AMBER parameters before residue renaming and
# OpenMM template assignment.  Atom ordering is unchanged by the preparation
# performed below, so these provide an independent congruence reference.
divcon_input_atom_charges = [float(atom.charge) for atom in prmtop.atoms]
divcon_input_atom_types = [str(atom.type) for atom in prmtop.atoms]

# Extract defined residues from the selected OpenMM water/ion force field.
defined_residues = set()
forcefield_xml_file = importlib.resources.files('openmm.app.data') / WATER_MODEL_XML
print(f"OpenMM water/ion template force field: {WATER_MODEL_XML}")

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
param_set.write(GENERATED_FORCEFIELD_XML)

# Load the force field XML file
tree = ET.parse(GENERATED_FORCEFIELD_XML)
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
tree.write(GENERATED_FORCEFIELD_XML, xml_declaration=True, encoding='utf-8', method="xml")

# Load the modified topology, but retain coordinates/box information from the
# original supplied inpcrd.  The residue-renaming/topology preparation above
# does not change atom count/order or coordinates, and rewriting an inpcrd via
# ParmEd may omit optional periodic box records.
prmtop = app.AmberPrmtopFile('modified_with_unique_residues.prmtop')
inpcrd = app.AmberInpcrdFile(inpcrdFile)

# Prefer the box carried by the authoritative input coordinate file.  If that
# file does not carry box vectors, fall back to the box retained in the
# prepared prmtop topology (AMBER prmtop files can also contain box data).
input_box_vectors = inpcrd.boxVectors
input_box_source = "input inpcrd"
if input_box_vectors is None:
    input_box_vectors = prmtop.topology.getPeriodicBoxVectors()
    input_box_source = "prepared prmtop"

# Extract positions from inpcrd
positions = inpcrd.getPositions()

# Convert the positions to a NumPy array for easy manipulation
positions_array = np.array([[pos[0].value_in_unit(unit.nanometers),
                              pos[1].value_in_unit(unit.nanometers),
                              pos[2].value_in_unit(unit.nanometers)] for pos in positions])

# Calculate coordinate bounds.  These are only needed when this script is
# responsible for creating a new solvent box.
x_coords = positions_array[:, 0]
y_coords = positions_array[:, 1]
z_coords = positions_array[:, 2]

min_x, max_x = min(x_coords), max(x_coords)
min_y, max_y = min(y_coords), max(y_coords)
min_z, max_z = min(z_coords), max(z_coords)

# The generated XML carries the DivCon-transferred AMBER atom classes and
# nonbonded parameters (including OW/HW), while the stock TIP3P file supplies
# the standard HOH/WAT and ion residue templates.
forcefield = app.ForceField(WATER_MODEL_XML, GENERATED_FORCEFIELD_XML)

# Use Modeller for the final topology/positions passed to OpenMM.
modeller = app.Modeller(prmtop.topology, inpcrd.positions)

if skip_waterbox:
    # Preserve an explicit periodic cell when the AMBER inputs carry one.
    # Some DivCon-generated parm7/inpcrd pairs do not currently serialize box
    # metadata.  In that case, mirror MDDriver's fallback: build an
    # orthorhombic box from the coordinate extrema with 1.0 A padding on each
    # side.  This replaces the old (and erroneous) 10 nm padding.
    if input_box_vectors is None:
        fallback_padding_angstrom = 1.0
        fallback_padding_nm = fallback_padding_angstrom * 0.1
        box_size_x = (max_x - min_x) + 2.0 * fallback_padding_nm
        box_size_y = (max_y - min_y) + 2.0 * fallback_padding_nm
        box_size_z = (max_z - min_z) + 2.0 * fallback_padding_nm
        input_box_vectors = (
            mm.Vec3(box_size_x, 0.0, 0.0) * unit.nanometer,
            mm.Vec3(0.0, box_size_y, 0.0) * unit.nanometer,
            mm.Vec3(0.0, 0.0, box_size_z) * unit.nanometer,
        )
        input_box_source = (
            f"coordinate bounds + {fallback_padding_angstrom:.1f} A padding/side"
        )

    modeller.topology.setPeriodicBoxVectors(input_box_vectors)
    box_dimensions = [
        input_box_vectors[0][0].value_in_unit(unit.nanometers),
        input_box_vectors[1][1].value_in_unit(unit.nanometers),
        input_box_vectors[2][2].value_in_unit(unit.nanometers),
    ]
    print(
        f"Skipping water box addition; using {input_box_source} periodic box (nm):",
        box_dimensions,
    )
else:
    # When Python creates the solvent box, retain the established 1 nm padding
    # behavior.  Modeller.addSolvent() sets the resulting periodic box.
    padding = 1.0 * unit.nanometers
    box_size_x = (max_x - min_x) + 2 * padding.value_in_unit(unit.nanometers)
    box_size_y = (max_y - min_y) + 2 * padding.value_in_unit(unit.nanometers)
    box_size_z = (max_z - min_z) + 2 * padding.value_in_unit(unit.nanometers)
    box_vector = mm.Vec3(box_size_x, box_size_y, box_size_z)
    modeller.addSolvent(
        forcefield,
        model='tip3p',
        boxSize=box_vector,
        ionicStrength=0.15*unit.molar,
    )
    final_box_vectors = modeller.topology.getPeriodicBoxVectors()
    if final_box_vectors is not None:
        print(
            "Generated periodic box (nm):",
            [
                final_box_vectors[0][0].value_in_unit(unit.nanometers),
                final_box_vectors[1][1].value_in_unit(unit.nanometers),
                final_box_vectors[2][2].value_in_unit(unit.nanometers),
            ],
        )
# Now create the system again after modifying the topology
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometers,
    constraints=app.HBonds,
    rigidWater=True,
    removeCMMotion=True,
    ewaldErrorTolerance=1.0e-4,
)

def _find_nonbonded_force(openmm_system):
    forces = [
        openmm_system.getForce(i)
        for i in range(openmm_system.getNumForces())
        if isinstance(openmm_system.getForce(i), mm.NonbondedForce)
    ]
    if len(forces) != 1:
        raise RuntimeError(
            f"Expected exactly one OpenMM NonbondedForce, found {len(forces)}"
        )
    return forces[0]

def _generated_nonbonded_class_parameters(xml_filename, wanted_classes):
    """Read sigma/epsilon (OpenMM XML units) for selected transferred classes."""
    xml_root = ET.parse(xml_filename).getroot()
    nb_element = xml_root.find("NonbondedForce")
    if nb_element is None:
        raise RuntimeError(f"No NonbondedForce found in {xml_filename}")

    parameters = {}
    for atom_element in nb_element.findall("Atom"):
        atom_class = atom_element.get("class")
        if atom_class in wanted_classes:
            parameters[atom_class] = (
                float(atom_element.get("sigma")),
                float(atom_element.get("epsilon")),
            )
    missing = set(wanted_classes) - set(parameters)
    if missing:
        raise RuntimeError(
            f"Missing transferred nonbonded classes in {xml_filename}: "
            + ", ".join(sorted(missing))
        )
    return parameters

def align_and_report_nonbonded_parameters(openmm_system, topology):
    """
    Align existing DivCon water particles to the DivCon-transferred OW/HW
    nonbonded parameters, then report the effective OpenMM nonbonded setup.

    The stock amber14/tip3p.xml file supplies the HOH residue template and the
    rigid TIP3P geometry.  The generated force-field XML contains the exact
    OW/HW Lennard-Jones values transferred from the DivCon input.  Applying
    those values here makes the Python validation Hamiltonian use the same
    per-particle water charge/LJ parameters as MDDriver while retaining the
    standard TIP3P topology/constraint geometry.
    """
    nonbonded = _find_nonbonded_force(openmm_system)
    transferred_water_lj = _generated_nonbonded_class_parameters(
        GENERATED_FORCEFIELD_XML, {"OW", "HW"}
    )

    topology_atoms = list(topology.atoms())
    can_map_to_input = (
        len(topology_atoms) == len(divcon_input_atom_charges)
        and openmm_system.getNumParticles() == len(divcon_input_atom_charges)
    )

    water_particle_count = 0
    if skip_waterbox and not can_map_to_input:
        raise RuntimeError(
            "Prepared topology atom count/order no longer matches the supplied "
            "DivCon topology; can not safely transfer water parameters by index."
        )

    # For the validation path (--skip-waterbox), transfer the authoritative
    # DivCon charge plus OW/HW LJ values to every existing water particle.
    # For Python-generated solvent, retain the selected stock TIP3P particles.
    if skip_waterbox:
        for atom in topology_atoms:
            if atom.residue.name not in WATER_RESIDUE_NAMES:
                continue

            if atom.element is not None and atom.element.symbol == "O":
                atom_class = "OW"
            elif atom.element is not None and atom.element.symbol == "H":
                atom_class = "HW"
            else:
                raise RuntimeError(
                    f"Unexpected atom {atom.name} in water residue {atom.residue.name}"
                )

            sigma_nm, epsilon_kj = transferred_water_lj[atom_class]
            charge_e = divcon_input_atom_charges[atom.index]
            nonbonded.setParticleParameters(
                atom.index,
                charge_e * unit.elementary_charge,
                sigma_nm * unit.nanometer,
                epsilon_kj * unit.kilojoule_per_mole,
            )
            water_particle_count += 1

    method_names = {
        mm.NonbondedForce.NoCutoff: "NoCutoff",
        mm.NonbondedForce.CutoffNonPeriodic: "CutoffNonPeriodic",
        mm.NonbondedForce.CutoffPeriodic: "CutoffPeriodic",
        mm.NonbondedForce.Ewald: "Ewald",
        mm.NonbondedForce.PME: "PME",
        mm.NonbondedForce.LJPME: "LJPME",
    }
    print("OpenMM nonbonded configuration:")
    print(f"  water template XML       : {WATER_MODEL_XML}")
    print(f"  transferred parameter XML: {GENERATED_FORCEFIELD_XML}")
    print(f"  method                   : {method_names.get(nonbonded.getNonbondedMethod(), nonbonded.getNonbondedMethod())}")
    print(f"  cutoff (nm)              : {nonbonded.getCutoffDistance().value_in_unit(unit.nanometer):.8f}")
    print(f"  Ewald error tolerance    : {nonbonded.getEwaldErrorTolerance():.8g}")
    print(f"  dispersion correction    : {nonbonded.getUseDispersionCorrection()}")
    print(f"  switching function       : {nonbonded.getUseSwitchingFunction()}")
    if nonbonded.getUseSwitchingFunction():
        print(f"  switching distance (nm)  : {nonbonded.getSwitchingDistance().value_in_unit(unit.nanometer):.8f}")
    print(f"  nonbonded particles      : {nonbonded.getNumParticles()}")
    print(f"  nonbonded exceptions     : {nonbonded.getNumExceptions()}")
    if skip_waterbox:
        print(f"  DivCon water particles aligned: {water_particle_count}")

    total_charge_e = 0.0
    for particle_index in range(nonbonded.getNumParticles()):
        charge, _, _ = nonbonded.getParticleParameters(particle_index)
        total_charge_e += charge.value_in_unit(unit.elementary_charge)
    print(f"  total particle charge (e): {total_charge_e:.10f}")

    print("DivCon-transferred OW/HW LJ parameters:")
    for atom_class in ("OW", "HW"):
        sigma_nm, epsilon_kj = transferred_water_lj[atom_class]
        print(
            f"  {atom_class}: sigma={sigma_nm:.12f} nm  "
            f"epsilon={epsilon_kj:.12f} kJ/mol"
        )

    # Report actual parameters on one representative water and the first Na/Cl
    # ions, if present.  This reports the final values OpenMM will actually use.
    representative_atoms = []
    first_water_residue = next(
        (r for r in topology.residues() if r.name in WATER_RESIDUE_NAMES), None
    )
    if first_water_residue is not None:
        representative_atoms.extend(list(first_water_residue.atoms()))

    seen_ion_residues = set()
    for atom in topology_atoms:
        residue_name = atom.residue.name.upper()
        if residue_name in {"NA", "CL"} and residue_name not in seen_ion_residues:
            representative_atoms.append(atom)
            seen_ion_residues.add(residue_name)
        if len(seen_ion_residues) == 2:
            break

    if representative_atoms:
        print("Representative effective OpenMM particle parameters:")
        for atom in representative_atoms:
            charge, sigma, epsilon = nonbonded.getParticleParameters(atom.index)
            input_type = (
                divcon_input_atom_types[atom.index] if atom.index < len(divcon_input_atom_types)
                else "<generated>"
            )
            input_charge = (
                divcon_input_atom_charges[atom.index] if atom.index < len(divcon_input_atom_charges)
                else float("nan")
            )
            print(
                f"  index={atom.index:5d} {atom.residue.name}:{atom.name:<4s} "
                f"DivConType={input_type:<8s} "
                f"DivConQ={input_charge:+.8f} "
                f"OpenMMQ={charge.value_in_unit(unit.elementary_charge):+.8f} "
                f"sigma={sigma.value_in_unit(unit.nanometer):.12f} nm "
                f"epsilon={epsilon.value_in_unit(unit.kilojoule_per_mole):.12f} kJ/mol"
            )

    if first_water_residue is not None:
        water_indices = {atom.index for atom in first_water_residue.atoms()}
        print("Representative water constraints:")
        for constraint_index in range(openmm_system.getNumConstraints()):
            atom1, atom2, distance = openmm_system.getConstraintParameters(constraint_index)
            if atom1 in water_indices and atom2 in water_indices:
                print(
                    f"  {atom1}-{atom2}: "
                    f"{distance.value_in_unit(unit.nanometer):.12f} nm"
                )

    return nonbonded

nonbonded_force = align_and_report_nonbonded_parameters(system, modeller.topology)

# Match the C++ MDDriver force-group layout so that OpenMM energies can be
# compared component-by-component without changing the Hamiltonian.  Any force
# type outside these four expected AMBER/OpenMM terms is isolated in group 4
# and reported separately as a diagnostic.
for force_index in range(system.getNumForces()):
    force = system.getForce(force_index)
    if isinstance(force, mm.HarmonicBondForce):
        force.setForceGroup(0)
    elif isinstance(force, mm.HarmonicAngleForce):
        force.setForceGroup(1)
    elif isinstance(force, mm.PeriodicTorsionForce):
        force.setForceGroup(2)
    elif isinstance(force, mm.NonbondedForce):
        force.setForceGroup(3)
    else:
        force.setForceGroup(4)

def report_openmm_energy_components(simulation, label):
    """Report OpenMM potential-energy components using the C++ force groups."""
    component_groups = (
        ("Bond", 0),
        ("Angle", 1),
        ("Torsion", 2),
        ("Nonbonded", 3),
    )

    print(f"OpenMM energy decomposition - {label} (kJ/mol):")
    component_sum = 0.0
    for component_name, group in component_groups:
        state = simulation.context.getState(getEnergy=True, groups=(1 << group))
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        component_sum += energy
        print(f"  {component_name:<12}: {energy:.6f}")

    other_state = simulation.context.getState(getEnergy=True, groups=(1 << 4))
    other_energy = other_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    total_state = simulation.context.getState(getEnergy=True)
    total_energy = total_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    if abs(other_energy) > 1.0e-8:
        print(f"  {'Other':<12}: {other_energy:.6f}")
    print(f"  {'4-part sum':<12}: {component_sum:.6f}")
    print(f"  {'Total':<12}: {total_energy:.6f}")
    print(f"  {'Residual':<12}: {(total_energy - component_sum - other_energy):.6e}")

# Set up the integrator
integrator = mm.LangevinIntegrator(
    298*unit.kelvin,
    1.0/unit.picoseconds,
    0.002*unit.picoseconds
)

# Set up the simulation
simulation = app.Simulation(modeller.topology, system, integrator)
#simulation = app.Simulation(prmtop.topology, system, integrator)

# Set initial positions from Amber coordinates
simulation.context.setPositions(modeller.positions)

# Record the actual OpenMM execution platform and the authoritative starting
# periodic box used by this Context.  These are high-value congruence checks.
platform = simulation.context.getPlatform()
print(f"OpenMM execution platform: {platform.getName()}")
print(f"OpenMM particles: {system.getNumParticles()}")
print(f"OpenMM constraints: {system.getNumConstraints()}")
initial_state = simulation.context.getState(getEnergy=True, getPositions=True)
initial_box = initial_state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometers)
print(
    "OpenMM initial box (nm):",
    [float(initial_box[0][0]), float(initial_box[1][1]), float(initial_box[2][2])],
)
print(
    "OpenMM initial potential energy (kJ/mol):",
    initial_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
)
report_openmm_energy_components(simulation, "initial")

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



# Minimize energy using the same OpenMM LocalEnergyMinimizer semantics as the
# C++ driver: 10 kJ/mol/nm RMS-force tolerance and a maximum iteration count.
minimization_tolerance = 10.0
print(
    f'Minimizing with OpenMM: tolerance={minimization_tolerance:.2f} kJ/mol/nm, '
    f'maxIterations={minimize_nsteps} ...',
    flush=True,
)
start_time = time.time()
mm.LocalEnergyMinimizer.minimize(
    simulation.context,
    minimization_tolerance,
    minimize_nsteps,
)
minimized_state = simulation.context.getState(getEnergy=True)
print(
    "Post-minimization Potential Energy: "
    f"{minimized_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.2f} kJ/mol"
)
report_openmm_energy_components(simulation, "post-minimization")
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.6f} seconds")

# Match the C++ protocol: initialize velocities once after minimization and
# carry them continuously through NVT, NPT, and production.
simulation.context.setVelocitiesToTemperature(298*unit.kelvin)

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

def get_equilibration_atom_indices(simulation):
    """Select atoms for equilibration RMSD (prefer solute heavy atoms)."""
    md_topology = md.Topology.from_openmm(simulation.topology)

    # Prefer protein heavy atoms for biomolecular systems.
    atom_indices = md_topology.select("protein and not element H")
    if len(atom_indices) > 0:
        return atom_indices

    # Fallback to non-water heavy atoms (covers ligand-only and mixed systems).
    atom_indices = md_topology.select("not water and not element H")
    if len(atom_indices) > 0:
        return atom_indices

    # Final fallback if the system has no heavy-atom selection.
    atom_indices = md_topology.select("not element H")
    if len(atom_indices) > 0:
        return atom_indices

    return np.arange(md_topology.n_atoms)


def monitor_rmsd_equilibration(simulation, reference_xml, threshold=0.05, max_steps=100000, report_interval=1000, atom_indices=None):
    """Run exactly max_steps while monitoring aligned RMSD on a meaningful atom subset."""
    rmsd_values = []
    convergence_reported = False

    md_topology = md.Topology.from_openmm(simulation.topology)
    if atom_indices is None:
        atom_indices = get_equilibration_atom_indices(simulation)

    with open(reference_xml, "r") as f:
        reference_state = mm.XmlSerializer.deserialize(f.read())

    ref_positions = reference_state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
    ref_traj = md.Trajectory(np.array(ref_positions)[np.newaxis, :, :], md_topology)

    completed_steps = 0
    while completed_steps < max_steps:
        chunk = min(report_interval, max_steps - completed_steps)
        simulation.step(chunk)
        completed_steps += chunk

        # Extract current positions and compute aligned RMSD against the fixed reference.
        state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
        cur_positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
        cur_traj = md.Trajectory(np.array(cur_positions)[np.newaxis, :, :], md_topology)
        cur_traj.superpose(ref_traj, atom_indices=atom_indices)
        rmsd = md.rmsd(cur_traj, ref_traj, atom_indices=atom_indices)[0]
        rmsd_values.append(rmsd)

        print(f"Step {completed_steps}: RMSD = {rmsd:.3f} nm", flush=True)

        if (
            not convergence_reported
            and len(rmsd_values) > 5
            and all(r < threshold for r in rmsd_values[-5:])
        ):
            print(
                f'RMSD convergence criterion satisfied at step {completed_steps} '
                f'with RMSD {rmsd:.3f} nm; continuing to requested equilibration length.'
            )
            convergence_reported = True

    if not convergence_reported:
        print(
            f"Completed requested equilibration length of {max_steps} steps; "
            "RMSD convergence criterion was not reached."
        )


def compute_target_pressure_bar(barostat=None):
    """Return barostat target pressure in bar."""
    if barostat is None:
        return float("nan")
    return barostat.getDefaultPressure().value_in_unit(unit.bar)


def compute_instantaneous_pressure_bar(simulation):
    """Return instantaneous pressure in bar when available, otherwise NaN."""
    if not hasattr(mm.MonteCarloBarostat, "computeCurrentPressure"):
        return float("nan")
    try:
        pressure = mm.MonteCarloBarostat.computeCurrentPressure(simulation.context)
        return pressure.value_in_unit(unit.bar)
    except Exception:
        return float("nan")


def compute_structural_metrics(simulation, reference_positions_nm, atom_indices):
    """Compute aligned RMSD and radius of gyration (nm) on selected atoms."""
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    cur_positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)

    idx = np.asarray(atom_indices, dtype=np.int64)
    ref_subset = np.asarray(reference_positions_nm)[idx]
    cur_subset = np.asarray(cur_positions)[idx]

    # Center both subsets.
    ref_centered = ref_subset - ref_subset.mean(axis=0)
    cur_centered = cur_subset - cur_subset.mean(axis=0)

    # Kabsch alignment: rotate current subset onto reference subset.
    covariance = cur_centered.T @ ref_centered
    u_mat, _, v_t = np.linalg.svd(covariance)
    rotation = v_t.T @ u_mat.T
    if np.linalg.det(rotation) < 0:
        v_t[-1, :] *= -1
        rotation = v_t.T @ u_mat.T

    cur_aligned = cur_centered @ rotation
    diff = cur_aligned - ref_centered
    rmsd_nm = float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))

    # Radius of gyration from current coordinates (translation invariant).
    rg_nm = float(np.sqrt(np.mean(np.sum(cur_centered * cur_centered, axis=1))))
    return rmsd_nm, rg_nm


def compute_degrees_of_freedom(system):
    """Compute mechanical degrees of freedom for instantaneous temperature estimates."""
    dof = 0
    for i in range(system.getNumParticles()):
        if system.getParticleMass(i) > 0 * unit.dalton:
            dof += 3

    dof -= system.getNumConstraints()

    for i in range(system.getNumForces()):
        if isinstance(system.getForce(i), mm.CMMotionRemover):
            dof -= 3
            break

    return max(dof, 1)


def compute_density_g_per_ml(simulation, state):
    """Compute density from topology mass and periodic volume."""
    total_mass = 0.0 * unit.dalton
    for atom in simulation.topology.atoms():
        if atom.element is not None:
            total_mass += atom.element.mass

    # Atom masses are molar masses; convert to absolute mass via Avogadro's constant.
    absolute_mass = total_mass / unit.AVOGADRO_CONSTANT_NA
    density = absolute_mass / state.getPeriodicBoxVolume()
    return density.value_in_unit(unit.gram / unit.milliliter)


def compute_temperature_k(kinetic_energy, dof):
    """Compute instantaneous temperature from kinetic energy and DOF."""
    return (2 * kinetic_energy / (dof * unit.MOLAR_GAS_CONSTANT_R)).value_in_unit(unit.kelvin)

# NVT Equilibration with RMSD monitoring
print('Equilibrating (NVT) with RMSD monitoring...', flush=True)
start_time = time.time()
equilibration_atom_indices = get_equilibration_atom_indices(simulation)
print(f"Tracking RMSD over {len(equilibration_atom_indices)} selected atoms.", flush=True)
monitor_rmsd_equilibration(simulation, "initial_state.xml", 0.05, nvt_equil_nsteps, atom_indices=equilibration_atom_indices)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.6f} seconds")

# Save the state after NVT equilibration
with open("nvt_equilibrated.xml", "w") as f:
    f.write(mm.XmlSerializer.serialize(simulation.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)))
    
# NPT Equilibration with RMSD monitoring
print('Equilibrating (NPT) with RMSD monitoring...', flush=True)
barostat = mm.MonteCarloBarostat(1 * unit.bar, 298 * unit.kelvin, 25)
system.addForce(barostat)
simulation.context.reinitialize(preserveState=True)
start_time = time.time()
monitor_rmsd_equilibration(simulation, "nvt_equilibrated.xml", 0.05, ntp_equil_nsteps, atom_indices=equilibration_atom_indices)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.6f} seconds")

# Save final equilibrated positions
save_imaged_pdb(simulation,"equilibrated_with_NVT+NPT.pdb")

# Continue directly from the equilibrated NPT velocities into production.
# Reset step count so production logs/reporters are production-relative.
simulation.currentStep = 0
#simulation.reporters.append(app.PDBReporter('output.pdb', preport_interval))
simulation.reporters.append(app.StateDataReporter(
    'energies.csv',
    preport_interval,
    step=True,
    potentialEnergy=True,
    kineticEnergy=True,
    temperature=True,
    volume=True,
    density=True,
    elapsedTime=True,
))
# Keep a handle to the DCD reporter so the equilibrated production starting
# structure can be written explicitly as frame 0.  Subsequent frames are still
# written automatically at preport_interval (for example 0, 500, 1000, ...).
dcd_reporter = app.DCDReporter('output.dcd', preport_interval)
simulation.reporters.append(dcd_reporter)
print (f'Running Production NPT Simulation - {production_nsteps * 0.002} ps ....', flush=True)

# Capture the production reference before taking any production MD steps.  Use
# this exact same State both as the structural RMSD reference and as DCD frame 0
# so the first trajectory frame is the authoritative production starting point.
production_reference_state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
production_reference_positions = production_reference_state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
dcd_reporter.report(simulation, production_reference_state)
print("Saved production step 0 frame to output.dcd", flush=True)
production_dof = compute_degrees_of_freedom(system)
instantaneous_pressure_supported = hasattr(mm.MonteCarloBarostat, "computeCurrentPressure")

if not instantaneous_pressure_supported:
    print(
        "Note: instantaneous pressure not supported on this OpenMM build; reporting NaN for instantaneous_pressure_bar.",
        flush=True,
    )

print(
    "Production metrics (step, pe_kj_per_mol, ke_kj_per_mol, temp_k, volume_nm3, density_g_per_ml, target_pressure_bar, instantaneous_pressure_bar, rmsd_nm, rg_nm, wall_time_s):",
    flush=True,
)
start_time = time.time()
steps_remaining = production_nsteps
with open("production_metrics.csv", "w", newline="") as metrics_csv:
    writer = csv.writer(metrics_csv)
    writer.writerow([
        "step",
        "potential_energy_kj_per_mol",
        "kinetic_energy_kj_per_mol",
        "temperature_k",
        "volume_nm3",
        "density_g_per_ml",
        "target_pressure_bar",
        "instantaneous_pressure_bar",
        "pressure_source",
        "rmsd_nm",
        "rg_nm",
        "wall_time_s",
    ])

    while steps_remaining > 0:
        try:
            chunk = min(preport_interval, steps_remaining)
            simulation.step(chunk)
            steps_remaining -= chunk

            report_state = simulation.context.getState(getEnergy=True, enforcePeriodicBox=True)
            pe_kj_per_mol = report_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            ke_quantity = report_state.getKineticEnergy()
            ke_kj_per_mol = ke_quantity.value_in_unit(unit.kilojoule_per_mole)
            temp_k = compute_temperature_k(ke_quantity, production_dof)
            volume_nm3 = report_state.getPeriodicBoxVolume().value_in_unit(unit.nanometer**3)
            density_g_per_ml = compute_density_g_per_ml(simulation, report_state)
            target_pressure_bar = compute_target_pressure_bar(barostat)
            instantaneous_pressure_bar = compute_instantaneous_pressure_bar(simulation)
            pressure_source = "instantaneous" if not np.isnan(instantaneous_pressure_bar) else "target_only"
            rmsd_nm, rg_nm = compute_structural_metrics(
                simulation,
                production_reference_positions,
                equilibration_atom_indices,
            )
            wall_time_s = time.time() - start_time

            writer.writerow([
                simulation.currentStep,
                pe_kj_per_mol,
                ke_kj_per_mol,
                temp_k,
                volume_nm3,
                density_g_per_ml,
                target_pressure_bar,
                instantaneous_pressure_bar,
                pressure_source,
                rmsd_nm,
                rg_nm,
                wall_time_s,
            ])

            print(
                f"  {simulation.currentStep}, {pe_kj_per_mol:.6f}, {ke_kj_per_mol:.6f}, {temp_k:.3f}, {volume_nm3:.6f}, {density_g_per_ml:.6f}, {target_pressure_bar:.6f}, {instantaneous_pressure_bar:.6f}, {rmsd_nm:.6f}, {rg_nm:.6f}, {wall_time_s:.3f}",
                flush=True,
            )
        except Exception as exc:
            print(
                f"ERROR during production metrics at step {simulation.currentStep}: {exc}",
                flush=True,
            )
            raise
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

summary_state = simulation.context.getState(getEnergy=True, enforcePeriodicBox=True)
summary_values = {
    "step": simulation.currentStep,
    "potential_energy_kj_per_mol": summary_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
    "kinetic_energy_kj_per_mol": summary_state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole),
    "temperature_k": compute_temperature_k(summary_state.getKineticEnergy(), production_dof),
    "volume_nm3": summary_state.getPeriodicBoxVolume().value_in_unit(unit.nanometer**3),
    "density_g_per_ml": compute_density_g_per_ml(simulation, summary_state),
    "target_pressure_bar": compute_target_pressure_bar(barostat),
    "instantaneous_pressure_bar": compute_instantaneous_pressure_bar(simulation),
    "rmsd_nm": compute_structural_metrics(simulation, production_reference_positions, equilibration_atom_indices)[0],
    "rg_nm": compute_structural_metrics(simulation, production_reference_positions, equilibration_atom_indices)[1],
    "wall_time_s": elapsed_time,
}

print("Simulation summary:")
for label, value in summary_values.items():
    print(f"  {label}: {value}")

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

