#!/bin/bash
set -e

#  // BEGIN COPYRIGHT
#  /***********************************************************************
#     Copyright (c) 2021 QuantumBio Inc. and/or its affiliates.
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
#  // END COPYRIGHT


if [ -z "${QBHOME}" ]; then
    echo "ERROR: QBHOME is not set! You MUST source /path/to/DivConSuite/etc/qbenv.sh THEN call this function"
    exit
fi

if [ -z "${PHENIX}" ]; then
    echo "ERROR: PHENIX is not set! You MUST source /path/to/phenix-1.XX.XXX/phenix_env.sh THEN call this function"
    exit
fi

DIVCON_BIN=${QBHOME}/bin/qmechanic

if [ ! -x "${DIVCON_BIN}" ]; then
    echo "ERROR: ${DIVCON_BIN} is not an executable! Set this path to the qmechanic executable"
    exit
fi

currentDate=`date`
echo "BEGIN PHENIX/DivCon Tutorial Test at ${currentDate} using ${DIVCON_BIN}"

WORKDIR=$PWD
PDBID=2WOR
echo "Tutorial #1 (Running Phenix/DivCon without the qbphenix wrapper): $PDBID"
rm -rf ${WORKDIR}/NakedPHENIX ; mkdir -p ${WORKDIR}/NakedPHENIX ; cd ${WORKDIR}/NakedPHENIX

wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/2WOR.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/2WOR-sf.cif

${PHENIX}/build/bin/phenix.ready_set ${PDBID}.pdb add_h_to_water=true

${PHENIX}/build/bin/phenix.refine ${PDBID}.updated.pdb ${PDBID}-sf.cif ${PDBID}.ligands.cif output.write_geo_file=False  output.write_def_file=False refinement.refine.strategy=individual_sites+individual_adp main.number_of_macro_cycles=2 qblib=True qblib_method=pm6 qblib_region_selection="chain A and resname 2AN and resid 1098" qblib_region_radius=3.0 qblib_buffer_radius=2.5 qblib_np=4 qblib_target_method=3 qblib_gradient_target=3

echo "Tutorial #2a (using qbphenix execution script and MOE on single ligand): 3IX1"
rm -rf ${WORKDIR}/qbphenix_moe ; mkdir -p ${WORKDIR}/qbphenix_moe ; cd ${WORKDIR}/qbphenix_moe

wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/3ix1.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/3ix1-sf.cif

$QBHOME/bin/qbphenix --dataFile 3ix1-sf.cif --pdbFile 3ix1.pdb --selection "chain A resname NFM resid 401"  --phenixOptions "main.number_of_macro_cycles=2" --protonation MOE --qmMethod pm6 --region-radius 3.0 --buffer-radius 2.5 --Nproc 4

echo "Tutorial #2b (using qbphenix execution script and phenix.ready_set on single ligand): 3IX1"
rm -rf ${WORKDIR}/qbphenix_readyset ; mkdir -p ${WORKDIR}/qbphenix_readyset ; cd ${WORKDIR}/qbphenix_readyset

$QBHOME/bin/qbphenix --pdbID 3ix1 --selection "chain A resname NFM resid 401" --phenixOptions "main.number_of_macro_cycles=2" --protonation ReadySet  --qmMethod pm6 --region-radius 3.0 --buffer-radius 2.5  --Nproc 4

echo "Tutorial #3 (using qbphenix execution script and MOE on all ligands): 3IX1"
rm -rf ${WORKDIR}/qbphenix_moe_all ; mkdir -p ${WORKDIR}/qbphenix_moe_all ; cd ${WORKDIR}/qbphenix_moe_all

wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/3ix1.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/3ix1-sf.cif

$QBHOME/bin/qbphenix --dataFile 3ix1-sf.cif --pdbFile 3ix1.pdb --selection "resname NFM" --phenixOptions "main.number_of_macro_cycles=1" --protonation MOE  --qmMethod pm6 --region-radius 3.0 --buffer-radius 2.5  --Nproc 4  

echo "Tutorial #4 (using qbphenix execution script and MOE): 1LRI "
rm -rf ${WORKDIR}/qbphenix_moe_1lri ; mkdir -p ${WORKDIR}/qbphenix_moe_1lri ; cd ${WORKDIR}/qbphenix_moe_1lri

wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/1lri.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/1lri-sf.cif

$QBHOME/bin/qbphenix --dataFile 1lri-sf.cif --pdbFile 1lri.pdb --selection "chain A resname CLR resid 99"  --phenixOptions "main.number_of_macro_cycles=2" --protonation MOE --qmMethod pm6 --mmMethod amberff14sb --region-radius 3.0 --buffer-radius 2.5 --Nproc 4 

echo "Tutorial #5: ONIOM (QM/MM) Based X-ray Refinement and Phenix Clash Score Analysis"
rm -rf ${WORKDIR}/Protonation ; mkdir -p ${WORKDIR}/Protonation ; cd ${WORKDIR}/Protonation

$QBHOME/bin/qbphenix --pdbID 1NAV --mmMethod amberff14sb --qmMethod pm6 --selection "resname IH5" --np 4 --region-radius 3.0 --buffer-radius 0.0 --protonation MOE >& OUT.screen

echo "Tutorial #6 (using qbphenix execution script and MOE on covalently bound ligand): 3NCK"
rm -rf ${WORKDIR}/Protonation ; mkdir -p ${WORKDIR}/Protonation ; cd ${WORKDIR}/Protonation

$QBHOME/bin/qbphenix --pdbID 3NCK --selection "resname NFF" --phenixOptions "main.number_of_macro_cycles=2" --protonation MOE  --qmMethod pm6 --mmMethod amberff14sb --region-radius 3.0 --buffer-radius 2.5  --Nproc 4  --scriptName run

grep -v "H1  NFF A" 3NCK.pdb > tmp.$$ ; mv tmp.$$ 3NCK.pdb
grep -v "H3  NFF A" 3NCK.pdb > tmp.$$ ; mv tmp.$$ 3NCK.pdb

./run


currentDate=`date`
echo "END Tutorial Test at ${currentDate} using ${DIVCON_BIN}"
