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


DIVCON_BIN=${QBHOME}/bin/qmechanic

if [ ! -x "${DIVCON_BIN}" ]; then
    echo "ERROR: ${DIVCON_BIN} is not an executable! Set this path to the qmechanic executable"
    exit
fi

currentDate=`date`
echo "BEGIN DivCon Xrya Tutorial Test at ${currentDate} using ${DIVCON_BIN}"

WORKDIR=$PWD

echo "Tutorial #1: All Atom refinement with ONIOM region using protonated PDBid:3e34"
tutorFolder=3e34-All_Atoms
rm -rf ${WORKDIR}/$tutorFolder ; mkdir -p ${WORKDIR}/$tutorFolder ; cd ${WORKDIR}/$tutorFolder

wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/3e34+H.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/3e34.mtz

$QBHOME/bin/qmechanic 3e34+H.pdb 3e34.mtz --opt all 50 0.01 --qm-region /B/ED1/1003// 3.0 0 -h pm6 amberff14sb   --np 4 -O -p 3e34_refined.pdb 3e34_refined.mtz

echo "Tutorial #2: All Atom refinement wtih ONIOM region on the structure downloaded from PDB: PDBid:1lri"
tutorFolder=1lri-allAtom
rm -rf ${WORKDIR}/$tutorFolder ; mkdir -p ${WORKDIR}/$tutorFolder ; cd ${WORKDIR}/$tutorFolder
export DEV_ENABLE_FIRST_ALTERNATE=1

$QBHOME/bin/qmechanic 1lri --prepare --opt all 50 0.01 --qm-region /A/CLR/99// 3.0 0 -h pm6 amberff14sb   --np 4 -O -p 1lri_refined.pdb 1lri_refined.mtz

echo "Tutorial #3: XModeScore executed on protonated PDBid:4b72"
tutorFolder=xmodeScore_4b72
rm -rf ${WORKDIR}/$tutorFolder ; mkdir -p ${WORKDIR}/$tutorFolder ; cd ${WORKDIR}/$tutorFolder
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4b72+H.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4b72.mtz

$QBHOME/bin/qbphenix --pdbFile 4b72+H.pdb --sfFile 4b72.mtz --XmodeScore --protonation Skip --protomers "-1..1" --engine divcon --protonateTautomers divcon --qmMethod pm6 --mmMethod amberff14sb --nSmallCycles 50 --region-radius 3.0 --buffer-radius 0.0 --selection "resname 2FB and resid 1503 and chain A" --dir testXmode --Nproc 16

echo "Tutorial #4: Cryo_EM Refinement on protonated PDBid:7jsy"
tutorFolder=cryoEM_7jsy
rm -rf ${WORKDIR}/$tutorFolder ; mkdir -p ${WORKDIR}/$tutorFolder ; cd ${WORKDIR}/$tutorFolder
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/7jsy+H.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/emd_22463.map

$QBHOME/bin/qmechanic 7jsy+H.pdb emd_22463.map --opt all 35 0.01 --qm-region /A/I3C/501// 0.0 0.0 -h pm6 amberff14sb -O -p 7jsy+H_refined.pdb 7jsy+H_refined.mtz --np 4 -v 2 --nb-cutoff 25


currentDate=`date`
echo "END Tutorial Test at ${currentDate} using ${DIVCON_BIN}"