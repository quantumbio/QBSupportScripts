#!/bin/bash
set -e

#  // BEGIN COPYRIGHT
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
echo "BEGIN DivCon Xray Tutorial Test at ${currentDate} using ${DIVCON_BIN}"

WORKDIR=$PWD

echo "Tutorial #1: All Atom refinement with ONIOM region using protonated PDBid:3e34"
tutorFolder=3e34-All_Atoms
rm -rf ${WORKDIR}/$tutorFolder ; mkdir -p ${WORKDIR}/$tutorFolder ; cd ${WORKDIR}/$tutorFolder

wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/3e34+H.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/3e34.mtz

$QBHOME/bin/qmechanic 3e34+H.pdb 3e34.mtz --opt all 50 0.01 --qm-region /B/ED1/1003// 3.0 0 -h pm6 amberff14sb   --np 4 -O -p 3e34_refined.pdb 3e34_refined.mtz


echo "Tutorial #2: All Atom MM refinement with MTZ labels set by user PDBid:4jp4"
tutorFolder=3e34-All_Atoms
rm -rf ${WORKDIR}/$tutorFolder ; mkdir -p ${WORKDIR}/$tutorFolder ; cd ${WORKDIR}/$tutorFolder

wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4jp4+H.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4jp4.mtz

$QBHOME/bin/qmechanic 4jp4+H.pdb "4jp4.mtz#/*/*[I-obs(+),SIGI-obs(+),I-obs(-),SIGI-obs(-),R-free-flags]" --opt all 25 0.01  -h pm6 amberff14sb   --np 4 -O -p 4jp4_refined.pdb 4jp4_refined.mtz



echo "Tutorial #3: All Atom refinement wtih ONIOM region on the structure downloaded from PDB: PDBid:1lri"
tutorFolder=1lri-allAtom
rm -rf ${WORKDIR}/$tutorFolder ; mkdir -p ${WORKDIR}/$tutorFolder ; cd ${WORKDIR}/$tutorFolder

$QBHOME/bin/qmechanic 1lri --prepare --opt all 50 0.01 --qm-region /A/CLR/99// 3.0 0 -h pm6 amberff14sb   --np 4 -O -p 1lri_refined.pdb 1lri_refined.mtz

echo "Tutorial #4: XModeScore executed on protonated PDBid:4b72"
tutorFolder=xmodeScore_4b72
rm -rf ${WORKDIR}/$tutorFolder ; mkdir -p ${WORKDIR}/$tutorFolder ; cd ${WORKDIR}/$tutorFolder
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4b72+H.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4b72.mtz

$QBHOME/bin/qbphenix --pdbFile 4b72+H.pdb --sfFile 4b72.mtz --XmodeScore --protonation Skip --protomers "-1..1" --engine divcon --protonateTautomers divcon --qmMethod pm6 --mmMethod amberff14sb --nSmallCycles 50 --region-radius 3.0 --buffer-radius 0.0 --selection "resname 2FB and resid 1503 and chain A" --dir testXmode --Nproc 16

echo "Tutorial #5: XModeScore with Dock executed on protonated PDBid:1bzc"
tutorFolder=xmodeScore_1bzc
rm -rf ${WORKDIR}/$tutorFolder ; mkdir -p ${WORKDIR}/$tutorFolder ; cd ${WORKDIR}/$tutorFolder
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/1bzc+H.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/1bzc.mtz

$QBHOME/bin/qbphenix --pdbFile 1bzc+H.pdb --sfFile 1bzc.mtz --XmodeScore --protonation Skip --protomers "0" --exploreFlip --exploreChiral --exploreDocking --engine divcon --protonateTautomers divcon --qmMethod pm6 --mmMethod amberff14sb --nSmallCycles 40 --region-radius 3.0 --selection "resname TPI and resid 902 and chain A" --dir testXmode --Nproc 16

echo "Tutorial #6: Multi-Blob Docking PDBid:4o9s"
tutorFolder=multiBlob_4o9s
rm -rf ${WORKDIR}/$tutorFolder ; mkdir -p ${WORKDIR}/$tutorFolder ; cd ${WORKDIR}/$tutorFolder
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Dock/4o9s+HnL5.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Dock/4o9s_lig.mol2
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Dock/4o9s.mtz
export QB_FIX_METAL_CHARGES=1

$QBHOME/bin/qmechanic 4o9s+HnL5.pdb 4o9s.mtz --ligand blobs 4o9s_lig.mol2 -h amberff14sb --np 2 -v2 --dock opt rigid >test2.log

currentDate=`date`
echo "END Tutorial Test at ${currentDate} using ${DIVCON_BIN}"
