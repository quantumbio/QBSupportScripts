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
echo "BEGIN DivCon CryoEM Tutorial Test at ${currentDate} using ${DIVCON_BIN}"

WORKDIR=$PWD

echo "Tutorial #1: Cryo_EM Refinement on protonated PDBid:7jsy"
tutorFolder=cryoEM_7jsy
rm -rf ${WORKDIR}/$tutorFolder ; mkdir -p ${WORKDIR}/$tutorFolder ; cd ${WORKDIR}/$tutorFolder
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/7jsy+H.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/emd_22463.map

$QBHOME/bin/qmechanic 7jsy+H.pdb emd_22463.map --resolution 1.8 --experiment cryoEM --opt all 35 0.01 --qm-region /A/I3C/501// 0.0 0.0 -h pm6 amberff14sb -O -p 7jsy+H_refined.pdb 7jsy+H_refined.mtz --np 4 -v 2 --nb-cutoff 25

#
echo "Tutorial #2: Cryo_EM ONIOM Refinement with qbdivcon: PDBid:7efc"
tutorFolder=cryoEM_qbdivcon_7efc
rm -rf ${WORKDIR}/$tutorFolder ; mkdir -p ${WORKDIR}/$tutorFolder ; cd ${WORKDIR}/$tutorFolder

wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/7efc+H.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/emd_31083.map

$QBHOME/bin/qbdivcon --pdbfile 7efc+H.pdb --sfFile emd_31083.map --experiment cryoEM --resolution 1.7 --protonation skip  --engine divcon --qmMethod pm6 --mmMethod amberff14sb --resname BTN --chain A --resid 5100 --np 4 --region-radius 3.0 --nSmallCycles 40 
#

echo "Tutorial #3: Cryo_EM XModeScore with qbdivcon: PDBid:7jsy"
tutorFolder=cryoEM_xmodescore_7jsy
rm -rf ${WORKDIR}/$tutorFolder ; mkdir -p ${WORKDIR}/$tutorFolder ; cd ${WORKDIR}/$tutorFolder

wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/7jsy+H.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/emd_22463.map

$QBHOME/bin/qbdivcon --pdbfile 7jsy+H.pdb --sfFile emd_22463.map --experiment cryoEM --resolution 1.8 --XmodeScore --protomers "-1..1" --exploreFlip --protonation skip  --engine divcon --qmMethod pm6 --mmMethod amberff14sb --resname I3C --chain A --resid 501 --np 16 --region-radius 3.0 --nSmallCycles 40 --dir testXmode
#

currentDate=`date`
echo "END Tutorial Test at ${currentDate} using ${DIVCON_BIN}"

