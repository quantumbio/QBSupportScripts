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
echo "BEGIN Workshop 2003 Tutorials at ${currentDate} using ${DIVCON_BIN}"

qbExec=qbphenix
cloud=""
WORKDIR=$PWD
if [[ -v GRID_MARKET ]];  
then 
qbExec=qbdivcon
cloud=" --cloud gridmarkets"
else
ENGINE_DIVCON=""
fi

echo "Tutorial #1:ONIOM  Refinement of PDBid:1LRI with Protonation of the input PDB"
dir1=Tutorial1_1LRI
rm -rf ${WORKDIR}/$dir1 ; mkdir -p ${WORKDIR}/$dir1 ; cd ${WORKDIR}/$dir1
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/1LRI.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/1LRI.mtz

$QBHOME/bin/$qbExec --pdbfile 1LRI.pdb --sfFile 1LRI.mtz --protonation divcon  --engine divcon --qmMethod pm6 --mmMethod amberff14sb --resname CLR --np 4 --region-radius 3.0 --nSmallCycles 25 --ncycles 2 --dir qm_results $cloud
#
echo "Tutorial #1a:ONIOM PHENIX Refinement of PDBid:1LRI"
dir1=Tutorial1_1LRI_PHENIX
rm -rf ${WORKDIR}/$dir1 ; mkdir -p ${WORKDIR}/$dir1 ; cd ${WORKDIR}/$dir1
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/1LRI.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/1LRI.mtz

$QBHOME/bin/qbphenix --pdbfile 1LRI.pdb --sfFile 1LRI.mtz --protonation divcon --makeCIF divcon --qmMethod pm6 --mmMethod amberff14sb --resname CLR --np 4 --region-radius 3.0 --nSmallCycles 25 --ncycles 2 --dir qm_results
#
echo "Tutorial #2:ONIOM  Refinement of Structure with FAD: PDBid:1SIQ"
dir1=Tutorial2_1SIQ
rm -rf ${WORKDIR}/$dir1 ; mkdir -p ${WORKDIR}/$dir1 ; cd ${WORKDIR}/$dir1
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/1SIQ+H.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/1SIQ.mtz

$QBHOME/bin/$qbExec --pdbfile 1SIQ+H.pdb --sfFile 1SIQ.mtz --protonation skip  --engine divcon --qmMethod pm6 --mmMethod amberff14sb --resname FAD --chain A --resid 399 --np 4 --region-radius 3.0 --nSmallCycles 40 --ncycles 2 --dir qm_results $cloud
#
echo "Tutorial #3:ONIOM  XModeScore on AZM : PDBid:3HS4"
dir1=Tutorial3_AZM_XmodeScore
rm -rf ${WORKDIR}/$dir1 ; mkdir -p ${WORKDIR}/$dir1 ; cd ${WORKDIR}/$dir1
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4+H.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4.mtz

$QBHOME/bin/$qbExec --pdbfile 3HS4+H.pdb --sfFile 3HS4.mtz --Xmodescore --protomers -1..0  --protonation skip  --engine divcon --qmMethod pm6 --mmMethod amberff14sb --resname AZM --chain A --resid 701 --np 10 --region-radius 3.0 --nSmallCycles 35 --dir XmodeScore_results $cloud
#



echo "Tutorial 4: XModeScore on user-provided AZM tautomer files "
rm -rf ${WORKDIR}/xmodeScore_3HS4 ; mkdir -p ${WORKDIR}/xmodeScore_3HS4 ; cd ${WORKDIR}/xmodeScore_3HS4
dir1=Tutorial4_AZM_XmodeScore_ExternalFiles
rm -rf ${WORKDIR}/$dir1 ; mkdir -p ${WORKDIR}/$dir1 ; cd ${WORKDIR}/$dir1
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4+H.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4.mtz

mkdir 3HS4-tautomers; cd 3HS4-tautomers
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4+H_0_0_0_0_-1_C1_refined.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4+H_1_1_0_0_-1_C1_refined.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4+H_3_1_0_2_-1_C1_refined.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4+H_4_0_0_0_0_C1_refined.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4+H_6_1_0_1_0_C1_refined.pdb
cd ../
$QBHOME/bin/$qbExec --pdbFile 3HS4+H.pdb --dataFile  3HS4.mtz --XModeScore 3HS4-tautomers --protonation skip  --protonateTautomers skip --engine divcon --mmMethod amberff14sb --qmMethod pm6 --selection "chain A resname AZM resid 701" --nproc 10 --region-radius 3.0 --nSmallCycles 35 --dir XmodeScore_results $cloud 
#
echo "Tutorial #5: Automated MOE-based Docking – using qbphenix: REQUIRES MOE "
dir4=Tutorial5_MOE_AutoDock_Xmode
rm -rf ${WORKDIR}/$dir4 ; mkdir -p ${WORKDIR}/$dir4 ; cd ${WORKDIR}/$dir4
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K.mtz
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.mtz
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4XF_A_402.pdb

$QBHOME/bin/qbphenix --pdbFile 5C3K-H_refine_001.pdb --dataFile 5C3K.mtz --densityFile 5C3K-H_refine_001.mtz --ligandFile 4XF_A_402.pdb --dock dockFolderResults --protonation skip --Xmodescore --protonateTautomers MOE --mmMethod amberff14sb --qmMethod pm6 --region-radius 3.0 --np 22 --dir Dock_Scored_results
    
echo "Tutorial #6: XModeScore with the new features to explore flip and chiral states "
dir3=Tutorial6_xmodeScore_4wq6
rm -rf ${WORKDIR}/$dir3 ; mkdir -p ${WORKDIR}/$dir3 ; cd ${WORKDIR}/$dir3
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4wq6+H.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4wq6.mtz
$QBHOME/bin/$qbExec --pdbFile 4wq6+H.pdb --sfFile 4wq6.mtz --XModeScore --protomers "-1..1" --exploreFlip --exploreChiral --protonation skip  --mmMethod amberff14sb --qmMethod pm6 --engine divcon --nSmallCycles 40 --resname 3TQ --chain A --resid 601 --np 12 --region-radius 3.0 --dir XmodeScore_results $cloud


currentDate=`date`
echo "END Tutorial Test at ${currentDate} using ${DIVCON_BIN}"
