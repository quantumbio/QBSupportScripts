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

echo "Tutorial #1: XModeScore executed on PDBid:1OPK"
rm -rf ${WORKDIR}/xmodeScore_1OPK ; mkdir -p ${WORKDIR}/xmodeScore_1OPK ; cd ${WORKDIR}/xmodeScore_1OPK

$QBHOME/bin/qbphenix --pdbID 4YJR --XModeScore --protomers "-1..1" --mmMethod amberff14sb --qmMethod pm6 --protonation MOE --selection "chain A resname 4DJ resid 701" --np 20 --dir test1

echo "Tutorial 2a: XModeScore on user-provided AZM 3HS4.pdb file"
rm -rf ${WORKDIR}/xmodeScore_3HS4 ; mkdir -p ${WORKDIR}/xmodeScore_3HS4 ; cd ${WORKDIR}/xmodeScore_3HS4
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/3HS4-tautomers/3HS4+H.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/3HS4-tautomers/3HS4.mtz
$QBHOME/bin/qbphenix --pdbFile 3HS4+H.pdb --dataFile  3HS4.mtz --XModeScore --protomers "-1..0" --mmMethod amberff14sb --qmMethod pm6 --protonation skip --nproc 20 --dir xmodeData --selection "chain A resname AZM resid 701" --protPH 9.0 --generateTautomers divcon 
    
echo "Tutorial 2b: XModeScore on pre-generated files"
rm -rf ${WORKDIR}/xmodeScore_pregen3HS4 ; mkdir -p ${WORKDIR}/xmodeScore_pregen3HS4 ; cd ${WORKDIR}/xmodeScore_pregen3HS4
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/3HS4-tautomers/3HS4+H.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/3HS4-tautomers/3HS4.mtz
mkdir 3HS4-tautomers; cd 3HS4-tautomers
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/3HS4-tautomers/0_0_0_0_-1.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/3HS4-tautomers/1_1_0_0_-1.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/3HS4-tautomers/3_1_0_2_-1.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/3HS4-tautomers/4_0_0_0_0.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/3HS4-tautomers/6_1_0_1_0.pdb
cd ../
$QBHOME/bin/qbphenix --pdbFile 3HS4+H.pdb --dataFile  3HS4.mtz --XModeScore 3HS4-tautomers --selection "chain A resname AZM resid 701" --mmMethod amberff14sb --qmMethod pm6 --protonation skip --protonateTautomers skip --nproc 20 --dir xmodeData --v 1 

echo "Tutorial #3: XModeScore on Buster structures"
rm -rf ${WORKDIR}/xmodeScore_external_2BSM ; mkdir -p ${WORKDIR}/xmodeScore_external_2BSM ; cd ${WORKDIR}/xmodeScore_external_2BSM
source /share/apps/GlobalPhasing/linux-x86_64/BUSTER_snapshot_20180515/setup.sh
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/2BSM+H.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/2BSM.mtz
qbbuster --pdbFile 2BSM+H.pdb --sfFile 2BSM.mtz --XModeScore --protomers "-1..1" --protonation skip --makeCIF grade --mmMethod amberff14sb --qmMethod pm6 --qmWeight 5.0 --ncycles 1 --resname BSM --chain A --resid 1224 --np 20 --buffer-radius 0.0 --region-radius 3.0

echo "Tutorial #4: Automated MOE-based Docking – using qbphenix"
dir4=xmodeScore_autoDock
rm -rf ${WORKDIR}/$dir4 ; mkdir -p ${WORKDIR}/$dir4 ; cd ${WORKDIR}/$dir4

$QBHOME/bin/qbphenix --pdbFile 5C3K-H_refine_001.pdb --dataFile 5C3K.mtz --densityFile 5C3K-H_refine_001.mtz --ligandFile 4XF_A_402.pdb --dock dockFolderResults --protonation skip --Xmodescore --protonateTautomers MOE --mmMethod amberff14 --qmMethod pm6 --region-radius 3.0 --nproc 22
    
echo "Tutorial #5: Real-Space Ligand Statistics – using qmechanic"
dir5=RealSpaceStatictics
rm -rf ${WORKDIR}/$dir5 ; mkdir -p ${WORKDIR}/$dir5 ; cd ${WORKDIR}/$dir5
wget http://downloads.quantumbioinc.com/media/tutorials/RealSpaceStats/0_0_0_0_-1_refine_001.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/RealSpaceStats/0_0_0_0_-1_refine_001.mtz
$QBHOME/bin/qmechanic 0_0_0_0_-1_refine_001.pdb 0_0_0_0_-1_refine_001.mtz --xstats "/A/AZM/701//,/A/HIS/96//"


currentDate=`date`
echo "END Tutorial Test at ${currentDate} using ${DIVCON_BIN}"
