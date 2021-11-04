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

$QBHOME/bin/qbphenix --pdbID 1OPK --XModeScore --protonation MOE --selection "chain A resname 18K resid 299" --np 20 --dir test1

echo "Tutorial 2a: XModeScore on AZM PDBid:3HS4"
rm -rf ${WORKDIR}/xmodeScore_3HS4 ; mkdir -p ${WORKDIR}/xmodeScore_3HS4 ; cd ${WORKDIR}/xmodeScore_3HS4

$QBHOME/bin/qbphenix --pdbID 3HS4 --XModeScore --protonation MOE --nproc 20 --dir xmodeData --selection "chain A resname AZM resid 701" --protPH 9.0 --generateTautomers divcon --protomers "-1..0"
    
echo "Tutorial 2b: XModeScore on user-provided AZM 3HS4.pdb file"
rm -rf ${WORKDIR}/xmodeScore_user3HS4 ; mkdir -p ${WORKDIR}/xmodeScore_user3HS4 ; cd ${WORKDIR}/xmodeScore_user3HS4

wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/3HS4.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/3HS4.sf

$QBHOME/bin/qbphenix --pdbFile 3HS4.pdb --dataFile  3HS4.mtz --cifFile  3HS4.cif --XModeScore --protonation skip --nproc 20 --dir xmodeData --selection "chain A resname AZM resid 701" --protPH 9.0 --generateTautomers divcon --protomers "-1..0"

echo "Tutorial 2c: XModeScore on pre-generated files"
rm -rf ${WORKDIR}/xmodeScore_pregen3HS4 ; mkdir -p ${WORKDIR}/xmodeScore_pregen3HS4 ; cd ${WORKDIR}/xmodeScore_pregen3HS4

#INSERT EXAMPLE HERE

echo "Tutorial #3: XModeScore on Refmac/Buster structures"
rm -rf ${WORKDIR}/xmodeScore_external1TOW ; mkdir -p ${WORKDIR}/xmodeScore_external1TOW ; cd ${WORKDIR}/xmodeScore_external1TOW

#INSERT EXAMPLE HERE

echo "Tutorial #4: ONIOM (QM/MM) XModeScore on PDBid:3HS4"
rm -rf ${WORKDIR}/xmodeScore_ONIOM_3HS4 ; mkdir -p ${WORKDIR}/xmodeScore_ONIOM_3HS4 ; cd ${WORKDIR}/xmodeScore_ONIOM_3HS4

$QBHOME/bin/qbphenix --pdbid 3HS4 --XmodeScore --protonation moe --nproc 20 --dir xmodeData --selection "chain A resname AZM resid 701" --protPH 9.0 --generateTautomers divcon --protomers "-1,0" --mmMethod amberff14 --qmMethod pm6 --buffer-radius 0.0 --region-radius 3.0
    
echo "Automated MOE-based Docking – using qbphenix"
rm -rf ${WORKDIR}/xmodeScore_MOEDOCK_3HS4 ; mkdir -p ${WORKDIR}/Protonation ; cd ${WORKDIR}/Protonation

# INSERT EXAMPLE HERE


currentDate=`date`
echo "END Tutorial Test at ${currentDate} using ${DIVCON_BIN}"
