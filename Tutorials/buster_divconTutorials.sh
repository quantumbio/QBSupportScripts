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

if [ -z "${BDG_home}" ]; then
    echo "ERROR: BDG_home is not set! You MUST source /path/to/BUSTER_snapshot_XXXXXXXX/setup.sh THEN call this function"
    exit
fi

DIVCON_BIN=${QBHOME}/bin/qmechanic

if [ ! -x "${DIVCON_BIN}" ]; then
    echo "ERROR: ${DIVCON_BIN} is not an executable! Set this path to the qmechanic executable"
    exit
fi

currentDate=`date`
echo "BEGIN BUSTER/DivCon Tutorial Test at ${currentDate} using ${DIVCON_BIN}"

WORKDIR=$PWD

echo "Tutorial #1 (Running BUSTER/DivCon without the qbbuster wrapper): 2BSM"
rm -rf ${WORKDIR}/NakedBUSTER ; mkdir -p ${WORKDIR}/NakedBUSTER ; cd ${WORKDIR}/NakedBUSTER

tee divcon.ini <<EOF
qm-region = /A/BSM/1224// 3 0
hamiltonian = pm6 amberff14sb
np = 4
EOF

tee 2BSM.qm <<EOF
NOTE BUSTER_SET QM01 =  ALL
NOTE BUSTER_QM_CHARGE_01 0
NOTE BUSTER_QM_MULTIP_01 1
NOTE BUSTER_QM_WRESTR_01 0.0
NOTE BUSTER_QM_WRITE_PDB_FOR_HELPER ON
NOTE BUSTER_QM_MAXDISP 1.0E05
NOTE BUSTER_QM_HELPER $QBHOME/scripts/BusterQMHelper.sh
NOTE BUSTER_QM_METHOD PM6
NOTE BUSTER_QM_WEIGHT 3.0
EOF
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/2BSM+H.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/2BSM.mtz
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/BSM.cif

#refine -Gelly 2BSM.qm -p 2BSM+H.pdb -m 2BSM.mtz -l BSM.cif -nbig 2  RunGellySanityCheck=no RunGellyScreen=no -qm_weight 3.0


echo "Tutorial #2 (using qbbuster execution script and Grade & DivCon protonator): 2BSM"
rm -rf ${WORKDIR}/qbbuster_divcon ; mkdir -p ${WORKDIR}/qbbuster_divcon ; cd ${WORKDIR}/qbbuster_divcon

wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/2BSM.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/2BSM.mtz

$QBHOME/bin/qbbuster --pdbFile 2BSM.pdb --sfFile 2BSM.mtz --protonation divcon --makeCIF divcon --mmMethod amberff14sb --qmMethod pm6 --qmWeight 5.0 --ncycles 1  --nSmallCycles 15 --selection "resname BSM" --region-radius 3.0 --np 4 --dir qmRun

echo "Tutorial #3 (running XModeScore with qbbuster execution script on the structure with multiple ligand copies): 4ntk"
rm -rf ${WORKDIR}/qbbuster_4ntk ; mkdir -p ${WORKDIR}/qbbuster_4ntk ; cd ${WORKDIR}/qbbuster_4ntk

wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/4ntk.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/4ntk.mtz

qbbuster --pdbFile 4ntk.pdb --sfFile 4ntk.mtz --XModeScore --protomers "-1" --protonation moe --protonateTautomers MOE --makeCIF moe --mmMethod amberff14sb --qmMethod pm6 --qmWeight 3.0 --engine buster --ncycles 1 --nSmallCycles 20 --selection "resname ZSP and resid 202 and chain A" --np 32 --buffer-radius 0.0 --region-radius 3.0 --dir xmodeScore_results
currentDate=`date`
echo "END Tutorial Test at ${currentDate} using ${DIVCON_BIN}"


echo "Tutorial #4 (running XModeScore with qbbuster execution script on the structure with multiple chiral centers): 4gr0"
rm -rf ${WORKDIR}/qbbuster_4gr0 ; mkdir -p ${WORKDIR}/qbbuster_4gr0 ; cd ${WORKDIR}/qbbuster_4gr0

wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/4gr0+H.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/4gr0.mtz

# C31 chiral atom name in the selection R4B chain A 

qbbuster --pdbFile 4gr0+H.pdb --sfFile 4gr0.mtz --XModeScore --protomers "0" --exploreChiral "C31" --protonation skip --protonateTautomers MOE --makeCIF moe --mmMethod amberff14sb --qmMethod pm6 --qmWeight 3.0 --engine buster --ncycles 1 --nSmallCycles 20 --selection "resname R4B and resid 306 and chain A" --np 32 --region-radius 3.0 --dir xmodeScore_results
currentDate=`date`
echo "END Tutorial Test at ${currentDate} using ${DIVCON_BIN}"

echo "Tutorial #5 (running XModeScore with qbbuster execution script on the structure with multiple flips): 3r4p"
rm -rf ${WORKDIR}/qbbuster_3r4p ; mkdir -p ${WORKDIR}/qbbuster_3r4p ; cd ${WORKDIR}/qbbuster_3r4p

wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/3r4p+H.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/3r4p.mtz

# C31 chiral atom name in the selection R4B chain A 

qbbuster --pdbFile 3r4p+H.pdb --sfFile 3r4p.mtz --XModeScore --protomers "0" --exploreFLip "N6" --protonation skip --protonateTautomers MOE --makeCIF moe --mmMethod amberff14sb --qmMethod pm6 --qmWeight 3.0 --engine buster --ncycles 1 --nSmallCycles 20 --selection "resname FU7 and resid 901 and chain A" --np 32 --region-radius 3.0 --dir xmodeScore_results
currentDate=`date`
echo "END Tutorial Test at ${currentDate} using ${DIVCON_BIN}"