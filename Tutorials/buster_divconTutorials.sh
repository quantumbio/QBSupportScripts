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
region-selection = /A/BSM/1224// 3 0
hamiltonian = pm6 amberff14
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

refine -Gelly 2BSM.qm -p 2BSM+H.pdb -m 2BSM.mtz -l BSM.cif -nbig 2  RunGellySanityCheck=no RunGellyScreen=no -qm_weight 3.0


echo "Tutorial #2 (using qbbuster execution script and Grade & DivCon protonator): 2BSM"
rm -rf ${WORKDIR}/qbbuster_divcon ; mkdir -p ${WORKDIR}/qbbuster_divcon ; cd ${WORKDIR}/qbbuster_divcon

wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/2bsm.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/2bsm-sf.cif

$QBHOME/bin/bin/qbbuster --pdbFile 2bsm.pdb --sfFile 2bsm-sf.cif --protonation divcon --makeCIF grade --mmMethod amberff14 --qmMethod pm6 --qmWeight 5.0 --ncycles 2 --selection "resname BSM" --region-radius 3.0 --np 4 --dir qmRun

echo "Tutorial #3 (qbbuster and XModeScore): 1TOW"
rm -rf ${WORKDIR}/qbbuster_xmodescore ; mkdir -p ${WORKDIR}/qbbuster_xmodescore ; cd ${WORKDIR}/qbbuster_xmodescore

wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/1tow.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/XModeScore/1tow-sf.cif

/path/to/DivConSuite/bin/qbbuster --pdbFile 1tow.pdb --sfFile 1tow-sf.cif --XModeScore --protonation divcon --makeCIF grade --mmMethod amberff14 --qmMethod pm6 --qmWeight 2.0 --selection "resname CRZ chain A resid 501" --region-radius 3.0 --protomers "-1..1" --np 24 --dir XModeScoreRun

currentDate=`date`
echo "END Tutorial Test at ${currentDate} using ${DIVCON_BIN}"
