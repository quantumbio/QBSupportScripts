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

echo "Tutorial #1: XModeScore executed on PDBid:1OPK"
rm -rf ${WORKDIR}/xmodeScore_1OPK ; mkdir -p ${WORKDIR}/xmodeScore_1OPK ; cd ${WORKDIR}/xmodeScore_1OPK

$QBHOME/bin/qbphenix --pdbID 4YJR --XModeScore --protomers "-1..1" --mmMethod amberff14sb --qmMethod pm6 --protonation MOE --selection "chain A resname 4DJ resid 701" --np 20 --dir test1
