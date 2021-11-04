#!/bin/bash
#set -e

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
    echo "ERROR: QBHOME is not set! Must set the QBHOME environment variable THEN call this function"
    exit
fi

DIVCON_BIN=${QBHOME}/bin/qmechanic

if [ ! -x "${DIVCON_BIN}" ]; then
    echo "ERROR: ${DIVCON_BIN} is not an executable! Set this path to the qmechanic executable"
    exit
fi

currentDate=`date`
echo "BEGIN Tutorial Test at ${currentDate} using ${DIVCON_BIN}"

WORKDIR=$PWD

echo "First Tutorial: Single Point Calculation"
rm -rf ${WORKDIR}/SinglePointCalculation ; mkdir -p ${WORKDIR}/SinglePointCalculation ; cd ${WORKDIR}/SinglePointCalculation

wget http://downloads.quantumbioinc.com/media/tutorials/cli/4EK4_out.pdb
mv 4EK4_out.pdb 4EK4.pdb
${DIVCON_BIN} 4EK4.pdb --np 8 -h pm6 -v 2 -p pdb 

echo "Second Tutorial: Interaction Energy Decomposition"
rm -rf ${WORKDIR}/InteractionEnergyDecomposition ; mkdir -p ${WORKDIR}/InteractionEnergyDecomposition ; cd ${WORKDIR}/InteractionEnergyDecomposition

wget http://downloads.quantumbioinc.com/media/tutorials/cli/1LRI-addH.pdb
${DIVCON_BIN} 1LRI-addH.pdb -h pm6 --np 8 -v 2 --decompose "/A/CLR/99//" --dc

echo "Third Tutorial: Structure Optimization"
rm -rf ${WORKDIR}/StructureOptimization ; mkdir -p ${WORKDIR}/StructureOptimization ; cd ${WORKDIR}/StructureOptimization

wget http://downloads.quantumbioinc.com/media/tutorials/cli/1LRI-addH.pdb
${DIVCON_BIN} 1LRI-addH.pdb -h pm6 --np 8 -v 2 --opt all 3 0.01 --symmetry -p pdb

echo "Forth Tutorial: Active Site Structure Optimization"
rm -rf ${WORKDIR}/ActiveSiteStructureOptimization ; mkdir -p ${WORKDIR}/ActiveSiteStructureOptimization ; cd ${WORKDIR}/ActiveSiteStructureOptimization

wget http://downloads.quantumbioinc.com/media/tutorials/cli/1P2Y-addH.pdb
${DIVCON_BIN} 1P2Y-addH.pdb --opt /*/HEM/*//,/*/NCT/*// 5.0 0.0 --np 8 -h pm6 -v 2 -p pdb

echo "Fifth Tutorial: ONIOM (mixed-QM/MM) Simulations"
rm -rf ${WORKDIR}/ONIOM_Simulations ; mkdir -p ${WORKDIR}/ONIOM_Simulations ; cd ${WORKDIR}/ONIOM_Simulations

wget http://downloads.quantumbioinc.com/media/tutorials/cli/1LRI-addH.pdb
${DIVCON_BIN} 1LRI-addH.pdb -h pm6 amberff14sb --opt 100 0.01 --qm-region /A/CLR/99// 3.0 0.0 --np 8 -v 2 -p pdb

echo "Sixth Tutorial: Protonation"
rm -rf ${WORKDIR}/Protonation ; mkdir -p ${WORKDIR}/Protonation ; cd ${WORKDIR}/Protonation

${DIVCON_BIN} 4EK4 --prepare --np 8 -v 2 -p pdb

currentDate=`date`
echo "END Tutorial Test at ${currentDate} using ${DIVCON_BIN}"
