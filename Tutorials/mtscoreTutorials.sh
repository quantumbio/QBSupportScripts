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

echo "First A Tutorial: Command Line: MTConfSearch"
rm -rf ${WORKDIR}/MTCS-SMILES ; mkdir -p ${WORKDIR}/MTCS-SMILES ; cd ${WORKDIR}/MTCS-SMILES

${DIVCON_BIN} c1ccccc1CCCC --mtcs -p sdf

echo "First B Tutorial: Command Line: MTConfSearch"
rm -rf ${WORKDIR}/MTCS-MOL2 ; mkdir -p ${WORKDIR}/MTCS-MOL2 ; cd ${WORKDIR}/MTCS-MOL2


echo "Second Tutorial: Command Line: MTScore (Endstate)"
rm -rf ${WORKDIR}/MTScoreES ; mkdir -p ${WORKDIR}/MTScoreES ; cd ${WORKDIR}/MTScoreES

wget http://downloads.quantumbioinc.com/media/tutorials/MT/4w7t_protein.pdb
wget http://downloads.quantumbioinc.com/media/tutorials/MT/4w7t_ligand.mol2

${DIVCON_BIN} 4w7t_protein.pdb --ligand 4w7t_ligand.mol2 --mtscore endstate

echo "Third Tutorial: 3Step MTCS+MTScoreE (Ensemble) with an External Docker"
rm -rf ${WORKDIR}/3step ; mkdir -p ${WORKDIR}/3step ; cd ${WORKDIR}/3step

wget http://downloads.quantumbioinc.com/media/tutorials/MT/Bace_030215_CAT_4p.tar.gz
tar xvf Bace_030215_CAT_4p.tar.gz
cp Bace_030215_CAT_4p/Bace_030215.pdb .
cp Bace_030215_CAT_4p/Bace_030215_CAT_4p.mol2 .

${DIVCON_BIN} Bace_030215.pdb --ligand Bace_030215_CAT_4p.mol2 --mtcs -p sdf --np 4 -v 2

moebatch -licwait -run "${QBHOME}/svl/run/qbDockPair.svl" -rec Bace_030215.pdb -lig Bace_030215_CAT_4p.mol2 -conf Bace_030215_CAT_4p_conf.sdf -protonate -delwat

${DIVCON_BIN} Bace_030215.pdb --ligand Bace_030215_CAT_4p.mol2 -h garf --mtdock Bace_030215_CAT_4p_dock.sdf --mtscore ensemble --np 4 -v 2

echo "Forth Tutorial: 2Step MTScoreE with an External Docker (no-MTCS)"
rm -rf ${WORKDIR}/2Step ; mkdir -p ${WORKDIR}/2Step ; cd ${WORKDIR}/2Step

wget http://downloads.quantumbioinc.com/media/tutorials/MT/Bace_030215_CAT_4p.tar.gz
tar xvf Bace_030215_CAT_4p.tar.gz
cp Bace_030215_CAT_4p/Bace_030215.pdb .
cp Bace_030215_CAT_4p/Bace_030215_CAT_4p.mol2 .

moebatch -licwait -run "${QBHOME}/svl/run/qbDockPair.svl" -rec Bace_030215.pdb -lig Bace_030215_CAT_4p.mol2 -delwat -protonate

${DIVCON_BIN} pro_Bace_030215_CAT_4p_predock.pdb --ligand lig_Bace_030215_CAT_4p_predock.mol2 -h garf --mtdock Bace_030215_CAT_4p_dock.sdf --mtscore ensemble --np 4 -v 2


echo "Fifth Tutorial: MTScoreE with induced fit protein:ligand docking"
rm -rf ${WORKDIR}/InducedFit ; mkdir -p ${WORKDIR}/InducedFit ; cd ${WORKDIR}/InducedFit

wget http://downloads.quantumbioinc.com/media/tutorials/MT/Bace_030215_CAT_4p.tar.gz
tar xvf Bace_030215_CAT_4p.tar.gz
cp Bace_030215_CAT_4p/Bace_030215.pdb .
cp Bace_030215_CAT_4p/Bace_030215_CAT_4p.mol2 .

moebatch -licwait -run "${QBHOME}/svl/run/qbDockPair.svl" -rec Bace_030215.pdb -lig Bace_030215_CAT_4p.mol2 -inducedfit -delwat -protonate

${DIVCON_BIN} pro_Bace_030215_CAT_4p_predock.pdb --ligand lig_Bace_030215_CAT_4p_predock.mol2 -h amberff14sb --mtdock *_dock.pdb --mtscore ensemble --np 4 -v 2

echo "Sixth Tutorial: Manipulating cutoffs to maximize predictions"
rm -rf ${WORKDIR}/Cutoffs ; mkdir -p ${WORKDIR}/Cutoffs ; cd ${WORKDIR}/Cutoffs

wget http://downloads.quantumbioinc.com/media/tutorials/MT/Bace_030215_CAT_4p.tar.gz
tar xvf Bace_030215_CAT_4p.tar.gz
cp Bace_030215_CAT_4p/Bace_030215.pdb .
cp Bace_030215_CAT_4p/Bace_030215_CAT_4p.mol2 .

moebatch -licwait -run "${QBHOME}/svl/run/qbDockPair.svl" -rec Bace_030215.pdb -lig Bace_030215_CAT_4p.mol2 -delwat -protonate

${DIVCON_BIN} pro_Bace_030215_CAT_4p_predock.pdb --ligand lig_Bace_030215_CAT_4p_predock.mol2 -h garf --mtdock Bace_030215_CAT_4p_dock.sdf --mtscore ensemble --nb-cutoff 8.0  --np 4 -v 2

echo "Sixth Tutorial: Manipulating cutoffs to maximize predictions"



currentDate=`date`
echo "END Tutorial Test at ${currentDate} using ${DIVCON_BIN}"
