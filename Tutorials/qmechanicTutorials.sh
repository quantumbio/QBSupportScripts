#!/bin/bash
set -e

currentDate=`date`
echo "BEGIN Tutorial Test at ${currentDate}"

#UPDATE DIVCON_BIN to the proper path to the qmechanic binary:
DIVCON_BIN=/home/lance/Release/DivConDiscoverySuite-DEV.895/bin/qmechanic
WORKDIR=$PWD

echo "First Tutorial: Single Point Calculation"
rm -rf ${WORKDIR}/SinglePointCalculation ; mkdir -p ${WORKDIR}/SinglePointCalculation ; cd ${WORKDIR}/SinglePointCalculation

${DIVCON_BIN} 4EK4 --prepare --np 8 -v 2 -p pdb 
mv 4EK4_out.pdb 4EK4.pdb
${DIVCON_BIN} 4EK4.pdb --np 8 -h pm6 -v 2 -p pdb 

echo "Second Tutorial: Interaction Energy Decomposition"
rm -rf ${WORKDIR}/InteractionEnergyDecomposition ; mkdir -p ${WORKDIR}/InteractionEnergyDecomposition ; cd ${WORKDIR}/InteractionEnergyDecomposition

wget http://downloads.quantumbioinc.com/media/tutorials/cli/1LRI-addH.pdb
${DIVCON_BIN} 1LRI-addH.pdb -h pm6 --np 8 -v 2 --decompose "/A/CLR/99//" --dc

echo "Third Tutorial: Structure Optimization"
rm -rf ${WORKDIR}/StructureOptimization ; mkdir -p ${WORKDIR}/StructureOptimization ; cd ${WORKDIR}/StructureOptimization

wget http://downloads.quantumbioinc.com/media/tutorials/cli/1LRI-addH.pdb
${DIVCON_BIN} 1LRI-addH.pdb -h pm6 --np 8 -v 2 --opt all 5 0.01 --symmetry -p pdb

echo "Forth Tutorial: Active Site Structure Optimization"
rm -rf ${WORKDIR}/ActiveSiteStructureOptimization ; mkdir -p ${WORKDIR}/ActiveSiteStructureOptimization ; cd ${WORKDIR}/ActiveSiteStructureOptimization

wget http://downloads.quantumbioinc.com/media/tutorials/cli/1P2Y-addH.pdb
${DIVCON_BIN} 1P2Y-addH.pdb --opt /*/HEM/*//,/*/NCT/*// 5.0 0.0 --np 8 -h pm6 -v 2 -p pdb

echo "Fifth Tutorial: ONIOM (mixed-QM/MM) Simulations"
rm -rf ${WORKDIR}/ONIOM_Simulations ; mkdir -p ${WORKDIR}/ONIOM_Simulations ; cd ${WORKDIR}/ONIOM_Simulations

${DIVCON_BIN} 1LRI --prepare -h pm6 amberff14 --opt 100 0.01 --region-selection /A/CLR/99// 3.0 0.0 --np 8 -v 2 -p pdb

echo "Sixth Tutorial: Protonation"
rm -rf ${WORKDIR}/Protonation ; mkdir -p ${WORKDIR}/Protonation ; cd ${WORKDIR}/Protonation

${DIVCON_BIN} 4EK4 --prepare --np 8 -v 2 -p pdb

currentDate=`date`
echo "END Tutorial Test at ${currentDate}"
