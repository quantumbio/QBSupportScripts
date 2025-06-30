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

function tutorial1_action {
echo "First Tutorial: Single Point Calculation"
rm -rf ${WORKDIR}/SinglePointCalculation ; mkdir -p ${WORKDIR}/SinglePointCalculation ; cd ${WORKDIR}/SinglePointCalculation

wget http://downloads.quantumbioinc.com/media/tutorials/cli/1TOW-H.pdb
mv  1TOW-H.pdb 1TOW.pdb
${DIVCON_BIN} 1TOW.pdb --np 8 -h pm6 -v 2 
}

function tutorial2_action {
echo "Second Tutorial: Interaction Energy Decomposition"
rm -rf ${WORKDIR}/InteractionEnergyDecomposition ; mkdir -p ${WORKDIR}/InteractionEnergyDecomposition ; cd ${WORKDIR}/InteractionEnergyDecomposition

wget http://downloads.quantumbioinc.com/media/tutorials/cli/1LRI-addH.pdb
${DIVCON_BIN} 1LRI-addH.pdb -h pm6 --np 8 -v 2 --decompose "/A/CLR/99//" --dc
}

function tutorial3_action {
echo "Third Tutorial: Structure Optimization"
rm -rf ${WORKDIR}/StructureOptimization ; mkdir -p ${WORKDIR}/StructureOptimization ; cd ${WORKDIR}/StructureOptimization

wget http://downloads.quantumbioinc.com/media/tutorials/cli/1LRI-addH.pdb
${DIVCON_BIN} 1LRI-addH.pdb -h pm6 --np 8 -v 2 --opt all 3 0.01 --symmetry -p pdb
}

function tutorial4_action {
echo "Forth Tutorial: Active Site Structure Optimization"
rm -rf ${WORKDIR}/ActiveSiteStructureOptimization ; mkdir -p ${WORKDIR}/ActiveSiteStructureOptimization ; cd ${WORKDIR}/ActiveSiteStructureOptimization

wget http://downloads.quantumbioinc.com/media/tutorials/cli/1TOW-H.pdb
${DIVCON_BIN} 1TOW-H.pdb --opt /*/CRZ/*// 3.0 0.0 --np 8 -h pm6 -v 2 -p pdb
}

function tutorial5_action {
echo "Fifth Tutorial: ONIOM (mixed-QM/MM) Simulations"
rm -rf ${WORKDIR}/ONIOM_Simulations ; mkdir -p ${WORKDIR}/ONIOM_Simulations ; cd ${WORKDIR}/ONIOM_Simulations

wget http://downloads.quantumbioinc.com/media/tutorials/cli/1LRI-addH.pdb
${DIVCON_BIN} 1LRI-addH.pdb -h pm6 amberff14sb --opt 25 0.01 --qm-region /A/CLR/99// 3.0 0.0 --np 8 -v 1 -p pdb
}

function tutorial6_action {
echo "Sixth Tutorial: Protonation"
rm -rf ${WORKDIR}/Protonation ; mkdir -p ${WORKDIR}/Protonation ; cd ${WORKDIR}/Protonation

${DIVCON_BIN} 4EK4 --prepare --np 8 -v 2 -p pdb -h amberff14sb

currentDate=`date`
echo "END Tutorial Test at ${currentDate} using ${DIVCON_BIN}"
}

function tutorial7_action {
echo "Seventh Tutorial: Protonation"
rm -rf ${WORKDIR}/macrocycle ; mkdir -p ${WORKDIR}/macrocycle ; cd ${WORKDIR}/macrocycle

cd macrocycle
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/cli/3p72_macrocycle_bond.json
${DIVCON_BIN} 3p72  --standards 3p72_macrocycle_bond.json -h amberff14sb --prepare all --np 4 -v2 -p 3p72+H.pdb 3p72.mtz -O 2>&1 
cd ..
grep 
currentDate=`date`
echo "END Tutorial Test at ${currentDate} using ${DIVCON_BIN}"
}

function tutorial8_action {
echo "Eighth Tutorial: Covalently Bound Ligand"
rm -rf ${WORKDIR}/bound_ligand ; mkdir -p ${WORKDIR}/bound_ligand ; mkdir bound_ligand
cd bound_ligand
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/cli/5y41+H.pdb
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/cli/5y41_covalent_bond.json
${DIVCON_BIN} 5y41+H.pdb  --standards 5y41_covalent_bond.json -h amberff14sb --np 4 -v2  -O 2>&1 
cd ..
grep 
currentDate=`date`

echo "END Tutorial Test at ${currentDate} using ${DIVCON_BIN}"
}

if [ "$1" == "menu" ]; then
PS3="Please choose an option: " # Prompt displayed to the user
options=("Tut 1" "Tut 2" "Tut 3" "Tut 4" "Tut 5" "Tut 6" "Tut 7" "Tut 8" "Exit")

select opt in "${options[@]}"; do
    case $opt in
        "Tut 1")
            tutorial1_action
            ;;
        "Tut 2")
            tutorial2_action
            ;;
        "Tut 3")
            tutorial3_action
            ;;
        "Tut 4")
            tutorial4_action
            ;;
        "Tut 5")
            tutorial5_action
            ;;
        "Tut 6")
            tutorial6_action
            ;;
        "Tut 7")
            tutorial7_action
            ;;
        "Tut 8")
            tutorial8_action
            ;;
        "Exit")
            echo "Exiting menu."
            break # Exit the select loop
            ;;
        *)
            echo "Invalid option. Please try again."
            ;;
    esac
done
else
    tutorial1_action
    tutorial2_action
    tutorial3_action
    tutorial4_action
    tutorial5_action
    tutorial6_action
    tutorial7_action
    tutorial8_action
fi
