#!/bin/bash -e

export MT_ALLOW_MTDOCK=1 

# This script will test all of the various supported command line options available in the 
# EnumerationDriver/MTCS -> Dock/MTDock -> XModeScore/MTScore protocol

# WARNING: currently only support for options is being tested. Proper functional tests should be designed to test underlying features

# First, create the work environment for this test. We'll use uuidgen and the first hexadecimal group.

jobid=`uuidgen | awk -F "-" '{print $1}'`

if [ -z ${QBHOME+x} ]; then 
    echo "Set export QBHOME=/path/to/DivConSuite_to_test and re-run"
    exit
else 
    echo "Testing DivCon Suite available at ${QBHOME} in directory $jobid"
fi

mkdir ${jobid} ; cd ${jobid}
WORKDIR="${PWD}"

wget https://github.com/quantumbio/QBSupportScripts/raw/master/Tutorials/data/MT/4w7t_ligand.mol2
wget https://github.com/quantumbio/QBSupportScripts/raw/master/Tutorials/data/MT/4w7t_protein.pdb
wget https://github.com/quantumbio/QBSupportScripts/raw/master/Tutorials/data/MT/4w7t-sf.cif
wget https://github.com/quantumbio/QBSupportScripts/raw/master/Tutorials/data/MT/4w7t_phases.mtz

wget https://github.com/quantumbio/QBSupportScripts/raw/master/Tutorials/data/XModeScore/1bzc%2BH.pdb
wget https://github.com/quantumbio/QBSupportScripts/raw/master/Tutorials/data/XModeScore/1bzc.mtz
wget https://github.com/quantumbio/QBSupportScripts/raw/master/Tutorials/data/XModeScore/1bzc-sf.cif


echo "Default MTCS+MTDock"
mkdir -p "${WORKDIR}/default-MTCS_MTDOCK_MTSCOREE"
cd "${WORKDIR}/default-MTCS_MTDOCK_MTSCOREE"
ln -s ../4w7t_ligand.mol2
ln -s ../4w7t_protein.pdb
CLI_ARGS="4w7t_protein.pdb        \
    --ligand 4w7t_ligand.mol2                   \
    --mtcs 2 preopt opt torsion 100 1e-3      \
    --mtdock 2,3 opt torsion 100 0.01  \
    --mtscore ensemble                          \
    -p pdb -v 2 -O -h amberff14sb --np 10"
#    --mtdock 5,10 opt torsion ligand 100 0.01  \

echo ${QBHOME}/bin/qmechanic ${CLI_ARGS}
${QBHOME}/bin/qmechanic ${CLI_ARGS}

echo "Default MTCS+MTDock+Density (MTZ) test"
mkdir -p "${WORKDIR}/default-MTCS_MTDOCK_MTSCOREE-mtz"
cd "${WORKDIR}/default-MTCS_MTDOCK_MTSCOREE-mtz"
ln -s ../4w7t_ligand.mol2
ln -s ../4w7t_protein.pdb
ln -s ../4w7t_phases.mtz
CLI_ARGS="4w7t_protein.pdb 4w7t_phases.mtz       \
    --ligand 4w7t_ligand.mol2                   \
    --mtcs 2 preopt opt torsion 100 1e-3      \
    --mtdock 2,3 opt torsion 100 0.01  \
    --mtscore ensemble                          \
    -p pdb -v 2 -O -h amberff14sb --np 10"
#    --mtdock 2,3 opt torsion ligand 100 0.01  \

echo ${QBHOME}/bin/qmechanic ${CLI_ARGS}
${QBHOME}/bin/qmechanic ${CLI_ARGS}

echo "Default MTCS+MTDock+Density (SF-CIF) test"
mkdir -p "${WORKDIR}/default-MTCS_MTDOCK_MTSCOREE-sf"
cd "${WORKDIR}/default-MTCS_MTDOCK_MTSCOREE-sf"
ln -s ../4w7t_ligand.mol2
ln -s ../4w7t_protein.pdb
ln -s ../4w7t-sf.cif
CLI_ARGS="4w7t_protein.pdb 4w7t-sf.cif       \
    --ligand 4w7t_ligand.mol2                   \
    --mtcs 2 preopt opt torsion 100 1e-3      \
    --mtdock 2,3 opt torsion 100 0.01  \
    --mtscore ensemble                          \
    -p pdb -v 2 -O -h amberff14sb --np 10"
#    --mtdock 5,10 opt torsion ligand 100 0.01  \

echo ${QBHOME}/bin/qmechanic ${CLI_ARGS}
${QBHOME}/bin/qmechanic ${CLI_ARGS}

echo "Default ED+MTDock+Density (MTZ) test"
mkdir -p "${WORKDIR}/default-ED_MTDOCK_MTSCOREE-mtz"
cd "${WORKDIR}/default-ED_MTDOCK_MTSCOREE-mtz"
ln -s ../1bzc+H.pdb
ln -s ../1bzc.mtz
CLI_ARGS="1bzc+H.pdb 1bzc.mtz      \
    --structures [-1..1] /A/TPI/902// \
   --mtdock opt torsion 15 0.01  \
    -p pdb -v 2 -O -h amberff14sb --np 10"
#    --mtdock 2,3 opt torsion pocket 100 0.01  \

echo ${QBHOME}/bin/qmechanic ${CLI_ARGS}
${QBHOME}/bin/qmechanic ${CLI_ARGS}

echo "Default ED+MTDock+Density (SF-CIF) test"
mkdir -p "${WORKDIR}/default-ED_MTDOCK_MTSCOREE-sf"
cd "${WORKDIR}/default-ED_MTDOCK_MTSCOREE-sf"
ln -s ../1bzc+H.pdb
ln -s ../1bzc-sf.cif
CLI_ARGS="1bzc+H.pdb 1bzc-sf.cif      \
    --structures [-1..1] /A/TPI/902// \
   --mtdock opt torsion 15 0.01  \
    -p pdb -v 2 -O -h amberff14sb --np 10"
#    --mtdock 2,3 opt torsion pocket 100 0.01  \

echo ${QBHOME}/bin/qmechanic ${CLI_ARGS}
${QBHOME}/bin/qmechanic ${CLI_ARGS}

# now run the ED tests without density to make sure they work

echo "Default ED+MTDock+noDensity test"
mkdir -p "${WORKDIR}/default-ED_MTDOCK_MTSCOREE"
cd "${WORKDIR}/default-ED_MTDOCK_MTSCOREE"
ln -s ../1bzc+H.pdb
CLI_ARGS="1bzc+H.pdb     \
    --structures [-1..1] /A/TPI/902// \
   --mtdock opt torsion 15 0.01  \
    -p pdb -v 2 -O -h amberff14sb --np 10"
#    --mtdock 2,3 opt torsion pocket 100 0.01  \

echo ${QBHOME}/bin/qmechanic ${CLI_ARGS}
${QBHOME}/bin/qmechanic ${CLI_ARGS}

# testing the ED+MTDock+Density cases with density but with --ligand instead of selections

echo "Default ED+MTDock+Density (MTZ) test"
mkdir -p "${WORKDIR}/default-ED_MTDOCK_MTSCOREE-mtz"
cd "${WORKDIR}/default-ED_MTDOCK_MTSCOREE-mtz"
ln -s ../1bzc+H.pdb
ln -s ../1bzc.mtz
CLI_ARGS="1bzc+H.pdb 1bzc.mtz      \
    --ligand /A/TPI/902//       \
    --structures [-1..1]  \
   --mtdock opt torsion 15 0.01  \
    -p pdb -v 2 -O -h amberff14sb --np 10"
#    --mtdock 2,3 opt torsion pocket 100 0.01  \

echo ${QBHOME}/bin/qmechanic ${CLI_ARGS}
${QBHOME}/bin/qmechanic ${CLI_ARGS}

echo "Default ED+MTDock+Density (SF-CIF) test"
mkdir -p "${WORKDIR}/default-ED_MTDOCK_MTSCOREE-sf"
cd "${WORKDIR}/default-ED_MTDOCK_MTSCOREE-sf"
ln -s ../1bzc+H.pdb
ln -s ../1bzc-sf.cif
CLI_ARGS="1bzc+H.pdb 1bzc-sf.cif      \
    --ligand /A/TPI/902//       \
    --structures [-1..1]     \
   --mtdock opt torsion 15 0.01  \
    -p pdb -v 2 -O -h amberff14sb --np 10"
#    --mtdock 2,3 opt torsion pocket 100 0.01  \

echo ${QBHOME}/bin/qmechanic ${CLI_ARGS}
${QBHOME}/bin/qmechanic ${CLI_ARGS}

# now run the ED tests without density to make sure they work

echo "Default ED+MTDock+noDensity test"
mkdir -p "${WORKDIR}/default-ED_MTDOCK_MTSCOREE"
cd "${WORKDIR}/default-ED_MTDOCK_MTSCOREE"
ln -s ../1bzc+H.pdb
CLI_ARGS="1bzc+H.pdb     \
    --ligand /A/TPI/902//       \
    --structures [-1..1]     \
   --mtdock opt torsion 15 0.01  \
    -p pdb -v 2 -O -h amberff14sb --np 10"
#    --mtdock 2,3 opt torsion pocket 100 0.01  \

echo ${QBHOME}/bin/qmechanic ${CLI_ARGS}
${QBHOME}/bin/qmechanic ${CLI_ARGS}


