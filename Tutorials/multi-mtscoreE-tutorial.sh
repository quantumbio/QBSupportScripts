#!/bin/bash

# if not running in PBS, need to set these env variables:
if [ -z ${PBS_JOBID} ]; then export PBS_JOBID=`date +%H%M%S%m%d%y`; fi
if [ -z ${PBS_O_WORKDIR} ]; then export PBS_O_WORKDIR=${PWD}; fi
if [ -z ${PBS_NUM_PPN} ]; then export PBS_NUM_PPN=10; fi

# Multi-ligand example: CMET from MCompChem benchmark https://github.com/MCompChem/fep-benchmark/tree/master/cmet

export MOE_INSTALL=/share/apps/MOE/moe2020.0901
export DIVCON_INSTALL=/home/lance/Release/DivConDiscoverySuite-2022-b4882
export MOE_SVL=${DIVCON_INSTALL}/svl/run/qbDockPair.svl

# support function
function max_bg_procs {
    if [[ $# -eq 0 ]] ; then
            echo "Usage: max_bg_procs NUM_PROCS.  Will wait until the number of background (&)"
            echo "           bash processes (as determined by 'jobs -pr') falls below NUM_PROCS"
            return
    fi
    local max_number=$((0 + ${1:-0}))
    while true; do
            local current_number=$(jobs -pr | wc -l)
            if [[ $current_number -lt $max_number ]]; then
                    break
            fi
            sleep 5
    done
}


# STEP 1: mimic what is done on the CLIENT computer (this is in the MOE GUI for example)

# input files from the MCompChem project:
#wget https://raw.githubusercontent.com/MCompChem/fep-benchmark/master/cmet/ligands.sdf
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Examples/5-cmet-ligands.sdf
wget https://raw.githubusercontent.com/MCompChem/fep-benchmark/master/cmet/4r1y_prepared.pdb

# SVL execution files from the QBSupportScripts project:
wget https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/svl/run/sdf2sdf.svl


targetfile="4r1y_prepared.pdb"
targetbasename=`basename "$targetfile" .pdb`

${MOE_INSTALL}/bin/moebatch -licwait -run ./sdf2sdf.svl -docksdf 5-cmet-ligands.sdf -split -outf mol2

# STEP 2: run the MOE/dock tasks (one on each ligand) - this can be done completely in parallel

for ligandfile in `ls 5-cmet-ligands-MOE-*.mol2` ; do
    max_bg_procs $( echo "scale=0;((${PBS_NUM_PPN})+0.5)/2" | bc )
    ligandbasename=`basename "$ligandfile" .mol2`
    echo "${MOE_INSTALL}/bin/moebatch -licwait -run "${MOE_SVL}" -rec ${targetbasename}.pdb -lig ${ligandfile} -delwat -protonate"
    ${MOE_INSTALL}/bin/moebatch -licwait -run "${MOE_SVL}" -rec ${targetbasename}.pdb -lig ${ligandfile} -delwat -protonate >> OUT.${ligandbasename}_proteinE 2>&1  &
done
wait    # wait until all ligand tasks (docking) have been run

# STEP 3: run the DivCon/MTScore tasks (one on each ligand) - this can be done completely in parallel

for ligandfile in `ls 5-cmet-ligands-MOE-*.sdf` ; do
    max_bg_procs $( echo "scale=0;((${PBS_NUM_PPN})+0.5)/2" | bc )
    ligandbasename=`basename "$ligandfile" .sdf`
    DOCKFILES="${ligandbasename}_dock.sdf"
    echo "${DIVCON_INSTALL}/bin/qmechanic pro_${ligandbasename}_predock.pdb --ligand lig_${ligandbasename}_predock.mol2 -O -h garf --mtdock ${DOCKFILES} --mtscore endstate --np 2 -v 2 -p sdf"
    ${DIVCON_INSTALL}/bin/qmechanic pro_${ligandbasename}_predock.pdb --ligand lig_${ligandbasename}_predock.mol2 -O -h garf --mtdock ${DOCKFILES} --mtscore endstate --np 2 -v 2 -p sdf >> OUT.${ligcorename}_proteinE-garf 2>&1 &
done
wait    # wait until all ligand tasks (MTScore) have been run
