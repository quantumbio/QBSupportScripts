#!/bin/bash

#
#   Documentation:
#       The MTScoreE (ensemble) workflow is a 3 step process. This script encapulates these steps 
#           and demonstrates the use of the method with qmechanic and MOE from CCG. 
#           The script assumes that it is to run the pdb file (target) and a mol2 file representing 
#           the PLACED or DOCKED ligand which is also the "novel" ligand. Generally, the PLACED ligand is
#           used to define the pocket into which the novel ligand will be docked. In this case, the PLACED ligand 
#           and the NOVEL ligand are the same species / file.
#
#       INPUT Assumptions:
#           (1) The script takes two input
#
#           (2) A number of 
#

POSITIONAL=()
while [[ $# -gt 0 ]] ; do
    key="$1"
    case $key in
        -t|--target)
        targetfile="$2"
        shift # past argument
        shift # past value
        ;;
        -l|--ligand)
        ligandfile="$2"
        shift # past argument
        shift # past value
        ;;
        *)    # unknown option
        POSITIONAL+=("$1") # save it in an array for later
        shift # past argument
        ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters




if [ -z ${DIVCON_INSTALL} ]; then echo "ERROR: must set DIVCON_INSTALL env variable to QBHOME path"; exit; fi
if [ -z ${MOE_INSTALL} ]; then echo "ERROR: must set MOE_INSTALL env variable to QBHOME path"; exit; fi

if [ -z ${PREV_DIVCON_INSTALL} ]; then "WARNING: PREV_DIVCON_INSTALL not set. No test on previous version will be run."; fi     # Generally PREV_DIVCON_INSTALL is to test a previous version. If not included then this test will be skipped.
if [ -z ${MOE_SVL} ]; then MOE_SVL=${DIVCON_INSTALL}/svl fi

if [ -z ${MOE_MTSCOREES} ]; then MOEOPT_MTSES=""; else MOEOPT_MTSES="-mtscorees"; fi
if [ -z ${MTCSCONF} ]; then MTCSCONF=5; fi
if [ -z ${MAXPOSE} ]; then MAXPOSE=50; fi
if [ -z ${REMAXPOSE} ]; then REMAXPOSE=5; fi

# default MT options
if [ -z ${QM_OVERWRITE} ]; then QM_OVERWRITE="-O"; fi           # overwrite all output files with new files
if [ -z ${MT_HAM_TYPE} ]; then MT_HAM_TYPE="-h garf"; fi        # use the GARF statistical potential
if [ -z ${MT_MTCS_TYPE} ]; then MT_MTCS_TYPE="--mtcs input"; fi        # use the input bond lengths/angles in MTCS
if [ -z ${MTCSOPT} ]; then MTCSOPT=""; fi        # do not optimize the ligand conformation using the potential
if [ -z ${PBS_NUM_PPN} ]; then PBS_NUM_PPN=4; fi        # use 4 cores when we can 

echo "MOEOPT_MTSES: ${MOEOPT_MTSES}"
echo "MTCSCONF:     ${MTCSCONF}"
echo "MAXPOSE:      ${MAXPOSE}"
echo "REMAXPOSE:    ${REMAXPOSE}"
echo "QM_OVERWRITE: ${QM_OVERWRITE}"
echo "MT_HAM_TYPE:  ${MT_HAM_TYPE}"
echo "MT_MTCS_TYPE: ${MT_MTCS_TYPE}"
echo "MTCSOPT:      ${MTCSOPT}"
echo "PBS_NUM_PPN:  ${PBS_NUM_PPN}"

targetbasename=$1

# Go through all mol2 files in the previous directory and run each one as a 
for ligandfile in `ls ../*.mol2` ; do
    ligandbasename=`basename "$ligandfile" .mol2`
    rm -f OUT.${ligandbasename}_proteinE
    rm -f OUT.${ligandbasename}_proteinES
    
    if [ -e ${targetbasename}.pdb ] ; then
        ln -s -f ../${targetbasename}.pdb
    else
        echo "ERROR: ${PWD}/../${targetbasename}.pdb does not exist"
    fi
    ln -s -f ../${ligandbasename}.mol2

    starttime=`date`
    echo "$ligandbasename STARTTIME: $starttime" &> OUT.${ligandbasename}_proteinE
    
    env >> OUT.${ligandbasename}_proteinE 2>&1

    timeout 60m ${DIVCON_INSTALL}/bin/qmechanic ${targetbasename}.pdb --ligand ${ligandbasename}.mol2 ${QM_OVERWRITE} ${MT_HAM_TYPE} ${MT_MTCS_TYPE} ${MTCSOPT} -p sdf --np ${PBS_NUM_PPN} -v 2 >> OUT.${ligandbasename}_proteinE 2>&1
    if [ $? -eq 124 ] ; then
        echo "TIMEOUT ERROR (60min): ${targetbasename}.pdb ${ligandbasename}.mol2 --mtcs - CHECK STRUCTURE FOR ERRORS/ODDITIES" >> OUT.${ligandbasename}_proteinE
    else
        echo "QMECHANIC RUN COMPLETE" >> OUT.${ligandbasename}_proteinE
    fi
    
    if [ -e ${targetbasename}_conf.sdf ] ; then
        mv ${targetbasename}_conf.sdf ${ligandbasename}_conf.sdf
    fi

    ${MOE_INSTALL} -licwait -run "${MOE_SVL}" -rec ${targetbasename}.pdb -lig ${ligandbasename}.mol2 -conf ${ligandbasename}_conf.sdf -delwat $MOEOPT_MTSES --maxpose ${MAXPOSE} --mtcsconf ${MTCSCONF} --remaxpose ${REMAXPOSE} >> OUT.${ligandbasename}_proteinE 2>&1
    echo "MOE RUN COMPLETE" >> OUT.${ligandbasename}_proteinE

    timeout 60m ${DIVCON_INSTALL}/bin/qmechanic ${targetbasename}.pdb --ligand ${ligandbasename}.mol2 ${QM_OVERWRITE} ${MT_HAM_TYPE} --mtdock ${ligandbasename}_dock.sdf opt off --mtscore --np ${PBS_NUM_PPN} -v 2 >> OUT.${ligandbasename}_proteinE 2>&1
    if [ $? -eq 124 ] ; then
        echo "TIMEOUT ERROR (60min): ${targetbasename}.pdb ${ligandbasename}.mol2 --mtscore" >> OUT.${ligandbasename}_proteinE
    else
        echo "QMECHANIC RUN COMPLETE" >> OUT.${ligandbasename}_proteinE
    fi
    
    endtime=`date`
    echo "$ligandbasename ENDTIME: $endtime" >> OUT.${ligandbasename}_proteinE 2>&1
   
    rm -f *.h5

    if [ ! -z ${PREV_DIVCON_INSTALL} ]; then
        jobname=`basename ${PREV_DIVCON_INSTALL} | sed 's/DivConDiscoverySuite-//'`
        timeout 60m ${PREV_DIVCON_INSTALL}/bin/qmechanic ${targetbasename}.pdb --ligand ${ligandbasename}.mol2 ${QM_OVERWRITE} ${MT_HAM_TYPE} --mtdock ${ligandbasename}_dock.sdf opt off --mtscore --np ${PBS_NUM_PPN} -v 2 > OUT.${ligandbasename}_proteinE-${jobname} 2>&1
        if [ $? -eq 124 ] ; then
            echo "TIMEOUT ERROR (60min): ${targetbasename}.pdb ${ligandbasename}.mol2 --mtscore" >> OUT.${ligandbasename}_proteinE-${jobname}
        else
            echo "QMECHANIC RUN COMPLETE" >> OUT.${ligandbasename}_proteinE-${jobname}
        fi

        rm -f *.h5
    endif

    mdbfile=${ligandbasename}_dock.mdb
    dockedbasename=`basename "${mdbfile}" .mdb`
    timeout 60m  ${DIVCON_INSTALL}/bin/qmechanic ${targetbasename}.pdb --ligand ${dockedbasename}.mol2 ${QM_OVERWRITE} ${MT_HAM_TYPE} --mtscore endstate --np ${PBS_NUM_PPN} -v 2 > OUT.${ligandbasename}_proteinES
    if [ $? -eq 124 ] ; then
        echo "TIMEOUT ERROR (60min): ${targetbasename}.pdb ${ligandbasename}.mol2 --mtscore endstate" >> OUT.${ligandbasename}_proteinES 2>&1
    else
        echo "QMECHANIC RUN COMPLETE" >> OUT.${ligandbasename}_proteinES
    fi

    rm -f *.h5 *.log
    
done
