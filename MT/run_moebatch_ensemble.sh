#!/bin/bash

#  // BEGIN COPYRIGHT
#  /***********************************************************************
#     Copyright (c) 2020-2021 QuantumBio Inc. and/or its affiliates.
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


# START VARIABLES
bNOVEL=1        # assume that a novel and a placed ligand will be provided (set to 0 later if not)
bTHREESTEP=0    # assume that the two step process (skipping MTCS) will be performed.
bFORCEEXE=""     # assume the protein/ligand is not ready to roll and therefore throw an error if MOE finds a problem.

# SET DEFAULTS:
if [ -z "${MOE_MTSCOREES}" ]; then MOEOPT_MTSES=""; else MOEOPT_MTSES="-mtscorees"; fi
if [ -z "${MTCSCONF}" ]; then MTCSCONF=0; fi
if [ -z "${MAXPOSE}" ]; then MAXPOSE=125; fi
if [ -z "${REMAXPOSE}" ]; then REMAXPOSE=25; fi

# DEFAULT MT options (these may be set in the environment prior to calling the script if you wish to change them)
if [ -z "${QM_OVERWRITE}" ]; then QM_OVERWRITE="-O"; fi             # overwrite all output files with new files
if [ -z "${MT_HAM_TYPE}" ]; then MT_HAM_TYPE="-h garf"; fi          # use the GARF statistical potential
if [ -z "${MT_MTCS_TYPE}" ]; then MT_MTCS_TYPE="--mtcs input"; fi   # use the input bond lengths/angles in MTCS
if [ -z "${MTCSOPT}" ]; then MTCSOPT=""; fi                         # do not optimize the ligand conformation using the potential
if [ -z "${PBS_NUM_PPN}" ]; then PBS_NUM_PPN=4; fi                  # use 4 cores when we can 

# CHECK FOR EXISTENCE OF DIVCON / QMECHANIC (see ERROR TRAP BELOW)
if [ ! -z "${DIVCON_INSTALL}" ]; then
    QMECHANIC_BIN="${DIVCON_INSTALL}/bin/qmechanic"
elif [ ! -z "${QBHOME}" ]; then
    QMECHANIC_BIN="${QBHOME}/bin/qmechanic"
else
    QMECHANIC_BIN=`which qmechanic`
fi
if [ -f "${QMECHANIC_BIN}" ] ; then
    progdir=`dirname "$QMECHANIC_BIN"`
    DIVCON_INSTALL=`dirname "$progdir"`
else
    DIVCON_INSTALL=""
    QMECHANIC_BIN=""
fi

# CHECK FOR EXISTENCE OF PREVIOUS DIVCON / QMECHANIC (see ERROR TRAP BELOW)
if [ ! -z "${PREV_DIVCON_INSTALL}" ]; then
    PREV_QMECHANIC_BIN="${PREV_DIVCON_INSTALL}/bin/qmechanic"
else
    PREV_QMECHANIC_BIN=""
fi
if [ -f "${PREV_QMECHANIC_BIN}" ] ; then
    progdir=`dirname "$PREV_QMECHANIC_BIN"`
    PREV_DIVCON_INSTALL=`dirname "$progdir"`
else
    PREV_DIVCON_INSTALL=""
    PREV_QMECHANIC_BIN=""
fi

# CHECK FOR EXISTENCE OF MOE / moebatch (see ERROR TRAP BELOW)
if [ ! -z "${MOE_INSTALL}" ]; then
    MOE_BIN="${MOE_INSTALL}/bin/moebatch"
else
    MOE_BIN=`which moebatch`
fi
if [ -f "${MOE_BIN}" ] ; then
    progdir=`dirname "$MOE_BIN"`
    MOE_INSTALL=`dirname "$progdir"`
else
    MOE_INSTALL=""
    MOE_BIN=""
fi

# CHECK FOR EXISTENCE OF DC_SVL/svl/run/qbDockPair.svl
if [ ! -z "${DC_SVL}" ]; then
    DC_SVL_DP="${DC_SVL}/run/qbDockPair.svl"
else
    DC_SVL_DP="${DIVCON_INSTALL}/svl/run/qbDockPair.svl"
fi
if [ -f "${DC_SVL_DP}" ] ; then
    progdir=`dirname "$DC_SVL_DP"`
    DC_SVL=`dirname "$progdir"`
else
    DC_SVL=""
    DC_SVL_DP=""
fi

# SCIPT FUNCTION DEFINITIONS

# Report the current settings as of the point the function is executed
current_settings()
{
    echo "--------------------------------------------------------------------"
    echo "Current Job Settings: "
    echo " *DIVCON_INSTALL:      ${DIVCON_INSTALL}"
    echo " *MOE_INSTALL:         ${MOE_INSTALL}"
    echo " *DC_SVL:              ${DC_SVL}"
    echo "  PREV_DIVCON_INSTALL: ${PREV_DIVCON_INSTALL}"
    echo ""
    echo "  MOEOPT_MTSES:        ${MOEOPT_MTSES}"
    echo "  MTCSCONF:            ${MTCSCONF}"
    echo "  MAXPOSE:             ${MAXPOSE}"
    echo "  REMAXPOSE:           ${REMAXPOSE}"
    echo "  QM_OVERWRITE:        ${QM_OVERWRITE}"
    echo "  MT_HAM_TYPE:         ${MT_HAM_TYPE}"
    echo "  MT_MTCS_TYPE:        ${MT_MTCS_TYPE}"
    echo "  MTCSOPT:             ${MTCSOPT}"
    echo "  PBS_NUM_PPN:         ${PBS_NUM_PPN}"
    echo "  THREE_STEP:          ${bTHREESTEP}  // 0 = OFF"
    echo "  FORCE_EXE:           ${bFORCEEXE}    // "" = OFF"
    echo ""
    echo "* = MUST BE SET OR ERROR WILL THROW"
    echo "--------------------------------------------------------------------"
}

# Error/help text used to provide information to the user about why his/her script died
error_exit()
{
    echo ""
	echo "$1"
	echo ""
	echo "Usage: $0 --target protein.pdb --ligand placed_ligand.mol2 [novel_ligand.mol2] [--tlimit 60]"
	echo ""
	echo "Documentation:"
    echo "  The MTScoreE (ensemble) workflow is a 2 or 3 step process. This script encapulates these steps"
    echo "  and demonstrates the use of the method with qmechanic for scoring and MOE from CCG for docking."
    echo "  The script requires that you provide a fully prepared PDB file for the target alone and "
    echo "  a fully prepared mol2 file for the PLACED [wildtype] ligand to define the active site."
    echo "  You may also provide a NOVEL ligand which will be put in place of the PLACED ligand."
    echo "  If the NOVEL ligand is not provided, then the PLACED ligand will be used instead."
	echo ""
	echo "Command Line Options:"
	echo "  --help     Print this help and exit."
	echo "  --target   [required] A PDB file corresponding to a protonated PROTEIN (no ligand)"
	echo "  --ligand   [required] A MOL2 file corresponding to a protonated LIGAND placed in active site".
	echo "             - This placed_ligand is used to define the active site residues for docking."
	echo "             - You may use a wildcard '*' (such as *.mol2) in order to repeat this process for"
	echo "             - many placed ligands (if --ligand has a '*' then --novel may NOT be provided)."
	echo "  --novel    [optional] A MOL2 file corresponding to a protonated, novel ligand."
	echo "             - (Script will use --ligand if --novel is not defined.)"
	echo "             - You may use a wildcard '*' (such as *.mol2) in order to repeat this process for"
	echo "             - many placed ligands (if --ligand has a '*' then --novel may NOT be provided)."
	echo "  --tlimit   [optional] Time limit in minutes for each qmechanic (MTCS, MTScore) run."
	echo "             - (Script will use 60 minutes if --tlimit is not defined.)"
	echo "             - Note: generally these calculations should take a few seconds/minutes so "
	echo "               you should never see this limit reached."
	echo "  --3step    [optional] By default, currently MTCS is only used to generate the unbound ZL."
	echo "             - Using the 3step option, these MTCS conformers will also be used in MOE/Dock."
	echo "  --forceexe [optional] By default, MOE will throw an error if it finds any trouble in the input."
	echo "             - Using the forceexe option, the calculation will continue."
	echo ""
	echo "Useful env variables include:"
    echo "  DIVCON_INSTALL  = path to QBHOME to be used"
    echo "  MOE_INSTALL     = path to MOE/batch"
	echo ""
    echo "Optional env variables include:"
    echo "  MT_HAM_TYPE     = pair potential to be used"
    echo "  MT_MTCS_TYPE    = ligand (MTCS) lengths/angles to use"
    echo "  MTCSOPT         = ligand (MTCS) optimization options"
    echo "  PBS_NUM_PPN     = number of processors cores to use"
	echo ""
    echo "For more guidance, visit: http://www.quantumbioinc.com/support/manual/movabletype"
    echo "                Or email: support@quantumbioinc.com"
	echo ""
	current_settings
	exit 1
}

# Execution function: used to wrap any call to qmechanic / moebatch for timeout and other reporting
execute_binary()
{
    echo "INFO: RUNNING: $1"
    echo "INFO: SEE ${PWD}/$2 to observe progress / errors / success / results"
    # check for timeout and use it if available to run the job, otherwise just move forward (not all distributions have timeout)
    timeoutbin=`which timeout`
    if [ -z "${timeoutbin}" ] ; then
        $1 >> $2 2>&1
    else
        timeout --foreground ${timeoutlimit}m $1 >> $2 2>&1
        if [ $? -eq 124 ] ; then
            echo "TIMEOUT ERROR (${timeoutlimit} min): CHECK STRUCTURE FOR ERRORS/ODDITIES" >> "$2"
            echo "TIMEOUT ERROR (${timeoutlimit} min): CHECK STRUCTURE FOR ERRORS/ODDITIES"
            exit
        fi
    fi
}

# MAIN CODE BEGINS

# Script argument processing:
POSITIONAL=()
while [[ $# -gt 0 ]] ; do
    key="$1"
    case $key in
        -t|--target)
            targetfile="$2"
            if [ ! -f "${targetfile}" ]; then error_exit "ERROR: ${targetfile} does not exist"; fi
            shift # past argument
            shift # past value
        ;;
        -l|--ligand)
            unset inligandfile
            shift # past argument
            filearray=( "$@" )
            i=0
            for file in "${filearray[@]}" ; do
                if [[ "${file}" == "--"* ]] ; then break; fi    # this is where next argument (if exists) hits.
                if [ ! -f "${file}" ]; then error_exit "ERROR: ${file} does not exist"; fi
                inligandfile=( "${inligandfile[@]}" "$file" )
                shift # past value
            done
        ;;
        -n|--novel)
            unset innovelfile
            shift # past argument
            filearray=( "$@" )
            i=0
            for file in "${filearray[@]}" ; do
                if [[ "${file}" == "--"* ]] ; then break; fi    # this is where next argument (if exists) hits.
                if [ ! -f "${file}" ]; then error_exit "ERROR: ${file} does not exist"; fi
                innovelfile=( "${innovelfile[@]}" "$file" )
                shift # past value
            done
        ;;
        --tlimit)
            timeoutlimit="$2"
            shift # past argument
            shift # past value
        ;;
        --3step)
            bTHREESTEP=1
            shift # past argument
        ;;
        --forceexe)
            bFORCEEXE="-forceexe"
            shift # past argument
        ;;
        -h|--hamiltonian)
            hamiltonian="$2"
            shift # past argument
            shift # past value
        ;;
        --help)
        error_exit "Help:"         
        ;;
        *)    # unknown option
            
            error_exit "ERROR: unknown command line argument $1" 
        ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

# ERROR TRAPS
if [ -z "${targetfile}" ] ; then error_exit "ERROR: you must define --targetfile as the protein as a separate pdb file."; fi
if [ -z "${inligandfile}" ] ; then error_exit "ERROR: you must define --ligand as the PLACED [wildtype] ligand as a separate mol2 file."; fi

if [ ! -z "${innovelfile}" ] && [ "${#inligandfile[@]}" -gt 1 ] ; then error_exit "ERROR: PLACED [WT] ligand may not include wildcard '*' when NOVEL ligand is provided."; fi

if [ -z "${innovelfile}" ] ; then
    echo "NOTE: No NOVEL ligand provided - assuming self-docking/scoring of PLACED [WT] ligand."
    echo "          Use --novel to define alternate ligand for docking/scoring."
    b=("${a[@]}") 
    innovelfile=("${inligandfile[@]}") 
    bNOVEL=0
fi

if [ -z "${timeoutlimit}" ] ; then timeoutlimit=60; fi

if [ -z "${QMECHANIC_BIN}" ]; then error_exit "ERROR: qmechanic executable not found!\n    Set DIVCON_INSTALL env variable or qmechanic must be in the system PATH"; fi
if [ -z "${MOE_BIN}" ]; then error_exit "ERROR: moebatch executable not found!\n    Set MOE_INSTALL env variable or moebatch must be in the system PATH" ; fi
if [ -z "${DC_SVL_DP}" ]; then error_exit "ERROR: DivCon Suite SVL NOT FOUND!\n    Set DC_SVL env variable or svl directory must be in DIVCON_INSTALL directory" ; fi

if [ -z "${PREV_DIVCON_INSTALL}" ]; then echo "NOTE: PREV_DIVCON_INSTALL not set. No test on previous version will be run." ; fi     # Generally PREV_DIVCON_INSTALL is to test a previous version. If not included then this test will be skipped.

# BEGIN JOB

current_settings

targetbasename=`basename "$targetfile" .pdb`
if [ ! -e ${targetbasename}.pdb ] ; then error_exit "ERROR: ${PWD}/${targetbasename}.pdb does not exist" ; fi

# Go through all mol2 files in the previous directory and run each one as a 
for novelfile in "${innovelfile[@]}" ; do
    novelbasename=`basename "$novelfile" .mol2`
    if [ "$bNOVEL" -eq 0 ]; then
        ligandbasename=${novelbasename}
    else
        ligandbasename=`basename "$inligandfile" .mol2`
    fi
    echo "Running: $novelfile $ligandbasename"

    rm -f OUT.${novelbasename}_proteinE
    rm -f OUT.${novelbasename}_proteinES
    
    if [ ! -e ${novelbasename}.mol2 ] ; then echo "ERROR: ${PWD}/${novelbasename}.mol2 does not exist"; continue; fi
    
    starttime=`date`
    echo "$novelbasename STARTTIME: $starttime" &> OUT.${novelbasename}_proteinE
    
    env >> OUT.${novelbasename}_proteinE 2>&1

#   STEP #1: EXECUTE MTCS (conformational search) to generate conformers which match the chosen potential. (for 3 step)
    if [ "$bTHREESTEP" -eq 1 ]; then
        execute_binary "${DIVCON_INSTALL}/bin/qmechanic ${targetbasename}.pdb --ligand ${novelbasename}.mol2 ${QM_OVERWRITE} ${MT_HAM_TYPE} ${MT_MTCS_TYPE} ${MTCSOPT} -p sdf --np ${PBS_NUM_PPN} -v 2" "OUT.${novelbasename}_proteinE"
        echo "QMECHANIC RUN COMPLETE" >> OUT.${novelbasename}_proteinE
        if [ -e ${targetbasename}_conf.sdf ] ; then mv ${targetbasename}_conf.sdf ${novelbasename}_conf.sdf;  fi
        MOE_CONF_FILE="-conf ${novelbasename}_conf.sdf --mtcsconf ${MTCSCONF}"
    fi

#   STEP #2: EXECUTE MOE (docking) to generate poses which fit within the active site of the placed_ligand.mol2
    execute_binary "${MOE_BIN} -licwait -run ${DC_SVL}/run/qbDockPair.svl -rec ${targetbasename}.pdb -lig ${ligandbasename}.mol2 ${MOE_CONF_FILE} -o ${novelbasename}_dock -delwat $MOEOPT_MTSES -maxpose ${MAXPOSE} -remaxpose ${REMAXPOSE} ${bFORCEEXE}" "OUT.${novelbasename}_proteinE"
    echo "MOE RUN COMPLETE" >> OUT.${novelbasename}_proteinE

#   STEP #3: EXECUTE MTScore (Ensemble scoring) to score MOE-generated poses which fit within the active site of the placed_ligand.mol2
    execute_binary "${DIVCON_INSTALL}/bin/qmechanic pro_${ligandbasename}_predock.pdb --ligand lig_${ligandbasename}_predock.mol2 ${QM_OVERWRITE} ${MT_HAM_TYPE} --mtdock ${novelbasename}_dock.sdf opt off --mtscore --np ${PBS_NUM_PPN} -v 2" "OUT.${novelbasename}_proteinE"
    echo "QMECHANIC RUN COMPLETE" >> OUT.${novelbasename}_proteinE
    
    endtime=`date`
    echo "$novelbasename ENDTIME: $endtime" >> OUT.${novelbasename}_proteinE 2>&1
   
    rm -f *.h5

#   OPTIONAL STEP #3: EXECUTE MTScore (Ensemble scoring) using an alternate chosen qmechanic to score MOE-generated poses which fit within the active site of the placed_ligand.mol2
    if [ ! -z "${PREV_DIVCON_INSTALL}" ]; then
        jobname=`basename ${PREV_DIVCON_INSTALL} | sed 's/DivConDiscoverySuite-//'`
        execute_binary "${PREV_DIVCON_INSTALL}/bin/qmechanic pro_${targetbasename}_predock.pdb --ligand lig_${ligandbasename}_predock.mol2 ${QM_OVERWRITE} ${MT_HAM_TYPE} --mtdock ${novelbasename}_dock.sdf opt off --mtscore --np ${PBS_NUM_PPN} -v 2" "OUT.${novelbasename}_proteinE-${jobname}"

        rm -f *.h5
    fi

#   OPTION STEP #4: EXECUTE MTScore (EndState scoring) to score the "winning" or "best" MOE-generated pose which fit within the active site of the placed_ligand.mol2
    mdbfile=${novelbasename}_dock.mdb
    dockedbasename=`basename "${mdbfile}" .mdb`
    execute_binary "${DIVCON_INSTALL}/bin/qmechanic ${targetbasename}.pdb --ligand ${dockedbasename}.mol2 ${QM_OVERWRITE} ${MT_HAM_TYPE} --mtscore endstate --np ${PBS_NUM_PPN} -v 2"  "OUT.${novelbasename}_proteinES"

    rm -f *.h5 *.log
    
done
