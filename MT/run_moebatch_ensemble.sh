#!/bin/bash

#
#   Documentation:
#       This script runs qmechanic-mtcs / moedock (mtcs confomers) / qmechanic-mtscoreE
#


if [ -z ${MOE_MTSCOREES} ]; then MOEOPT_MTSES=""; else MOEOPT_MTSES="-mtscorees"; fi
if [ -z ${MTCSCONF} ]; then MTCSCONF=5; fi
if [ -z ${MAXPOSE} ]; then MAXPOSE=50; fi
if [ -z ${REMAXPOSE} ]; then REMAXPOSE=5; fi

echo "MOEOPT_MTSES: ${MOEOPT_MTSES}"
echo "MTCSCONF:     ${MTCSCONF}"
echo "MAXPOSE:      ${MAXPOSE}"
echo "REMAXPOSE:    ${REMAXPOSE}"

targetbasename=$1

for ligandfile in `ls ../*.mol2` ; do
    ligandbasename=`basename "$ligandfile" .mol2`
    rm -f OUT.${ligandbasename}_proteinE
    rm -f OUT.${ligandbasename}_proteinES

    ln -s -f ../${targetbasename}.pdb
    ln -s -f ../${ligandbasename}.mol2

    starttime=`date`
    echo "$ligandbasename STARTTIME: $starttime" &> OUT.${ligandbasename}_proteinE
    
    env >> OUT.${ligandbasename}_proteinE 2>&1

    timeout 60m ${DIVCON_INSTALL}/bin/qmechanic ${targetbasename}.pdb --ligand ${ligandbasename}.mol2 ${QM_OVERWRITE} ${MT_HAM_TYPE} ${MT_MTCS_TYPE} ${MTCSOPT} -p sdf --np ${PBS_NUM_PPN} -v 2 >> OUT.${ligandbasename}_proteinE 2>&1
    if [ $? -eq 124 ] ; then
        echo "TIMEOUT ERROR (60min): ${targetbasename}.pdb ${ligandbasename}.mol2 --mtcs" >> OUT.${ligandbasename}_proteinE
    else
        echo "QMECHANIC RUN COMPLETE" >> OUT.${ligandbasename}_proteinE
    fi
    
    if [ -e ${targetbasename}_conf.sdf ] ; then
        mv ${targetbasename}_conf.sdf ${ligandbasename}_conf.sdf
    fi
#    /share/apps/MOE/moe2016/bin/moebatch -licwait -run "/home/lance/Release/current/svl/run/qbDockPair.svl" -rec ${targetbasename}.pdb -lig ${ligandbasename}.mol2 -conf ${ligandbasename}_conf.sdf -delwat $MOEOPT_MTSES --maxpose ${MAXPOSE} --mtcsconf ${MTCSCONF} --remaxpose ${REMAXPOSE} >> OUT.${ligandbasename}_proteinE 2>&1
#    /share/apps/MOE/moe2018/bin/moebatch -licwait -run "${DIVCON_INSTALL}/svl/run/qbDockPair.svl" -rec ${targetbasename}.pdb -lig ${ligandbasename}.mol2 -conf ${ligandbasename}_conf.sdf -delwat $MOEOPT_MTSES --maxpose ${MAXPOSE} --mtcsconf ${MTCSCONF} --remaxpose ${REMAXPOSE} >> OUT.${ligandbasename}_proteinE 2>&1
#    /share/apps/MOE/moe2018/bin/moebatch -licwait -run "/home/lance/Dropbox/Developer/MOEDivCon/svl/run/qbDockPair.svl" -rec ${targetbasename}.pdb -lig ${ligandbasename}.mol2 -conf ${ligandbasename}_conf.sdf -delwat $MOEOPT_MTSES --maxpose ${MAXPOSE} --mtcsconf ${MTCSCONF} --remaxpose ${REMAXPOSE} >> OUT.${ligandbasename}_proteinE 2>&1
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

    jobname=`basename ${PREV_DIVCON_INSTALL} | sed 's/DivConDiscoverySuite-//'`
    
    timeout 60m ${PREV_DIVCON_INSTALL}/bin/qmechanic ${targetbasename}.pdb --ligand ${ligandbasename}.mol2 ${QM_OVERWRITE} ${MT_HAM_TYPE} --mtdock ${ligandbasename}_dock.sdf opt off --mtscore --np ${PBS_NUM_PPN} -v 2 > OUT.${ligandbasename}_proteinE-${jobname} 2>&1
    if [ $? -eq 124 ] ; then
        echo "TIMEOUT ERROR (60min): ${targetbasename}.pdb ${ligandbasename}.mol2 --mtscore" >> OUT.${ligandbasename}_proteinE-${jobname}
    else
        echo "QMECHANIC RUN COMPLETE" >> OUT.${ligandbasename}_proteinE-${jobname}
    fi

    rm -f *.h5

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
