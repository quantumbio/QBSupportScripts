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


targetbasename=$1

for ligandfile in `ls ../*.mol2` ; do
    ligandbasename=`basename "$ligandfile" .mol2`
    rm -f OUT.${ligandbasename}_proteinE
    rm -f OUT.${ligandbasename}_proteinES

    ln -s -f ../${targetbasename}.pdb
    ln -s -f ../${ligandbasename}.mol2
    
    ${DIVCON_INSTALL}/bin/qmechanic ${targetbasename}.pdb --ligand ${ligandbasename}.mol2 ${QM_OVERWRITE} ${MT_HAM_TYPE} --mtscore endstate -v 2 &> OUT.${ligandbasename}_proteinES
    echo "QMECHANIC RUN COMPLETE" >> OUT.${ligandbasename}_proteinES
    rm -f *.h5

    starttime=`date`
    echo "$ligandbasename STARTTIME: $starttime" &> OUT.${ligandbasename}_proteinE
    
    env >> OUT.${ligandbasename}_proteinE 2>&1

    timeout 240m ${DIVCON_INSTALL}/bin/qmechanic ${targetbasename}.pdb --ligand ${ligandbasename}.mol2 ${QM_OVERWRITE} ${MT_HAM_TYPE} ${MT_MTCS_TYPE} ${MTCSOPT} ${MTDOCK_TYPE} ${MTDOCK_OPT} --mtscore -p sdf --np 4 -v 2 >& OUT.${ligandbasename}_proteinE
    if [ $? -eq 124 ] ; then
        echo "TIMEOUT ERROR (240m): ${targetbasename}.pdb ${ligandbasename}.mol2 --mtscore" >> OUT.${ligandbasename}_proteinE
    else
        echo "QMECHANIC RUN COMPLETE" >> OUT.${ligandbasename}_proteinE
    fi
    rm -f *.h5 *.log
    
    endtime=`date`
    echo "$ligandbasename ENDTIME: $endtime" >> OUT.${ligandbasename}_proteinE 2>&1

done
