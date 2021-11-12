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


POSITIONAL=()

fileName=${1}

while [[ $# -gt 0 ]] ; do
    key="$1"
    case $key in
        -s|--scale)
        SCALE=YES
        shift # past argument
        ;;
        -d|--details)
        PRINTDETAILS=YES
        shift # past argument
        ;;
        *)    # unknown option
        POSITIONAL+=("$1") # save it in an array for later
        shift # past argument
        ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

echo "SCALE: ${SCALE}"
echo "PRINTDETAILS: ${PRINTDETAILS}"

PROCFILE="proc.out"

# defaut scaling is no scaling
mES=1
bES=0
mE=1
bE=0

# default intra-calculation scale for calls to bc
intscale="scale=10"

checkNumber ()
{
    re='^[+-]?[0-9]+([.][0-9]+)?$'
    if ! [[ $1 =~ $re ]] ; then
       return 1
    else
       return 0
    fi
}

WORKDIR=$PWD

dirlist=`ls -d */`
if test -z "$dirlist" ; then
    dirlist="./"
fi

for dir in $dirlist ; do
    dir=${dir%/}

    cd ${WORKDIR}/$dir

    rm -f ${PROCFILE}
    grandtotaltime=0
    jobcount=0
    for outfile in `ls OUT.*protein${fileName}` ; do
        ligand=`echo $outfile | awk -F '.' '{print $2}' | awk -F "_proteinE" '{print $1}'`
        scoreline=`grep $ligand $outfile | grep -v TIME | tail -1`
        echo $scoreline  | sed 's/"//g' >> ${PROCFILE}
        
    # tally the time
        # MTCS
        mtcstime=`grep "Total Computation Time (Seconds)" ${outfile} | awk '{print $5}' | head -1`
        # MTScoreE or the entirety of MTCS+MTDOCK+MTSCOREE
        mtscoretime=`grep "Total Computation Time" ${outfile} | tail -1 | awk '{print $5}'`
        
        # MOEDOCK
        startdock=`grep "db_ImportSD" ${outfile} | head -1 | awk -F "db_ImportSD:" '{print $2}'`
        enddock=`grep "Docking finished in" ${outfile} | tail -1 | awk '{print $4}'`
        if [ -z "$startdock" ] || [ -z "$enddock" ] ; then
            totaltime=${mtscoretime}
        else
            startdock=`date -d"${startdock}" +%s`
            enddock=${enddock#-}
            totaltime=`echo "${mtcstime} + ${enddock} + ${mtscoretime}" | bc`
        fi
        grandtotaltime=`echo "${grandtotaltime} + ${totaltime}" | bc`
        
        jobcount=`echo "1+${jobcount}" | bc`
    done
    
    averagetime=`echo "scale=1;${grandtotaltime}/${jobcount}/60" | bc`

    sed -i 's/_dock//g' ${PROCFILE}

# calculate scaling factors
    sESCount=0
    sECount=0
    sESTotal=0
    sETotal=0
    sESAvg=0
    sEAvg=0

    expCount=0
    expTotal=0
    expAvg=0
    
    while read line ; do
        firstone=`echo $line | awk '{print $1}' | sed 's/ //g'| sed 's/\n//g'`
        exp_dG=`echo $line | awk '{print $2}' | sed 's/ //g'| sed 's/\n//g'`
                
        firstName=`grep "$firstone\b" ${PROCFILE} | awk '{print $1}'| sed 's/\n//g'`
        firstmtE=`grep "$firstone\b" ${PROCFILE} | awk '{print $6}'| sed 's/\n//g'`
        firstmtES=`grep "$firstone\b" ${PROCFILE} | awk '{print $2}'| sed 's/\n//g'`

        if checkNumber $exp_dG ; then
            expCount=`echo "$expCount + 1" | bc`
            expTotal=`echo "$expTotal + $exp_dG" | bc`
        fi
        
        if checkNumber $firstmtES ; then
            sESCount=`echo "$sESCount + 1" | bc`
            sESTotal=`echo "${intscale}; $sESTotal + $firstmtES" | bc`
        fi

        if checkNumber $firstmtE ; then
            sECount=`echo "$sECount + 1" | bc`
            sETotal=`echo "${intscale}; $sETotal + $firstmtE" | bc`
        fi  
    done < "../../singles.txt"
    
    echo "$sESCount | $sECount | ${expCount}"
    sESAvg=`echo "${intscale}; $sESTotal / $sESCount" | bc`
    sEAvg=`echo "${intscale}; $sETotal / $sECount" | bc`
    expAvg=`echo "${intscale}; $expTotal / $expCount" | bc`

    tmpESNumr=0
    tmpESDenr=0
    tmpENumr=0
    tmpEDenr=0
    while read line ; do
        firstone=`echo $line | awk '{print $1}' | sed 's/ //g'| sed 's/\n//g'`
        exp_dG=`echo $line | awk '{print $2}' | sed 's/ //g'| sed 's/\n//g'`

        firstName=`grep "$firstone\b" ${PROCFILE} | awk '{print $1}'| sed 's/\n//g'`
        firstmtE=`grep "$firstone\b" ${PROCFILE} | awk '{print $6}'| sed 's/\n//g'`
        firstmtES=`grep "$firstone\b" ${PROCFILE} | awk '{print $2}'| sed 's/\n//g'`

        tmpESNumr=`echo "${intscale}; $tmpESNumr + (( $firstmtES - $sESAvg ) * ( $exp_dG - $expAvg ))" | bc`
        tmpESDenr=`echo "${intscale}; $tmpESDenr + (( $firstmtES - $sESAvg ) * ( $firstmtES - $sESAvg ))" | bc`

        tmpENumr=`echo "${intscale}; $tmpENumr + (( $firstmtE - $sEAvg ) * ( $exp_dG - $expAvg ))" | bc`
        tmpEDenr=`echo "${intscale}; $tmpEDenr + (( $firstmtE - $sEAvg ) * ( $firstmtE - $sEAvg ))" | bc`
        
        tmpExpDenr=`echo "${intscale}; $tmpEDenr + (( $exp_dG - $expAvg ) * ( $exp_dG - $expAvg ))" | bc`

    done < "../../singles.txt"
    
    Ses=`echo "${intscale}; sqrt ( $tmpESDenr / ( $sESCount - 1 ) )" | bc`
    Se=`echo "${intscale}; sqrt ( $tmpEDenr / ( $sECount - 1 ) )" | bc`
    Sexp=`echo "${intscale}; sqrt ( $tmpExpDenr / ( $expCount - 1 ) )" | bc`

    if [ ! -z ${SCALE+x} ] ; then
        mES=`echo "${intscale}; $tmpESNumr / $tmpESDenr" | bc`
        mE=`echo "${intscale}; $tmpENumr / $tmpEDenr" | bc`
        bES=`echo "${intscale}; $expAvg - ( $mES * $sESAvg )" | bc`
        bE=`echo "${intscale}; $expAvg - ( $mE * $sEAvg )" | bc`
    fi

    echo "$dir : Scaling Used: E: y = $mE x + $bE | ES: y = $mES x + $bES "

# process the singles (i.e. dG's)
    sESCount=0
    sECount=0
    sES_RunningABSdiff=0
    sE_RunningABSdiff=0
    sES_RunningSQdiff=0
    sE_RunningSQdiff=0
    if [ ! -z ${PRINTDETAILS+x} ] ; then
        echo "Lig1 exp_dG E_${dir} ES_${dir}"
    fi
    while read line ; do

        firstone=`echo $line | awk '{print $1}' | sed 's/ //g'| sed 's/\n//g'`
        exp_dG=`echo $line | awk '{print $2}' | sed 's/ //g'| sed 's/\n//g'`
                
        firstName=`grep "$firstone\b" ${PROCFILE} | awk '{print $1}'| sed 's/\n//g'`
        firstmtE=`grep "$firstone\b" ${PROCFILE} | awk '{print $6}'| sed 's/\n//g'`
        firstmtES=`grep "$firstone\b" ${PROCFILE} | awk '{print $2}'| sed 's/\n//g'`
# scale according to input
        firstmtE=`echo "((${mE} * $firstmtE)+(${bE}))" | bc`
        firstmtES=`echo "((${mES} * $firstmtES)+(${bES}))" | bc`

        printline=0
        mtESdG=$firstmtES
        absMtESdG_exp_dG="N/A"
        sqAbsMtESdG_exp_dG="N/A"
        if checkNumber $firstmtES ; then
            mtESdG_exp_dG=`echo "$mtESdG - $exp_dG" | bc`
            absMtESdG_exp_dG=${mtESdG_exp_dG#-}
            sqAbsMtESdG_exp_dG=`echo "$absMtESdG_exp_dG * $absMtESdG_exp_dG" | bc`
            sESCount=`echo "$sESCount + 1" | bc`
            sES_RunningABSdiff=`echo "$sES_RunningABSdiff + $absMtESdG_exp_dG" | bc`
            sES_RunningSQdiff=`echo "$sES_RunningSQdiff + $sqAbsMtESdG_exp_dG" | bc`
            sES_MUE=`echo "scale=5; $sES_RunningABSdiff / $sESCount" | bc`
            sES_RMSE=`echo "scale=5; sqrt($sES_RunningSQdiff / $sESCount)" | bc`
            printline=1
            
            sESTotal=`echo "${intscale}; $sESTotal + $mtESdG" | bc`                
        else
            echo "ERROR (MTScoreES) IN: $dir / $firstone"      
        fi
    
        mtEdG=$firstmtE
        absMtEdG_exp_dG="N/A"
        sqAbsMtEdG_exp_dG="N/A"
        if checkNumber $firstmtE ; then
            mtEdG_exp_dG=`echo "$mtEdG - $exp_dG" | bc`
            absMtEdG_exp_dG=${mtEdG_exp_dG#-}
            sqAbsMtEdG_exp_dG=`echo "$absMtEdG_exp_dG * $absMtEdG_exp_dG" | bc`
            sECount=`echo "$sECount + 1" | bc`
            sE_RunningABSdiff=`echo "$sE_RunningABSdiff + $absMtEdG_exp_dG" | bc`
            sE_RunningSQdiff=`echo "$sE_RunningSQdiff + $sqAbsMtEdG_exp_dG" | bc`
            sE_MUE=`echo "scale=5; $sE_RunningABSdiff / $sECount" | bc`
            sE_RMSE=`echo "scale=5; sqrt($sE_RunningSQdiff / $sECount)" | bc`
            printline=1
        else
            echo "ERROR (MTScoreE) IN: $dir / $firstone"      
        fi
        
        if [ ! -z ${PRINTDETAILS+x} ] ; then
            echo "${firstName} ${exp_dG} ${mtEdG} ${mtESdG}"
        fi
    done < "../../singles.txt"
    
#     echo "sE_MUE: $sE_MUE  sES_MUE: $sES_MUE"
#     echo "sE_RMSE: $sE_RMSE sES_RMSE: $sES_RMSE"

# process the pairs (i.e. ddG's)
    ESCount=0
    ECount=0
    ES_RunningABSdiff=0
    E_RunningABSdiff=0
    ES_RunningSQdiff=0
    E_RunningSQdiff=0

    if [ ! -z ${PRINTDETAILS+x} ] ; then
        echo "Lig1 Lig2 exp_ddG MTScoreE_ddG MTScoreES_ddG"
    fi
    
    while read line ; do

        firstone=`echo $line | awk '{print $1}' | sed 's/ //g'| sed 's/\n//g'`
        secondone=`echo $line | awk '{print $2}' | sed 's/ //g'| sed 's/\n//g'`
        exp_ddG=`echo $line | awk '{print $3}' | sed 's/ //g'| sed 's/\n//g'`
    
    #    grep "$firstone\b" ${PROCFILE}
    
        firstName=`grep "$firstone\b" ${PROCFILE} | awk '{print $1}'| sed 's/\n//g'`
        firstmtE=`grep "$firstone\b" ${PROCFILE} | awk '{print $6}'| sed 's/\n//g'`
        firstmtES=`grep "$firstone\b" ${PROCFILE} | awk '{print $2}'| sed 's/\n//g'`
# scale according to input
        firstmtE=`echo "((${mE} * $firstmtE)+(${bE}))" | bc`
        firstmtES=`echo "((${mES} * $firstmtES)+(${bES}))" | bc`

        secondName=`grep "$secondone\b" ${PROCFILE} | awk '{print $1}'| sed 's/\n//g'`
        secondmtE=`grep "$secondone\b" ${PROCFILE} | awk '{print $6}'| sed 's/\n//g'`
        secondmtES=`grep "$secondone\b" ${PROCFILE} | awk '{print $2}'| sed 's/\n//g'`
# scale according to input
        secondmtE=`echo "((${mE} * $secondmtE) + (${bE}))" | bc`
        secondmtES=`echo "((${mES} * $secondmtES) + (${bES}))" | bc`

        printline=0
        mtESddG="N/A"
        absMtESddG_exp_ddG="N/A"
        sqAbsMtESddG_exp_ddG="N/A"
        if checkNumber $firstmtES ; then
            if checkNumber $secondmtES ; then
                mtESddG=`echo "$secondmtES - $firstmtES" | bc`
                mtESddG_exp_ddG=`echo "$mtESddG - $exp_ddG" | bc`
                absMtESddG_exp_ddG=${mtESddG_exp_ddG#-}
                sqAbsMtESddG_exp_ddG=`echo "$absMtESddG_exp_ddG * $absMtESddG_exp_ddG" | bc`
                ESCount=`echo "$ESCount + 1" | bc`
                ES_RunningABSdiff=`echo "$ES_RunningABSdiff + $absMtESddG_exp_ddG" | bc`
                ES_RunningSQdiff=`echo "$ES_RunningSQdiff + $sqAbsMtESddG_exp_ddG" | bc`
                ES_MUE=`echo "scale=5; $ES_RunningABSdiff / $ESCount" | bc`
                ES_RMSE=`echo "scale=5; sqrt($ES_RunningSQdiff / $ESCount)" | bc`
                printline=1
            else
                echo "ERROR (MTScoreES) IN: $dir / $secondone"      
            fi
        else
            echo "ERROR (MTScoreES) IN: $dir / $firstone"      
        fi
    
        mtEddG="N/A"
        absMtEddG_exp_ddG="N/A"
        sqAbsMtEddG_exp_ddG="N/A"
        if checkNumber $firstmtE ; then
            if checkNumber $secondmtE ; then
                mtEddG=`echo "$secondmtE - $firstmtE" | bc`
                mtEddG_exp_ddG=`echo "$mtEddG - $exp_ddG" | bc`
                absMtEddG_exp_ddG=${mtEddG_exp_ddG#-}
                sqAbsMtEddG_exp_ddG=`echo "$absMtEddG_exp_ddG * $absMtEddG_exp_ddG" | bc`
                ECount=`echo "$ECount + 1" | bc`
                E_RunningABSdiff=`echo "$E_RunningABSdiff + $absMtEddG_exp_ddG" | bc`
                E_RunningSQdiff=`echo "$E_RunningSQdiff + $sqAbsMtEddG_exp_ddG" | bc`
                E_MUE=`echo "scale=5; $E_RunningABSdiff / $ECount" | bc`
                E_RMSE=`echo "scale=5; sqrt($E_RunningSQdiff / $ECount)" | bc`
                printline=1
            else
                echo "ERROR (MTScoreE) IN: $dir / $secondone"      
            fi
        else
            echo "ERROR (MTScoreE) IN: $dir / $firstone"      
        fi
            
        if [ ! -z ${PRINTDETAILS+x} ] ; then
            echo "${secondName} ${firstName} ${exp_ddG} ${mtEddG} ${mtEddG_exp_ddG} ${mtESddG} ${mtESddG_exp_ddG} "
        fi
        
    done < "../../pairs.txt"

#    echo "$dir : E_MUE/RMSE: $ECount $E_MUE $E_RMSE | ES_MUE/RMSE: $ESCount $ES_MUE $ES_RMSE | avgtime(min): $averagetime"

#     echo "E_MUE:  $E_MUE   ES_MUE:  $ES_MUE"
#     echo "E_RMSE: $E_RMSE   ES_RMSE: $ES_RMSE"
#     echo "avgtime(min): $averagetime"

    baseversion=`echo $fileName | awk -F "." '{print $2}'`
    if [ -e OUT.${baseversion} ] ; then
        calccount=`grep RUNNING OUT.${baseversion} | wc -l`
        taveragetime=`cat OUT.${baseversion} | grep user | awk '{print $2}' | awk -F "m" '{print ($1*60+$2)/60}' | paste -sd+ | awk -v calccount=${calccount} '{print "scale=2; ("$1") / "calccount}' | bc -l`
#        echo "testing: $calccount $sESCount | $sECount | ${expCount} | $taveragetime"
    fi

    echo "FILE, averagetime, averagetime sE_MUE, sE_RMSE, sES_MUE, sES_RMSE, E_MUE, E_RMSE, ES_MUE, ES_RMSE"
    echo "$fileName, $averagetime, $taveragetime, $sE_MUE, $sE_RMSE, $sES_MUE, $sES_RMSE, $E_MUE, $E_RMSE, $ES_MUE, $ES_RMSE"

done

