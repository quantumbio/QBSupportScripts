#!/bin/bash

    export ENVOY_URL='http://$(hostname).local:8090'
    export ENVOY_UI_URL='http://localhost:8091'
    export DATA_HEADER='Content-Type: application/json'
        
    declare -a COMM_LIST
    TESTDATA=$( curl -isS ${ENVOY_URL} 2>&1 >/dev/null )
    if [[ $TESTDATA == *"Connection refused"* ]]; then
        echo "ERROR: $TESTDATA"
        echo "ERROR: Envoy is not running at $ENVOY_URL. Please start Envoy and try again."
        echo "    Visit https://www.pharma.gridmarkets.com/setup for more information or to download Envoy."
        exit 1
    else
        isFail=0
        AUTHDATA=$( curl -s  ${ENVOY_URL}/auth )
        echo ${AUTHDATA}
#        username=`echo ${AUTHDATA} | jq '.Username'`
        username=$( echo ${AUTHDATA} | jq '.Username' || echo "ERROR" )
        [[ ${username} == "ERROR" ]] && isFail=1 && echo "ERROR USERNAME: ${AUTHDATA}"
        
        AUTHDATA=$( curl -s  ${ENVOY_URL}/user-info )
        echo ${AUTHDATA}
#        credits_available=`echo ${AUTHDATA} | jq '.data.credits_available'`
        credits_available=$(echo ${AUTHDATA} | jq '.data.credits_available' || echo "ERROR" )
        [[ ${credits_available} == "ERROR" ]] && isFail=1 && echo "ERROR CREDITS: ${AUTHDATA}"
                
        AUTHDATA=$( curl -s  ${ENVOY_URL}/info )
        echo ${AUTHDATA}
#        envoy_version=`echo ${AUTHDATA} | jq '.version'`
        envoy_version=$( echo ${AUTHDATA} | jq '.version' || echo "ERROR" )
        [[ ${envoy_version} == "ERROR" ]] && isFail=1 && echo "ERROR ENVOY: ${AUTHDATA}"

        AUTHDATA=$( curl -s  ${ENVOY_URL}/products )
        echo ${AUTHDATA}
#        export divcon_version=`echo ${AUTHDATA} | jq '[.[] | select(.app_type=="divcon")] | .[-1] | .version'`
        export divcon_version=0.DEV
        [[ ${divcon_version} == "ERROR" ]] && isFail=1 && echo "ERROR DIVCON: ${AUTHDATA}"

        AUTHDATA=$( curl -s -X POST ${ENVOY_URL}/machines -H "${DATA_HEADER}" -d "{\"operation\":\"simulation\",\"app\":\"DivCon\"}" )
        echo ${AUTHDATA}
#        export machine=`echo ${AUTHDATA} | jq '.data[0] | .id'`
        export machine=$( echo ${AUTHDATA} | jq '.data[0] | .id' || echo "ERROR" )
        [[ ${machine} == "ERROR" ]] && isFail=1 && echo "ERROR MACHINE: ${AUTHDATA}"

        if [[ ${isFail} == 1 ]]; then
            echo ""
            echo "ERROR: You are not logged into GridMarkets or there is a communication issue."
            echo "    Visit ${ENVOY_UI_URL} to verify your connection."
            echo "    Visit https://www.pharma.gridmarkets.com/ for support."
            echo ""
            exit 1
        fi
        
    fi
