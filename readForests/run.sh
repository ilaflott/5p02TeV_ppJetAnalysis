#!/bin/bash

## general note; 
## busy cluster (check w/ 'condor_q --allusers') --> use larger # of files per job --> longer run time but overall kinder during busy hours
## not-busy cluster --> can use smaller # of files per job --> shorter run time but occupies more cluster machines

subm_ppData_jetPlots=0
subm_ppData_jetTrigEff=0
subm_ppData_makeNTuple=0

subm_ppMC_jetPlots=1
subm_ppMC_JERS=0
subm_ppMC_unf=0
subm_ppMC_MCEff=0

sleep_between_subm=0
Nmin=60

#rarely used
subm_ppData_jetMult=0
subm_ppMC_jetMult=0
subm_ppData_findEvt=0


echo ""
echo "compiling code(s) + submitting jobs"
echo ""

if [[ $subm_ppData_jetPlots -eq 1 ]]
then
    echo "...ppData jetPlots..."
    rootcompile readForests_ppData_jetPlots.C
    source run_readForests_jetPlots.sh 4 "0.0" "2.5" "ppData"
    if [[ $sleep_between_subm -eq 1 ]]
    then
	#sleep ${Nmin}m
	askCondor 360 10
    fi
fi


if [[ $subm_ppData_jetTrigEff -eq 1 ]]
then
    echo "...ppData jetTrigEff..."
    rootcompile readForests_ppData_jetTrigEff.C
    source run_readForests_jetTrigEff.sh 4 "0.0" "5.1" "ppData"
    if [[ $sleep_between_subm -eq 1 ]]
	then
	sleep ${Nmin}m
    fi
fi

if [[ $subm_ppData_makeNTuple -eq 1 ]]
then
    echo "...ppData makeNTuple..."
    rootcompile readForests_ppData_makeNTuple.C
    source run_readForests_makeNTuple.sh 4 "0.0" "2.0"
    if [[ $sleep_between_subm -eq 1 ]]
	then
	sleep ${Nmin}m
    fi
fi

if [[ $subm_ppMC_jetPlots -eq 1 ]]
then
    echo "...ppMC jetPlots..."
    rootcompile readForests_ppMC_jetPlots.C
    source run_readForests_jetPlots.sh 4 "0.0" "2.5" "ppMC"
    if [[ $sleep_between_subm -eq 1 ]]
	then
	sleep ${Nmin}m
    fi
fi


if [[ $subm_ppMC_JERS -eq 1 ]]
then
    echo "...ppMC JERS..."
    rootcompile readForests_ppMC_JERS.C
    source run_readForests_JERS.sh 4 "0.0" "2.5"
    if [[ $sleep_between_subm -eq 1 ]]
	then
	sleep ${Nmin}m
    fi
fi


if [[ $subm_ppMC_unf -eq 1 ]]
then
    echo "...ppMC unf..."
    rootcompile readForests_ppMC_unf.C
    source run_readForests_unf.sh 4 "0.0" "2.5"
    if [[ $sleep_between_subm -eq 1 ]]
	then
	sleep ${Nmin}m
    fi
fi


if [[ $subm_ppMC_MCEff -eq 1 ]]
then
    echo "...ppMC MCEff..."
    rootcompile readForests_ppMC_MCEff.C
    source run_readForests_MCEff.sh 4 "0.0" "2.0"
    if [[ $sleep_between_subm -eq 1 ]]
	then
	sleep ${Nmin}m
    fi
fi





if [[ $subm_ppData_jetMult -eq 1 ]]
then
    echo "...ppData jetMult..."
    rootcompile readForests_ppData_jetMult.C
    source run_readForests_jetMult.sh 4 "0.0" "2.0" "ppData"
#    source run_readForests_jetMult.sh 4 "0.0" "0.5" "ppData"
#    source run_readForests_jetMult.sh 4 "0.5" "1.0" "ppData"
#    source run_readForests_jetMult.sh 4 "1.0" "1.5" "ppData"
#    source run_readForests_jetMult.sh 4 "1.5" "2.0" "ppData"
fi

if [[ $subm_ppMC_jetMult -eq 1 ]]
then
    echo "...ppMC jetMult..."
    rootcompile readForests_ppMC_jetMult.C
    source run_readForests_jetMult.sh 4 "0.0" "2.0" "ppMC"
#    source run_readForests_jetMult.sh 4 "0.0" "0.5" "ppMC"
#    source run_readForests_jetMult.sh 4 "0.5" "1.0" "ppMC"
#    source run_readForests_jetMult.sh 4 "1.0" "1.5" "ppMC"
#    source run_readForests_jetMult.sh 4 "1.5" "2.0" "ppMC"
fi
if [[ $subm_ppData_findEvt -eq 1 ]]
then
    echo "...ppData findEvt..."
    rootcompile readForests_ppData_findEvt.C
    source run_readForests_findEvt.sh 4 "0.0" "0.1"  # "ppData"
fi

askCondor 1000 10


return

