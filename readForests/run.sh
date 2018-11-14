#!/bin/bash

subm_ppData_jetPlots=0
subm_ppMC_jetPlots=0
subm_ppMC_JERS=1
subm_ppData_jetMult=0
subm_ppMC_jetMult=0


echo ""
echo "compiling code(s) + submitting jobs"
echo ""

if [[ $subm_ppData_jetPlots -eq 1 ]]
then
    echo "...ppData jetPlots..."
    rootcompile readForests_ppData_jetPlots.C
    source run_readForests_jetPlots.sh 4 "0.0" "2.0" "ppData"
#    source run_readForests_jetPlots.sh 4 "0.0" "0.5" "ppData"
#    source run_readForests_jetPlots.sh 4 "0.5" "1.0" "ppData"
#    source run_readForests_jetPlots.sh 4 "1.0" "1.5" "ppData"
#    source run_readForests_jetPlots.sh 4 "1.5" "2.0" "ppData"
fi

if [[ $subm_ppMC_jetPlots -eq 1 ]]
then
    echo "...ppMC jetPlots..."
    rootcompile readForests_ppMC_jetPlots.C
    source run_readForests_jetPlots.sh 4 "0.0" "2.0" "ppMC"
#    source run_readForests_jetPlots.sh 4 "0.0" "0.5" "ppMC"
#    source run_readForests_jetPlots.sh 4 "0.5" "1.0" "ppMC"
#    source run_readForests_jetPlots.sh 4 "1.0" "1.5" "ppMC"
#    source run_readForests_jetPlots.sh 4 "1.5" "2.0" "ppMC"
fi

if [[ $subm_ppMC_JERS -eq 1 ]]
then
    echo "...ppMC JERS..."
    rootcompile readForests_ppMC_JERS.C
    source run_readForests_JERS.sh 4 "0.0" "2.0"
#    source run_readForests_JERS.sh 4 "0.0" "0.5"
#    source run_readForests_JERS.sh 4 "0.5" "1.0"
#    source run_readForests_JERS.sh 4 "1.0" "1.5"
#    source run_readForests_JERS.sh 4 "1.5" "2.0"
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

askCondor 100 100


return

