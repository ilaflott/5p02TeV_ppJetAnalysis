#!/bin/bash

echo ""
echo "compiling code(s)"
echo ""

#echo "...ppMC JERS..."
#rootcompile readForests_ppMC_JERS.C

echo "...ppMC jetPlots..."
rootcompile readForests_ppMC_jetPlots.C

echo "...ppData jetPlots..."
rootcompile readForests_ppData_jetPlots.C

echo ""
echo "done compiling code(s), submitting."
echo ""


## ppData and ppMC jetPlots
source run_readForests_jetPlots.sh 4

source askCondor.sh 1 1
sleep 30s

source askCondor.sh 11 60

### ppMC JERS
#source run_readForests_ppMC_JERS.sh 4
#source askCondor.sh 1 1
#sleep 0.5s
#
#source askCondor.sh 11 60


return
