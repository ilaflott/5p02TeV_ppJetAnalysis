#!/bin/bash

echo ""
echo "compiling code(s)"
echo ""

#echo "...ppMC JERS..."
#rootcompile readForests_ppMC_JERS.C

#echo "...ppMC jetPlots..."
#rootcompile readForests_ppMC_jetPlots.C

echo "...ppData jetPlots..."
rootcompile readForests_ppData_jetPlots.C

echo ""
echo "done compiling code(s), submitting."
echo ""

## ppData and ppMC jetPlots
source run_readForests_jetPlots.sh 4 "0.0" "2.0"
#sleep 60s
#source run_readForests_jetPlots.sh 4 "2.0" "3.0"
#sleep 60s
#source run_readForests_jetPlots.sh 4 "3.2" "4.7"
#sleep 60s

echo ""
echo "done running code(s)."
echo ""

source askCondor.sh 100 60

return







