#!/bin/bash

echo ""
echo "compiling code(s)"
echo ""
#
#echo "...ppMC JERS..."
#rootcompile readForests_ppMC_JERS.C
#
#echo "...ppMC jetPlots..."
#rootcompile readForests_ppMC_jetPlots.C
#
echo "...ppData jetPlots..."
rootcompile readForests_ppData_jetPlots.C

echo ""
echo "done compiling code(s), submitting."
echo ""

#source run_readForests_jetPlots.sh 3
#source askCondor.sh 1 1
#sleep 60s
#
#source run_readForests_jetPlots.sh 4
#source askCondor.sh 1 1
#sleep 60s

source run_readForests_jetPlots.sh 5
source askCondor.sh 1 1
sleep 60s


#source run_readForests_ppMC_JERS.sh 3
#source askCondor.sh 1 1
#sleep 60s

echo ""
echo "submission/compilation done. askingCondor."
echo ""

source askCondor.sh 24 5

return
