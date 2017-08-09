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
#source run_readForests_jetPlots_noHLT100.sh 4

echo ""
echo "done running code(s)."
echo ""

return







