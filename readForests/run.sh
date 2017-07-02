#!/bin/bash

#echo ""
#echo "compiling code(s)"
#echo ""
#
#echo "...ppMC JERS..."
#rootcompile readForests_ppMC_JERS.C
#
#echo "...ppMC jetPlots..."
#rootcompile readForests_ppMC_jetPlots.C
#
#echo "...ppData jetPlots..."
#rootcompile readForests_ppData_jetPlots.C

#echo ""
#echo "done compiling code(s), submitting."
#echo ""


echo ""
echo "submitting R=3 now."
echo ""

## ppData and ppMC jetPlots
source run_readForests_jetPlots.sh 3
source askCondor.sh 1 1
sleep 0.5s

## ppMC JERS
source run_readForests_ppMC_JERS.sh 3
source askCondor.sh 1 1
sleep 0.5s


echo ""
echo "submission for R=3 done. askingCondor."
echo ""

source askCondor.sh 11 60



echo ""
echo "submitting R=4 now."
echo ""

## ppData and ppMC jetPlots
source run_readForests_jetPlots.sh 4
source askCondor.sh 1 1
sleep 0.5s

## ppMC JERS
source run_readForests_ppMC_JERS.sh 4
source askCondor.sh 1 1
sleep 0.5s


echo ""
echo "submission for R=4 done. askingCondor."
echo ""

source askCondor.sh 11 60


echo ""
echo "submitting R=5 now."
echo ""

## ppData and ppMC jetPlots
source run_readForests_jetPlots.sh 5
source askCondor.sh 1 1
sleep 0.5s

## ppMC JERS
source run_readForests_ppMC_JERS.sh 5
source askCondor.sh 1 1
sleep 0.5s


echo ""
echo "submission for R=5 done. askingCondor."
echo ""

source askCondor.sh 11 60


return
