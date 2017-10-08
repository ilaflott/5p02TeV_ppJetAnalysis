#!/bin/bash

echo ""
echo "compiling unfoldMCSpectra"
echo ""
rooUnfoldCompile unfoldMCSpectra.C

echo ""
echo "done compiling unfoldMCSpectra"
echo ""
sleep 1s

#source run_unfoldMCSpectra.sh <R> <condorDir>          <date_output> <etaBin>    <outputTag> <kReg> <useSimpleBinning>
#source run_unfoldMCSpectra.sh "4" "10.02.17_outputCondor" "10-02-17" "0.0eta2.0" "00eta20_10.2.17" 10 0
#source run_unfoldMCSpectra.sh "4" "10.02.17_outputCondor" "10-02-17" "0.0eta2.0" "00eta20_10.2.17" 6 1
source run_unfoldMCSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.3.17" 6 0
source run_unfoldMCSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.3.17" 6 1

echo ""
echo "done unfolding MC."
echo ""

return

