#!/bin/bash

R=$1
date_output="06-07-17"

#echo ""
#echo "compiling printPlots_jetPlots"
#echo ""
#rootcompile printPlots_jetPlots_combine.C
#
#echo ""
#echo "compilation done..."
#echo ""
#sleep 1s

echo ""
echo "Running..."
echo ""
./printPlots_jetPlots_combine.exe   ppData_HighPtJetTrig_ak${R}PFJets_${date_output}_jetPlots_0.0eta2.0/ ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_jetPlots_0.0eta2.0/   HPtJetTrig_semiOffPy8_0eta20

sleep 1s
