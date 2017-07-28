#!/bin/bash

R=$1
date_output="07-25-17"

echo ""
echo "compiling printPlots_jetPlots"
echo ""
rootcompile printPlots_jetPlots_combine.C

echo ""
echo "compilation done..."
echo ""
sleep 1s

echo ""
echo "Running w/ Jet80+LowerJets PDs..."
echo ""
./printPlots_jetPlots_combine.exe   ppData_HighPtJetTrig_ak${R}PFJets_${date_output}_jetPlots_0.0eta0.5/ ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_jetPlots_0.0eta0.5/   HPtJetTrig_semiOffPy8_0eta05

#echo ""
#echo "Running w/ Jet80 PDs..."
#echo ""
#./printPlots_jetPlots_combine.exe   ppData_HighPtJetTrig_ak${R}PFJets_${date_output}_jetPlots_0.0eta0.5/ppData_HighPtJet80_ak${R}PFJets_${date_output}_jetPlots_0.0eta0.5/ ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_jetPlots_0.0eta0.5/   HPtJet80_semiOffPy8_0eta05



#echo ""
#echo "Running w/ LowerJets PDs..."
#echo ""
#./printPlots_jetPlots_combine.exe   ppData_HighPtJetTrig_ak${R}PFJets_${date_output}_jetPlots_0.0eta0.5/ppData_HighPtLowerJets_ak${R}PFJets_${date_output}_jetPlots_0.0eta0.5/ ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_jetPlots_0.0eta0.5/   HPtLowerJets_semiOffPy8_0eta05



sleep 1s
