#!/bin/bash

R=$1
dataCondorDir=$2
dataDate=$3
dataEtaBin=$4
MCCondorDir=$5
MCDate=$6
MCEtaBin=$7
outtag=$8

echo ""
echo "Running w/ Jet80+LowerJets PDs..."
echo ""


#./printPlots_jetMultPlots.exe   ${condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${dateDataIn}_jetMult_${dataEtaBin}/ HPtJetTrig_semiOffPy8_${outtag}

./printPlots_jetMultPlots.exe   ${dataCondorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${dataDate}_jetMult_${dataEtaBin}/ ${MCCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${MCDate}_jetMult_${MCEtaBin}/ HPtJetTrig_semiOffPy8_${outtag}


return
