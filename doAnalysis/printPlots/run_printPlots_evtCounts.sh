#!/bin/bash

R=$1
dataCondorDir=$2
dateDataIn=$3
dataEtaBin=$4
MCCondorDir=$5
dateMCIn=$6
MCEtaBin=$7
outtag=$8

echo ""
echo "Running w/ Jet80+LowerJets+PY8 PDs..."
echo ""
./printPlots_evtCounts.exe   ${dataCondorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${dateDataIn}_jetPlots_${dataEtaBin}/ ${MCCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${dateMCIn}_unf_${MCEtaBin}/   HPtJetTrig_semiOffPy8_${outtag}


return
