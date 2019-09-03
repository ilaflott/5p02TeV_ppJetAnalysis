#!/bin/bash

R=$1
condorDir=$2
dateDataIn=$3
etabin=$4
filename=$5
outtag=$6
targtrig=$7
reftrig=$8
#filenameLowJets=$5
#filenameJet80=$6
#outtag=$7
#targtrig=$8
#reftrig=$9

echo ""
echo "Running printPlots_jetTrigEff..."
echo ""

#./printPlots_jetTrigEff.exe ${condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${dateDataIn}_jetTrigEff_${etabin}/ ${filenameLowJets} ${filenameJet80} HPtJetTrig_${outtag} ${targtrig} ${reftrig}
./printPlots_jetTrigEff.exe ${condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${dateDataIn}_jetTrigEff_${etabin}/ ${filename} HPtJetTrig_${outtag} ${targtrig} ${reftrig}


return
