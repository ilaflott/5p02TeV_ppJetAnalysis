#!/bin/bash

R=$1
condorDir=$2
dateDataIn=$3
etabin=$4
outtag=$5

echo ""
echo "Running w/ Jet80+LowerJets PDs..."
echo ""
./printPlots_jetTrig.exe ${condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${dateDataIn}_jetPlots_${etabin}/ HPtJetTrig_${outtag}


return
