#!/bin/bash

R=$1
condorDir=$2
date_output=$3
etabin=$4
etabinOut=$5


echo ""
echo "Running ppMC JERS for R=${R}, date=${date_output}, and etabin ${etabin}"
echo ""

./printPlots_ppMC_JERS.exe ${condorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_JERS_${etabin} ${etabinOut}

return