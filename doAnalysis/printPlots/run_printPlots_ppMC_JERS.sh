#!/bin/bash

R=$1
date_output=$2
etabin=$3
etabinOut=$4


echo ""
echo "Running ppMC JERS for R=${R}, date=${date_output}, and etabin ${etabin}"
echo ""

./printPlots_ppMC_JERS.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_JERS_${etabin} ${etabinOut}

return