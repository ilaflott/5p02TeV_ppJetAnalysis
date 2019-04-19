#!/bin/bash

if [[ $# -eq 9 ]]
then 

    R=$1
    condorDir=$2
    date_output=$3
    descIn=$4
    descOut=$5
    JERplots=$6
    MCEffplots=$7
    absetabin=$8
    fitquant=$9


else
    echo ""
    echo "use error."
    echo "source run_printPlots_ppMC_JERS_etabin.sh <R> <condorDir> <date_output> <descIn> <descOut> <JERplots> <MCeffplots> <absetabin> <fitquant>"
    echo ""
    return
fi

echo ""
echo "Running ppMC JERS for R=${R}, date=${date_output}, and descIn ${descIn}"
echo ""

#./printPlots_ppMC_JERS.exe ${condorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_JERS_${descIn} ${descOut} $JERplots $MCEffplots
./printPlots_ppMC_JERS_etabin.exe ${condorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_JERS_${descIn} ${descOut} $JERplots $MCEffplots $absetabin $fitquant

return