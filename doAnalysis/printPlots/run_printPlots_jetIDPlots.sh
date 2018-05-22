#!/bin/bash

R=$1
condorDir=$2
dateIn=$3
etabin=$4
outtag=$5
type=$6

echo ""
echo "Running w/ Jet80+LowerJets PDs..."
echo ""

if [ "${type}" = "ppData" ]
then
    ./printPlots_jetIDPlots.exe   ${condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${dateIn}_jetPlots_${etabin}/                 HPtJetTrig_${outtag}
elif [ "${type}" = "ppMC" ]
then
    ./printPlots_jetIDPlots.exe   ${condorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${dateIn}_jetPlots_${etabin}/      Py8_CUETP8M1_${outtag}
else
    echo ""
    echo "data type ${type} not known. exit"
    echo ""
fi

return

















#echo ""
#echo "Running w/ Jet80 PDs..."
#echo ""
#./printPlots_jetPlots.exe   ppData_HighPtJetTrig_ak${R}PFJets_${date_output}_jetPlots_${etabin}/ppData_HighPtJet80_ak${R}PFJets_${date_output}_jetPlots_${etabin}/ ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_jetPlots_${etabin}/   HPtJet80_semiOffPy8_${etabinOut}

#echo ""
#echo "Running w/ LowerJets PDs..."
#echo ""
#./printPlots_jetPlots.exe   ppData_HighPtJetTrig_ak${R}PFJets_${date_output}_jetPlots_${etabin}/ppData_HighPtLowerJets_ak${R}PFJets_${date_output}_jetPlots_${etabin}/ ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_jetPlots_${etabin}/   HPtLowerJets_semiOffPy8_${etabinOut}




