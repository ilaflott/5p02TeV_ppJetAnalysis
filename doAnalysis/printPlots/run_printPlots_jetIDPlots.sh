#!/bin/bash

R=$1
condorDir=$2
dateIn=$3
etabin_jetID=$4
etabin_nojetID=$5
outtag=$6
type=$7

echo ""
echo "Running w/ Jet80+LowerJets PDs..."
echo ""

if [ "${type}" = "ppData" ]
then
#    ./printPlots_jetIDPlots.exe   ${condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${dateIn}_jetPlots_${etabin}/       HPtJetTrig_${outtag}
    ./printPlots_jetIDPlots.exe   ${condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${dateIn}_jetPlots_${etabin_jetID}/  ${condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${dateIn}_jetPlots_${etabin_nojetID}/     HPtJetTrig_${outtag}
elif [ "${type}" = "ppMC" ]
then
#    ./printPlots_jetIDPlots.exe   ${condorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${dateIn}_jetPlots_${etabin}/      Py8_CUETP8M1_${outtag}
    ./printPlots_jetIDPlots.exe   ${condorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${dateIn}_jetPlots_${etabin_jetID}/ ${condorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${dateIn}_jetPlots_${etabin_nojetID}/ Py8_CUETP8M1_${outtag}
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




