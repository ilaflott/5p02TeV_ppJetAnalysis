#!/bin/bash

if [ $# -eq 0 ]
then
    echo "give me a radius (3,4... etc.)"
    return
fi

R=$1

#### Jet80, then LowerJets
echo ""
echo "submitting jet80/lowerJets ppData job(s)"
echo ""

source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 20 0 filelists/5p02TeV_HighPtJet80_forests.txt     ${R} PF 0  0.0 4.7

echo ""
echo "done submitting ppData jetPlots Jet80 job(s)"
echo ""
source askCondor.sh 10 30

source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 20 0 filelists/5p02TeV_HighPtLowerJets_forests.txt ${R} PF 0  0.0 4.7

echo ""
echo "done submitting ppData jetPlots LowerJets job(s)"
echo ""
source askCondor.sh 10 30






##### semiOfficial Py8
echo ""
echo "submitting Py8 ppMC job(s)"
echo ""

source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 60 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0  0.0  4.7

echo ""
echo "done submitting ppMC jetPlots job(s)"
echo ""
source askCondor.sh 10 60

return