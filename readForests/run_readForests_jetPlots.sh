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

source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 7 0 filelists/5p02TeV_HighPtJet80_forests.txt     ${R} PF 0  0.0  0.5
sleep 0.5													                      
source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 7 0 filelists/5p02TeV_HighPtLowerJets_forests.txt ${R} PF 0  0.0  0.5
sleep 0.5s													                      


#### semiOfficial Py8
echo ""
echo "submitting Py8 ppMC job(s)"
echo ""

source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 15 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0  0.0  0.5
sleep 0.5s

echo ""
echo "done submitting jetPlots job(s)"
echo ""

return