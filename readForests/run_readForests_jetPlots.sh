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

source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt     ${R} PF 0  0.0  0.5
source askCondor.sh 2 10
sleep 60s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 5 0 filelists/5p02TeV_HighPtLowerJets_forests.txt ${R} PF 0  0.0  0.5
#source askCondor.sh 2 10
#sleep 60s													                      

###### semiOfficial Py8
#echo ""
#echo "submitting Py8 ppMC job(s)"
#echo ""
#
#source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 13 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0  0.0  0.5
#source askCondor.sh 2 10
#sleep 60s



echo ""
echo "done submitting jetPlots job(s)"
echo ""

source askCondor.sh 20 30

return