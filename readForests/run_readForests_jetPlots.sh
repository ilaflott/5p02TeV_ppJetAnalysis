#!/bin/bash

echo ""
echo "compiling readForests_ppData_jetPlots.C"
echo ""
rootcompile readForests_ppData_jetPlots.C

##### TEST ppData jetPlots
#
#echo ""
#echo "submitting test job(s)"
#echo ""
#source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 2 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0     0.0 0.5
#source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 2 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 0.0 0.5
#
#source askCondor.sh 1 3

#### Jet80, then LowerJets
echo ""
echo "submitting job(s) for real!"
echo ""

source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 7 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0     0.0 0.5
sleep 1s
source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 7 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 0.0 0.5
sleep 1s
source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 7 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0     0.5 1.0
sleep 1s
source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 7 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 0.5 1.0
sleep 1s
source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 7 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0     1.0 1.5
sleep 1s
source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 7 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 1.0 1.5
sleep 1s
source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 7 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0     1.5 2.0
sleep 1s
source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 7 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 1.5 2.0
sleep 1s


echo ""
echo "condor says..."
echo ""
source askCondor.sh 2 5



### TEST Official Py 8

#source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 0.0 1.0
#source askCondor.sh

### Official Py 8

#source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 0.0 1.0
#
#source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 1.0 2.0
#
#source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 2.0 3.0
#
#source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 3.2 4.7
#


echo ""
echo "compiling readForests_ppMC_jetPlots.C"
echo ""
rootcompile readForests_ppMC_jetPlots.C

### TEST semiOfficial Py8
#echo ""
#echo "submitting test job"
#echo ""

#source condorSubmit_readForests.sh readForests_ppMC_jetPlots 1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0.0 0.5
#source askCondor.sh 1 3

### semiOfficial Py8
echo ""
echo "submitting job(s) for real"
echo ""

source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 15 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0.0 0.5
sleep 1s
source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 15 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0.5 1.0
sleep 1s
source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 15 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 1.0 1.5
sleep 1s
source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 15 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 1.5 2.0
sleep 1s

echo ""
echo "done submitting, asking condor..."
echo ""

source askCondor.sh 2 5

echo ""
#echo "done."
echo "done. Submitting JERS jobs now too!"
echo ""

source run_readForests_ppMC_JERS.sh

return