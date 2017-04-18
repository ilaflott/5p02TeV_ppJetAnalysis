#!/bin/bash

rootcompile readForests_ppData_jetPlots.C

#### TEST ppData jetPlots
#source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 2 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0  0.0 1.0
#source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 2 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 0.0 1.0
#source askCondor.sh

### Jet80, then LowerJets
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 8 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0  0.0 1.0
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 8 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 0.0 1.0
#sleep 1s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 8 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0  1.0 2.0
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 8 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 1.0 2.0
#sleep 1s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 8 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0  2.0 3.0
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 8 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 2.0 3.0
#sleep 1s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 8 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0  3.2 4.7
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 8 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 3.2 4.7
#sleep 1s
#source askCondor.sh

### TEST Official Py 8
#source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 0.0 1.0
#source askCondor.sh

### Official Py 8
#source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 0.0 1.0
#sleep 5s
#source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 1.0 2.0
#sleep 5s
#source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 2.0 3.0
#sleep 5s
#source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 3.2 4.7
#sleep 5s


## TEST semiOfficial Py8
source condorSubmit_readForests.sh readForests_ppMC_jetPlots 1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0.0 1.0

### semiOfficial Py8
#source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0.0 1.0
#sleep 5s
#source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 1.0 2.0
#sleep 5s
#source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 2.0 3.0
#sleep 5s
#source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 3.2 4.7
#sleep 5s


return