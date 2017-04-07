#!/bin/bash


##### R=4, up to 2.4
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0  0.0 2.4
#sleep 5s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 0.0 2.4
#sleep 5s
#
##### R=4, 2.4 to  2.7
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0  2.4 2.7
#sleep 5s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 2.4 2.7
#sleep 5s
#
##### R=4, 2.7 to 3.0
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0  2.7 3.0
#sleep 5s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 2.7 3.0
#sleep 5s
#
##### R=4, 3.0 to 3.2
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0  3.0 3.2
#sleep 5s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 3.0 3.2
#sleep 5s
#
##### R=4, 3.2 to 4.7
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0  3.2 4.7
#sleep 5s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 3.2 4.7
##sleep 5s
#
##sleep 30s



##### R=4, up to 2.4
source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 0.0 2.4
sleep 5s

##### R=4, 2.4 to  2.7
source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 2.4 2.7
sleep 5s

##### R=4, 2.7 to  3.0
source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 2.7 3.0
sleep 5s

##### R=4, 2.7 to  3.0
source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 3.0 3.2
sleep 5s

##### R=4, 3.2 to  4.7
source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 3.2 4.7
sleep 5s


return