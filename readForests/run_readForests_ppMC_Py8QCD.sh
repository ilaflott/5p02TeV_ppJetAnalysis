#!/bin/bash

#### TEST READFORESTS
#./readForests_ppMC_JetQA_RapBins.exe
#./readForests_ppMC_Unfold.exe
#./readForests_ppMC_MCEff.exe
#./readForests_ppMC_MCJEC_JetID0.exe
#./readForests_ppMC_MCJEC_JetID1.exe

#### TEST CONDOR SUBMIT
## read ak3PFjetAnalyzer
#source condorSubmit_readForests.sh readForests_ppMC_JetQA_RapBins 1 1 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_Unfold        1 1 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_MCEff         1 1 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_MCJEC_JetID0  1 1 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
source condorSubmit_readForests.sh readForests_ppMC_MCJEC_JetID1  1 1 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0


##### CONDOR SUBMIT
#### read ak3PFjetAnalyzer
#source condorSubmit_readForests.sh readForests_ppMC_JetQA_RapBins -1 10 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_Unfold        -1 10 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_MCEff         -1 10 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
##source condorSubmit_readForests.sh readForests_ppMC_MCJEC_JetID0  -1 10 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
##source condorSubmit_readForests.sh readForests_ppMC_MCJEC_JetID1  -1 10 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
##
### read ak4PFjetAnalyzer
#source condorSubmit_readForests.sh readForests_ppMC_JetQA_RapBins -1 10 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_Unfold        -1 10 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_MCEff         -1 10 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
##source condorSubmit_readForests.sh readForests_ppMC_MCJEC_JetID0  -1 10 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
##source condorSubmit_readForests.sh readForests_ppMC_MCJEC_JetID1  -1 10 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
##
### read ak5PFjetAnalyzer
##
