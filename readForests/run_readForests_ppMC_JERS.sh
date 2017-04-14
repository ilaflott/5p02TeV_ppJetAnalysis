#!/bin/bash



#source condorSubmit_readForests.sh readForests_ppMC_JERS_v0 1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 0 2.4
#sleep 1s
#source askCondor.sh

##for jetRECO meeting
#rootcompile readForests_ppMC_JERS_v0.C
#
##### OFFICIAL MC
#source condorSubmit_readForests.sh readForests_ppMC_JERS_v0 -1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 0 1.0
#sleep 60s
#
#source condorSubmit_readForests.sh readForests_ppMC_JERS_v0 -1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 0 2.0
#sleep 60s
##source askCondor.sh
#
###for analysis meeting
rootcompile readForests_ppMC_JERS.C
#
##### OFFICIAL MC
source condorSubmit_readForests.sh readForests_ppMC_JERS -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 0.0 1.0
#source condorSubmit_readForests.sh readForests_ppMC_JERS_v1 -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 0.0 1.0
#sleep 60s
sleep 1s
#
#source condorSubmit_readForests.sh readForests_ppMC_JERS_v1 -1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 0.0 2.0
#sleep 60s
source askCondor.sh


#### SEMI OFFICIAL MC
#source condorSubmit_readForests.sh readForests_ppMC_JERS_v1 -1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0.0 1.0
#source condorSubmit_readForests.sh readForests_ppMC_JERS_v1 -1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0.0 2.0



