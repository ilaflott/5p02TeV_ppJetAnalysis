#!/bin/bash

#### MC TEST SUBMIT
#
#rootcompile readForests_ppMC_JERS.C
#source condorSubmit_readForests.sh readForests_ppMC_JERS 1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0 2.0
#sleep 1s
#source askCondor.sh


## -------------------------------------------- #
#### OFFICIAL MC for jetRECO/analysis meeting ### NOTE: NO JET ID, NO UNFOLDING HISTS, GEN/RECO JET PT > 50,60
## -------------------------------------------- #
#
#rootcompile readForests_ppMC_JERS.C
#
#source condorSubmit_readForests.sh readForests_ppMC_JERS -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 0 1.0
#sleep 10s
#source condorSubmit_readForests.sh readForests_ppMC_JERS -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 1.0 2.0
#sleep 10s
#source condorSubmit_readForests.sh readForests_ppMC_JERS -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 2.0 3.0
#sleep 10s
#source condorSubmit_readForests.sh readForests_ppMC_JERS -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 3.2 4.7
##sleep 10s
#
#source askCondor.sh



# ---------------------------------------- #
### SEMIOFFICIAL MC for analysis meeting ### NOTE: NO JET ID, WITH UNFOLDING HISTS, GEN/RECO JET PT > 50,60
# ---------------------------------------- #

rootcompile readForests_ppMC_JERS.C

source condorSubmit_readForests.sh readForests_ppMC_JERS -1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0 1.0
sleep 10s
source condorSubmit_readForests.sh readForests_ppMC_JERS -1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 1.0 2.0
sleep 10s
source condorSubmit_readForests.sh readForests_ppMC_JERS -1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 2.0 3.0
sleep 10s
source condorSubmit_readForests.sh readForests_ppMC_JERS -1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 3.2 4.7
#sleep 10s

source askCondor.sh

