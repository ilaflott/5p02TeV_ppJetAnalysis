#!/bin/bash


##### SEMI PRIVATE MC
#
#### test job
##source condorSubmit_readForests.sh readForests_ppMC_JERS 1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0 2.0
#
#### all regions
##source condorSubmit_readForests.sh readForests_ppMC_JERS -1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0 4.7
#
### jetID eta bins
#source condorSubmit_readForests.sh readForests_ppMC_JERS -1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0 2.4
#source condorSubmit_readForests.sh readForests_ppMC_JERS -1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 2.4 2.7
#source condorSubmit_readForests.sh readForests_ppMC_JERS -1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 2.7 3.0
#source condorSubmit_readForests.sh readForests_ppMC_JERS -1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 3.0 3.2
#source condorSubmit_readForests.sh readForests_ppMC_JERS -1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 3.2 4.7
#
#### barrel/endcap only
##source condorSubmit_readForests.sh readForests_ppMC_JERS -1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0 2.0






#### OFFICIAL MC

### test job
#source condorSubmit_readForests.sh readForests_ppMC_JERS 1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 0 1.0

### jetID eta bins
source condorSubmit_readForests.sh readForests_ppMC_JERS -1 50 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 0 2.4
source condorSubmit_readForests.sh readForests_ppMC_JERS -1 50 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 2.4 2.7
source condorSubmit_readForests.sh readForests_ppMC_JERS -1 50 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 2.7 3.0
source condorSubmit_readForests.sh readForests_ppMC_JERS -1 50 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 3.0 3.2
source condorSubmit_readForests.sh readForests_ppMC_JERS -1 50 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 3.2 4.7

### 0-1, central barrel
source condorSubmit_readForests.sh readForests_ppMC_JERS -1 50 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt 4 PF 0 0 1.0




