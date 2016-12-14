#!/bin/bash

#source condorSubmit_readForests.sh readForests_ppMC_MCJEC_jtID0 -1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
#sleep 60s
#source condorSubmit_readForests.sh readForests_ppMC_MCJEC_jtIDL -1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
#sleep 60s
#source condorSubmit_readForests.sh readForests_ppMC_MCJEC_jtIDT -1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
#sleep 60s
sleep 120s

source condorSubmit_readForests.sh readForests_ppMC_MCJEC_jtID0 -1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
sleep 120s
source condorSubmit_readForests.sh readForests_ppMC_MCJEC_jtIDL -1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
sleep 120s
source condorSubmit_readForests.sh readForests_ppMC_MCJEC_jtIDT -1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
sleep 120s
