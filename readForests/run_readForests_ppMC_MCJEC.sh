#!/bin/bash

source condorSubmit_readForests.sh readForests_ppMC_MCJEC_jtIDTv2_eta50 -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
sleep 30s

source condorSubmit_readForests.sh readForests_ppMC_MCJEC_jtID0v2_eta50 -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
sleep 30s

#source condorSubmit_readForests.sh readForests_ppMC_MCJEC_jtID0 -1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_MCJEC_jtIDL -1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_MCJEC_jtIDT -1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0

