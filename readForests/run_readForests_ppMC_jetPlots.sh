#!/bin/bash

source condorSubmit_readForests.sh readForests_ppMC_jetPlots_jtIDTv2_eta27to50 -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
sleep 30s

source condorSubmit_readForests.sh readForests_ppMC_jetPlots_jtIDTv2_eta30 -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
sleep 30s

source condorSubmit_readForests.sh readForests_ppMC_jetPlots_jtIDTv2_eta20 -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
sleep 30s