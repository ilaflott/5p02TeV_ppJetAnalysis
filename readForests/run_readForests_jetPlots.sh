#!/bin/bash



#### eta20
#rootcompile_ver readForests_ppMC_jetPlots.C jtIDv3L_eta20
#
#rootcompile_ver readForests_ppData_jetPlots.C jtIDv3L_eta20
#
#
#### R=3
#source condorSubmit_readForests.sh readForests_ppMC_jetPlots_jtIDv3L_eta20 1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
#sleep 30s
#
source condorSubmit_readForests.sh readForests_ppData_jetPlots_jtIDv3L_eta20 1 6 0 filelists/5p02TeV_HighPtJet80_forests.txt 3 PF 0
sleep 10s
source condorSubmit_readForests.sh readForests_ppData_jetPlots_jtIDv3L_eta20 1 6 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 3 PF 0
sleep 10s
#
#### R=4
#source condorSubmit_readForests.sh readForests_ppMC_jetPlots_jtIDv3L_eta20 1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
#sleep 30s
#
source condorSubmit_readForests.sh readForests_ppData_jetPlots_jtIDv3L_eta20 1 6 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0
sleep 30s
source condorSubmit_readForests.sh readForests_ppData_jetPlots_jtIDv3L_eta20 1 6 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0
sleep 30s




#### eta50
#rootcompile_ver readForests_ppMC_jetPlots.C jtIDv3L_eta50
#
#rootcompile_ver readForests_ppData_jetPlots.C jtIDv3L_eta50
#
#
#### R=3
#source condorSubmit_readForests.sh readForests_ppMC_jetPlots_jtIDv3L_eta50 1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
#sleep 30s
#
#source condorSubmit_readForests.sh readForests_ppData_jetPlots_jtIDv3L_eta50 1 6 0 filelists/5p02TeV_HighPtJet80_forests.txt 3 PF 0
#sleep 30s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots_jtIDv3L_eta50 1 6 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 3 PF 0
#sleep 30s
#
#### R=4
#source condorSubmit_readForests.sh readForests_ppMC_jetPlots_jtIDv3L_eta50 1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
#sleep 30s
#
#source condorSubmit_readForests.sh readForests_ppData_jetPlots_jtIDv3L_eta50 1 6 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0
#sleep 30s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots_jtIDv3L_eta50 1 6 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0
#sleep 30s
