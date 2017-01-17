#!/bin/bash



##### ppData, eta20, jetQAPlots 
##### R=3
#source condorSubmit_readForests.sh readForests_ppData_jetPlots_jtIDv3L_eta20  -1 6 0 filelists/5p02TeV_HighPtJet80_forests.txt 3 PF 0
#sleep 10s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots_jtIDv3L_eta20  -1 6 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 3 PF 0
#sleep 10s
#
##### R=4
#source condorSubmit_readForests.sh readForests_ppData_jetPlots_jtIDv3L_eta20  -1 6 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0
#sleep 10s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots_jtIDv3L_eta20  -1 6 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0
#sleep 10s
#
#
##### ppData, eta47, jetPlots + jetIDEff Plots
##### R=3
#source condorSubmit_readForests.sh readForests_ppData_jetPlots_jtIDv3L_eta47  -1 6 0 filelists/5p02TeV_HighPtJet80_forests.txt 3 PF 0
#sleep 10s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots_jtIDv3L_eta47  -1 6 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 3 PF 0
#sleep 10s
#
##### R=4
#source condorSubmit_readForests.sh readForests_ppData_jetPlots_jtIDv3L_eta47  -1 6 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0
#sleep 10s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots_jtIDv3L_eta47  -1 6 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0
#sleep 10s




#### ppMC, eta20, jetQA, jetIDEff
### R=3
source condorSubmit_readForests.sh readForests_ppMC_jetPlots_jtIDv3L_eta20 -1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
sleep 30s

### R=4
source condorSubmit_readForests.sh readForests_ppMC_jetPlots_jtIDv3L_eta20 -1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
sleep 30s

#### ppMC, eta47, jetQA, jetIDEff Plots
### R=3
source condorSubmit_readForests.sh readForests_ppMC_jetPlots_jtIDv3L_eta47 -1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
sleep 30s

### R=4
source condorSubmit_readForests.sh readForests_ppMC_jetPlots_jtIDv3L_eta47 -1 25 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
sleep 30s