#!/bin/bash

rootcompile readForests_ppData_jetPlots.C

#### TEST ppData jetPlots
#source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 2 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0  0.0 1.0
#source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 2 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 0.0 1.0
#source askCondor.sh

### Jet80, then LowerJets
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 8 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0  0.0 1.0
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 8 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 0.0 1.0
#sleep 1s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 8 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0  1.0 2.0
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 8 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 1.0 2.0
sleep 1s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 8 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0  2.0 3.0
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 8 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 2.0 3.0
#sleep 1s
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 8 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0  3.2 4.7
#source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 8 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 3.2 4.7
#sleep 1s
#source askCondor.sh

return