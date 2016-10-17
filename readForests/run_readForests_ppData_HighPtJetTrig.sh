#!/bin/bash

## read ak4CalojetAnalyzer
#source condorSubmit_readForests_ppData.sh -1 8 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 Calo 0
#source condorSubmit_readForests_ppData.sh -1 10 8 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 Calo 0



### read ak3PFjetAnalyzer
source condorSubmit_readForests_ppData.sh readForests_ppData -1 5 0 filelists/5p02TeV_HighPtJet80_forests.txt      3 PF 0
source condorSubmit_readForests_ppData.sh readForests_ppData -1 5 0 filelists/5p02TeV_HighPtLowerJets_forests.txt  3 PF 0
sleep 30s

### read ak4PFjetAnalyzer
source condorSubmit_readForests_ppData.sh readForests_ppData -1 5 0 filelists/5p02TeV_HighPtJet80_forests.txt      4 PF 0
source condorSubmit_readForests_ppData.sh readForests_ppData -1 5 0 filelists/5p02TeV_HighPtLowerJets_forests.txt  4 PF 0
sleep 30s

### read ak5PFjetAnalyzer
source condorSubmit_readForests_ppData.sh readForests_ppData -1 5 0 filelists/5p02TeV_HighPtJet80_forests.txt      5 PF 0
source condorSubmit_readForests_ppData.sh readForests_ppData -1 5 0 filelists/5p02TeV_HighPtLowerJets_forests.txt  5 PF 0
sleep 30s