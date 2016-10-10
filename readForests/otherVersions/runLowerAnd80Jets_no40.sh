#!/bin/bash

## Usage is...
## source condorSubmit_ppData.sh <NJobs> <NFilesPerJob> <startFilePos> <filelistIn> <radius> <jetType> <debug>

### read ak4CalojetAnalyzer
#source condorSubmit_readForests_ppData.sh -1 8 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 Calo 0
#source condorSubmit_readForests_ppData.sh -1 10 8 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 Calo 0
#
### read ak3PFjetAnalyzer
#source condorSubmit_readForests_ppData.sh -1 8 0 filelists/5p02TeV_HighPtJet80_forests.txt 3 PF 0
#source condorSubmit_readForests_ppData.sh -1 10 8 filelists/5p02TeV_HighPtLowerJets_forests.txt 3 PF 0

# read ak4PFjetAnalyzer
source condorSubmit_readForests_ppData.sh -1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0
source condorSubmit_readForests_ppData.sh -1 8 8 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0