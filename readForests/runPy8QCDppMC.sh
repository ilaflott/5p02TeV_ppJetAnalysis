#!/bin/bash

## Usage is...
## source condorSubmit_ppMC.sh <NJobs> <NFilesPerJob> <startFilePos> <filelistIn> <radius> <jetType> <debug>

## read ak4CalojetAnalyzer
source condorSubmit_readForests_ppMC.sh -1 15 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 Calo 0

## read ak3PFjetAnalyzer
source condorSubmit_readForests_ppMC.sh -1 15 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0

# read ak4PFjetAnalyzer
source condorSubmit_readForests_ppMC.sh -1 15 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
