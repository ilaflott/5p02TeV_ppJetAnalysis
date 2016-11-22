#!/bin/bash

##### CONDOR SUBMIT ppData TEST
##### read ak3PFjetAnalyzer
#source condorSubmit_readForests.sh readForests_ppData 1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt  3 PF 1
#
##### read ak4PFjetAnalyzer
#source condorSubmit_readForests.sh readForests_ppData 1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt  4 PF 1
#
###### read ak5PFjetAnalyzer
##source condorSubmit_readForests.sh readForests_ppData 1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt  5 PF 0
#
#
#
#### CONDOR SUBMIT ppMC TEST
##### read ak3PFjetAnalyzer
#source condorSubmit_readForests.sh readForests_ppMC 1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt      3 PF 1
#
##### read ak4PFjetAnalyzer
#source condorSubmit_readForests.sh readForests_ppMC 1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt      4 PF 1
#
####### read ak5PFjetAnalyzer
###source condorSubmit_readForests.sh readForests_ppMC 1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt      5 PF 0



#### CONDOR SUBMIT ppData
#### read ak3PFjetAnalyzer
source condorSubmit_readForests.sh readForests_ppData -1 8 0 filelists/5p02TeV_HighPtJet80_forests.txt      3 PF 0
source condorSubmit_readForests.sh readForests_ppData -1 8 0 filelists/5p02TeV_HighPtLowerJets_forests.txt  3 PF 0
echo "done submitting ak3PF ppData...."
echo ""

#### read ak4PFjetAnalyzer
source condorSubmit_readForests.sh readForests_ppData -1 8 0 filelists/5p02TeV_HighPtJet80_forests.txt      4 PF 0
source condorSubmit_readForests.sh readForests_ppData -1 8 0 filelists/5p02TeV_HighPtLowerJets_forests.txt  4 PF 0
echo "done submitting ak4PF ppData...."
echo ""

##### read ak5PFjetAnalyzer
#source condorSubmit_readForests.sh readForests_ppData -1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt      5 PF 0
#source condorSubmit_readForests.sh readForests_ppData -1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt  5 PF 0



#### CONDOR SUBMIT ppData
#### read ak3PFjetAnalyzer
source condorSubmit_readForests.sh readForests_ppMC -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt      3 PF 0
echo "done submitting ak3PF ppMC...."
echo ""

#### read ak4PFjetAnalyzer
source condorSubmit_readForests.sh readForests_ppMC -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt      4 PF 0
echo "done submitting ak4PF ppMC...."
echo ""

##### read ak5PFjetAnalyzer
#source condorSubmit_readForests.sh readForests_ppMC -1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt      5 PF 0


