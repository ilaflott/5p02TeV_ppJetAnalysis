#!/bin/bash

##### COMPILE READFORESTS
#echo "compiling readforests code(s)..."
#echo ""
#echo ""
#rootcompile readForests_ppData.C 


##### TEST READFORESTS
#echo "running readForests tests..."
#echo ""
#echo ""
#./readForests_ppData.exe


#### TEST CONDOR SUBMIT
#source condorSubmit_readForests.sh readForests_ppData 1 13 0 filelists/5p02TeV_HighPtJet80_forests.txt      3 PF 0
#source condorSubmit_readForests.sh readForests_ppData 1 13 0 filelists/5p02TeV_HighPtLowerJets_forests.txt      3 PF 0


#### CONDOR SUBMIT
#### read ak3PFjetAnalyzer
source condorSubmit_readForests.sh readForests_ppData -1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt      3 PF 0
source condorSubmit_readForests.sh readForests_ppData -1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt  3 PF 0



#### read ak4PFjetAnalyzer
source condorSubmit_readForests.sh readForests_ppData -1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt      4 PF 0
source condorSubmit_readForests.sh readForests_ppData -1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt  4 PF 0


#### read ak5PFjetAnalyzer
source condorSubmit_readForests.sh readForests_ppData -1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt      5 PF 0
source condorSubmit_readForests.sh readForests_ppData -1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt  5 PF 0


