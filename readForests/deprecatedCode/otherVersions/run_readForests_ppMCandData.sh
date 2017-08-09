#!/bin/bash

##### COMPILE READFORESTS
#echo "compiling readforests code(s)..."
#echo ""
#echo ""
#rootcompile readForests_ppData.C 
#rootcompile readForests_ppMC.C 

##### TEST READFORESTS
#echo "running readForests tests..."
#echo ""
#echo ""
#./readForests_ppData.exe
#./readForests_ppMC.exe
#./readForests_ppMC_JetQA_RapBins.exe  filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 0 1 3 PF 0 readForests_ppMC_JetQA_RapBins_defOut.root  
#./readForests_ppMC_Unfold.exe         filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 0 1 3 PF 0 readForests_ppMC_Unfold_defOut.root  
#./readForests_ppMC_MCEff.exe          filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 0 1 3 PF 0 readForests_ppMC_MCEff_defOut.root  
#./readForests_ppMC_MCJEC_JetID0.exe   filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 0 1 3 PF 0 readForests_ppMC_MCJEC_JetID0_defOut.root  
#./readForests_ppMC_MCJEC_JetID1.exe   filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 0 1 3 PF 0 readForests_ppMC_MCJEC_JetID1_defOut.root  



#### TEST CONDOR SUBMIT
#source condorSubmit_readForests.sh readForests_ppData 1 13 0 filelists/5p02TeV_HighPtJet80_forests.txt      3 PF 0
#source condorSubmit_readForests.sh readForests_ppData 1 13 0 filelists/5p02TeV_HighPtLowerJets_forests.txt      3 PF 0
#
#source condorSubmit_readForests.sh readForests_ppMC_JetQA_RapBins 1 5 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_Unfold        1 5 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_MCEff         1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_MCJEC_JetID0  1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_MCJEC_JetID1  1 40 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0



#### CONDOR SUBMIT
#### read ak3PFjetAnalyzer
#source condorSubmit_readForests.sh readForests_ppData -1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt      3 PF 0
#source condorSubmit_readForests.sh readForests_ppData -1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt  3 PF 0
#
#source condorSubmit_readForests.sh readForests_ppMC_JetQA_RapBins -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
#source condorSubmit_readForests.sh readForests_ppMC__JetQA_MCEff -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_Unfold        -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_MCEff         -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
source condorSubmit_readForests.sh readForests_ppMC_MCJEC_JetID0  -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_MCJEC_JetID1  -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 3 PF 0


#### read ak4PFjetAnalyzer
#source condorSubmit_readForests.sh readForests_ppData -1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt      4 PF 0
#source condorSubmit_readForests.sh readForests_ppData -1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt  4 PF 0
#
#source condorSubmit_readForests.sh readForests_ppMC_JetQA_RapBins -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 
#source condorSubmit_readForests.sh readForests_ppMC__JetQA_MCEff -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_Unfold        -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_MCEff         -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
source condorSubmit_readForests.sh readForests_ppMC_MCJEC_JetID0  -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_MCJEC_JetID1  -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0


#### read ak5PFjetAnalyzer
#source condorSubmit_readForests.sh readForests_ppData -1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt      5 PF 0
#source condorSubmit_readForests.sh readForests_ppData -1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt  5 PF 0

#source condorSubmit_readForests.sh readForests_ppMC_JetQA_RapBins -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 5 PF 0
#source condorSubmit_readForests.sh readForests_ppMC__JetQA_MCEff -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 5 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_Unfold        -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 5 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_MCEff         -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 5 PF 0
source condorSubmit_readForests.sh readForests_ppMC_MCJEC_JetID0  -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 5 PF 0
#source condorSubmit_readForests.sh readForests_ppMC_MCJEC_JetID1  -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 5 PF 0

