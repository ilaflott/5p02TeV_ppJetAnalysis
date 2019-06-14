#!/bin/bash

#make_rooUnfold

rooUnfoldCompile bayesUnfoldDataSpectra_etabin.C

############ FULL GEN + RECO PY8 ###############

#./bayesUnfoldDataSpectra_etabin.exe "05.17.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_05-17-19_jetPlots_0.0eta2.0_SMPbins" "Bayes_PY8_FullRECO_05.17.19_SMPbins"  0  "05.17.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-17-19_JERS_0.0eta2.0_unf/" "Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root" 
#
#./bayesUnfoldDataSpectra_etabin.exe "05.17.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_05-17-19_jetPlots_0.0eta2.0_SMPbins" "Bayes_PY8_FullRECO_05.17.19_SMPbins"  1  "05.17.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-17-19_JERS_0.0eta2.0_unf/" "Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root" 
#
#./bayesUnfoldDataSpectra_etabin.exe "05.17.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_05-17-19_jetPlots_0.0eta2.0_SMPbins" "Bayes_PY8_FullRECO_05.17.19_SMPbins"  2  "05.17.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-17-19_JERS_0.0eta2.0_unf/" "Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root" 
#
#./bayesUnfoldDataSpectra_etabin.exe "05.17.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_05-17-19_jetPlots_0.0eta2.0_SMPbins" "Bayes_PY8_FullRECO_05.17.19_SMPbins"  3  "05.17.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-17-19_JERS_0.0eta2.0_unf/" "Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root" 
#
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_PY8_FullRECO_05.17.19_SMPbins_LOMC_??eta??.*"


./bayesUnfoldDataSpectra_etabin.exe "06.07.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_06-07-19_jetPlots_0.0eta2.0_SMPbins_withJECsys" "Bayes_PY8_FullRECO_06.07.19_SMPbins_withJECsys"  0  "05.17.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-17-19_JERS_0.0eta2.0_unf/" "Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root" 
./bayesUnfoldDataSpectra_etabin.exe "06.07.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_06-07-19_jetPlots_0.0eta2.0_SMPbins_withJECsys" "Bayes_PY8_FullRECO_06.07.19_SMPbins_withJECsys"  1  "05.17.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-17-19_JERS_0.0eta2.0_unf/" "Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root" 
./bayesUnfoldDataSpectra_etabin.exe "06.07.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_06-07-19_jetPlots_0.0eta2.0_SMPbins_withJECsys" "Bayes_PY8_FullRECO_06.07.19_SMPbins_withJECsys"  2  "05.17.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-17-19_JERS_0.0eta2.0_unf/" "Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root" 
./bayesUnfoldDataSpectra_etabin.exe "06.07.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_06-07-19_jetPlots_0.0eta2.0_SMPbins_withJECsys" "Bayes_PY8_FullRECO_06.07.19_SMPbins_withJECsys"  3  "05.17.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-17-19_JERS_0.0eta2.0_unf/" "Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root" 

scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_PY8_FullRECO_06.07.19_SMPbins_withJECsys_LOMC_??eta??.*"




#./bayesUnfoldDataSpectra_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_PY8_FullRECO"  0  "02.18.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-18-19_JERS_0.0eta3.0_unf/" "Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root" 
##scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_PY8_FullRECO_LOMC_00eta05.*"
#./bayesUnfoldDataSpectra_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_PY8_FullRECO"  1  "02.18.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-18-19_JERS_0.0eta3.0_unf/" "Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root" 
##scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_PY8_FullRECO_LOMC_00eta05.*"
#./bayesUnfoldDataSpectra_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_PY8_FullRECO"  2  "02.18.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-18-19_JERS_0.0eta3.0_unf/" "Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root" 
##scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_PY8_FullRECO_LOMC_00eta05.*"
#./bayesUnfoldDataSpectra_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_PY8_FullRECO"  3  "02.18.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-18-19_JERS_0.0eta3.0_unf/" "Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root" 
##scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_PY8_FullRECO_LOMC_00eta05.*"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_PY8_FullRECO_LOMC_??eta??.*"


#------------------------------#


#rooUnfoldCompile SVDUnfoldDataSpectra_wNLO_etabin.C

#./SVDUnfoldDataSpectra_wNLO_etabin.exe "test"

#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_*_*_NLOMC_wNP_??eta??*.*"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_HLT60Thresh*_NLOMC_wNP_00eta05*"