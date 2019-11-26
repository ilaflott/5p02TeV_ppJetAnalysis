#!/bin/bash

#make_rooUnfold

rooUnfoldCompile bayesUnfoldDataSpectra_etabin.C

############ FULL GEN + RECO PY8 ###############


###WITH JEC SYS, WITH LUMI CORR, SMP BINS PT LO 56
#DATADIR="08.06.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_08-06-19_jetPlots_0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal"
#MCDIR="06.25.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_06-25-19_unf_0.0eta2.0_semifinal/"
#MCFILESTR="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
#OUTPUTTAG="Bayes_PY8_FullRECO_08.06.19_SMPbins_withJECsys_withLumiCorr"

###WITH JEC SYS, WITH LUMI CORR, SMP BINS PT LO 64
#DATADIR="08.06.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_08-06-19_jetPlots_0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal"
#MCDIR="06.25.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_06-25-19_unf_0.0eta2.0_semifinal/"
#MCFILESTR="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
#OUTPUTTAG="Bayes_PY8_FullRECO_08.06.19_SMPbins_withJECsys_SMPbins_ptLo64_withLumiCorr"

###WITH JEC SYS V2, WITH LUMI CORR, SMP BINS PT LO 56
DATADIR="08.06.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_08-06-19_jetPlots_0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinal"
MCDIR="06.25.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_06-25-19_unf_0.0eta2.0_semifinal/"
MCFILESTR="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
OUTPUTTAG="Bayes_PY8_FullRECO_08.06.19_SMPbins_withJECsysv2_SMPbins_ptLo56_withLumiCorr"

###WITH JEC SYS V2, WITH LUMI CORR, JOHN BINS
#DATADIR="08.06.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_08-06-19_jetPlots_0.0eta2.0_JOHNbins_withJECsys_withLumiCorr_semifinal"
#MCDIR="06.25.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_06-25-19_unf_0.0eta2.0_semifinal/"
#MCFILESTR="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
#OUTPUTTAG="Bayes_PY8_FullRECO_08.06.19_JOHNbins_withJECsys_withLumiCorr"

###WITH JEC SYS V2, WITH LUMI CORR, CHRIS BINS
#DATADIR="08.06.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_08-06-19_jetPlots_0.0eta2.0_CHRISbins2_JECsysv2_withLumiCorr_semifinal"
#MCDIR="06.25.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_06-25-19_unf_0.0eta2.0_semifinal/"
#MCFILESTR="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
#OUTPUTTAG="Bayes_PY8_FullRECO_08.06.19_CHRISbins2_withJECsysv2_SMPbins_withLumiCorr"

./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  0  ${MCDIR} ${MCFILESTR}
./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  1  ${MCDIR} ${MCFILESTR}
./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  2  ${MCDIR} ${MCFILESTR}
./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  3  ${MCDIR} ${MCFILESTR}

scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_LOMC_??eta??.*"
scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_LOMC_??eta??_respMat.*"
source run_unfoldDataSpectra_wNLO_etabin.sh





#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  0  ${MCDIR} ${MCFILESTR}
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  1  ${MCDIR} ${MCFILESTR}
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  2  ${MCDIR} ${MCFILESTR}
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  3  ${MCDIR} ${MCFILESTR}
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_LOMC_??eta??.*"







#------------------------------#

