#!/bin/bash

#make_rooUnfold

rooUnfoldCompile bayesUnfoldDataSpectra_etabin.C

############ FULL GEN + RECO PY8 ###############

###NO JEC SYS
#DATADIR="06.07.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_06-07-19_jetPlots_0.0eta2.0_SMPbins_withJECsys"
#MCDIR="05.17.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-17-19_JERS_0.0eta2.0_unf/"
#MCFILESTR="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
#OUTPUTTAG="Bayes_PY8_FullRECO_06.07.19_SMPbins"
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  0  ${MCDIR} ${MCFILESTR}
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  1  ${MCDIR} ${MCFILESTR}
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  2  ${MCDIR} ${MCFILESTR}
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  3  ${MCDIR} ${MCFILESTR}
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_LOMC_??eta??.*"

###WITH JEC SYS
#DATADIR="06.07.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_06-07-19_jetPlots_0.0eta2.0_SMPbins_withJECsys"
#MCDIR="05.17.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-17-19_JERS_0.0eta2.0_unf/"
#MCFILESTR="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
#OUTPUTTAG="Bayes_PY8_FullRECO_06.07.19_SMPbins_withJECsys"
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  0  ${MCDIR} ${MCFILESTR}
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  1  ${MCDIR} ${MCFILESTR}
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  2  ${MCDIR} ${MCFILESTR}
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  3  ${MCDIR} ${MCFILESTR}
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_LOMC_??eta??.*"

## NO JEC SYS + CHRIS MCGINNS BINNING
#DATADIR="06.07.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_06-07-19_jetPlots_0.0eta2.0_Chrisbins"
#MCDIR="05.17.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-17-19_JERS_0.0eta2.0_unf/"
#MCFILESTR="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
#OUTPUTTAG="Bayes_PY8_FullRECO_06.07.19_Chrisbins"
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  0  ${MCDIR} ${MCFILESTR}
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  1  ${MCDIR} ${MCFILESTR}
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  2  ${MCDIR} ${MCFILESTR}
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  3  ${MCDIR} ${MCFILESTR}
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_LOMC_??eta??.*"







#------------------------------#

