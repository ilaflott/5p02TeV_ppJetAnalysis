#!/bin/bash

#make_rooUnfold

rooUnfoldCompile bayesUnfoldMCSpectra_etabin.C

############ FULL GEN + RECO PY8 ###############

MCDIR="06.25.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_06-25-19_unf_0.0eta2.0_semifinal/"
MCFILESTR="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
OUTPUTTAG="Bayes_Closure_PY8_FullRECO_06.25.19_SMPbins"
./bayesUnfoldMCSpectra_etabin.exe  ${OUTPUTTAG}  0  ${MCDIR} ${MCFILESTR}
./bayesUnfoldMCSpectra_etabin.exe  ${OUTPUTTAG}  1  ${MCDIR} ${MCFILESTR}
./bayesUnfoldMCSpectra_etabin.exe  ${OUTPUTTAG}  2  ${MCDIR} ${MCFILESTR}
./bayesUnfoldMCSpectra_etabin.exe  ${OUTPUTTAG}  3  ${MCDIR} ${MCFILESTR}
scp2Serin "output/unfoldMCSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_LOMC_??eta??.*"

#------------------------------#

