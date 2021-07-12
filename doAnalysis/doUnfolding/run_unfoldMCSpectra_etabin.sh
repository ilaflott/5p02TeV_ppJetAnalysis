#!/bin/bash

#make_rooUnfold

rooUnfoldCompile bayesUnfoldMCSpectra_etabin.C

############ FULL GEN + RECO PY8 ###############
#MCDIR="03.18.20_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_03-18-20_unf_0.0eta2.5_SMPbins_withjety_semifinalv3/"
#MCDIR="03.18.20_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_03-18-20_unf_0.0eta2.5_SMPbins_withgenjety_semifinalv3/"
#MCFILESTR="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
#OUTPUTTAG="03.18.20_Bayes_Closure_PY8_FullRECO_withgenjety_SMPbins_semifinalv3"
MCDIR="01.05.21_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_01-05-21_unf_0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4/"
MCFILESTR="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
OUTPUTTAG="01.05.21_Bayes_Closure_PY8_FullRECO_withgenjety_SMPbins_semifinalv4"
etabinStart=0
etabinEnd=1
for (( j=${etabinStart}; j<${etabinEnd}; j++ ));
do
    ./bayesUnfoldMCSpectra_etabin.exe ${OUTPUTTAG}  $j  ${MCDIR} ${MCFILESTR} "y"
done
return
#./bayesUnfoldMCSpectra_etabin.exe  ${OUTPUTTAG}  3  ${MCDIR} ${MCFILESTR}
#scp2Serin "output/unfoldMCSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_LOMC_??eta??.*"

#------------------------------#

