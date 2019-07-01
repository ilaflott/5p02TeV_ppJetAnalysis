#!/bin/bash

#make_rooUnfold

rooUnfoldCompile bayesUnfoldMCSpectra_wNLO_etabin.C


############ NNPDF NLO + HERWIG EE5C LO NPCs ###############

FITTYPE="sigmu"
THYDIR="smearTheory/NNPDF_03.14.19_${FITTYPE}_noJERscales_spl3wgts_gaussSmear_NLO_plots/"
THYFILESTR="NNPDF_03.14.19_${FITTYPE}_noJERscales_spl3wgts_gaussSmear_"
OUTPUTTAG="06.07.19_Bayes_Closure_NNPDF_NLO_${FITTYPE}_noJERscales"


##vanilla no sys
./bayesUnfoldMCSpectra_wNLO_etabin.exe  "${OUTPUTTAG}"  0    ${THYDIR} ${THYFILESTR} 
#./bayesUnfoldMCSpectra_wNLO_etabin.exe  "${OUTPUTTAG}"  1    ${THYDIR} ${THYFILESTR} 
#./bayesUnfoldMCSpectra_wNLO_etabin.exe  "${OUTPUTTAG}"  2    ${THYDIR} ${THYFILESTR} 
#./bayesUnfoldMCSpectra_wNLO_etabin.exe  "${OUTPUTTAG}"  3    ${THYDIR} ${THYFILESTR} 
scp2Serin "output/unfoldMCSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_NLOMC_wNP_??eta??.*"
return



############## SMEARED GEN PY 8 LO + HERWIG EE5C LO NPCs ###############
#
#FITTYPE="sigmu"
#THYDIR="smearTheory/PY8_03.14.19_${FITTYPE}_spl3wgts_gaussSmear_NLO_plots/"
#THYFILESTR="PY8_03.14.19_${FITTYPE}_spl3wgts_gaussSmear_"
#OUTPUTTAG="06.07.19_Bayes_Closure_PY8_LO_${FITTYPE}"
#
#
###vanilla no sys
##./bayesUnfoldMCSpectra_wNLO_etabin.exe  "${OUTPUTTAG}"  0   
##./bayesUnfoldMCSpectra_wNLO_etabin.exe  "${OUTPUTTAG}"  1    
##./bayesUnfoldMCSpectra_wNLO_etabin.exe  "${OUTPUTTAG}"  2    
##./bayesUnfoldMCSpectra_wNLO_etabin.exe  "${OUTPUTTAG}"  3    
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_nosys_NLOMC_wNP_??eta??.*"
##return
