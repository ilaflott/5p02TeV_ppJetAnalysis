#!/bin/bash

#make_rooUnfold

rooUnfoldCompile bayesUnfoldMCSpectra_wNLO_etabin.C


############ NNPDF NLO + HERWIG EE5C LO NPCs ###############

FITTYPE="sigmu"
THYDIR="smearTheory/CT14NNLO_NLO_v3_03.18.20_${FITTYPE}_geny_spl3wgts_gaussSmear_murmufpt_plots/"
THYFILESTR="CT14NNLO_NLO_v3_03.18.20_${FITTYPE}_geny_spl3wgts_gaussSmear_murmufpt_"
OUTPUTTAG="03.18.20_Bayes_Closure_CT14NNLO_NLO_v3_${FITTYPE}_SMPbins_withjety_semifinalv3_REBIN_NODATACOVMAT_murmufpt"
etabinStart=0
etabinEnd=4
for (( j=${etabinStart}; j<${etabinEnd}; j++ ));
do
    ./bayesUnfoldMCSpectra_wNLO_etabin.exe ${OUTPUTTAG}  $j  ${THYDIR} ${THYFILESTR} 
done
return
#scp2Serin "output/unfoldMCSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_NLOMC_wNP_??eta??.*"
#return



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
