#!/bin/bash

#make_rooUnfold

rooUnfoldCompile bayesUnfoldDataSpectra_wNLO_etabin.C


############# NNPDF NLO + HERWIG EE5C LO NPCs ###############
####NNPDF_03.14.19_sigmu_noJERscales_spl3wgts_gaussSmear_NLO_plots
#
#FITTYPE="sigmu"
#DATADIR="08.06.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_08-06-19_jetPlots_0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal"
#THYDIR="smearTheory/NNPDF_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_NLO_plots/"
#THYFILESTR="NNPDF_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_"
#OUTPUTTAG="08.06.19_Bayes_NNPDF_NLO_${FITTYPE}_noJERscales_withLumiCorr_semifinal"

FITTYPE="sigmu"
DATADIR="08.06.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_08-06-19_jetPlots_0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal"
THYDIR="smearTheory/NNPDF_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_NLO_plots/"
THYFILESTR="NNPDF_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_"
OUTPUTTAG="08.06.19_Bayes_NNPDF_NLO_${FITTYPE}_noJERscales_withLumiCorr_SMPbins_ptLo64_semifinal"

#FITTYPE="sigmu"
#DATADIR="08.06.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_08-06-19_jetPlots_0.0eta2.0_JOHNbins_withJECsys_withLumiCorr_semifinal"
#THYDIR="smearTheory/NNPDF_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_NLO_plots/"
#THYFILESTR="NNPDF_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_"
#OUTPUTTAG="08.06.19_Bayes_NNPDF_NLO_${FITTYPE}_noJERscales_withLumiCorr_JOHNbins_semifinal"

###vanilla no sys
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_nosys"  0    ${THYDIR} ${THYFILESTR}  0 "" ""
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_nosys"  1    ${THYDIR} ${THYFILESTR}  0 "" ""
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_nosys"  2    ${THYDIR} ${THYFILESTR}  0 "" ""
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_nosys"  3    ${THYDIR} ${THYFILESTR}  0 "" ""
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_nosys_NLOMC_wNP_??eta??.*"
#return
#
##JEC sys
./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JECsys"  0    ${THYDIR} ${THYFILESTR}  1 "JEC" ""
./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JECsys"  1    ${THYDIR} ${THYFILESTR}  1 "JEC" ""
./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JECsys"  2    ${THYDIR} ${THYFILESTR}  1 "JEC" ""
./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JECsys"  3    ${THYDIR} ${THYFILESTR}  1 "JEC" ""
scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_JECsys_NLOMC_wNP_??eta??.*"
return
##
#####JER sys
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JERsys"  0    ${THYDIR} ${THYFILESTR}  1 "JER" ""
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JERsys"  1    ${THYDIR} ${THYFILESTR}  1 "JER" ""
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JERsys"  2    ${THYDIR} ${THYFILESTR}  1 "JER" ""
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JERsys"  3    ${THYDIR} ${THYFILESTR}  1 "JER" ""
###scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_JERsys_NLOMC_wNP_??eta??.*"
###return
##
####NP sys12, alt NP choices, HERWIG EE4C/PYTHIA8
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsys12" 0   ${THYDIR} ${THYFILESTR}  1 "NP" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsys12" 1   ${THYDIR} ${THYFILESTR}  1 "NP" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsys12" 2   ${THYDIR} ${THYFILESTR}  1 "NP" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsys12" 3   ${THYDIR} ${THYFILESTR}  1 "NP" "12"
###scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_NPsys12_NLOMC_wNP_??eta??.*"
###return
##
####NP sysupdown, shifted HERWIG EE5C NP up/down, 1 sigma
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsysupdown" 0 ${THYDIR} ${THYFILESTR}  1 "NP" "updown"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsysupdown" 1 ${THYDIR} ${THYFILESTR}  1 "NP" "updown"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsysupdown" 2 ${THYDIR} ${THYFILESTR}  1 "NP" "updown"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsysupdown" 3 ${THYDIR} ${THYFILESTR}  1 "NP" "updown"
###scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_NPsysupdown_NLOMC_wNP_??eta??.*"
###return
##
####PDF sys12, alt PDF choices, CT14/HERA
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsys12" 0   ${THYDIR} ${THYFILESTR}  1 "PDF" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsys12" 1   ${THYDIR} ${THYFILESTR}  1 "PDF" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsys12" 2   ${THYDIR} ${THYFILESTR}  1 "PDF" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsys12" 3   ${THYDIR} ${THYFILESTR}  1 "PDF" "12"
###scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_PDFsys12_NLOMC_wNP_??eta??.*"
###return
##
####PDF sysupdown, shifted NNPDF PDF up/down, 6 Pt scale + PDF errs, 1 sigma
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsysupdown" 0 ${THYDIR} ${THYFILESTR}  1 "PDF" "updown"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsysupdown" 1 ${THYDIR} ${THYFILESTR}  1 "PDF" "updown"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsysupdown" 2 ${THYDIR} ${THYFILESTR}  1 "PDF" "updown"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsysupdown" 3 ${THYDIR} ${THYFILESTR}  1 "PDF" "updown"
##scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_PDFsysupdown_NLOMC_wNP_??eta??.*"
##return

#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_*sys*_NLOMC_wNP_??eta??.*"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_*sys*_NLOMC_wNP_??eta??_respMat.*"
#return
###-----







############## SMEARED GEN PY 8 LO + HERWIG EE5C LO NPCs ###############
####PY8_03.14.19_sigmu_noJERscales_spl3wgts_gaussSmear_NLO_plots
#
FITTYPE="sigmu"
DATADIR="06.25.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_06-25-19_jetPlots_0.0eta2.0_SMPbins_withJECsys_semifinal"
THYDIR="smearTheory/PY8_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_NLO_plots/"
THYFILESTR="PY8_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_"
OUTPUTTAG="06.25.19_Bayes_PY8_LO_${FITTYPE}_noJERscales_semifinal"
#
#
###vanilla no sys
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_nosys"  0    ${THYDIR} ${THYFILESTR}  0 "" ""
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_nosys"  1    ${THYDIR} ${THYFILESTR}  0 "" ""
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_nosys"  2    ${THYDIR} ${THYFILESTR}  0 "" ""
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_nosys"  3    ${THYDIR} ${THYFILESTR}  0 "" ""
##scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_nosys_NLOMC_wNP_??eta??.*"
##return
#
###JEC sys
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JECsys"  0    ${THYDIR} ${THYFILESTR}  1 "JEC" ""
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JECsys"  1    ${THYDIR} ${THYFILESTR}  1 "JEC" ""
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JECsys"  2    ${THYDIR} ${THYFILESTR}  1 "JEC" ""
./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JECsys"  3    ${THYDIR} ${THYFILESTR}  1 "JEC" ""
##scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_JECsys_NLOMC_wNP_??eta??.*"
##return
#
###JER sys
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JERsys"  0    ${THYDIR} ${THYFILESTR}  1 "JER" ""
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JERsys"  1    ${THYDIR} ${THYFILESTR}  1 "JER" ""
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JERsys"  2    ${THYDIR} ${THYFILESTR}  1 "JER" ""
./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JERsys"  3    ${THYDIR} ${THYFILESTR}  1 "JER" ""
##scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_JERsys_NLOMC_wNP_??eta??.*"
##return
#
###NP sys12, alt NP choices, HERWIG EE4C/PYTHIA8
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsys12" 0   ${THYDIR} ${THYFILESTR}  1 "NP" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsys12" 1   ${THYDIR} ${THYFILESTR}  1 "NP" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsys12" 2   ${THYDIR} ${THYFILESTR}  1 "NP" "12"
./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsys12" 3   ${THYDIR} ${THYFILESTR}  1 "NP" "12"
##scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_NPsys12_NLOMC_wNP_??eta??.*"
##return
#
###NP sysupdown, shifted HERWIG EE5C NP up/down, 1 sigma
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsysupdown" 0 ${THYDIR} ${THYFILESTR}  1 "NP" "updown"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsysupdown" 1 ${THYDIR} ${THYFILESTR}  1 "NP" "updown"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsysupdown" 2 ${THYDIR} ${THYFILESTR}  1 "NP" "updown"
./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsysupdown" 3 ${THYDIR} ${THYFILESTR}  1 "NP" "updown"
##scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_NPsysupdown_NLOMC_wNP_??eta??.*"
##return
#
###PDF sys12, alt PDF choices, CT14/HERA
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsys12" 0   ${THYDIR} ${THYFILESTR}  1 "PDF" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsys12" 1   ${THYDIR} ${THYFILESTR}  1 "PDF" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsys12" 2   ${THYDIR} ${THYFILESTR}  1 "PDF" "12"
./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsys12" 3   ${THYDIR} ${THYFILESTR}  1 "PDF" "12"
##scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_PDFsys12_NLOMC_wNP_??eta??.*"
##return
#
scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_*sys*_NLOMC_wNP_??eta??.*"
return
###-----
