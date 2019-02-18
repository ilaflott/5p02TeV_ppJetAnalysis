#!/bin/bash

#make_rooUnfold

rooUnfoldCompile bayesUnfoldDataSpectra_wNLO_etabin.C

#./bayesUnfoldDataSpectra_wNLO_etabin.exe "test"
##OLD JOHNS THRESHS FOR |Y|<2.0, WANTED 3.0
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "01.06.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_01-06-19_jetPlots_0.0eta2.0" "Bayes_HLT60Thresh75"

#./bayesUnfoldDataSpectra_wNLO_etabin.exe "01.31.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_01-31-19_jetPlots_0.0eta3.0" "Bayes_HLT60Thresh75"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "01.28.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_01-28-19_jetPlots_0.0eta3.0" "Bayes_HLT60Thresh80"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "01.29.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_01-29-19_jetPlots_0.0eta3.0" "Bayes_HLT60Thresh85"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "01.29.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_01-29-19_jetPlots__0.0eta3.0_1" "Bayes_HLT60Thresh90"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "01.30.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_01-30-19_jetPlots_0.0eta3.0" "Bayes_HLT60Thresh90HLT80Thresh110"

#### DO ME NEXT
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.04.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-04-19_jetPlots_0.0eta3.0" "Bayes_noDataCovMat_kNoErrs"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.04.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-04-19_jetPlots_0.0eta3.0" "Bayes_noDataCovMat_kHistErrs"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.04.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-04-19_jetPlots_0.0eta3.0" "Bayes_noDataCovMat_kCovErrs"


#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.04.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-04-19_jetPlots_0.0eta3.0" "Bayes_wDataCovMat_kNoErrs"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.04.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-04-19_jetPlots_0.0eta3.0" "Bayes_wDataCovMat_kHistErrs"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.04.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-04-19_jetPlots_0.0eta3.0" "Bayes_wDataCovMat_sqrtlumiscale_kCovErrs"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.04.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-04-19_jetPlots_0.0eta3.0" "Bayes_wDataCovMat_lumisqrd_kCovErrs"

#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.04.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-04-19_jetPlots_0.0eta3.0" "Bayes_wDataCovMat_lumisqrd_measErrsSet_kCovErrs"

#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.07.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-07-19_jetPlots_0.0eta3.0" "Bayes_wDataCovMat_psclsqrd_lumisqrd_measErrsSet_kCovErrs_LONNPDF" "smearTheory/NNPDF_NNLO_01.10.19_spl3wgts_gaussSmear_plots/" "NNPDF_NNLO_01.10.19_spl3wgts_gaussSmear_"

#scp2Serin output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_wDataCovMat_psclsqrd_lumisqrd_measErrsSet_kCovErrs_LONNPDF_NLOMC_wNP_00eta05.pdf

./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.07.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-07-19_jetPlots_0.0eta3.0" "Bayes_wDataCovMat_psclsqrd_lumisqrd_measErrsSet_kCovErrs_NLONNPDF" "smearTheory/NNPDF_NNLO_02.14.19_v2_spl3wgts_gaussSmear_plots/" "NNPDF_NNLO_02.14.19_v2_spl3wgts_gaussSmear_"

scp2Serin output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_wDataCovMat_psclsqrd_lumisqrd_measErrsSet_kCovErrs_NLONNPDF_NLOMC_wNP_00eta05.pdf

#rooUnfoldCompile SVDUnfoldDataSpectra_wNLO_etabin.C

#./SVDUnfoldDataSpectra_wNLO_etabin.exe "test"

#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_*_*_NLOMC_wNP_??eta??*.*"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_HLT60Thresh*_NLOMC_wNP_00eta05*"