#!/bin/bash

#make_rooUnfold

rooUnfoldCompile bayesUnfoldDataSpectra_etabin.C

############ FULL GEN + RECO PY8 ###############


#####THIS WAS USED FOR PREAPPROVAL
#DATADIR="01.05.21_outputCondor/ppData_HighPtJetTrig_ak4PFJets_01-05-21_jetPlots_0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"
#MCDIR="01.05.21_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_01-05-21_unf_0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4/"
#MCFILESTR="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
#OUTPUTTAG="01.05.21_Bayes_PY8_FullRECO_SMPbins_withgenjety_wJERSFs_SMPbins_wPFMET0p3_semifinalv4"

#### use crystal ball fit mean to apply scale factor. misses and fakes are included!
#DATADIR="01.05.21_outputCondor/ppData_HighPtJetTrig_ak4PFJets_01-05-21_jetPlots_0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"
#MCDIR="07.13.21_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_07-13-21_unf_0.0eta2.5_wFakesAndBkgJets_JERSFwCrysballMu_semifinalv4/"
#MCFILESTR="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
#OUTPUTTAG="07.13.21_Bayes_PY8_FullRECO_SMPbins_JERSFwCrysballMu_wUnmatchedJets_semifinalv4"
#OUTPUTTAG="07.13.21_Bayes_PY8_FullRECO_SMPbins_JERSFwCrysballMu_semifinalv4"

#### use crystal ball fit mean to apply scale factor. misses and fakes are included!
DATADIR="01.05.21_outputCondor/ppData_HighPtJetTrig_ak4PFJets_01-05-21_jetPlots_0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"
MCDIR="07.22.21_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_07-22-21_unf_0.0eta2.5_wUnMatchedAndMatchedFakesAndBkgJets_JERSFwCrysballMu_semifinalv4/"
MCFILESTR="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
#OUTPUTTAG="07.22.21_Bayes_PY8_FullRECO_SMPbins_JERSFwCrysballMu_noExtraFakeNorMiss_semifinalv4"
#OUTPUTTAG="07.22.21_Bayes_PY8_FullRECO_SMPbins_JERSFwCrysballMu_wUnmatchedFakeAndMiss_semifinalv4"
#OUTPUTTAG="07.22.21_Bayes_PY8_FullRECO_SMPbins_JERSFwCrysballMu_wUnmatchedFakeAndMiss_wMatchedFake_semifinalv4"
#OUTPUTTAG="07.22.21_Bayes_PY8_FullRECO_SMPbins_JERSFwCrysballMu_wUnmatchedFakeAndMiss_wMatchedMiss_semifinalv4"
OUTPUTTAG="07.22.21_Bayes_PY8_FullRECO_SMPbins_JERSFwCrysballMu_wUnmatchedFakeAndMiss_wMatchedFakeAndMiss_semifinalv4"


dosystunf=0
#declare -a sysarr=("JEC" "JER")
declare -a sysarr=("JER")
etabinStart=0
etabinEnd=4
sysarrlen=${#sysarr[@]}
for (( j=${etabinStart}; j<${etabinEnd}; j++ ));
do
    #./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  $j  ${MCDIR} ${MCFILESTR} "y"
    for (( i=0; i<${sysarrlen}; i++ ));
    do
	echo ""
	echo "j=$j"
	echo "${OUTPUTTAG}_${sysarr[i]}sys"
	echo "${sysarr[i]}"
	./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}_${sysarr[i]}sys $j ${MCDIR} ${MCFILESTR} "y" ${dosystunf} ${sysarr[i]}
    done
done
return



#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_LOMC_??y??.*"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_LOMC_??y??_respMat.*"
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  0  ${MCDIR} ${MCFILESTR} "eta"
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  1  ${MCDIR} ${MCFILESTR} "eta"
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  2  ${MCDIR} ${MCFILESTR} "eta"
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  3  ${MCDIR} ${MCFILESTR} "eta"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_LOMC_??eta??.*"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_LOMC_??eta??_respMat.*"

####WITH JEC SYS V2, WITH LUMI CORR, SMP BINS PT LO 56
#DATADIR="11.26.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_11-26-19_jetPlots_0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinalv2"
#MCDIR="06.25.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_06-25-19_unf_0.0eta2.0_semifinal/"
#MCFILESTR="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
#OUTPUTTAG="Bayes_PY8_FullRECO_11.26.19_SMPbins_withJECsysv2_SMPbins_ptLo56_withLumiCorr_semifinalv2"
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  0  ${MCDIR} ${MCFILESTR} "eta"
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  1  ${MCDIR} ${MCFILESTR} "eta"
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  2  ${MCDIR} ${MCFILESTR} "eta"
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  3  ${MCDIR} ${MCFILESTR} "eta"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_LOMC_??eta??.*"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_LOMC_??eta??_respMat.*"


####WITH JEC SYS V2, WITH LUMI CORR, SMP BINS PT LO 56, LOWHLT40THRESH
#DATADIR="11.26.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_11-26-19_jetPlots_0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh"
#MCDIR="06.25.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_06-25-19_unf_0.0eta2.0_semifinal/"
#MCFILESTR="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
#OUTPUTTAG="Bayes_PY8_FullRECO_11.26.19_SMPbins_withJECsysv2_SMPbins_ptLo56_withLumiCorr_semifinalv2_LowHLT40Thresh"
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  0  ${MCDIR} ${MCFILESTR} "eta"
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  1  ${MCDIR} ${MCFILESTR} "eta"
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  2  ${MCDIR} ${MCFILESTR} "eta"
#./bayesUnfoldDataSpectra_etabin.exe ${DATADIR} ${OUTPUTTAG}  3  ${MCDIR} ${MCFILESTR} "eta"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_LOMC_??eta??.*"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_LOMC_??eta??_respMat.*"







#------------------------------#
#-------------OLD--------------#
#------------------------------#


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
