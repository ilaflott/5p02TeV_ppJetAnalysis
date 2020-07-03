#!/bin/bash

#make_rooUnfold

rooUnfoldCompile bayesUnfoldDataSpectra_etabin.C

############ FULL GEN + RECO PY8 ###############


####WITH JEC SYS V2, WITH LUMI CORR, SMP BINS PT LO 56, LOWHLT40THRESH
DATADIR="03.18.20_outputCondor/ppData_HighPtJetTrig_ak4PFJets_03-18-20_jetPlots_0.0eta2.5_SMPbins_withjety_semifinalv3"
#MCDIR="03.18.20_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_03-18-20_unf_0.0eta2.5_SMPbins_withjety_semifinalv3/"
#MCDIR="03.18.20_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_03-18-20_unf_0.0eta2.5_SMPbins_withgenjety_semifinalv3/"
MCDIR="03.18.20_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_03-18-20_unf_0.0eta2.5_SMPbins_withgenjety_wJERSFs_semifinalv3/"
MCFILESTR="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
#OUTPUTTAG="Bayes_PY8_FullRECO_03.18.20_SMPbins_withJECsysv2_SMPbins_ptLo56_withLumiCorr_semifinalv3"
#OUTPUTTAG="03.18.20_Bayes_PY8_FullRECO_SMPbins_withJECsysv2_withgenjety_SMPbins_semifinalv3"
OUTPUTTAG="03.18.20_Bayes_PY8_FullRECO_SMPbins_withgenjety_wJERSFs_SMPbins_semifinalv3"

dosystunf=1
declare -a sysarr=("JEC" "JER")
#declare -a sysarr=("JER")
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
