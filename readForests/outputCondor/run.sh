#!/bin/bash

echo ""
echo "rename/hadd'ing..."
echo ""

#localDir=/Users/ilaflott/Working/Dropbox/localAnalysis/src/5p02TeV_ppJetAnalysis-master/doAnalysis/printPlots/outputCondor

###SEMIFINAL DATE DO NOT REMOVE
#MM="06"
#DD="25"
#YY="19"

#MM="08"
#DD="06"
#YY="19"

#MM="11"
#DD="26"
#YY="19"

#MM="03"
#DD="18"
#YY="20"

#MM="12"
#DD="02"
#YY="20"

#MM="01"
#DD="05"
#YY="21"

#MM="05"
#DD="08"
#YY="21"

MM="07"
DD="22"
YY="21"




#mv                ppData_HighPtJet80_ak4PFJets_01-06-21_jetTrigEff_0.0eta5.1                 ppData_HighPtJet80_ak4PFJets_01-05-21_jetTrigEff_0.0eta5.1
#mv            ppData_HighPtLowerJets_ak4PFJets_01-06-21_jetTrigEff_0.0eta5.1             ppData_HighPtLowerJets_ak4PFJets_01-05-21_jetTrigEff_0.0eta5.1
#mv ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_01-07-21_JERS_0.0eta2.5        ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_01-05-21_JERS_0.0eta2.5      




#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.5" "ppData" "SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.5" "ppData" "SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_PFMETplotsOnly_semifinalv4"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetTrigEff_0.0eta5.1" "ppData" "wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.5" "ppMC" "SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_PFMETplotsOnly_semifinalv4"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.5" "ppMC" "SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.5" "ppMC" "SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_tightJetID_semifinalv4"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.5" "ppMC" "SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "unf_0.0eta2.5" "ppMC" "SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "unf_0.0eta3.0" "ppMC" "SMPbins_withjetyrespmat_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.5" "ppData" "SMPbins_wtype1PFMETFrac0p3_semifinalv4"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.5" "ppMC" "SMPbins_wtype1PFMETFrac0p3_semifinalv4"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.5" "ppMC" "SMPbins_wPFMETFrac_genbinsonly_semifinalv4"

mv ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_08-16-21_JERS_0.0eta2.5 ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_07-22-21_JERS_0.0eta2.5
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.5" "ppMC" "SMPbins_wPFMETFrac0p5_wgenDR0p2_genbinsonly_semifinalv4"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.5" "ppMC" "SMPbins_wPFMET_GENRECOmatchREDOwDR0p2_genbinsonly_semifinalv4"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.5" "ppMC" "SMPbins_wPFMET_Nvtxeq1_evtvr0p1_GENRECOmatchREDOwDR0p2_genbinsonly_semifinalv4"
source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.5" "ppMC" "SMPbins_wPFMET_novzwgt_Nvtxeq1_evtvr0p1_GENRECOmatchREDOwDR0p2_genbinsonly_semifinalv4"
return














#mv ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_03-26-20_JERS_0.0eta2.5/ ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_03-18-20_JERS_0.0eta2.5/ 
#mv ppData_HighPtJet80_ak4PFJets_11-20-20_jetPlots_0.0eta2.5_1     ppData_HighPtJet80_ak4PFJets_03-18-20_jetPlots_0.0eta2.5    
#mv ppData_HighPtLowerJets_ak4PFJets_11-20-20_jetPlots_0.0eta2.5_1 ppData_HighPtLowerJets_ak4PFJets_03-18-20_jetPlots_0.0eta2.5
#mv ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-20-20_jetPlots_0.0eta2.0 ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-19-20_jetPlots_0.0eta2.0
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta1.0" "ppMC" "semifinal_novzweights"





#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "useTupel_withLumiCorr"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "makeNTupel_0.0eta2.0" "ppData" "useTupel_withLumiCorr"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "SMPbins_JECsysv2_withLumiCorr_semifinalv2"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "SMPbins_JECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "SMPbins_JECsysv2_withLumiCorr_semifinal"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "CHRISbins_JECsysv2_withLumiCorr_semifinal"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetTrigEff_0.0eta2.0" "ppData" "forHLTTurnOns_secondAttempt"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetTrigEff_0.0eta5.1" "ppData" "semifinalv2"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.5" "ppData" "SMPbins_withjety_semifinalv4"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.5" "ppData" "SMPbins_withjety_eMaxJetID_semifinalv4"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.5" "ppData" "SMPbins_withjety_eSumJetID_semifinalv4"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.5" "ppData" "SMPbins_withjety_hotspotsv3_2k17tightJetID_semifinalv4"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.5" "ppData" "SMPbins_withjety_wHotSpotMask_JetID_semifinalv4"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.5" "ppData" "SMPbins_withjety_wMETFrac0p3Cut_JetID_semifinalv4"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.5" "ppData" "SMPbins_withjety_wHBHEIsoNoiseCut_JetID_semifinalv4"


#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "JOHNbins_withJECsys_withLumiCorr_semifinal"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta1.0" "ppData" "lumiRunStudy_diffDupeSkip"


####source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "Chrisbins"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "SMPbins_withJECsys_semifinal_doublecheck"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0"      "ppMC" "noDRCut_tightJetID_semifinal"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "MCEff_0.0eta2.0"      "ppMC" "noDRCut_semifinal"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "unf_0.0eta2.5"      "ppMC" "genetabins_semifinalv3"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.5"      "ppMC" "semifinalv3"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.5"      "ppMC" "SMPbins_withjety_semifinalv3"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.5"      "ppMC" "SMPbins_withjety_2k17tightJetID_semifinalv4"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "unf_0.0eta2.5"      "ppMC" "SMPbins_withjety_semifinalv3"

##mv ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_06-26-19_MCEff_0.0eta2.0 ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_06-25-19_MCEff_0.0eta2.0
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "MCEff_0.0eta2.0"    "ppMC" "semifinal"










echo ""
echo "done rename/hadd'ing."
echo ""

return
