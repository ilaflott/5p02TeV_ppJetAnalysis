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
MM="11"
DD="26"
YY="19"


#mv ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_07-31-19_jetPlots_0.0eta1.0 ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_07-29-19_jetPlots_0.0eta1.0
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta1.0" "ppMC" "semifinal_novzweights"




#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "useTupel_withLumiCorr"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "makeNTupel_0.0eta2.0" "ppData" "useTupel_withLumiCorr"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "SMPbins_JECsysv2_withLumiCorr_semifinalv2"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "SMPbins_JECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "SMPbins_JECsysv2_withLumiCorr_semifinal"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "CHRISbins_JECsysv2_withLumiCorr_semifinal"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetTrigEff_0.0eta2.0" "ppData" "forHLTTurnOns_secondAttempt"
source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetTrigEff_0.0eta5.1" "ppData" "semifinalv2"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetTrigEff_0.0eta5.1" "ppData" "semifinalv2_LowHLT40Thresh"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "JOHNbins_withJECsys_withLumiCorr_semifinal"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta1.0" "ppData" "lumiRunStudy_diffDupeSkip"


####source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "Chrisbins"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "SMPbins_withJECsys_semifinal_doublecheck"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0"      "ppMC" "noDRCut_tightJetID_semifinal"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "MCEff_0.0eta2.0"      "ppMC" "noDRCut_semifinal"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "unf_0.0eta2.0"      "ppMC" "noDRCut_semifinal"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.0"      "ppMC" "noDRCut_semifinal"

##mv ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_06-26-19_MCEff_0.0eta2.0 ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_06-25-19_MCEff_0.0eta2.0
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "MCEff_0.0eta2.0"    "ppMC" "semifinal"










echo ""
echo "done rename/hadd'ing."
echo ""

return
