#!/bin/bash

echo ""
echo "rename/hadd'ing..."
echo ""

#localDir=/Users/ilaflott/Working/Dropbox/localAnalysis/src/5p02TeV_ppJetAnalysis-master/doAnalysis/printPlots/outputCondor

MM="06"
DD="07"
YY="19"
#mv ppData_HighPtJet80_ak4PFJets_05-20-19_jetPlots_0.0eta2.0 ppData_HighPtJet80_ak4PFJets_05-17-19_jetPlots_0.0eta2.0
#mv ppData_HighPtLowerJets_ak4PFJets_05-20-19_jetPlots_0.0eta2.0 ppData_HighPtLowerJets_ak4PFJets_05-17-19_jetPlots_0.0eta2.0
#mv ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-20-19_JERS_0.0eta2.0 ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-17-19_JERS_0.0eta2.0

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "JohnSbins"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "Chrisbins"
source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "SMPbins_withJECsys"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppMC" "duplicate"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.0" "ppMC" "unf"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.0" "ppMC" "JER"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta3.0" "ppData" "JEC_final"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta3.0" "ppData" "wSpecialTrgPt_hists_classicDupeSkip"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta3.0" "ppMC" "JEC_final"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta3.0" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta3.0_JER" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta3.0_unf" "ppMC"


echo ""
echo "done rename/hadd'ing."
echo ""

return




