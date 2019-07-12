#!/bin/bash

echo ""
echo "rename/hadd'ing..."
echo ""

#localDir=/Users/ilaflott/Working/Dropbox/localAnalysis/src/5p02TeV_ppJetAnalysis-master/doAnalysis/printPlots/outputCondor

MM="06"
DD="25"
YY="19"


#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppMC" "semifinal_doublecheck"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.0"     "ppMC" "semifinal_doublecheck"



####source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "JohnSbins"
####source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "Chrisbins"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData" "SMPbins_withJECsys_semifinal_doublecheck"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "unf_0.0eta2.0"      "ppMC" "semifinal"

##mv ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_06-26-19_MCEff_0.0eta2.0 ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_06-25-19_MCEff_0.0eta2.0
source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "MCEff_0.0eta2.0"    "ppMC" "semifinal"










echo ""
echo "done rename/hadd'ing."
echo ""

return
