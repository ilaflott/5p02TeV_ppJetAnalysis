#!/bin/bash

echo ""
echo "rename/hadd'ing..."
echo ""

#localDir=/Users/ilaflott/Working/Dropbox/localAnalysis/src/5p02TeV_ppJetAnalysis-master/doAnalysis/printPlots/outputCondor

MM="02"
DD="07"
YY="19"


############################################################

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta0.5" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.5eta1.0" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_1.0eta1.5" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_1.5eta2.0" "ppMC"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta1.0" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_1.0eta2.0" "ppMC"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppMC"

############################################################

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta0.5_gendr0p3_recoJetID" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.5eta1.0_gendr0p3_recoJetID" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_1.0eta1.5_gendr0p3_recoJetID" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_1.5eta2.0_gendr0p3_recoJetID" "ppMC"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta1.0" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_1.0eta2.0" "ppMC"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.0_2ndBin" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS__0.0eta2.0_1" "ppMC"

############################################################

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta0.5" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.5eta1.0" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_1.0eta1.5" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_1.5eta2.0" "ppData"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta1.0" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_1.0eta2.0" "ppData"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppData"
source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta3.0" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots__0.0eta3.0_1" "ppData"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta0.5" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.5eta1.0" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_1.0eta1.5" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_1.5eta2.0" "ppMC"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta1.0" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_1.0eta2.0" "ppMC"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppMC"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta0.5" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.5eta1.0" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_1.0eta1.5" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_1.5eta2.0" "ppMC"


#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0_trgPtNoRes_NoJetID" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0_trgPtL2L3Res" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0_trgPtL2L3Res_NoJetID" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0_NoJetID" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta3.0" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.0_unf" "ppMC"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots__0.0eta2.0_1" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots__0.0eta2.0_1" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "findEvt_0.0eta0.1" "ppData"

echo ""
echo "done rename/hadd'ing."
echo ""

return








































#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta2.0" "ppMC"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.0eta0.5" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_0.5eta1.0" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_1.0eta1.5" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "jetPlots_1.5eta2.0" "ppMC"



#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "0.0eta0.5" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "0.5eta1.0" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "1.0eta1.5" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "1.5eta2.0" "ppData"



#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "2.0eta2.5" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "2.5eta3.0" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "3.2eta4.7" "ppMC"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.0" "ppMC"



## ppDataMC w/ larger eta regions
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "0.0eta2.0_HLTCalo" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "0.0eta2.0_HLTPF" "ppData"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "2.0eta3.0" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "2.0eta3.0" "ppMC"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "3.2eta4.7" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "3.2eta4.7" "ppMC"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "0.0eta4.7" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "0.0eta4.7" "ppMC"

## ppData analysis eta bins
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "0.0eta0.5" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "0.5eta1.0" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "1.0eta1.5" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "1.5eta2.0" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "2.0eta2.5" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "2.5eta3.0" "ppData"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "3.2eta4.7" "ppData"



#scp -r 4scp/ppMC* ${localCPU}:${localDir}/${MM}.${DD}.${YY}_outputCondor/
#scp -r 4scp/ppData* ${localCPU}:${localDir}/${date_scratch}_outputCondor/ppData_CaloHLT
#scp -r 4scp/ppData* ${localCPU}:${localDir}/${MM}.${DD}.${YY}_outputCondor/




echo ""
echo "done rename/hadd'ing."
echo ""

return