#!/bin/bash

echo ""
echo "rename/hadd'ing..."
echo ""

localDir=/Users/ilaflott/Working/Dropbox/localAnalysis/src/5p02TeV_ppJetAnalysis-master/doAnalysis/printPlots/outputCondor

MM="01"
DD="09"
YY="18"


#rm -rf 4scp/pp*
#rm -rf backup/pp*
#cp -r pp* backup/.




## ppMC analysis eta bins

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.0_rec56" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "0.0eta2.0_rec56_HLTCalo" "ppData"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.0_rec64" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "0.0eta2.0_rec64_HLTCalo" "ppData"
#
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "0.0eta2.0_rec56_HLTPF" "ppData"
#
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "0.0eta2.0_rec64_HLTPF" "ppData"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.0_gendr0p2_HLTCalo" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.0_gendr0p1_HLTCalo" "ppMC"

#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.0_rec56_gen32" "ppMC"
#source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "JERS_0.0eta2.0_rec74_gen43" "ppMC"










# high mult
source renameNhadd.sh 4 "${MM}-${DD}-${YY}" "${MM}.${DD}.${YY}" "0.0eta2.0" "ppData"





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