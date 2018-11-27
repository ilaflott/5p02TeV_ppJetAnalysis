#list="filelists/5p02TeV_HighPtLowerJets_forests.txt"

#!/bin/bash

rootcompile readForests_ppData_findEvt.C

#JET80 JOBS
list="filelists/5p02TeV_HighPtJet80_forests.txt"
#JET80 JOBS

dir="ppData_HighPtJet80_ak4PFJets_11-19-18_findEvt_0.0eta0.1"
R=4
etaLo=0.0
etaHi=0.1
ARRAY=(190 210 270)
for i in ${ARRAY[@]}; do
    echo "i=$i"
    source condorSubmit_readForests.sh readForests_ppData_findEvt 1 10 $i $list  ${R} PF 0 ${etaLo} ${etaHi} ${dir}
    sleep 30s
done


#return

#LOWERJETS JOBS
list="filelists/5p02TeV_HighPtLowerJets_forests.txt"
#LOWERJETS JOBS

dir="ppData_HighPtLowerJets_ak4PFJets_11-19-18_findEvt_0.0eta0.1"
R=4
etaLo=0.0
etaHi=0.1
ARRAY=(120 360 380 420 60)
for i in ${ARRAY[@]}; do
    echo "i=$i"
    source condorSubmit_readForests.sh readForests_ppData_findEvt 1 20 $i $list  ${R} PF 0 ${etaLo} ${etaHi} ${dir}
    sleep 60s
done
return

