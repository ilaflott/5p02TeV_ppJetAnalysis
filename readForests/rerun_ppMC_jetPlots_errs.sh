#!/bin/bash


list="filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt"

dir="ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_11-14-18_JERS_0.0eta2.0"
R=4
etaLo=0.0
etaHi=2.0
ARRAY=(6360 6720 6840)
for i in ${ARRAY[@]}; do
    echo "i=$i"
    source condorSubmit_readForests.sh readForests_ppMC_JERS 1 120 $i $list  ${R} PF 0 ${etaLo} ${etaHi} ${dir}
done






return


#dir="ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_11-06-18_JERS_0.0eta0.5"
#R=4
#etaLo=0.0
#etaHi=0.5
#ARRAY=(5775 6475)
#for i in ${ARRAY[@]}; do
#    echo "i=$i"
#    source condorSubmit_readForests.sh readForests_ppMC_JERS 1 80 $i $list  ${R} PF 0 ${etaLo} ${etaHi} ${dir}
#done


dir="ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_11-06-18_JERS_0.5eta1.0"
R=4
etaLo=0.5
etaHi=1.0
ARRAY=(320 560 720 960)
for i in ${ARRAY[@]}; do
    echo "i=$i"
    source condorSubmit_readForests.sh readForests_ppMC_JERS 1 80 $i $list  ${R} PF 0 ${etaLo} ${etaHi} ${dir}
done

sleep 10m


dir="ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_11-06-18_JERS_1.0eta1.5"
R=4
etaLo=1.0
etaHi=1.5
ARRAY=(4000 4080 4480 4640 4960 5120 5360 5760 5840 5920 6000 6080 6240 6320 6480 6720 6800)
for i in ${ARRAY[@]}; do
    echo "i=$i"
    source condorSubmit_readForests.sh readForests_ppMC_JERS 1 80 $i $list  ${R} PF 0 ${etaLo} ${etaHi} ${dir}
done

sleep 10m

dir="ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_11-06-18_JERS_1.5eta2.0"
R=4
etaLo=1.5
etaHi=2.0
ARRAY=(1120 1200 1280 1440 1520 1680 1760 1920 2160 2240 2320 2400 2480 2560 2720 2880 3040 3120 3280 4400 4480 4560 4640 4720 4880 4960 5120 5200 5360 5520 720 800 960)
for i in ${ARRAY[@]}; do
    echo "i=$i"
    source condorSubmit_readForests.sh readForests_ppMC_JERS 1 80 $i $list  ${R} PF 0 ${etaLo} ${etaHi} ${dir}
done