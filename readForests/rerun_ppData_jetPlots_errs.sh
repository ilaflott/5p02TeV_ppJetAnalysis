#!/bin/bash

rootcompile readForests_ppData_jetPlots.C

#JET80 JOBS
list="filelists/5p02TeV_HighPtLowerJets_forests.txt"
#JET80 JOBS


dir="ppData_HighPtLowerJets_ak4PFJets_11-14-18_jetPlots_0.0eta2.0_trgPtL2L3Res_NoJetID"
R=4
etaLo=0.0
etaHi=2.0
ARRAY=(430)
for i in ${ARRAY[@]}; do
    echo "i=$i"
    source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 10 $i $list  ${R} PF 0 ${etaLo} ${etaHi} ${dir}
done
return

#
#dir="ppData_HighPtJet80_ak4PFJets_11-07-18_jetPlots_0.0eta0.5"
#R=4
#etaLo=0.0
#etaHi=0.5
#ARRAY=(140 150 180 200 220 250 280 300 310 320 340 360 380 390)
#for i in ${ARRAY[@]}; do
#    echo "i=$i"
#    source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 10 $i $list  ${R} PF 0 ${etaLo} ${etaHi} ${dir}
#done

#echo""
#echo""
#echo""
#
#dir="ppData_HighPtJet80_ak4PFJets_11-07-18_jetPlots_1.0eta1.5"
#R=4
#etaLo=1.0
#etaHi=1.5
#ARRAY=(100 160 20 220 240 350 370)
#for i in ${ARRAY[@]}; do
#    echo "i=$i"
#    source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 10 $i $list  ${R} PF 0 ${etaLo} ${etaHi} ${dir}
#done
#
#echo""
#echo""
#echo""

#dir="ppData_HighPtJet80_ak4PFJets_11-07-18_jetPlots_1.5eta2.0"
#R=4
#etaLo=1.5
#etaHi=2.0
#ARRAY=(330 350 380)
#for i in ${ARRAY[@]}; do
#    echo "i=$i"
#    source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 10 $i $list  ${R} PF 0 ${etaLo} ${etaHi} ${dir}
#done
#
#echo""
#echo""
#echo""
#
#
#
##LOWERJETS JOBS
#list="filelists/5p02TeV_HighPtLowerJets_forests.txt"
##LOWERJETS JOBS
#
#dir="ppData_HighPtLowerJets_ak4PFJets_11-07-18_jetPlots_0.0eta0.5"
#R=4
#etaLo=0.0
#etaHi=0.5
#ARRAY=(10 120 130 140 150 170 180 200 230 250 260 280 290 30 310 330 340 360 380 400 430 70 80 90)
#for i in ${ARRAY[@]}; do
#    echo "i=$i"
#    source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 10 $i $list  ${R} PF 0 ${etaLo} ${etaHi} ${dir}
#done
#
#dir="ppData_HighPtLowerJets_ak4PFJets_11-07-18_jetPlots_0.5eta1.0"
#R=4
#etaLo=0.5
#etaHi=1.0
#ARRAY=(0 110 150 190 200 230 280 290 300 30 310 320 330 400 410)
#for i in ${ARRAY[@]}; do
#    echo "i=$i"
#    source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 10 $i $list  ${R} PF 0 ${etaLo} ${etaHi} ${dir}
#done
#
#dir="ppData_HighPtLowerJets_ak4PFJets_11-07-18_jetPlots_1.0eta1.5"
#R=4
#etaLo=1.0
#etaHi=1.5
#ARRAY=(10 130 140 160 20 240 250 270 330 340 360 40)
#for i in ${ARRAY[@]}; do
#    echo "i=$i"
#    source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 10 $i $list  ${R} PF 0 ${etaLo} ${etaHi} ${dir}
#done
#
#dir="ppData_HighPtLowerJets_ak4PFJets_11-07-18_jetPlots_1.5eta2.0"
#R=4
#etaLo=1.5
#etaHi=2.0
#ARRAY=(0 100 110 130 160 180 20 30 370 380 390 40 410 430 60 70 80)
#for i in ${ARRAY[@]}; do
#    echo "i=$i"
#    source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 10 $i $list  ${R} PF 0 ${etaLo} ${etaHi} ${dir}
#done