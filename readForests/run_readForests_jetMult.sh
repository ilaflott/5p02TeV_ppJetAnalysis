#!/bin/bash

if [ $# -eq 3 ]
then
    R=$1
    etaLo=$2
    etaHi=$3
    echo "running jobs for R=${R} w/ ${etaLo} < |eta| < ${etaHi}"
else
    echo "usage:"
    echo "source run_readForests_jetMult.sh [R=3,4] [etaLo=\"0.0\"] [etaHi=\"2.0\"]"
    return
fi

##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################

##### Jet80, then LowerJets
echo ""
echo "submitting jet80/lowerJets ppData job(s)"
echo ""

source condorSubmit_readForests.sh readForests_ppData_jetMult -1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt     ${R} PF 0 ${etaLo} ${etaHi}

echo ""
echo "done submitting ppData jetMult Jet80 job(s)"
echo ""
#sleep 5s

source condorSubmit_readForests.sh readForests_ppData_jetMult -1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt ${R} PF 0 ${etaLo} ${etaHi}

echo ""
echo "done submitting ppData jetMult LowerJets job(s)"
echo ""
#sleep 5s


####### semiOfficial Py8
#echo ""
#echo "submitting Py8 ppMC job(s)"
#echo ""
#
#source condorSubmit_readForests.sh readForests_ppMC_jetMult -1 200 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0 ${etaLo} ${etaHi}
#
#echo ""
#echo "done submitting ppMC jetMult job(s)"
#echo ""
#sleep 5s

return