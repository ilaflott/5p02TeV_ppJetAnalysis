#!/bin/bash

if [ $# -eq 3 ]
then
    R=$1
    etaLo=$2
    etaHi=$3
    #datatype=$4
    echo "running jobs for ${type}, R=${R} w/ ${etaLo} < |eta| < ${etaHi}"
else
    echo "usage:"
    echo "source run_readForests_findEvt.sh [R=3,4] [etaLo=\"0.0\"] [etaHi=\"2.0\"]"
    return
fi

##########################################################################################################################################################################

##### Jet80, then LowerJets
echo ""
echo "submitting jet80/lowerJets ppData job(s)"
echo ""

source condorSubmit_readForests.sh readForests_ppData_findEvt -1 20 0 filelists/5p02TeV_HighPtJet80_forests.txt     ${R} PF 0 ${etaLo} ${etaHi}

#source condorSubmit_readForests.sh readForests_ppData_findEvt 1 5 0 filelists/5p02TeV_HighPtJet80_forests.txt     ${R} PF 0 ${etaLo} ${etaHi} #DEBUG
#return #DEBUG

echo ""
echo "done submitting ppData jetPlots Jet80 job(s)"
echo ""



source condorSubmit_readForests.sh readForests_ppData_findEvt -1 20 0 filelists/5p02TeV_HighPtLowerJets_forests.txt ${R} PF 0 ${etaLo} ${etaHi}

echo ""
echo "done submitting ppData jetPlots LowerJets job(s)"
echo ""

return