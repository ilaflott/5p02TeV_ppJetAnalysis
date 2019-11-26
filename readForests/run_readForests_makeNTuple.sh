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
    echo "source run_readForests_makeNTuple.sh [R=3,4] [etaLo=\"0.0\"] [etaHi=\"2.0\"] [type=\"ppData\"]"
    return
fi

##########################################################################################################################################################################

#if [[ $datatype == "ppData" ]]
#then
##### Jet80, then LowerJets
echo ""
echo "submitting ppData Jet80 ppData job(s)"
echo ""

#source condorSubmit_readForests.sh readForests_ppData_makeNTuple 2 2 0 filelists/5p02TeV_HighPtJet80_forests.txt     ${R} PF 1 ${etaLo} ${etaHi} ## DEBUG
#source condorSubmit_readForests.sh readForests_ppData_makeNTuple -1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt     ${R} PF 0 ${etaLo} ${etaHi}
#source condorSubmit_readForests.sh readForests_ppData_makeNTuple -1 15 0 filelists/5p02TeV_HighPtJet80_forests.txt     ${R} PF 0 ${etaLo} ${etaHi}
#source condorSubmit_readForests.sh readForests_ppData_makeNTuple -1 20 0 filelists/5p02TeV_HighPtJet80_forests.txt     ${R} PF 0 ${etaLo} ${etaHi}
#source condorSubmit_readForests.sh readForests_ppData_makeNTuple -1 50 0 filelists/5p02TeV_HighPtJet80_forests.txt     ${R} PF 0 ${etaLo} ${etaHi}

#echo ""
#echo "submitting ppData LowerJets job(s)"
#echo ""

#source condorSubmit_readForests.sh readForests_ppData_makeNTuple 2 2 0 filelists/5p02TeV_HighPtLowerJets_forests.txt ${R} PF 1 ${etaLo} ${etaHi} ## DEBUG
#source condorSubmit_readForests.sh readForests_ppData_makeNTuple -1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt ${R} PF 0 ${etaLo} ${etaHi}
source condorSubmit_readForests.sh readForests_ppData_makeNTuple -1 20 0 filelists/5p02TeV_HighPtLowerJets_forests.txt ${R} PF 0 ${etaLo} ${etaHi}
#source condorSubmit_readForests.sh readForests_ppData_makeNTuple 10 20 0 filelists/5p02TeV_HighPtLowerJets_forests.txt ${R} PF 0 ${etaLo} ${etaHi}
#source condorSubmit_readForests.sh readForests_ppData_makeNTuple 20 20 200 filelists/5p02TeV_HighPtLowerJets_forests.txt ${R} PF 0 ${etaLo} ${etaHi}
#source condorSubmit_readForests.sh readForests_ppData_makeNTuple -1 50 0 filelists/5p02TeV_HighPtLowerJets_forests.txt ${R} PF 0 ${etaLo} ${etaHi}

#echo ""
#echo "done submitting ppData MinBias6 job(s)"
#echo ""

#source condorSubmit_readForests.sh readForests_ppData_makeNTuple 1 1 0 filelists/5p02TeV_MinBias6_758p3_forests.txt ${R} PF 1 ${etaLo} ${etaHi} ## DEBUG
#source condorSubmit_readForests.sh readForests_ppData_makeNTuple 1 3 0 filelists/5p02TeV_MinBias6_758p3_forests.txt ${R} PF 0 ${etaLo} ${etaHi} ## TRY IT?!
#source condorSubmit_readForests.sh readForests_ppData_makeNTuple -1 1 0 filelists/5p02TeV_MinBias6_758p3_forests.txt ${R} PF 0 ${etaLo} ${etaHi}

echo ""
echo "done submitting ppData makeNTuple job(s)"
echo ""

#elif [[ $datatype == "ppMC" ]]
#then
#    ####### semiOfficial Py8
#    echo ""
#    echo "submitting Py8 ppMC job(s)"
#    echo ""
#    
#    #source condorSubmit_readForests.sh readForests_ppMC_makeNTuple 2 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0 ${etaLo} ${etaHi}
#    
#    source condorSubmit_readForests.sh readForests_ppMC_makeNTuple -1 200 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0 ${etaLo} ${etaHi}
#    #source condorSubmit_readForests.sh readForests_ppMC_makeNTuple -1 300 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0 ${etaLo} ${etaHi}
#
#    
#    
#    echo ""
#    echo "done submitting ppMC makeNTuple job(s)"
#    echo ""
#fi


return
