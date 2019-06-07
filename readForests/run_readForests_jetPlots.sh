#!/bin/bash

if [ $# -eq 4 ]
then
    R=$1
    etaLo=$2
    etaHi=$3
    datatype=$4
    echo "running jobs for ${type}, R=${R} w/ ${etaLo} < |eta| < ${etaHi}"
else
    echo "usage:"
    echo "source run_readForests_jetPlots.sh [R=3,4] [etaLo=\"0.0\"] [etaHi=\"2.0\"] [type=\"ppData\"]"
    return
fi

##########################################################################################################################################################################

if [[ $datatype == "ppData" ]]
then
##### Jet80, then LowerJets
    echo ""
    echo "submitting jet80+lowerJets and/or MinBias ppData job(s)"
    echo ""

    #source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 1 0 filelists/5p02TeV_HighPtJet80_forests.txt     ${R} PF 1 ${etaLo} ${etaHi} ## DEBUG
    #source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 15 0 filelists/5p02TeV_HighPtJet80_forests.txt     ${R} PF 0 ${etaLo} ${etaHi}
    source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 25 0 filelists/5p02TeV_HighPtJet80_forests.txt     ${R} PF 0 ${etaLo} ${etaHi}
    #source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 50 0 filelists/5p02TeV_HighPtJet80_forests.txt     ${R} PF 0 ${etaLo} ${etaHi}
    
    echo ""
    echo "done submitting ppData jetPlots Jet80 job(s)"
    echo ""
    
    #source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt ${R} PF 1 ${etaLo} ${etaHi} ## DEBUG
    #source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 15 0 filelists/5p02TeV_HighPtLowerJets_forests.txt ${R} PF 0 ${etaLo} ${etaHi}
    source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 25 0 filelists/5p02TeV_HighPtLowerJets_forests.txt ${R} PF 0 ${etaLo} ${etaHi}
    #source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 50 0 filelists/5p02TeV_HighPtLowerJets_forests.txt ${R} PF 0 ${etaLo} ${etaHi}
    
    echo ""
    echo "done submitting ppData jetPlots LowerJets job(s)"
    echo ""

    #source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 1 0 filelists/5p02TeV_MinBias6_758p3_forests.txt ${R} PF 1 ${etaLo} ${etaHi} ## DEBUG
    #source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 3 0 filelists/5p02TeV_MinBias6_758p3_forests.txt ${R} PF 0 ${etaLo} ${etaHi} ## TRY IT?!
    #source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 1 0 filelists/5p02TeV_MinBias6_758p3_forests.txt ${R} PF 0 ${etaLo} ${etaHi}
    
    echo ""
    echo "done submitting ppData jetPlots MinBias job(s)"
    echo ""

elif [[ $datatype == "ppMC" ]]
then
    ####### semiOfficial Py8
    echo ""
    echo "submitting Py8 ppMC job(s)"
    echo ""
    
    source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 175 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0 ${etaLo} ${etaHi}
    ##100 files per job seems to produce errors in about 1-2 dozen files each time
    #source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 100 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0 ${etaLo} ${etaHi}
    
    echo ""
    echo "done submitting ppMC jetPlots job(s)"
    echo ""
fi


return