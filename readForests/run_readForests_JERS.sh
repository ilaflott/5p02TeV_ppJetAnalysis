#!/bin/bash

if [ $# -eq 3 ]
then
    R=$1
    etaLo=$2
    etaHi=$3
    echo "running jobs for R=${R} w/ ${etaLo} < |eta| < ${etaHi}"
else
    echo "usage:"
    echo "source run_readForests_JERS.sh [R=3,4] [etaLo=\"0.0\"] [etaHi=\"2.0\"]"
    return
fi

### semiOfficial Py8

echo ""
echo "submitting Py8 ppMC JERS job(s)"
echo ""

#source condorSubmit_readForests.sh     <readForestsCode>  <NJobs> <NFilesPerJob> <startFilePos> <filelistIn> <radius> <jetType> <debug> [<etaCutLo> <etaCutHi>]
#source  condorSubmit_readForests.sh readForests_ppMC_JERS -1 80 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0  ${etaLo} ${etaHi}
source  condorSubmit_readForests.sh readForests_ppMC_JERS -1 120 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0  ${etaLo} ${etaHi}

echo ""
echo "done submitting JERS job(s)"
echo ""

condor_q ilaflott

return

















