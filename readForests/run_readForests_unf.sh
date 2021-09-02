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
echo "submitting Py8 ppMC unf job(s)"
echo ""

#source condorSubmit_readForests.sh     <readForestsCode>  <NJobs> <NFilesPerJob> <startFilePos> <filelistIn> <radius> <jetType> <debug> [<etaCutLo> <etaCutHi>]

### test/debug job, only 1 and small.
#source  condorSubmit_readForests.sh readForests_ppMC_unf 1 50 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0  ${etaLo} ${etaHi}

### run jobs in one go
#source  condorSubmit_readForests.sh readForests_ppMC_unf -1 200 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0  ${etaLo} ${etaHi}

### run jobs in two parts
#source  condorSubmit_readForests.sh readForests_ppMC_unf 17 200 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0  ${etaLo} ${etaHi}
source  condorSubmit_readForests.sh readForests_ppMC_unf 18 200 3400 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0  ${etaLo} ${etaHi}

echo ""
echo "done submitting unf job(s)"
echo ""

condor_q ilaflott

return

















