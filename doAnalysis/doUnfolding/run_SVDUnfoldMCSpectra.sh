#!/bin/bash


##CONST
jtID=1
##CONST

if [ $# -eq 7 ]
then
    R=$1
    condorDir=$2
    date_output=$3
    etaBin=$4
    etaBinOut=$5
    kReg=$6
    useSimpleBinning=$7
else
    echo "usage:"
    echo "source run_SVDUnfoldMCSpectra.sh [R=3,4] [condorDir] [date_output] [etaBin] [etaBinOut] [kReg] [useSimpleBinning]"
    return
fi

#### unfoldMCSpectra
echo ""
echo "running unfoldMCSpectra"
echo ""

./SVDUnfoldMCSpectra.exe ${condorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_JERS_${etaBin}  Py8_closureTest_${etaBinOut} ${jtID} ${useSimpleBinning} ${kReg} 

return