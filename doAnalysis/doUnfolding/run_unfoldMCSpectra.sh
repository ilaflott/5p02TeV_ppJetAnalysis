#!/bin/bash


##CONST
jtID=1
doSVD=1
doBayes=1
##CONST

R=$1
condorDir=$2
date_output=$3
etaBin=$4
etaBinOut=$5
kReg=$6
useSimpleBinning=$7



#### unfoldMCSpectra
echo ""
echo "running unfoldMCSpectra"
echo ""

./unfoldMCSpectra.exe ${condorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_JERS_${etaBin}  Py8_closureTest_${etaBinOut} ${jtID} ${kReg} ${doBayes} ${doSVD} ${useSimpleBinning}

return