#!/bin/bash


##CONST
jtID=1
##CONST

R=$1
data_condorDir=$2
data_date=$3
MC_condorDir=$4
MC_date=$5
etaBin=$6
etaBinOut=$7
useSimpleBinning=$8
kReg=$9



#### unfoldDataSpectra
echo ""
echo "running SVDUnfoldDataSpectra"
echo ""

./SVDUnfoldDataSpectra.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${etaBin} ${MC_condorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${MC_date}_JERS_${etaBin} unf_HPtJetTrig_w_Py8_${etaBinOut} ${jtID} ${useSimpleBinning} ${kReg} 

return