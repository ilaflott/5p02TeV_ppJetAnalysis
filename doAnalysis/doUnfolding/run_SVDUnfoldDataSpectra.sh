#!/bin/bash


##CONST
jtID=1
##CONST


R=$1
data_condorDir=$2
data_date=$3
dataetaBin=$4
MC_condorDir=$5
MC_date=$6
MCetaBin=$7
etaBinOut=$8
kReg=$9
useSimpleBinning=${10}

#### unfoldDataSpectra
echo ""
echo "running SVDUnfoldDataSpectra"
echo ""

./SVDUnfoldDataSpectra.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${MC_condorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${MC_date}_JERS_${MCetaBin} unf_HPtJetTrig_w_Py8_${etaBinOut} ${jtID} ${useSimpleBinning} ${kReg} 

return