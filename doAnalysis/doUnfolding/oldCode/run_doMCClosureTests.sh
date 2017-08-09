#!/bin/bash


##CONST
date_data="05-25-17"
date_condor="5.25.17"
outCondorDir="${date_condor}_outputCondor"
##CONST


if [ $# -ne 5 ]
then
    echo "bad usage."
    echo "source run_unfoldMCSpectra <R> <jtID> <compile> <SVDkReg> <etaBin#>"
    echo ""

    echo "e.g."
    echo "source run_unfoldMCSpectra 3 1 1 20 1"
    echo ""

    return
fi

R=$1
jtID=$2
compile=$3
kReg=$4
etaBin=$5

if [ $compile -eq 1 ]
then
    echo ""
    echo "compiling unfoldMCSpectra"
    echo ""
    rooUnfoldCompile unfoldMCSpectra.C

    echo ""
    echo "done compiling unfoldMCSpectra"
    echo ""
    sleep 0.5s
else
    echo ""
    echo "skipping compilation..."
    echo ""
fi

#### unfoldMCSpectra
echo ""
echo "running unfoldMCSpectra"
echo ""

if [ $etaBin -eq 1 ]
then
    ./unfoldMCSpectra.exe   ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_0.0eta0.5   Py8_closureTest_00eta05_wjtID  1
elif [ $etaBin -eq 2 ]
then
    ./unfoldMCSpectra.exe   ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_0.5eta1.0   Py8_closureTest_05eta10_wjtID  1
else
    echo ""
    echo "etaBin not found; exit"
    echo ""
    return
fi




#echo ""
#echo "done running unfoldDataSpectra"
#echo ""
#sleep 0.5s
#
#
#
#
#
##### unfoldMCSpectra
##echo ""
##echo "compiling unfoldMCSpectra"
##echo ""
##rooUnfoldCompile unfoldMCSpectra.C
##
##echo ""
##echo "done compiling unfoldMCSpectra"
##echo ""
##sleep 0.5s
##
###echo ""
###echo "running unfoldMCSpectra"
###echo ""
##
###./unfoldMCSpectra.exe   ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_0.0eta0.5     Py8_closureTest_00eta05_0jtID  0
###./unfoldMCSpectra.exe   ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_0.5eta1.0     Py8_closureTest_05eta10_0jtID  0
###./unfoldMCSpectra.exe   ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_1.0eta1.5     Py8_closureTest_10eta15_0jtID  0
###./unfoldMCSpectra.exe   ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_1.5eta2.0     Py8_closureTest_15eta20_0jtID  0
#

###./unfoldMCSpectra.exe   ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_1.0eta1.5     Py8_closureTest_10eta15_wjtID  1
###./unfoldMCSpectra.exe   ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_1.5eta2.0     Py8_closureTest_15eta20_wjtID  1


echo ""
echo "done!"
echo ""

return