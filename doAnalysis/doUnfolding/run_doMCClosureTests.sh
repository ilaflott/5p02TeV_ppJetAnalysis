#!/bin/bash


##CONST
date_data="05-25-17"
date_condor="5.25.17"
outCondorDir="${date_condor}_outputCondor"
##CONST


if [ $# -ne 5 ]
then
    echo "bad usage."
    echo "source run_doMCClosureTests <R> <jtID> <compile> <SVDkReg> <etaBin#>"
    echo ""

    echo "e.g."
    echo "source run_doMCClosureTests 3 1 1 20 1"
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
    echo "compiling doMCClosureTests"
    echo ""
    rooUnfoldCompile doMCClosureTests.C

    echo ""
    echo "done compiling doMCClosureTests"
    echo ""
    sleep 0.5s
else
    echo ""
    echo "skipping compilation..."
    echo ""
fi

#### doMCClosureTests
echo ""
echo "running doMCClosureTests"
echo ""

if [ $etaBin -eq 1 ]
then
    ./doMCClosureTests.exe   ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_0.0eta0.5     Py8_closureTest_00eta05_wjtID  1
elif [ $etaBin -eq 2 ]
then
    ./doMCClosureTests.exe   ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_0.5eta1.0     Py8_closureTest_05eta10_wjtID  1
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
##### doMCClosureTests
##echo ""
##echo "compiling doMCClosureTests"
##echo ""
##rooUnfoldCompile doMCClosureTests.C
##
##echo ""
##echo "done compiling doMCClosureTests"
##echo ""
##sleep 0.5s
##
###echo ""
###echo "running doMCClosureTests"
###echo ""
##
###./doMCClosureTests.exe   ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_0.0eta0.5     Py8_closureTest_00eta05_0jtID  0
###./doMCClosureTests.exe   ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_0.5eta1.0     Py8_closureTest_05eta10_0jtID  0
###./doMCClosureTests.exe   ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_1.0eta1.5     Py8_closureTest_10eta15_0jtID  0
###./doMCClosureTests.exe   ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_1.5eta2.0     Py8_closureTest_15eta20_0jtID  0
#

###./doMCClosureTests.exe   ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_1.0eta1.5     Py8_closureTest_10eta15_wjtID  1
###./doMCClosureTests.exe   ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_1.5eta2.0     Py8_closureTest_15eta20_wjtID  1


echo ""
echo "done!"
echo ""

return