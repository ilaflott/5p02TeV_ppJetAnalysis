#!/bin/bash

if [ $# -eq 0 ]
then
    echo ""
    echo "give me a radius (3,4... etc.), and flags for jtID and compilation (0 or 1 only)"
    echo ""

    echo "e.g...."
    echo "source run_unfoldSpectra 3 1 1"
    echo ""

    return
fi

R=$1
jtID=$2
compile=$3

if [ $compile -eq 1 ]
then
    echo ""
    echo "compiling unfoldDataSpectra"
    echo ""
    rooUnfoldCompile unfoldDataSpectra.C

    echo ""
    echo "done compiling unfoldDataSpectra"
    echo ""
    sleep 0.5s
else
    echo ""
    echo "skipping compilation..."
    echo ""
fi





#### unfoldDataSpectra

echo ""
echo "running unfoldDataSpectra"
echo ""

./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak${R}PFJets_05-25-17_jetPlots_0.0eta0.5  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_05-25-17_JERS_0.0eta0.5     HPtJetTrig_Py8_unfData_00eta05   ${jtID}
./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak${R}PFJets_05-25-17_jetPlots_0.5eta1.0  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_05-25-17_JERS_0.5eta1.0     HPtJetTrig_Py8_unfData_05eta10   ${jtID}
#./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak${R}PFJets_05-25-17_jetPlots_1.0eta1.5  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_05-25-17_JERS_1.0eta1.5     HPtJetTrig_Py8_unfData_10eta15  ${jtID}
#./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak${R}PFJets_05-25-17_jetPlots_1.5eta2.0  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_05-25-17_JERS_1.5eta2.0     HPtJetTrig_Py8_unfData_15eta20  ${jtID}

echo ""
echo "done running unfoldDataSpectra"
echo ""
sleep 0.5s





#### doMCClosureTests
#echo ""
#echo "compiling doMCClosureTests"
#echo ""
#rooUnfoldCompile doMCClosureTests.C
#
#echo ""
#echo "done compiling doMCClosureTests"
#echo ""
#sleep 0.5s
#
##echo ""
##echo "running doMCClosureTests"
##echo ""
#
##./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_05-25-17_JERS_0.0eta0.5     Py8_closureTest_00eta05_0jtID  0
##./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_05-25-17_JERS_0.5eta1.0     Py8_closureTest_05eta10_0jtID  0
##./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_05-25-17_JERS_1.0eta1.5     Py8_closureTest_10eta15_0jtID  0
##./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_05-25-17_JERS_1.5eta2.0     Py8_closureTest_15eta20_0jtID  0

##./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_05-25-17_JERS_0.0eta0.5     Py8_closureTest_00eta05_wjtID  1
##./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_05-25-17_JERS_0.5eta1.0     Py8_closureTest_05eta10_wjtID  1
##./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_05-25-17_JERS_1.0eta1.5     Py8_closureTest_10eta15_wjtID  1
##./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_05-25-17_JERS_1.5eta2.0     Py8_closureTest_15eta20_wjtID  1


echo ""
echo "done!"
echo ""

return