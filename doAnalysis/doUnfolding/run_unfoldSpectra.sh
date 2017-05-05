#!/bin/bash

### unfoldDataSpectra
echo ""
echo "compiling unfoldDataSpectra"
echo ""
rooUnfoldCompile unfoldDataSpectra.C

echo ""
echo "done compiling. sleep."
echo ""
sleep 3s

echo ""
echo "running unfoldDataSpectra"
echo ""

./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak4PFJets_05-04-17_jetPlots_  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_     HPtJetTrig_Py8_unfoldedData_eta0to10   0
./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak4PFJets_05-04-17_jetPlots_  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_     HPtJetTrig_Py8_unfoldedData_eta10to20  0
./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak4PFJets_05-04-17_jetPlots_  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_     HPtJetTrig_Py8_unfoldedData_eta20to30  0
./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak4PFJets_05-04-17_jetPlots_  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_     HPtJetTrig_Py8_unfoldedData_eta32to47  0

echo ""
echo "done running unfoldDataSpectra. sleep."
echo ""
sleep 3s


### doMCClosureTests
echo ""
echo "compiling doMCClosureTests"
echo ""
rooUnfoldCompile doMCClosureTests.C

echo ""
echo "done compiling. sleep."
echo ""
sleep 3s


echo ""
echo "running do MCClosureTests"
echo ""

./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_eta0to10      Py8_closureTest_eta0to10   1
./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_eta10to20     Py8_closureTest_eta10to20  1
./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_eta20to30     Py8_closureTest_eta20to30  1
./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_eta32to47     Py8_closureTest_eta32to47  1

echo ""
echo "done!"
echo ""
