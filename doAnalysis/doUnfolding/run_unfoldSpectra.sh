#!/bin/bash

### unfoldDataSpectra
echo ""
echo "compiling unfoldDataSpectra"
echo ""
rooUnfoldCompile unfoldDataSpectra.C

echo ""
echo "done compiling. sleep."
echo ""
sleep 0.5s

echo ""
echo "running unfoldDataSpectra"
echo ""

#./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak4PFJets_05-04-17_jetPlots_1.5eta2.0  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_1.5eta2.0     HPtJetTrig_Py8_unfoldedData_15eta20  0
#./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak4PFJets_05-04-17_jetPlots_1.0eta1.5  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_1.0eta1.5     HPtJetTrig_Py8_unfoldedData_10eta15  0
#./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak4PFJets_05-04-17_jetPlots_0.5eta1.0  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_0.5eta1.0     HPtJetTrig_Py8_unfoldedData_05eta10  0
./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak4PFJets_05-04-17_jetPlots_0.0eta0.5  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_0.0eta0.5     HPtJetTrig_Py8_unfoldedData_00eta05  0
#
#./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak4PFJets_05-04-17_jetPlots_1.5eta2.0  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_1.5eta2.0     HPtJetTrig_Py8_unfoldedData_15eta20  1
#./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak4PFJets_05-04-17_jetPlots_1.0eta1.5  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_1.0eta1.5     HPtJetTrig_Py8_unfoldedData_10eta15  1
#./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak4PFJets_05-04-17_jetPlots_0.5eta1.0  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_0.5eta1.0     HPtJetTrig_Py8_unfoldedData_05eta10  1
#./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak4PFJets_05-04-17_jetPlots_0.0eta0.5  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_0.0eta0.5     HPtJetTrig_Py8_unfoldedData_00eta05  1

echo ""
echo "done running unfoldDataSpectra. sleep."
echo ""
#sleep 1s


#### doMCClosureTests
#echo ""
#echo "compiling doMCClosureTests"
#echo ""
#rooUnfoldCompile doMCClosureTests.C
#
#echo ""
#echo "done compiling. sleep."
#echo ""
#sleep 3s
#
#
#echo ""
#echo "running do MCClosureTests"
#echo ""
#
#./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_1.5eta2.0     Py8_closureTest_15eta20  1
##./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_1.0eta1.5     Py8_closureTest_10eta15  1
##./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_0.5eta1.0     Py8_closureTest_05eta10  1
##./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_05-04-17_JERS_0.0eta0.5     Py8_closureTest_00eta05  1
##

echo ""
echo "done!"
echo ""

return