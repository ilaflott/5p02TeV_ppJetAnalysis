#!/bin/bash

#### unfoldDataSpectra
#rooUnfoldCompile unfoldDataSpectra.C
#
#./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak4PFJets_04-18-17_jetPlots_eta0to10   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_04-17-17_JERS_eta0to10      HPtJetTrig_Py8_unfoldedData_eta0to10   0
#./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak4PFJets_04-18-17_jetPlots_eta10to20  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_04-17-17_JERS_eta10to20     HPtJetTrig_Py8_unfoldedData_eta10to20  0
#./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak4PFJets_04-18-17_jetPlots_eta20to30  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_04-17-17_JERS_eta20to30     HPtJetTrig_Py8_unfoldedData_eta20to30  0
#./unfoldDataSpectra.exe ppData_HighPtJetTrig_ak4PFJets_04-18-17_jetPlots_eta32to47  ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_04-17-17_JERS_eta32to47     HPtJetTrig_Py8_unfoldedData_eta32to47  0


### closure 
echo ""
echo "compiling doMCClosureTests"
echo ""
rooUnfoldCompile doMCClosureTests.C

echo ""
echo "running do MCClosureTests"
echo ""
./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_04-17-17_JERS_eta0to10      Py8_closureTest_eta0to10   0
#./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_04-17-17_JERS_eta10to20     Py8_closureTest_eta10to20  0
#./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_04-17-17_JERS_eta20to30     Py8_closureTest_eta20to30  0
#./doMCClosureTests.exe   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_04-17-17_JERS_eta32to47     Py8_closureTest_eta32to47  0
##./doMCClosureTests.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_   pp_MCPy8_closureTests 0
##./doMCClosureTests.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_   pp_MCPy8_closureTests 1

echo ""
echo "done!"
echo ""
