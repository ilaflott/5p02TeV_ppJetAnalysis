#!/bin/bash


./unfoldDataSpectra.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_10-18-16__Unfold   ppData_HighPtJetTrig_ak3PFJets_10-18-16   HPtJetTrig_Py8_unfoldedData 0
./unfoldDataSpectra.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_10-18-16__Unfold   ppData_HighPtJetTrig_ak4PFJets_10-18-16   HPtJetTrig_Py8_unfoldedData 0
./unfoldDataSpectra.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak5PFJets_10-18-16__Unfold   ppData_HighPtJetTrig_ak5PFJets_10-18-16   HPtJetTrig_Py8_unfoldedData 0
./unfoldDataSpectra.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_10-18-16__Unfold   ppData_HighPtJetTrig_ak3PFJets_10-18-16   HPtJetTrig_Py8_unfoldedData 1
./unfoldDataSpectra.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_10-18-16__Unfold   ppData_HighPtJetTrig_ak4PFJets_10-18-16   HPtJetTrig_Py8_unfoldedData 1
./unfoldDataSpectra.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak5PFJets_10-18-16__Unfold   ppData_HighPtJetTrig_ak5PFJets_10-18-16   HPtJetTrig_Py8_unfoldedData 1


./doMCClosureTests.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_10-18-16__Unfold   pp_MCPy8_closureTests 0
./doMCClosureTests.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_10-18-16__Unfold   pp_MCPy8_closureTests 0
./doMCClosureTests.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak5PFJets_10-18-16__Unfold   pp_MCPy8_closureTests 0
./doMCClosureTests.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_10-18-16__Unfold   pp_MCPy8_closureTests 1
./doMCClosureTests.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_10-18-16__Unfold   pp_MCPy8_closureTests 1
./doMCClosureTests.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak5PFJets_10-18-16__Unfold   pp_MCPy8_closureTests 1
