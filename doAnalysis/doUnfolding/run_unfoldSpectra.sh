#!/bin/bash

rooUnfoldCompile unfoldDataSpectra.C

## unfold the data
./unfoldDataSpectra.exe    ppData_HighPtJetTrig_ak4PFJets_03-06-17_jetPlots_eta0to24   ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_03-23-17_JERS_eta0to24 HPtJetTrig_Py8_unfoldedData 1

## closure 
#./doMCClosureTests.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_   pp_MCPy8_closureTests 0
#./doMCClosureTests.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_   pp_MCPy8_closureTests 1
