#!/bin/bash

rootcompile printPlots_jetPlots.C
./printPlots_jetPlots.exe   11.22.16_outputCondor/ppData_HighPtJetTrig_ak3PFJets_11-22-16__loose  11.22.16_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_11-22-16__loose_jetPlots TEST_jtIDL_ak3PF_jetPlots_11.28.16 0
./printPlots_jetPlots.exe   11.22.16_outputCondor/ppData_HighPtJetTrig_ak3PFJets_11-22-16__loose  11.22.16_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_11-22-16__loose_jetPlots TEST_jtIDL_ak3PF_jetPlots_11.28.16 1

#./printPlots_jetPlots.exe   11.22.16_outputCondor/ppData_HighPtJetTrig_ak3PFJets_11-22-16__loose  11.22.16_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_11-22-16__loose_jetPlots jtIDL_ak3PF_jetPlots_11.28.16 0
#./printPlots_jetPlots.exe   11.22.16_outputCondor/ppData_HighPtJetTrig_ak3PFJets_11-22-16__tight  11.22.16_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_11-22-16__tight_jetPlots jtIDT_ak3PF_jetPlots_11.28.16 0
#./printPlots_jetPlots.exe   11.22.16_outputCondor/ppData_HighPtJetTrig_ak4PFJets_11-22-16__loose  11.22.16_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_11-22-16__loose_jetPlots jtIDL_ak4PF_jetPlots_11.28.16 0
#./printPlots_jetPlots.exe   11.22.16_outputCondor/ppData_HighPtJetTrig_ak4PFJets_11-22-16__tight  11.22.16_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_11-22-16__tight_jetPlots jtIDT_ak4PF_jetPlots_11.28.16 0
#											                                                                                                   
#./printPlots_jetPlots.exe   11.22.16_outputCondor/ppData_HighPtJetTrig_ak3PFJets_11-22-16__loose  11.22.16_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_11-22-16__loose_jetPlots jtIDL_ak3PF_jetPlots_11.28.16 1
#./printPlots_jetPlots.exe   11.22.16_outputCondor/ppData_HighPtJetTrig_ak3PFJets_11-22-16__tight  11.22.16_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_11-22-16__tight_jetPlots jtIDT_ak3PF_jetPlots_11.28.16 1
#./printPlots_jetPlots.exe   11.22.16_outputCondor/ppData_HighPtJetTrig_ak4PFJets_11-22-16__loose  11.22.16_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_11-22-16__loose_jetPlots jtIDL_ak4PF_jetPlots_11.28.16 1
#./printPlots_jetPlots.exe   11.22.16_outputCondor/ppData_HighPtJetTrig_ak4PFJets_11-22-16__tight  11.22.16_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_11-22-16__tight_jetPlots jtIDT_ak4PF_jetPlots_11.28.16 1
