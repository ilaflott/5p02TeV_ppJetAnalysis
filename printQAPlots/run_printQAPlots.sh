#!/bin/bash

echo ""
echo ""
echo ""
echo "compiling..."
echo ""
echo ""
echo ""
echo ""
rootcompile printQAPlots.C

echo ""
echo "R=3..."
echo ""
./printQAPlots.exe 11.17.16_outputCondor/ppData_HighPtJetTrig_ak3PFJets_11-17-16     11.17.16_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_11-17-16 jetPlots_wResid_jetID1 1

./printQAPlots.exe 11.17.16_outputCondor/ppData_HighPtJetTrig_ak3PFJets_11-17-16     11.17.16_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_11-17-16 jetPlots_wResid_jetID1 0


echo ""
echo "R=4..."
echo ""
./printQAPlots.exe 11.17.16_outputCondor/ppData_HighPtJetTrig_ak4PFJets_11-17-16     11.17.16_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_11-17-16 jetPlots_wResid_jetID1 1

./printQAPlots.exe 11.17.16_outputCondor/ppData_HighPtJetTrig_ak4PFJets_11-17-16     11.17.16_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_11-17-16 jetPlots_wResid_jetID1 0




