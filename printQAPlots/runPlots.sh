#!/bin/bash

#//int printQAPlots(const std::string input_ppData_condorDir, const std::string input_ppData_Filename,
#//                 const std::string input_ppMC_condorDir, const std::string input_ppMC_Filename,
#//                 const std::string output_PDFname_base_full, const bool drawDataMCOverlaysInput=true )

echo "compiling plot script"
echo ""
rootcompile printQAPlots.C


echo "running plots for ak3PFJets_HLT.ak4CaloJets"
echo ""
./printQAPlots.exe ppData_HighPtJetTrig_ak3PFJets_HLT.ak4CaloJets_09-19-16__wFull13TeVJetID/ HighPtJetTrig_ak3PF-allFiles.root ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_09-19-16__wFull13TeVJetID/ Py8_CUETP8M1_QCDjetAllPtBins_ak3PF-allFiles.root printQAplots_9.21.16_13TeVJetID_ak3PFCaloJets_HLT.ak4CaloJets_QAPlots_ 0
return
sleep 5.0s
./printQAPlots.exe ppData_HighPtJetTrig_ak3PFJets_HLT.ak4CaloJets_09-19-16__wFull13TeVJetID/ HighPtJetTrig_ak3PF-allFiles.root ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_09-19-16__wFull13TeVJetID/ Py8_CUETP8M1_QCDjetAllPtBins_ak3PF-allFiles.root printQAplots_9.21.16_13TeVJetID_ak3PFCaloJets_HLT.ak4CaloJets_QAPlots_ 1
sleep 5.0s


echo "running plots for ak4PFJets_HLT.ak4CaloJets"
echo ""
./printQAPlots.exe ppData_HighPtJetTrig_ak4PFJets_HLT.ak4CaloJets_09-19-16__wFull13TeVJetID/ HighPtJetTrig_ak4PF-allFiles.root ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_09-19-16__wFull13TeVJetID/ Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root printQAplots_9.21.16_13TeVJetID_ak4PFCaloJets_HLT.ak4CaloJets_QAPlots_ 0
sleep 5.0s
./printQAPlots.exe ppData_HighPtJetTrig_ak4PFJets_HLT.ak4CaloJets_09-19-16__wFull13TeVJetID/ HighPtJetTrig_ak4PF-allFiles.root ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_09-19-16__wFull13TeVJetID/ Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root printQAplots_9.21.16_13TeVJetID_ak4PFCaloJets_HLT.ak4CaloJets_QAPlots_ 1
sleep 5.0s


echo "running plots for ak4CaloJets_HLT.ak4CaloJets"
echo ""
./printQAPlots.exe ppData_HighPtJetTrig_ak4CaloJets_HLT.ak4CaloJets_09-19-16__wFull13TeVJetID/ HighPtJetTrig_ak4Calo-allFiles.root ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4CaloJets_09-19-16__wFull13TeVJetID/ Py8_CUETP8M1_QCDjetAllPtBins_ak4Calo-allFiles.root printQAplots_9.21.16_13TeVJetID_ak4PFCaloJets_HLT.ak4CaloJets_QAPlots_ 0
sleep 5.0s
./printQAPlots.exe ppData_HighPtJetTrig_ak4CaloJets_HLT.ak4CaloJets_09-19-16__wFull13TeVJetID/ HighPtJetTrig_ak4Calo-allFiles.root ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4CaloJets_09-19-16__wFull13TeVJetID/ Py8_CUETP8M1_QCDjetAllPtBins_ak4Calo-allFiles.root printQAplots_9.21.16_13TeVJetID_ak4PFCaloJets_HLT.ak4CaloJets_QAPlots_ 1
sleep 5.0s

return