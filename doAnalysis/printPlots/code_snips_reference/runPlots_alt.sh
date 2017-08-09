#!/bin/bash

#//int printQAPlots(const std::string input_ppData_condorDir, const std::string input_ppData_Filename,
#//                 const std::string input_ppMC_condorDir, const std::string input_ppMC_Filename,
#//                 const std::string output_PDFname_base_full, const bool drawDataMCOverlaysInput=true )

theCode=printQAPlots_altLumiEff

echo "compiling plot script"
echo ""
rootcompile ${theCode}.C


echo "running plots for ak3PFJets_HLT.ak4CaloJets"
echo ""
./${theCode}.exe ppData_HighPtJetTrig_ak3PFJets_HLT.ak4CaloJets_09-19-16__wFull13TeVJetID/ HighPtJet80_ak3PF-allFiles.root ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_09-19-16__wFull13TeVJetID/ Py8_CUETP8M1_QCDjetAllPtBins_ak3PF-allFiles.root printQAplots_alt_9.21.16_13TeVJetID_ak3PFJets_HLT.ak4CaloJets_QAPlots_ 0
sleep 2.0s

./${theCode}.exe ppData_HighPtJetTrig_ak3PFJets_HLT.ak4CaloJets_09-19-16__wFull13TeVJetID/ HighPtJet80_ak3PF-allFiles.root ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_09-19-16__wFull13TeVJetID/ Py8_CUETP8M1_QCDjetAllPtBins_ak3PF-allFiles.root printQAplots_alt_9.21.16_13TeVJetID_ak3PFJets_HLT.ak4CaloJets_QAPlots_ 1
sleep 2.0s

echo "running plots for ak4PFJets_HLT.ak4CaloJets"
echo ""
./${theCode}.exe ppData_HighPtJetTrig_ak4PFJets_HLT.ak4CaloJets_09-19-16__wFull13TeVJetID/ HighPtJet80_ak4PF-allFiles.root ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_09-19-16__wFull13TeVJetID/ Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root printQAplots_alt_9.21.16_13TeVJetID_ak4PFJets_HLT.ak4CaloJets_QAPlots_ 0
sleep 2.0s

./${theCode}.exe ppData_HighPtJetTrig_ak4PFJets_HLT.ak4CaloJets_09-19-16__wFull13TeVJetID/ HighPtJet80_ak4PF-allFiles.root ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_09-19-16__wFull13TeVJetID/ Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root printQAplots_alt_9.21.16_13TeVJetID_ak4PFJets_HLT.ak4CaloJets_QAPlots_ 1
sleep 2.0s


echo "running plots for ak4CaloJets_HLT.ak4CaloJets"
echo ""
./${theCode}.exe ppData_HighPtJetTrig_ak4CaloJets_HLT.ak4CaloJets_09-19-16__wFull13TeVJetID/ HighPtJet80_ak4Calo-allFiles.root ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4CaloJets_09-19-16__wFull13TeVJetID/ Py8_CUETP8M1_QCDjetAllPtBins_ak4Calo-allFiles.root printQAplots_alt_9.21.16_13TeVJetID_ak4CaloJets_HLT.ak4CaloJets_QAPlots_ 0
sleep 2.0s

./${theCode}.exe ppData_HighPtJetTrig_ak4CaloJets_HLT.ak4CaloJets_09-19-16__wFull13TeVJetID/ HighPtJet80_ak4Calo-allFiles.root ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4CaloJets_09-19-16__wFull13TeVJetID/ Py8_CUETP8M1_QCDjetAllPtBins_ak4Calo-allFiles.root printQAplots_alt_9.21.16_13TeVJetID_ak4CaloJets_HLT.ak4CaloJets_QAPlots_ 1
sleep 2.0s

return