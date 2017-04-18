#!/bin/bash

#./printPlots_jetPlots_wRapBins.exe  <ppData dir> <ppMC dir <overlay flag>
#./printPlots_jetPlots_combine.exe   <ppData dir> <ppMC dir> <outputNameBase>

#################################

# DATA
# 4.17.17_outputCondor/ppData_HighPtJetTrig_ak4PFJets_04-18-17_jetPlots_eta0to10/
# 4.17.17_outputCondor/ppData_HighPtJetTrig_ak4PFJets_04-18-17_jetPlots_eta10to20/
# 4.17.17_outputCondor/ppData_HighPtJetTrig_ak4PFJets_04-18-17_jetPlots_eta20to30/
# 4.17.17_outputCondor/ppData_HighPtJetTrig_ak4PFJets_04-18-17_jetPlots_eta32to47/

# MC
# 4.17.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_04-18-17_jetPlots_eta0to10/
# 4.17.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_04-18-17_jetPlots_eta10to20/
# 4.17.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_04-18-17_jetPlots_eta20to30/
# 4.17.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_04-18-17_jetPlots_eta32to47/
# 4.17.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_Official_ak4PFJets_04-18-17_jetPlots_eta0to10/
# 4.17.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_Official_ak4PFJets_04-18-17_jetPlots_eta10to20/
# 4.17.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_Official_ak4PFJets_04-18-17_jetPlots_eta20to30/
# 4.17.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_Official_ak4PFJets_04-18-17_jetPlots_eta32to47/


echo ""
echo ""
echo ""
echo ""
echo "compiling printPlots_jetPlots"
rootcompile printPlots_jetPlots_combine.C

echo ""
echo ""
echo ""
echo ""
echo "compilation done..."
sleep 2s


echo ""
echo ""
echo "Running eta0to10..."
echo ""
echo ""
./printPlots_jetPlots_combine.exe   ppData_HighPtJetTrig_ak4PFJets_04-18-17_jetPlots_eta0to10/ ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_04-18-17_jetPlots_eta0to10/   eta0to10_semiOffPy8_ppMC
sleep 1s

echo ""
echo ""
echo "Running eta10to20..."
echo ""
echo ""
./printPlots_jetPlots_combine.exe   ppData_HighPtJetTrig_ak4PFJets_04-18-17_jetPlots_eta10to20/ ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_04-18-17_jetPlots_eta10to20/   eta10to20_semiOffPy8_ppMC
sleep 1s

echo ""
echo ""
echo "Running eta20to30..."
echo ""
echo ""
./printPlots_jetPlots_combine.exe   ppData_HighPtJetTrig_ak4PFJets_04-18-17_jetPlots_eta20to30/ ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_04-18-17_jetPlots_eta20to30/   eta20to30_semiOffPy8_ppMC
sleep 1s

echo ""
echo ""
echo "Running eta32to47..."
echo ""
echo ""
./printPlots_jetPlots_combine.exe   ppData_HighPtJetTrig_ak4PFJets_04-18-17_jetPlots_eta32to47/ ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_04-18-17_jetPlots_eta32to47/   eta32to47_semiOffPy8_ppMC
sleep 1s