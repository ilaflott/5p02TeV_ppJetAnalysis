#!/bin/bash

#./printPlots_jetPlots_wRapBins.exe  <ppData dir> <ppMC dir <overlay flag>
#./printPlots_jetPlots_combine.exe   <ppData dir> <ppMC dir> <outputNameBase>

#################################

# in /export/d00/scratch/ilaflott/5p02TeV_ppJetAnalysis/readForests/

# DATA
# 2.20.17_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-20-17_jetPlots_eta0to24/
# 2.20.17_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-20-17_jetPlots_eta24to27/
# 2.20.17_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-20-17_jetPlots_eta27to30/
# 2.20.17_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-20-17_jetPlots_eta32to47/

# MC
# 2.20.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-20-17_jetPlots_eta0to24/
# 2.20.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-20-17_jetPlots_eta24to27/
# 2.20.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-20-17_jetPlots_eta27to30/
# 2.20.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-20-17_jetPlots_eta32to47/


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
echo ""
echo ""
echo "Running eta0to24..."
./printPlots_jetPlots_combine.exe  2.20.17_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-20-17_jetPlots_eta0to24/ 2.20.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-20-17_jetPlots_eta0to24/ eta0to24
sleep 1s







echo ""
echo ""
echo ""
echo ""
echo "Running eta24to27..."
./printPlots_jetPlots_combine.exe  2.20.17_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-20-17_jetPlots_eta24to27/ 2.20.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-20-17_jetPlots_eta24to27/ eta24to27
sleep 1s




echo ""
echo ""
echo ""
echo ""
echo "Running eta27to30..."
./printPlots_jetPlots_combine.exe  2.20.17_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-20-17_jetPlots_eta27to30/ 2.20.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-20-17_jetPlots_eta27to30/ eta27to30
sleep 1s



echo ""
echo ""
echo ""
echo ""
echo "Running eta32to47..."
./printPlots_jetPlots_combine.exe  2.20.17_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-20-17_jetPlots_eta32to47/ 2.20.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-20-17_jetPlots_eta32to47/ eta32to47
sleep 1s