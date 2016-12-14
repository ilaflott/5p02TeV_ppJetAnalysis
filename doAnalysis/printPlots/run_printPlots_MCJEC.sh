#!/bin/bash


echo "compiling printPlots_MCJEC"
#rootcompile printPlots_MCJEC.C
echo ""
echo "done compiling. Running..."
echo ""
echo ""

./printPlots_MCJEC.exe outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_11-28-16_MCJEC_loose testingPrintPlots 



echo "compiling eta check, genpt 30 to 50..."
#rootcompile printPlots_MCJEC_genpt30to50.C
echo ""
echo "done compiling. Running..."
echo ""
echo ""

./printPlots_MCJEC_genpt30to50.exe outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_11-28-16_MCJEC_loose testingPrintPlots 


echo "compiling eta check, genpt 150 to 200..."
#rootcompile printPlots_MCJEC_genpt150to200.C
echo ""
echo "done compiling. Running..."
echo ""
echo ""

./printPlots_MCJEC_genpt150to200.exe outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_11-28-16_MCJEC_loose testingPrintPlots 



echo "compiling printPlots_MCJEC_MCEff"
#rootcompile printPlots_MCJEC_MCEff.C
echo ""
echo "done compiling. Running..."
echo ""
echo ""

./printPlots_MCJEC_MCEff.exe outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_11-28-16_MCJEC_loose testingPrintPlots 