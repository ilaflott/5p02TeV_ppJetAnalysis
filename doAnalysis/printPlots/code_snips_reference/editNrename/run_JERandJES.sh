#!/bin/bash

echo "compiling..."
rootcompile plotNfit.C
echo ""
echo "done compiling. Running..."
echo ""
echo ""

./plotNfit.exe readForests/outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_11-22-16__tight_MCJEC/ Py8_QCD_plotNfit_ak3jtIDtgt

#
#
#
echo "compiling eta check, genpt 30 to 50..."
rootcompile plotNfit_checkEta30to50.C
echo ""
echo "done compiling. Running..."
echo ""
echo ""

./plotNfit_checkEta30to50.exe readForests/outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_11-22-16__tight_MCJEC/ Py8_QCD_plotNfit_ak3_jtIDtgt_eta_30pt50


echo "compiling eta check, genpt 150 to 200..."
rootcompile plotNfit_checkEta150to200.C
echo ""
echo "done compiling. Running..."
echo ""
echo ""

./plotNfit_checkEta150to200.exe readForests/outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_11-22-16__tight_MCJEC/ Py8_QCD_plotNfit_ak3_jtIDtgt_eta_150pt200