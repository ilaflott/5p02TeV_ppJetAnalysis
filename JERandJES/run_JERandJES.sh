#!/bin/bash


#./deriveJERUnc.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_10-20-16__MCJEC_JetID0 Py8_QCD_deriveJERUnc
#./deriveJERUnc.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_10-20-16__MCJEC_JetID0 Py8_QCD_deriveJERUnc
#./deriveJERUnc.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak5PFJets_10-20-16__MCJEC_JetID0 Py8_QCD_deriveJERUnc
#./deriveJERUnc.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_10-20-16__MCJEC_JetID1 Py8_QCD_deriveJERUnc
#./deriveJERUnc.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_10-20-16__MCJEC_JetID1 Py8_QCD_deriveJERUnc
#./deriveJERUnc.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak5PFJets_10-20-16__MCJEC_JetID1 Py8_QCD_deriveJERUnc

echo "compiling..."
rootcompile plotNfit.C
echo ""
echo "done compiling. Running..."
echo ""
echo ""
#
#./plotNfit.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_10-20-16__MCJEC_JetID0 Py8_QCD_plotNfit
#./plotNfit.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_10-20-16__MCJEC_JetID0 Py8_QCD_plotNfit
#./plotNfit.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak5PFJets_10-20-16__MCJEC_JetID0 Py8_QCD_plotNfit
./plotNfit.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_10-20-16__MCJEC_JetID1 Py8_QCD_plotNfit
./plotNfit.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_10-20-16__MCJEC_JetID1 Py8_QCD_plotNfit
./plotNfit.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak5PFJets_10-20-16__MCJEC_JetID1 Py8_QCD_plotNfit
#
#
#
echo "compiling eta check, genpt 30 to 50..."
rootcompile plotNfit_checkEta30to50.C
echo ""
echo "done compiling. Running..."
echo ""
echo ""

#./plotNfit_checkEta30to50.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_10-20-16__MCJEC_JetID0 Py8_QCD_plotNfit_checkEta30to50
#./plotNfit_checkEta30to50.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_10-20-16__MCJEC_JetID0 Py8_QCD_plotNfit_checkEta30to50
#./plotNfit_checkEta30to50.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak5PFJets_10-20-16__MCJEC_JetID0 Py8_QCD_plotNfit_checkEta30to50
./plotNfit_checkEta30to50.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_10-20-16__MCJEC_JetID1 Py8_QCD_plotNfit_checkEta30to50
./plotNfit_checkEta30to50.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_10-20-16__MCJEC_JetID1 Py8_QCD_plotNfit_checkEta30to50
./plotNfit_checkEta30to50.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak5PFJets_10-20-16__MCJEC_JetID1 Py8_QCD_plotNfit_checkEta30to50



echo "compiling eta check, genpt 150 to 200..."
rootcompile plotNfit_checkEta150to200.C
echo ""
echo "done compiling. Running..."
echo ""
echo ""
#
#./plotNfit_checkEta150to200.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_10-20-16__MCJEC_JetID0 Py8_QCD_plotNfit_checkEta150to200
#./plotNfit_checkEta150to200.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_10-20-16__MCJEC_JetID0 Py8_QCD_plotNfit_checkEta150to200
#./plotNfit_checkEta150to200.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak5PFJets_10-20-16__MCJEC_JetID0 Py8_QCD_plotNfit_checkEta150to200
./plotNfit_checkEta150to200.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_10-20-16__MCJEC_JetID1 Py8_QCD_plotNfit_checkEta150to200
./plotNfit_checkEta150to200.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_10-20-16__MCJEC_JetID1 Py8_QCD_plotNfit_checkEta150to200
./plotNfit_checkEta150to200.exe ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak5PFJets_10-20-16__MCJEC_JetID1 Py8_QCD_plotNfit_checkEta150to200