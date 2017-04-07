#!/bin/bash


#[ppMC_pthatWeights $] ./ppMC_pthatWeights.exe HELP
#to use evtPthatWeights, you need a filelist. Then do...
#./evtPthatWeights.exe <inputFileList> <outputWeightFile> <radius> <jetType>
#or leave arguments blank for default tests
#rStatus=-1

echo ""
echo "running ppMC_pthatWeights in individual bins for later combination..."
echo ""
echo ""

echo ""
echo "running pthat370..."
echo ""
nice -n 19  ./ppMC_pthatWeights.exe 5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet370_Official_forests.txt   QCDjet370_pthatWeight.out 4 PF 370

echo ""
echo "done w/ pthat370. Sleeping to be kind."
echo ""
sleep 1 min

echo ""
echo "running pthat280..."
echo ""
nice -n 19  ./ppMC_pthatWeights.exe 5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet280_Official_forests.txt   QCDjet280_pthatWeight.out 4 PF 280

echo ""
echo "done w/ pthat280. Sleeping to be kind."
echo ""
sleep 1 min

echo ""
echo "running pthat220..."
echo ""
nice -n 19  ./ppMC_pthatWeights.exe 5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet220_Official_forests.txt   QCDjet220_pthatWeight.out 4 PF 220

echo ""
echo "done w/ pthat220. Sleeping to be kind."
echo ""
sleep 1 min

echo ""
echo "running pthat170..."
echo ""
nice -n 19  ./ppMC_pthatWeights.exe 5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet170_Official_forests.txt   QCDjet170_pthatWeight.out 4 PF 170

echo ""
echo "done w/ pthat170. Sleeping to be kind."
echo ""
sleep 1 min

echo ""
echo "running pthat120..."
echo ""
nice -n 19  ./ppMC_pthatWeights.exe 5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet120_Official_forests.txt   QCDjet120_pthatWeight.out 4 PF 120

echo ""
echo "done w/ pthat120. Sleeping to be kind."
echo ""
sleep 1 min

echo ""
echo "running pthat80..."
echo ""
nice -n 19  ./ppMC_pthatWeights.exe 5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet80_Official_forests.txt    QCDjet80_pthatWeight.out 4 PF 80
echo ""
echo "done w/ pthat80. Sleeping to be kind."
echo ""
sleep 1 min

echo ""
echo "running pthat50..."
echo ""
nice -n 19  ./ppMC_pthatWeights.exe 5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet50_Official_forests.txt    QCDjet50_pthatWeight.out 4 PF 50

echo ""
echo "done w/ pthat50. Sleeping to be kind."
echo ""
sleep 1 min

echo ""
echo "running pthat30..."
echo ""
nice -n 19  ./ppMC_pthatWeights.exe 5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet30_Official_forests.txt    QCDjet30_pthatWeight.out 4 PF 30

echo ""
echo "done w/ pthat30. Sleeping to be kind."
echo ""
sleep 1 min

echo ""
echo "running pthat15 (LAST BIN)..."
echo ""
nice -n 19  ./ppMC_pthatWeights.exe 5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet15_Official_forests.txt    QCDjet15_pthatWeight.out 4 PF 15

echo ""
echo "no sleeping... done!"
echo ""

exit 1









