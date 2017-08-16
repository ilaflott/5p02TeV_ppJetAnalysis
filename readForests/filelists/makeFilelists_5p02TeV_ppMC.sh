#!/bin/bash

echo ""
echo "making filelists for ppMC."
echo ""





### NOTES
# pp dijetMC w/ HI reco available, see 2015 HiForest twiki if interested
# pp dijetMC+MinBias w/ pp reco available for pileup study, see 2015 HiForest twiki if interested







touch 5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet15_RUT3_forests.txt 
touch 5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet30_RUT3_forests.txt 
touch 5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet50_RUT3_forests.txt 
touch 5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet80_RUT3_forests.txt 
touch 5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet120_RUT3_forests.txt
touch 5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet170_RUT3_forests.txt
touch 5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet220_RUT3_forests.txt
touch 5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet280_RUT3_forests.txt
touch 5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet370_RUT3_forests.txt
touch 5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet460_RUT3_forests.txt
touch 5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet540_RUT3_forests.txt

touch 5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_RUT3_forests.txt

touch 5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_RUT3_merged_forests.txt







### pp Reco Dijet MC Forests
### SemiOfficial sample, unmerged files, RUT3
### pthat bins {15,30,50,80,120,170,220,280,370
ls /cms/se/store/user/ilaflott/jobOutput2/Pythia8_Dijet_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/dijet15/HiForestAOD_*.root >>   5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet15_RUT3_forests.txt 
ls /cms/se/store/user/ilaflott/jobOutput2/Pythia8_Dijet_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/dijet30/HiForestAOD_*.root >>   5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet30_RUT3_forests.txt 
ls /cms/se/store/user/ilaflott/jobOutput2/Pythia8_Dijet_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/dijet50/HiForestAOD_*.root >>   5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet50_RUT3_forests.txt 
ls /cms/se/store/user/ilaflott/jobOutput2/Pythia8_Dijet_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/dijet80/HiForestAOD_*.root >>   5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet80_RUT3_forests.txt 
ls /cms/se/store/user/ilaflott/jobOutput2/Pythia8_Dijet_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/dijet120/HiForestAOD_*.root >>  5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet120_RUT3_forests.txt
ls /cms/se/store/user/ilaflott/jobOutput2/Pythia8_Dijet_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/dijet170/HiForestAOD_*.root >>  5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet170_RUT3_forests.txt
ls /cms/se/store/user/ilaflott/jobOutput2/Pythia8_Dijet_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/dijet220/HiForestAOD_*.root >>  5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet220_RUT3_forests.txt
ls /cms/se/store/user/ilaflott/jobOutput2/Pythia8_Dijet_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/dijet280/HiForestAOD_*.root >>  5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet280_RUT3_forests.txt
ls /cms/se/store/user/ilaflott/jobOutput2/Pythia8_Dijet_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/dijet370/HiForestAOD_*.root >>  5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet370_RUT3_forests.txt
ls /cms/se/store/user/ilaflott/jobOutput2/Pythia8_Dijet_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/dijet460/HiForestAOD_*.root >>  5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet460_RUT3_forests.txt
ls /cms/se/store/user/ilaflott/jobOutput2/Pythia8_Dijet_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/dijet540/HiForestAOD_*.root >>  5p02TeV_Py8_CUETP8M1_QCDjet_RUT3_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet540_RUT3_forests.txt

## combined list, pthat bins {15,30,50,80,120,170,220,280,370,460,540}, semiofficial
ls /cms/se/store/user/ilaflott/jobOutput2/Pythia8_Dijet_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/dijet*/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_RUT3_forests.txt

## merged forests, same events at above files, HcalRespCorr_v4
ls /cms/se/store/user/ilaflott/jobOutput2/Pythia8_Dijet_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/merged/*.root >> 5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_RUT3_merged_forests.txt


echo ""
echo "done. "
echo ""

return