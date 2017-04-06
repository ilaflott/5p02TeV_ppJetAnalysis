#!/bin/bash

echo ""
echo "making filelists for ppMC."
echo ""


### pp Reco Dijet MC Forests
## pp reconstruction, jec v12, unmerged files
### pthat bins {15,30,50,80,120,170,220,280,370,460,540}
#ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet15_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet15_forests.txt
#ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet30_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet30_forests.txt
#ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet50_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet50_forests.txt
#ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet80_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet80_forests.txt
#ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet120_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet120_forests.txt
#ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet170_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet170_forests.txt
#ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet220_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet220_forests.txt
#ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet280_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet280_forests.txt
#ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet370_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet370_forests.txt
#ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet460_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet460_forests.txt
#ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet540_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet540_forests.txt

## combined list, pthat bins {15,30,50,80,120,170,220,280,370,460,540}, semi private
#ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet*_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt

### merged forests, same events at above files, HcalRespCorr_v4
#ls /mnt/hadoop/cms/store/user/abaty/mergedForests/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV_ppSignal/Pythia8_Dijet*_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/0.root >> 5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_merged_forests.txt


### NOTES
# pp dijetMC w/ HI reco available, see 2015 HiForest twiki if interested
# pp dijetMC+MinBias w/ pp reco available for pileup study, see 2015 HiForest twiki if interested


### pp Reco Dijet MC Forests
### Official sample, unmerged files
### pthat bins {15,30,50,80,120,170,220,280,370
ls /mnt/hadoop/cms/store/user/dgulhan/Pythia8_Dijet15_pp502_TuneCUETP8M1/PYTHIA_QCD15_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_genpt0GeV/*/000?/HiForestAOD_*.root >>  5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet15_Official_forests.txt
ls /mnt/hadoop/cms/store/user/dgulhan/Pythia8_Dijet30_pp502_TuneCUETP8M1/PYTHIA_QCD30_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_genpt0GeV/*/000?/HiForestAOD_*.root >>  5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet30_Official_forests.txt
ls /mnt/hadoop/cms/store/user/dgulhan/Pythia8_Dijet50_pp502_TuneCUETP8M1/PYTHIA_QCD50_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_genpt0GeV/*/000?/HiForestAOD_*.root >>  5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet50_Official_forests.txt
ls /mnt/hadoop/cms/store/user/dgulhan/Pythia8_Dijet80_pp502_TuneCUETP8M1/PYTHIA_QCD80_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_genpt0GeV/*/000?/HiForestAOD_*.root >>  5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet80_Official_forests.txt
ls /mnt/hadoop/cms/store/user/dgulhan/Pythia8_Dijet120_pp502_TuneCUETP8M1/PYTHIA_QCD120_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_genpt0GeV/*/000?/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet120_Official_forests.txt
ls /mnt/hadoop/cms/store/user/dgulhan/Pythia8_Dijet170_pp502_TuneCUETP8M1/PYTHIA_QCD170_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_genpt0GeV/*/000?/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet170_Official_forests.txt
ls /mnt/hadoop/cms/store/user/dgulhan/Pythia8_Dijet220_pp502_TuneCUETP8M1/PYTHIA_QCD220_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_genpt0GeV/*/000?/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet220_Official_forests.txt
ls /mnt/hadoop/cms/store/user/dgulhan/Pythia8_Dijet280_pp502_TuneCUETP8M1/PYTHIA_QCD280_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_genpt0GeV/*/000?/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet280_Official_forests.txt
ls /mnt/hadoop/cms/store/user/dgulhan/Pythia8_Dijet370_pp502_TuneCUETP8M1/PYTHIA_QCD370_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_genpt0GeV/*/000?/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/5p02TeV_Py8_CUETP8M1_QCDjet370_Official_forests.txt

### combined list, pthat bins {15,30,50,80,120,170,220,280,370,460,540}, official
#ls /mnt/hadoop/cms/store/user/dgulhan/Pythia8_Dijet*_pp502_TuneCUETP8M1/PYTHIA_QCD*_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_genpt0GeV/*/000?/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt


echo ""
echo "done. "
echo ""

return