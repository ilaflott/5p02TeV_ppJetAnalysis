#!/bin/bash

# unmerged files
# comment: pp reconstruction, jec v12, PrivMC
ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet15_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet15_forests.txt
ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet30_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet30_forests.txt
ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet50_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet50_forests.txt
ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet80_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet80_forests.txt
ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet120_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet120_forests.txt
ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet170_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet170_forests.txt
ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet220_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet220_forests.txt
ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet280_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet280_forests.txt
ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet370_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet370_forests.txt
ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet460_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet460_forests.txt
ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet540_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjet540_forests.txt
# this list is a conglomerate of all above lists
ls /mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet*_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/*/0000/HiForestAOD_*.root >> 5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt
##

# same events as QCDJet???/QCDjetAllPtBins as above, but all file of specific pthat merged into one file
# comment: HcalRespCorr_v4
ls /mnt/hadoop/cms/store/user/abaty/mergedForests/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV_ppSignal/Pythia8_Dijet*_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/0.root >> 5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_merged_forests.txt
##

# Official MinBias, HINppWinter16DR-NoPU_75X_mcRun2_asymptotic_ppAt5TeV_forest_v2
ls /mnt/hadoop/cms/store/user/rbi/merged/MinBias_TuneCUETP8M1_5p02TeV-pythia8-HINppWinter16DR-NoPU_75X_mcRun2_asymptotic_ppAt5TeV_forest_v2/0.root >> 5p02TeV_Py8_CUETP8M1_MinBias_HINppWinter16_merged_forests.txt
##