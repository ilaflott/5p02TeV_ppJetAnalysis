#!/bin/bash

## pp Prompt Reco Forest (PD: /HighPt{Jet80,LowerJets}/Run2015E-PromptReco-v1/AOD)
# 262081 - 262328
## HighPtLowerJets, 99% 25.8/pb Golden JSON, v12 JECs
ls /mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/2015pp_HighPtLowerJets/*.root >> 5p02TeV_HighPtLowerJets_forests.txt 
## Jet80, (full) Golden JSON
ls /mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/2015pp_HighPtJet80/*.root >> 5p02TeV_HighPtJet80_forests.txt
ls /mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/2015pp_HighPtJet80/HiForestAOD_?.root >> 5p02TeV_HighPtJet80_9Files_debug_forests.txt
## combined Jet80, LowerJets 
ls /mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/2015pp_HighPtLowerJets/*.root >> 5p02TeV_HighPtJetTrig_forests.txt
ls /mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/2015pp_HighPtJet80/*.root     >> 5p02TeV_HighPtJetTrig_forests.txt
##

## pp Prompt Reco Forest (PD: /MinimumBias{1,2,6,10}/Run2015E-PromptReco-v1/AOD)
# 261395-262273
## MinimumBias{1,2,6,10}
ls   /mnt/hadoop/cms/store/user/dgulhan/MinimumBias1/HiForest_MinimumBias1_Run2015E_PromptReco_pp_5020GeV/151201_113508/0000/*.root >>  5p02TeV_MinBias1_forests.txt
ls   /mnt/hadoop/cms/store/user/dgulhan/MinimumBias2/HiForest_MinimumBias2_Run2015E_PromptReco_pp_5020GeV/151202_155519/0000/*.root >>  5p02TeV_MinBias2_forests.txt
ls   /mnt/hadoop/cms/store/user/dgulhan/MinimumBias6/HiForest_MinimumBias6_Run2015E_PromptReco_pp_5020GeV/151204_102344/000?/*.root >>  5p02TeV_MinBias6_forests.txt
ls   /mnt/hadoop/cms/store/user/dgulhan/MinimumBias10/HiForest_MinimumBias10_Run2015E_PromptReco_pp_5020GeV/151204_170246/000?/*.root >> 5p02TeV_MinBias10_forests.txt
## combined MinimumBias{1,2,6,10}
ls   /mnt/hadoop/cms/store/user/dgulhan/MinimumBias1/HiForest_MinimumBias1_Run2015E_PromptReco_pp_5020GeV/151201_113508/0000/*.root >> 5p02TeV_allMinBias_forests.txt
ls   /mnt/hadoop/cms/store/user/dgulhan/MinimumBias2/HiForest_MinimumBias2_Run2015E_PromptReco_pp_5020GeV/151202_155519/0000/*.root >> 5p02TeV_allMinBias_forests.txt
ls   /mnt/hadoop/cms/store/user/dgulhan/MinimumBias6/HiForest_MinimumBias6_Run2015E_PromptReco_pp_5020GeV/151204_102344/000?/*.root >>  5p02TeV_allMinBias_forests.txt
ls   /mnt/hadoop/cms/store/user/dgulhan/MinimumBias10/HiForest_MinimumBias10_Run2015E_PromptReco_pp_5020GeV/151204_170246/000?/*.root >> 5p02TeV_allMinBias_forests.txt
#ls   /mnt/hadoop/cms/store/user/dgulhan/MinimumBias6/HiForest_MinimumBias6_Run2015E_PromptReco_pp_5020GeV/151202_155519/0000/*.root >> 5p02TeV_allMinBias_forests.txt
#ls /mnt/hadoop/cms/store/user/dgulhan/MinimumBias10/HiForest_MinimumBias10_Run2015E_PromptReco_pp_5020GeV/151202_155519/0000/*.root >> 5p02TeV_allMinBias_forests.txt
##

