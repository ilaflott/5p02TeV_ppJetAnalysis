#!/bin/bash

#make_rooUnfold

rooUnfoldCompile bayesUnfoldDataSpectra_wNLO_etabin.C


############ SMEARED GEN PY 8 LO + HERWIG EE5C LO NPCs ###############

####use sigma/mu JER fit ---------------

##JER sys
./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_PY8_LO_sigmu_noJERscales_JERsys"       0          "smearTheory/PY8_ptLo28_drCut_geny_sigmu_noJERscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "PY8_ptLo28_drCut_geny_sigmu_noJERscales_spl3wgts_gaussSmear_"  1  "JER" ""

./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_PY8_LO_sigmu_noJERscales_JERsys"       1          "smearTheory/PY8_ptLo28_drCut_geny_sigmu_noJERscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "PY8_ptLo28_drCut_geny_sigmu_noJERscales_spl3wgts_gaussSmear_"  1  "JER" ""

./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_PY8_LO_sigmu_noJERscales_JERsys"       2          "smearTheory/PY8_ptLo28_drCut_geny_sigmu_noJERscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "PY8_ptLo28_drCut_geny_sigmu_noJERscales_spl3wgts_gaussSmear_"  1  "JER" ""

./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_PY8_LO_sigmu_noJERscales_JERsys"       3          "smearTheory/PY8_ptLo28_drCut_geny_sigmu_noJERscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "PY8_ptLo28_drCut_geny_sigmu_noJERscales_spl3wgts_gaussSmear_"  1  "JER" ""


scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_PY8_LO_sigmu_noJERscales_JERsys_NLOMC_wNP_??eta??.*"

return





############ NNPDF NLO + HERWIG EE5C LO NPCs ###############

####use sigma/mu JER fit ---------------

##vanilla no sys
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigmu_8TeVscales_nosys"                 "smearTheory/NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_spl3wgts_gaussSmear_"  0  "" ""
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigmu_nosys_NLOMC_wNP_00eta05.*"


#JER sys
./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigmu_noJERscales_JERsys" 0    "smearTheory/NNPDF_ptLo28_drCut_geny_sigmu_noJERscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigmu_noJERscales_spl3wgts_gaussSmear_"  1 "JER" ""
./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigmu_noJERscales_JERsys" 1    "smearTheory/NNPDF_ptLo28_drCut_geny_sigmu_noJERscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigmu_noJERscales_spl3wgts_gaussSmear_"  1 "JER" ""
./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigmu_noJERscales_JERsys" 2    "smearTheory/NNPDF_ptLo28_drCut_geny_sigmu_noJERscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigmu_noJERscales_spl3wgts_gaussSmear_"  1 "JER" ""
./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigmu_noJERscales_JERsys" 3    "smearTheory/NNPDF_ptLo28_drCut_geny_sigmu_noJERscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigmu_noJERscales_spl3wgts_gaussSmear_"  1 "JER" ""
scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigmu_noJERscales_JERsys_NLOMC_wNP_??eta??.*"
#return


##NP sys12, alt NP choices, HERWIG EE4C/PYTHIA8
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigmu_8TeVscales_NPsys12"    "smearTheory/NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_spl3wgts_gaussSmear_"  1 "NP" "12"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigmu_8TeVscales_NPsys12_NLOMC_wNP_00eta05.*"

##NP sysupdown, shifted HERWIG EE5C NP up/down, 1 sigma
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigmu_8TeVscales_NPsysupdown" "smearTheory/NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_spl3wgts_gaussSmear_"  1 "NP" "updown"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigmu_8TeVscales_NPsysupdown_NLOMC_wNP_00eta05.*"


##PDF sys12, alt PDF choices, CT14/HERA
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigmu_8TeVscales_PDFsys12"    "smearTheory/NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_spl3wgts_gaussSmear_"  1 "PDF" "12"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigmu_8TeVscales_PDFsys12_NLOMC_wNP_00eta05.*"
#return

##PDF sysupdown, shifted NNPDF PDF up/down, 6 Pt scale + PDF errs, 1 sigma
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigmu_8TeVscales_PDFsysupdown" "smearTheory/NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_spl3wgts_gaussSmear_"  1 "PDF" "updown"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigmu_8TeVscales_PDFsysupdown_NLOMC_wNP_00eta05.*"


##-----


##NP sys12, alt NP choices, POW+PY8/POW+PY8 CTEQ
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigmu_8TeVscales_v2_NPsys12"    "smearTheory/NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_v2_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_spl3wgts_gaussSmear_"  1 "NP" "12"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigmu_8TeVscales_v2_NPsys12_NLOMC_wNP_00eta05.*"

##NP sysupdown, shifted HERWIG EE5C NP up/down, 2 sigma
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigmu_8TeVscales_v2_NPsysupdown" "smearTheory/NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_v2_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_spl3wgts_gaussSmear_"  1 "NP" "updown"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigmu_8TeVscales_v2_NPsysupdown_NLOMC_wNP_00eta05.*"


##PDF sys12, alt PDF choices, CT10/MMHT
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigmu_8TeVscales_v2_PDFsys12"    "smearTheory/NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_v2_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_spl3wgts_gaussSmear_"  1 "PDF" "12"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigmu_8TeVscales_v2_PDFsys12_NLOMC_wNP_00eta05.*"

##PDF sysupdown, shifted NNPDF PDF up/down, 6 Pt scale + PDF errs, 2sigma
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigmu_8TeVscales_v2_PDFsysupdown" "smearTheory/NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_v2_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigmu_w8TeVscales_spl3wgts_gaussSmear_"  1 "PDF" "updown"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigmu_8TeVscales_v2_PDFsysupdown_NLOMC_wNP_00eta05.*"


####use sigma JER fit ---------------


##vanilla no sys
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigma_8TeVscales_nosys"                       "smearTheory/NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_spl3wgts_gaussSmear_"  0  "" ""
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigma_nosys_NLOMC_wNP_00eta05.*"


##JER sys
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigma_8TeVscales_JERsys"                       "smearTheory/NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_spl3wgts_gaussSmear_"  1 "JER" ""
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigma_8TeVscales_JERsys_NLOMC_wNP_00eta05.*"


##NP sys12, alt NP choices, HERWIG EE4C/PYTHIA8
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigma_8TeVscales_NPsys12"    "smearTheory/NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_spl3wgts_gaussSmear_"  1 "NP" "12"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigma_8TeVscales_NPsys12_NLOMC_wNP_00eta05.*"

##NP sysupdown, shifted HERWIG EE5C NP up/down, 1 sigma
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigma_8TeVscales_NPsysupdown" "smearTheory/NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_spl3wgts_gaussSmear_"  1 "NP" "updown"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigma_8TeVscales_NPsysupdown_NLOMC_wNP_00eta05.*"


##PDF sys12, alt PDF choices, CT14/HERA 
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigma_8TeVscales_PDFsys12"    "smearTheory/NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_spl3wgts_gaussSmear_"  1 "PDF" "12"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigma_8TeVscales_PDFsys12_NLOMC_wNP_00eta05.*"

##PDF sysupdown, shifted NNPDF PDF up/down, 6 Pt scale + PDF errs, 1 sigma
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigma_8TeVscales_PDFsysupdown" "smearTheory/NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_spl3wgts_gaussSmear_"  1 "PDF" "updown"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigma_8TeVscales_PDFsysupdown_NLOMC_wNP_00eta05.*"


##-----


##NP sys12, alt NP choices, POW+PY8/POW+PY8 CTEQ
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigma_8TeVscales_v2_NPsys12"    "smearTheory/NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_v2_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_spl3wgts_gaussSmear_"  1 "NP" "12"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigma_8TeVscales_v2_NPsys12_NLOMC_wNP_00eta05.*"

##NP sysupdown, shifted HERWIG EE5C NP up/down, 2 sigma
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigma_8TeVscales_v2_NPsysupdown" "smearTheory/NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_v2_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_spl3wgts_gaussSmear_"  1 "NP" "updown"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigma_8TeVscales_v2_NPsysupdown_NLOMC_wNP_00eta05.*"


##PDF sys12, alt PDF choices, CT10/MMHT
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigma_8TeVscales_v2_PDFsys12"    "smearTheory/NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_v2_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_spl3wgts_gaussSmear_"  1 "PDF" "12"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigma_8TeVscales_v2_PDFsys12_NLOMC_wNP_00eta05.*"

##PDF sysupdown, shifted NNPDF PDF up/down, 6 Pt scale + PDF errs, 2 sigma
#./bayesUnfoldDataSpectra_wNLO_etabin.exe "02.18.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-18-19_jetPlots_0.0eta3.0" "Bayes_NNPDF_NLO_sigma_8TeVscales_v2_PDFsysupdown" "smearTheory/NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_v2_fullstat_spl3wgts_gaussSmear_NLO_plots/" "NNPDF_ptLo28_drCut_geny_sigma_w8TeVscales_spl3wgts_gaussSmear_"  1 "PDF" "updown"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_NNPDF_NLO_sigma_8TeVscales_v2_PDFsysupdown_NLOMC_wNP_00eta05.*"









#------------------------------#


#rooUnfoldCompile SVDUnfoldDataSpectra_wNLO_etabin.C

#./SVDUnfoldDataSpectra_wNLO_etabin.exe "test"

#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_*_*_NLOMC_wNP_??eta??*.*"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_Bayes_HLT60Thresh*_NLOMC_wNP_00eta05*"