#!/bin/bash


dojetQAPlots=0
dojetTrigQAPlots=0   
dojetTrigEff=0
dojetIDQAPlots=1
dojetMultPlots=0   #old
doevtCounts=0
scpOutput=0




if [[ $dojetQAPlots -eq 1 ]]
then
    
    echo ""
    echo "compiling printPlots_jetPlots.C"
    echo ""
    
    rootcompile printPlots_jetPlots.C
    
    ##############################################
    
    echo ""
    echo "done compiling. Running!"
    echo ""    
    
    #####SPECIAL######
    #source run_printPlots_jetPlots.sh "4" "06.25.19_outputCondor" "06-25-19" "0.0eta2.0_SMPbins_withJECsys_semifinal"    "06.25.19_outputCondor" "06-25-19" "0.0eta2.0_semifinal"  "00eta20_06.25.19_noMB_SMPbins_semifinal"    
    #scp2Serin "output/ak4PF_HPtJetTrig_semiOffPy8_00eta20_06.25.19_noMB_SMPbins_semifinal_jetPlots.*"
    #return
    
    #####SPECIAL######
#    source run_printPlots_jetPlots.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal"    "06.25.19_outputCondor" "06-25-19" "0.0eta2.0_semifinal"  "00eta20_08.06.19_noMB_SMPbins_withLumiCorr_semifinal"    
 #   scp2Serin "output/ak4PF_HPtJetTrig_semiOffPy8_00eta20_08.06.19_noMB_SMPbins_withLumiCorr_semifinal_jetPlots.*"
    
    #source run_printPlots_jetPlots.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal"    "08.19.19_outputCondor" "08-19-19" "0.0eta2.0_noDRCut_semifinal"  "00eta20_08.06.19_noMB_SMPbins_withLumiCorr_noDRCut_semifinal"    
    #scp2Serin "output/ak4PF_HPtJetTrig_semiOffPy8_00eta20_08.06.19_noMB_SMPbins_withLumiCorr_noDRCut_semifinal_jetPlots.*"

    #source run_printPlots_jetPlots.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal"    "08.19.19_outputCondor" "08-19-19" "0.0eta1.0_noDRCut_semifinal_00eta10Only"  "00eta10_08.06.19_noMB_SMPbins_withLumiCorr_noDRCut_semifinal"    
    #scp2Serin "output/ak4PF_HPtJetTrig_semiOffPy8_00eta10_08.06.19_noMB_SMPbins_withLumiCorr_noDRCut_semifinal_jetPlots.*"

    
    ## BUGGED OR MISSING HISTS I WANT
    #source run_printPlots_jetPlots.sh "4" "11.21.19_outputCondor" "11-21-19" "0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinalv2"                   "06.25.19_outputCondor" "06-25-19" "0.0eta2.0_semifinal"  "00eta20_11.21.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2"    
    #scp2Serin "output/ak4PF_HPtJetTrig_semiOffPy8_00eta20_11.21.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2_jetPlots.*"
    ## BUGGED OR MISSING HISTS I WANT
    #source run_printPlots_jetPlots.sh "4" "11.21.19_outputCondor" "11-21-19" "0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh"    "06.25.19_outputCondor" "06-25-19" "0.0eta2.0_semifinal"  "00eta20_11.21.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh"    
    #scp2Serin "output/ak4PF_HPtJetTrig_semiOffPy8_00eta20_11.21.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh_jetPlots.*"

    #### SEMIFINAL V2
    #source run_printPlots_jetPlots.sh "4" "11.26.19_outputCondor" "11-26-19" "0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinalv2"                   "06.25.19_outputCondor" "06-25-19" "0.0eta2.0_semifinal"  "00eta20_11.26.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2"    
    #scp2Serin "output/ak4PF_HPtJetTrig_semiOffPy8_00eta20_11.26.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2_jetPlots.*"
    #### SEMIFINAL V2 LOWHLT40THRESH
    #source run_printPlots_jetPlots.sh "4" "11.26.19_outputCondor" "11-26-19" "0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh"    "06.25.19_outputCondor" "06-25-19" "0.0eta2.0_semifinal"  "00eta20_11.26.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh"    
    #scp2Serin "output/ak4PF_HPtJetTrig_semiOffPy8_00eta20_11.26.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh_jetPlots.*"


    #### SEMIFINAL V3, 03.18.20 outputCondor results #ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_03-18-20_jetPlots_0.0eta2.5_semifinalv3/ #ppData_HighPtJetTrig_ak4PFJets_03-18-20_jetPlots_0.0eta2.5_SMPbins_withjety_semifinalv3/
    #source run_printPlots_jetPlots.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_semifinalv3" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_semifinalv3"  "00eta20_03.18.20_semifinalv3"    

    ### SEMIFINAL V4, 01.05.21 outputCondor results 
    source run_printPlots_jetPlots.sh "4" "01.05.21_outputCondor" "01-05-21" "0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4" "01.05.21_outputCondor" "01-05-21" "0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"  "00eta20_01.05.21_semifinalv4"    

    ##############################################
    
    echo ""
    echo "done running jetPlots!"
    echo ""
    

fi




if [[ $dojetTrigQAPlots -eq 1 ]]
then
    echo ""
    echo "compiling printPlots_jetTrig.C"
    echo ""
    JETTRIGOUTDIR=/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output
    rootcompile printPlots_jetTrig.C
    
    ##############################################
    
    echo ""
    echo "done compiling. Running!"
    echo ""    

    
    #####SPECIAL
    #source run_printPlots_jetTrig.sh "4" "06.25.19_outputCondor" "06-25-19" "0.0eta2.0_SMPbins_withJECsys_semifinal"      "00eta20_06.25.19_noMB_SMPbins_semifinal"

    #source run_printPlots_jetTrig.sh "4" "07.29.19_outputCondor" "07-29-19" "0.0eta1.0_lumiRunStudy_usualDupeSkip"      "00eta10_07.29.19_noMB_lumiRunStudy_usualDupeSkip_semifinal"
    #source run_printPlots_jetTrig.sh "4" "07.29.19_outputCondor" "07-29-19" "0.0eta1.0_lumiRunStudy_diffDupeSkip"      "00eta10_07.29.19_noMB_lumiRunStudy_diffDupeSkip_semifinal"
    #source run_printPlots_jetTrig.sh "4" "08.05.19_outputCondor" "08-05-19" "0.0eta2.0_lumiRunStudy_lumieffcalcOnly"      "00eta20_08.05.19_noMB_SMPbins_semifinal_SMPTrigCombo"
    
    ####SPECIAL
    #source run_printPlots_jetTrig.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal"      "00eta20_08.06.19_noMB_SMPbins_withLumiCorr_semifinal"
    #scp2Serin "output/ak4PF_HPtJetTrig_00eta20_08.06.19_noMB_SMPbins_withLumiCorr_semifinal_jetTrig.*"


    ###SPECIAL
    source run_printPlots_jetTrig.sh "4" "11.26.19_outputCondor" "11-26-19" "0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinalv2"      "00eta20_11.26.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2"
    #scp2Serin "output/ak4PF_HPtJetTrig_00eta20_08.06.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2_jetTrig.*"

    #source run_printPlots_jetTrig.sh "4" "11.26.19_outputCondor" "11-26-19" "0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh"      "00eta20_11.26.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh"
    #scp2Serin "output/ak4PF_HPtJetTrig_00eta20_08.06.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh_jetTrig.*"

    
    ##############################################
    
    echo ""
    echo "done running jetTrig!"
    echo ""
    

    

fi




if [[ $dojetTrigEff -eq 1 ]]
then
    echo ""
    echo "compiling printPlots_jetTrigEff.C"
    echo ""
    rootcompile printPlots_jetTrigEff.C
    echo ""
    echo "done compiling. Running!"
    echo ""

    ##FIRST ATTEMPT, WORKS
    #source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal" "00eta20_08.06.19_HLT80Eff_useHLT60Ref_firstAttempt" "HLT80" "HLT60"
    #source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal" "00eta20_08.06.19_HLT60Eff_useHLT40Ref_firstAttempt" "HLT60" "HLT40"
    #scp2Serin "output/ak4PF_HPtJetTrig_00eta20_08.06.19_HLT??Eff_useHLT??Ref_firstAttempt_jetTrigEff.*"
    
    
    ##SECOND ATTEMPT, WORKS, OUTPUT IDENTICAL TO FIRST ATTEMPT (WHAT WAS THE DIFFERENCE?!)
    #source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_secondAttempt" "0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal" "00eta20_08.06.19_HLT80Eff_useHLT60Ref_secondAttempt" "HLT80" "HLT60"
    #source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_secondAttempt" "0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal" "00eta20_08.06.19_HLT60Eff_useHLT40Ref_secondAttempt" "HLT60" "HLT40"
    #scp2Serin "output/ak4PF_HPtJetTrig_00eta20_08.06.19_HLT??Eff_useHLT??Ref_secondAttempt_jetTrigEff.*"
    
    
    ##THIRD ATTEMPT, WORKS, COMBO HIGHER THAN INDIV LEAD TRIG JET SPECTRA, WHY?
    #source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta5.1_forHLTTurnOns_thirdAttempt" "0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal" "00eta51_08.06.19_HLT80Eff_useHLT60Ref_thirdAttempt" "HLT80" "HLT60"
    #source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta5.1_forHLTTurnOns_thirdAttempt" "0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal" "00eta51_08.06.19_HLT60Eff_useHLT40Ref_thirdAttempt" "HLT60" "HLT40"
    #scp2Serin "output/ak4PF_HPtJetTrig_00eta51_08.06.19_HLT??Eff_useHLT??Ref_thirdAttempt_jetTrigEff.*"


    
    #SEMIFINAL V2
    #source run_printPlots_jetTrigEff.sh "4" "11.26.19_outputCondor" "11-26-19" "0.0eta5.1_semifinalv2" "0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinalv2" "00eta51_11.26.19_HLT80Eff_useHLT60Ref_semifinalv2" "HLT80" "HLT60"
    #source run_printPlots_jetTrigEff.sh "4" "11.26.19_outputCondor" "11-26-19" "0.0eta5.1_semifinalv2" "0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinalv2" "00eta51_11.26.19_HLT60Eff_useHLT40Ref_semifinalv2" "HLT60" "HLT40"
    #scp2Serin "output/ak4PF_HPtJetTrig_00eta51_11.26.19_HLT??Eff_useHLT??Ref_semifinalv2_jetTrigEff.*"

    ##SEMIFINAL V2 LOWHLT40 THRESH
#    source run_printPlots_jetTrigEff.sh "4" "11.26.19_outputCondor" "11-26-19" "0.0eta5.1_semifinalv2_LowHLT40Thresh" "0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh" "00eta51_11.26.19_HLT80Eff_useHLT60Ref_semifinalv2_LowHLT40Thresh" "HLT80" "HLT60"
#    source run_printPlots_jetTrigEff.sh "4" "11.26.19_outputCondor" "11-26-19" "0.0eta5.1_semifinalv2_LowHLT40Thresh" "0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh" "00eta51_11.26.19_HLT60Eff_useHLT40Ref_semifinalv2_LowHLT40Thresh" "HLT60" "HLT40"
    #scp2Serin "output/ak4PF_HPtJetTrig_00eta51_11.26.19_HLT??Eff_useHLT??Ref_semifinalv2_LowHLT40Thresh_jetTrigEff.*"

    ##SEMIFINAL V4 TIGHT JET ID, HBHEISONOISEFILTER, PFMETFRAC CUT
    source run_printPlots_jetTrigEff.sh "4" "01.05.21_outputCondor" "01-05-21"        "0.0eta5.1_wHBHEIsoNoise_wPFMETFrac0p3_tightJetID_semifinalv4"      "0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_tightJetID_semifinalv4" "00eta51_01.05.21_HLT80Eff_useHLT60Ref_semifinalv4" "HLT80" "HLT60"
    source run_printPlots_jetTrigEff.sh "4" "01.05.21_outputCondor" "01-05-21"        "0.0eta5.1_wHBHEIsoNoise_wPFMETFrac0p3_tightJetID_semifinalv4"      "0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_tightJetID_semifinalv4" "00eta51_01.05.21_HLT60Eff_useHLT40Ref_semifinalv4" "HLT60" "HLT40"





    ##############################################
    
    echo ""
    echo "done running jetTrigEff!"
    echo ""    
    

fi





if [[ $dojetIDQAPlots -eq 1 ]]
then
    #####SPECIAL
    #####SPECIAL
    #####SPECIAL    
    #echo ""
    #echo "compiling printPlots_HINvSMPjetIDPlots.C"
    #echo ""
    #
    #rootcompile printPlots_HINvSMPjetIDPlots.C
    #
    ###############################################
    #
    #echo ""
    #echo "done compiling. Running!"
    #echo ""
    #
    #
    #rootcompile printPlots_HINvSMPjetIDPlots.C
    #source run_printPlots_HINvSMPjetIDPlots.sh "4" "08.06.19_outputCondor" "08-06-19"  "0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal"   "0.0eta2.0_useTupel_withLumiCorr" "00eta20_08.06.19_useIncJetAna_withLumiCorr_semifinal" "ppData"
    #scp2Serin "/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta20_08.06.19_useIncJetAna_withLumiCorr_semifinal_HINvSMPjetIDPlots.pdf"
    #####SPECIAL
    #####SPECIAL
    #####SPECIAL

    
    ####SPECIAL
    ####SPECIAL
    ####SPECIAL    
    #echo ""
    #echo "compiling printPlots_tightVloose_jetIDPlots.C"
    #echo ""
    #
    #rootcompile printPlots_tightVloose_jetIDPlots.C
    #
    ###############################################
    #
    #echo ""
    #echo "done compiling. Running!"
    #echo ""
    #
    #
    #source run_printPlots_tightVloose_jetIDPlots.sh "4" "08.06.19_outputCondor" "08-06-19"  "0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinal"   "0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_tightJetID_semifinal" "00eta20_08.06.19_tightVloose_withLumiCorr_semifinal" "ppData"
    #scp2Serin "/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta20_08.06.19_tightVloose_withLumiCorr_semifinal_tightVloose_jetIDPlots.pdf"
    ####SPECIAL
    ####SPECIAL
    ####SPECIAL





    echo ""
    echo "compiling printPlots_jetIDPlots.C"
    echo ""
    
    rootcompile printPlots_jetIDPlots.C
    
    ##############################################
    
    echo ""
    echo "done compiling. Running!"
    echo ""
    
    
    #source run_printPlots_jetIDPlots.sh "4" "08.06.19_outputCondor" "08-06-19"  "0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal"  "00eta20_08.06.19_useIncJetAna_withLumiCorr_semifinal" "ppData"
    #scp2Serin "output/ak4PF_HPtJetTrig_00eta20_08.06.19_useIncJetAna_withLumiCorr_semifinal_jetIDPlots.*"
    
    #source run_printPlots_jetIDPlots.sh "4" "08.06.19_outputCondor" "08-06-19"  "0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_tightJetID_semifinal"  "00eta20_08.06.19_useIncJetAna_withLumiCorr_tightJetID_semifinal" "ppData"
    #scp2Serin "output/ak4PF_HPtJetTrig_00eta20_08.06.19_useIncJetAna_withLumiCorr_tightJetID_semifinal_jetIDPlots.*"
    
    #source run_printPlots_jetIDPlots.sh "4" "08.06.19_outputCondor" "08-06-19"  "0.0eta2.0_useTupel_withLumiCorr"  "00eta20_08.06.19_useTupel_withLumiCorr_semifinal" "ppData"
    #scp2Serin "output/ak4PF_HPtJetTrig_00eta20_08.06.19_useTupel_withLumiCorr_semifinal_jetIDPlots.*"

    #source run_printPlots_jetIDPlots.sh "4" "11.26.19_outputCondor" "11-26-19"  "0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinalv2"  "00eta20_11.26.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2"     "ppData"
    #scp2Serin "output/ak4PF_HPtJetTrig_00eta20_11.26.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2_jetIDPlots.*"    
    #
    #source run_printPlots_jetIDPlots.sh "4" "11.26.19_outputCondor" "11-26-19"  "0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh"  "00eta20_11.26.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh"     "ppData"
    #scp2Serin "output/ak4PF_HPtJetTrig_00eta20_11.26.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh_jetIDPlots.*"    
    
#    source run_printPlots_jetIDPlots.sh "4" "03.18.20_outputCondor" "03-18-20"  "0.0eta2.5_SMPbins_withjety_semifinalv3"  "00eta20_03.18.20_semifinalv3"     "ppData"
#    source run_printPlots_jetIDPlots.sh "4" "03.18.20_outputCondor" "03-18-20"  "0.0eta2.5_SMPbins_withjety_semifinalv4"  "00eta20_03.18.20_semifinalv4"     "ppData"
#    source run_printPlots_jetIDPlots.sh "4" "03.18.20_outputCondor" "03-18-20"  "0.0eta2.5_SMPbins_withjety_semifinalv4"  "00eta20_03.18.20_jet80PD_semifinalv4"     "ppData"
#    source run_printPlots_jetIDPlots.sh "4" "03.18.20_outputCondor" "03-18-20"  "0.0eta2.5_SMPbins_withjety_semifinalv4"  "00eta20_03.18.20_lowjetsPD_semifinalv4"     "ppData"
#    source run_printPlots_jetIDPlots.sh "4" "03.18.20_outputCondor" "03-18-20"  "0.0eta2.5_semifinalv3"  "00eta20_03.18.20_semifinalv4"     "ppMC"
    
    #source run_printPlots_jetIDPlots.sh "4" "06.25.19_outputCondor" "06-25-19"  "0.0eta2.0_semifinal"  "00eta20_06.25.19_SMPbins_semifinal" "ppMC"
    #scp2Serin "output/ak4PF_Py8_CUETP8M1_00eta20_06.25.19_SMPbins_semifinal_jetIDPlots.*"

    #source run_printPlots_jetIDPlots.sh "4" "08.19.19_outputCondor" "08-19-19" "0.0eta2.0_noDRCut_semifinal"  "00eta20_08.19.19_SMPbins_noDRCut_semifinal"    "ppMC"
    #scp2Serin "output/ak4PF_Py8_CUETP8M1_00eta20_08.19.19_SMPbins_noDRCut_semifinal_jetIDPlots.*"

    #source run_printPlots_jetIDPlots.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_noDRCut_tightJetID_semifinal"  "00eta20_08.06.19_SMPbins_noDRCut_tightJetID_semifinal"    "ppMC"
    #scp2Serin "output/ak4PF_Py8_CUETP8M1_00eta20_08.06.19_SMPbins_noDRCut_tightJetID_semifinal_jetIDPlots.*"

    #source run_printPlots_jetIDPlots.sh "4" "08.19.19_outputCondor" "08-19-19" "0.0eta1.0_noDRCut_semifinal_00eta10Only"  "00eta10_08.19.19_SMPbins_noDRCut_semifinal"    "ppMC"
    #scp2Serin "output/ak4PF_Py8_CUETP8M1_00eta10_08.19.19_SMPbins_noDRCut_semifinal_jetIDPlots.*"

    source run_printPlots_jetIDPlots.sh "4" "01.05.21_outputCondor" "01-05-21"  "0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"  "00eta20_01.05.21_semifinalv4"     "ppData"
    #    source run_printPlots_jetIDPlots.sh "4" "01.05.21_outputCondor" "01-05-21"  "0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"  "00eta20_01.05.21_semifinalv4"     "ppMC"
    
    ##############################################
    
    echo ""
    echo "done running jetIDPlots!"
    echo ""

    #scp2Serin output/\*.\*

    ##############################################
    
    #return
fi






if [[ $doevtCounts -eq 1 ]]
then
    
    echo ""
    echo "compiling printPlots_evtCounts.C"
    echo ""
    
    rootcompile printPlots_evtCounts.C
    
    ##############################################
    
    echo ""
    echo "done compiling. Running!"
    echo ""    

    #source run_printPlots_evtCounts.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_JetID_semifinalv4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withgenjety_wJERSFs_semifinalv3"  "00eta20_03.18.20_semifinalv3_Jet80_LowJets_PY8"    
    source run_printPlots_evtCounts.sh "4" "01.05.21_outputCondor" "01-05-21" "0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4" "01.05.21_outputCondor" "01-05-21" "0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"  "00eta20_01.05.21_semifinalv4_Jet80_LowJets_PY8"    

    echo ""
    echo "done running evtCounts!"
    echo ""


fi



















if [[ $dojetMultPlots -eq 1 ]]
then
    
    echo ""
    echo "compiling printPlots_jetMult.C"
    echo ""
    
    rootcompile printPlots_jetMultPlots.C

    ##############################################
    
    echo ""
    echo "done compiling. Running!"
    echo ""
    
    #source run_printPlots_jetMultPlots.sh "4" "10.27.17_outputCondor" "10-27-17" "0.0eta2.0" "10-27-17" "0.0eta2.0" "00eta20_firstJetMultPlots"
    #source run_printPlots_jetMultPlots.sh "4" "10.27.17_outputCondor" "10-27-17" "0.0eta2.0" "00eta20_first"
    #source run_printPlots_jetMultPlots.sh "4" "12.19.17_outputCondor" "12-19-17" "0.0eta2.0" "00eta20"
    #source run_printPlots_jetMultPlots.sh "4" "01.22-23.18_outputCondor" "01-22-18" "0.0eta2.0" "01.22-23.18_outputCondor" "01-22-18" "0.0eta2.0" "00eta20_ppMCFirst" 
    source run_printPlots_jetMultPlots.sh "4" "01.22-23.18_outputCondor" "01-22-18" "0.0eta2.0" "01.22-23.18_outputCondor" "01-22-18" "0.0eta2.0" "00eta20_ppMC_tksNotMatched" 
    source run_printPlots_jetMultPlots.sh "4" "01.22-23.18_outputCondor" "01-22-18" "0.0eta2.0" "01.22-23.18_outputCondor" "01-22-18" "matchedTracks_0.0eta2.0" "00eta20_ppMC_tksMatched" 
    
    
    ##############################################
    
    echo ""
    echo "done running jetMultPlots!"
    echo ""
    
    ###############################################
    
    #return
    
    ###############################################
fi





return

















    #scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta30_04.12.19_wMB_jetTrig.pdf
    #scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta30_04.12.19_noMB_jetTrig.pdf
    #scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta30_04.12.19_wJetID_wMB_jetTrig.pdf
    #scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta30_04.12.19_wJetID_noMB_jetTrig.pdf
    #scp2Serin output/\*jetTrig\*
    #scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta20_11.06.18_good_jetTrig.pdf
    #scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta20_11.14.18_trgPt\*.pdf
    #scp2Serin "/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta?0_01.??.19_HLT60Thresh??_jetTrig.pdf"
    ##############################################
    
    #return



#if [[ $dojetTrigEff -eq 1 ]]
#then
#    echo ""
#    echo "compiling printPlots_jetTrigEff.C"
#    echo ""
#    rootcompile printPlots_jetTrigEff.C
#    echo ""
#    echo "done compiling. Running!"
#    echo ""
#    
#    ############################################### USE LOWERJETS PD ONLY
#    source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtLowerJets_ak4PF-allFiles.root" "00eta20_08.06.19_HLT80Eff_useHLT60Ref_LowerJetsPD" "HLT80" "HLT60"
#    source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtLowerJets_ak4PF-allFiles.root" "00eta20_08.06.19_HLT80Eff_useHLT40Ref_LowerJetsPD" "HLT80" "HLT40"
#    source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtLowerJets_ak4PF-allFiles.root" "00eta20_08.06.19_HLT60Eff_useHLT40Ref_LowerJetsPD" "HLT60" "HLT40"
#    
#    ############################################## USE JET80 PD ONLY
#    source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtJet80_ak4PF-allFiles.root" "00eta20_08.06.19_HLT80Eff_useHLT60Ref_Jet80PD" "HLT80" "HLT60"
#    source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtJet80_ak4PF-allFiles.root" "00eta20_08.06.19_HLT80Eff_useHLT40Ref_Jet80PD" "HLT80" "HLT40"
#    source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtJet80_ak4PF-allFiles.root" "00eta20_08.06.19_HLT60Eff_useHLT40Ref_Jet80PD" "HLT60" "HLT40"
#    
#
#    ################################################ USE BOTH (SHOULDNT WORK RIGHT NOW BUT I WANNA SEE)
#    #source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtJetTrig_noMB_ak4PF-allFiles.root" "00eta20_08.06.19_HLT80Eff_useHLT60Ref_JetTrigPD" "HLT80" "HLT60"
#    #source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtJetTrig_noMB_ak4PF-allFiles.root" "00eta20_08.06.19_HLT80Eff_useHLT40Ref_JetTrigPD" "HLT80" "HLT40"
#    #source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtJetTrig_noMB_ak4PF-allFiles.root" "00eta20_08.06.19_HLT60Eff_useHLT40Ref_JetTrigPD" "HLT60" "HLT40"
#
#    
#    scp2Serin "output/ak4PF_HPtJetTrig_00eta20_08.06.19_HLT??Eff_useHLT??Ref_*PD_jetTrigEff.*"
#    ##############################################
#    
#    echo ""
#    echo "done running jetTrigEff!"
#    echo ""    
#    
#
#if

#if [[ $scpOutput -eq 1 ]]
#then
#    #scp2Serin output/\*
#    if [[ $dojetQAPlots -eq 1 ]]
#    then
#	scp2Serin output/\*jetPlots\*
#    fi
#    if [[ $dojetTrigQAPlots -eq 1 ]]
#    then
#	scp2Serin output/\*jetTrig.\*
#    fi
##    if [[ $dojetIDQAPlots -eq 1 ]]
##    then
##	scp2Serin output/
##    fi
##    if [[ $dojetMultPlots -eq 1 ]]
##    then
##	scp2Serin output/
##    fi
#fi
