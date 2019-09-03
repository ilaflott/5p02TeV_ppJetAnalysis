#!/bin/bash


dojetQAPlots=0        #done for ppMC + (ppData from 08.06.19)
                      #looks good
dojetTrigQAPlots=0   #done for ppData, 08.06.19
dojetTrigEff=1                      #looks good
dojetIDQAPlots=0      #done for ppMC, DO FOR ppData, 08.06.19, look at ouput

dojetMultPlots=0
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
    
    source run_printPlots_jetPlots.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal"    "08.19.19_outputCondor" "08-19-19" "0.0eta2.0_noDRCut_semifinal"  "00eta20_08.06.19_noMB_SMPbins_withLumiCorr_noDRCut_semifinal"    
    scp2Serin "output/ak4PF_HPtJetTrig_semiOffPy8_00eta20_08.06.19_noMB_SMPbins_withLumiCorr_noDRCut_semifinal_jetPlots.*"




    #source run_printPlots_jetPlots.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal"    "08.19.19_outputCondor" "08-19-19" "0.0eta1.0_noDRCut_semifinal_00eta10Only"  "00eta10_08.06.19_noMB_SMPbins_withLumiCorr_noDRCut_semifinal"    
   # scp2Serin "output/ak4PF_HPtJetTrig_semiOffPy8_00eta10_08.06.19_noMB_SMPbins_withLumiCorr_noDRCut_semifinal_jetPlots.*"

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
    source run_printPlots_jetTrig.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal"      "00eta20_08.06.19_noMB_SMPbins_withLumiCorr_semifinal"
    scp2Serin "output/ak4PF_HPtJetTrig_00eta20_08.06.19_noMB_SMPbins_withLumiCorr_semifinal_jetTrig.*"
    
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
    
    ############################################## USE JET80 PD ONLY
    source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtJet80_ak4PF-allFiles.root" "00eta20_08.06.19_Jet80PD_HLT80Eff_useHLT60Ref" "HLT80" "HLT60"
    source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtJet80_ak4PF-allFiles.root" "00eta20_08.06.19_Jet80PD_HLT80Eff_useHLT40Ref" "HLT80" "HLT40"
    source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtJet80_ak4PF-allFiles.root" "00eta20_08.06.19_Jet80PD_HLT60Eff_useHLT40Ref" "HLT60" "HLT40"
    
    ############################################### USE LOWERJETS PD ONLY
    source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtLowerJets_ak4PF-allFiles.root" "00eta20_08.06.19_LowerJetsPD_HLT80Eff_useHLT60Ref" "HLT80" "HLT60"
    source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtLowerJets_ak4PF-allFiles.root" "00eta20_08.06.19_LowerJetsPD_HLT80Eff_useHLT40Ref" "HLT80" "HLT40"
    source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtLowerJets_ak4PF-allFiles.root" "00eta20_08.06.19_LowerJetsPD_HLT60Eff_useHLT40Ref" "HLT60" "HLT40"
    
    ############################################### USE BOTH (SHOULDNT WORK RIGHT NOW BUT I WANNA SEE)
    source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtJetTrig_noMB_ak4PF-allFiles.root" "00eta20_08.06.19_JetTrigPD_HLT80Eff_useHLT60Ref" "HLT80" "HLT60"
    source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtJetTrig_noMB_ak4PF-allFiles.root" "00eta20_08.06.19_JetTrigPD_HLT80Eff_useHLT40Ref" "HLT80" "HLT40"
    source run_printPlots_jetTrigEff.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_forHLTTurnOns_firstAttempt" "HighPtJetTrig_noMB_ak4PF-allFiles.root" "00eta20_08.06.19_JetTrigPD_HLT60Eff_useHLT40Ref" "HLT60" "HLT40"

    
    scp2Serin "output/ak4PF_HPtJetTrig_00eta20_08.06.19_*PD_HLT??Eff_useHLT??Ref_jetTrigEff.*"
    ##############################################
    
    echo ""
    echo "done running jetTrigEff!"
    echo ""    
    

fi





if [[ $dojetIDQAPlots -eq 1 ]]
then
    
    echo ""
#    echo "compiling printPlots_jetIDPlots.C"
    echo "compiling printPlots_jetIDPlots_v2.C"
    echo ""
    
    #rootcompile printPlots_jetIDPlots.C
    rootcompile printPlots_jetIDPlots_v2.C

    ##############################################
    
    echo ""
    echo "done compiling. Running!"
    echo ""
    

    #####SPECIAL
    #source run_printPlots_jetIDPlots_v2.sh "4" "06.25.19_outputCondor" "06-25-19"  "0.0eta2.0_SMPbins_withJECsys_semifinal"  "00eta20_06.25.19_noMB_SMPbins_semifinal" "ppData"
    #source run_printPlots_jetIDPlots_v2.sh "4" "08.06.19_outputCondor" "08-06-19"  "0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal"  "00eta20_08.06.19_noMB_SMPbins_withLumiCorr_semifinal" "ppData"
    #scp2Serin "output/ak4PF_HPtJetTrig_00eta20_08.06.19_noMB_SMPbins_withLumiCorr_semifinal_jetIDPlots.*"
    




    #source run_printPlots_jetIDPlots_v2.sh "4" "06.25.19_outputCondor" "06-25-19"  "0.0eta2.0_semifinal"  "00eta20_06.25.19_SMPbins_semifinal" "ppMC"
    
#   source run_printPlots_jetPlots.sh "4" "08.06.19_outputCondor" "08-06-19" "0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal"    "08.19.19_outputCondor" "08-19-19" "0.0eta2.0_noDRCut_semifinal"  "00eta20_08.06.19_noMB_SMPbins_withLumiCorr_noDRCut_semifinal"    
    
    source run_printPlots_jetIDPlots_v2.sh "4" "08.19.19_outputCondor" "08-19-19" "0.0eta2.0_noDRCut_semifinal"  "00eta20_08.19.19_SMPbins_noDRCut_semifinal"    "ppMC"
    scp2Serin "output/ak4PF_Py8_CUETP8M1_00eta20_08.19.19_SMPbins_noDRCut_semifinal_jetIDPlots.*"
    #source run_printPlots_jetIDPlots_v2.sh "4" "08.19.19_outputCondor" "08-19-19" "0.0eta1.0_noDRCut_semifinal_00eta10Only"  "00eta10_08.19.19_SMPbins_noDRCut_semifinal"    "ppMC"
    #scp2Serin "output/ak4PF_Py8_CUETP8M1_00eta10_08.19.19_SMPbins_noDRCut_semifinal_jetIDPlots.*"

    
    
    ##############################################
    
    echo ""
    echo "done running jetIDPlots!"
    echo ""

    #scp2Serin output/\*.\*

    ##############################################
    
    #return
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



if [[ $scpOutput -eq 1 ]]
then
    #scp2Serin output/\*
    if [[ $dojetQAPlots -eq 1 ]]
    then
	scp2Serin output/\*jetPlots\*
    fi
    if [[ $dojetTrigQAPlots -eq 1 ]]
    then
	scp2Serin output/\*jetTrig.\*
    fi
#    if [[ $dojetIDQAPlots -eq 1 ]]
#    then
#	scp2Serin output/
#    fi
#    if [[ $dojetMultPlots -eq 1 ]]
#    then
#	scp2Serin output/
#    fi
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



