#!/bin/bash


dojetQAPlots=0
dojetTrigQAPlots=1
dojetIDQAPlots=0
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
    
    #source run_printPlots_jetPlots.sh "4" "11.06.18_outputCondor" "11-06-18" "0.0eta0.5_wjetid" "11.06.18_outputCondor" "11-06-18" "0.0eta0.5_wjetid" "00eta05_wJetID"
    #source run_printPlots_jetPlots.sh "4" "11.06.18_outputCondor" "11-06-18" "0.5eta1.0_wjetid" "11.06.18_outputCondor" "11-06-18" "0.5eta1.0_wjetid" "05eta10_wJetID"
    #source run_printPlots_jetPlots.sh "4" "11.06.18_outputCondor" "11-06-18" "1.0eta1.5_wjetid" "11.06.18_outputCondor" "11-06-18" "1.0eta1.5_wjetid" "10eta15_wJetID"
    #source run_printPlots_jetPlots.sh "4" "11.06.18_outputCondor" "11-06-18" "1.5eta2.0_wjetid" "11.06.18_outputCondor" "11-06-18" "1.5eta2.0_wjetid" "15eta20_wJetID"    
    #source run_printPlots_jetPlots.sh "4" "11.06.18_outputCondor" "11-06-18" "0.0eta2.0_wjetid" "11.06.18_outputCondor" "11-06-18" "0.0eta2.0_wjetid" "00eta20_wJetID"
    
    #source run_printPlots_jetPlots.sh "4" "11.14.18_outputCondor" "11-14-18" "0.0eta2.0_trgPtNoRes"   "11.14.18_outputCondor" "11-14-18" "0.0eta2.0" "00eta20_trgPtNoRes"
    #source run_printPlots_jetPlots.sh "4" "11.14.18_outputCondor" "11-14-18" "0.0eta2.0_trgPtNoRes_NoJetID"   "11.14.18_outputCondor" "11-14-18" "0.0eta2.0_NoJetID" "00eta20_trgPtNoRes_NoJetID"

    #source run_printPlots_jetPlots.sh "4" "11.27.18_outputCondor" "11-27-18" "0.0eta2.0"   "11.27.18_outputCondor" "11-27-18" "0.0eta2.0" "00eta20"
    #source run_printPlots_jetPlots.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta3.0"   "11.27.18_outputCondor" "11-27-18" "0.0eta2.0" "00eta20"


    #source run_printPlots_jetPlots.sh "4" "04.12.19_outputCondor" "04-12-19" "0.0eta3.0_noMinBias"   "04.12.19_outputCondor" "04-12-19" "0.0eta3.0" "00eta30_noMB"
    #source run_printPlots_jetPlots.sh "4" "04.12.19_outputCondor" "04-12-19" "0.0eta3.0"   "04.12.19_outputCondor" "04-12-19" "0.0eta3.0" "00eta30_wMB"
    
    #source run_printPlots_jetPlots.sh "4" "04.12.19_outputCondor" "04-12-19" "0.0eta3.0_wJetID_noMinBias"   "04.12.19_outputCondor" "04-12-19" "0.0eta3.0_wJetID"   "00eta30_wJetID_noMB"
    #source run_printPlots_jetPlots.sh "4" "04.12.19_outputCondor" "04-12-19" "0.0eta3.0_wJetID"             "04.12.19_outputCondor" "04-12-19" "0.0eta3.0_wJetID"   "00eta30_wJetID_wMB"

    #source run_printPlots_jetPlots.sh "4" "04.16.19_outputCondor" "04-16-19" "0.0eta3.0_tupel"    "04.16.19_outputCondor" "04-16-19" "0.0eta3.0_wJetID"   "00eta30_tupel_wMB_noJetID" 
    #scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_semiOffPy8_00eta30_tupel_wMB_noJetID_jetPlots.pdf
    #source run_printPlots_jetPlots.sh "4" "04.16.19_outputCondor" "04-16-19" "0.0eta3.0_tupel"    "04.16.19_outputCondor" "04-16-19" "0.0eta3.0_wJetID"   "00eta30_tupel_noMB_noJetID"
    #scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_semiOffPy8_00eta30_tupel_noMB_noJetID_jetPlots.pdf
    #source run_printPlots_jetPlots.sh "4" "04.16.19_outputCondor" "04-16-19" "0.0eta3.0_tupel"    "04.16.19_outputCondor" "04-16-19" "0.0eta3.0_wJetID"   "00eta30_tupel_wMB_wJetID" 
    #scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_semiOffPy8_00eta30_tupel_wMB_wJetID_jetPlots.pdf
    #source run_printPlots_jetPlots.sh "4" "04.16.19_outputCondor" "04-16-19" "0.0eta3.0_tupel"    "04.16.19_outputCondor" "04-16-19" "0.0eta3.0_wJetID"   "00eta30_tupel_noMB_wJetID"
    #scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_semiOffPy8_00eta30_tupel_noMB_wJetID_jetPlots.pdf
    source run_printPlots_jetPlots.sh "4" "06.07.19_outputCondor" "06-07-19" "0.0eta2.0_SMPbins_withJECsys"    "06.07.19_outputCondor" "06-07-19" "0.0eta2.0"   "00eta20_noMB_wJetID"
    


    
    ##############################################
    
    echo ""
    echo "done running jetPlots!"
    echo ""
    
    #scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_semiOffPy8_00eta20_wJetID_jetPlots.\*
    #scp2Serin output/\*jetPlots\*.\*
#    scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_semiOffPy8_00eta30_wMB_jetPlots.pdf
#    scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_semiOffPy8_00eta30_noMB_jetPlots.pdf
    #scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_semiOffPy8_00eta30_wJetID_noMB_jetPlots.pdf
    #scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_semiOffPy8_00eta30_wJetID_wMB_jetPlots.pdf

    ###############################################
    
    #return    
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
    
    #source run_printPlots_jetTrig.sh "4" "10.11.17_outputCondor" "10-11-17" "0.0eta2.0_rec56_HLTCalo" "00eta20_APSDNP_10.17.17_HLTCalo"
    #source run_printPlots_jetTrig.sh "4" "10.11.17_outputCondor" "10-11-17" "0.0eta2.0_rec56_HLTPF" "00eta20_APSDNP_10.17.17_HLTPF"  #goodHLTPFthresh
    
    #source run_printPlots_jetTrig.sh "4" "07.02.18_outputCondor" "07-02-18" "0.0eta2.0_L2L3recoJetID_HLTPF" "00eta20_07.02.18"    #badHLTPFthresh
    #source run_printPlots_jetTrig.sh "4" "10.29.18_outputCondor" "10-29-18" "0.0eta2.0" "00eta20_10.29.18"  #badHLTPFthresh    
    #source run_printPlots_jetTrig.sh "4" "11.06.18_outputCondor" "11-06-18" "0.0eta2.0_wjetid" "00eta20_11.06.18"  #badHLTPFthresh    
    #source run_printPlots_jetTrig.sh "4" "11.06.18_outputCondor" "11-06-18" "0.0eta2.0_wjetid_badHLTPFthresh" "00eta20_11.06.18_bad"  #badHLTPFthresh
    #source run_printPlots_jetTrig.sh "4" "11.06.18_outputCondor" "11-06-18" "0.0eta2.0_wjetid_goodHLTPFthresh" "00eta20_11.06.18_good"  #goodHLTPFthresh
    
    #source run_printPlots_jetTrig.sh "4" "11.06.18_outputCondor" "11-06-18" "0.0eta0.5_wjetid" "00eta05"
    #source run_printPlots_jetTrig.sh "4" "11.06.18_outputCondor" "11-06-18" "0.5eta1.0_wjetid" "05eta10"
    #source run_printPlots_jetTrig.sh "4" "11.06.18_outputCondor" "11-06-18" "1.0eta1.5_wjetid" "10eta15"
    #source run_printPlots_jetTrig.sh "4" "11.06.18_outputCondor" "11-06-18" "1.5eta2.0_wjetid" "15eta20"    
    #source run_printPlots_jetTrig.sh "4" "11.06.18_outputCondor" "11-06-18" "0.0eta2.0_wjetid" "00eta20"
    
    #source run_printPlots_jetTrig.sh "4" "11.14.18_outputCondor" "11-14-18" "0.0eta2.0_trgPtNoRes" "00eta20_11.14.18_trgPtNoRes"  #goodHLTPFthresh
    #source run_printPlots_jetTrig.sh "4" "11.14.18_outputCondor" "11-14-18" "0.0eta2.0_trgPtNoRes_NoJetID" "00eta20_11.14.18_trgPtNoRes_NoJetID"  #goodHLTPFthresh
    #source run_printPlots_jetTrig.sh "4" "11.14.18_outputCondor" "11-14-18" "0.0eta2.0_trgPtL2L3Res" "00eta20_11.14.18_trgPtL2L3Res"  #goodHLTPFthresh
    
    #source run_printPlots_jetTrig.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta2.0" "00eta20_01.06.19_HLT60Thresh75"  #john's recc, HLT60 turn on slow #OLD

    #source run_printPlots_jetTrig.sh "4" "01.31.19_outputCondor" "01-31-19" "0.0eta3.0" "00eta30_01.31.19_HLT60Thresh75"  #john's recc, HLT60 turn on slow #NEW W |Y|<3.0
    #source run_printPlots_jetTrig.sh "4" "01.28.19_outputCondor" "01-28-19" "0.0eta3.0" "00eta30_01.28.19_HLT60Thresh80"  #my guess, looks slightly better
    #source run_printPlots_jetTrig.sh "4" "01.29.19_outputCondor" "01-29-19" "0.0eta3.0" "00eta30_01.29.19_HLT60Thresh85"  #my guess, looks 
    #source run_printPlots_jetTrig.sh "4" "01.29.19_outputCondor" "01-29-19" "_0.0eta3.0_1" "00eta30_01.29.19_HLT60Thresh90"  #my guess, looks 
#    source run_printPlots_jetTrig.sh "4" "01.30.19_outputCondor" "01-30-19" "0.0eta3.0" "00eta30_01.30.19_HLT60Thresh90HLT80Thresh110"  #my guess, looks 
 
    #source run_printPlots_jetTrig.sh "4" "02.07.19_outputCondor" "02-07-19" "0.0eta3.0" "00eta30_02.07.19_HLT60Thresh90HLT80Thresh110"  #my guess, looks     
 #   source run_printPlots_jetTrig.sh "4" "04.12.19_outputCondor" "04-12-19" "0.0eta3.0" "00eta30_04.12.19_wMB" 
 #   source run_printPlots_jetTrig.sh "4" "04.12.19_outputCondor" "04-12-19" "0.0eta3.0_noMinBias" "00eta30_04.12.19_noMB" 

#    source run_printPlots_jetTrig.sh "4" "04.12.19_outputCondor" "04-12-19" "0.0eta3.0_wJetID"           "00eta30_04.12.19_wJetID_wMB"
#    source run_printPlots_jetTrig.sh "4" "04.12.19_outputCondor" "04-12-19" "0.0eta3.0_wJetID_noMinBias" "00eta30_04.12.19_wJetID_noMB"  




    #source run_printPlots_jetTrig.sh "4" "04.16.19_outputCondor" "04-16-19" "0.0eta3.0_tupel"           "00eta30_04.16.19_wJetID_tupel_wMB"
    #source run_printPlots_jetTrig.sh "4" "04.16.19_outputCondor" "04-16-19" "0.0eta3.0_tupel"           "00eta30_04.16.19_wJetID_tupel_noMB"
    #source run_printPlots_jetTrig.sh "4" "04.16.19_outputCondor" "04-16-19" "0.0eta3.0_incjetana"       "00eta30_04.16.19_wJetID_incjetana_wMB"

    #source run_printPlots_jetTrig.sh "4" "04.16.19_outputCondor" "04-16-19" "0.0eta3.0_incjetana"       "00eta30_04.16.19_wJetID_incjetana_noMB"
    #source run_printPlots_jetTrig.sh "4" "04.16.19_outputCondor" "04-16-19" "0.0eta3.0_incjetana"       "00eta30_04.16.19_wJetID_incjetana_noMB_jet80Only"
    #source run_printPlots_jetTrig.sh "4" "04.16.19_outputCondor" "04-16-19" "0.0eta3.0_incjetana"       "00eta30_04.16.19_wJetID_incjetana_noMB_lowjetsOnly"


    #source run_printPlots_jetTrig.sh "4" "04.18.19_outputCondor" "04-18-19" "0.0eta3.0"       "00eta30_04.18.19_wJetID_incjetana_noMB"
    #source run_printPlots_jetTrig.sh "4" "04.18.19_outputCondor" "04-18-19" "0.0eta3.0"       "00eta30_04.18.19_wJetID_incjetana_wMB_wExclusionIssue"

    #scp2Serin "/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta30_04.16.19_wJetID_*_*MB_jetTrig.pdf"
    #scp2Serin "/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta30_04.16.19_wJetID_incjetana_noMB_jet80Only_jetTrig.pdf"    
    #scp2Serin "/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta30_04.18.19_wJetID_incjetana_*MB_jetTrig.*"    
#    scp2Serin "/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta30_04.18.19_wJetID_incjetana_wMB_noMissedHLT40_jetTrig_pt55.*"
    #scp2Serin "/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta30_04.18.19_wJetID_incjetana_wMB_noMissedHLT40_pt??_jetTrig.*"



    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_noL1JetsMB"        "00eta30_04.19.19_wJetID_incjetana_wMB_nolepveto_noL1JetsMB"
    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_noL1JetsMB"       "00eta30_04.19.19_wJetID_incjetana_noMB_nolepveto_noL1JetsMB"
    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_noL1JetsMB"       "00eta30_04.19.19_wJetID_incjetana_wMB_nolepveto_noL1JetsMB_noJetID"

    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_lepveto"        "00eta30_04.19.19_wJetID_incjetana_wMB_lepveto_noL1JetsMB"
    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_lepveto"       "00eta30_04.19.19_wJetID_incjetana_noMB_lepveto_noL1JetsMB"
    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_lepveto"       "00eta30_04.19.19_wJetID_incjetana_wMB_lepveto_noL1JetsMB_noJetID"
    
    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_lepveto_MBwHLT40L1seed"        "00eta30_04.19.19_wJetID_incjetana_wMB_lepveto_MBwHLT40L1seed"
    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_lepveto_MBwHLT40L1seed"       "00eta30_04.19.19_wJetID_incjetana_noMB_lepveto_MBwHLT40L1seed"
    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_lepveto_MBwHLT40L1seed"       "00eta30_04.19.19_wJetID_incjetana_wMB_lepveto_MBwHLT40L1seed_noJetID"

    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_lepveto_MBwHLT4060L1seed"        "00eta30_04.19.19_wJetID_incjetana_wMB_lepveto_MBwHLT4060L1seed"
    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_lepveto_MBwHLT4060L1seed"       "00eta30_04.19.19_wJetID_incjetana_noMB_lepveto_MBwHLT4060L1seed"
    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_lepveto_MBwHLT4060L1seed"       "00eta30_04.19.19_wJetID_incjetana_wMB_lepveto_MBwHLT4060L1seed_noJetID"

    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_lepveto_MBwHLT406080L1seed"        "00eta30_04.19.19_wJetID_incjetana_wMB_lepveto_MBwHLT406080L1seed"
    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_lepveto_MBwHLT406080L1seed"       "00eta30_04.19.19_wJetID_incjetana_noMB_lepveto_MBwHLT406080L1seed"
    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_lepveto_MBwHLT406080L1seed"       "00eta30_04.19.19_wJetID_incjetana_wMB_lepveto_MBwHLT406080L1seed_noJetID"

    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_lepveto_MBwAllL1seed"        "00eta30_04.19.19_wJetID_incjetana_wMB_lepveto_MBwAllL1seed"
    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_lepveto_MBwAllL1seed"       "00eta30_04.19.19_wJetID_incjetana_noMB_lepveto_MBwAllL1seed"
    #source run_printPlots_jetTrig.sh "4" "04.19.19_outputCondor" "04-19-19" "0.0eta3.0_lepveto_MBwAllL1seed"       "00eta30_04.19.19_wJetID_incjetana_wMB_lepveto_MBwAllL1seed_noJetID"

    #scp2Serin "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.18.19_wJetID_incjetana_wMB_wExclusionIssue_jetTrig.*"
    #scp2Serin "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.19.19_wJetID_incjetana_wMB_nolepveto_noL1JetsMB_jetTrig.*"
    #scp2Serin "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.19.19_wJetID_incjetana_wMB_lepveto_noL1JetsMB_jetTrig.*"
    #scp2Serin "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.19.19_wJetID_incjetana_wMB_lepveto_MBwHLT40L1seed_jetTrig.*"
    #scp2Serin "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.19.19_wJetID_incjetana_wMB_lepveto_MBwHLT4060L1seed_jetTrig.*"
    #scp2Serin "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.19.19_wJetID_incjetana_wMB_lepveto_MBwHLT406080L1seed_jetTrig.*"
    #scp2Serin "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.19.19_wJetID_incjetana_wMB_lepveto_MBwAllL1seed_jetTrig.*"
    #scp2Serin "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.19.19_wJetID_incjetana_wMB_lepveto_*_jetTrig.*"
    #scp2Serin "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.19.19_wJetID_incjetana_wMB_nolepveto_noL1JetsMB_noJetID_jetTrig.*"
    #scp2Serin "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.19.19_wJetID_incjetana_wMB_lepveto_*_noJetID_jetTrig.*"

    ### SKIPS HLTPF AND HLTCALO EVTS FOR MB DECISION
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBwAllL1seed"       "00eta30_04.23.19_wJetID_incjetana_noMB_lepveto_MBwAllL1seed"
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBwAllL1seed"        "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_MBwAllL1seed"
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBwAllL1seed"        "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_noL1Jet68"
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBwAllL1seed"        "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_noL1Jet60"
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBwAllL1seed"        "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_noL1Jet52"
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBwAllL1seed"        "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_noL1Jet48"
    
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_MBwAllL1seed_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_noL1Jet68_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_noL1Jet60_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_noL1Jet52_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_noL1Jet48_jetTrig.*"
    

    ### SKIPS ONLY HLTPF EVTS FOR MB DECISION

    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBwAllL1seed_skipHLTPFOnly"      "00eta30_04.23.19_wJetID_incjetana_noMB_lepveto_skipHLTPFOnly_MBwAllL1seed"
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBwAllL1seed_skipHLTPFOnly"       "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_skipHLTPFOnly_MBwAllL1seed"
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBwAllL1seed_skipHLTPFOnly"       "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_skipHLTPFOnly_noL1Jet68"
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBwAllL1seed_skipHLTPFOnly"       "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_skipHLTPFOnly_noL1Jet60"
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBwAllL1seed_skipHLTPFOnly"       "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_skipHLTPFOnly_noL1Jet52"
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBwAllL1seed_skipHLTPFOnly"       "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_skipHLTPFOnly_noL1Jet48"

    
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_skipHLTPFOnly_MBwAllL1seed_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_skipHLTPFOnly_noL1Jet68_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_skipHLTPFOnly_noL1Jet60_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_skipHLTPFOnly_noL1Jet52_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_skipHLTPFOnly_noL1Jet48_jetTrig.*"



    ### SKIPS NO EVTS FOR MB DECISION

    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBallHLTPF"      "00eta30_04.23.19_wJetID_incjetana_noMB_lepveto_MBallHLTPF"
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBallHLTPF"      "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_MBallHLTPF"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_MBallHLTPF_jetTrig.*"

    ##subtract HLTJet && isMB one by one
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBallHLTPF"       "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_noHLTJet80"
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBallHLTPF"       "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_noHLTJet60"
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBallHLTPF"       "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_noHLTJet40"

    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_noHLTJet80_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_noHLTJet60_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_noHLTJet40_jetTrig.*"

    ##subtract L1Jet && isMBt one by one
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBallHLTPF"       "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_allHLT_noL1Jet68"
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBallHLTPF"       "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_allHLT_noL1Jet60"
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBallHLTPF"       "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_allHLT_noL1Jet52"
    #source run_printPlots_jetTrig.sh "4" "04.23.19_outputCondor" "04-23-19" "0.0eta3.0_lepveto_MBallHLTPF"       "00eta30_04.23.19_wJetID_incjetana_wMB_lepveto_allHLT_noL1Jet48"

    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_allHLT_noL1Jet68_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_allHLT_noL1Jet60_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_allHLT_noL1Jet52_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.23.19_wJetID_incjetana_*MB_lepveto_allHLT_noL1Jet48_jetTrig.*"



    ### SIMPLE TRIGPT EXCLUSION W/ HLT40's LOWER THRESHOLD (IF NO TRIGGER JET OBJECT FROM OTHERS, TRGPT==0)
    #source run_printPlots_jetTrig.sh "4" "04.24.19_outputCondor" "04-24-19" "0.0eta3.0_lepveto_MBexc_wTrgPt"      "00eta30_04.24.19_wJetID_incjetana_noMB_lepveto_MBexc_wTrgPt"
    #source run_printPlots_jetTrig.sh "4" "04.24.19_outputCondor" "04-24-19" "0.0eta3.0_lepveto_MBexc_wTrgPt"       "00eta30_04.24.19_wJetID_incjetana_wMB_lepveto_MBexc_wTrgPt"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.24.19_wJetID_incjetana_*MB_lepveto_MBexc_wTrgPt_jetTrig.*"

    ##subtract HLTJet && isMB one by one
    #source run_printPlots_jetTrig.sh "4" "04.24.19_outputCondor" "04-24-19" "0.0eta3.0_lepveto_MBexc_wTrgPt"       "00eta30_04.24.19_wJetID_incjetana_wMB_lepveto_MBexc_wTrgPt_noHLTJet80"
    #source run_printPlots_jetTrig.sh "4" "04.24.19_outputCondor" "04-24-19" "0.0eta3.0_lepveto_MBexc_wTrgPt"       "00eta30_04.24.19_wJetID_incjetana_wMB_lepveto_MBexc_wTrgPt_noHLTJet60"
    #source run_printPlots_jetTrig.sh "4" "04.24.19_outputCondor" "04-24-19" "0.0eta3.0_lepveto_MBexc_wTrgPt"       "00eta30_04.24.19_wJetID_incjetana_wMB_lepveto_MBexc_wTrgPt_noHLTJet40"

    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.24.19_wJetID_incjetana_*MB_lepveto_MBexc_wTrgPt_noHLTJet80_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.24.19_wJetID_incjetana_*MB_lepveto_MBexc_wTrgPt_noHLTJet60_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.24.19_wJetID_incjetana_*MB_lepveto_MBexc_wTrgPt_noHLTJet40_jetTrig.*"

    ##subtract L1Jet && isMB one by one
    #source run_printPlots_jetTrig.sh "4" "04.24.19_outputCondor" "04-24-19" "0.0eta3.0_lepveto_MBexc_wTrgPt"       "00eta30_04.24.19_wJetID_incjetana_wMB_lepveto_MBexc_wTrgPt_noL1Jet68"
    #source run_printPlots_jetTrig.sh "4" "04.24.19_outputCondor" "04-24-19" "0.0eta3.0_lepveto_MBexc_wTrgPt"       "00eta30_04.24.19_wJetID_incjetana_wMB_lepveto_MBexc_wTrgPt_noL1Jet60"
    #source run_printPlots_jetTrig.sh "4" "04.24.19_outputCondor" "04-24-19" "0.0eta3.0_lepveto_MBexc_wTrgPt"       "00eta30_04.24.19_wJetID_incjetana_wMB_lepveto_MBexc_wTrgPt_noL1Jet52"

    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.24.19_wJetID_incjetana_*MB_lepveto_MBexc_wTrgPt_noL1Jet68_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.24.19_wJetID_incjetana_*MB_lepveto_MBexc_wTrgPt_noL1Jet60_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.24.19_wJetID_incjetana_*MB_lepveto_MBexc_wTrgPt_noL1Jet52_jetTrig.*"


    ### SIMPLE TRIGPT EXCLUSION W/ HLT40's LOWER THRESHOLD (IF NO TRIGGER JET OBJECT FROM OTHERS, TRGPT==0), W/ is/HLT40 && NOT HLTMB plots, w/ isMB && NOT L1_HF_COINCIDENCE plots
    #source run_printPlots_jetTrig.sh "4" "04.25.19_outputCondor" "04-25-19" "0.0eta3.0_lepveto_MBexc_wTrgPt"      "00eta30_04.25.19_wJetID_incjetana_noMB_lepveto_MBexc_wTrgPt"
    #source run_printPlots_jetTrig.sh "4" "04.25.19_outputCondor" "04-25-19" "0.0eta3.0_lepveto_MBexc_wTrgPt"      "00eta30_04.25.19_wJetID_incjetana_wMB_lepveto_MBexc_wTrgPt"
    #source run_printPlots_jetTrig.sh "4" "04.25.19_outputCondor" "04-25-19" "0.0eta3.0_lepveto_MBexc_wTrgPt"      "00eta30_04.25.19_wJetID_incjetana_wMB_lepveto_MBexc_wTrgPt_wL1HFcoinc"


    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.25.19_wJetID_incjetana_*MB_lepveto_MBexc_wTrgPt_jetTrig.*"
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.25.19_wJetID_incjetana_*MB_lepveto_MBexc_wTrgPt_wL1HFcoinc_jetTrig.*"



    ### 
    #source run_printPlots_jetTrig.sh "4" "04.29.19_outputCondor" "04-29-19" "0.0eta3.0_L1MBforexclu"      "00eta30_04.29.19_wJetID_incjetana_noMB_isMBexc_wTrgPt_L1MBexclu"
    #source run_printPlots_jetTrig.sh "4" "04.29.19_outputCondor" "04-29-19" "0.0eta3.0_L1MBforexclu"       "00eta30_04.29.19_wJetID_incjetana_wMB_isMBexc_wTrgPt_L1MBexclu"
    #
    #source run_printPlots_jetTrig.sh "4" "04.29.19_outputCondor" "04-29-19" "0.0eta3.0_HLTMBforexclu"      "00eta30_04.29.19_wJetID_incjetana_noMB_isMBexc_wTrgPt_HLTMBexclu"
    #source run_printPlots_jetTrig.sh "4" "04.29.19_outputCondor" "04-29-19" "0.0eta3.0_HLTMBforexclu"       "00eta30_04.29.19_wJetID_incjetana_wMB_isMBexc_wTrgPt_HLTMBexclu"
    #
    #source run_printPlots_jetTrig.sh "4" "04.29.19_outputCondor" "04-29-19" "0.0eta3.0_HLTMBforexclu_johnHLTthresh"      "00eta30_04.29.19_wJetID_incjetana_noMB_isMBexc_wTrgPt_HLTMBexclu_johnHLTthresh"
    #source run_printPlots_jetTrig.sh "4" "04.29.19_outputCondor" "04-29-19" "0.0eta3.0_HLTMBforexclu_johnHLTthresh"       "00eta30_04.29.19_wJetID_incjetana_wMB_isMBexc_wTrgPt_HLTMBexclu_johnHLTthresh"
    #
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_04.29.19_wJetID_incjetana_*MB_isMBexc_wTrgPt_*MBexclu*_jetTrig.*"


    
    ###
    #source run_printPlots_jetTrig.sh "4" "05.01.19_outputCondor" "05-01-19" "0.0eta3.0_isHLTforDupeSkip"      "00eta30_05.01.19_wJetID_incjetana_noMB_isHLTforDupeSkip_wTrgPt_HLTMBexclu"
    #source run_printPlots_jetTrig.sh "4" "05.01.19_outputCondor" "05-01-19" "0.0eta3.0_isHLTforDupeSkip"       "00eta30_05.01.19_wJetID_incjetana_wMB_isHLTforDupeSkip_wTrgPt_HLTMBexclu"
    
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_05.01.19_wJetID_incjetana_*MB_isHLTforDupeSkip_wTrgPt_HLTMBexclu_jetTrig.*"

    ##
    source run_printPlots_jetTrig.sh "4" "06.07.19_outputCondor" "06-07-19" "0.0eta2.0_SMPbins_withJECsys"      "00eta20_06.07.19_wJetID_incjetana_noMB_SMPbins_withJECsys"
    #source run_printPlots_jetTrig.sh "4" "05.01.19_outputCondor" "05-01-19" "0.0eta3.0_isHLTforDupeSkip"       "00eta30_05.01.19_wJetID_incjetana_wMB_isHLTforDupeSkip_wTrgPt_HLTMBexclu"
    
    #scp2Serin  "${JETTRIGOUTDIR}/ak4PF_HPtJetTrig_00eta30_05.01.19_wJetID_incjetana_*MB_isHLTforDupeSkip_wTrgPt_HLTMBexclu_jetTrig.*"


    ##############################################
    
    echo ""
    echo "done running jetTrig!"
    echo ""

    return
    

    

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
    
    #source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "0.0eta0.5_wjetid"  "0.0eta2.0_NOjetid" "00eta05" "ppData"
    #source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "0.5eta1.0_wjetid"  "0.0eta2.0_NOjetid" "05eta10" "ppData"
    #source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "1.0eta1.5_wjetid"  "0.0eta2.0_NOjetid" "10eta15" "ppData"
    #source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "1.5eta2.0_wjetid"  "0.0eta2.0_NOjetid" "15eta20" "ppData"
    #source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "0.0eta0.5_wjetid"  "0.0eta2.0_NOjetid" "00eta05" "ppMC"
    #source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "0.5eta1.0_wjetid"  "0.0eta2.0_NOjetid" "05eta10" "ppMC"
    #source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "1.0eta1.5_wjetid"  "0.0eta2.0_NOjetid" "10eta15" "ppMC"
    #source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "1.5eta2.0_wjetid"  "0.0eta2.0_NOjetid" "15eta20" "ppMC"

    #source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "0.0eta2.0_wjetid"  "0.0eta2.0_NOjetid" "00eta20" "ppData"
    #source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "0.0eta2.0_wjetid"  "0.0eta2.0_NOjetid" "00eta20" "ppMC"
    
    
    #source run_printPlots_jetIDPlots.sh "4" "11.14.18_outputCondor" "11-14-18"  "0.0eta2.0_trgPtNoRes"  "0.0eta2.0_trgPtNoRes_NoJetID" "00eta20_trgPtNoRes" "ppData"
    #source run_printPlots_jetIDPlots.sh "4" "11.14.18_outputCondor" "11-14-18"  "0.0eta2.0_trgPtL2L3Res"  "0.0eta2.0_trgPtL2L3Res_NoJetID" "00eta20_trgPtL2L3Res" "ppData"
    
    #source run_printPlots_jetIDPlots.sh "4" "11.14.18_outputCondor" "11-14-18"  "0.0eta2.0"  "0.0eta2.0_NoJetID" "00eta20" "ppMC"
    
    source run_printPlots_jetIDPlots_v2.sh "4" "04.16.19_outputCondor" "04-16-19"  "0.0eta3.0_tupel"  "00eta30_tupel_wMB" "ppData"
    source run_printPlots_jetIDPlots_v2.sh "4" "04.16.19_outputCondor" "04-16-19"  "0.0eta3.0_tupel"  "00eta30_tupel_noMB" "ppData"
    
    #source run_printPlots_jetIDPlots_v2.sh "4" "04.16.19_outputCondor" "04-16-19"  "0.0eta3.0_incjetana"  "00eta30_incjetana_wMB" "ppData"
    #source run_printPlots_jetIDPlots_v2.sh "4" "04.16.19_outputCondor" "04-16-19"  "0.0eta3.0_incjetana"  "00eta30_incjetana_noMB" "ppData"
    
    
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



