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
    source run_printPlots_jetPlots.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta3.0"   "11.27.18_outputCondor" "11-27-18" "0.0eta2.0" "00eta20"

    
    ##############################################
    
    echo ""
    echo "done running jetPlots!"
    echo ""
    
    #scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_semiOffPy8_00eta20_wJetID_jetPlots.\*
    #scp2Serin output/\*jetPlots\*.\*
    
    ###############################################
    
    #return    
fi




if [[ $dojetTrigQAPlots -eq 1 ]]
then
    echo ""
    echo "compiling printPlots_jetTrig.C"
    echo ""
    
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
 
    source run_printPlots_jetTrig.sh "4" "02.07.19_outputCondor" "02-07-19" "0.0eta3.0" "00eta30_02.07.19_HLT60Thresh90HLT80Thresh110"  #my guess, looks     
    
    ##############################################
    
    echo ""
    echo "done running jetTrig!"
    echo ""
    
    #scp2Serin output/\*jetTrig\*
    #scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta20_11.06.18_good_jetTrig.pdf
    #scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta20_11.14.18_trgPt\*.pdf
    #scp2Serin "/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta?0_01.??.19_HLT60Thresh??_jetTrig.pdf"
    ##############################################
    
    #return
fi




if [[ $dojetIDQAPlots -eq 1 ]]
then
    
    echo ""
    echo "compiling printPlots_jetIDPlots.C"
    echo ""
    
    rootcompile printPlots_jetIDPlots.C

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
    
    
    source run_printPlots_jetIDPlots.sh "4" "11.14.18_outputCondor" "11-14-18"  "0.0eta2.0_trgPtNoRes"  "0.0eta2.0_trgPtNoRes_NoJetID" "00eta20_trgPtNoRes" "ppData"
    #source run_printPlots_jetIDPlots.sh "4" "11.14.18_outputCondor" "11-14-18"  "0.0eta2.0_trgPtL2L3Res"  "0.0eta2.0_trgPtL2L3Res_NoJetID" "00eta20_trgPtL2L3Res" "ppData"
    
    source run_printPlots_jetIDPlots.sh "4" "11.14.18_outputCondor" "11-14-18"  "0.0eta2.0"  "0.0eta2.0_NoJetID" "00eta20" "ppMC"
    
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










