#!/bin/bash


dojetIDQAPlots=0
dojetQAPlots=0
dojetTrigQAPlots=1


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
    
    #source run_printPlots_jetIDPlots.sh "4" "10.29.18_outputCondor" "10-29-18"  "0.0eta2.0"  "_0.0eta2.0_1" "00eta20" "ppData"
    #source run_printPlots_jetIDPlots.sh "4" "10.29.18_outputCondor" "10-29-18"  "0.0eta2.0"  "_0.0eta2.0_1" "00eta20" "ppMC"
    
    source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "0.0eta0.5_wjetid"  "0.0eta2.0_NOjetid" "00eta05" "ppData"
    source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "0.5eta1.0_wjetid"  "0.0eta2.0_NOjetid" "05eta10" "ppData"
    source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "1.0eta1.5_wjetid"  "0.0eta2.0_NOjetid" "10eta15" "ppData"
    source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "1.5eta2.0_wjetid"  "0.0eta2.0_NOjetid" "15eta20" "ppData"

    source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "0.0eta2.0_wjetid"  "0.0eta2.0_NOjetid" "00eta20" "ppData"

    source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "0.0eta0.5_wjetid"  "0.0eta2.0_NOjetid" "00eta05" "ppMC"
    source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "0.5eta1.0_wjetid"  "0.0eta2.0_NOjetid" "05eta10" "ppMC"
    source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "1.0eta1.5_wjetid"  "0.0eta2.0_NOjetid" "10eta15" "ppMC"
    source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "1.5eta2.0_wjetid"  "0.0eta2.0_NOjetid" "15eta20" "ppMC"

    source run_printPlots_jetIDPlots.sh "4" "11.06.18_outputCondor" "11-06-18"  "0.0eta2.0_wjetid"  "0.0eta2.0_NOjetid" "00eta20" "ppMC"
    
    ##############################################
    
    echo ""
    echo "done running jetIDPlots!"
    echo ""

    scp2Serin output/\*.\*

    ##############################################
    
    return
fi






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
    
    #source run_printPlots_jetPlots.sh "4" "10.29.18_outputCondor" "10-29-18" "0.0eta2.0" "10.29.18_outputCondor" "10-29-18" "0.0eta2.0" "00eta20_wJetID"
    
    source run_printPlots_jetPlots.sh "4" "11.06.18_outputCondor" "11-06-18" "0.0eta0.5_wjetid" "11.06.18_outputCondor" "11-06-18" "0.0eta0.5_wjetid" "00eta05_wJetID"
    source run_printPlots_jetPlots.sh "4" "11.06.18_outputCondor" "11-06-18" "0.5eta1.0_wjetid" "11.06.18_outputCondor" "11-06-18" "0.5eta1.0_wjetid" "05eta10_wJetID"
    source run_printPlots_jetPlots.sh "4" "11.06.18_outputCondor" "11-06-18" "1.0eta1.5_wjetid" "11.06.18_outputCondor" "11-06-18" "1.0eta1.5_wjetid" "10eta15_wJetID"
    source run_printPlots_jetPlots.sh "4" "11.06.18_outputCondor" "11-06-18" "1.5eta2.0_wjetid" "11.06.18_outputCondor" "11-06-18" "1.5eta2.0_wjetid" "15eta20_wJetID"
    
    source run_printPlots_jetPlots.sh "4" "11.06.18_outputCondor" "11-06-18" "0.0eta2.0_wjetid" "11.06.18_outputCondor" "11-06-18" "0.0eta2.0_wjetid" "00eta20_wJetID"
    
    ##############################################
    
    echo ""
    echo "done running jetPlots!"
    echo ""
    
    #scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_semiOffPy8_00eta20_wJetID_jetPlots.\*
    scp2Serin output/\*jetPlots\*.\*
    
    ###############################################
    
    return    
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
    source run_printPlots_jetTrig.sh "4" "11.06.18_outputCondor" "11-06-18" "0.0eta2.0_wjetid_goodHLTPFthresh" "00eta20_11.06.18_good"  #goodHLTPFthresh

    #source run_printPlots_jetTrig.sh "4" "11.06.18_outputCondor" "11-06-18" "0.0eta0.5_wjetid" "00eta05"
    #source run_printPlots_jetTrig.sh "4" "11.06.18_outputCondor" "11-06-18" "0.5eta1.0_wjetid" "05eta10"
    #source run_printPlots_jetTrig.sh "4" "11.06.18_outputCondor" "11-06-18" "1.0eta1.5_wjetid" "10eta15"
    #source run_printPlots_jetTrig.sh "4" "11.06.18_outputCondor" "11-06-18" "1.5eta2.0_wjetid" "15eta20"    
    #source run_printPlots_jetTrig.sh "4" "11.06.18_outputCondor" "11-06-18" "0.0eta2.0_wjetid" "00eta20"
    
    ##############################################
    
    echo ""
    echo "done running jetTrig!"
    echo ""
    
    #scp2Serin output/\*jetTrig\*
    scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_00eta20_11.06.18_good_jetTrig.pdf
    ##############################################
    
    return
fi






return






























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

return
###############################################


















#source run_printPlots_jetTrig.sh "4" "8.28.17_outputCondor" "08-28-17" "0.0eta4.7" "00eta47_8.28.17"

#source run_printPlots_jetTrig.sh "4" "8.17-22.17_outputCondor/genpt49"  "08-28-17" "0.0eta0.5" "00eta05_8.23.17_genpt49"
#source run_printPlots_jetTrig.sh "4" "8.17-22.17_outputCondor/genpt49"  "08-28-17" "0.5eta1.0" "05eta10_8.23.17_genpt49"
#source run_printPlots_jetTrig.sh "4" "8.17-22.17_outputCondor/genpt49"  "08-28-17" "1.0eta1.5" "10eta15_8.23.17_genpt49"
#source run_printPlots_jetTrig.sh "4" "8.17-22.17_outputCondor/genpt49"  "08-28-17" "1.5eta2.0" "15eta20_8.23.17_genpt49"
#source run_printPlots_jetTrig.sh "4" "8.17-22.17_outputCondor/genpt49"  "08-28-17" "2.0eta2.5" "20eta25_8.23.17_genpt49"
#source run_printPlots_jetTrig.sh "4" "8.17-22.17_outputCondor/genpt49"  "08-28-17" "2.5eta3.0" "25eta30_8.23.17_genpt49"
#source run_printPlots_jetTrig.sh "4" "8.17-22.17_outputCondor/genpt49"  "08-28-17" "3.2eta4.7" "32eta47_8.23.17_genpt49"

#source run_printPlots_jetTrig.sh "4" "08.28.17_outputCondor/ppData_CaloHLT" "08-28-17" "0.0eta2.0" "00eta20_8.28.17_HLTCalo"
#source run_printPlots_jetTrig.sh "4" "08.28.17_outputCondor/ppData_CaloHLT" "08-28-17" "2.0eta3.0" "20eta30_8.28.17_HLTCalo"
#source run_printPlots_jetTrig.sh "4" "08.28.17_outputCondor/ppData_CaloHLT" "08-28-17" "3.2eta4.7" "32eta47_8.28.17_HLTCalo"

#source run_printPlots_jetTrig.sh "4" "08.28.17_outputCondor/ppData_PFHLT" "08-28-17" "0.0eta2.0" "00eta20_8.28.17_HLTPF"
#source run_printPlots_jetTrig.sh "4" "08.28.17_outputCondor/ppData_PFHLT" "08-28-17" "2.0eta3.0" "20eta30_8.28.17_HLTPF"
#source run_printPlots_jetTrig.sh "4" "08.28.17_outputCondor/ppData_PFHLT" "08-28-17" "3.2eta4.7" "32eta47_8.28.17_HLTPF"

#source run_printPlots_jetTrig.sh "4" "09.21.17_outputCondor/HLTCalo" "09-21-17" "0.0eta2.0" "00eta20_9.21.17_HLTCalo"
