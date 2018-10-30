#!/bin/bash

echo ""
echo "compiling printPlots_jetIDPlots.C"
echo ""

rootcompile printPlots_jetIDPlots.C
##############################################

echo ""
echo "done compiling. Running!"
echo ""


source run_printPlots_jetIDPlots.sh "4" "10.29.18_outputCondor" "10-29-18"  "0.0eta2.0"  "_0.0eta2.0_1" "00eta20" "ppData"
source run_printPlots_jetIDPlots.sh "4" "10.29.18_outputCondor" "10-29-18"  "0.0eta2.0"  "_0.0eta2.0_1" "00eta20" "ppMC"
#scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_\*.\*

##############################################

echo ""
echo "done running jetIDPlots!"
echo ""





##############################################

#return
##############################################





echo ""
echo "compiling printPlots_jetPlots.C"
echo ""

rootcompile printPlots_jetPlots.C
##############################################

echo ""
echo "done compiling. Running!"
echo ""


source run_printPlots_jetPlots.sh "4" "10.29.18_outputCondor" "10-29-18" "0.0eta2.0" "10.29.18_outputCondor" "10-29-18" "0.0eta2.0" "00eta20_wJetID"

#scp2Serin /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/output/ak4PF_HPtJetTrig_semiOffPy8_00eta20_wJetID_jetPlots.\*
scp2Serin output/\*.\*
##############################################


echo ""
echo "done running jetPlots!"
echo ""
###############################################

return
################################################









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







echo ""
echo "compiling printPlots_jetTrig.C"
echo ""

#rootcompile printPlots_jetTrig.C
##############################################

echo ""
echo "done compiling. Running!"
echo ""

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




#source run_printPlots_jetTrig.sh "4" "10.11.17_outputCondor" "10-11-17" "0.0eta2.0_rec56_HLTCalo" "00eta20_APSDNP_10.17.17"
##############################################

echo ""
echo "done running jetTrig!"
echo ""
##############################################

return
##############################################















