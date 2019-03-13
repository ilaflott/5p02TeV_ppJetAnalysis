#!/bin/bash

#echo ""
#echo "compiling printPlots_ppMC_JERS.C"
#echo ""
#
#rootcompile printPlots_ppMC_JERS.C
#
#echo ""
#echo "done compiling. Running!"
#echo ""

#source run_printPlots_ppMC_JERS.sh "4" "08.07.18_outputCondor" "08-07-18" "0.0eta2.0" "00eta20"

#source run_printPlots_ppMC_JERS.sh "4" "11.06.18_outputCondor" "11-06-18" "0.0eta0.5" "00eta05"
#source run_printPlots_ppMC_JERS.sh "4" "11.06.18_outputCondor" "11-06-18" "0.5eta1.0" "05eta10"
#source run_printPlots_ppMC_JERS.sh "4" "11.06.18_outputCondor" "11-06-18" "1.0eta1.5" "10eta15"
#source run_printPlots_ppMC_JERS.sh "4" "11.06.18_outputCondor" "11-06-18" "1.5eta2.0" "15eta20"


#source run_printPlots_ppMC_JERS.sh "4" "11.14.18_outputCondor" "11-14-18" "0.0eta2.0" "00eta20" "1" "0"         #JER only
#source run_printPlots_ppMC_JERS.sh "4" "11.14.18_outputCondor" "11-14-18" "0.0eta2.0_unf" "00eta20_unf" "0" "1" #MCEff only

## usage
## source run_printPlots_ppMC_JERS.sh <radius> <condordir> <date input> <input etabin tag> <output etabin desc> <doJERplots> <doMCEffplots> <abs eta bin>

#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta2.0" "00eta20" "1" "0" "-1"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta2.0" "00eta05" "1" "0" "0"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta2.0" "05eta10" "1" "0" "1"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta2.0" "10eta15" "1" "0" "2"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta2.0" "15eta20" "1" "0" "3"        #JER 0-2only


#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta3.0" "00eta05" "1" "0" "0"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta3.0" "05eta10" "1" "0" "1"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta3.0" "10eta15" "1" "0" "2"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta3.0" "15eta20" "1" "0" "3"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta3.0" "20eta25" "1" "0" "4"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta3.0" "25eta30" "1" "0" "5"        #JER 0-2only

#echo ""
#echo "done running printPlots_ppMC_JERS"
#echo ""

#echo ""
#echo "scp2Serin"
#echo ""

#scp2Serin output/ak4PF\*


echo ""
echo "compiling printPlots_ppMC_JERS_etabin.C"
echo ""

rootcompile printPlots_ppMC_JERS_etabin.C

echo ""
echo "done compiling. Running!"
echo ""

source run_printPlots_ppMC_JERS_etabin.sh "4" "02.18.19_outputCondor" "02-18-19" "0.0eta3.0_JER" "00eta05_semifinal" "1" "0" "0"         #JER any |y|bin
source run_printPlots_ppMC_JERS_etabin.sh "4" "02.18.19_outputCondor" "02-18-19" "0.0eta3.0_JER" "05eta10_semifinal" "1" "0" "1"        #JER any |y|bin
source run_printPlots_ppMC_JERS_etabin.sh "4" "02.18.19_outputCondor" "02-18-19" "0.0eta3.0_JER" "10eta15_semifinal" "1" "0" "2"        #JER any |y|bin
source run_printPlots_ppMC_JERS_etabin.sh "4" "02.18.19_outputCondor" "02-18-19" "0.0eta3.0_JER" "15eta20_semifinal" "1" "0" "3"        #JER any |y|bin
source run_printPlots_ppMC_JERS_etabin.sh "4" "02.18.19_outputCondor" "02-18-19" "0.0eta3.0_JER" "20eta25_semifinal" "1" "0" "4"        #JER any |y|bin
source run_printPlots_ppMC_JERS_etabin.sh "4" "02.18.19_outputCondor" "02-18-19" "0.0eta3.0_JER" "25eta30_semifinal" "1" "0" "5"        #JER any |y|bin

echo ""
echo "done running printPlots_ppMC_JERS_etabin"
echo ""

#echo ""
#echo "scp2Serin"
#echo ""

scp2Serin "output/ak4PF_PY8JER_??eta??_semifinal.*"

##############################################
return
##############################################

























































##### DEPRECATED

##echo ""
##echo "compiling printPlots_jetPlots.C"
##echo ""
##
##rootcompile printPlots_jetPlots.C
##
##echo ""
##echo "done compiling. Running!"
##echo ""
##
##source run_printPlots_jetPlots.sh 4 "08-15-17" "3.2eta4.7" "32eta47_8.15.17"
##source run_printPlots_jetPlots.sh 4 "08-15-17" "2.0eta3.0" "20eta30_8.15.17"
##source run_printPlots_jetPlots.sh 4 "08-15-17" "0.0eta2.0" "00eta20_8.15.17"
#
#echo ""
#echo "done running jetPlots!"
#echo ""
#
###############################################
#return
###############################################
#
#
#
#
#
#
##echo ""
##echo "compiling printPlots_jetIDPlots.C"
##echo ""
##
##rootcompile printPlots_jetIDPlots.C
##
##echo ""
##echo "done compiling. Running!"
##echo ""
##
##source run_printPlots_jetIDPlots.sh 4 "08-15-17" "3.2eta4.7" "32eta47_8.15.17"
##source run_printPlots_jetIDPlots.sh 4 "08-15-17" "2.0eta3.0" "20eta30_8.15.17"
##source run_printPlots_jetIDPlots.sh 4 "08-15-17" "0.0eta2.0" "00eta20_8.15.17"
#
#echo ""
#echo "done running jetIDPlots!"
#echo ""
###############################################
#return
###############################################
#
#
##echo ""
##echo "compiling printPlots_jetTrig.C"
##echo ""
##
##rootcompile printPlots_jetTrig.C
##
##echo ""
##echo "done compiling. Running!"
##echo ""
##
##source run_printPlots_jetTrig.sh 4 "08-15-17" "3.2eta4.7" "32eta47_8.15.17"
##source run_printPlots_jetTrig.sh 4 "08-15-17" "2.0eta3.0" "20eta30_8.15.17"
##source run_printPlots_jetTrig.sh 4 "08-15-17" "0.0eta2.0" "00eta20_8.15.17"
#
#echo ""
#echo "done running jetTrig!"
#echo ""
###############################################
#return
###############################################









#source run_printPlots_ppMC_JERS.sh "4" "08.16.17_outputCondor" "08-16-17" "0.0eta0.5" "00eta05_9.7.17_genpt32"
#source run_printPlots_ppMC_JERS.sh "4" "08.16.17_outputCondor" "08-16-17" "0.5eta1.0" "05eta10_9.7.17_genpt32"
#source run_printPlots_ppMC_JERS.sh "4" "08.16.17_outputCondor" "08-16-17" "1.0eta1.5" "10eta15_9.7.17_genpt32"
#source run_printPlots_ppMC_JERS.sh "4" "08.16.17_outputCondor" "08-16-17" "1.5eta2.0" "15eta20_9.7.17_genpt32"

#source run_printPlots_ppMC_JERS.sh "4" "09.07.17_outputCondor/recpt56_genpt49_ptbins56" "09-07-17" "0.0eta0.5" "00eta05_9.7.17_gen49rec56bin56"
#source run_printPlots_ppMC_JERS.sh "4" "09.07.17_outputCondor/recpt56_genpt49_ptbins56" "09-07-17" "0.5eta1.0" "05eta10_9.7.17_gen49rec56bin56"
#source run_printPlots_ppMC_JERS.sh "4" "09.07.17_outputCondor/recpt56_genpt49_ptbins56" "09-07-17" "1.0eta1.5" "10eta15_9.7.17_gen49rec56bin56"
#source run_printPlots_ppMC_JERS.sh "4" "09.07.17_outputCondor/recpt56_genpt49_ptbins56" "09-07-17" "1.5eta2.0" "15eta20_9.7.17_gen49rec56bin56"

#source run_printPlots_ppMC_JERS.sh "4" "09.07.17_outputCondor/recpt37_genpt49_ptbins49" "09-07-17" "0.0eta0.5" "00eta05_9.7.17_gen49rec37bin49"
#source run_printPlots_ppMC_JERS.sh "4" "09.07.17_outputCondor/recpt37_genpt49_ptbins49" "09-07-17" "0.5eta1.0" "05eta10_9.7.17_gen49rec37bin49"
#source run_printPlots_ppMC_JERS.sh "4" "09.07.17_outputCondor/recpt37_genpt49_ptbins49" "09-07-17" "1.0eta1.5" "10eta15_9.7.17_gen49rec37bin49"
#source run_printPlots_ppMC_JERS.sh "4" "09.07.17_outputCondor/recpt37_genpt49_ptbins49" "09-07-17" "1.5eta2.0" "15eta20_9.7.17_gen49rec37bin49"

#source run_printPlots_ppMC_JERS.sh "4" "09.13.17_outputCondor/recpt10_genpt10_ptbins32" "09-13-17" "0.0eta0.5" "00eta05_9.13.17_recgenpt_10_10_ptbins32"
#source run_printPlots_ppMC_JERS.sh "4" "09.13.17_outputCondor/recpt43_genpt32_ptbins32" "09-13-17" "0.0eta0.5" "00eta05_9.13.17_recgenpt_43_32_ptbins32"
#source run_printPlots_ppMC_JERS.sh "4" "09.13.17_outputCondor/recpt56_genpt49_ptbins49" "09-13-17" "0.0eta0.5" "00eta05_9.13.17_recgenpt_56_49_ptbins49"

#source run_printPlots_ppMC_JERS.sh "4" "09.13.17_outputCondor/recpt20_genpt20_ptbins32" "09-13-17" "0.0eta0.5" "00eta05_9.13.17_recgenpt_20_20_ptbins32"
#source run_printPlots_ppMC_JERS.sh "4" "09.13.17_outputCondor/recpt20_genpt20_ptbins32" "09-13-17" "0.5eta1.0" "05eta10_9.13.17_recgenpt_20_20_ptbins32"
#source run_printPlots_ppMC_JERS.sh "4" "09.13.17_outputCondor/recpt20_genpt20_ptbins32" "09-13-17" "1.0eta1.5" "10eta15_9.13.17_recgenpt_20_20_ptbins32"
#source run_printPlots_ppMC_JERS.sh "4" "09.13.17_outputCondor/recpt20_genpt20_ptbins32" "09-13-17" "1.5eta2.0" "15eta20_9.13.17_recgenpt_20_20_ptbins32"

#source run_printPlots_ppMC_JERS.sh "4" "09.18.17_outputCondor/ppMC_JERS_gendrjt0p2" "09-18-17" "0.0eta0.5" "00eta05_9.18.17_gendrjt0p2"
#source run_printPlots_ppMC_JERS.sh "4" "09.18.17_outputCondor/ppMC_JERS_gendrjt0p2" "09-18-17" "0.5eta1.0" "05eta10_9.18.17_gendrjt0p2"
#source run_printPlots_ppMC_JERS.sh "4" "09.18.17_outputCondor/ppMC_JERS_gendrjt0p2" "09-18-17" "1.0eta1.5" "10eta15_9.18.17_gendrjt0p2"
#source run_printPlots_ppMC_JERS.sh "4" "09.18.17_outputCondor/ppMC_JERS_gendrjt0p2" "09-18-17" "1.5eta2.0" "15eta20_9.18.17_gendrjt0p2"

#source run_printPlots_ppMC_JERS.sh "4" "09.18.17_outputCondor/ppMC_JERS_gendrjt0p1" "09-18-17" "0.0eta0.5" "00eta05_9.18.17_gendrjt0p1"
#source run_printPlots_ppMC_JERS.sh "4" "09.18.17_outputCondor/ppMC_JERS_gendrjt0p1" "09-18-17" "0.5eta1.0" "05eta10_9.18.17_gendrjt0p1"
#source run_printPlots_ppMC_JERS.sh "4" "09.18.17_outputCondor/ppMC_JERS_gendrjt0p1" "09-18-17" "1.0eta1.5" "10eta15_9.18.17_gendrjt0p1"
#source run_printPlots_ppMC_JERS.sh "4" "09.18.17_outputCondor/ppMC_JERS_gendrjt0p1" "09-18-17" "1.5eta2.0" "15eta20_9.18.17_gendrjt0p1"


#source run_printPlots_ppMC_JERS.sh "4" "09.19.17_outputCondor" "09-19-17" "0.0eta2.0" "00eta20_9.19.17"

#source run_printPlots_ppMC_JERS.sh "4" "09.19.17_outputCondor" "09-19-17" "0.0eta0.5" "00eta05_9.19.17"
#source run_printPlots_ppMC_JERS.sh "4" "09.19.17_outputCondor" "09-19-17" "0.5eta1.0" "05eta10_9.19.17"
#source run_printPlots_ppMC_JERS.sh "4" "09.19.17_outputCondor" "09-19-17" "1.0eta1.5" "10eta15_9.19.17"
#source run_printPlots_ppMC_JERS.sh "4" "09.19.17_outputCondor" "09-19-17" "1.5eta2.0" "15eta20_9.19.17"

#source run_printPlots_ppMC_JERS.sh "4" "10.02.17_outputCondor" "10-02-17" "0.0eta2.0" "00eta20_10.2.17"

#source run_printPlots_ppMC_JERS.sh "4" "10.02.17_outputCondor" "10-02-17" "0.0eta0.5" "00eta05_10.2.17_badTF1"
#source run_printPlots_ppMC_JERS.sh "4" "10.02.17_outputCondor" "10-02-17" "0.5eta1.0" "05eta10_10.2.17_badTF1"
#source run_printPlots_ppMC_JERS.sh "4" "10.02.17_outputCondor" "10-02-17" "1.0eta1.5" "10eta15_10.2.17_badTF1"
#source run_printPlots_ppMC_JERS.sh "4" "10.02.17_outputCondor" "10-02-17" "1.5eta2.0" "15eta20_10.2.17_badTF1"

#source run_printPlots_ppMC_JERS.sh "4" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.3.17"
