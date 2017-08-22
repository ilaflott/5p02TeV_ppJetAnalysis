#!/bin/bash

##############################################
##############################################





#echo ""
#echo "compiling printPlots_jetPlots.C"
#echo ""
#
#rootcompile printPlots_jetPlots.C
#
#echo ""
#echo "done compiling. Running!"
#echo ""
#
#source run_printPlots_jetPlots.sh 4 "08-15-17" "3.2eta4.7" "32eta47_8.15.17"
#source run_printPlots_jetPlots.sh 4 "08-15-17" "2.0eta3.0" "20eta30_8.15.17"
#source run_printPlots_jetPlots.sh 4 "08-15-17" "0.0eta2.0" "00eta20_8.15.17"

echo ""
echo "done running jetPlots!"
echo ""

##############################################
##############################################






#echo ""
#echo "compiling printPlots_jetIDPlots.C"
#echo ""
#
#rootcompile printPlots_jetIDPlots.C
#
#echo ""
#echo "done compiling. Running!"
#echo ""
#
#source run_printPlots_jetIDPlots.sh 4 "08-15-17" "3.2eta4.7" "32eta47_8.15.17"
#source run_printPlots_jetIDPlots.sh 4 "08-15-17" "2.0eta3.0" "20eta30_8.15.17"
#source run_printPlots_jetIDPlots.sh 4 "08-15-17" "0.0eta2.0" "00eta20_8.15.17"

echo ""
echo "done running jetIDPlots!"
echo ""

##############################################
##############################################


#echo ""
#echo "compiling printPlots_jetTrig.C"
#echo ""
#
#rootcompile printPlots_jetTrig.C
#
#echo ""
#echo "done compiling. Running!"
#echo ""
#
#source run_printPlots_jetTrig.sh 4 "08-15-17" "3.2eta4.7" "32eta47_8.15.17"
#source run_printPlots_jetTrig.sh 4 "08-15-17" "2.0eta3.0" "20eta30_8.15.17"
#source run_printPlots_jetTrig.sh 4 "08-15-17" "0.0eta2.0" "00eta20_8.15.17"

echo ""
echo "done running jetTrig!"
echo ""

##############################################
##############################################



echo ""
echo "compiling printPlots_ppMC_JERS.C"
echo ""

rootcompile printPlots_ppMC_JERS.C

echo ""
echo "done compiling. Running!"
echo ""

source run_printPlots_ppMC_JERS.sh 4 "08-16-17" "0.0eta0.5" "00eta05_8.16.17"
source run_printPlots_ppMC_JERS.sh 4 "08-16-17" "0.5eta1.0" "05eta10_8.16.17"
source run_printPlots_ppMC_JERS.sh 4 "08-16-17" "1.0eta1.5" "10eta15_8.16.17"
source run_printPlots_ppMC_JERS.sh 4 "08-16-17" "1.5eta2.0" "15eta20_8.16.17"

source run_printPlots_ppMC_JERS.sh 4 "08-16-17" "2.0eta2.5" "20eta25_8.16.17"
source run_printPlots_ppMC_JERS.sh 4 "08-16-17" "2.5eta3.0" "25eta30_8.16.17"

source run_printPlots_ppMC_JERS.sh 4 "08-16-17" "3.2eta4.7" "32eta47_8.16.17"

echo ""
echo "done running ppMC JERS!"
echo ""

##############################################
##############################################

return
