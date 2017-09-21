#!/bin/bash






####################################################
##################### JETPLOTS #####################
####################################################

echo ""
echo "compiling code(s)"
echo ""
#
echo "...ppMC jetPlots..."
rootcompile readForests_ppMC_jetPlots.C
#
echo "...ppData jetPlots..."
rootcompile readForests_ppData_jetPlots.C
#
echo ""
echo "done compiling code(s), submitting."
echo ""

### INCL. JET SPECTRA, ENTIRE |y| RANGE
#source run_readForests_jetPlots.sh 4 "0.0" "4.7"
#sleep 60s

### INCL. JET SPECTRA, LARGER |y| REGIONS
source run_readForests_jetPlots.sh 4 "0.0" "2.0"
#source run_readForests_jetPlots.sh 4 "2.0" "3.0"
#source run_readForests_jetPlots.sh 4 "3.2" "4.7"
#sleep 60s

#### INCL. JET SPECTRA |y| bins
#source run_readForests_jetPlots.sh 4 "0.0" "0.5"
#sleep 60s
#source run_readForests_jetPlots.sh 4 "0.5" "1.0"
#sleep 60s
#source run_readForests_jetPlots.sh 4 "1.0" "1.5"
#sleep 60s
#source run_readForests_jetPlots.sh 4 "1.5" "2.0"
#sleep 60s
#source run_readForests_jetPlots.sh 4 "2.0" "2.5"
#sleep 60s
#source run_readForests_jetPlots.sh 4 "2.5" "3.0"
#sleep 60s
#source run_readForests_jetPlots.sh 4 "3.2" "4.7"
#sleep 60s


echo ""
echo "done running code(s)."
echo ""

source askCondor.sh 100 30

return








#################################################
###################### JERS #####################
#################################################
echo ""
echo "compiling code(s)"
echo ""

echo "...ppMC JERS..."
rootcompile readForests_ppMC_JERS.C
#
##echo "...ppData JERS..."
##rootcompile readForests_ppData_JERS.C
#
echo ""
echo "done compiling code(s), submitting."
echo ""

source run_readForests_JERS.sh 4 "0.0" "0.5"
#sleep 60s
source run_readForests_JERS.sh 4 "0.5" "1.0"
#sleep 60s
source run_readForests_JERS.sh 4 "1.0" "1.5"
#sleep 60s
source run_readForests_JERS.sh 4 "1.5" "2.0"
#sleep 60s
#source run_readForests_JERS.sh 4 "2.0" "2.5"
#sleep 60s
#source run_readForests_JERS.sh 4 "2.5" "3.0"
#sleep 60s
#source run_readForests_JERS.sh 4 "3.2" "4.7"
#sleep 60s


echo ""
echo "done running code(s)."
echo ""

source askCondor.sh 100 60

return



