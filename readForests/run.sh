#!/bin/bash



echo ""
echo "compiling code(s)"
echo ""

#echo "...ppMC jetPlots..."
#rootcompile readForests_ppMC_jetPlots.C

echo "...ppMC JERS..."
rootcompile readForests_ppMC_JERS.C

#echo "...ppData jetPlots..."
#rootcompile readForests_ppData_jetPlots.C

#echo "...ppData jetMult..."
#rootcompile readForests_ppData_jetMult.C

#echo "...ppMC jetMult..."
#rootcompile readForests_ppMC_jetMult.C

echo ""
echo "done compiling code(s), submitting."
echo ""


####################################################

#source run_readForests_JERS.sh 4 "0.0" "0.5" 
#source run_readForests_JERS.sh 4 "0.5" "1.0" 
#source run_readForests_JERS.sh 4 "1.0" "1.5" 
#source run_readForests_JERS.sh 4 "1.5" "2.0" 

#source run_readForests_JERS.sh 4 "0.0" "1.0"
#source run_readForests_JERS.sh 4 "1.0" "2.0"

source run_readForests_JERS.sh 4 "0.0" "2.0"

####################################################

#source run_readForests_jetPlots.sh 4 "0.0" "0.5" "ppData"
#source run_readForests_jetPlots.sh 4 "0.5" "1.0" "ppData"
#source run_readForests_jetPlots.sh 4 "1.0" "1.5" "ppData"
#source run_readForests_jetPlots.sh 4 "1.5" "2.0" "ppData"

#source run_readForests_jetPlots.sh 4 "0.0" "1.0" "ppData"
#source run_readForests_jetPlots.sh 4 "1.0" "2.0" "ppData"

#source run_readForests_jetPlots.sh 4 "0.0" "2.0" "ppData"

####################################################

#source run_readForests_jetPlots.sh 4 "0.0" "0.5" "ppMC"
#source run_readForests_jetPlots.sh 4 "0.5" "1.0" "ppMC"
#source run_readForests_jetPlots.sh 4 "1.0" "1.5" "ppMC"
#source run_readForests_jetPlots.sh 4 "1.5" "2.0" "ppMC"

#source run_readForests_jetPlots.sh 4 "0.0" "1.0" "ppMC"
#source run_readForests_jetPlots.sh 4 "1.0" "2.0" "ppMC"

#source run_readForests_jetPlots.sh 4 "0.0" "2.0" "ppMC"

#askCondor  1 1
#sleep 120s


echo ""
echo "done running code(s)."
echo ""

#askCondor  100 30
askCondor  100 100

return













## INCL. JET SPECTRA, ENTIRE |y| RANGE
#source run_readForests_jetPlots.sh 4 "0.0" "4.7"
#source run_readForests_jetPlots.sh 4 "0.0" "2.0"

#source run_readForests_jetPlots.sh 4 "0.0" "0.5" "ppData"
#source run_readForests_jetPlots.sh 4 "0.5" "1.0" "ppData"
#source run_readForests_jetPlots.sh 4 "0.0" "0.5" "ppMC"
#source run_readForests_jetPlots.sh 4 "0.5" "1.0" "ppMC"

#source run_readForests_jetPlots.sh 3 "0.0" "0.5" "ppData"
#source run_readForests_jetPlots.sh 3 "0.5" "1.0" "ppData"

#source run_readForests_jetPlots.sh 4 "0.0" "0.5" "ppData"
#source run_readForests_jetPlots.sh 4 "0.5" "1.0" "ppData"
#source run_readForests_jetPlots.sh 4 "0.0" "1.0" "ppData"

#source run_readForests_jetPlots.sh 3 "0.0" "0.5" "ppMC"
#source run_readForests_jetPlots.sh 3 "0.5" "1.0" "ppMC"

#source run_readForests_jetPlots.sh 4 "0.0" "0.5" "ppMC"
#source run_readForests_jetPlots.sh 4 "0.5" "1.0" "ppMC"
#source run_readForests_jetPlots.sh 4 "0.0" "1.0" "ppMC"

#### INCL. JET SPECTRA |y| bins
#source run_readForests_jetPlots.sh 4 "0.0" "0.5"
#askCondor  1 1
#sleep 120s
#source run_readForests_jetPlots.sh 4 "0.5" "1.0"
#askCondor  1 1
#sleep 120s
#source run_readForests_jetPlots.sh 4 "1.0" "1.5"
#askCondor  1 1
#sleep 120s
#source run_readForests_jetPlots.sh 4 "1.5" "2.0"
#askCondor  1 1
#sleep 120s

#source run_readForests_jetPlots.sh 4 "2.0" "2.5"
#askCondor  1 1
#sleep 60s
#source run_readForests_jetPlots.sh 4 "2.5" "3.0"
#askCondor  1 1
#sleep 60s
#source run_readForests_jetPlots.sh 4 "3.2" "4.7"
#askCondor  1 1
#sleep 60s


### INCL. JET SPECTRA, LARGER |y| REGIONS
#source run_readForests_jetPlots.sh 4 "0.0" "2.0"

#source run_readForests_jetMult.sh 4 "0.0" "2.0"












##temp
source run_readForests_JERS.sh 4 "0.0" "2.0"
askCondor  100 60
return
##temp




#################################################
###################### JERS #####################
#################################################
echo ""
echo "compiling code(s)"
echo ""

echo "...ppMC JERS..."
rootcompile readForests_ppMC_JERS.C

##echo "...ppData JERS..."
##rootcompile readForests_ppData_JERS.C

echo ""
echo "done compiling code(s), submitting."
echo ""

source run_readForests_JERS.sh 4 "0.0" "2.0"
askCondor  1 1
sleep 60s

source run_readForests_JERS.sh 4 "0.0" "0.5"
askCondor  1 1
sleep 60s
source run_readForests_JERS.sh 4 "0.5" "1.0"
askCondor  1 1
sleep 60s
source run_readForests_JERS.sh 4 "1.0" "1.5"
askCondor  1 1
sleep 60s
source run_readForests_JERS.sh 4 "1.5" "2.0"
askCondor  1 1
sleep 60s

#source run_readForests_JERS.sh 4 "2.0" "2.5"
#sleep 60s
#source run_readForests_JERS.sh 4 "2.5" "3.0"
#sleep 60s
#source run_readForests_JERS.sh 4 "3.2" "4.7"
#sleep 60s


echo ""
echo "done running code(s)."
echo ""

askCondor  1 1


return

