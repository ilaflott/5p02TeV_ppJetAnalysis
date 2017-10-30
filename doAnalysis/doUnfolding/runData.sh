#!/bin/bash

### ----------------------------------------------------------------------------------
datadirTag=$1
MCdirTag=$2
MM=$3
DD=$4
YY=$5
runBayes=$6
runSVD=$7

condorDate="${MM}.${DD}.${YY}"
sampleDate="${MM}-${DD}-${YY}"


### ----------------------------------------------------------------------------------
if [ $runBayes -eq 1 ]
then
    
    source run_bayesUnfoldDataSpectra.sh  "4" "${condorDate}_outputCondor" "${sampleDate}" "0.0eta2.0${datadirTag}" "${condorDate}_outputCondor" "${sampleDate}" "0.0eta2.0${MCdirTag}" "Bayes_00eta20_${condorDate}${MCdirTag}" 0
    #source run_bayesUnfoldDataSpectra.sh  "4" "${condorDate}_outputCondor" "${sampleDate}" "0.0eta2.0${datadirTag}" "${condorDate}_outputCondor" "${sampleDate}" "0.0eta2.0${MCdirTag}" "Bayes_00eta20_${condorDate}${MCdirTag}" 1
    
    echo ""
    echo "Bayes Data Unfolding Done."
    echo ""
    
fi

### ----------------------------------------------------------------------------------
if [ $runSVD -eq 1 ]
then
    
    source run_SVDUnfoldDataSpectra.sh "4" "${condorDate}_outputCondor" "${sampleDate}" "0.0eta2.0${datadirTag}" "${condorDate}_outputCondor" "${sampleDate}" "0.0eta2.0${MCdirTag}" "SVD_00eta20_${condorDate}${MCdirTag}" 6 0
    source run_SVDUnfoldDataSpectra.sh "4" "${condorDate}_outputCondor" "${sampleDate}" "0.0eta2.0${datadirTag}" "${condorDate}_outputCondor" "${sampleDate}" "0.0eta2.0${MCdirTag}" "SVD_00eta20_${condorDate}${MCdirTag}" 6 1
    
    echo ""
    echo "SVD Data Unfolding Done."
    echo ""
    
fi


echo ""
echo "done unfolding MC."
echo ""

return
















































### ----------------------------------------------------------------------------------
### ----------------------------------------------------------------------------------
### OLD
### ----------------------------------------------------------------------------------
### ----------------------------------------------------------------------------------

echo ""
echo "compiling SVDUnfoldDataSpectra"
echo ""
rooUnfoldCompile SVDUnfoldDataSpectra.C

echo ""
echo "done compiling unfoldDataSpectra"
echo ""
#sleep 1s

source run_SVDUnfoldDataSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.9.17_alt"  0 5 
source run_SVDUnfoldDataSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.9.17_alt"  1 15

echo ""
echo "done unfolding Data w/ SVD."
echo ""

#return












### ----------------------------------------------------------------------------------

echo ""
echo "compiling bayesUnfoldDataSpectra"
echo ""
rooUnfoldCompile bayesUnfoldDataSpectra.C

echo ""
echo "done compiling unfoldDataSpectra"
echo ""
#sleep 1s

source run_bayesUnfoldDataSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.9.17_alt"  0  
source run_bayesUnfoldDataSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.9.17_alt"  1 

echo ""
echo "done unfolding Data w/ Bayes."
echo ""

return


