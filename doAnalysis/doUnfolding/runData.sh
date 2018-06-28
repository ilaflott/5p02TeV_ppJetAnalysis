#!/bin/bash

### ----------------------------------------------------------------------------------
if [ $# -eq 8 ]
then
    datadirTag=$1
    MCdirTag=$2
    MM=$3
    DD=$4
    YY=$5
    runBayes=$6
    runSVD=$7
    simpbins=$8
else
    echo "usage:"
    echo "source runData.sh [datadirTag] [MCdirTag] [MM] [DD] [YY] [runBayes] [runSVD] [useSimpBins]"
    return
fi

#echo "in runData.sh, datadirTag=${datadirTag}"

condorDate="${MM}.${DD}.${YY}"
sampleDate="${MM}-${DD}-${YY}"


### ----------------------------------------------------------------------------------
if [ $runBayes -eq 1 ]
then
    echo ""
    echo "compiling bayesUnfoldDataSpectra.C"
    echo ""
    rooUnfoldCompile bayesUnfoldDataSpectra.C
    
# "source run_bayesUnfoldDataSpectra.sh [R=3,4] [dataCondorDir] [data_date] [MCCondorDir] [MC_date] [MCetaBin] [etaBinOut] [useSimpleBinning]"    
    
    source run_bayesUnfoldDataSpectra.sh  "4" "${condorDate}_outputCondor" "${sampleDate}" "${datadirTag}" "${condorDate}_outputCondor" "${sampleDate}" "${MCdirTag}" "Bayes_${condorDate}_${MCdirTag}" "${simpbins}"
    
    echo ""
    echo "Bayes Data Unfolding Done."
    echo ""
    
fi

### ----------------------------------------------------------------------------------
if [ $runSVD -eq 1 ]
then
    if [ ${simpbins} -eq 1 ]
    then
	kReg=9
    else
	kReg=8
    fi

    echo ""
    echo "compiling SVDUnfoldDataSpectra.C"
    echo ""
    echo "kReg=${kReg}"
    rooUnfoldCompile SVDUnfoldDataSpectra.C

# "source run_SVDUnfoldDataSpectra.sh [R=3,4] [dataCondorDir] [data_date] [dataEtaBin] [MCCondorDir] [MC_date] [MCetaBin] [etaBinOut] [kReg] [useSimpleBinning]"    

    source run_SVDUnfoldDataSpectra.sh  "4" "${condorDate}_outputCondor" "${sampleDate}" "${datadirTag}" "${condorDate}_outputCondor" "${sampleDate}" "${MCdirTag}" "SVD_${condorDate}_${MCdirTag}" "${kReg}" "${simpbins}" 
    

    echo ""
    echo "SVD Data Unfolding Done."
    echo ""
    
fi


echo ""
echo "done unfolding MC."
echo ""

return










































#
#
#
#
#
#
#### ----------------------------------------------------------------------------------
#### ----------------------------------------------------------------------------------
#### OLD
#### ----------------------------------------------------------------------------------
#### ----------------------------------------------------------------------------------
#
#echo ""
#echo "compiling SVDUnfoldDataSpectra"
#echo ""
#rooUnfoldCompile SVDUnfoldDataSpectra.C
#
#echo ""
#echo "done compiling unfoldDataSpectra"
#echo ""
##sleep 1s
#
#source run_SVDUnfoldDataSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.9.17_alt"  0 5 
#source run_SVDUnfoldDataSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.9.17_alt"  1 15
#
#echo ""
#echo "done unfolding Data w/ SVD."
#echo ""
#
##return
#
#
#
#
#
#
#
#
#
#
#
#
#### ----------------------------------------------------------------------------------
#
#echo ""
#echo "compiling bayesUnfoldDataSpectra"
#echo ""
#rooUnfoldCompile bayesUnfoldDataSpectra.C
#
#echo ""
#echo "done compiling unfoldDataSpectra"
#echo ""
##sleep 1s
#
#source run_bayesUnfoldDataSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.9.17_alt"  0  
#source run_bayesUnfoldDataSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.9.17_alt"  1 
#
#echo ""
#echo "done unfolding Data w/ Bayes."
#echo ""
#
#return
#
#
#
#
#
#





    #source run_bayesUnfoldDataSpectra.sh  "4" "${condorDate}_outputCondor" "${sampleDate}" "0.0eta2.0${datadirTag}" "${condorDate}_outputCondor" "${sampleDate}" "0.0eta2.0${MCdirTag}" "Bayes_00eta20_${condorDate}${MCdirTag}" 0
    #source run_bayesUnfoldDataSpectra.sh  "4" "${condorDate}_outputCondor" "${sampleDate}" "0.0eta2.0${datadirTag}" "${condorDate}_outputCondor" "${sampleDate}" "0.0eta2.0${MCdirTag}" "Bayes_00eta20_${condorDate}${MCdirTag}" 1
    

#source run_SVDUnfoldDataSpectra.sh "4" "${condorDate}_outputCondor" "${sampleDate}" "{datadirTag}" "${condorDate}_outputCondor" "${sampleDate}" "${MCdirTag}" "SVD_${condorDate}_${MCdirTag}" "${kReg}" "${simpbins}"



#source run_SVDUnfoldDataSpectra.sh "4" "${condorDate}_outputCondor" "${sampleDate}" "0.0eta2.0${datadirTag}" "${condorDate}_outputCondor" "${sampleDate}" "0.0eta2.0${MCdirTag}" "SVD_00eta20_${condorDate}${MCdirTag}" 6 "${simpbins}"
    

