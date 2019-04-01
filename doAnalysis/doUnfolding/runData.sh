#!/bin/bash

### ----------------------------------------------------------------------------------
if [ $# -eq 11 ]
then
    datadirTag=$1
    MCdirTag=$2
    MM=$3
    DD=$4
    YY=$5
    simpbins=$6
    runBayes=$7
    runSVD=$8
    SVDkReg=$9
    BayeskIter=${10}
    etabinint=${11}
else
    echo "usage:"
    echo "source runData.sh [datadirTag] [MCdirTag] [MM] [DD] [YY] [useSimpBins] [runBayes] [runSVD] [SVDkReg] [BayeskIter] [etabinint]"
    return
fi

#echo "in runData.sh, datadirTag=${datadirTag}"

condorDate="${MM}.${DD}.${YY}"
sampleDate="${MM}-${DD}-${YY}"


### ----------------------------------------------------------------------------------
if [ $runBayes -eq 1 ]
then
    echo ""
    echo "kIter=${BayeskIter}"
    echo ""
    
    ## echo "source run_bayesUnfoldDataSpectra.sh [R=3,4] [dataCondorDir] [data_date] [MCCondorDir] [MC_date] [MCetaBin] [tagOut] [useSimpleBinning] [BayeskIter] [etabinint]"
    source run_bayesUnfoldDataSpectra.sh  "4" "${condorDate}_outputCondor" "${sampleDate}" "${datadirTag}" "${condorDate}_outputCondor" "${sampleDate}" "${MCdirTag}" "Bayes_${condorDate}_${MCdirTag}" "${simpbins}" "${BayeskIter}" "${etabinint}"
    
    echo ""
    echo "Bayes Data Unfolding Done."
    echo ""
    
fi

### ----------------------------------------------------------------------------------
if [ $runSVD -eq 1 ]
then
    echo ""
    echo "kReg=${SVDkReg}"
    echo ""
    
# "source run_SVDUnfoldDataSpectra.sh [R=3,4] [dataCondorDir] [data_date] [dataEtaBin] [MCCondorDir] [MC_date] [MCetaBin] [etaBinOut] [kReg] [useSimpleBinning]"    
    echo "kReg setting = ${SVDkReg}"
    source run_SVDUnfoldDataSpectra.sh  "4" "${condorDate}_outputCondor" "${sampleDate}" "${datadirTag}" "${condorDate}_outputCondor" "${sampleDate}" "${MCdirTag}" "SVD_${condorDate}_${MCdirTag}" "${SVDkReg}" "${simpbins}" 

    echo ""
    echo "SVD Data Unfolding Done."
    echo ""
    
fi


echo ""
echo "done unfolding Data."
echo ""

return






#source run_SVDUnfoldDataSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.9.17_alt"  0 5 
#source run_SVDUnfoldDataSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.9.17_alt"  1 15
#source run_bayesUnfoldDataSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.9.17_alt"  0  
#source run_bayesUnfoldDataSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.9.17_alt"  1 
