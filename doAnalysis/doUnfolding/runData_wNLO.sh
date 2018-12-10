#!/bin/bash

### ----------------------------------------------------------------------------------
if [ $# -eq 10 ]
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
else
    echo "usage:"
    echo "source runData_wNLO.sh [datadirTag] [MCdirTag] [MM] [DD] [YY] [useSimpBins] [runBayes] [runSVD] [SVDkReg] [BayeskIter]"
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

    #"source run_bayesUnfoldDataSpectra_wNLO.sh [R=3,4] [dataCondorDir] [data_date] [MCCondorDir] [MC_date] [MCetaBin] [etaBinOut] [useSimpleBinning] [kIter]"        
    source run_bayesUnfoldDataSpectra_wNLO.sh  "4" "${condorDate}_outputCondor" "${sampleDate}" "${datadirTag}" "${condorDate}_outputCondor" "${sampleDate}" "${MCdirTag}" "Bayes_${condorDate}" "${simpbins}" "${BayeskIter}"
    
    echo ""
    echo "Bayes Data Unfolding Done."
    echo ""
    
fi

### ----------------------------------------------------------------------------------
if [ $runSVD -eq 1 ]
then
    
    echo "kReg=${SVDkReg}"
    
# "source run_SVDUnfoldDataSpectra_wNLO.sh [R=3,4] [dataCondorDir] [data_date] [dataEtaBin] [MCCondorDir] [MC_date] [MCetaBin] [etaBinOut] [kReg] [useSimpleBinning]"        
    source run_SVDUnfoldDataSpectra_wNLO.sh  "4" "${condorDate}_outputCondor" "${sampleDate}" "${datadirTag}" "${condorDate}_outputCondor" "${sampleDate}" "${MCdirTag}" "SVD_${condorDate}_${MCdirTag}" "${SVDkReg}" "${simpbins}" 
    
    
    echo ""
    echo "SVD Data Unfolding Done."
    echo ""
    
fi


echo ""
echo "done unfolding Data."
echo ""

return



#source run_SVDUnfoldDataSpectra_wNLO.sh "4" "10.03.17_outputCondor" "10-03-17" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.9.17_alt"  0 5 
#source run_SVDUnfoldDataSpectra_wNLO.sh "4" "10.03.17_outputCondor" "10-03-17" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.9.17_alt"  1 15
#source run_bayesUnfoldDataSpectra_wNLO.sh "4" "10.03.17_outputCondor" "10-03-17" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.9.17_alt"  0  
#source run_bayesUnfoldDataSpectra_wNLO.sh "4" "10.03.17_outputCondor" "10-03-17" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.9.17_alt"  1 
