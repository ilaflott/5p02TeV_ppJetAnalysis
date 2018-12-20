#!/bin/bash

### ----------------------------------------------------------------------------------
if [ $# -eq 10 ]
then
    datadirTag=$1
    MM=$2
    DD=$3
    YY=$4
    simpbins=$5
    runBayes=$6
    runSVD=$7
    SVDkReg=$8
    BayeskIter=$9
    useNPCorr=${10}
else
    echo "usage:"
    echo "source runData_wNLO.sh [datadirTag] [MM] [DD] [YY] [useSimpBins] [runBayes] [runSVD] [SVDkReg] [BayeskIter] [useNPCorr]"
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
    
    source run_bayesUnfoldDataSpectra_wNLO.sh  "4" "${condorDate}_outputCondor" "${sampleDate}" "${datadirTag}" "Bayes_${condorDate}" "${simpbins}" "${BayeskIter}" "${useNPCorr}"
    
    echo ""
    echo "Bayes Data Unfolding Done."
    echo ""
    
fi

### ----------------------------------------------------------------------------------
if [ $runSVD -eq 1 ]
then
    
    echo "kReg=${SVDkReg}"    
    source run_SVDUnfoldDataSpectra_wNLO.sh  "4" "${condorDate}_outputCondor" "${sampleDate}" "${datadirTag}" "SVD_${condorDate}" "${SVDkReg}" "${simpbins}" "${useNPCorr}"
    
    
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
