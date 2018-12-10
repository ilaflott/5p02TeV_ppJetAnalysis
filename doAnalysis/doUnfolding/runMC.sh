#!/bin/bash

### ----------------------------------------------------------------------------------
### ----------------------------------------------------------------------------------
if [ $# -eq 9 ]
then
    
    dirTag=$1
    MM=$2
    DD=$3
    YY=$4
    simpbins=$5
    runBayes=$6
    runSVD=$7
    SVDkReg=$8
    BayeskIter=$9
    
else
    echo "usage:"
    echo "source runMC.sh [dirTag] [MM] [DD] [YY] [useSimpBins] [runBayes] [runSVD] [SVDkReg] [BayeskIter]"
    return
fi

condorDate="${MM}.${DD}.${YY}"
sampleDate="${MM}-${DD}-${YY}"


### ----------------------------------------------------------------------------------
if [ $runBayes -eq 1 ]
then
    echo ""
    echo "kIter=${BayeskIter}"
    echo ""
    
    #USE: "source run_bayesUnfoldMCSpectra.sh [R=3,4] [condorDir] [date_output] [etaBin] [etaBinOut] [useSimpleBinning]"
    source run_bayesUnfoldMCSpectra.sh  "4" "${condorDate}_outputCondor" "${sampleDate}" "${dirTag}" "Bayes_${condorDate}_${dirTag}" "${simpbins}" "${BayeskIter}"
    
    echo ""
    echo "Bayes MC Unfolding Done."
    echo ""
    
fi

### ----------------------------------------------------------------------------------
if [ $runSVD -eq 1 ]
then
    echo ""
    echo "kReg=${SVDkReg}"
    echo ""

    #USE: "source run_SVDUnfoldMCSpectra.sh [R=3,4] [condorDir] [date_output] [etaBin] [etaBinOut] [kReg] [useSimpleBinning]"    
    source run_SVDUnfoldMCSpectra.sh "4" "${condorDate}_outputCondor" "${sampleDate}" "${dirTag}" "SVD_${condorDate}_${dirTag}" "${SVDkReg}" "${simpbins}"
    
    echo ""
    echo "SVD MC Unfolding Done."
    echo ""
    
fi

echo ""
echo "done unfolding MC."
echo ""

return


#source run_unfoldMCSpectra.sh "4" "10.02.17_outputCondor" "10-02-17" "0.0eta2.0" "00eta20_10.2.17" 10 0
#source run_unfoldMCSpectra.sh "4" "10.02.17_outputCondor" "10-02-17" "0.0eta2.0" "00eta20_10.2.17" 6 1

#source run_unfoldMCSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.3.17" 6 0
#source run_unfoldMCSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "00eta20_10.3.17" 6 1

#source run_bayesUnfoldMCSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "Bayes_00eta20_10.9.17" 6 0
#source run_bayesUnfoldMCSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "Bayes_00eta20_10.9.17_alt2" 6 0
#source run_bayesUnfoldMCSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "Bayes_00eta20_10.9.17" 6 1

#source run_SVDUnfoldMCSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "SVD_00eta20_10.9.17" 6 0
#source run_SVDUnfoldMCSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "SVD_00eta20_10.9.17_alt2" 6 0
#source run_SVDUnfoldMCSpectra.sh "4" "10.03.17_outputCondor" "10-03-17" "0.0eta2.0" "SVD_00eta20_10.9.17" 6 1
