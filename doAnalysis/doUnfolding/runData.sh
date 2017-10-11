#!/bin/bash



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


