#!/bin/bash



### ----------------------------------------------------------------------
echo ""
echo "compiling SVDUnfoldDataSpectra"
echo ""
rooUnfoldCompile SVDUnfoldDataSpectra.C

echo ""
echo "done compiling unfoldDataSpectra"
echo ""

echo ""
echo "compiling bayesUnfoldDataSpectra"
echo ""
rooUnfoldCompile bayesUnfoldDataSpectra.C

echo ""
echo "done compiling unfoldDataSpectra"
echo ""

echo ""
echo " ------------------------------------------------------------------"
echo " ------------------------------------------------------------------"
echo " ------------------------------------------------------------------"
echo ""

echo ""
echo "running Data."
echo ""

source runData.sh "_rec56_HLTCalo" "10" "10" "17" 1 1
source runData.sh "_rec64_HLTCalo" "10" "10" "17" 1 1

source runData.sh "_rec56_HLTPF" "10" "10" "17" 1 1
source runData.sh "_rec64_HLTPF" "10" "10" "17" 1 1


echo ""
echo "done w/ Data."
echo ""


#return 







### ----------------------------------------------------------------------
echo ""
echo "compiling SVDUnfoldMCSpectra"
echo ""
rooUnfoldCompile SVDUnfoldMCSpectra.C

echo ""
echo "done compiling unfoldMCSpectra"
echo ""

echo ""
echo "compiling bayesUnfoldMCSpectra"
echo ""
rooUnfoldCompile bayesUnfoldMCSpectra.C

echo ""
echo "done compiling unfoldMCSpectra"
echo ""


echo ""
echo " ------------------------------------------------------------------"
echo " ------------------------------------------------------------------"
echo " ------------------------------------------------------------------"
echo ""

echo ""
echo "running MC."
echo ""

#note; _HLTCalo/HLTPF added to JERS dirs for 10.10.17 sample for convenience of I/O.
# there is no trigger emulation in the MC samples the label is superficial and does not matter.
source runMC.sh "_rec56_HLTCalo" "10" "10" "17" 1 1
source runMC.sh "_rec64_HLTCalo" "10" "10" "17" 1 1

echo ""
echo "done w/ MC."
echo ""


return

