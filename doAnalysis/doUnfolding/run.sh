#!/bin/bash



### ----------------------------------------------------------------------
echo ""
echo "compiling SVDUnfoldDataSpectra"
echo ""
rooUnfoldCompile SVDUnfoldDataSpectra.C

echo ""
echo "compiling bayesUnfoldDataSpectra"
echo ""
rooUnfoldCompile bayesUnfoldDataSpectra.C


echo ""
echo ""
echo "------------------------------------------------------------------"
echo ""
echo ""


echo ""
echo "running Data."
echo ""

source runData.sh "_rec56_HLTCalo" "_gendr0p1_HLTCalo" "10" "10" "17" 0 1
#source runData.sh "_rec56_HLTCalo" "_gendr0p2_HLTCalo" "10" "10" "17" 1 1
#source runData.sh "_rec56_HLTPF" "_gendr0p1_HLTPF" "10" "10" "17" 1 1
#source runData.sh "_rec56_HLTPF" "_gendr0p2_HLTPF" "10" "10" "17" 1 1

#source runData.sh "_rec64_HLTCalo" "_gendr0p1_HLTCalo" "10" "10" "17" 1 1
#source runData.sh "_rec64_HLTCalo" "_gendr0p2_HLTCalo" "10" "10" "17" 1 1
#source runData.sh "_rec64_HLTPF" "_gendr0p1_HLTPF" "10" "10" "17" 1 1
#source runData.sh "_rec64_HLTPF" "_gendr0p2_HLTPF" "10" "10" "17" 1 1

#source runData.sh "_rec56_HLTCalo" "_rec56_gen32_HLTCalo" "10" "11" "17" 1 1
#source runData.sh "_rec56_HLTPF" "_rec56_gen32_HLTPF" "10" "11" "17" 1 1
#source runData.sh "_rec56_HLTCalo" "_rec74_gen43_HLTCalo" "10" "11" "17" 1 1
#source runData.sh "_rec56_HLTPF" "_rec74_gen43_HLTPF" "10" "11" "17" 1 1

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

#source runMC.sh "_rec56_HLTCalo" "10" "10" "17" 0 1
#source runMC.sh "_rec64_HLTCalo" "10" "10" "17" 0 1

##prefer these two right here 01/23/18
source runMC.sh "_gendr0p1_HLTCalo" "10" "10" "17"  0 1
#source runMC.sh "_gendr0p2_HLTCalo" "10" "10" "17" 1 0

#source runMC.sh   "_rec56_gen32_HLTCalo" "10" "11" "17" 1 1
#source runMC.sh  "_rec56_gen32_HLTPF" "10" "11" "17" 1 1
#source runMC.sh  "_rec74_gen43_HLTCalo" "10" "11" "17" 1 1
#source runMC.sh  "_rec74_gen43_HLTPF" "10" "11" "17" 1 1

echo ""
echo "done w/ MC."
echo ""


return

