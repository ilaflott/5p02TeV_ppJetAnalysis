#!/bin/bash


#### USAGE
## "source runData.sh [datadirTag] [MCdirTag] [MM] [DD] [YY] [runBayes] [runSVD] [useSimpBins]"
#
#echo ""
#echo "running Data."
#echo ""
#
#mv output/unfoldDataSpectra/* output/unfoldDataSpectra/temp/.
#
####UPDATE runData.sh ALA runMC.sh WHEN YOU GET A CHANCE
###prefer these right here, gendr0p1 slightly better than gendr0p2
##source runData.sh "_rec56_HLTCalo" "_gendr0p1_HLTCalo" "10" "10" "17" 1 1
##source runData.sh "_rec56_HLTCalo" "_gendr0p2_HLTCalo" "10" "10" "17" 1 1
##source runData.sh "_rec56_HLTPF" "_gendr0p1_HLTPF" "10" "10" "17" 1 1
##source runData.sh "_rec56_HLTPF" "_gendr0p2_HLTPF" "10" "10" "17" 1 1
#
#
#source runData.sh "0.0eta1.0" "0.0eta1.0" "04" "24" "18" 1 0 1
#source runData.sh "0.0eta1.0" "0.0eta1.0" "04" "24" "18" 1 0 0
#
#echo ""
#echo "done w/ Data."
#echo ""
#
#scp2Serin "output/unfoldDataSpectra/*.pdf"






### USAGE
#"source runMC.sh [dirTag] [MM] [DD] [YY] [runBayes] [runSVD] [useSimpBins]"

echo ""
echo "running MC."
echo ""

mv output/unfoldMCSpectra/* output/unfoldMCSpectra/temp/.

#note; _HLTCalo/HLTPF added to JERS dirs for 10.10.17 sample for convenience of I/O.
# there is no trigger emulation in the MC samples the label is superficial and does not matter.
##prefer these two right here 01/23/18
#source runMC.sh "_gendr0p1_HLTCalo" "10" "10" "17"  1 1
#source runMC.sh "_gendr0p2_HLTCalo" "10" "10" "17" 1 0 1

##working with these two right now
#source runMC.sh "0.0eta1.0" "04" "29" "18" 0 1 1
#source runMC.sh "0.0eta1.0" "04" "29" "18" 0 1 0
source runMC.sh "0.0eta2.0" "04" "29" "18" 1 1 1
source runMC.sh "0.0eta2.0" "04" "29" "18" 1 1 0
#source runMC.sh "1.0eta2.0" "04" "29" "18" 0 1 1
#source runMC.sh "1.0eta2.0" "04" "29" "18" 0 1 0


#source runMC.sh "0.0eta1.0" "04" "27" "18" 0 1 1
#source runMC.sh "0.0eta1.0" "04" "27" "18" 0 1 0

## had issue where the unf reco/genpt cuts were the same as the JER reco/genpt cuts (they go to much lower pt for a reliable resolution measurement in our low pt range, pre 100 GeV)
#source runMC.sh "0.0eta1.0" "04" "24" "18" 1 0 1
#source runMC.sh "0.0eta1.0" "04" "24" "18" 1 0 0

echo ""
echo "done w/ MC."
echo ""

echo ""
echo "moving MC pdf output to local cpu"
echo ""

##this command needs to have input in quotes to work for some reason
scp2Serin "output/unfoldMCSpectra/*.pdf"

return













### OLD MC RUNS NOT PREFERRED
#source runMC.sh "_rec56_HLTCalo" "10" "10" "17" 0 1
#source runMC.sh "_rec64_HLTCalo" "10" "10" "17" 0 1
#source runMC.sh   "_rec56_gen32_HLTCalo" "10" "11" "17" 1 1
#source runMC.sh  "_rec56_gen32_HLTPF" "10" "11" "17" 1 1
#source runMC.sh  "_rec74_gen43_HLTCalo" "10" "11" "17" 1 1
#source runMC.sh  "_rec74_gen43_HLTPF" "10" "11" "17" 1 1


###OLD DATA RUNS NOT PREFERRED
#larger recpt cut did not help, if anything, worse.
#source runData.sh "_rec64_HLTCalo" "_gendr0p1_HLTCalo" "10" "10" "17" 1 1
#source runData.sh "_rec64_HLTCalo" "_gendr0p2_HLTCalo" "10" "10" "17" 1 1
#source runData.sh "_rec64_HLTPF" "_gendr0p1_HLTPF" "10" "10" "17" 1 1
#source runData.sh "_rec64_HLTPF" "_gendr0p2_HLTPF" "10" "10" "17" 1 1

#meh
#source runData.sh "_rec56_HLTCalo" "_rec56_gen32_HLTCalo" "10" "11" "17" 1 1
#source runData.sh "_rec56_HLTPF" "_rec56_gen32_HLTPF" "10" "11" "17" 1 1
#source runData.sh "_rec56_HLTCalo" "_rec74_gen43_HLTCalo" "10" "11" "17" 1 1
#source runData.sh "_rec56_HLTPF" "_rec74_gen43_HLTPF" "10" "11" "17" 1 1





##### ----------------------------------------------------------------------
#echo ""
#echo "compiling SVDUnfoldDataSpectra"
#echo ""
#rooUnfoldCompile SVDUnfoldDataSpectra.C
#
#echo ""
#echo "compiling bayesUnfoldDataSpectra"
#echo ""
#rooUnfoldCompile bayesUnfoldDataSpectra.C
#
#
#echo ""
#echo ""
#echo "------------------------------------------------------------------"
#echo ""
#echo ""
#
#
### ----------------------------------------------------------------------
#echo ""
#echo "compiling SVDUnfoldMCSpectra"
#echo ""
#rooUnfoldCompile SVDUnfoldMCSpectra.C
#
#echo ""
#echo "done compiling unfoldMCSpectra"
#echo ""
#
#echo ""
#echo "compiling bayesUnfoldMCSpectra"
#echo ""
#rooUnfoldCompile bayesUnfoldMCSpectra.C
#
#echo ""
#echo "done compiling unfoldMCSpectra"
#echo ""
#
#
#echo ""
#echo " ------------------------------------------------------------------"
#echo " ------------------------------------------------------------------"
#echo " ------------------------------------------------------------------"
#echo ""

