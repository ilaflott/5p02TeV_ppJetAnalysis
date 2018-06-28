#!/bin/bash






make_rooUnfold



#### USAGE
## "source runData.sh [datadirTag] [MCdirTag] [MM] [DD] [YY] [runBayes] [runSVD] [useSimpBins]"
#



mv output/unfoldDataSpectra/*.pdf output/unfoldDataSpectra/temp/.
mv output/unfoldDataSpectra/*.root output/unfoldDataSpectra/temp/.

echo ""
echo "running Data...."
echo ""


echo ""
echo "// ------------------------------------------------------------------------------------ //"
echo "// ------------------------------------------------------------------------------------ //"
echo "// ------------------------------------- BEGIN ---------------------------------------- //"
echo "// ------------------------------------------------------------------------------------ //"
echo "// ------------------------------------------------------------------------------------ //"
echo ""

####UPDATE runData.sh ALA runMC.sh WHEN YOU GET A CHANCE
###prefer these right here, gendr0p1 slightly better than gendr0p2
##source runData.sh "_rec56_HLTCalo" "_gendr0p1_HLTCalo" "10" "10" "17" 1 1
##source runData.sh "_rec56_HLTCalo" "_gendr0p2_HLTCalo" "10" "10" "17" 1 1
##source runData.sh "_rec56_HLTPF" "_gendr0p1_HLTPF" "10" "10" "17" 1 1
##source runData.sh "_rec56_HLTPF" "_gendr0p2_HLTPF" "10" "10" "17" 1 1


#source runData.sh "0.0eta2.0" "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49" "05" "01" "18" 1 1 1
#source runData.sh "0.0eta2.0" "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49" "05" "01" "18" 1 1 0

#source runData.sh "0.0eta2.0_ptHi1588_ptLoRec56" "0.0eta2.0_ptHi1032_ptLoRec56_ptLoGen56" "05" "09" "18" 1 0 0
#source runData.sh "0.0eta2.0_ptHi1032_ptLoRec56" "0.0eta2.0_ptHi1032_ptLoRec56_ptLoGen56" "05" "09" "18" 1 0 0
#source runData.sh "0.0eta2.0_ptHi1032_ptLoRec32" "0.0eta2.0_ptHi1032_ptLoRec32_ptLoGen32" "05" "09" "18" 1 0 0
#source runData.sh "0.0eta2.0_ptHi967_ptLoRec49" "0.0eta2.0_ptHi1032_ptLoRec32_ptLoGen32" "05" "09" "18" 1 0 0
#source runData.sh "0.0eta2.0_ptHi967_ptLoRec49"  "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49" "05" "09" "18" 1 0 0
#source runData.sh "0.0eta2.0_ptHi967_ptLoRec49"  "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49" "05" "09" "18" 1 1 1
#source runData.sh "0.0eta2.0_ptHi967_ptLoRec49"  "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49" "05" "09" "18" 1 0 0
#source runData.sh "0.0eta2.0_ptHi967_ptLoRec49" "0.0eta2.0_ptHi1032_ptLoRec56_ptLoGen56" "05" "09" "18" 1 1 0
#source runData.sh "0.0eta2.0_ptHi1588_ptLoRec56"  "0.0eta2.0_ptHi1032_ptLoRec32_ptLoGen32" "05" "09" "18" 1 0 0

####best unfolding (between 64-967 GeV) of 05/09/18 samples. 
#source runData.sh "0.0eta2.0_ptHi1588_ptLoRec56"  "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49" "05" "09" "18" 1 1 0


#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTCalo"  "0.0eta2.0_ptHi2116_ptLoRec49_ptLoGen49_HLTCalo" "05" "29" "18" 1 1 0
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTPF"    "0.0eta2.0_ptHi2116_ptLoRec49_ptLoGen49_HLTPF"   "05" "29" "18" 0 1 0

#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTCalo"    "0.0eta2.0_ptHi2116_ptLoRec49_ptLoGen49_HLTCalo"   "05" "29" "18" 1 1 1
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTPF"    "0.0eta2.0_ptHi2116_ptLoRec49_ptLoGen49_HLTPF"   "05" "29" "18" 1 1 1

#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTPF"    "0.0eta2.0_ptHi1248_ptLoRec49_ptLoGen49_HLTPF"   "06" "18" "18" 1 1 0
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTPF"    "0.0eta2.0_ptHi1172_ptLoRec49_ptLoGen49_HLTPF"   "06" "18" "18" 1 1 0
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTPF"    "0.0eta2.0_ptHi1101_ptLoRec49_ptLoGen49_HLTPF"   "06" "18" "18" 1 1 0
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTPF"    "0.0eta2.0_ptHi1032_ptLoRec49_ptLoGen49_HLTPF"   "06" "18" "18" 1 1 0
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTPF"     "0.0eta2.0_ptHi967_ptLoRec64_ptLoGen64_HLTPF"   "06" "18" "18" 1 1 0


#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTCalo"    "0.0eta2.0_ptHi1248_ptLoRec49_ptLoGen49_HLTCalo"   "06" "18" "18" 1 1 0
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTCalo"    "0.0eta2.0_ptHi1172_ptLoRec49_ptLoGen49_HLTCalo"   "06" "18" "18" 1 1 0
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTCalo"    "0.0eta2.0_ptHi1101_ptLoRec49_ptLoGen49_HLTCalo"   "06" "18" "18" 1 1 0
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTCalo"    "0.0eta2.0_ptHi1032_ptLoRec49_ptLoGen49_HLTCalo"   "06" "18" "18" 1 1 0
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTCalo"    "0.0eta2.0_ptHi967_ptLoRec64_ptLoGen64_HLTCalo"   "06" "18" "18" 1 1 0


source runData.sh "0.0eta2.0"    "0.0eta2.0"   "06" "21" "18" 1 1 0

#source runData.sh "0.0eta2.0_ptHi1588_ptLoRec56"  "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49" "05" "09" "18" 0 1 0


echo ""
echo "// ------------------------------------------------------------------------------------ //"
echo "// ------------------------------------------------------------------------------------ //"
echo "// -------------------------------------- END ----------------------------------------- //"
echo "// ------------------------------------------------------------------------------------ //"
echo "// ------------------------------------------------------------------------------------ //"
echo ""

echo ""
echo "done w/ Data."
echo ""




echo ""
echo "moving data pdf output to local cpu"
echo ""
scp2Serin "output/unfoldDataSpectra/*.pdf"
scp2Serin "output/unfoldDataSpectra/*.root"






return















### USAGE
#"source runMC.sh [dirTag] [MM] [DD] [YY] [runBayes] [runSVD] [useSimpBins]"

echo ""
echo "running MC."
echo ""

mv output/unfoldMCSpectra/*.pdf output/unfoldMCSpectra/temp/.
mv output/unfoldMCSpectra/*.root output/unfoldMCSpectra/temp/.

source runMC.sh "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49"   "05" "09" "18" 1 1 1
source runMC.sh "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49"   "05" "09" "18" 1 1 0

#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec56_ptLoGen56"  "05" "09" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec56_ptLoGen56"  "05" "09" "18" 1 1 0

#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec32_ptLoGen32"  "05" "09" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec32_ptLoGen32"  "05" "09" "18" 1 1 0

#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec21_ptLoGen21"  "05" "09" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec21_ptLoGen21"  "05" "09" "18" 1 1 0


echo ""
echo "done w/ MC."
echo ""




#echo ""
#echo "moving MC pdf output to local cpu"
#echo ""
###this command needs to have input in quotes to work for some reason
#scp2Serin "output/unfoldMCSpectra/*.pdf"

echo ""
echo "moving MC + data pdf output to local cpu"
echo ""
scp2Serin "output/unfold*Spectra/*.pdf"

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









############################################
#### WINNER WINNER CHICKEN DINNER 49/49 ####
############################################

#### Symm gen/recpt cuts, 49/49 wins, both binnings
#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec37_ptLoGen37"   "05" "04" "18" 1 1 0
#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec43_ptLoGen43"   "05" "04" "18" 1 1 0
source runMC.sh "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49"   "05" "04" "18" 0 1 0

#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec37_ptLoGen37"   "05" "04" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec43_ptLoGen43"   "05" "04" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49"   "05" "04" "18" 1 1 1


##### Asymm gen/recpt cuts w/ genpt cut lower than recpt cut, 56/49 wins, both binnings
#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec56_ptLoGen37"   "05" "04" "18" 1 1 0
#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec56_ptLoGen43"   "05" "04" "18" 1 1 0
source runMC.sh "0.0eta2.0_ptHi967_ptLoRec56_ptLoGen49"   "05" "04" "18" 0 1 0

#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec56_ptLoGen37"   "05" "04" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec56_ptLoGen43"   "05" "04" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec56_ptLoGen49"   "05" "04" "18" 1 1 1


##### Asymm gen/pt cuts w/ recpt cut lower than genpt cut, approx tie for 43/49 + 49/56, YYY bins
#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec37_ptLoGen43"   "05" "04" "18" 1 1 0
source runMC.sh "0.0eta2.0_ptHi967_ptLoRec43_ptLoGen49"   "05" "04" "18" 0 1 0
source runMC.sh "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen56"   "05" "04" "18" 0 1 0
#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec56_ptLoGen64"   "05" "04" "18" 1 1 0

#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec37_ptLoGen43"   "05" "04" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec43_ptLoGen49"   "05" "04" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen56"   "05" "04" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec56_ptLoGen64"   "05" "04" "18" 1 1 1

#source runMC.sh "0.0eta2.0" "05" "01" "18" 1 1 1
#source runMC.sh "0.0eta2.0" "05" "01" "18" 1 1 0

#source runMC.sh "0.0eta2.0_ptHi1588" "05" "02" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi1588" "05" "02" "18" 1 1 0
#source runMC.sh "0.0eta2.0_ptHi1497" "05" "02" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi1497" "05" "02" "18" 1 1 0
#source runMC.sh "0.0eta2.0_ptHi1410" "05" "02" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi1410" "05" "02" "18" 1 1 0

#source runMC.sh "0.0eta2.0_ptHi1327" "05" "02" "18" 1 0 1
#source runMC.sh "0.0eta2.0_ptHi1327" "05" "02" "18" 1 0 0
#source runMC.sh "0.0eta2.0_ptHi1327" "05" "02" "18" 0 1 1
#source runMC.sh "0.0eta2.0_ptHi1327" "05" "02" "18" 0 1 0

#source runMC.sh "0.0eta2.0_ptHi1248" "05" "02" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi1248" "05" "02" "18" 1 1 0

#source runMC.sh "0.0eta1.0" "04" "29" "18" 0 1 1
#source runMC.sh "0.0eta1.0" "04" "29" "18" 0 1 0
#source runMC.sh "0.0eta2.0" "04" "29" "18" 1 1 1
#source runMC.sh "0.0eta2.0" "04" "29" "18" 1 1 0
#source runMC.sh "1.0eta2.0" "04" "29" "18" 0 1 1
#source runMC.sh "1.0eta2.0" "04" "29" "18" 0 1 0


#source runMC.sh "0.0eta1.0" "04" "27" "18" 0 1 1
#source runMC.sh "0.0eta1.0" "04" "27" "18" 0 1 0

## had issue where the unf reco/genpt cuts were the same as the JER reco/genpt cuts (they go to much lower pt for a reliable resolution measurement in our low pt range, pre 100 GeV)
#source runMC.sh "0.0eta1.0" "04" "24" "18" 1 0 1
#source runMC.sh "0.0eta1.0" "04" "24" "18" 1 0 0

#note; _HLTCalo/HLTPF added to JERS dirs for 10.10.17 sample for convenience of I/O.
# there is no trigger emulation in the MC samples the label is superficial and does not matter.
##prefer these two right here 01/23/18
#source runMC.sh "_gendr0p1_HLTCalo" "10" "10" "17"  1 1
#source runMC.sh "_gendr0p2_HLTCalo" "10" "10" "17" 1 0 1
