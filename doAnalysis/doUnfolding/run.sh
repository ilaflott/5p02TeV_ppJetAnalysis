#!/bin/bash



###unfolding with PYTHIA8
#source run_unfoldDataSpectra_etabin.sh
#source   run_unfoldMCSpectra_etabin.sh

####smeared NLO unfolding
source run_unfoldDataSpectra_wNLO_etabin.sh
#source   run_unfoldMCSpectra_wNLO_etabin.sh


































#### dates of interesting samples
## 10.10.17 --> APS DNP results
#source runData.sh "_rec56_HLTCalo" "_gendr0p1_HLTCalo" "10" "10" "17" 1 1
#source runData.sh "_rec56_HLTPF" "_gendr0p1_HLTPF" "10" "10" "17" 1 1
    
## 05.09.18 --> figuring out pt ranges post bin width division fix
#source runData.sh "0.0eta2.0_ptHi1588_ptLoRec56"  "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49" "05" "09" "18" $doBayes $doSVD $useSimpBins
    
## 05.29.18 --> fakes issue fix, unfolding starts to look much better. SVD Unfolding looks good except divector dist suggests weird reg. strength 
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTCalo"  "0.0eta2.0_ptHi2116_ptLoRec49_ptLoGen49_HLTCalo" "05" "29" "18" $doBayes $doSVD $useSimpBins
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTPF"    "0.0eta2.0_ptHi2116_ptLoRec49_ptLoGen49_HLTPF"   "05" "29" "18" 0 1 0
    
## 06.18.18 --> HLTPF v HLTCalo in ppData does not make divectors better
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTPF"      "0.0eta2.0_ptHi1248_ptLoRec49_ptLoGen49_HLTPF"   "06" "18" "18" $doBayes $doSVD $useSimpBins
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTCalo"    "0.0eta2.0_ptHi1248_ptLoRec49_ptLoGen49_HLTCalo"   "06" "18" "18" $doBayes $doSVD $useSimpBins
    
## 06.28.18 --> TH1D->TH1F fix, different jetIDs do not greatly affect divector distribution issue (raw/reco/L2L3recoJetID, HLTPF/HLTCalo)
#source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTPF"      "0.0eta2.0_ptHi1248_ptLoRec49_ptLoGen49_L2L3recoJetID"   "06" "28" "18"   $doBayes $doSVD $useSimpBins
#source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTCalo"    "0.0eta2.0_ptHi1248_ptLoRec49_ptLoGen49_L2L3recoJetID"   "06" "28" "18"   $doBayes $doSVD $useSimpBins
    

    #source runData.sh    "0.0eta0.5_L2L3recoJetID_HLTPF"      "0.0eta0.5_gendr0p2_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
    #source runData.sh    "0.5eta1.0_L2L3recoJetID_HLTPF"      "0.5eta1.0_gendr0p2_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
    #source runData.sh    "1.0eta1.5_L2L3recoJetID_HLTPF"      "1.0eta1.5_gendr0p2_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
    #source runData.sh    "1.5eta2.0_L2L3recoJetID_HLTPF"      "1.5eta2.0_gendr0p2_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
    #source runData.sh    "0.0eta1.0_L2L3recoJetID_HLTPF"      "0.0eta1.0_gendr0p2_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
    #source runData.sh    "1.0eta2.0_L2L3recoJetID_HLTPF"      "1.0eta2.0_gendr0p2_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
    #source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTPF"      "0.0eta2.0_gendr0p2_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
    #dir3=output/unfoldDataSpectra/L2L3recoJetID_HLTPF_gendr0p2
    #rm -r $dir3/
    #mkdir $dir3
    #mv output/unfoldDataSpectra/*.pdf   $dir3
    #mv output/unfoldDataSpectra/*.root  $dir3
    #mv $dir3 output/unfoldDataSpectra/temp
    
    #source runData.sh    "0.0eta0.5_L2L3recoJetID_HLTCalo"      "0.0eta0.5_gendr0p2_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
    #source runData.sh    "0.5eta1.0_L2L3recoJetID_HLTCalo"      "0.5eta1.0_gendr0p2_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
    #source runData.sh    "1.0eta1.5_L2L3recoJetID_HLTCalo"      "1.0eta1.5_gendr0p2_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
    #source runData.sh    "1.5eta2.0_L2L3recoJetID_HLTCalo"      "1.5eta2.0_gendr0p2_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
    #source runData.sh    "0.0eta1.0_L2L3recoJetID_HLTCalo"      "0.0eta1.0_gendr0p2_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
    #source runData.sh    "1.0eta2.0_L2L3recoJetID_HLTCalo"      "1.0eta2.0_gendr0p2_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
    #source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTCalo"      "0.0eta2.0_gendr0p2_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
    #dir4=output/unfoldDataSpectra/L2L3recoJetID_HLTCalo_gendr0p2
    #mkdir $dir4
    ##rm $dir4/*
    #mv output/unfoldDataSpectra/*.pdf   $dir4
    #mv output/unfoldDataSpectra/*.root  $dir4
    #mv -f $dir4 output/unfoldDataSpectra/temp

    #source runData.sh    "0.0eta0.5_L2L3recoJetID_HLTPF"      "0.0eta0.5_gendr0p3_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
    #source runData.sh    "0.5eta1.0_L2L3recoJetID_HLTPF"      "0.5eta1.0_gendr0p3_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
    #source runData.sh    "1.0eta1.5_L2L3recoJetID_HLTPF"      "1.0eta1.5_gendr0p3_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
    #source runData.sh    "1.5eta2.0_L2L3recoJetID_HLTPF"      "1.5eta2.0_gendr0p3_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
    #source runData.sh    "0.0eta1.0_L2L3recoJetID_HLTPF"      "0.0eta1.0_gendr0p3_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
    #source runData.sh    "1.0eta2.0_L2L3recoJetID_HLTPF"      "1.0eta2.0_gendr0p3_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
    #source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTPF"      "0.0eta2.0_gendr0p3_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
    #dir5=output/unfoldDataSpectra/L2L3recoJetID_HLTPF_gendr0p3
    #rm -r $dir5/
    #mkdir $dir5
    #mv output/unfoldDataSpectra/*.pdf   $dir5
    #mv output/unfoldDataSpectra/*.root  $dir5
    #mv $dir5 output/unfoldDataSpectra/temp
    
    #source runData.sh    "0.0eta0.5_L2L3recoJetID_HLTCalo"      "0.0eta0.5_gendr0p3_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
    #source runData.sh    "0.5eta1.0_L2L3recoJetID_HLTCalo"      "0.5eta1.0_gendr0p3_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
    #source runData.sh    "1.0eta1.5_L2L3recoJetID_HLTCalo"      "1.0eta1.5_gendr0p3_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
    #source runData.sh    "1.5eta2.0_L2L3recoJetID_HLTCalo"      "1.5eta2.0_gendr0p3_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
    #source runData.sh    "0.0eta1.0_L2L3recoJetID_HLTCalo"      "0.0eta1.0_gendr0p3_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
    #source runData.sh    "1.0eta2.0_L2L3recoJetID_HLTCalo"      "1.0eta2.0_gendr0p3_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
    #source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTCalo"      "0.0eta2.0_gendr0p3_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
    #dir6=output/unfoldDataSpectra/L2L3recoJetID_HLTCalo_gendr0p3
    #mkdir $dir6
    ##rm $dir6/*
    #mv output/unfoldDataSpectra/*.pdf   $dir6
    #mv output/unfoldDataSpectra/*.root  $dir6
    #mv -f $dir6 output/unfoldDataSpectra/temp


#    source runData.sh    "0.0eta0.5_L2L3recoJetID_HLTPF"      "0.0eta0.5_gendr0p4_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
#    source runData.sh    "0.5eta1.0_L2L3recoJetID_HLTPF"      "0.5eta1.0_gendr0p4_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
#    source runData.sh    "1.0eta1.5_L2L3recoJetID_HLTPF"      "1.0eta1.5_gendr0p4_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
#    source runData.sh    "1.5eta2.0_L2L3recoJetID_HLTPF"      "1.5eta2.0_gendr0p4_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
#    source runData.sh    "0.0eta1.0_L2L3recoJetID_HLTPF"      "0.0eta1.0_gendr0p4_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
#    source runData.sh    "1.0eta2.0_L2L3recoJetID_HLTPF"      "1.0eta2.0_gendr0p4_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
#    source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTPF"      "0.0eta2.0_gendr0p4_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
#    dir7=output/unfoldDataSpectra/L2L3recoJetID_HLTPF_gendr0p4
#    mkdir $dir7
#    rm $dir7/*
#    mv output/unfoldDataSpectra/*.pdf   $dir7
#    mv output/unfoldDataSpectra/*.root  $dir7
#    mv $dir7 output/unfoldDataSpectra/temp
#    
#    source runData.sh    "0.0eta0.5_L2L3recoJetID_HLTCalo"      "0.0eta0.5_gendr0p4_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
#    source runData.sh    "0.5eta1.0_L2L3recoJetID_HLTCalo"      "0.5eta1.0_gendr0p4_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
#    source runData.sh    "1.0eta1.5_L2L3recoJetID_HLTCalo"      "1.0eta1.5_gendr0p4_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
#    source runData.sh    "1.5eta2.0_L2L3recoJetID_HLTCalo"      "1.5eta2.0_gendr0p4_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
#    source runData.sh    "0.0eta1.0_L2L3recoJetID_HLTCalo"      "0.0eta1.0_gendr0p4_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
#    source runData.sh    "1.0eta2.0_L2L3recoJetID_HLTCalo"      "1.0eta2.0_gendr0p4_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
#    source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTCalo"      "0.0eta2.0_gendr0p4_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
#    dir8=output/unfoldDataSpectra/L2L3recoJetID_HLTCalo_gendr0p4
#    mkdir $dir8
#    rm $dir8/*
#    mv output/unfoldDataSpectra/*.pdf   $dir8
#    mv output/unfoldDataSpectra/*.root  $dir8
#    mv $dir8 output/unfoldDataSpectra/temp




#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49"   "05" "09" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49"   "05" "09" "18" $doBayes $doSVD $useSimpBins 
    
#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec56_ptLoGen56"  "05" "09" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec56_ptLoGen56"  "05" "09" "18" $doBayes $doSVD $useSimpBins 
    
#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec32_ptLoGen32"  "05" "09" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec32_ptLoGen32"  "05" "09" "18" $doBayes $doSVD $useSimpBins 
    
#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec21_ptLoGen21"  "05" "09" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec21_ptLoGen21"  "05" "09" "18" $doBayes $doSVD $useSimpBins 

