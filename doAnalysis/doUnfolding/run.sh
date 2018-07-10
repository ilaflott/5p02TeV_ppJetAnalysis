#!/bin/bash


################################################################################################################################
################################################################################################################################


function unfoldData(){

    #echo ""
    #echo "... moving unfoldDataSpectra output ..."
    #echo ""
    
    #mv output/unfoldDataSpectra/*.pdf output/unfoldDataSpectra/temp/.
    #mv output/unfoldDataSpectra/*.root output/unfoldDataSpectra/temp/.
    
    
    echo ""
    echo "// -------------------------------------------------------------------------------------------------- //"
    echo "// -------------------------------------------------------------------------------------------------- //"
    echo "// ------------------------------------- RUNNING DATA UNFOLDING ------------------------------------- //"
    echo "// -------------------------------------------------------------------------------------------------- //"
    echo "// -------------------------------------------------------------------------------------------------- //"
    echo ""
    
#### USAGE
## "source runData.sh [datadirTag] [MCdirTag] [MM] [DD] [YY] [runBayes] [runSVD] [useSimpBins]"
    
#### dates of interesting samples
## 10.10.17 --> APS DNP results
#source runData.sh "_rec56_HLTCalo" "_gendr0p1_HLTCalo" "10" "10" "17" 1 1
#source runData.sh "_rec56_HLTPF" "_gendr0p1_HLTPF" "10" "10" "17" 1 1
    
## 05.09.18 --> figuring out pt ranges post bin width division fix
#source runData.sh "0.0eta2.0_ptHi1588_ptLoRec56"  "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49" "05" "09" "18" 1 1 0
    
## 05.29.18 --> fakes issue fix, unfolding starts to look much better. SVD Unfolding looks good except divector dist suggests weird reg. strength 
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTCalo"  "0.0eta2.0_ptHi2116_ptLoRec49_ptLoGen49_HLTCalo" "05" "29" "18" 1 1 0
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTPF"    "0.0eta2.0_ptHi2116_ptLoRec49_ptLoGen49_HLTPF"   "05" "29" "18" 0 1 0
    
## 06.18.18 --> HLTPF v HLTCalo in ppData does not make divectors better
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTPF"      "0.0eta2.0_ptHi1248_ptLoRec49_ptLoGen49_HLTPF"   "06" "18" "18" 1 1 0
#source runData.sh "0.0eta2.0_ptHi2116_ptLoRec49_HLTCalo"    "0.0eta2.0_ptHi1248_ptLoRec49_ptLoGen49_HLTCalo"   "06" "18" "18" 1 1 0
    
## 06.28.18 --> TH1D->TH1F fix, different jetIDs do not greatly affect divector distribution issue (raw/reco/L2L3recoJetID, HLTPF/HLTCalo)
#source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTPF"      "0.0eta2.0_ptHi1248_ptLoRec49_ptLoGen49_L2L3recoJetID"   "06" "28" "18"   1 1 0
#source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTCalo"    "0.0eta2.0_ptHi1248_ptLoRec49_ptLoGen49_L2L3recoJetID"   "06" "28" "18"   1 1 0
    

    ### 07.02.18 --> TH1D's, eta bins of 0.5, two ppMC JER's w/ diff gendrjt cuts (0.1, 0.2)    
    #source runData.sh    "0.0eta0.5_L2L3recoJetID_HLTPF"      "0.0eta0.5_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "0.5eta1.0_L2L3recoJetID_HLTPF"      "0.5eta1.0_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "1.0eta1.5_L2L3recoJetID_HLTPF"      "1.0eta1.5_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "1.5eta2.0_L2L3recoJetID_HLTPF"      "1.5eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0     

    #source runData.sh    "0.0eta1.0_L2L3recoJetID_HLTPF"      "0.0eta1.0_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0     
    #source runData.sh    "1.0eta2.0_L2L3recoJetID_HLTPF"      "1.0eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0     

    source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTPF"      "0.0eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0     

    dir1=output/unfoldDataSpectra/L2L3recoJetID_HLTPF_gendr0p1_FavoredBinning

    rm -r $dir1
    mkdir $dir1
    mv output/unfoldDataSpectra/*.pdf   $dir1
    mv output/unfoldDataSpectra/*.root  $dir1

    echo ""
    echo "moving data output to local cpu"
    echo ""

    #scp2Serin $dir1
    mv -f $dir1 output/unfoldDataSpectra/temp



    
    #source runData.sh    "0.0eta0.5_L2L3recoJetID_HLTCalo"      "0.0eta0.5_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "0.5eta1.0_L2L3recoJetID_HLTCalo"      "0.5eta1.0_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "1.0eta1.5_L2L3recoJetID_HLTCalo"      "1.0eta1.5_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "1.5eta2.0_L2L3recoJetID_HLTCalo"      "1.5eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0         
    
    #source runData.sh    "0.0eta1.0_L2L3recoJetID_HLTCalo"      "0.0eta1.0_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0     
    #source runData.sh    "1.0eta2.0_L2L3recoJetID_HLTCalo"      "1.0eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0     
    
    #source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTCalo"      "0.0eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0     

    #dir2=output/unfoldDataSpectra/L2L3recoJetID_HLTCalo_gendr0p1_FavoredBinning

    #rm -r $dir2
    #mkdir $dir2
    #mv output/unfoldDataSpectra/*.pdf   $dir2
    #mv output/unfoldDataSpectra/*.root  $dir2
    #
    #echo ""
    #echo "moving data output to local cpu"
    #echo ""
    #
    ##scp2Serin $dir2
    #mv -f $dir2 output/unfoldDataSpectra/temp
    
    #source runData.sh    "0.0eta0.5_L2L3recoJetID_HLTPF"      "0.0eta0.5_gendr0p2_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "0.5eta1.0_L2L3recoJetID_HLTPF"      "0.5eta1.0_gendr0p2_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "1.0eta1.5_L2L3recoJetID_HLTPF"      "1.0eta1.5_gendr0p2_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "1.5eta2.0_L2L3recoJetID_HLTPF"      "1.5eta2.0_gendr0p2_recoJetID"   "07" "02" "18"   1 1 0     
    #source runData.sh    "0.0eta1.0_L2L3recoJetID_HLTPF"      "0.0eta1.0_gendr0p2_recoJetID"   "07" "02" "18"   1 1 0     
    #source runData.sh    "1.0eta2.0_L2L3recoJetID_HLTPF"      "1.0eta2.0_gendr0p2_recoJetID"   "07" "02" "18"   1 1 0     
    #source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTPF"      "0.0eta2.0_gendr0p2_recoJetID"   "07" "02" "18"   1 1 0     
    #dir3=output/unfoldDataSpectra/L2L3recoJetID_HLTPF_gendr0p2
    #rm -r $dir3/
    #mkdir $dir3
    #mv output/unfoldDataSpectra/*.pdf   $dir3
    #mv output/unfoldDataSpectra/*.root  $dir3
    #mv $dir3 output/unfoldDataSpectra/temp
    
    #source runData.sh    "0.0eta0.5_L2L3recoJetID_HLTCalo"      "0.0eta0.5_gendr0p2_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "0.5eta1.0_L2L3recoJetID_HLTCalo"      "0.5eta1.0_gendr0p2_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "1.0eta1.5_L2L3recoJetID_HLTCalo"      "1.0eta1.5_gendr0p2_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "1.5eta2.0_L2L3recoJetID_HLTCalo"      "1.5eta2.0_gendr0p2_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "0.0eta1.0_L2L3recoJetID_HLTCalo"      "0.0eta1.0_gendr0p2_recoJetID"   "07" "02" "18"   1 1 0     
    #source runData.sh    "1.0eta2.0_L2L3recoJetID_HLTCalo"      "1.0eta2.0_gendr0p2_recoJetID"   "07" "02" "18"   1 1 0     
    #source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTCalo"      "0.0eta2.0_gendr0p2_recoJetID"   "07" "02" "18"   1 1 0     
    #dir4=output/unfoldDataSpectra/L2L3recoJetID_HLTCalo_gendr0p2
    #mkdir $dir4
    ##rm $dir4/*
    #mv output/unfoldDataSpectra/*.pdf   $dir4
    #mv output/unfoldDataSpectra/*.root  $dir4
    #mv -f $dir4 output/unfoldDataSpectra/temp

    #source runData.sh    "0.0eta0.5_L2L3recoJetID_HLTPF"      "0.0eta0.5_gendr0p3_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "0.5eta1.0_L2L3recoJetID_HLTPF"      "0.5eta1.0_gendr0p3_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "1.0eta1.5_L2L3recoJetID_HLTPF"      "1.0eta1.5_gendr0p3_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "1.5eta2.0_L2L3recoJetID_HLTPF"      "1.5eta2.0_gendr0p3_recoJetID"   "07" "02" "18"   1 1 0     
    #source runData.sh    "0.0eta1.0_L2L3recoJetID_HLTPF"      "0.0eta1.0_gendr0p3_recoJetID"   "07" "02" "18"   1 1 0     
    #source runData.sh    "1.0eta2.0_L2L3recoJetID_HLTPF"      "1.0eta2.0_gendr0p3_recoJetID"   "07" "02" "18"   1 1 0     
    #source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTPF"      "0.0eta2.0_gendr0p3_recoJetID"   "07" "02" "18"   1 1 0     
    #dir5=output/unfoldDataSpectra/L2L3recoJetID_HLTPF_gendr0p3
    #rm -r $dir5/
    #mkdir $dir5
    #mv output/unfoldDataSpectra/*.pdf   $dir5
    #mv output/unfoldDataSpectra/*.root  $dir5
    #mv $dir5 output/unfoldDataSpectra/temp
    
    #source runData.sh    "0.0eta0.5_L2L3recoJetID_HLTCalo"      "0.0eta0.5_gendr0p3_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "0.5eta1.0_L2L3recoJetID_HLTCalo"      "0.5eta1.0_gendr0p3_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "1.0eta1.5_L2L3recoJetID_HLTCalo"      "1.0eta1.5_gendr0p3_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "1.5eta2.0_L2L3recoJetID_HLTCalo"      "1.5eta2.0_gendr0p3_recoJetID"   "07" "02" "18"   1 1 0 
    #source runData.sh    "0.0eta1.0_L2L3recoJetID_HLTCalo"      "0.0eta1.0_gendr0p3_recoJetID"   "07" "02" "18"   1 1 0     
    #source runData.sh    "1.0eta2.0_L2L3recoJetID_HLTCalo"      "1.0eta2.0_gendr0p3_recoJetID"   "07" "02" "18"   1 1 0     
    #source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTCalo"      "0.0eta2.0_gendr0p3_recoJetID"   "07" "02" "18"   1 1 0     
    #dir6=output/unfoldDataSpectra/L2L3recoJetID_HLTCalo_gendr0p3
    #mkdir $dir6
    ##rm $dir6/*
    #mv output/unfoldDataSpectra/*.pdf   $dir6
    #mv output/unfoldDataSpectra/*.root  $dir6
    #mv -f $dir6 output/unfoldDataSpectra/temp


#    source runData.sh    "0.0eta0.5_L2L3recoJetID_HLTPF"      "0.0eta0.5_gendr0p4_recoJetID"   "07" "02" "18"   1 1 0 
#    source runData.sh    "0.5eta1.0_L2L3recoJetID_HLTPF"      "0.5eta1.0_gendr0p4_recoJetID"   "07" "02" "18"   1 1 0 
#    source runData.sh    "1.0eta1.5_L2L3recoJetID_HLTPF"      "1.0eta1.5_gendr0p4_recoJetID"   "07" "02" "18"   1 1 0 
#    source runData.sh    "1.5eta2.0_L2L3recoJetID_HLTPF"      "1.5eta2.0_gendr0p4_recoJetID"   "07" "02" "18"   1 1 0     
#    source runData.sh    "0.0eta1.0_L2L3recoJetID_HLTPF"      "0.0eta1.0_gendr0p4_recoJetID"   "07" "02" "18"   1 1 0     
#    source runData.sh    "1.0eta2.0_L2L3recoJetID_HLTPF"      "1.0eta2.0_gendr0p4_recoJetID"   "07" "02" "18"   1 1 0     
#    source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTPF"      "0.0eta2.0_gendr0p4_recoJetID"   "07" "02" "18"   1 1 0     
#    dir7=output/unfoldDataSpectra/L2L3recoJetID_HLTPF_gendr0p4
#    mkdir $dir7
#    rm $dir7/*
#    mv output/unfoldDataSpectra/*.pdf   $dir7
#    mv output/unfoldDataSpectra/*.root  $dir7
#    mv $dir7 output/unfoldDataSpectra/temp
#    
#    source runData.sh    "0.0eta0.5_L2L3recoJetID_HLTCalo"      "0.0eta0.5_gendr0p4_recoJetID"   "07" "02" "18"   1 1 0 
#    source runData.sh    "0.5eta1.0_L2L3recoJetID_HLTCalo"      "0.5eta1.0_gendr0p4_recoJetID"   "07" "02" "18"   1 1 0 
#    source runData.sh    "1.0eta1.5_L2L3recoJetID_HLTCalo"      "1.0eta1.5_gendr0p4_recoJetID"   "07" "02" "18"   1 1 0 
#    source runData.sh    "1.5eta2.0_L2L3recoJetID_HLTCalo"      "1.5eta2.0_gendr0p4_recoJetID"   "07" "02" "18"   1 1 0 
#    source runData.sh    "0.0eta1.0_L2L3recoJetID_HLTCalo"      "0.0eta1.0_gendr0p4_recoJetID"   "07" "02" "18"   1 1 0     
#    source runData.sh    "1.0eta2.0_L2L3recoJetID_HLTCalo"      "1.0eta2.0_gendr0p4_recoJetID"   "07" "02" "18"   1 1 0     
#    source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTCalo"      "0.0eta2.0_gendr0p4_recoJetID"   "07" "02" "18"   1 1 0     
#    dir8=output/unfoldDataSpectra/L2L3recoJetID_HLTCalo_gendr0p4
#    mkdir $dir8
#    rm $dir8/*
#    mv output/unfoldDataSpectra/*.pdf   $dir8
#    mv output/unfoldDataSpectra/*.root  $dir8
#    mv $dir8 output/unfoldDataSpectra/temp




    #scp2Serin $dir1
    #rm -r $dir1
    ##mv $dir1 output/unfoldDataSpectra/temp    
    #scp2Serin $dir2
    #rm -r $dir2
    ##mv $dir2 output/unfoldDataSpectra/temp    
    #scp2Serin $dir3
    #rm -r $dir3
    ##mv $dir3 output/unfoldDataSpectra/temp    
    #scp2Serin $dir4
    #rm -r $dir4
    ##mv $dir4 output/unfoldDataSpectra/temp    
    #scp2Serin $dir5
    #rm -r $dir5
    ##mv $dir4 output/unfoldDataSpectra/temp    
    #scp2Serin $dir6
    #rm -r $dir6
    ##mv $dir4 output/unfoldDataSpectra/temp    



### idea; always give name to run.sh script for dir to conglomerate all results i'm running?


### default way i generally deal with output
    #echo ""
    #echo "moving data output to local cpu"
    #echo ""
    #scp2Serin "output/unfoldDataSpectra/*.pdf"
    #scp2Serin "output/unfoldDataSpectra/*.root"
    
    
    echo ""
    echo "// ---------------------------------------------------------------------------------------------- //"
    echo "// ---------------------------------------------------------------------------------------------- //"
    echo "// ------------------------------------- END DATA UNFOLDING ------------------------------------- //"
    echo "// ---------------------------------------------------------------------------------------------- //"
    echo "// ---------------------------------------------------------------------------------------------- //"
    echo ""
    
    return
    
}


################################################################################################################################
################################################################################################################################




function unfoldMC(){
    
    #echo ""
    #echo "... moving unfoldMCSpectra output ..."
    #echo ""
    
    #mv output/unfoldMCSpectra/*.pdf  output/unfoldMCSpectra/temp/.
    #mv output/unfoldMCSpectra/*.root output/unfoldMCSpectra/temp/.
    
    echo ""
    echo "// ------------------------------------------------------------------------------------------------ //"
    echo "// ------------------------------------------------------------------------------------------------ //"
    echo "// ------------------------------------- RUNNING MC UNFOLDING ------------------------------------- //"
    echo "// ------------------------------------------------------------------------------------------------ //"
    echo "// ------------------------------------------------------------------------------------------------ //"
    echo ""
    
### USAGE
#"source runMC.sh [dirTag] [MM] [DD] [YY] [runBayes] [runSVD] [useSimpBins]"
    
#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49"   "05" "09" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49"   "05" "09" "18" 1 1 0 
    
#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec56_ptLoGen56"  "05" "09" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec56_ptLoGen56"  "05" "09" "18" 1 1 0 
    
#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec32_ptLoGen32"  "05" "09" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec32_ptLoGen32"  "05" "09" "18" 1 1 0 
    
#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec21_ptLoGen21"  "05" "09" "18" 1 1 1
#source runMC.sh "0.0eta2.0_ptHi1032_ptLoRec21_ptLoGen21"  "05" "09" "18" 1 1 0 

    #source runMC.sh   "0.0eta0.5_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0 
    #source runMC.sh   "0.5eta1.0_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0 
    #source runMC.sh   "1.0eta1.5_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0 
    #source runMC.sh   "1.5eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0         

    #source runMC.sh   "0.0eta1.0_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0     
    #source runMC.sh   "1.0eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0     

    source runMC.sh   "0.0eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   1 1 0     

    MCdir1=output/unfoldMCSpectra/L2L3recoJetID_MCClosure_gendr0p1_FavoredBinning

    rm -r $MCdir1
    mkdir $MCdir1
    mv output/unfoldMCSpectra/*.pdf   $MCdir1
    mv output/unfoldMCSpectra/*.root  $MCdir1

    echo ""
    echo "moving MC output to local cpu"
    echo ""

    #scp2Serin $MCdir1
    mv -f $MCdir1 output/unfoldMCSpectra/temp    
    
    #scp2Serin "output/unfoldMCSpectra/*.pdf"
    #scp2Serin "output/unfoldMCSpectra/*.root"
    
    
    echo ""
    echo "// -------------------------------------------------------------------------------------------- //"
    echo "// -------------------------------------------------------------------------------------------- //"
    echo "// ------------------------------------- END MC UNFOLDING ------------------------------------- //"
    echo "// -------------------------------------------------------------------------------------------- //"
    echo "// -------------------------------------------------------------------------------------------- //"
    echo ""
    
    return

}

################################################################################################################################
################################################################################################################################










echo ""
echo "... making rooUnfold ..."
echo ""

#make_rooUnfold

unfoldData

unfoldMC

return





