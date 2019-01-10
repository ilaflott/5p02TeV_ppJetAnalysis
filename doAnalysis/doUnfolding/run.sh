#!/bin/bash


#Data and SVD Unfolding settings
dataHLTPF=1
dataHLTCalo=0
doBayes=1
doSVD=1
useSimpBins=0
makeRooUnf=0

###LO Unf ###
#MC
doLOMCUnfolding=0
LOMCBayeskIter=6
LOMCSVDkReg=6
#Data
doDataUnfoldingwLO=1
dataLOBayeskIter=4
dataLOSVDkReg=6

###NLO Unf###
useNPCorr=1
#MC
doNLOMCUnfolding=0
NLOMCBayeskIter=6
NLOMCSVDkReg=6
#Data
doDataUnfoldingwNLO=0
dataNLOBayeskIter=6
dataNLOSVDkReg=6








################################################################################################################################
FINALOUTDIR=/cms/heavyion/ilaflott/temp/4scp
DATAOUTDIR=output/unfoldDataSpectra
HLTPFdir=L2L3recoJetID_HLTPF_w_Py8LOMC_NLOMC_NLONPMC
HLTCalodir=L2L3recoJetID_HLTCalo_w_Py8LOMC_NLOMC_NLONPMC
################################################################################################################################
function unfoldData(){

    echo ""
    echo "... cleaning up previous unfoldDataSpectra output ..."
    echo ""
    
    mv ${DATAOUTDIR}/*.pdf  ${DATAOUTDIR}/temp/.
    mv ${DATAOUTDIR}/*.root ${DATAOUTDIR}/temp/.
    
    
    echo ""
    echo "// -------------------------------------------------------------------------------------------------- //"
    echo "// -------------------------------------------------------------------------------------------------- //"
    echo "// ------------------------------------- RUNNING DATA UNFOLDING ------------------------------------- //"
    echo "// -------------------------------------------------------------------------------------------------- //"
    echo "// -------------------------------------------------------------------------------------------------- //"
    echo ""

    
    if [[ $# -eq 1 ]]
    then
	#HLTPFdir=L2L3recoJetID_HLTPF_${1}
	#HLTCalodir=L2L3recoJetID_HLTCalo_${1}
	HLTPFdir=${HLTPFdir}_${1}
	HLTCalodir=${HLTCalodir}_${1}
    #else
	#HLTPFdir=L2L3recoJetID_HLTPF
	#HLTCalodir=L2L3recoJetID_HLTCalo
    fi
    
    
    
    
    if [[ dataHLTPF -eq 1 ]]
    then
	
	if [[ doDataUnfoldingwLO -eq 1 ]]
	then
        ### 07.02.18 --> TH1D's, eta bins of 0.5, two ppMC JER's w/ diff gendrjt cuts (0.1, 0.2)    
        #source runData.sh    "0.0eta0.5_L2L3recoJetID_HLTPF"      "0.0eta0.5_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
        #source runData.sh    "0.5eta1.0_L2L3recoJetID_HLTPF"      "0.5eta1.0_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
        #source runData.sh    "1.0eta1.5_L2L3recoJetID_HLTPF"      "1.0eta1.5_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
        #source runData.sh    "1.5eta2.0_L2L3recoJetID_HLTPF"      "1.5eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
    	    
        #source runData.sh    "0.0eta1.0_L2L3recoJetID_HLTPF"      "0.0eta1.0_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
        #source runData.sh    "1.0eta2.0_L2L3recoJetID_HLTPF"      "1.0eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
	    
	#source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTPF"      "0.0eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins         
	#source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTPF"      "0.0eta2.0_gendr0p1_genPtLo43_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins  
	#source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTPF"      "0.0eta2.0_gendr0p4_genPtLo43_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins  
	    
	    #source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTPF"      "0.0eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   $useSimpBins $doBayes $doSVD $dataLOSVDkReg $dataLOBayeskIter
	    #source runData.sh    "0.0eta2.0_trgPtNoRes"      "0.0eta2.0_unf"   "11" "14" "18"   $useSimpBins $doBayes $doSVD $dataLOSVDkReg $dataLOBayeskIter
	    source runData.sh    "0.0eta2.0"      "0.0eta2.0_PY8unf"   "01" "06" "19"   $useSimpBins $doBayes $doSVD $dataLOSVDkReg $dataLOBayeskIter
	fi
	
	if [[ doDataUnfoldingwNLO -eq 1 ]]
	then
	    #source runData_wNLO.sh "0.0eta2.0_L2L3recoJetID_HLTPF"  "0.0eta2.0_gendr0p1_recoJetID" "07" "02" "18"   $useSimpBins $doBayes $doSVD $dataNLOSVDkReg $dataNLOBayeskIter
	    #source runData_wNLO.sh "0.0eta2.0_trgPtNoRes"  "0.0eta2.0_unf" "11" "14" "18"   $useSimpBins $doBayes $doSVD $dataNLOSVDkReg $dataNLOBayeskIter
	    source runData_wNLO.sh "0.0eta2.0_trgPtNoRes"  "11" "14" "18"   $useSimpBins $doBayes $doSVD $dataNLOSVDkReg $dataNLOBayeskIter $useNPCorr
	fi
	
	#echo ""
	#echo "cleaning up previous DATA output in ${DATAOUTDIR}/${HLTPFdir}"
	#echo ""
	#
	#rm -r ${DATAOUTDIR}/${HLTPFdir}
	#mkdir ${DATAOUTDIR}/${HLTPFdir}
	#
	#echo ""
	#echo "moving DATA output to ${DATAOUTDIR}/${HLTPFdir}"
	#echo ""
	#
	#mv ${DATAOUTDIR}/*.pdf   ${DATAOUTDIR}/${HLTPFdir}
	#mv ${DATAOUTDIR}/*.root  ${DATAOUTDIR}/${HLTPFdir}
	#
	#echo ""
	#echo "moving DATA output to ${FINALOUTDIR}/${HLTPFdir}"
	#echo ""
	#
        ##scp2Serin $HLTPFdir
        ##mv -f $HLTPFdir output/unfoldDataSpectra/temp    
        #rm -r ${FINALOUTDIR}/${HLTPFdir}
	#mv  ${DATAOUTDIR}/${HLTPFdir} ${FINALOUTDIR}/${HLTPFdir}
	#scp2Serin ${FINALOUTDIR}/${HLTPFdir}
    fi
    
    
    if [[ dataHLTCalo -eq 1 ]]
    then
        #source runData.sh    "0.0eta0.5_L2L3recoJetID_HLTCalo"      "0.0eta0.5_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
        #source runData.sh    "0.5eta1.0_L2L3recoJetID_HLTCalo"      "0.5eta1.0_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
        #source runData.sh    "1.0eta1.5_L2L3recoJetID_HLTCalo"      "1.0eta1.5_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
        #source runData.sh    "1.5eta2.0_L2L3recoJetID_HLTCalo"      "1.5eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins         
    	
        #source runData.sh    "0.0eta1.0_L2L3recoJetID_HLTCalo"      "0.0eta1.0_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
        #source runData.sh    "1.0eta2.0_L2L3recoJetID_HLTCalo"      "1.0eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
	
	source runData.sh    "0.0eta2.0_L2L3recoJetID_HLTCalo"      "0.0eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
	
	#echo ""
	#echo "cleaning up previous output in ${DATAOUTDIR}/${HLTCalodir}"
	#echo ""
	#
	#rm -r ${DATAOUTDIR}/${HLTCalodir}
	#mkdir ${DATAOUTDIR}/${HLTCalodir}
	#
	#echo ""
	#echo "moving output to ${DATAOUTDIR}/${HLTCalodir}"
	#echo ""
	#
	#mv ${DATAOUTDIR}/*.pdf   ${DATAOUTDIR}/${HLTCalodir}
	#mv ${DATAOUTDIR}/*.root  ${DATAOUTDIR}/${HLTCalodir}
	#
	#echo ""
	#echo "moving DATA output to ${FINALOUTDIR}/${HLTCalodir}"
	#echo ""
	#
        ##scp2Serin $HLTCalodir
        ##mv -f $HLTCalodir output/unfoldDataSpectra/temp    
        #rm -r ${FINALOUTDIR}/${HLTCalodir}
	#mv  ${DATAOUTDIR}/${HLTCalodir} ${FINALOUTDIR}/${HLTCalodir}
	
    fi
    
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
#FINALOUTDIR=/cms/heavyion/ilaflott/temp #defined above, here for convenience
MCdir=L2L3recoJetID_MCClosure_Py8LOMC_NLOMC_NLONPMC
MCOUTDIR=output/unfoldMCSpectra
################################################################################################################################
function unfoldMC(){
    
    echo ""
    echo "... cleaning up previous unfoldMCSpectra output ..."
    echo ""
    
    mv ${MCOUTDIR}/*.pdf  ${MCOUTDIR}/temp/.
    mv ${MCOUTDIR}/*.root ${MCOUTDIR}/temp/.
    
    echo ""
    echo "// ------------------------------------------------------------------------------------------------ //"
    echo "// ------------------------------------------------------------------------------------------------ //"
    echo "// ------------------------------------- RUNNING MC UNFOLDING ------------------------------------- //"
    echo "// ------------------------------------------------------------------------------------------------ //"
    echo "// ------------------------------------------------------------------------------------------------ //"
    

    if [[ $# -eq 1 ]]
    then
	#MCdir=L2L3recoJetID_MCClosure_gendr0p1_${1}
	MCdir=${MCdir}_${1}
#    else
	#MCdir=L2L3recoJetID_MCClosure_gendr0p1
	#MCdir=0.0eta2.0_toyMC_NLO_wANDwoNP
    fi

    if [[ doLOMCUnfolding -eq 1 ]]
    then
        #source runMC.sh   "0.0eta0.5_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
        #source runMC.sh   "0.5eta1.0_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
        #source runMC.sh   "1.0eta1.5_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins 
        #source runMC.sh   "1.5eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins         
        	
        #source runMC.sh   "0.0eta1.0_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
        #source runMC.sh   "1.0eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
        	
        #source runMC.sh   "0.0eta2.0_gendr0p1_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
	#source runMC.sh   "0.0eta2.0_gendr0p1_recoJetID"   "07" "02" "18" $useSimpBins  $doBayes $doSVD $LOMCSVDkReg $LOMCBayeskIter
        #source runMC.sh   "0.0eta2.0_gendr0p1_genPtLo43_recoJetID"   "07" "02" "18"   $doBayes $doSVD $useSimpBins     
	#source runMC.sh   "0.0eta2.0_unf"   "11" "14" "18" $useSimpBins  $doBayes $doSVD $LOMCSVDkReg $LOMCBayeskIter
	source runMC.sh   "0.0eta2.0_PY8unf"   "01" "06" "19" $useSimpBins  $doBayes $doSVD $LOMCSVDkReg $LOMCBayeskIter

    fi
    
    if [[ doNLOMCUnfolding -eq 1 ]]
    then
	#source runNLOMC.sh   "0.0eta2.0_gendr0p1_recoJetID" "07" "02" "18" $useSimpBins $doBayes $doSVD $NLOMCSVDkReg $NLOMCBayeskIter
	source runNLOMC.sh   "0.0eta2.0_unf" "11" "14" "18" $useSimpBins $doBayes $doSVD $NLOMCSVDkReg $NLOMCBayeskIter $useNPCorr
    fi


    
    #echo ""
    #echo "cleaning up previous MC output in ${MCOUTDIR}/${MCdir}"
    #echo ""
    #
    #rm -r ${MCOUTDIR}/${MCdir}
    #mkdir ${MCOUTDIR}/${MCdir}
    #
    #echo ""
    #echo "moving MC output to ${MCOUTDIR}/${MCdir}"
    #echo ""
    #
    #mv ${MCOUTDIR}/*.pdf   ${MCOUTDIR}/${MCdir}
    #mv ${MCOUTDIR}/*.root  ${MCOUTDIR}/${MCdir}
    #
    #
    #echo ""
    #echo "moving MC output to ${FINALOUTDIR}/${MCdir}"
    #echo ""
    #
    #rm -r ${FINALOUTDIR}/${MCdir}
    #mv  ${MCOUTDIR}/${MCdir} ${FINALOUTDIR}/${MCdir}
    #scp2Serin ${FINALOUTDIR}/${MCdir}
    
    
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










if [[ $makeRooUnf -eq 1 ]]
then
    echo ""
    echo "... making rooUnfold ..."
    echo ""
    make_rooUnfold
fi


##i'm sorry scripting gods.
echo ""
echo "compiling..."
echo ""
if [[ $doBayes -eq 1 ]]
then
    if [[ $doDataUnfoldingwLO -eq 1 ]]
    then
	rooUnfoldCompile bayesUnfoldDataSpectra.C
    fi
    if [[ $doDataUnfoldingwNLO -eq 1 ]]
    then
	rooUnfoldCompile bayesUnfoldDataSpectra_wNLO.C
    fi
    if [[ $doLOMCUnfolding -eq 1 ]]
    then
	rooUnfoldCompile bayesUnfoldMCSpectra.C
    fi
    if [[ $doNLOMCUnfolding -eq 1 ]]
    then
#	rooUnfoldCompile bayesUnfoldNLOMCSpectra.C #TODO
	rooUnfoldCompile bayesUnfoldNLOMCSpectra_v2.C
    fi
fi

if [[ $doSVD -eq 1 ]]
then
    if [[ $doDataUnfoldingwLO -eq 1 ]]
    then
	rooUnfoldCompile SVDUnfoldDataSpectra.C
    fi
    if [[ $doDataUnfoldingwNLO -eq 1 ]]
    then
	rooUnfoldCompile SVDUnfoldDataSpectra_wNLO.C
    fi
    if [[ $doLOMCUnfolding -eq 1 ]]
    then
	rooUnfoldCompile SVDUnfoldMCSpectra.C
    fi
    if [[ $doNLOMCUnfolding -eq 1 ]]
    then
	rooUnfoldCompile SVDUnfoldNLOMCSpectra.C
#	rooUnfoldCompile SVDUnfoldNLOMCSpectra_v2.C #TODO
    fi    
fi


if [[ $# -eq 1 ]]
then
    dirName=$1
    if ([ $doDataUnfoldingwLO == 1 ] || [ $doDataUnfoldingwNLO == 1 ])
    then
	unfoldData $1    
    fi
    if ([ $doLOMCUnfolding -eq 1 ] || [ $doNLOMCUnfolding == 1 ])
    then
	unfoldMC $1
    fi
else
    if ([ $doDataUnfoldingwLO == 1 ] || [ $doDataUnfoldingwNLO == 1 ])
    then
	unfoldData    
    fi
    if ([ $doLOMCUnfolding -eq 1 ] || [ $doNLOMCUnfolding == 1 ])
    then
	unfoldMC    
    fi
fi

echo ""
echo "cleaning up previous unfolding output in ${FINALOUTDIR}/current_output"
echo ""
rm -r  ${FINALOUTDIR}/old_output/*
mv -f ${FINALOUTDIR}/current_output/*  ${FINALOUTDIR}/old_output/.


if ([ $doDataUnfoldingwLO == 1 ] || [ $doDataUnfoldingwNLO == 1 ])
then
    
    echo ""
    echo "cleaning up previous DATA output in ${DATAOUTDIR}/${HLTPFdir}"
    echo ""
    
    rm -r ${DATAOUTDIR}/${HLTPFdir}
    mkdir ${DATAOUTDIR}/${HLTPFdir}
    
    echo ""
    echo "moving DATA output to ${DATAOUTDIR}/${HLTPFdir}"
    echo ""
    
    mv ${DATAOUTDIR}/*.pdf   ${DATAOUTDIR}/${HLTPFdir}
    mv ${DATAOUTDIR}/*.root  ${DATAOUTDIR}/${HLTPFdir}
    
    echo ""
    echo "moving DATA output to ${FINALOUTDIR}/current_output/${HLTPFdir}"
    echo ""
    
    #rm -r ${FINALOUTDIR}/current_output/${HLTPFdir}
    mv ${DATAOUTDIR}/${HLTPFdir} ${FINALOUTDIR}/current_output/${HLTPFdir}

fi


if ([ $doLOMCUnfolding -eq 1 ] || [ $doNLOMCUnfolding == 1 ])
then
    
    echo ""
    echo "cleaning up previous MC output in ${MCOUTDIR}/${MCdir}"
    echo ""
    
    rm -r ${MCOUTDIR}/${MCdir}
    mkdir ${MCOUTDIR}/${MCdir}
    
    echo ""
    echo "moving MC output to ${MCOUTDIR}/${MCdir}"
    echo ""
    
    mv ${MCOUTDIR}/*.pdf   ${MCOUTDIR}/${MCdir}
    mv ${MCOUTDIR}/*.root  ${MCOUTDIR}/${MCdir}
    
    
    echo ""
    echo "moving MC output to ${FINALOUTDIR}/current_output/${MCdir}"
    echo ""
    
    #rm -r ${FINALOUTDIR}/current_output/${MCdir}
    mv  ${MCOUTDIR}/${MCdir} ${FINALOUTDIR}/current_output/${MCdir}

fi

scp2Serin ${FINALOUTDIR}/current_output




return













































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

