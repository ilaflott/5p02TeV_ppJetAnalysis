#!/bin/bash


smeartheory=1
smearPY8=0



if [[ $smeartheory -eq 1 ]]
then
    
    ##### --------------------------------------------------------------------------- #####
    # NEW SMEARTHEORY SCRIPT; 0.5 BIN WIDTHS UP TO 3.0 (pending available output of course)
    ##### --------------------------------------------------------------------------- #####
    
    echo ""
    echo "compiling smearTheorySpectra_gaussCoreJER_etabin"
    echo ""
    
    rooUnfoldCompile smearTheorySpectra_gaussCoreJER_etabin.C
    
    
### INPUTS FOR run_smearTheory_etabin
    #nlofstr="fnl5020_LO2_R04Jets_modify_CT10nlo_HS"
    #nlofstr="fnl5020_LO2_R04Jets_modify_CT14nlo_HS"
    #nlofstr="fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS"
    #nlofstr="fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS"
    nlofstr="fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC"
    
    #nlofshortstr="CT10"
    #nlofshortstr="CT14"
    #nlofshortstr="HERAPDF"
    #nlofshortstr="MMHT"
    nlofshortstr="NNPDF"
    
    jerfstr="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots_JERS/output/ak4PF_PY8JER_"
    #descstr="02.18.19_semifinal_v2_w8TeVJERscales" #_ptmax1000"
    #descstr="ptLo28_drCut_geny_sigma"
    #descstr="ptLo28_drCut_geny_sigmu"
    #descstr="03.14.19_sigmu"
    descstr="03.14.19_sigma"
    
    #fittypestr="modLog"
    #fittypestr="7TeV"
    fittypestr="spl3wgts"
    #fittypestr="spl3wgts_extv2"
    ##fittypestr="spl3wgts_extv1"
    ##fittypestr="spl3wgts_extv3"
    
    scpoutput=0
    scpafterlast=1
    etabin_i=0
    #Netabins=1 ##debug; runs first etabin only
    Netabins=4 
    while [ $etabin_i -lt $Netabins ]
    do
	if [[ $scpafterlast -eq 1 ]]
	then
	    if [[ $etabin_i -eq $(($Netabins - 1)) ]]
	    then
    		scpoutput=1
	    fi
	fi
	
	echo "etabin_i=${etabin_i}"
	source run_smearTheory_etabin.sh $etabin_i   $nlofstr $nlofshortstr $jerfstr $descstr $fittypestr $scpoutput
	etabin_i=$(($etabin_i + 1))
	
    done

fi

if [[ $smearPY8 -eq 1 ]]
then
    
    ##### --------------------------------------------------------------------------- #####
    # NEW SMEARTHEORY SCRIPT; 0.5 BIN WIDTHS UP TO 3.0 (pending available output of course)
    ##### --------------------------------------------------------------------------- #####
    
    echo ""
    echo "compiling smearPY8Spectra_gaussCoreJER_etabin"
    echo ""
    
    rooUnfoldCompile smearPY8Spectra_gaussCoreJER_etabin.C
    
    
    ### INPUTS FOR run_smearTheory_etabin
    #nlofstr="fnl5020_LO2_R04Jets_modify_CT10nlo_HS"
    #nlofstr="fnl5020_LO2_R04Jets_modify_CT14nlo_HS"
    #nlofstr="fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS"
    #nlofstr="fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS"
    #nlofstr="fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC"
    py8fstr="02.18.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-18-19_JERS_0.0eta3.0_unf/Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root"
    
    #nlofshortstr="CT10"
    #nlofshortstr="CT14"
    #nlofshortstr="HERAPDF"
    #nlofshortstr="MMHT"
    #nlofshortstr="NNPDF"
    py8fshortstr="PY8"
    
    jerfstr="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots_JERS/output/ak4PF_PY8JER_"
    #descstr="02.18.19_semifinal_v2_w8TeVJERscales" #_ptmax1000"
    #descstr="ptLo28_drCut_geny_sigma"
    #descstr="ptLo28_drCut_geny_sigmu"
    descstr="ptLo28_drCut_geny_sigmu"
    
    #fittypestr="modLog"
    #fittypestr="7TeV"
    fittypestr="spl3wgts"
    #fittypestr="spl3wgts_extv2"
    ##fittypestr="spl3wgts_extv1"
    ##fittypestr="spl3wgts_extv3"
    
    scpoutput=0
    scpafterlast=1
    #etabin_i=0
    #Netabins=6 ##full set of etabins 
    etabin_i=0 
    Netabins=4 ##debug; runs first etabin only
    while [ $etabin_i -lt $Netabins ]
    do
	if [[ $scpafterlast -eq 1 ]]
	then
	    if [[ $etabin_i -eq $(($Netabins - 1)) ]]
	    then
    		scpoutput=1
	    fi
	fi
	
	echo "etabin_i=${etabin_i}"
	source run_smearPY8_etabin.sh $etabin_i   $py8fstr $py8fshortstr $jerfstr $descstr $fittypestr $scpoutput
	etabin_i=$(($etabin_i + 1))
	
    done
    
fi

return
    
    


















#### OLD    
##### ------------------------------ #####
# DEFAULT SMEARTHEORY SCRIPT; 00eta20 ONLY
##### ------------------------------ #####
#echo ""
#echo "compiling smearTheorySpectra_gaussCoreJER"
#echo ""
#rooUnfoldCompile smearTheorySpectra_gaussCoreJER.C
#source run_smearTheory.sh   ##00eta20

