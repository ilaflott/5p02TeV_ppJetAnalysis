#!/bin/bash





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

#nlofshortstr="CT10_NLO"
#nlofshortstr="CT14_NLO"
#nlofshortstr="HERAPDF_NLO"
#nlofshortstr="MMHT_NLO"
nlofshortstr="NNPDF_NNLO"

jerfstr="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots_JERS/output/ak4PF_MCJEC_"
descstr="1.10.19_ptmax1000"

#fittypestr="modLog"
#fittypestr="7TeV"
#fittypestr="spl3wgts"
fittypestr="spl3wgts_extv2"
##fittypestr="spl3wgts_extv1"
##fittypestr="spl3wgts_extv3"

scpoutput=0
etabin_i=0
Netabins=4
#Netabins=1 ##debug; runs one bin only

while [ $etabin_i -lt $Netabins ]
do
    if [[ $etabin_i -eq $(($Netabins - 1)) ]]
    then
	scpoutput=1
    fi
    
    echo "etabin_i=${etabin_i}"
    source run_smearTheory_etabin.sh $etabin_i   $nlofstr $nlofshortstr $jerfstr $descstr $fittypestr $scpoutput
    etabin_i=$(($etabin_i + 1))
done





##### ------------------------------ #####
# DEFAULT SMEARTHEORY SCRIPT; 00eta20 ONLY
##### ------------------------------ #####
#echo ""
#echo "compiling smearTheorySpectra_gaussCoreJER"
#echo ""
#rooUnfoldCompile smearTheorySpectra_gaussCoreJER.C
#source run_smearTheory.sh   ##00eta20

