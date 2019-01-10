#!/bin/bash

#    NLOfilestr="fnl5020_LO2_R04Jets_modify_CT10nlo_HS"
#    NLOfilestr="fnl5020_LO2_R04Jets_modify_CT14nlo_HS"
#    NLOfilestr="fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS"
#    NLOfilestr="fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS"
NLOfilestr="fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC"
NLOfileshortstr="NNPDF_NNLO"
JERfilestr="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots_JERS/output/ak4PF_MCJEC_"    #"00eta05.root"
#JERfilestr="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots_JERS/output/ak4PF_MCJEC_00eta20.root"
Descstr="spl3wgts_ptmax1000_extv2.2"
Fittypestr="modLog"

#### RIG THIS FOR smearTheory_gaussCoreJER_etabin.exe #####
if [[ "$#" -eq 1 ]]
then
    
    etaBin=$1
    echo ""
    echo "etaBin=$etaBin"
    echo ""

    if [[ $etaBin -eq 0 ]]
    then
	etabinstr="00eta05"
    elif [[ $etaBin -eq 1 ]]
    then
	etabinstr="05eta10"
    elif [[ $etaBin -eq 2 ]]
    then
	etabinstr="10eta15"
    elif [[ $etaBin -eq 3 ]]
    then
	etabinstr="15eta20"
    fi


    echo ""
    echo "etabinstr=${etabinstr}"
    echo ""
    
    OUTDIR=${NLOfileshortstr}_${etabinstr}_${descstr}_gausSmThy_plots/
    OUTFILE=${NLOfilestr}_${descstr}_gaussSmear_${etabinstr}.root     
    JERfilestr=${JERfilestr}${etabinstr}.root

    
    echo ""
    echo "OUTDIR     =${OUTDIR}"
    echo "OUTFILE=${OUTFILE}"
    echo "JERfilestr=${JERfilestr}"
    echo ""
    
    echo ""
    echo "RUN NLO smearing"
    echo ""
    
    ## USAGE
    ## ./smearTheorySpectra_gaussCoreJER_etabin.exe <NLOfilestr> <useSplineWeights> <fitType_str> <JERfilestr> <etabin> <outputtag>
    ./smearTheorySpectra_gaussCoreJER_etabin.exe "${NLOfilestr}" "1"  "${Fittypestr}"  "${JERfilestr}"  "${etaBin}"  "${etabinstr}"    
    
    echo "" 
    echo "DONE NLO smearing."
    echo ""
    
    rm -rf $OUTDIR
    mkdir  $OUTDIR
    cp $OUTFILE $OUTDIR/.
    
    #root -l -b -q quickCheck_gausSmThy.C++ $OUTDIR $OUTFILE
    
    echo ""
    echo "scp2Serin"
    echo ""
    
    scp2Serin ${OUTDIR}    
    
    
elif [[ "$#" -eq 0 ]]
then
    
    SPLINEDIR=NNPDF_NNLO_00eta20_spl3wgts_ptmax1000_extv2.2_gausSmThy_plots/
    SPLINEFILENAME=fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_spl3wgts_ptmax1000_extv2.2_gaussSmear_00eta20.root 
    ##DEBUG
    #rootcompile smearTheorySpectra_gaussCoreJER.C    
    
    echo ""
    echo "RUN NLO smearing"
    echo ""

    #./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_CT10nlo_HS
    #./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_CT14nlo_HS
    #./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS 
    #./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS
    
    #./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC 0    #fitweights
    ./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC 1    #splineweights
    
    #rm NNPDF_NNLO_00eta20_fitwgts_gausSmThy_plots/*
    rm NNPDF_NNLO_00eta20_spl3wgts_gausSmThy_plots/*
    
    root -l -b -q quickCheck_gausSmThy.C++
    
    rm -rf $SPLINEDIR
    mkdir $SPLINEDIR
    mv NNPDF_NNLO_00eta20_spl3wgts_gausSmThy_plots/*png $SPLINEDIR/.
    
    echo ""
    echo "DONE NLO smearing."
    echo "DONE NLO smearing. scp2Serin"
    echo ""
    
    ##scp2Serin fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_00eta20.root
    ##scp2Serin fnl5020_LO2_R04Jets_modify_CT14nlo_HS_gaussSmear_00eta20.root
    ##scp2Serin fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS_gaussSmear_00eta20.root
    ##scp2Serin fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS_gaussSmear_00eta20.root
    
    #cp fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_gaussSmear_00eta20.root NNPDF_NNLO_00eta20_fitwgts_gausSmThy_plots/.
    #cp fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_gaussSmear_00eta20.root fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_ptmax686_gaussSmear_00eta20.root


    cp fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_spl3wgts_gaussSmear_00eta20.root $SPLINEFILENAME    
    cp $SPLINEFILENAME $SPLINEDIR
    
    #scp2Serin "NNPDF_NNLO_00eta20_*_gausSmThy_plots"
    #scp2Serin NNPDF_NNLO_00eta20_fitwgts_gausSmThy_plots    
    #scp2Serin NNPDF_NNLO_00eta20_spl3wgts_gausSmThy_plots        
    scp2Serin $SPLINEDIR



else 
    echo "bad input. see macro's insides."
fi

echo ""
echo "done w/ run_smearTheory.sh"
echo ""

return











##################OLD




#
#
#elif [[ "$#" -eq 1 ]]
#then
#    etaBin=$1
#    echo ""
#    echo "compiling JER_gaussCore_smearTheoryCurves_singleYbin"
#    echo ""
#    
#    rooUnfoldCompile smearTheorySpectra_gaussCoreJER_singleYbin.C
#    
#    echo ""
#    echo "RUN NLO smearing"
#    echo ""
#    
#    ./smearTheorySpectra_gaussCoreJER_singleYbin.exe fnl5020_LO2_R04Jets_modify_CT10nlo_HS $etaBin
##./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_CT14nlo_HS $etaBin
##./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS $etaBin
##./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS $etaBin
##./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC $etaBin
#    
#    echo ""
#    echo "DONE NLO smearing."
#    echo "DONE NLO smearing. scp2Serin"
#    echo ""
#    
#    outfile="fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_etaBin${etaBin}.root"
#    scp2Serin $outfile
#elif [[ "$#" -eq 2 ]]
#then
#    
#    echo ""
#    echo "compiling JER_gaussCore_smearTheoryCurves_singleYbin"
#    echo ""
#    
#    rooUnfoldCompile smearTheorySpectra_gaussCoreJER_singleYbin.C
#    
#    echo ""
#    echo "RUN NLO smearing"
#    echo ""
#    
#    ./smearTheorySpectra_gaussCoreJER_singleYbin.exe fnl5020_LO2_R04Jets_modify_CT10nlo_HS 1
#    ./smearTheorySpectra_gaussCoreJER_singleYbin.exe fnl5020_LO2_R04Jets_modify_CT10nlo_HS 2
#    ./smearTheorySpectra_gaussCoreJER_singleYbin.exe fnl5020_LO2_R04Jets_modify_CT10nlo_HS 3
#    ./smearTheorySpectra_gaussCoreJER_singleYbin.exe fnl5020_LO2_R04Jets_modify_CT10nlo_HS 4
##./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_CT14nlo_HS $etaBin
##./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS $etaBin
##./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS $etaBin
##./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC $etaBin
#    
#    echo ""
#    echo "DONE NLO smearing."
#    echo "DONE NLO smearing. scp2Serin"
#    echo ""
#    
#    scp2Serin fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_etaBin\?.root