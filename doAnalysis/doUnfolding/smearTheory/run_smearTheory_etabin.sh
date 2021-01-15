#!/bin/bash

NARGS=8
if [[ "$#" -eq $NARGS ]]
then
    etaBin=$1
    NLOfilestr=$2
    NLOfileshortstr=$3
    JERfilestr=$4
    jerdescstr=$5
    Fittypestr=$6
    scp2Serin=$7
    outdescst=$8
else
    echo ""
    echo "error! usage is:"
    echo "source run_smearTheory_etabin.sh <etaBin>  <NLOfilestr> <NLOfileshortstr> <JERfilestr> <jerdescstr> <Fittypestr> <scp2Serin> <outdescstr> "
    echo ""
    return
fi



#### RIG THIS FOR smearTheory_gaussCoreJER_etabin.exe #####
if [[ "$#" -eq $NARGS ]]
then
    

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
    elif [[ $etaBin -eq 4 ]]
    then
	etabinstr="20eta25"
    elif [[ $etaBin -eq 5 ]]
    then
	etabinstr="25eta30"
    fi
    
    
    echo ""
    echo "etaBin=$etaBin"
    echo "etabinstr=${etabinstr}"
    echo ""
    
    #OUTDIR=${NLOfileshortstr}_NLO_v2_${descstr}_noJERscales_${Fittypestr}_gaussSmear_TEST_plots/
    #OUTFILE=${NLOfileshortstr}_NLO_v2_${descstr}_noJERscales_${Fittypestr}_gaussSmear_TEST_${etabinstr}.root    

    ##OUTDIR=${NLOfileshortstr}_NLO_v3_${descstr}_noJERscales_${Fittypestr}_gaussSmear_TEST_murmufpt_plots/
    ##OUTFILE=${NLOfileshortstr}_NLO_v3_${descstr}_noJERscales_${Fittypestr}_gaussSmear_TEST_murmufpt_${etabinstr}.root    
    ##OUTDIR=${NLOfileshortstr}_NLO_v3_${descstr}_noJERscales_${Fittypestr}_gaussSmear_TEST_murmufpt1_plots/
    ##OUTFILE=${NLOfileshortstr}_NLO_v3_${descstr}_noJERscales_${Fittypestr}_gaussSmear_TEST_murmufpt1_${etabinstr}.root    
    ##OUTDIR=${NLOfileshortstr}_NLO_v3_${descstr}_${Fittypestr}_gaussSmear_murmufpt_DEBUG_JERSF1.0_wREDRAW_plots/
    ##OUTFILE=${NLOfileshortstr}_NLO_v3_${descstr}_${Fittypestr}_gaussSmear_murmufpt_DEBUG_JERSF1.0_wREDRAW_${etabinstr}.root    
    #OUTDIR=${NLOfileshortstr}_NLO_v3_${descstr}_${Fittypestr}_gaussSmear_murmufpt1_plots/
    #OUTFILE=${NLOfileshortstr}_NLO_v3_${descstr}_${Fittypestr}_gaussSmear_murmufpt1_${etabinstr}.root    
    #OUTDIR=${NLOfileshortstr}_NLO_v3_${descstr}_${Fittypestr}_gaussSmear_murmufpt_plots/
    #OUTFILE=${NLOfileshortstr}_NLO_v3_${descstr}_${Fittypestr}_gaussSmear_murmufpt_${etabinstr}.root    
    OUTDIR=${NLOfileshortstr}_NLO_${jerdescstr}_${Fittypestr}_${outdescstr}_plots/
    OUTFILE=${NLOfileshortstr}_NLO_${jerdescstr}_${Fittypestr}_${outdescstr}_${etabinstr}.root    
    JERDIR=/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots_JERS/output/
    JERFILE=${JERfilestr}${etabinstr}_${jerdescstr}.root
    FULLJERFILEPATH=${JERDIR}${JERFILE}
    JERPDFFILE=${JERfilestr}${etabinstr}_${jerdescstr}.pdf    
    FULLJERPDFFILEPATH=${JERDIR}${JERPDFFILE}

    #return     #DEBUG
    
    echo ""
    echo "RUN NLO smearing"
    echo "OUTFILE=${OUTFILE}"
    echo "FULLJERFILEPATH=${FULLJERFILEPATH}"
    echo ""

    ## USAGE
    # ./smearTheorySpectra_gaussCoreJER_etabin.exe <NLOfilestr> <useSplineWeights> <fitType_str> <JERFILE> <etabin> <outputfile>
    ./smearTheorySpectra_gaussCoreJER_etabin.exe "${NLOfilestr}"  "${Fittypestr}"  "${FULLJERFILEPATH}"  "${etaBin}"  "${OUTFILE}"

    echo "" 
    echo "DONE NLO smearing."
    echo ""

    if [[ -d "${OUTDIR}" ]]
    then	
	echo "WARNING $OUTDIR exists!"
	#echo "recreating directory!"
	#mkdir $OUTDIR
    else 
	mkdir  $OUTDIR
    fi
    




    echo ""
    echo "copying code to $OUTDIR"
    echo ""    

    cp smearTheorySpectra_gaussCoreJER_etabin.C $OUTDIR/.
    cp smearTheorySpectra.h $OUTDIR/.
    
    echo ""
    echo "moving $OUTFILE to $OUTDIR"
    echo ""

    mv $OUTFILE $OUTDIR/.
    
    
    echo ""
    echo "JERDIR=${JERDIR}"
    echo "JERFILE=${JERFILE}"
    echo "FULLJERFILEPATH=${FULLJERFILEPATH}"
    cp $FULLJERFILEPATH $OUTDIR/.
    echo ""
    echo "JERPDFFILE=${JERPDFFILE}"
    echo "FULLJERPDFFILEPATH=${FULLJERPDFFILEPATH}"
    cp $FULLJERPDFFILEPATH $OUTDIR/.
    echo ""

    #return    #DEBUG
    
    echo ""
    echo "writing plots in $OUTFILE to $OUTDIR"
    echo ""
    NLOROOTEXC="quickCheck_gausSmThy_etabin.C++( \"${OUTDIR}\", \"${OUTFILE}\", ${etaBin})"
    root -l -b -q "${NLOROOTEXC}"

    #echo ""
    #echo "writing plots in $JERFILE to $OUTDIR"
    #echo ""
    #JERROOTEXC="quickCheck_gausSmThy_etabin.C++( \"${OUTDIR}\", \"${JERFILE}\", ${etaBin})"
    #root -l -b -q "${JERROOTEXC}"

    echo ""
    echo "scp2Serin"
    echo ""

    if [[ $scp2Serin -eq 1 ]]
    then
	echo "OUTDIR=$OUTDIR"
	scp2Serin ${OUTDIR}    
    fi
    
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
