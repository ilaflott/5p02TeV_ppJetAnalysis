#!/bin/bash

NARGS=7
if [[ "$#" -eq $NARGS ]]
then
    etaBin=$1
    NLOfilestr=$2
    NLOfileshortstr=$3
    JERfilestr=$4
    Descstr=$5
    Fittypestr=$6
    scp2Serin=$7
else
    echo ""
    echo "error! usage is:"
    echo "source run_smearTheory_etabin.sh <etaBin>  <NLOfilestr> <NLOfileshortstr> <JERfilestr> <Descstr> <Fittypestr> <scp2Serin>"
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
    
    OUTDIR=${NLOfileshortstr}_${descstr}_noJERscales_${Fittypestr}_gaussSmear_NLO_plots/
    OUTFILE=${NLOfileshortstr}_${descstr}_noJERscales_${Fittypestr}_gaussSmear_${etabinstr}.root     
    JERFILE=${JERfilestr}${etabinstr}_${descstr}.root
    JERPDFFILE=${JERfilestr}${etabinstr}_${descstr}.pdf
    
    #JERFILE=${JERfilestr}${etabinstr}_semifinal.root
    #JERPDFFILE=${JERfilestr}${etabinstr}_semifinal.pdf

    echo ""
    echo "OUTDIR =${OUTDIR}"
    echo "OUTFILE=${OUTFILE}"
    echo "JERFILE=${JERFILE}"
    echo "JERPDFFILE=${JERPDFFILE}"
    echo ""    
    
    #return     #DEBUG
    
    echo ""
    echo "RUN NLO smearing"
    echo ""
    
    ## USAGE
    ## ./smearTheorySpectra_gaussCoreJER_etabin.exe <NLOfilestr> <useSplineWeights> <fitType_str> <JERFILE> <etabin> <outputfile>
    ./smearTheorySpectra_gaussCoreJER_etabin.exe "${NLOfilestr}"  "${Fittypestr}"  "${JERFILE}"  "${etaBin}"  "${OUTFILE}"    
    
    echo "" 
    echo "DONE NLO smearing."
    echo ""

    if [[ -d "${OUTDIR}" ]]
    then	
	echo "$OUTDIR exists!"
	#echo "removing directory!"
	#rm -rf $OUTDIR
	#echo "recreating directory!"
	#mkdir $OUTDIR
    else 
	mkdir  $OUTDIR
    fi
    
    echo ""
    echo "moving $OUTFILE to $OUTDIR"
    echo ""
    
    cp smearTheorySpectra_gaussCoreJER_etabin.C $OUTDIR/.
    cp smearTheorySpectra.h $OUTDIR/.
    cp $OUTFILE $OUTDIR/.
    cp $JERFILE $OUTDIR/.
    cp $JERPDFFILE $OUTDIR/.

    #return    #DEBUG
    
    echo ""
    echo "writing plots in $OUTFILE to $OUTDIR"
    echo ""
    ROOTEXC="quickCheck_gausSmThy_etabin.C++( \"${OUTDIR}\", \"${OUTFILE}\", ${etaBin})"
    ##TEMP COMMENT OUT
    root -l -b -q "${ROOTEXC}"
    #return
    ##TEMP COMMENT OUT

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