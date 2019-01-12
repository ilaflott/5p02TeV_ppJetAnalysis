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
    cp2Serin=$7
else
    echo ""
    echo "error! usage is:"
    echo "source run_smearTheory_etabin.sh <etaBin>  <NLOfilestr> <NLOfileshortstr> <JERfilestr> <Descstr> <Fittypestr> <cp2Serin>"
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
    
    OUTDIR=${NLOfileshortstr}_${descstr}_${Fittypestr}_gaussSmear_plots/
    OUTFILE=${NLOfileshortstr}_${descstr}_${Fittypestr}_gaussSmear_${etabinstr}.root     
    JERFILE=${JERfilestr}${etabinstr}.root
    
    
    echo ""
    echo "OUTDIR =${OUTDIR}"
    echo "OUTFILE=${OUTFILE}"
    echo "JERFILE=${JERFILE}"
    echo ""    
    

    echo ""
    echo "RUN NLO smearing"
    echo ""
    
    ## USAGE
    ## ./smearTheorySpectra_gaussCoreJER_etabin.exe <NLOfilestr> <useSplineWeights> <fitType_str> <JERFILE> <etabin> <outputfile>
    ./smearTheorySpectra_gaussCoreJER_etabin.exe "${NLOfilestr}"  "${Fittypestr}"  "${JERFILE}"  "${etaBin}"  "${OUTFILE}"    
    
    #return

    echo "" 
    echo "DONE NLO smearing."
    echo ""
    
    if [[ -d "${OUTDIR}" ]]
    then
	echo "$OUTDIR exists!"
    else 
	mkdir  $OUTDIR
    fi
    
    echo ""
    echo "moving $OUTFILE to $OUTDIR"
    echo ""
    
    cp $OUTFILE $OUTDIR/.
    
    echo ""
    echo "writing plots in $OUTFILE to $OUTDIR"
    echo ""
    ROOTEXC="quickCheck_gausSmThy_etabin.C++( \"${OUTDIR}\", \"${OUTFILE}\", ${etaBin})"
    root -l -b -q "${ROOTEXC}"
    
    echo ""
    echo "scp2Serin"
    echo ""

    if [[ $cp2Serin -eq 1 ]]
    then
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