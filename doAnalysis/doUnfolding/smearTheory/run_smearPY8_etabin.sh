#!/bin/bash

NARGS=8
if [[ "$#" -eq $NARGS ]]
then
    etaBin=$1
    PY8filestr=$2
    PY8fileshortstr=$3
    JERfilestr=$4
    jerdescstr=$5
    Fittypestr=$6
    scp2Serin=$7
    outdescstr=$8
else
    echo ""
    echo "error! usage is:"
    echo "source run_smearPY8_etabin.sh <etaBin>  <PY8filestr> <PY8fileshortstr> <JERfilestr> <jerdescstr> <Fittypestr> <csp2Serin> <outdescstr>"
    echo ""
    return
fi



#### RIG THIS FOR smearPY8_gaussCoreJER_etabin.exe #####
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
    
    #OUTDIR=${NLOfileshortstr}_NLO_v3_${jerdescstr}_${Fittypestr}_${outdescstr}_plots/
    #OUTFILE=${NLOfileshortstr}_NLO_v3_${jerdescstr}_${Fittypestr}_${outdescstr}_${etabinstr}.root    
    OUTDIR=${PY8fileshortstr}_${jerdescstr}_noJERscales_${Fittypestr}_${outdescstr}_plots/
    OUTFILE=${PY8fileshortstr}_${jerdescstr}_noJERscales_${Fittypestr}_${outdescstr}_${etabinstr}.root     
    JERDIR=/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots_JERS/output/
    JERFILE=${JERfilestr}${etabinstr}_${jerdescstr}.root
    FULLJERFILEPATH=${JERDIR}${JERFILE}
    JERPDFFILE=${JERfilestr}${etabinstr}_${jerdescstr}.pdf    
    FULLJERPDFFILEPATH=${JERDIR}${JERPDFFILE}
    
    echo ""
    echo "RUN PY8 smearing"
    echo "OUTDIR =${OUTDIR}"
    echo "OUTFILE=${OUTFILE}"
    #echo "JERFILE=${JERFILE}"
    #echo "JERPDFFILE=${JERPDFFILE}"
    echo ""        
    
    ## USAGE
    ## ./smearPY8Spectra_gaussCoreJER_etabin.exe <PY8filestr> <useSplineWeights> <fitType_str> <JERFILE> <etabin> <outputfile>
    ./smearPY8Spectra_gaussCoreJER_etabin.exe "${PY8filestr}"  "${Fittypestr}"  "${FULLJERFILEPATH}"  "${etaBin}"  "${OUTFILE}"    
    
    echo "" 
    echo "DONE PY8 smearing."
    echo ""

    #return

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
    echo "copying code to $OUTDIR"
    echo ""    
    
    cp smearPY8Spectra_gaussCoreJER_etabin.C $OUTDIR/.
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
	echo "OUTDIR=$OUTDIR"
	scp2Serin ${OUTDIR}    
    fi
    
else 
    echo "bad input. see macro's insides."
fi

echo ""
echo "done w/ run_smearPY8.sh"
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
