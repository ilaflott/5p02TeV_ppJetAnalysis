#!/bin/bash

#make_rooUnfold

rooUnfoldCompile bayesUnfoldMCSpectra_wPY8_etabin.C

############# NNPDF NLO + HERWIG EE5C LO NPCs ###############


####CT14NNLO_03.18.20_sigmu_geny_noJERscales_spl3wgts_gaussSmear_TEST_plots

FITTYPE="sigmu"
#DATADIR="03.18.20_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_03-18-20_unf_0.0eta2.5_SMPbins_withgenjety_semifinalv3"
DATADIR="01.05.21_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_01-05-21_unf_0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4/"
THYDIR="smearTheory/PY8_03.18.20_sigmu_geny_noJERscales_spl3wgts_gaussSmear_semifinalv4_plots/"
THYFILESTR="PY8_03.18.20_sigmu_geny_noJERscales_spl3wgts_gaussSmear_semifinalv4_"
OUTPUTTAG="01.05.21_Bayes_RECOPY8_smearPY8_v4_${FITTYPE}_mergedSMPbins_semifinalv4"
#OUTPUTTAG="03.18.20_Bayes_CT14NNLO_NLO_v3_${FITTYPE}_johnbins2_semifinalv3_murmufHTp_JohnNPs"

applyNPcorrs=0
dosystunf=0
declare -a sysarr=("JEC")
declare -a sysarropt=("")
etabinStart=0
etabinEnd=4
sysarrlen=${#sysarr[@]}
for (( j=${etabinStart}; j<${etabinEnd}; j++ ));
do
    for (( i=0; i<${sysarrlen}; i++ ));
    do
	echo ""
	echo "j=$j"
	echo "${OUTPUTTAG}_${sysarr[i]}${sysarropt[i]}sys"
	echo "${sysarr[i]}"
	echo "${sysarropt[i]}"
	./bayesUnfoldMCSpectra_wPY8_etabin.exe ${DATADIR} ${OUTPUTTAG}_${sysarr[i]}sys${sysarropt[i]} $j ${THYDIR} ${THYFILESTR} ${dosystunf} ${sysarr[i]} ${sysarropt[i]} ${applyNPcorrs}
	#./bayesUnfoldMCSpectra_wPY8_etabin.exe ${DATADIR} ${OUTPUTTAG}_${sysarr[i]}sys${sysarropt[i]} $j ${THYDIR} ${THYFILESTR} ${dosystunf} ${sysarr[i]} ${sysarropt[i]} 
    done
done
return



#./bayesUnfoldMCSpectra_wPY8_etabin.exe      argv[0]
#${DATADIR} 				     argv[1]
#${OUTPUTTAG}_${sysarr[i]}sys${sysarropt[i]} argv[2]
#$j 					     argv[3]
#${THYDIR} 				     argv[4]
#${THYFILESTR} 				     argv[5]
#${dosystunf} 				     argv[6]
#${sysarr[i]} 				     argv[7]
#${sysarropt[i]} 			     argv[8]
#${applyNPcorrs}                             argv[9]
