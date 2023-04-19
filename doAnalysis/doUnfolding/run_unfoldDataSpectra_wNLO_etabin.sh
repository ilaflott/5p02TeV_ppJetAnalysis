#!/bin/bash

FITTYPE="sigma"
#FITTYPE="sigmu"

#make_rooUnfold

rooUnfoldCompile bayesUnfoldDataSpectra_wNLO_etabin.C


DATADIR="01.05.21_outputCondor/ppData_HighPtJetTrig_ak4PFJets_01-05-21_jetPlots_0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"

#####NEW CRYSTAL BALL FIT SMEARING
THYDIR="smearTheory/CT14NNLO_NLO_03.18.20_sigmu_geny_crysball_ptLo56_spl3wgts_murmufHTp_JohnNPs_crysballsmear_plots/"
THYFILESTR="CT14NNLO_NLO_03.18.20_sigmu_geny_crysball_ptLo56_spl3wgts_murmufHTp_JohnNPs_crysballsmear_"
#OUTPUTTAG="09.27.21_Bayes_CT14NNLO_NLO_sigmu_murmufHTp_SMPbins_wJERSFs_semifinalv4_crysballsmear_wPY8MatchedFakesv4"
OUTPUTTAG="09.27.21_Bayes_CT14NNLO_NLO_sigmu_murmufHTp_SMPbins_wJERSFs_semifinalv4_crysballsmear_NOPY8MatchedFakes"

##__________________________________________________________________________________##
###my gauss core-extracted params+ fits of them
#THYDIR="smearTheory/CT14NNLO_NLO_09.27.21_${FITTYPE}_wPFMET_RECOptLo10_GENptLo3_genbinsonly_TESTFIT_JOHNDSCB_GENYbinEQYbin_semifinalv4_spl3wgts_murmufHTp_JohnNPs_gaussmear_plots/"
#THYFILESTR="CT14NNLO_NLO_09.27.21_${FITTYPE}_wPFMET_RECOptLo10_GENptLo3_genbinsonly_TESTFIT_JOHNDSCB_GENYbinEQYbin_semifinalv4_spl3wgts_murmufHTp_JohnNPs_gaussmear_"
#OUTPUTTAG="09.27.21_Bayes_CT14NNLO_NLO_${FITTYPE}_murmufHTp_SMPbins_wJERSFs_semifinalv4_gaussmear_wPY8MatchedFakesv4"

##__________________________________________________________________________________##
#johns DSCB-extracted params+ fits of them
#THYDIR="smearTheory/CT14NNLO_NLO_09.27.21_${FITTYPE}_wPFMET_RECOptLo10_GENptLo3_genbinsonly_TESTFIT_JOHNDSCB_GENYbinEQYbin_semifinalv4_spl3wgts_murmufHTp_JohnNPs_gauscrysballsmear_plots/"
#THYFILESTR="CT14NNLO_NLO_09.27.21_${FITTYPE}_wPFMET_RECOptLo10_GENptLo3_genbinsonly_TESTFIT_JOHNDSCB_GENYbinEQYbin_semifinalv4_spl3wgts_murmufHTp_JohnNPs_gauscrysballsmear_"
#OUTPUTTAG="09.27.21_Bayes_CT14NNLO_NLO_${FITTYPE}_murmufHTp_SMPbins_wJERSFs_semifinalv4_gauscrysballsmear_wPY8MatchedFakesv4"

##__________________________________________________________________________________##
###my gauss core-extracted params+ fits of them, mean==1
#THYDIR="smearTheory/CT14NNLO_NLO_09.27.21_${FITTYPE}_wPFMET_RECOptLo10_GENptLo3_genbinsonly_TESTFIT_JOHNDSCB_GENYbinEQYbin_semifinalv4_spl3wgts_murmufHTp_JohnNPs_gausAt1smear_plots/"
#THYFILESTR="CT14NNLO_NLO_09.27.21_${FITTYPE}_wPFMET_RECOptLo10_GENptLo3_genbinsonly_TESTFIT_JOHNDSCB_GENYbinEQYbin_semifinalv4_spl3wgts_murmufHTp_JohnNPs_gausAt1smear_"
#OUTPUTTAG="09.27.21_Bayes_CT14NNLO_NLO_${FITTYPE}_murmufHTp_SMPbins_wJERSFs_semifinalv4_gausAt1smear_wPY8MatchedFakesv4"

##__________________________________________________________________________________##
##johns DSCB core-extracted params+ fits of them, mean==1
#THYDIR="smearTheory/CT14NNLO_NLO_09.27.21_${FITTYPE}_wPFMET_RECOptLo10_GENptLo3_genbinsonly_TESTFIT_JOHNDSCB_GENYbinEQYbin_semifinalv4_spl3wgts_murmufHTp_JohnNPs_gausAt1crysballsmear_plots/"
#THYFILESTR="CT14NNLO_NLO_09.27.21_${FITTYPE}_wPFMET_RECOptLo10_GENptLo3_genbinsonly_TESTFIT_JOHNDSCB_GENYbinEQYbin_semifinalv4_spl3wgts_murmufHTp_JohnNPs_gausAt1crysballsmear_"
#OUTPUTTAG="09.27.21_Bayes_CT14NNLO_NLO_${FITTYPE}_murmufHTp_SMPbins_wJERSFs_semifinalv4_gausAt1crysballsmear_wPY8MatchedFakesv4"

##__________________________________________________________________________________##
####this will be for when john is done with improving the DSCB fits; hopefully soonish
#THYDIR="smearTheory/CT14NNLO_NLO_09.27.21_${FITTYPE}_wPFMET_RECOptLo10_GENptLo3_genbinsonly_TESTFIT2_JOHNDSCB_GENYbinEQYbin_semifinalv4_spl3wgts_murmufHTp_JohnNPs_gauscrysballsmear_plots/"
#THYFILESTR="CT14NNLO_NLO_09.27.21_${FITTYPE}_wPFMET_RECOptLo10_GENptLo3_genbinsonly_TESTFIT2_JOHNDSCB_GENYbinEQYbin_semifinalv4_spl3wgts_murmufHTp_JohnNPs_gauscrysballsmear_"
#OUTPUTTAG="09.27.21_Bayes_CT14NNLO_NLO_${FITTYPE}_murmufHTp_SMPbins_wJERSFs_semifinalv4_TESTFIT2_gauscrysballsmear_wPY8MatchedFakesv4"

etabinStart=0
etabinEnd=4
applyNPcorrs=1
dosystunf=1

declare -a sysarr=("JER" "JEC")
declare -a sysarropt=("" "")
#declare -a sysarr=("JEC")
#declare -a sysarropt=("")
#declare -a sysarr=("JEC" "JER" "PDF" "NP")
#declare -a sysarropt=("" "" "updown" "updown")

sysarrlen=${#sysarr[@]}
#sysarrlen=1

#etabinStart=0
#etabinEnd=1
#applyNPcorrs=1
#dosystunf=1

#sysarrlen=${#sysarr[@]}

for (( j=${etabinStart}; j<${etabinEnd}; j++ ));
do
    for (( i=0; i<${sysarrlen}; i++ ));
    do
	echo ""
	echo "j=$j"
	echo "${OUTPUTTAG}_${sysarr[i]}${sysarropt[i]}sys"
	echo "${sysarr[i]}"
	echo "${sysarropt[i]}"
	./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} ${OUTPUTTAG}_${sysarr[i]}sys${sysarropt[i]} $j ${THYDIR} ${THYFILESTR} ${dosystunf} ${sysarr[i]} ${sysarropt[i]} ${applyNPcorrs}
    done
done

return





############ ------------- OLD *** BUT IMPORTANT *** ###############

#####GAUSS CORE BASED ON CRYSTAL BALL FIT PARAMS SMEARING [GAUSS ALWAYS CENTERED AT 1]
#DATADIR="01.05.21_outputCondor/ppData_HighPtJetTrig_ak4PFJets_01-05-21_jetPlots_0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"
#THYDIR="smearTheory/CT14NNLO_NLO_03.18.20_sigmu_geny_crysball_ptLo56_spl3wgts_gaussSmear_murmufHTp_JohnNPs_semifinalv4_gausAt1crysballsmear_wJERhistsmear_wresid_plots/"
#THYFILESTR="CT14NNLO_NLO_03.18.20_sigmu_geny_crysball_ptLo56_spl3wgts_gaussSmear_murmufHTp_JohnNPs_semifinalv4_gausAt1crysballsmear_wJERhistsmear_wresid_"
#OUTPUTTAG="01.05.21_Bayes_CT14NNLO_NLO_${FITTYPE}_murmufHTp_SMPbins_withgenjety_wJERSFs_SMPbins_wPFMET0p3_semifinalv4_gausAt1crysballsmear_wJERhistsmear_wresid"

#####GAUSS CORE BASED ON CRYSTAL BALL FIT PARAMS SMEARING
#DATADIR="01.05.21_outputCondor/ppData_HighPtJetTrig_ak4PFJets_01-05-21_jetPlots_0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"
#THYDIR="smearTheory/CT14NNLO_NLO_03.18.20_sigmu_geny_crysball_ptLo56_spl3wgts_gaussSmear_murmufHTp_JohnNPs_semifinalv4_gauscrysballsmear_wJERhistsmear_wresid_plots/"
#THYFILESTR="CT14NNLO_NLO_03.18.20_sigmu_geny_crysball_ptLo56_spl3wgts_gaussSmear_murmufHTp_JohnNPs_semifinalv4_gauscrysballsmear_wJERhistsmear_wresid_"
#OUTPUTTAG="01.05.21_Bayes_CT14NNLO_NLO_${FITTYPE}_murmufHTp_SMPbins_withgenjety_wJERSFs_SMPbins_wPFMET0p3_semifinalv4_gauscrysballsmear_wJERhistsmear_wresid"

#### WITH RESIDUALS NO JER HIST SMEARING (OLD WAY, LATEST NLO IN MY BINS)
#DATADIR="01.05.21_outputCondor/ppData_HighPtJetTrig_ak4PFJets_01-05-21_jetPlots_0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"
#THYDIR="smearTheory/CT14NNLO_NLO_03.18.20_sigmu_geny_spl3wgts_gaussSmear_murmufHTp_JohnNPs_semifinalv4_smeartests_noJERhistsmear_wresid_plots/"
#THYFILESTR="CT14NNLO_NLO_03.18.20_sigmu_geny_spl3wgts_gaussSmear_murmufHTp_JohnNPs_semifinalv4_smeartests_noJERhistsmear_wresid_"
#OUTPUTTAG="01.05.21_Bayes_CT14NNLO_NLO_${FITTYPE}_murmufHTp_SMPbins_withgenjety_wJERSFs_SMPbins_wPFMET0p3_semifinalv4_smeartests_noJERhistsmear_wresid"

#### WITH RESIDUALS WITH JER HIST SMEARING (NEWER WAY, LATEST NLO IN MY BINS)
#DATADIR="01.05.21_outputCondor/ppData_HighPtJetTrig_ak4PFJets_01-05-21_jetPlots_0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"
#THYDIR="smearTheory/CT14NNLO_NLO_03.18.20_sigmu_geny_spl3wgts_gaussSmear_murmufHTp_JohnNPs_semifinalv4_smeartests_wJERhistsmear_wresid_plots/"
#THYFILESTR="CT14NNLO_NLO_03.18.20_sigmu_geny_spl3wgts_gaussSmear_murmufHTp_JohnNPs_semifinalv4_smeartests_wJERhistsmear_wresid_"
#OUTPUTTAG="01.05.21_Bayes_CT14NNLO_NLO_${FITTYPE}_murmufHTp_SMPbins_withgenjety_wJERSFs_SMPbins_wPFMET0p3_semifinalv4_smeartests_wJERhistsmear_wresid"


######USED FOR PREAPPROVAL
#DATADIR="01.05.21_outputCondor/ppData_HighPtJetTrig_ak4PFJets_01-05-21_jetPlots_0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_semifinalv4"
###THYDIR="smearTheory/CT14NNLO_NLO_03.18.20_${FITTYPE}_geny_spl3wgts_gaussSmear_murmufHTp_JohnNPs_semifinalv4_plots/"
###THYFILESTR="CT14NNLO_NLO_03.18.20_${FITTYPE}_geny_spl3wgts_gaussSmear_murmufHTp_JohnNPs_semifinalv4_"
###OUTPUTTAG="01.05.21_Bayes_CT14NNLO_NLO_${FITTYPE}_murmufHTp_SMPbins_withgenjety_wJERSFs_SMPbins_wPFMET0p3_semifinalv4"
#THYDIR="smearTheory/CT14NNLO_NLO_03.18.20_${FITTYPE}_geny_spl3wgts_gaussSmear_murmufpt_JohnNPs_semifinalv4_plots/"
#THYFILESTR="CT14NNLO_NLO_03.18.20_${FITTYPE}_geny_spl3wgts_gaussSmear_murmufpt_JohnNPs_semifinalv4_"
#OUTPUTTAG="01.05.21_Bayes_CT14NNLO_NLO_${FITTYPE}_murmufpt_SMPbins_withgenjety_wJERSFs_SMPbins_wPFMET0p3_semifinalv4"









































################### OLD OLD BELOW HERE; LIKELY IRRELEVANT AND/OR DEPRECATED




#####NNPDF_NLO_06.25.19_sigmu_noJERscales_spl3wgts_gaussSmear_NLO_plots
#
#FITTYPE="sigmu"
##DATADIR="11.26.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_11-26-19_jetPlots_0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinalv2"
#DATADIR="11.26.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_11-26-19_jetPlots_0.0eta2.0_SMPbins_JECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh"
##THYDIR="smearTheory/NNPDF_LO_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_plots/"
##THYFILESTR="NNPDF_LO_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_"
#THYDIR="smearTheory/NNPDF_NLO_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_plots/"
#THYFILESTR="NNPDF_NLO_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_"
##OUTPUTTAG="11.26.19_Bayes_NNPDF_NLO_${FITTYPE}_noJERscales_withLumiCorr_SMPbins_withJECsysv2_ptLo56_semifinalv2"
#OUTPUTTAG="11.26.19_Bayes_NNPDF_NLO_${FITTYPE}_noJERscales_withLumiCorr_SMPbins_withJECsysv2_ptLo56_semifinalv2_LowHLT40Thresh"



###vanilla no sys
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_nosys"  0    ${THYDIR} ${THYFILESTR}  0 "" "" 1
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_nosys"  1    ${THYDIR} ${THYFILESTR}  0 "" "" 1
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_nosys"  2    ${THYDIR} ${THYFILESTR}  0 "" "" 1
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_nosys"  3    ${THYDIR} ${THYFILESTR}  0 "" "" 1 
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_nosys_NLOMC_wNP_??eta??.*"
#return
#
##JEC sys
./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JECsysv2"  0    ${THYDIR} ${THYFILESTR}  1 "JEC" "" 1
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JECsysv2"  1    ${THYDIR} ${THYFILESTR}  1 "JEC" "" 1
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JECsysv2"  2    ${THYDIR} ${THYFILESTR}  1 "JEC" "" 1
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JECsysv2"  3    ${THYDIR} ${THYFILESTR}  1 "JEC" "" 1
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_JECsys_NLOMC_wNP_??eta??.*"
#return
##
#####JER sys
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JERsys"  0    ${THYDIR} ${THYFILESTR}  1 "JER" "" 1
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JERsys"  1    ${THYDIR} ${THYFILESTR}  1 "JER" "" 1
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JERsys"  2    ${THYDIR} ${THYFILESTR}  1 "JER" "" 1
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JERsys"  3    ${THYDIR} ${THYFILESTR}  1 "JER" "" 1
###scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_JERsys_NLOMC_wNP_??eta??.*"
###return
##
#
#####NP sys12, alt NP choices, HERWIG EE4C/PYTHIA8
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsys12" 0   ${THYDIR} ${THYFILESTR}  1 "NP" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsys12" 1   ${THYDIR} ${THYFILESTR}  1 "NP" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsys12" 2   ${THYDIR} ${THYFILESTR}  1 "NP" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsys12" 3   ${THYDIR} ${THYFILESTR}  1 "NP" "12"
####scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_NPsys12_NLOMC_wNP_??eta??.*"
####return
###
#####NP sysupdown, shifted HERWIG EE5C NP up/down, 1 sigma
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsysupdown" 0 ${THYDIR} ${THYFILESTR}  1 "NP" "updown"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsysupdown" 1 ${THYDIR} ${THYFILESTR}  1 "NP" "updown"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsysupdown" 2 ${THYDIR} ${THYFILESTR}  1 "NP" "updown"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsysupdown" 3 ${THYDIR} ${THYFILESTR}  1 "NP" "updown"
####scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_NPsysupdown_NLOMC_wNP_??eta??.*"
####return
###
#####PDF sys12, alt PDF choices, CT14/HERA
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsys12" 0   ${THYDIR} ${THYFILESTR}  1 "PDF" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsys12" 1   ${THYDIR} ${THYFILESTR}  1 "PDF" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsys12" 2   ${THYDIR} ${THYFILESTR}  1 "PDF" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsys12" 3   ${THYDIR} ${THYFILESTR}  1 "PDF" "12"
####scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_PDFsys12_NLOMC_wNP_??eta??.*"
####return
###
#####PDF sysupdown, shifted NNPDF PDF up/down, 6 Pt scale + PDF errs, 1 sigma
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsysupdown" 0 ${THYDIR} ${THYFILESTR}  1 "PDF" "updown" 1
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsysupdown" 1 ${THYDIR} ${THYFILESTR}  1 "PDF" "updown" 1
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsysupdown" 2 ${THYDIR} ${THYFILESTR}  1 "PDF" "updown" 1
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsysupdown" 3 ${THYDIR} ${THYFILESTR}  1 "PDF" "updown" 1
###scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_PDFsysupdown_NLOMC_wNP_??eta??.*"
###return

#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_*sys*_NLOMC_wNP_??eta??.*"
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_*sys*_NLOMC_wNP_??eta??_respMat.*"
#return
###-----

#
#
#
#
#
#
############### SMEARED GEN PY 8 LO + HERWIG EE5C LO NPCs ###############
#####PY8_03.14.19_sigmu_noJERscales_spl3wgts_gaussSmear_NLO_plots
##
#FITTYPE="sigmu"
#DATADIR="06.25.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_06-25-19_jetPlots_0.0eta2.0_SMPbins_withJECsys_semifinal"
#THYDIR="smearTheory/PY8_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_NLO_plots/"
#THYFILESTR="PY8_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_"
#OUTPUTTAG="06.25.19_Bayes_PY8_LO_${FITTYPE}_noJERscales_semifinal"
##
##
####vanilla no sys
###./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_nosys"  0    ${THYDIR} ${THYFILESTR}  0 "" ""
###./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_nosys"  1    ${THYDIR} ${THYFILESTR}  0 "" ""
###./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_nosys"  2    ${THYDIR} ${THYFILESTR}  0 "" ""
###./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_nosys"  3    ${THYDIR} ${THYFILESTR}  0 "" ""
###scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_nosys_NLOMC_wNP_??eta??.*"
###return
##
####JEC sys
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JECsys"  0    ${THYDIR} ${THYFILESTR}  1 "JEC" ""
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JECsys"  1    ${THYDIR} ${THYFILESTR}  1 "JEC" ""
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JECsys"  2    ${THYDIR} ${THYFILESTR}  1 "JEC" ""
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JECsys"  3    ${THYDIR} ${THYFILESTR}  1 "JEC" ""
###scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_JECsys_NLOMC_wNP_??eta??.*"
###return
##
####JER sys
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JERsys"  0    ${THYDIR} ${THYFILESTR}  1 "JER" ""
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JERsys"  1    ${THYDIR} ${THYFILESTR}  1 "JER" ""
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JERsys"  2    ${THYDIR} ${THYFILESTR}  1 "JER" ""
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_JERsys"  3    ${THYDIR} ${THYFILESTR}  1 "JER" ""
###scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_JERsys_NLOMC_wNP_??eta??.*"
###return
##
####NP sys12, alt NP choices, HERWIG EE4C/PYTHIA8
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsys12" 0   ${THYDIR} ${THYFILESTR}  1 "NP" "12"
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsys12" 1   ${THYDIR} ${THYFILESTR}  1 "NP" "12"
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsys12" 2   ${THYDIR} ${THYFILESTR}  1 "NP" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsys12" 3   ${THYDIR} ${THYFILESTR}  1 "NP" "12"
###scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_NPsys12_NLOMC_wNP_??eta??.*"
###return
##
####NP sysupdown, shifted HERWIG EE5C NP up/down, 1 sigma
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsysupdown" 0 ${THYDIR} ${THYFILESTR}  1 "NP" "updown"
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsysupdown" 1 ${THYDIR} ${THYFILESTR}  1 "NP" "updown"
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsysupdown" 2 ${THYDIR} ${THYFILESTR}  1 "NP" "updown"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_NPsysupdown" 3 ${THYDIR} ${THYFILESTR}  1 "NP" "updown"
###scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_NPsysupdown_NLOMC_wNP_??eta??.*"
###return
##
####PDF sys12, alt PDF choices, CT14/HERA
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsys12" 0   ${THYDIR} ${THYFILESTR}  1 "PDF" "12"
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsys12" 1   ${THYDIR} ${THYFILESTR}  1 "PDF" "12"
##./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsys12" 2   ${THYDIR} ${THYFILESTR}  1 "PDF" "12"
#./bayesUnfoldDataSpectra_wNLO_etabin.exe ${DATADIR} "${OUTPUTTAG}_PDFsys12" 3   ${THYDIR} ${THYFILESTR}  1 "PDF" "12"
###scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_PDFsys12_NLOMC_wNP_??eta??.*"
###return
##
#scp2Serin "output/unfoldDataSpectra/ak4PFJets_wjtID_anabins_${OUTPUTTAG}_*sys*_NLOMC_wNP_??eta??.*"
#return
####-----































############# NNPDF NLO + HERWIG EE5C LO NPCs ###############
####NNPDF_06.25.19_sigmu_noJERscales_spl3wgts_gaussSmear_NLO_plots

#FITTYPE="sigmu"
#DATADIR="08.06.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_08-06-19_jetPlots_0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal"
#THYDIR="smearTheory/NNPDF_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_NLO_plots/"
#THYFILESTR="NNPDF_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_"
#OUTPUTTAG="08.06.19_Bayes_NNPDF_NLO_${FITTYPE}_noJERscales_withLumiCorr_semifinal"

#FITTYPE="sigmu"
#DATADIR="08.06.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_08-06-19_jetPlots_0.0eta2.0_SMPbins_withJECsys_withLumiCorr_semifinal"
#THYDIR="smearTheory/NNPDF_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_NLO_plots/"
#THYFILESTR="NNPDF_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_"
#OUTPUTTAG="08.06.19_Bayes_NNPDF_NLO_${FITTYPE}_noJERscales_withLumiCorr_SMPbins_ptLo64_semifinal"

#FITTYPE="sigmu"
#DATADIR="08.06.19_outputCondor/ppData_HighPtJetTrig_ak4PFJets_08-06-19_jetPlots_0.0eta2.0_JOHNbins_withJECsys_withLumiCorr_semifinal"
#THYDIR="smearTheory/NNPDF_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_NLO_plots/"
#THYFILESTR="NNPDF_06.25.19_${FITTYPE}_semifinal_noJERscales_spl3wgts_gaussSmear_"
#OUTPUTTAG="08.06.19_Bayes_NNPDF_NLO_${FITTYPE}_noJERscales_withLumiCorr_JOHNbins_semifinal"




#THYDIR="smearTheory/CT14NNLO_NLO_v3_03.18.20_${FITTYPE}_geny_noJERscales_spl3wgts_gaussSmear_TEST_murmufpt_plots/"
#THYFILESTR="CT14NNLO_NLO_v3_03.18.20_${FITTYPE}_geny_noJERscales_spl3wgts_gaussSmear_TEST_murmufpt_"
#OUTPUTTAG="03.18.20_Bayes_CT14NNLO_NLO_v3_${FITTYPE}_noJERscales_SMPbins_withjety_ptLo56_semifinalv3_REBIN_NODATACOVMAT_murmufpt_TEST"
#THYDIR="smearTheory/CT14NNLO_NLO_v3_03.18.20_${FITTYPE}_geny_noJERscales_spl3wgts_gaussSmear_TEST_murmufpt1_plots/"
#THYFILESTR="CT14NNLO_NLO_v3_03.18.20_${FITTYPE}_geny_noJERscales_spl3wgts_gaussSmear_TEST_murmufpt1_"
#OUTPUTTAG="03.18.20_Bayes_CT14NNLO_NLO_v3_${FITTYPE}_noJERscales_SMPbins_withjety_ptLo56_semifinalv3_REBIN_NODATACOVMAT_murmufpt1_TEST"

#THYDIR="smearTheory/CT14NNLO_NLO_v3_03.18.20_${FITTYPE}_geny_spl3wgts_gaussSmear_murmufpt_plots/"
#THYFILESTR="CT14NNLO_NLO_v3_03.18.20_${FITTYPE}_geny_spl3wgts_gaussSmear_murmufpt_"
#OUTPUTTAG="03.18.20_Bayes_CT14NNLO_NLO_v3_${FITTYPE}_SMPbins_withjety_semifinalv3_REBIN_NODATACOVMAT_murmufpt"
#OUTPUTTAG="03.18.20_Bayes_CT14NNLO_NLO_v3_${FITTYPE}_defSMPbins_withjety_semifinalv3_wDATACOVMAT_murmufpt"
#OUTPUTTAG="03.18.20_Bayes_CT14NNLO_NLO_v3_${FITTYPE}_defSMPbins2_withjety_semifinalv3_wDATACOVMAT_murmufpt"
#THYDIR="smearTheory/CT14NNLO_NLO_v3_03.18.20_${FITTYPE}_geny_spl3wgts_gaussSmear_murmufpt1_plots/"
#THYFILESTR="CT14NNLO_NLO_v3_03.18.20_${FITTYPE}_geny_spl3wgts_gaussSmear_murmufpt1_"
#OUTPUTTAG="03.18.20_Bayes_CT14NNLO_NLO_v3_${FITTYPE}_SMPbins_withjety_semifinalv3_REBIN_NODATACOVMAT_murmufpt1"

