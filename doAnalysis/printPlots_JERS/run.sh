#!/bin/bash



echo ""
echo "compiling printPlots_ppMC_JERS_etabin.C"
echo ""

rootcompile printPlots_ppMC_JERS_etabin.C

echo ""
echo "done compiling. Running!"
echo ""

declare -a ybinstrARR=( "00eta05" "05eta10" "10eta15" "15eta20" )
ybinstrARRlen=${#ybinstrARR[@]}

##### echo "source run_printPlots_ppMC_JERS_etabin.sh <R> <condorDir> <date_output> <descIn> <descOut> <JERplots> <MCeffplots> <absetabin> <fitquant>"

##FITTYPE="sigmu"
#FITTYPE="sigma"
#JETR="4"
#MM="07"
#DD="22"
#YY="21"
#MMDDYY="${MM}-${DD}-${YY}"
#CONDORDIR="${MM}.${DD}.${YY}_outputCondor"
####### ak4PF_PY8JER_00eta05_07.22.21_sigmu_noPFMET_minlopTcuts_genbinsonly_semifinalv4.root
##INDESC="0.0eta2.5_SMPbins_noPFMETFrac_genbinsonly_semifinalv4"
#INDESC="0.0eta2.5_SMPbins_wPFMETFrac_genbinsonly_semifinalv4"
##INDESC="0.0eta2.5_SMPbins_wPFMETFrac0p5_genbinsonly_semifinalv4"
##INDESC="0.0eta2.5_SMPbins_wPFMETFrac0p5_genDR0p2_genbinsonly_semifinalv4"
##INDESC="0.0eta2.5_SMPbins_wPFMETFrac_genDR0p2_genbinsonly_semifinalv4"
##INDESC="0.0eta2.5_SMPbins_wPFMETFrac_genDR0p1_genbinsonly_semifinalv4"
##INDESC="0.0eta2.5_SMPbins_wPFMET_GENRECOmatchREDO_genbinsonly_semifinalv4"
##INDESC="0.0eta2.5_SMPbins_wPFMET_GENRECOmatchREDOwDR0p2_genbinsonly_semifinalv4"
##INDESC="0.0eta2.5_SMPbins_wPFMET_evtvr0p1_GENRECOmatchREDOwDR0p2_genbinsonly_semifinalv4"
##INDESC="0.0eta2.5_SMPbins_wPFMET_Nvtxeq1_evtvr0p1_GENRECOmatchREDOwDR0p2_genbinsonly_semifinalv4"
##INDESC="0.0eta2.5_SMPbins_wPFMET_novzwgt_Nvtxeq1_evtvr0p1_GENRECOmatchREDOwDR0p2_genbinsonly_semifinalv4"
##OUTDESC="${MM}.${DD}.${YY}_${FITTYPE}_noPFMET_minlopTcuts_genbinsonly_semifinalv4"
##OUTDESC="${MM}.${DD}.${YY}_${FITTYPE}_wPFMET_minlopTcuts_genbinsonly_semifinalv4"
#OUTDESC="${MM}.${DD}.${YY}_${FITTYPE}_wPFMET_minlopTcuts_genbinsonly_TESTFIT_JOHNDSCB_semifinalv4"
##OUTDESC="${MM}.${DD}.${YY}_${FITTYPE}_wPFMETFrac0p5_minlopTcuts_genbinsonly_semifinalv4"
##OUTDESC="${MM}.${DD}.${YY}_${FITTYPE}_wPFMETFrac0p5_genDR0p2_minlopTcuts_genbinsonly_semifinalv4"
##OUTDESC="${MM}.${DD}.${YY}_${FITTYPE}_wPFMET_minlopTcuts_genDR0p2_genbinsonly_semifinalv4"
##OUTDESC="${MM}.${DD}.${YY}_${FITTYPE}_wPFMET_minlopTcuts_genDR0p1_genbinsonly_semifinalv4"
##OUTDESC="${MM}.${DD}.${YY}_${FITTYPE}_wPFMET_GENRECOmatchREDO_minlopTcuts_genbinsonly_semifinalv4"
##OUTDESC="${MM}.${DD}.${YY}_${FITTYPE}_wPFMET_GENRECOmatchREDOwDR0p2_minlopTcuts_genbinsonly_semifinalv4"
##OUTDESC="${MM}.${DD}.${YY}_${FITTYPE}_wPFMET_evtvr0p1_GENRECOmatchREDOwDR0p2_minlopTcuts_genbinsonly_semifinalv4"
##OUTDESC="${MM}.${DD}.${YY}_${FITTYPE}_wPFMET_Nvtxeq1_evtvr0p1_GENRECOmatchREDOwDR0p2_minlopTcuts_genbinsonly_semifinalv4"
##OUTDESC="${MM}.${DD}.${YY}_${FITTYPE}_wPFMET_novzwgt_Nvtxeq1_evtvr0p1_GENRECOmatchREDOwDR0p2_minlopTcuts_genbinsonly_semifinalv4"


#FITTYPE="sigmu"
FITTYPE="sigma"
JETR="4"
MM="09"
DD="27"
YY="21"
MMDDYY="${MM}-${DD}-${YY}"
CONDORDIR="${MM}.${DD}.${YY}_outputCondor"
#INDESC="0.0eta2.5_SMPbins_wPFMET_RECOptLo10_GENptLo3_genbinsonly_semifinalv4"
INDESC="0.0eta2.5_SMPbins_wPFMET_RECOptLo10_GENptLo3_genbinsonly_noMatchedYbinFakesAndMiss_semifinalv4"
#OUTDESC="${MM}.${DD}.${YY}_${FITTYPE}_wPFMET_RECOptLo10_GENptLo3_genbinsonly_TESTFIT_JOHNDSCB_semifinalv4"
OUTDESC="${MM}.${DD}.${YY}_${FITTYPE}_wPFMET_RECOptLo10_GENptLo3_genbinsonly_TESTFIT_JOHNDSCB_GENYbinEQYbin_LOWPTFIX_semifinalv4"
#OUTDESC="${MM}.${DD}.${YY}_${FITTYPE}_wPFMET_RECOptLo10_GENptLo3_genbinsonly_TESTFIT_JOHNDSCB_GENYbinEQYbin_semifinalv4_DATAPOINTS_ONLY"
### NOTE for outdesc listed below this line; first ybin is diff from TESTFIT jobs above, next three are the same.
#OUTDESC="${MM}.${DD}.${YY}_${FITTYPE}_wPFMET_RECOptLo10_GENptLo3_genbinsonly_TESTFIT2_JOHNDSCB_GENYbinEQYbin_semifinalv4"
BINTYPE="geny"

ybinStart=0
ybinEnd=4
for (( YBIN=${ybinStart}; YBIN<${ybinEnd}; YBIN++ ));
do
    source run_printPlots_ppMC_JERS_etabin.sh $JETR $CONDORDIR $MMDDYY $INDESC ${ybinstrARR[YBIN]}_$OUTDESC "1" "0" "$YBIN" $FITTYPE $BINTYPE
done


return

#source run_printPlots_ppMC_JERS_etabin.sh "4" "01.05.21_outputCondor" "01-05-21" "0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_tightJetID_semifinalv4" "00eta05_01.05.21_${FITTYPE}_geny" "1" "0" "0" "$FITTYPE" "geny"
#source run_printPlots_ppMC_JERS_etabin.sh "4" "01.05.21_outputCondor" "01-05-21" "0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_tightJetID_semifinalv4" "05eta10_01.05.21_${FITTYPE}_geny" "1" "0" "1" "$FITTYPE" "geny"
#source run_printPlots_ppMC_JERS_etabin.sh "4" "01.05.21_outputCondor" "01-05-21" "0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_tightJetID_semifinalv4" "10eta15_01.05.21_${FITTYPE}_geny" "1" "0" "2" "$FITTYPE" "geny"
#source run_printPlots_ppMC_JERS_etabin.sh "4" "01.05.21_outputCondor" "01-05-21" "0.0eta2.5_SMPbins_withjety_wHBHEIsoNoise_wPFMETFrac0p3_tightJetID_semifinalv4" "15eta20_01.05.21_${FITTYPE}_geny" "1" "0" "3" "$FITTYPE" "geny"


###ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_03-18-20_JERS_0.0eta2.5_SMPbins_withjety_semifinalv3
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_semifinalv3" "00eta05_03.18.20_${FITTYPE}_recy" "1" "0" "0" "$FITTYPE" "y"
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_semifinalv3" "05eta10_03.18.20_${FITTYPE}_recy" "1" "0" "1" "$FITTYPE" "y"
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_semifinalv3" "10eta15_03.18.20_${FITTYPE}_recy" "1" "0" "2" "$FITTYPE" "y"
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_semifinalv3" "15eta20_03.18.20_${FITTYPE}_recy" "1" "0" "3" "$FITTYPE" "y"
																					
source run_printPlots_ppMC_JERS_etabin.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_semifinalv3" "00eta05_03.18.20_${FITTYPE}_geny" "1" "0" "0" "$FITTYPE" "geny"
source run_printPlots_ppMC_JERS_etabin.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_semifinalv3" "05eta10_03.18.20_${FITTYPE}_geny" "1" "0" "1" "$FITTYPE" "geny"
source run_printPlots_ppMC_JERS_etabin.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_semifinalv3" "10eta15_03.18.20_${FITTYPE}_geny" "1" "0" "2" "$FITTYPE" "geny"
source run_printPlots_ppMC_JERS_etabin.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_semifinalv3" "15eta20_03.18.20_${FITTYPE}_geny" "1" "0" "3" "$FITTYPE" "geny"

#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_semifinalv3" "00eta05_03.18.20_${FITTYPE}_receta" "1" "0" "0" "$FITTYPE" "eta"
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_semifinalv3" "05eta10_03.18.20_${FITTYPE}_receta" "1" "0" "1" "$FITTYPE" "eta"
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_semifinalv3" "10eta15_03.18.20_${FITTYPE}_receta" "1" "0" "2" "$FITTYPE" "eta"
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_semifinalv3" "15eta20_03.18.20_${FITTYPE}_receta" "1" "0" "3" "$FITTYPE" "eta"
																					       
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_semifinalv3" "00eta05_03.18.20_${FITTYPE}_geneta" "1" "0" "0" "$FITTYPE" "geneta"
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_semifinalv3" "05eta10_03.18.20_${FITTYPE}_geneta" "1" "0" "1" "$FITTYPE" "geneta"
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_semifinalv3" "10eta15_03.18.20_${FITTYPE}_geneta" "1" "0" "2" "$FITTYPE" "geneta"
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.18.20_outputCondor" "03-18-20" "0.0eta2.5_SMPbins_withjety_semifinalv3" "15eta20_03.18.20_${FITTYPE}_geneta" "1" "0" "3" "$FITTYPE" "geneta"

#echo "scp2Serin"																             
#scp2Serin "output/ak4PF_PY8JER_??eta??_03.18.20_sigm*.*"  







echo ""
echo "done running printPlots_ppMC_JERS_etabin"
echo ""

##############################################
return
##############################################































#############################################################################################################
########################################### NOT DEPRECATED BUT OLD ##########################################
#############################################################################################################



##ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_03-14-19_JERS_0.0eta3.0_ptLo28_JER

#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "00eta05_03.14.19_sigma" "1" "0" "0" "sigma" 
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "05eta10_03.14.19_sigma" "1" "0" "1" "sigma" 
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "10eta15_03.14.19_sigma" "1" "0" "2" "sigma" 
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "15eta20_03.14.19_sigma" "1" "0" "3" "sigma" 

#source run_printPlots_ppMC_JERS_etabin.sh "4" "06.25.19_outputCondor" "06-25-19" "0.0eta2.0_semifinal" "00eta05_06.25.19_sigmu_semifinal" "1" "0" "0" "sigmu" 
#source run_printPlots_ppMC_JERS_etabin.sh "4" "06.25.19_outputCondor" "06-25-19" "0.0eta2.0_semifinal" "05eta10_06.25.19_sigmu_semifinal" "1" "0" "1" "sigmu" 
#source run_printPlots_ppMC_JERS_etabin.sh "4" "06.25.19_outputCondor" "06-25-19" "0.0eta2.0_semifinal" "10eta15_06.25.19_sigmu_semifinal" "1" "0" "2" "sigmu" 
source run_printPlots_ppMC_JERS_etabin.sh "4" "06.25.19_outputCondor" "06-25-19" "0.0eta2.0_semifinal" "15eta20_06.25.19_sigmu_semifinal" "1" "0" "3" "sigmu" 

###ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_06-07-19_JERS_0.0eta2.0_JER
#
#source run_printPlots_ppMC_JERS_etabin.sh "4" "06.07.19_outputCondor" "06-07-19" "0.0eta2.0_JER" "00eta05_06.07.19_sigma" "1" "0" "0" "sigma" 
#source run_printPlots_ppMC_JERS_etabin.sh "4" "06.07.19_outputCondor" "06-07-19" "0.0eta2.0_JER" "05eta10_06.07.19_sigma" "1" "0" "1" "sigma" 
#source run_printPlots_ppMC_JERS_etabin.sh "4" "06.07.19_outputCondor" "06-07-19" "0.0eta2.0_JER" "10eta15_06.07.19_sigma" "1" "0" "2" "sigma" 
#source run_printPlots_ppMC_JERS_etabin.sh "4" "06.07.19_outputCondor" "06-07-19" "0.0eta2.0_JER" "15eta20_06.07.19_sigma" "1" "0" "3" "sigma" 
#
#
#source run_printPlots_ppMC_JERS_etabin.sh "4" "06.07.19_outputCondor" "06-07-19" "0.0eta2.0_JER" "00eta05_06.07.19_sigmu" "1" "0" "0" "sigmu" 
#source run_printPlots_ppMC_JERS_etabin.sh "4" "06.07.19_outputCondor" "06-07-19" "0.0eta2.0_JER" "05eta10_06.07.19_sigmu" "1" "0" "1" "sigmu" 
#source run_printPlots_ppMC_JERS_etabin.sh "4" "06.07.19_outputCondor" "06-07-19" "0.0eta2.0_JER" "10eta15_06.07.19_sigmu" "1" "0" "2" "sigmu" 
#source run_printPlots_ppMC_JERS_etabin.sh "4" "06.07.19_outputCondor" "06-07-19" "0.0eta2.0_JER" "15eta20_06.07.19_sigmu" "1" "0" "3" "sigmu" 
#
#
#echo "scp2Serin"																             
#scp2Serin "output/ak4PF_PY8JER_??eta??_06.07.19_sigm*.*"										             



###ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_03-14-19_JERS_0.0eta3.0_ptLo28_JER
#
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "00eta05_03.14.19_sigma" "1" "0" "0" "sigma" 
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "05eta10_03.14.19_sigma" "1" "0" "1" "sigma" 
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "10eta15_03.14.19_sigma" "1" "0" "2" "sigma" 
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "15eta20_03.14.19_sigma" "1" "0" "3" "sigma" 
#
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "00eta05_03.14.19_sigmu" "1" "0" "0" "sigmu" 
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "05eta10_03.14.19_sigmu" "1" "0" "1" "sigmu" 
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "10eta15_03.14.19_sigmu" "1" "0" "2" "sigmu" 
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "15eta20_03.14.19_sigmu" "1" "0" "3" "sigmu" 
#
#
#echo "scp2Serin"																             
#scp2Serin "output/ak4PF_PY8JER_??eta??_03.14.19_sigm*.*"										             







##ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_03-14-19_JERS_0.0eta3.0_ptLo28_JER

#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "00eta05_ptLo28_drCut_geny_sigma" "1" "0" "0" "sigma" 
#echo "scp2Serin"																             
#scp2Serin "output/ak4PF_PY8JER_00eta05_ptLo28_drCut_geny_sigma.*"										             

																		             
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "05eta10_ptLo28_drCut_geny_sigma" "1" "0" "1" "sigma"
#echo "scp2Serin"																             
#scp2Serin "output/ak4PF_PY8JER_05eta10_ptLo28_drCut_geny_sigma.*"										             
																		             
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "10eta15_ptLo28_drCut_geny_sigma" "1" "0" "2" "sigma"  
#echo "scp2Serin"																             
#scp2Serin "output/ak4PF_PY8JER_10eta15_ptLo28_drCut_geny_sigma.*"										             
																		             
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "15eta20_ptLo28_drCut_geny_sigma" "1" "0" "3" "sigma" 
#echo "scp2Serin"
#scp2Serin "output/ak4PF_PY8JER_15eta20_ptLo28_drCut_geny_sigma.*"

#echo "scp2Serin"
#scp2Serin  "output/ak4PF_PY8JER_??eta??_ptLo28_drCut_geny_sigma.*"

#return

#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "00eta05_ptLo28_drCut_geny_sigmu" "1" "0" "0" "sigmu"
#echo "scp2Serin"																             
#scp2Serin "output/ak4PF_PY8JER_00eta05_ptLo28_drCut_geny_sigmu.*"										             
																		             
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "05eta10_ptLo28_drCut_geny_sigmu" "1" "0" "1" "sigmu" 
#echo "scp2Serin"																             
#scp2Serin "output/ak4PF_PY8JER_05eta10_ptLo28_drCut_geny_sigmu.*"										             
																		             
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "10eta15_ptLo28_drCut_geny_sigmu" "1" "0" "2" "sigmu" 
#echo "scp2Serin"																             
#scp2Serin "output/ak4PF_PY8JER_10eta15_ptLo28_drCut_geny_sigmu.*"										             
																		             
#source run_printPlots_ppMC_JERS_etabin.sh "4" "03.14.19_outputCondor" "03-14-19" "0.0eta3.0_ptLo28_JER" "15eta20_ptLo28_drCut_geny_sigmu" "1" "0" "3" "sigmu" 
#echo "scp2Serin"
#scp2Serin "output/ak4PF_PY8JER_15eta20_ptLo28_drCut_geny_sigmu.*"

#echo "scp2Serin"
#scp2Serin  "output/ak4PF_PY8JER_??eta??_ptLo28_drCut_geny_sigmu.*"

#echo "scp2Serin"
#scp2Serin  "output/ak4PF_PY8JER_??eta??_ptLo28_drCut_geny_sigm*.*"


#ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-18-19_JERS_0.0eta3.0_JER
#source run_printPlots_ppMC_JERS_etabin.sh "4" "02.18.19_outputCondor" "02-18-19" "0.0eta3.0_JER" "00eta05_semifinal" "1" "0" "0"         #JER any |y|bin
#source run_printPlots_ppMC_JERS_etabin.sh "4" "02.18.19_outputCondor" "02-18-19" "0.0eta3.0_JER" "05eta10_semifinal" "1" "0" "1"        #JER any |y|bin
#source run_printPlots_ppMC_JERS_etabin.sh "4" "02.18.19_outputCondor" "02-18-19" "0.0eta3.0_JER" "10eta15_semifinal" "1" "0" "2"        #JER any |y|bin
#source run_printPlots_ppMC_JERS_etabin.sh "4" "02.18.19_outputCondor" "02-18-19" "0.0eta3.0_JER" "15eta20_semifinal" "1" "0" "3"        #JER any |y|bin
#source run_printPlots_ppMC_JERS_etabin.sh "4" "02.18.19_outputCondor" "02-18-19" "0.0eta3.0_JER" "20eta25_semifinal" "1" "0" "4"        #JER any |y|bin
#source run_printPlots_ppMC_JERS_etabin.sh "4" "02.18.19_outputCondor" "02-18-19" "0.0eta3.0_JER" "25eta30_semifinal" "1" "0" "5"        #JER any |y|bin

#echo "scp2Serin"
#scp2Serin "output/ak4PF_PY8JER_00eta05_semifinal.*"
#scp2Serin "output/ak4PF_PY8JER_??eta??_semifinal.*"























#####################################################################################
############################### DEPRECATED ##########################################
#####################################################################################
#echo ""
#echo "compiling printPlots_ppMC_JERS.C"
#echo ""
#
#rootcompile printPlots_ppMC_JERS.C
#
#echo ""
#echo "done compiling. Running!"
#echo ""

#source run_printPlots_ppMC_JERS.sh "4" "08.07.18_outputCondor" "08-07-18" "0.0eta2.0" "00eta20"

#source run_printPlots_ppMC_JERS.sh "4" "11.06.18_outputCondor" "11-06-18" "0.0eta0.5" "00eta05"
#source run_printPlots_ppMC_JERS.sh "4" "11.06.18_outputCondor" "11-06-18" "0.5eta1.0" "05eta10"
#source run_printPlots_ppMC_JERS.sh "4" "11.06.18_outputCondor" "11-06-18" "1.0eta1.5" "10eta15"
#source run_printPlots_ppMC_JERS.sh "4" "11.06.18_outputCondor" "11-06-18" "1.5eta2.0" "15eta20"


#source run_printPlots_ppMC_JERS.sh "4" "11.14.18_outputCondor" "11-14-18" "0.0eta2.0" "00eta20" "1" "0"         #JER only
#source run_printPlots_ppMC_JERS.sh "4" "11.14.18_outputCondor" "11-14-18" "0.0eta2.0_unf" "00eta20_unf" "0" "1" #MCEff only

## usage
## source run_printPlots_ppMC_JERS.sh <radius> <condordir> <date input> <input etabin tag> <output etabin desc> <doJERplots> <doMCEffplots> <abs eta bin>

#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta2.0" "00eta20" "1" "0" "-1"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta2.0" "00eta05" "1" "0" "0"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta2.0" "05eta10" "1" "0" "1"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta2.0" "10eta15" "1" "0" "2"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta2.0" "15eta20" "1" "0" "3"        #JER 0-2only


#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta3.0" "00eta05" "1" "0" "0"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta3.0" "05eta10" "1" "0" "1"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta3.0" "10eta15" "1" "0" "2"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta3.0" "15eta20" "1" "0" "3"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta3.0" "20eta25" "1" "0" "4"        #JER 0-2only
#source run_printPlots_ppMC_JERS.sh "4" "01.06.19_outputCondor" "01-06-19" "0.0eta3.0" "25eta30" "1" "0" "5"        #JER 0-2only

#echo ""
#echo "done running printPlots_ppMC_JERS"
#echo ""

#echo ""
#echo "scp2Serin"
#echo ""

#scp2Serin output/ak4PF\*
