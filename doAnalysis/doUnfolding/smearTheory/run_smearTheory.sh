#!/bin/bash

echo ""
echo "compiling JER_gaussCore_smearTheoryCurves"
echo ""

rooUnfoldCompile smearTheorySpectra_gaussCoreJER.C

echo ""
echo "RUN NLO smearing"
echo ""

./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_CT10nlo_HS
#./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_CT14nlo_HS
#./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS
#./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS
#./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC

echo ""
echo "DONE NLO smearing."
echo "DONE NLO smearing. scp2Serin"
echo ""

scp2Serin fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_00eta20.root



return

