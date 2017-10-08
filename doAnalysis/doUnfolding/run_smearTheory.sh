#!/bin/bash

echo ""
echo "compiling JER_gaussCore_smearTheoryCurves"
echo ""

rooUnfoldCompile smearTheorySpectra_gaussCoreJER.C

echo ""
echo "running smearing"
echo ""

./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_CT10nlo_HS

echo ""
echo "done."
echo ""

return












./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_CT14nlo_HS
./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS
./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS
./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC

echo ""
echo "done."
echo ""

return