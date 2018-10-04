#!/bin/bash

echo ""
echo "compiling smearTheorySpectra_gaussCoreJER"
echo ""
rooUnfoldCompile smearTheorySpectra_gaussCoreJER.C

source run_smearTheory.sh


#echo ""
#echo "compiling smearTheorySpectra_gaussCoreJER_singleYbin"
#echo ""
#rooUnfoldCompile smearTheorySpectra_gaussCoreJER_singleYbin.C
#source run_smearTheory.sh  1
#source run_smearTheory.sh  1 all