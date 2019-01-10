#!/bin/bash

#echo ""
#echo "compiling smearTheorySpectra_gaussCoreJER"
#echo ""
#rooUnfoldCompile smearTheorySpectra_gaussCoreJER.C
#source run_smearTheory.sh   ##00eta20



echo ""
echo "compiling smearTheorySpectra_gaussCoreJER_etabin"
echo ""

rooUnfoldCompile smearTheorySpectra_gaussCoreJER_etabin.C

#source run_smearTheory.sh 0 

source run_smearTheory.sh 1  

#source run_smearTheory.sh 2  

#source run_smearTheory.sh 3  
