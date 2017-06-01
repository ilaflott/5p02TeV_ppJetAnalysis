#!/bin/bash

if [ $# -eq 0 ]
then
    echo "give me a radius (3,4... etc.)"
    return
fi

R=$1





echo ""
echo "compiling readForests_ppMC_JERS"
echo ""
rootcompile readForests_ppMC_JERS.C

echo ""
echo "submitting job(s)"
echo ""

## -------------------------------------------- #
#### OFFICIAL MC for jetRECO/analysis meeting ### NOTE: NO JET ID, NO UNFOLDING HISTS, GEN/RECO JET PT > 56
## -------------------------------------------- #

source condorSubmit_readForests.sh readForests_ppMC_JERS -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0 0.0 5.0
#source condorSubmit_readForests.sh readForests_ppMC_JERS -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt ${R} PF 0 0.0 5.0
sleep 10s									       


## ---------------------------------------- #
#### SEMIOFFICIAL MC for analysis meeting ### NOTE: wJET ID, WITH UNFOLDING HISTS, GEN/RECO JET PT > 56
## ---------------------------------------- #
#source condorSubmit_readForests.sh readForests_ppMC_JERS -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0  0.0  0.5
#sleep 1s 															      	      
#source condorSubmit_readForests.sh readForests_ppMC_JERS -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0  0.5  1.0
#sleep 1s 															      	      
#source condorSubmit_readForests.sh readForests_ppMC_JERS -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0  1.0  1.5
#sleep 1s 															      	      
#source condorSubmit_readForests.sh readForests_ppMC_JERS -1 20 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0  1.5  2.0
#sleep 1s 

echo ""
echo "asking condor"
echo ""

source askCondor.sh 20 20

echo ""
echo "done"
echo ""

return
