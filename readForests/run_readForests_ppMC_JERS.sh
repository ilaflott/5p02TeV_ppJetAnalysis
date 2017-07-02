#!/bin/bash

if [ $# -eq 0 ]
then
    echo "give me a radius (3,4... etc.)"
    return
fi

R=$1

### semiOfficial Py8

echo ""
echo "submitting Py8 ppMC JERS job(s)"
echo ""

source condorSubmit_readForests.sh readForests_ppMC_JERS -1 15 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt ${R} PF 0  0.0  0.5

echo ""
echo "done submitting JERS job(s)"
echo ""

sleep 0.5s

return
