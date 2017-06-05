#!/bin/bash

echo ""
echo "compiling readForests_ppMC_JERS.C"
echo ""

rootcompile readForests_ppMC_JERS.C

echo ""
echo "done compiling."
echo ""

echo ""
echo "running R=3"
echo ""
source run_readForests_ppMC_JERS.sh 3

source askCondor.sh 12 5

echo ""
echo "running R=4"
echo ""
source run_readForests_ppMC_JERS.sh 4

source askCondor.sh 12 5

echo ""
echo "running R=5"
echo ""
source run_readForests_ppMC_JERS.sh 5

source askCondor.sh 12 5

echo ""
echo "done submitting. Condor says..."
echo ""

source askCondor.sh 20 30

echo ""
echo "done!"
echo ""


return
