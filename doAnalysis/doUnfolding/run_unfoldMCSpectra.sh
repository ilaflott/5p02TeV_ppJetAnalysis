#!/bin/bash


##CONST
#data_date="05-25-17"
#data_date_condor="5.25.17"
#data_outCondorDir="${data_date_condor}_outputCondor"

MC_date="06-05-17"
MC_date_condor="6.05.17"
MC_outCondorDir="${MC_date_condor}_outputCondor"
##CONST


if [ $# -ne 4 ]
then
    echo "bad usage."
    echo "source run_unfoldMCSpectra.sh <R> <jtID> <SVDkReg> <etaBin#>"
    echo ""

    echo "e.g."
    echo "source run_unfoldMCSpectra.sh 3 1 20 1"
    echo ""

    return
fi

R=$1
jtID=$2
#compile=$3
kReg=$3
etaBin=$4

#if [ $compile -eq 1 ]
#then
#    echo ""
#    echo "compiling unfoldMCSpectra"
#    echo ""
#    rooUnfoldCompile unfoldMCSpectra.C
#
#    echo ""
#    echo "done compiling unfoldMCSpectra"
#    echo ""
#    sleep 0.5s
#else
#    echo ""
#    echo "skipping compilation..."
#    echo ""
#fi

#### unfoldMCSpectra
echo ""
echo "running unfoldMCSpectra"
echo ""

if [ $etaBin -eq 1 ]
then
    ./unfoldMCSpectra.exe   ${MC_outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${MC_date}_JERS_0.0eta0.5   Py8_closureTest_00eta05  ${jtID}  ${kReg}
elif [ $etaBin -eq 2 ]
then
    ./unfoldMCSpectra.exe   ${MC_outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${MC_date}_JERS_0.5eta1.0   Py8_closureTest_05eta10  ${jtID}  ${kReg}
else
    echo ""
    echo "etaBin not found; exit"
    echo ""
    return
fi

echo ""
echo "unfoldMCSpectra done!"
echo ""

return