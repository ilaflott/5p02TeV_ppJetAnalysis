#!/bin/bash


##CONST
date_data="05-25-17"
date_condor="5.25.17"
outCondorDir="${date_condor}_outputCondor"
##CONST


if [ $# -ne 5 ]
then
    echo "bad usage."
    echo "source run_unfoldSpectra <R> <jtID> <compile> <SVDkReg> <etaBin#>"
    echo ""

    echo "e.g."
    echo "source run_unfoldSpectra 3 1 1 20 1"
    echo ""

    return
fi

R=$1
jtID=$2
compile=$3
kReg=$4
etaBin=$5

if [ $compile -eq 1 ]
then
    echo ""
    echo "compiling unfoldDataSpectra"
    echo ""
    rooUnfoldCompile unfoldDataSpectra.C

    echo ""
    echo "done compiling unfoldDataSpectra"
    echo ""
    sleep 0.5s
else
    echo ""
    echo "skipping compilation..."
    echo ""
fi




#### unfoldDataSpectra
echo ""
echo "running unfoldDataSpectra"
echo ""

if [ $etaBin -eq 1 ]
then
    ./unfoldDataSpectra.exe ${outCondorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${date_data}_jetPlots_0.0eta0.5  ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_0.0eta0.5     HPtJetTrig_Py8_unfData_00eta05   ${jtID} ${kReg}
elif [ $etaBin -eq 2 ]
then
    ./unfoldDataSpectra.exe ${outCondorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${date_data}_jetPlots_0.5eta1.0  ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_0.5eta1.0     HPtJetTrig_Py8_unfData_05eta10   ${jtID} ${kReg}
else
    echo ""
    echo "etaBin not found; exit"
    echo ""
    return
fi

#./unfoldDataSpectra.exe ${outCondorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${date_data}_jetPlots_1.0eta1.5  ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_1.0eta1.5     HPtJetTrig_Py8_unfData_10eta15  ${jtID}
#./unfoldDataSpectra.exe ${outCondorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${date_data}_jetPlots_1.5eta2.0  ${outCondorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_data}_JERS_1.5eta2.0     HPtJetTrig_Py8_unfData_15eta20  ${jtID}

echo ""
echo "unfoldDataSpectra done."
echo ""

return