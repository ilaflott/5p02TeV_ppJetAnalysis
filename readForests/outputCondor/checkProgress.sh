#!/bin/bash

function checkProgress
{
    dir=$1
    etabin=$2
    jobsCompleted=`ls -lt ${dir}/*.root | wc -l`
    totalJobs=`ls -lt ${dir}/*.err | wc -l`
    errLines=`cat ${dir}/*.err | wc -l`
    percentDone=$jobsCompleted
    percentDone=$(( $percentDone * 100 ))
    percentDone=$(( $percentDone / $totalJobs ))
    echo "dir: ${dir}"
    echo "====================================="
    echo " ${percentDone}% complete, ${jobsCompleted}/${totalJobs} jobs done"
    if [[ $errLines -gt 0 ]]
    then
	echo "WARNING!"
	echo "errors int jobs detected!"
	echo "inspect me further!"
    else
	echo " no errors!"
    fi
    echo ""
    echo ""
    return
}

Nasks=$1
T=$2
Ncount=0
date="10-02-17"

#note: apparently ls'ing into an array is not a favored approach, generally. 
#if funny behavior, suspect this.
dirarray=($(ls -d pp*/))

while [[ $Ncount -lt $Nasks ]]
do
    
    echo ""
    echo "checking job progress!"
    echo ""
    
    for i in "${dirarray[@]}"
    do
	checkProgress ${i} eta
    done
    




    #echo "_______________________"
    #echo ""
    #echo "bin: 0.0 < |eta| < 2.0"
    #echo ""
    #checkProgress ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_${date}_jetPlots_0.0eta2.0 0.0eta2.0
    #checkProgress                ppData_HighPtJet80_ak4PFJets_${date}_jetPlots_0.0eta2.0 0.0eta2.0
    #checkProgress            ppData_HighPtLowerJets_ak4PFJets_${date}_jetPlots_0.0eta2.0 0.0eta2.0
    #
    #
    #echo "_______________________"
    #echo ""
    #echo "bin: 0.0 < |eta| < 0.5"
    #echo ""
    #checkProgress ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_${date}_jetPlots_0.0eta0.5 0.0eta0.5
    #checkProgress                ppData_HighPtJet80_ak4PFJets_${date}_jetPlots_0.0eta0.5 0.0eta0.5
    #checkProgress            ppData_HighPtLowerJets_ak4PFJets_${date}_jetPlots_0.0eta0.5 0.0eta0.5
    #
    #
    #echo "_______________________"
    #echo ""
    #echo "bin: 0.5 < |eta| < 1.0"
    #echo ""    
    #checkProgress ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_${date}_jetPlots_0.5eta1.0 0.5eta1.0
    #checkProgress                ppData_HighPtJet80_ak4PFJets_${date}_jetPlots_0.5eta1.0 0.5eta1.0
    #checkProgress            ppData_HighPtLowerJets_ak4PFJets_${date}_jetPlots_0.5eta1.0 0.5eta1.0
    #
    #
    #echo "_______________________"
    #echo ""
    #echo "bin: 1.0 < |eta| < 1.5"
    #echo ""
    #checkProgress ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_${date}_jetPlots_1.0eta1.5 1.0eta1.5
    #checkProgress                ppData_HighPtJet80_ak4PFJets_${date}_jetPlots_1.0eta1.5 1.0eta1.5
    #checkProgress            ppData_HighPtLowerJets_ak4PFJets_${date}_jetPlots_1.0eta1.5 1.0eta1.5
    #
    #echo "_______________________"
    #echo ""
    #echo "bin: 1.5 < |eta| < 2.0"
    #echo ""    
    #checkProgress ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_${date}_jetPlots_1.5eta2.0 1.5eta2.0
    #checkProgress                ppData_HighPtJet80_ak4PFJets_${date}_jetPlots_1.5eta2.0 1.5eta2.0
    #checkProgress            ppData_HighPtLowerJets_ak4PFJets_${date}_jetPlots_1.5eta2.0 1.5eta2.0

    



    sleep ${T}s
    Ncount=$(( $Ncount + 1 ))
done

echo ""
echo "done."
echo ""
