#!/bin/bash


# hadd's root output in dir together; for 5 TeV ppJetAnalysis/readForests/outputCondor mostly
function hadd_ppMC_allInDir(){
    if [[ $# -ne 1 ]]
    then
        echo "usage is ..."
        echo "hadd_ppMC_allInDir <radiusInt> "
        return
    fi
    echo "in dir $PWD"
    R=$1
    hadd Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PF-allFiles.root ./*.root
    return
}

function hadd_ppData_allInDir(){
    if [[ $# -ne 1 ]]
    then
        echo "usage is ..."
        echo "hadd_ppData_allInDir <radiusInt> "
        return
    fi
    R=$1
    echo "in dir $PWD"
    hadd HighPtJetTrig_ak${R}PF-allFiles.root ./ppData*/*.root
    return
}





if [ $# -ne 5 ]
then
    echo ""
    echo "Useage is..."
    echo "source renameNhadd.sh [R=3,4...] [date_output=\"MM-DD-YY\"] [date_scratch=\"MM.DD.YY\"] [etabin=\"?.?eta?.?\"] [type=\"ppType\"]"
    echo ""
    
    echo ""
    echo "e.g."
    echo "source renameNhadd.sh 4 01-01-17 1.01.17 0.0eta0.5 ppData"
    echo ""
    
    echo ""
    echo "note, date_output is the date in the name of the dir (typically on scratch) to move the output to"
    echo "note, date_scratch is the date in the name of target dir of interest containing the files to hadd"
    echo ""
    
    echo ""
    echo "TLDR; bad input. Try again. Exit."
    echo "" 
    return
fi

R=$1
date_output=$2
date_scratch=$3
etabin=$4
type=$5



#theSCRATCHdir=${SCRATCH}/5p02TeV_ppJetAnalysis/readForests/${date_scratch}_outputCondor/
#theSCRATCHdir=/cms/heavyion/ilaflott/T2_US_MIT_SCRATCH/5p02TeV_ppJetAnalysis/readForests/${date_scratch}_outputCondor/
#theSCRATCHdir=/cms/heavyion/ilaflott/T2_US_MIT_SCRATCH/5p02TeV_ppJetAnalysis/readForests/${date_scratch}_outputCondor/ppData_CaloHLT
theSCRATCHdir=/cms/heavyion/ilaflott/T2_US_MIT_SCRATCH/5p02TeV_ppJetAnalysis/readForests/${date_scratch}_outputCondor
workingDir=$PWD


if [[ -d "${theSCRATCHdir}"  ]]
then
    echo "${date_scratch} condor dir on SCRATCH exists! moving on..."
else
    echo "creating directory  in ${SCRATCH}/5p02TeV_ppJetAnalysis/readForests/"
    mkdir ${theSCRATCHdir}
fi





##########################################################################################
#### HADD N MOVE ppMC JOBS
##########################################################################################
if [ $type == "ppMC" ]
then
    echo ""
    echo "hadding ppMC job outputs for etabin=${etabin} and date=${date_output}"
    echo ""
    sleep 1s
    
    #ppMCjetPlotsDir=ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_jetPlots_${etabin}/
    ppMCjetPlotsDir=ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_JERS_${etabin}
    sleep 1s
    
    
    cd  ${ppMCjetPlotsDir}
    hadd_ppMC_allInDir ${R}
    
    echo ""
    echo "done hadding ppMC"
    echo ""

    cd ${workingDir}
    sleep 1s
        
    cp --parents $ppMCjetPlotsDir/*allFiles*.root 4scp/
    
    echo ""
    echo "moving ppMC to ${theSCRATCHdir}"
    echo ""
    mv ${ppMCjetPlotsDir} ${theSCRATCHdir}
    sleep 1s
    
    
    
    

#########################################################################################
### HADD N MOVE ppData JOBS
#########################################################################################
elif [ $type == "ppData" ]
then
    echo ""
    echo "making ppData folders"
    echo ""
    sleep 1s
    
    ppDataHPtDir=ppData_HighPtJetTrig_ak${R}PFJets_${date_output}_jetPlots_${etabin}/    
    mkdir ${ppDataHPtDir}
    sleep 1s

    echo ""
    echo "moving ppData folders"
    echo ""
    sleep 1s
    
    ppDatajet80=ppData_HighPtJet80_ak${R}PFJets_${date_output}_jetPlots_${etabin}/
    ppDatalowJets=ppData_HighPtLowerJets_ak${R}PFJets_${date_output}_jetPlots_${etabin}/
    mv ${ppDatajet80} ${ppDataHPtDir}
    mv ${ppDatalowJets} ${ppDataHPtDir}
    cd ${ppDataHPtDir}
    sleep 1s
    
    echo ""
    echo "hadding ppData"
    echo ""
    sleep 1s
    
    hadd_ppData_allInDir ${R}
    sleep 1s
    
    echo ""
    echo "done hadd'ing ppData"
    echo ""
    cd ${workingDir}
    sleep 1s
    
    cp --parents $ppDataHPtDir/*allFiles*.root 4scp/
    
    
    echo ""
    echo "moving to SCRATCH space"
    echo ""
    mv ${ppDataHPtDir} ${theSCRATCHdir}
    sleep 1s
    
fi


echo ""
echo "done"
echo ""

return 