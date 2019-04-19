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
    R=$1
    jet80dir=$2
    lowjetsdir=$3    
    
    if [[ $# -eq 4 ]]
    then
	minbiasdir=$4	
	echo "hadd'ing MinBias, Jet80, and LowerJets output"
	#return
	hadd MinBias_ak${R}PF-allFiles.root         ./${minbiasdir}*.root
	hadd HighPtJet80_ak${R}PF-allFiles.root     ./${jet80dir}*.root
	hadd HighPtLowerJets_ak${R}PF-allFiles.root ./${lowjetsdir}*.root
	hadd HighPtJetTrig_noMB_ak${R}PF-allFiles.root HighPtJet80_ak${R}PF-allFiles.root HighPtLowerJets_ak${R}PF-allFiles.root
	hadd HighPtJetTrig_ak${R}PF-allFiles.root MinBias_ak${R}PF-allFiles.root HighPtLowerJets_ak${R}PF-allFiles.root HighPtJet80_ak${R}PF-allFiles.root
    elif [[ $# -eq 3 ]] 
    then
	echo "hadd'ing Jet80 and LowerJets output"
	#return
	hadd HighPtJet80_ak${R}PF-allFiles.root     ./${jet80dir}*.root
	hadd HighPtLowerJets_ak${R}PF-allFiles.root ./${lowjetsdir}*.root
	hadd HighPtJetTrig_noMB_ak${R}PF-allFiles.root HighPtJet80_ak${R}PF-allFiles.root HighPtLowerJets_ak${R}PF-allFiles.root
    else
        echo "error. usage is ..."
        echo "hadd_ppData_allInDir <radiusInt> <jet80dir> <lowjetsdir> <[opt]minbiasdir>"
	echo "exit."
        return
    fi
    echo ""
    echo "done hadd'ing output."
    echo "in dir $PWD"

    return
}





if [ $# -ne 5 ]
then
    echo ""
    echo "Useage is..."
    #echo "source renameNhadd.sh [R=3,4...] [date_output=\"MM-DD-YY\"] [date_scratch=\"MM.DD.YY\"] [etabin=\"?.?eta?.?\"] [type=\"ppType\"]"
    echo "source renameNhadd.sh [R=3,4...] [date_output=\"MM-DD-YY\"] [date_scratch=\"MM.DD.YY\"] [desc=\"jobtype_?.?eta?.?\"] [type=\"ppType\"]"
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
#etabin=$4
#desc argument *typically* of form; jetPlots_0.0eta2.0
desc=$4
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
    echo "hadding ppMC job outputs for desc=${desc} and date=${date_output}"
    echo ""
    sleep 1s
    
    ppMCjetPlotsDir=ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_${desc}/
    #ppMCjetPlotsDir=ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_jetMult_${etabin}/
    #ppMCjetPlotsDir=ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_jetPlots_${etabin}/
    #ppMCjetPlotsDir=ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_JERS_${etabin}
    sleep 1s
    
    if [[ -d "${ppMCjetPlotsDir}" ]]
    then
	cd  ${ppMCjetPlotsDir}
	hadd_ppMC_allInDir ${R}
    else
	echo "${ppMCjetPlotsDir} not found."
	echo "exit."
	return
    fi
    
    echo ""
    echo "done hadding ppMC"
    echo ""

    cd ${workingDir}
    sleep 1s
        
    #cp --parents $ppMCjetPlotsDir/*allFiles*.root 4scp/
    
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
    
    ppDataHPtDir=ppData_HighPtJetTrig_ak${R}PFJets_${date_output}_${desc}
    mkdir ${ppDataHPtDir}
    sleep 1s
    
    echo ""
    echo "moving ppData folders"
    echo ""
    sleep 1s
    
    ppDatajet80=ppData_HighPtJet80_ak${R}PFJets_${date_output}_${desc}/
    if [[ -d "${ppDatajet80}" ]]
    then
	echo "${ppDatajet80} dir found."
	mv ${ppDatajet80} ${ppDataHPtDir}
    else
	echo "${ppDatajet80} not found."
	#echo "exit."
	#return
    fi
    ppDatalowJets=ppData_HighPtLowerJets_ak${R}PFJets_${date_output}_${desc}/
    if [[ -d "${ppDatalowJets}" ]]
    then
	echo "${ppDatalowJets} dir found."
	mv ${ppDatalowJets} ${ppDataHPtDir}
    else
	echo "${ppDatalowJets} not found."
	#echo "exit."
	#return
    fi
    ppDataMinBias=ppData_MinBias6_758p3_ak${R}PFJets_${date_output}_${desc}/
    if [[ -d "${ppDataMinBias}" ]]
    then
	echo "${ppDataMinBias} dir found."
	mv ${ppDataMinBias} ${ppDataHPtDir}
    else
	echo "${ppDataMinBias} not found."
	#echo "exit."
	#return
    fi

    cd ${ppDataHPtDir}
    sleep 1s
    
    echo ""
    echo "hadding ppData"
    echo ""
    sleep 1s

    if [[ -d "${ppDataMinBias}" && -d "${ppDatalowJets}" && -d "${ppDatajet80}" ]]
    then	
	hadd_ppData_allInDir ${R} ${ppDatajet80} ${ppDatalowJets} ${ppDataMinBias}
    elif [[ -d "${ppDatalowJets}" && -d "${ppDatajet80}" ]]
    then
	hadd_ppData_allInDir ${R} ${ppDatajet80} ${ppDatalowJets} ${ppDataMinBias}
    fi
	
    cd ${workingDir}
    sleep 1s
    
    #cp --parents $ppDataHPtDir/*allFiles*.root 4scp/    
    
    echo ""
    echo "moving to SCRATCH space"
    echo ""
    mv ${ppDataHPtDir} ${theSCRATCHdir}
    sleep 1s
    
fi


echo ""
echo "done. output in"
echo "${theSCRATCHdir}"
echo ""

return 