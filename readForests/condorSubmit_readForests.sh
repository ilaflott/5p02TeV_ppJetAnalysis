#!/bin/bash
##const##
usetargdir=false
##const##
echo ""

## error conditions + I/O
if [[ $# -eq 11 ]]
then
    echo "target directory specified."
    usetargdir=true
elif [[ $# -ne 10 ]] # not enough arguments
then
    echo "Usage is... "
    echo "source condorSubmit_readForests.sh <readForestsCode> <NJobs> <NFilesPerJob> <startFilePos> <filelistIn> <radius> <jetType> <debug> [<etaCutLo> <etaCutHi>]"
    echo "to run over all files with <NFilesPerJob>, set <NJobs> to -1"
    return 1
fi


## input arguments to submit script
readForestsCode=$1 
NJobs=$2
NFilesPerJob=$3
filelistIn=$5  
startFilePos=$4
radius=$6
jetType=$7
debug=$8
etaCutLo=$9
etaCutHi=${10}
targdir=${11}


echo "etaCutLo=$etaCutLo"
echo "etaCutHi=$etaCutHi"

## NJobs=-1 case
nFiles=`wc -l < $filelistIn`
if [[ $NJobs -eq -1 ]]
then    
    NJobs=$(( $nFiles / $NFilesPerJob ))
    if [[ $(( $nFiles % $NFilesPerJob ))  -gt 0 ]]
    then
	NJobs=$(( $NJobs + 1 ))
    fi
    startFilePos=0
fi


## some debug info, just in case
NFilesRequested=$(( $NJobs * $NFilesPerJob ))
echo "running ${readForestsCode}"
echo "requesting ${NFilesRequested} files for ${NJobs} jobs" #, starting at ${startFilePos}"
echo ""
echo "number of files in ${filelistIn}..."
echo "=${nFiles}"
echo ""


## grab strings from filelistIn, readForestsCode input
filelist=${filelistIn##*/} #gets rid of "filelists/ on the left of filelistIn"

filelistTitle=${filelist%_*} #get rid of stuff to the right of the filelist title _forests.txt
energy=${filelistTitle%%_*} #get rid of everything to right of 5p02Tev
trig=${filelistTitle#*_} #get rid of everything to left of trigName

dataTypeAndVer=${readForestsCode#*_} #get rid of stuff to left, basically readForests_ 
dataType=${dataTypeAndVer%%_*} # get rid of _ver from ppDataMC_ver

readForestsVer=${dataTypeAndVer#*_} # get rid of ppDataMC_ from ppDataMC_ver
if [[ "$dataType" == "$readForestsVer" ]]
then
    readForestsVer=""
else
    readForestsVer="_${readForestsVer}"
fi


## make strings from input 
outName="${dataType}_${trig}_ak${radius}${jetType}Jets" #echo "outName is ${outName}" #debug
dirName="${outName}_$(date +"%m-%d-%y")${readForestsVer}_${etaCutLo}eta${etaCutHi}"
logFileDir="${PWD}/outputCondor/${dirName}"

## create output directory for condor job

AltCounter=0
while [[ -d "${logFileDir}"  ]]
do
    AltCounter=$(( $AltCounter + 1 ))
    echo "dir exists!"
    dirName="${outName}_$(date +"%m-%d-%y")${readForestsVer}_${etaCutLo}eta${etaCutHi}_${AltCounter}"
    logFileDir="${PWD}/outputCondor/${dirName}"    
done

# if targdir specified exists and there's output in it already, condor will overwrite it with the new output
if [[ "$usetargdir" = true ]]
then
    echo "target directory specified. "
    dirName=${targdir}
    logFileDir="${PWD}/outputCondor/${dirName}"
fi

echo "output in outputCondor/${dirName}"
#return ## DEBUG

mkdir $logFileDir

if [[ -d "${logFileDir}"  ]]
then
    echo "logFileDir exists."
else
    echo "logFileDir not created. exit."
    return
fi

### uncomment me to work on T2 US MIT
### cmsenv for condor
#echo "cmsenv'ing on CVMFS..."
#cd ${CVMFS_758}
#cmsenv
#cd -
### /uncomment me


## copy over code used for job running/submitting for archival purposes
#rootcompile ${readForestsCode}.C
cp ${readForestsCode}.* "${logFileDir}"
cp readForests_*.h "${logFileDir}"
#cp L2L3ResidualWFits.h "${logFileDir}"
#cp MCTruthResidual.h "${logFileDir}"
cp readForests_JEC_v12_ak4PF_74X_dataRun2_HLT_ppAt5TeV_v0.h "${logFileDir}"
cp readForests_JEC_v12_ak4PF_75X_data_Run2_v13.h "${logFileDir}"
cp condorRun_readForests.sh "${logFileDir}"
cp ${filelistIn} "${logFileDir}"

readForestsExe="${readForestsCode}.exe"

## if ppMC, no JEC tar ball needed (takes up to 2 MB each time)
readForestsCode2=${readForestsCode#*_}
#echo "readForestsCode2 = $readForestsCode2"
readForestsCode3=${readForestsCode2%_*}
#echo "readForestsCode3 = $readForestsCode3"
#if [[ "$readForestsCode3" == "ppData" ]]
#then
#    transferInputFileList="${filelist},${readForestsExe},JECDataDriven.tar.gz,residualMCTruth.tar.gz"
#    cp JECDataDriven.tar.gz "${logFileDir}"   
#    cp residualMCTruth.tar.gz "${logFileDir}"   
#else
transferInputFileList="${filelist},${readForestsExe}"
#fi
#echo "transferInputFileList=$transferInputFileList"
cd ${logFileDir}


### CREATE NAMES AND FILES, THEN SUBMIT ###
lastFilePos=$(( $nFiles-1 ))
NthJob=0
startfile=0
endfile=0
while [ $NthJob -lt $NJobs ]
do 
    JobNum=$(( $NthJob + 1 ))
    echo ""
    echo "SPLITTING FILES FOR JOB # ${JobNum} of ${NJobs}"

    ## start/end file 
    if [[ $NthJob -le 0 ]]
    then
	startfile=$startFilePos
	endfile=$(( $startfile + $NFilesPerJob ))
	endfile=$(( $endfile - 1 ))
    else 
	startfile=$(( $endfile + 1 ))
	endfile=$(( $startfile + $NFilesPerJob ))
	endfile=$(( $endfile - 1 ))
    fi
    #echo "startfile is ${startfile}" 

    ## check; end of filelist
    if [[ $endfile -ge $lastFilePos ]] 
    then 
	echo "end of filelist!"
	let endfile=$lastFilePos
	let NthJob=$(( $NJobs - 1 )) 
    fi
    #echo "endfile is ${endfile}"     

    
    ## for next job
    NthJob=$(($NthJob + 1))

    ## define output names for job submission
    fileRange="${startfile}to${endfile}"
    Error="${outName}-${fileRange}.err"
    Output="${outName}-${fileRange}.out"
    Log="${outName}-${fileRange}.log"
    outfile="${outName}-${fileRange}.root"
    
    ## create the condor submit file
    cat > ${logFileDir}/subfile <<EOF
Universe       = vanilla
Executable     = condorRun_readForests.sh
Arguments      = $readForestsExe $startfile $endfile $filelist $outfile $radius $jetType $debug $etaCutLo $etaCutHi
Input          = /dev/null
Error          = ${logFileDir}/$Error
Output         = ${logFileDir}/$Output
Log            = ${logFileDir}/$Log
# get the environment (path, etc.)
GetEnv         = True
# prefer to run on fast, 64 bit computers
Rank           = kflops
Requirements   = Arch == "X86_64"
should_transfer_files   = YES
transfer_input_files = ${transferInputFileList}
when_to_transfer_output = ON_EXIT
Notification  =  never
Queue
EOF
    
    ## submit the job defined in the above submit file
    echo "running ${readForestsCode} on files #${startfile} to #${endfile}"
    condor_submit ${logFileDir}/subfile    
    #sleep 1.0s  #my way of being nicer to condor, not sure it really matters but i'm paranoid
done

cd -
echo "done."
condor_q $USER
return
















#### MIT T2 SUBMIT FILE


#### Universe       = vanilla
#### Environment = "HOSTNAME=$HOSTNAME"
#### Executable     = condorRun_readForests.sh
#### +AccountingGroup = "group_cmshi.ilaflott"
#### Arguments      = $readForestsExe $startfile $endfile $filelist $outfile $radius $jetType $debug $etaCutLo $etaCutHi
#### Input          = /dev/null
#### Error          = ${logFileDir}/$Error
#### Output         = ${logFileDir}/$Output
#### Log            = ${logFileDir}/$Log
#### # get the environment (path, etc.)
#### GetEnv         = True
#### # prefer to run on fast, 64 bit computers
#### Rank           = kflops
#### Requirements   = Arch == "X86_64"
#### should_transfer_files   = YES
#### transfer_input_files = ${filelist},${readForestsExe},JECDataDriven.tar.gz
#### when_to_transfer_output = ON_EXIT
#### match_list_length = 5
#### Queue















### grab strings from filelistIn, readForestsCode input
#filelist=${filelistIn##*/} #gets rid of "filelists/ on the left of filelistIn"
##echo "filelist is ${filelist}" #debug
##echo ""
#
#filelistTitle=${filelist%_*} #get rid of stuff to the right of the filelist title _forests.txt
##echo "filelistTitle is ${filelistTitle}" #debug
##echo ""
#
#energy=${filelistTitle%%_*} #get rid of everything to right of 5p02Tev
##echo "energy is ${energy}" #debug
##echo ""
#
#trig=${filelistTitle#*_} #get rid of everything to left of trigName
##echo "trig is ${trig}" #debug
##echo ""
##echo ""
#
#
#dataTypeAndVer=${readForestsCode#*_} #get rid of stuff to left, basically readForests_ 
##echo "dataTypeAndVer is ${dataTypeAndVer}" #debug
##echo ""
#
#dataType=${dataTypeAndVer%_*} # get rid of _ver from ppDataMC_ver
##echo "dataType is ${dataType}" #debug
##echo ""
#
#readForestsVer=${dataTypeAndVer#*_} # get rid of ppDataMC_ from ppDataMC_ver
#echo "readForestsVer is ${readForestsVer}" #debug
##echo ""
##echo ""