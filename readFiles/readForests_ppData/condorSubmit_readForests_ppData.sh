#!/bin/bash


echo ""

# error conditions + I/O
if [[ $# -ne 6 ]] # not enough arguments
then
    echo "Usage is... "
    echo "source condorSubmit_ppData.sh <inputVersion> <NJobs> <NFilesPerJob> <startFilePos> <filelistIn> <debug>"
    echo "where <readForests Ver.> = HighPtLowerJets or Vanilla"
    return 1
fi
#if [[ ! $3 =~ ^-?[0-9]+$ ]] # check integer input against text reg ex.

# input arguments to submit script
inputVersion=$1
NJobs=$2
NFilesPerJob=$3
startFilePos=$4
filelistIn=$5  #echo "filelistIn is ${filelistIn}" #debug
debug=$6

# additional inputs to the code, may be input to this script in the near future
radius=4
jetType="PF"

## checking version input,
readFilesScript=""
if [ "$inputVersion" = "Vanilla" ]
then
    echo "picked vanilla version"
    readFilesScript="readForests_ppData"
    readFilesScriptExe="${readFilesScript}.exe"
elif [ "$inputVersion" = "HighPtLowerJets" ]
then
    echo "picked HighPtLowerJets version"
    readFilesScript="readForests_ppData_${inputVersion}"
    readFilesScriptExe="${readFilesScript}.exe"
    filelistIn="filelists/5p02TeV_HighPtLowerJets_forests.txt"
else
    echo "${inputVersion} is not a readForests version. exit."
    return
fi

#elif[[ "$inputVersion" = "HighPtJet80" ]]
#    readFilesScript="readForests_ppData_${inputVersion}"
#    readFilesScriptExe="${readFilesScript}.exe"
#    filelistIn="filelists/5p02TeV_HighPtJet80_forests.txt"


# NJobs=-1 case
nFiles=`wc -l < $filelistIn`
if [[ $NJobs -eq -1 ]]
then
    echo "submitting jobs for all files in list"
    NJobs=$(( $nFiles / $NFilesPerJob ))
    if [[ $(( $nFiles % $NFilesPerJob ))  -gt 0 ]]
    then
	NJobs=$(( $NJobs + 1 ))
    fi
    startFilePos=0
fi

# simple error cases for startFilePos
if [[ $startFilePos -ge $nFiles ]]
then
    echo "<startFilePos> larger than filelist, exit"
    #echo "setting it to 0..."
    #startFilePos=0
    return
elif [[ $startFilePos -lt 0 ]]
then
    echo "bad <startFilePos>"
    echo "setting it to 0..."
    startFilePos=0
fi

# some debug info, just in case
NFilesRequested=$(( $NJobs * $NFilesPerJob ))
echo "require ${NFilesRequested} files for ${NJobs} jobs"
echo "# of files in list: ${nFiles}"
echo "starting at file position ${startFilePos}..."
echo ""


# create output folder/logfileNames with name based on filelist
filelist=${filelistIn##*/} #echo "filelist is ${filelist}"
filelistTitle=${filelist%_*} #echo "filelistTitle is ${filelistTitle}"
energy=${filelistTitle%%_*} #echo "energy is ${energy}"
trig=${filelistTitle#*_} #echo "trig is ${trig}"

#dirName="readForests_ppData_${energy}_${trig}_$(date +"%Y-%m-%d__%H_%M")"
dirName="${readFilesScript}_${energy}_${trig}_$(date +"%m-%d-%y__%H_%M")"
outName="${trig}_ak${radius}${jetType}"

logFileDir="${PWD}/outputCondor/${dirName}"
if [ -d "${logFileDir}" ]; then
    rm -rf "${logFileDir}"
fi
mkdir $logFileDir
echo "log files in outputCondor/${dirName}"


# cmsenv for condor
echo "cmsenv'ing..."
#cd /cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_5_8/src
cd ${CVMFS_758}
cmsenv
cd -


# compile code executable, same as rootcompile in my .bashrc
echo "compiling..."
#g++ readForests_ppData.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readForests_ppData.exe || return 
rootcompile "${readFilesScript}.C"


# copy over code used for job running/submitting for archival purposes
cp ${readFilesScript}.* "${logFileDir}"
cp condorRun_readForests_ppData.sh "${logFileDir}"
cp ${filelistIn} "${logFileDir}"
cd ${logFileDir}


### CREATE NAMES AND FILES, THEN SUBMIT ###
NthJob=0
filelistLength=$nFiles
lastFilePos=$(( $filelistLength-1 ))
startfile=0
endfile=0
while [ $NthJob -lt $NJobs ]
do 
    JobNum=$(( $NthJob + 1 ))
    echo ""
    echo "SPLITTING FILES FOR JOB # ${JobNum} of ${NJobs}"

    # start/end file 
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
    # check; end of filelist
    if [[ $endfile -ge $lastFilePos ]] 
    then 
	echo "end of filelist!"
	let endfile=$lastFilePos
	let NthJob=$(( $NJobs - 1 )) 
    fi
    #echo "endfile is ${endfile}"     

    
    # for next job
    NthJob=$(($NthJob + 1))

    # define output names for job submission
    fileRange="${startfile}to${endfile}"
    Error="${outName}-${fileRange}.err"
    Output="${outName}-${fileRange}.out"
    Log="${outName}-${fileRange}.log"
    outfile="${outName}-${fileRange}.root"
    
    # create the condor submit file
    cat > ${logFileDir}/subfile <<EOF

Universe       = vanilla
Environment = "HOSTNAME=$HOSTNAME"
Executable     = condorRun_readForests_ppData.sh
+AccountingGroup = "group_cmshi.ilaflott"
Arguments      = $readFilesScriptExe $startfile $endfile $filelist $outfile $radius $jetType $debug
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
transfer_input_files = ${filelist},${readFilesScriptExe}
when_to_transfer_output = ON_EXIT
Queue
EOF
    
    # submit the job defined in the above submit file
    echo "running ${readFilesVersion} on files #${startfile} to #${endfile}"
    condor_submit ${logFileDir}/subfile    
    sleep 0.2s #my way of being nicer to condor, not sure it really matters but i'm paranoid
done

cd -
echo "done."
return