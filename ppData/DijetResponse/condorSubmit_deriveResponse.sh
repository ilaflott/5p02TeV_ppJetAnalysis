#!/bin/bash

# error conditions 
if [[ $# -ne 5 ]] # not enough arguments
then
    echo "Usage is... "
    echo "source condorSubmit_deriveResponse.sh <NJobs> <NFilesPerJob> <startFilePos> <filelistIn> <debug>"
    return 1
elif [[ ! $3 =~ ^-?[0-9]+$ ]] # check integer input
then
    echo "non integer <startFilePos>, exit"
    return 1
elif [[ $3 -lt 0 ]] # check for valid startFilePos
then
    echo "bad <startFilePos>, exit"
    echo "0 <= <startFilePos> <= nFiles-1"
    return 1
fi

# input arguments to submit script
NJobs=$1
NFilesPerJob=$2
startFilePos=$3
filelistIn=$4  #echo "filelistIn is ${filelistIn}" #debug
debug=$5

# one condor job submit per NFilesPerJob until we submit NJobs
nFiles=`wc -l < $filelistIn`
if [[ $startFilePos -ge $nFiles ]]
then
    echo "<startFilePos> too big, exit"
    echo "0 <= <startFilePos> < ${nFiles}"
    return 1
fi
echo ""
echo "# of files in list: ${nFiles}"
echo "you requested ${NJobs} jobs with ${NFilesPerJob} files per job"
NFilesRequested=$(( $NJobs * $NFilesPerJob ))
echo "${NFilesRequested} for jobs to run"
echo "starting at file position ${startFilePos}..."
echo ""

# additional inputs to the run script and .exe, these don't change too much
destination="/mnt/hadoop/cms/store/user/ilaflott/5p02TeV_ppJetAnalysis/ppData/DijetResponse"
radius=4
isMC=0
jetType="PF"

# create output folder/logfileNames with name based on filelist
filelist=${filelistIn##*/} #echo "filelist is ${filelist}"
filelistTitle=${filelist%_*} #echo "filelistTitle is ${filelistTitle}"
energy=${filelistTitle%%_*} #echo "energy is ${energy}"
trig=${filelistTitle#*_} #echo "trig is ${trig}"
dirName="deriveResponse_${energy}_${trig}_$(date +"%Y-%m-%d__%H_%M")"
outName="${trig}_ak${radius}${jetType}"
logFileDir="${PWD}/outputCondor/${dirName}"
mkdir $logFileDir
echo "log files in outputCondor/${dirName}"

# cmsenv for condor
echo "cmsenv'ing..."
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_5_8/src
cmsenv
cd -

# compile code executable, same as rootcompile in my .bashrc
echo "compiling..."
g++ DijetResponse.C $(root-config --cflags --libs) -Werror -Wall -O2 -o DijetResponse.exe || return 1
cp DijetResponse.* "${logFileDir}"

### CREATE NAMES AND FILES, THEN SUBMIT ###

NthJob=0
filelistLength=$nFiles
lastFilePos=$(( $filelistLength-1 ))
startfile=0
endfile=0
#echo "startfile is ${startfile}" 
#echo "endfile is ${endfile}"     
while [ $NthJob -lt $NJobs ]
do 
    echo ""
    echo "SPLITTING FILES FOR JOB #${NthJob}"

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
    echo "startfile is ${startfile}" 
    # check; end of filelist
    if [[ $endfile -ge $lastFilePos ]] 
    then 
	echo "end of filelist!"
	let endfile=$lastFilePos
	let NthJob=$(( $NJobs - 1 )) 
    fi
    echo "endfile is ${endfile}"     

    
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
Executable     = condorRun_deriveResponse.sh
+AccountingGroup = "group_cmshi.ilaflott"
Arguments      = $startfile $endfile $filelist $outfile $radius $isMC $debug
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
transfer_input_files = ${filelistIn},DijetResponse.exe
when_to_transfer_output = ON_EXIT
Queue
EOF

    # submit the job defined in the above submit file
    echo "running files #${startfile} to #${endfile}"
    echo "submitting deriveResponse job..."
    condor_submit ${logFileDir}/subfile    
    sleep 0.5s #my way of being nicer to condor, not sure it really matters but i'm paranoid
done