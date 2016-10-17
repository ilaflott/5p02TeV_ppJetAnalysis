#!/bin/bash

echo ""

## error conditions + I/O
if [[ $# -ne 8 ]] # not enough arguments
then
    echo "Usage is... "
    echo "source condorSubmit_ppMC.sh <readFilesScript> <NJobs> <NFilesPerJob> <startFilePos> <filelistIn> <radius> <jetType> <debug>"
    return 1
fi
#if [[ ! $3 =~ ^-?[0-9]+$ ]] # check integer input against text reg ex.

## input arguments to submit script
readFilesScript=$1 # "readForests_ppMC"
NJobs=$2
NFilesPerJob=$3
startFilePos=$4
filelistIn=$5  #echo "filelistIn is ${filelistIn}" 
radius=$6
jetType=$7
debug=$8

## additional inputs to the code, may be input to this script in the near future
readFilesScriptExe="${readFilesScript}.exe"

filelist=${filelistIn##*/} #echo "filelist is ${filelist}"
nFiles=`wc -l < $filelistIn`

## NJobs=-1 case
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

## simple error cases for startFilePos
if [[ $startFilePos -ge $nFiles ]]
then
    echo "<startFilePos> larger than filelist, exit"
    return
elif [[ $startFilePos -lt 0 ]]
then
    echo "bad <startFilePos>"
    echo "setting it to 0..."
    startFilePos=0
fi

## some debug info, just in case
NFilesRequested=$(( $NJobs * $NFilesPerJob ))
echo "the readFilesScript is ${readFilesScript}"
echo "require ${NFilesRequested} files for ${NJobs} jobs"
echo "# of files in list ${filelist}: ${nFiles}"
echo "starting at file position ${startFilePos}..."
echo ""

## grab some usefule string based on filelist name
filelistTitle=${filelist%_*} #echo "filelistTitle is ${filelistTitle}"
energy=${filelistTitle%%_*} #echo "energy is ${energy}"
trig=${filelistTitle#*_} #echo "trig is ${trig}"
outName="${trig}_ak${radius}${jetType}"

## create output directory for condor job
dirName="ppMC_${trig}_ak${radius}${jetType}Jets_$(date +"%m-%d-%y")"
logFileDir="${PWD}/outputCondor/${dirName}"
AltCounter=0
while [[ -d "${logFileDir}"  ]]
  do
  AltCounter=$(( $AltCounter + 1 ))
  echo "dir exists!"
  dirName="ppMC_${trig}_ak${radius}${jetType}Jets_$(date +"%m-%d-%y")_${AltCounter}"
  logFileDir="${PWD}/outputCondor/${dirName}"    
done
echo "output in outputCondor/${dirName}"
mkdir $logFileDir


## cmsenv for condor
echo "cmsenv'ing..."
cd ${CVMFS_758}
cmsenv
cd -


### compile code executable, same as rootcompile in my .bashrc
#echo "compiling..."
#rootcompile "${readFilesScript}.C"
echo "no compile needed, moving on...."

## copy over code used for job running/submitting for archival purposes
cp ${readFilesScript}.* "${logFileDir}"
cp readForests.h "${logFileDir}"
cp condorRun_readForests.sh "${logFileDir}"
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
Environment = "HOSTNAME=$HOSTNAME"
Executable     = condorRun_readForests.sh
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
    
    ## submit the job defined in the above submit file
    echo "running ${readFilesVersion} on files #${startfile} to #${endfile}"
    condor_submit ${logFileDir}/subfile    
    echo "sleeping..."
    sleep 3s  #my way of being nicer to condor, not sure it really matters but i'm paranoid
done

cd -
echo "done."
return