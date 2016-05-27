#!/bin/bash

# error condition
if [[ $# -ne 4 ]]
then
    echo "Usage is... "
    echo "source condor_makeTarAndSubmit.sh <NJobs> <NFilesPerJob> <filelistIn> <debug>"
    return 1
fi

# input arguments to submit script
NJobs=$1
NFilesPerJob=$2
filelistIn=$3  #echo "filelistIn is ${filelistIn}" #debug
debug=$4

# additional inputs to the run script and .exe, these don't change too much
#destination="/mnt/hadoop/cms/store/user/ilaflott/5p02TeV_ppJetAnalysis/PP_Data/DijetResponse"
radius=4
isMC=0

# create output folder/logfileNames with name based on filelist
filelist=${filelistIn##*/} #../filelists/5p02TeV_trig_forests.txt -> 5p02TeV_trig_forests.txt
filelistTitle=${filelist%_*} #5p02TeV_trig_forests.txt -> 5p02TeV_trig
energy=${listSubStr%_*}
trig=${listSubStr#*_}
dirName="deriveResponse_ppData_${energy}_${trig}_$(date +"%Y-%m-%d__%H_%M")"
outName="${trig}_ak${radius}${jetType}"
logFileDir="${PWD}/outputCondor/${dirName}"
mkdir $logFileDir
echo "log files in ${logFileDir}"

# cmsenv for condor
echo "cmsenv'ing..."
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_5_8/src
cmsenv
cd -

# compile code executable
echo "compiling..."
g++ DijetResponse.C $(root-config --cflags --libs) -Werror -Wall -O2 -o DijetResponse.exe || return 1
cp DijetResponse.* "${logFileDir}"

# one condor job submit per NFilesPerJob until we submit NJobs
nFiles=`wc -l < $filelistIn`
echo "nFiles in list: ${nFiles}"

####################################
## CREATE SUBFILE(S) + SUBMIT ##
####################################

NthJob=0
while [ $NthJob -lt $NJobs ]
do
    # job splitting
    startfile=$(( $NthJob * $NFilesPerJob ))
    endfile=$(( ($NthJob + 1) * $NFilesPerJob ))

    # check for end of file list
    if [ $endfile -gt $nFiles ]; then
        let endfile=$nFiles #last file for this submission is the last one in the list
        let NthJob=$(($NJobs - 1)) #last job submission for this loop
    fi
    
    # define output names for job submission
    fileRange="${startfile}to${endfile}"
    Error="${outName}-${fileRange}.err"
    Output="${outName}-${fileRange}.out"
    Log="${outName}-${fileRange}.log"
    outfile="${outName}-${fileRange}.root"

    # create the condor submit file
    #Arguments      = $startfile $endfile $filelist $radius $isMC $outfile $debug $destination
    cat > ${logFileDir}/subfile <<EOF

Universe         = vanilla
Environment      = "HOSTNAME=$HOSTNAME"
Executable       = condorRun_derive_DijetResponse.sh
+AccountingGroup = "group_cmshi.ilaflott"
Arguments        = $startfile $endfile $filelist $outfile $radius $isMC $debug
Input            = /dev/null
Error            = ${logFileDir}/$Error
Output           = ${logFileDir}/$Output
Log              = ${logFileDir}/$Log
# get the environment (path, etc.)
GetEnv           = True
# prefer to run on fast, 64 bit computers
Rank             = kflops
Requirements     = Arch == "X86_64"
should_transfer_files   = YES
transfer_input_files    = ${filelistIn},DijetResponse.exe
when_to_transfer_output = ON_EXIT
Queue
EOF

    # submit the job defined in the above submit file
    NthJob=$(($NthJob + 1))
    echo "submitting derive DijetResponse job ${NthJob} of ${NJobs}"
    condor_submit ${logFileDir}/subfile    
    sleep 1s
done