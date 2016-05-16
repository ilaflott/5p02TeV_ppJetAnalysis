#!/bin/bash

# error condition
if [[ $# -ne 4 ]]
then
    echo "Usage is... "
    echo "source condor_makeTarAndSubmit.sh <NJobs> <NFilesPerJob> <inputFileList> <debug>"
    return 1
fi

# make log folder for this submission
now="submit_$(date +"%Y-%m-%d__%H_%M_%S")"
logFileDir="${PWD}/condorLogs/${now}"
mkdir $logFileDir
echo "log files in ${logFileDir}"

# input arguments to submit script
NJobs=$1
NFilesPerJob=$2
filelist=$3
debug=$4

# cmsenv for condor
#export SCRAM_ARCH=slc6_amd64_gcc491
#source /osg/app/cmssoft/cms/cmsset_default.sh
echo "cmsenv'ing..."
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_5_8/src
cmsenv
cd -

# compile code executable
echo "compiling..."
g++ RAA_read_data_pp.C $(root-config --cflags --libs) -Werror -Wall -O2 -o RAA_read_data_pp.exe || return 1
cp RAA_read_data_pp.exe "$logFileDir"

# additional inputs to the run script and .exe, these don't change too much
destination="/mnt/hadoop/cms/store/user/ilaflott/5p02TeV_ppJetAnalysis/PP_Data/RAA_read_data_pp"
radius=4
jetType="PF"

# one condor job submit per NFilesPerJob until we submit NJobs
nFiles=`wc -l < $filelist`
echo "nFiles in list: $nFiles"


####################################
## CREATE SUBMIT FILE(S) + SUBMIT ##
####################################

NthJob=0
while [ $NthJob -lt $NJobs ]
do
    # job splitting
    startfile=$(( $NthJob * $NFilesPerJob ))
    endfile=$(( ($NthJob + 1) * $NFilesPerJob ))

    # check for end of filelist
    if [ $endfile -gt $nFiles ]; then
        let endfile=$nFiles
        let NthJob=$NJobs
    fi

    # define output names for job submission
    Error="ak$radius-ppData-${startfile}_to_${endfile}.err"
    Output="ak$radius-ppData-${startfile}_to_${endfile}.out"
    Log="ak$radius-ppData-${startfile}_to_${endfile}.log"
    outfile="5p02TeV_pp_data_ak${radius}${jetType}_20_eta_20_fileNum${startfile}to${endfile}.root"
    
    # create the condor submit file
    cat > ${logFileDir}/subfile <<EOF

Universe       = vanilla
Environment = "HOSTNAME=$HOSTNAME"
Executable     = condorRun_RAA_read_data_pp.sh
+AccountingGroup = "group_cmshi.ilaflott"
Arguments      = $startfile $endfile $filelist $radius $jetType $outfile $debug $destination
Input          = /dev/null
Error          = ${logFileDir}/$Error
Output         = ${logFileDir}/$Output
Log            = ${logFileDir}/$Log
# get the environment (path, etc.)
GetEnv         = True
# prefer to run on fast computers
Rank           = kflops
# only run on 64 bit computers
Requirements   = Arch == "X86_64"
should_transfer_files   = YES
transfer_input_files = ${filelist},RAA_read_data_pp.exe
when_to_transfer_output = ON_EXIT
# specify any extra input files (for example, an orcarc file)
Queue
EOF

    # submit the job defined in the above submit file
    NthJob=$(($NthJob + 1))
    echo "submitting RAA_read_data_pp job ${NthJob} of ${NJobs}"
    condor_submit ${logFileDir}/subfile
    
done