#!/bin/bash

echo ""

# error conditions 
if [[ $# -ne 1 ]] # not enough arguments
then
    echo "Usage is... "
    echo "source condorSubmit_generic.sh someScript.C"
    return 1
fi

# input arguments to submit script
theCode=$1
echo "theCode is = ${theCode}"

theCodeNoC=${theCode/%.C/}
echo "theCodeNoC is = ${theCodeNoC}"

theExe="${theCodeNoC}.exe"
echo "theExe is = ${theExe}"

dirName="${theCodeNoC}_$(date +"%m-%d-%y__%H_%M")"
echo "dirName is = ${dirName}"
echo "output will be in outputCondor/${dirName}, creating directory..."

logFileDir="${PWD}/outputCondor/${dirName}"
if [ -d "${logFileDir}" ]; then
    rm -rf "${logFileDir}"
fi
mkdir $logFileDir

# cmsenv for condor
echo "cmsenv'ing..."
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_5_8/src
cmsenv
cd -

echo "compiling..."
rootcompile ${theCode}
#echo "compiling..."

cp "${theCodeNoC}".* "${logFileDir}"
cp condorRun_generic.sh "${logFileDir}"

cd ${logFileDir}

### CREATE NAMES AND FILES, THEN SUBMIT ###
Error="${theCodeNoC}_condorRun_generic.err"
Output="${theCodeNoC}_condorRun_generic.out"
Log="${theCodeNoC}_condorRun_generic.log"
#outfile="${theCodeNoC}_condorRun_generic.root"

# create the condor submit file
cat > ${logFileDir}/subfile <<EOF

Universe       = vanilla
Environment = "HOSTNAME=$HOSTNAME"
Executable     = condorRun_generic.sh
+AccountingGroup = "group_cmshi.ilaflott"
Arguments      = $theExe
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
transfer_input_files = $theExe
when_to_transfer_output = ON_EXIT
Queue
EOF

# submit the job defined in the above submit file
echo "running ${theCode} on condor, submitting..."
condor_submit ${logFileDir}/subfile    

cd -
echo "done."
return