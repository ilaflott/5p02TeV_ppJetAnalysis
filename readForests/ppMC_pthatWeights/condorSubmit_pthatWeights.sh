#!/bin/bash

echo ""

## error conditions + I/O
if [[ $# -ne 1 ]] # not enough arguments
then
    echo "Usage is... "
    echo "source condorSubmit_pthatWeights.sh <minPthatBin>"

    return 1
fi

minPthatBin=$1

## cmsenv for condor
echo "cmsenv'ing on CVMFS..."
cd ${CVMFS_758}
cmsenv
cd -


# input files
exe="ppMC_pthatWeights.exe"
filelistTarBall="${PWD}/filelists.tar.gz"

# name output
outName="QCDjet${minPthatBin}_pthatWeight"
Error="${outName}.err"
Output="${outName}.out"
Log="${outName}.log"
outfile="${outName}.txt"
    
# output dir
logFileDir="${PWD}/pthatWeightsOut"

# subfile
cat > ${logFileDir}/subfile <<EOF

Universe       = vanilla
Environment = "HOSTNAME=$HOSTNAME"
Executable     = condorRun_ppMC_pthatWeights.sh
+AccountingGroup = "group_cmshi.ilaflott"
Arguments      = $minPthatBin
Input          = /dev/null
Error          = ${logFileDir}/${Error}
Output         = ${logFileDir}/${Output}
Log            = ${logFileDir}/${Log}
# get the environment (path, etc.)
GetEnv         = True
# prefer to run on fast, 64 bit computers
Rank           = kflops
Requirements   = Arch == "X86_64"
should_transfer_files   = YES
transfer_input_files = ${filelistTarBall},${exe}
when_to_transfer_output = ON_EXIT
Queue
EOF
    
    ## submit the job defined in the above submit file
echo "submitting ppMC_pthatWeights job"

condor_submit ${logFileDir}/subfile    

echo "done submitting pthatWeights code."

return