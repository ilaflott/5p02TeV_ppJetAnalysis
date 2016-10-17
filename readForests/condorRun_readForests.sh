#!/bin/bash

if [[ $# -ne 8 ]] # not enough arguments
then
    echo "not enough input arguments to condorRun script, exit"
    return 1
fi

# env info on condor CPU
echo ""
echo "----------------------------------------------------"
echo "Job started on `date` at WN: `hostname` "
echo "Job is running on `uname -a`"
echo "----------------------------------------------------"
echo ""
echo "free disk space is..."
df -h
echo ""
echo "root directory is..."
echo $ROOTSYS
echo ""
echo "working directory is..."
pwd
echo ""
echo "contents before job processing is..."
ls -al
echo ""

# inputs for run+root script
readFilesExe=$1
startfile=$2
endfile=$3
filelist=$4
outfile=$5
radius=$6
jetType=$7 
debug=$8

echo "Processing.."

./${readFilesExe} "${filelist}" ${startfile} ${endfile} ${radius} "${jetType}" ${debug} "${outfile}"

echo "Job done!"
