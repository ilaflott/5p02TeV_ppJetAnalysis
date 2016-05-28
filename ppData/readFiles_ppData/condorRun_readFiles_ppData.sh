#!/bin/bash

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
startfile=$1
endfile=$2
filelist=$3
outfile=$6 
radius=$4 
jetType=$5 
debug=$7
#destination=$8

echo "Processing..."
./readFiles_ppData.exe $startfile $endfile "${filelist}" "${outfile}" $radius "${jetType}" $debug

#echo "Copying file ${outfile} to ${destination}..."
#mv ${outfile} ${destination}

echo "Job done!"
