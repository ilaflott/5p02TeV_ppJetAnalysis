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
outfile=$4
radius=$5 
jetType=$6 
debug=$7

echo "Processing..."
./readForests_ppData.exe $startfile $endfile "${filelist}" "${outfile}" $radius "${jetType}" $debug

echo "Job done!"
