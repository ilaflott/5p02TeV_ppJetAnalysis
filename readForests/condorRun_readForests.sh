#!/bin/bash

if [[ $# -ne 10 ]] # not enough arguments
then
    echo "not enough input arguments to condorRun script, exit"
    return 1
fi

## env info on condor CPU
echo ""
echo "----------------------------------------------------"
echo "Job started on `date` at WN: `hostname` "
echo "Job is running on `uname -a`"
echo "----------------------------------------------------"
echo ""
#echo "free disk space is..."
#df -h
#echo ""
echo "root directory is..."
echo $ROOTSYS
echo ""
echo "working directory is..."
pwd
echo ""
echo "contents before job processing is..."
ls -al
echo ""

#if [ ! -f ./JECDataDriven.tar.gz ]
#then
#    sleep 0.01s
#    #echo "JEC file not found. Residual corrections will not be applied."
#else
#    #echo "JEC file found. Residual corrections may be applied."
#    #tar -zxvf JECDataDriven.tar.gz
#    tar -zxf JECDataDriven.tar.gz
#fi
#
#if [ ! -f ./residualMCTruth.tar.gz ]
#then
#    sleep 0.01s
#    #echo "JEC file not found. Residual corrections will not be applied."
#else
#    #echo "JEC file found. Residual corrections may be applied."
#    #tar -zxvf JECDataDriven.tar.gz
#    tar -zxf residualMCTruth.tar.gz
#fi


# inputs for run+root script
readFilesExe=${1}
startfile=${2}
endfile=${3}
filelist=${4}
outfile=${5}
radius=${6}
jetType=${7} 
debug=${8}
etaCutLo=${9}
etaCutHi=${10}

#echo ""
#echo "etaCutLo=$etaCutLo"
#echo "etaCutHi=$etaCutHi"
#echo ""
#
#echo "Processing.."
./${readFilesExe} "${filelist}" ${startfile} ${endfile} ${radius} "${jetType}" ${debug} "${outfile}" ${etaCutLo} ${etaCutHi}

#echo "Job done!"
