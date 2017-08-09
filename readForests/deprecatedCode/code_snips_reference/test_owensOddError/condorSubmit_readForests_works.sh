#!/bin/bash

echo ""

NJobs=$1

filelistIn=$2


## NJobs=-1 case
if [[ $NJobs -eq 1 ]]
then    
    echo "submitting jobs for all files in list"
fi

filelist=${filelistIn##*/}
nFiles=`wc -l < $filelistIn`