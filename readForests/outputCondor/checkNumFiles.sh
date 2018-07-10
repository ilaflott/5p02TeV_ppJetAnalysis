#!/bin/bash


if [[ $# -eq 1 ]]
then
    
    
    targdir=$1
    
    dirarray=($(ls -d ${targdir}/))
    
    for i in "${dirarray[@]}"
    do
	#cat $i | wc -l
	numLines=($(ls ${i}/*.root | wc -l))
	echo "dir = $i"
	echo "has ${numLines} files"	
    done
else
    echo "error, usage! do:"
    echo "source checkNumFiles.sh <targdir>"
    echo ""    
fi

return