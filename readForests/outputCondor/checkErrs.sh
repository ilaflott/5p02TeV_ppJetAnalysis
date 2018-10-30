#!/bin/bash


if [[ $# -eq 1 ]]
then
    
    targdir=$1
    
    filearray=($(ls ${targdir}/*.err))
    
    for i in "${filearray[@]}"
    do
	#cat $i | wc -l
	numLines=($(cat $i | wc -l))
	if [[ $numLines -gt 0 ]]
	then
	    echo "ERRORS for file $i"
	    #echo "has errors!"
	    #echo "numLines=${numLines}"
	fi
    done
else
    echo "error, usage! do:"
    echo "source checkErrs.sh <targdir>"
    echo ""    
fi

return