#!/bin/bash


if [[ $# -eq 1 ]]
then
    
    targdir=$1    
    filearray=($(ls ${targdir}/*.err))
    errLineCount=0
    for i in "${filearray[@]}"
    do
	#cat $i | wc -l
	numLines=($(cat $i | wc -l))
	if [[ $numLines -gt 0 ]]
	then
	    echo "ERRORS for file $i"
	    errLineCount=$(($errLineCount + $numLines))
	fi	
    done
    
    echo "errLineCount=${errLineCount}" ##DEBUG
    
    if [[ $errLineCount -gt 0 ]]
    then
	echo ""
	echo "********* WARNING **********"
	echo "directory ${targdir} has errors!!"
	echo "********* WARNING **********"
	echo ""
    else
	echo ""
	echo "directory ${targdir} has no errors."	
	echo ""
    fi
    
else
    echo ""
    echo "error of usage! do:"
    echo "source checkErrs.sh <targdir>"
    echo ""    
fi

return