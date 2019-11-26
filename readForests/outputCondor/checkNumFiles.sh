#!/bin/bash


if [[ $# -eq 1 ]]
then
        
    targdir=$1
    
    #dirarray=($(ls -d ${targdir}/))
    
    #for i in "${dirarray[@]}"
    #do
	#cat $i | wc -l
#	numLines=($(ls ${i}/*.root | wc -l))
    numLines=($(ls ${targdir}/*.root | wc -l))
    errnumLines=($(ls ${targdir}/*.err | wc -l))
    echo "dir = $targdir"
    echo "has ${numLines} root files"	
    echo "and ${errnumLines} err files"	
    if [[ $numLines -lt $errnumLines ]]
    then
	echo "!!!WARNINGWARNINGWARNINGWARNINGWARNING!!! -- looks like some jobs did not finish correctly! -- !!!WARNINGWARNINGWARNINGWARNINGWARNING!!!"
	echo "!!!WARNINGWARNINGWARNINGWARNINGWARNING!!! --   look at the output logs and figure it out!   -- !!!WARNINGWARNINGWARNINGWARNINGWARNING!!!"
    else
	echo "all is well!"
    fi

	#done
else
    echo "error, usage! do:"
    echo "source checkNumFiles.sh <targdir>"
    echo ""    
fi

return
