#!/bin/bash

#dirarray=($(ls -d ${PWD}/*/))
if [[ $# -eq 1 ]]
then
    expression=${1}
    dirarray=($(ls -d ${expression}/))
else
    dirarray=($(ls -d pp*/))    
fi

for i in "${dirarray[@]}"
do
    echo ""
    echo ""
    echo "_____________________"
    echo "checking directory $i"
    echo ""
    #echo "num files"
    source checkNumFiles.sh $i
    echo""
    #echo "errs"
    source checkErrs.sh $i
done
