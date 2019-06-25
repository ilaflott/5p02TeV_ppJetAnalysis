#!/bin/bash

expression=${1}
#dirarray=($(ls -d ${PWD}/*/))
if [[ $# -eq 0 ]]
then
    dirarray=($(ls -d */))
elif [[ $# -eq 1 ]]
then
    dirarray=($(ls -d ${expression}/))
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
