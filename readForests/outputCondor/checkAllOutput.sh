#!/bin/bash


#dirarray=($(ls -d ${PWD}/*/))
dirarray=($(ls -d */))
for i in "${dirarray[@]}"
do
    echo "checking directory $i"
    echo ""
    echo "num files"
    source checkNumFiles.sh $i
    echo "errs"
    echo ""
    source checkErrs.sh $i
done
