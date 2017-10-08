#!/bin/bash

Nasks=$1
T=$2
Ncount=0

#secondCount=0
while [[ $Ncount -lt $Nasks ]]
do
    cd outputCondor
    ./checkProgress.sh 1 1
    cd ..

    condor_q ilaflott
    sleep ${T}s

#    secondCount=$(( $secondCount + 5 ))

    Ncount=$(( $Ncount + 1 ))
done
#echo "hour of condor-querying done. "
echo ""
echo "done asking condor."
echo ""
return