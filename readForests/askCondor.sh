#!/bin/bash


secondCount=0
while [[ $secondCount -lt 3000 ]]
do
    condor_q ilaflott
    sleep 30s
    secondCount=$(( $secondCount + 5 ))
done
echo "hour of condor-querying done. "
return