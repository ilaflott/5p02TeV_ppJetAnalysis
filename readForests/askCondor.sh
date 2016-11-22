#!/bin/bash


secondCount=0
while [[ $secondCount -lt 3600 ]]
do
    condor_q ilaflott
    sleep 5s
    secondCount=$(( $secondCount + 5 ))
done
echo "hour of condor-querying done. "
return