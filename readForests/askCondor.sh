#!/bin/bash


secondCount=0
while [[ $secondCount -lt 1800 ]]
do
    condor_q ilaflott
    sleep 3s
    secondCount=$(( $secondCount + 3 ))
done
echo "half-hour of condor-querying done. "
return