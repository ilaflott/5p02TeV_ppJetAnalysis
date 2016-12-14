#!/bin/bash


secondCount=0
while [[ $secondCount -lt 18000 ]]
do
    condor_q ilaflott
    sleep 10s
    secondCount=$(( $secondCount + 5 ))
done
echo "hour of condor-querying done. "
return