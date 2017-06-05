#!/bin/bash

echo ""
echo "compiling unfoldDataSpectra"
echo ""
rooUnfoldCompile unfoldDataSpectra.C

echo ""
echo "done compiling unfoldDataSpectra"
echo ""
sleep 1s


## R = 5
#source run_unfoldDataSpectra.sh 5 0 0 9 1
#sleep 0.5s
#source run_unfoldDataSpectra.sh 5 1 0 9 1
#sleep 0.5s
#source run_unfoldDataSpectra.sh 4 0 0 10 2
#sleep 0.5s
#source run_unfoldDataSpectra.sh 4 1 0 10 2
#sleep 0.5s

## R = 4
source run_unfoldDataSpectra.sh 4 0 0 9 1
sleep 0.5s
source run_unfoldDataSpectra.sh 4 1 0 9 1
sleep 0.5s
#source run_unfoldDataSpectra.sh 4 0 0 10 2
#sleep 0.5s
#source run_unfoldDataSpectra.sh 4 1 0 10 2
#sleep 0.5s

## R = 3
source run_unfoldDataSpectra.sh 3 0 0 9 1
sleep 0.5s
source run_unfoldDataSpectra.sh 3 1 0 9 1
sleep 0.5s
#source run_unfoldDataSpectra.sh 3 0 0 10 2
#sleep 0.5s
#source run_unfoldDataSpectra.sh 3 1 0 10 2
#sleep 0.5s





echo ""
echo "done."
echo ""


return 