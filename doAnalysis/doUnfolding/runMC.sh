#!/bin/bash

echo ""
echo "compiling unfoldMCSpectra"
echo ""
rooUnfoldCompile unfoldMCSpectra.C

echo ""
echo "done compiling unfoldMCSpectra"
echo ""
sleep 1s
#if [[ rooUnfoldCompileResult -eq ]]



#### anabin

##### R = 3
##source run_unfoldMCSpectra.sh 3 1 27 1
##sleep 0.5s
##
#### R = 4
##source run_unfoldMCSpectra.sh 4 1 25 1   
##sleep 0.5s
##
##### R = 5
##source run_unfoldMCSpectra.sh 5 1 30 1
##sleep 0.5s


### simpbin

## R = 3
source run_unfoldMCSpectra.sh 3 1 21 1
sleep 0.5s

## R = 4
source run_unfoldMCSpectra.sh 4 1 20 1   
sleep 0.5s

## R = 5
source run_unfoldMCSpectra.sh 5 1 20 1
sleep 0.5s


echo ""
echo "done."
echo ""


return 