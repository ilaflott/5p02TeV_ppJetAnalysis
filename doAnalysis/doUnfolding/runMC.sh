#!/bin/bash

echo ""
echo "compiling unfoldMCSpectra"
echo ""
rooUnfoldCompile unfoldMCSpectra.C

echo ""
echo "done compiling unfoldMCSpectra"
echo ""
sleep 1s


doSVD=0
doBayes=1
useSimpleBinning=0


#### anabins, 30/30 for rec/genpt
### R = 3
source run_unfoldMCSpectra.sh 3 1 20 1 ${doBayes} ${doSVD} ${useSimpleBinning}
sleep 0.5s
### R = 4
source run_unfoldMCSpectra.sh 4 1 20 1 ${doBayes} ${doSVD} ${useSimpleBinning}
sleep 0.5s
#### R = 5
#source run_unfoldMCSpectra.sh 5 1 20 1 ${doBayes} ${doSVD} ${useSimpleBinning}
#sleep 0.5s

useSimpleBinning=1

#### simpbins, 30/30 for rec/genpt
### R = 3
source run_unfoldMCSpectra.sh 3 1 21 1 ${doBayes} ${doSVD} ${useSimpleBinning}
sleep 0.5s
### R = 4
source run_unfoldMCSpectra.sh 4 1 20 1 ${doBayes} ${doSVD} ${useSimpleBinning}
sleep 0.5s
#### R = 5
#source run_unfoldMCSpectra.sh 5 1 20 1 ${doBayes} ${doSVD} ${useSimpleBinning}
#sleep 0.5s




##### anabins, 40/40 for rec/genpt
#### R = 3
#source run_unfoldMCSpectra.sh 3 1 25 1 ${doBayes} ${doSVD} ${useSimpleBinning}
#sleep 0.5s
#### R = 4
#source run_unfoldMCSpectra.sh 4 1 25 1 ${doBayes} ${doSVD} ${useSimpleBinning}
#sleep 0.5s
#### R = 5
#source run_unfoldMCSpectra.sh 5 1 21 1 ${doBayes} ${doSVD} ${useSimpleBinning}
#sleep 0.5s
#
#useSimpleBinning=1

#### simpbins, 40/40 for rec/genpt
#### R = 3
#source run_unfoldMCSpectra.sh 3 1 19 1 ${doBayes} ${doSVD} ${useSimpleBinning}
#sleep 0.5s
#### R = 4
#source run_unfoldMCSpectra.sh 4 1 18 1 ${doBayes} ${doSVD} ${useSimpleBinning}
#sleep 0.5s
#### R = 5
#source run_unfoldMCSpectra.sh 5 1 21 1 ${doBayes} ${doSVD} ${useSimpleBinning}
#sleep 0.5s
#echo ""
#echo "done."
#echo ""




###### anabins, 40/40 for low kReg
#### R = 3
#source run_unfoldMCSpectra.sh 3 1 13 1 ${doBayes} ${doSVD} ${useSimpleBinning}
#sleep 0.5s
#### R = 4
#source run_unfoldMCSpectra.sh 4 1 13 1 ${doBayes} ${doSVD} ${useSimpleBinning}
#sleep 0.5s
#### R = 5
#source run_unfoldMCSpectra.sh 5 1 13 1 ${doBayes} ${doSVD} ${useSimpleBinning}
#sleep 0.5s




echo ""
echo "done."
echo ""


return 