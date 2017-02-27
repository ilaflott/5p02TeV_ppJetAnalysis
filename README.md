
ian w laflotte, started on 6/15/2016, last updated on 2.27.2017






// ABOUT THE REPO --------------------------------------------------------
------------------------

The primary analysis folder for 5p02TeV_ppJetAnalysis github repository. This repository is built on code inherited from Raghav Kunnawalkam's RAA_5p02_code repository. 

See below for notes regarding the use of scripts/code, goals, progress, etc. Some directories are split into data/MC versions (e.g. readForests, filelists), some not (e.g. DijetResponse). 

Some of this may be more readable in github via the RAW format button.





// ABOUT readForests --------------------------------------------------------
------------------------

purpose; 
to run over HiForests on hidsk, either locally or via condor batch, and create root files with histograms of interest for later plotting and analysis.

------------------------
// primaryCODE(s) -----
------------------------

readForests_ppMC_jetPlots.C
readForests_ppData_jetPlots.C
readForests_jetPlots.h
readForests_ppMC_MCJEC.Cr
readForests_MCJEC.h


// secondaryCODE(s) -----
------------------------

L2L3ResidualWFits.h
MCTruthResidual.h
JECDataDriven.tar.gz
residualMCTruth.tar.gz


// steering/condorMACRO(s) -----
------------------------

askCondor.sh
condorRun_readForests.sh
condorSubmit_readForests
run_readForests_jetPlots.sh
run_readForests_ppMC_MCJEC.sh



----------------------------------
// HOW TO SET UP readForests -----
----------------------------------

note; $CMSSW_BASE/src is a SCRAM-based area of your choice, no particular CMSSW required (for reference, 7_5_8 generally used)

```
cd $CMSSW_BASE/src
source /osg/app/cmssoft/cms/cmsset_default.sh
cmsenv
voms-proxy-init -voms cms
export CVMFS_758=/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_5_8/src/
cd readForests
tar -zxvf JECDataDriven.tar.gz
tar -zxvf residualMCTruth.tar.gz
```


----------------------------------
// HOW TO COMPILE readForests -----
----------------------------------

note; readForests_ppData_jetPlots.C used for example, but other .C's follow the same scheme

```
g++ readForest_ppData_jetPlots.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readForests_ppData_jetPlots.exe
```


----------------------------------------
// HOW TO TEST readForests LOCALLY -----
----------------------------------------

note; to test readForests on a file that resides on hadoop, with default parameters (see header for details)
```
./readForests_ppData_jetPlots.exe
```


note; to test on a file locally, edit the filelists "filelists/test_readForests_*_local.txt" to point to a local test file of your choice, then run the .exe and point to the local-test-file filelist of interest
```
./readForests_ppData_jetPlots.exe filelists/test_readForests_ppData_Jet80_local.txt 0 0 4 PF 1 test_readForests_ppData_jetPlots.root
```


note; general usage of the executable is output to screen by the executable itself if an incorrect # of arguments is given to the executable, e.g.
```
./readForests_ppData_jetPlots.exe HELP
```


note; the output fo the above command is
```
for actually running, do...
./readForests_ppData_jetPlots.exe <inputFileList> <startFile> <endFile> <jetRadius> <jetType> <debugMode> <outputFilename> [<absEtaCutLo> <absEtaCutHi>]

where args in [] are optional.


rStatus=-1
```

note; local tests and debugging of readForests should happen before one wishes to test on condor. 


---------------------------------------
// HOW TO TEST readForests ON CONDOR -----
---------------------------------------

note; so lets say you've set up readForests and you've acquired a local test file or two for quick, local debugging, and now you actually want to run an analysis job to learn something about your data. Though possible in theory, locally running one big job can work. But it can also take a long time, hog a lot of cluster memory (despite "nice"ing the command), and require starting from scratch upon failure after hours of running. Thus, running one job locally is a recipe for wasted time, avoid at all costs!


note; This is what the "steering/condor macros" are for; they can cut up a giant local job into N batch jobs of size M, for N and M of your choice. before running readForests with these macros, create and outputCondor folder 
```
cd $CMSSW_BASE/src/readForests
mkdir outputCondor
```


note; now, assuming you've compiled the readForests script you want to use, let's try submitting a single job to condor. The following 
~ calls the script/macro (input 0: source condorSubmit_readForests.sh) 
~ points to readForests_ppData_jetPlots.C (input 1: readForests_ppData_jetPlots)
~ submits one job (input 2: 1)
~ running over 10 files from the filelist (input 3: 10)
~ starting at the 0'th file (input 4: 0)
~ using HighPtJet80's filelist (input 5: filelists/filelists/5p02TeV_HighPtJet80_forests.txt)
~ reading the ak4PF Jet collection (inputs 6 and 7: 4 PF)
~ debug mode on (input 8: 1)
~ keeping jets with abs(eta) between 0 and 4.7 (optional inputs 9 and 10: 0 4.7)
```
source condorSubmit_readForests.sh readForests_ppData_jetPlots 1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 1  0.0 4.7
```


note; you should see some output to the screen confirming your submission to condor. to monitor progress, check the outputcondor folder for output. You can also ask condor directly. to monitor jobs without constantly querying eithe folder or condor yourself by hand, edit the askCondor script to use your MIT username and source it to look at your condor_q every few seconds or so. One job with 10 files shouldn't take more than 20 minutes to run, and on average, runs much more quickly than that, pending demand on the cluser
```
ls outputCondor/                       # check folder for output
condor_q <your MIT username>           # ask condor
source askCondor.sh                    # ask condor every few seconds, use ctrl-c to quit
```


note; everytime edits to readForests pass local tests, one should also test them on condor with a single, small test job. There are issues that will pop up on condor that won't pop up locally, and so this should not be skipped. condorSubmit_readForests.sh copies *.C and *.h files to directory created in outputCondor, for tracing down bugs after-the-fact. To see test job output from above test job submission, do
```
cat outputCondor/ppData_HighPtJet80_ak4PFJets_02-27-17_jetPlots/*.out     # output from code itself, often helpful
cat outputCondor/ppData_HighPtJet80_ak4PFJets_02-27-17_jetPlots/*.err     # error output from code itself, often helpful
cat outputCondor/ppData_HighPtJet80_ak4PFJets_02-27-17_jetPlots/*.log     # condor status during running of code, rarely helpful
root -l outputCondor/ppData_HighPtJet80_ak4PFJets_02-27-17_jetPlots/*.root # open root output file from job, check for empty histograms, etc., often helpful
```




-------------------------------------------------------
// RUNNING ONE SET OF readForests JOBS ON CONDOR -----
-------------------------------------------------------

note; So your code has passed local tests, and a condor test. Nothing in the *.err, the test-root file has stuff in (at least some of the) histograms. Now you want to run over the entire filelist. Take the condorSubmit command you used before, and turn the first numerical argument (input 2) from '1' to '-1', this tells the script to keep submitting 10-file-large jobs until it hits the end of the specified dataset. All the output from this one condorSubmit is put under one directory.


note; usually ppMC just needs one condorSubmit job. I usually submit something like;
```
source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0 4.7
```


note; usually ppData results are created combining output from Jet80 and LowerJets filelists. So often one "set" of readForests ppData jobs takes two condorSubmit jobs
```
source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 10 0 filelists/5p02TeV_HighPtJet80_forests.txt 4 PF 0 0.0 4.7
source condorSubmit_readForests.sh readForests_ppData_jetPlots -1 10 0 filelists/5p02TeV_HighPtLowerJets_forests.txt 4 PF 0 0.0 4.7
```


note; assuming a normal exit status for the jobs, and no run-time issues. You'll usually want to hadd the output in the folders in outputCondor. in readForests/outputCondor, for ppMC, do something like
```
hadd -f ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-27-17_jetPlots/Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root      ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-27-17_jetPlots/*.root
```

note; for ppData, assuming you're combining jet80 and LowerJets results, do something like
```
mkdir ppData_HighPtJetTrig_ak4PFJets_02-27-17_jetPlots/
mv ppData_HighPtJet80_ak4PFJets_02-27-17_jetPlots/ ppData_HighPtJetTrig_ak4PFJets_02-27-17_jetPlots/.
mv ppData_HighPtLowerJets_ak4PFJets_02-27-17_jetPlots/ ppData_HighPtJetTrig_ak4PFJets_02-27-17_jetPlots/.
hadd -f ppData_HighPtJetTrig_ak4PFJets_02-27-17_jetPlots/HighPtJetTrig_ak4PF-allFiles.root ppData_HighPtJetTrig_ak4PFJets_02-27-17_jetPlots/ppData_HighPt*/*.root
```


note; the hadd'd output file is what you want to point your printPlots macros to (combine, jetID eff code from raghav, +others, see doAnalysis/printPlots). After awhile, the condor output can get big, so i recommend utlizing space under your username on scratch for longer term storage/later reference. The directory is /export/d00/scratch/<your username>



-------------------------------------------------------
// RUNNING MANY SETs OF readForests JOBS ON CONDOR -----
-------------------------------------------------------

******WARNING!****** 
CAN GENERATE LOTS OF OUTPUT AND USE A LOT OF RESOURCES! 
ONLY USE IF COMFORTABLE WITH PREVIOUS WORKFLOWS DESCRIBED. 
ONLY USE IF YOUR CODE HAS BEEN TESTED LOCALLY AND ON CONDOR.

note; now lets say you're fine with running one set of condor jobs. Let's say you want to run a bunch of sets at once to get them running before you leave the office for the night, so you can look at fresh output in the morning. take a look at the macro run_readForests_jetPlots.sh. You'll notice it looks like many commands, each running one set of jobs. When I know what I want to look at and I trust my code, I make a bunch of slightly different versions of the readForests code under slightly different names, or use one version I want to use to look at different eta ranges, and then write up each submission command in the run_readForest_jetPlots.sh script. Then, to submit all at once, i do
```
source run_readForests_jetPlots.sh
```


note; there's already stuff in run_readForests_jetPlots.sh from the last time I had to submit a lot of jobs, but there's two main forms I find useful. Below is one version of run_readForests_jetPlots.sh running the ppMC version of readForests over the Py8 filelist using different eta ranges. The "sleep" command asks the terminal to just sit for a specified amount of time, I put it there to be nicer to the condor client that handles your job submissions and whatnot, or to give me another opportunity to "ctrl-C" and stop whatever I'm doing to re-consider.
```
#!/bin/bash
source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0.0 2.4
sleep 2s
source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 2.4 2.7
sleep 2s
source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 2.7 3.0
sleep 2s
source condorSubmit_readForests.sh readForests_ppMC_jetPlots -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 3.2 4.7
```


note; here's a version where I edited the code into three different versions that I know will work, and I want to compile + submit on the fly
```
#!/bin/bash
rootcompile readForests_ppMC_jetPlots_v1.C
source condorSubmit_readForests.sh readForests_ppMC_jetPlots_v1 -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0.0 4.7
sleep 2s

rootcompile readForests_ppMC_jetPlots_v2.C
source condorSubmit_readForests.sh readForests_ppMC_jetPlots_v2 -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0.0 4.7
sleep 2s

rootcompile readForests_ppMC_jetPlots_v3.C
source condorSubmit_readForests.sh readForests_ppMC_jetPlots_v3 -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0.0 4.7
sleep 2s

rootcompile readForests_ppMC_jetPlots_v4.C
source condorSubmit_readForests.sh readForests_ppMC_jetPlots_v4 -1 30 0 filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt 4 PF 0 0.0 4.7
```










-------------------------------------------------------
// UNDER CONSTRUCTION                             -----
-------------------------------------------------------

// ABOUT printPlots --------------------------------------------------------

purpose;
to take in hadd'd output from readForests code and conduct analysis/create graphical output

 




-------------------------------------------------------
// OLD, IGNORE FOR NOW BELOW THIS LINE            -----
-------------------------------------------------------


**********
// unfoldSpectra --------------------------------------------------------

FILE(s):
unfoldDataSpectra.C
doMCClosureTests.C
unfoldSpectra.h
utilities.h

GOAL(s) OF CODE: 
remove detector distortions+effects w/ Bayesian+SVD unfolding techniques and a spectra "training" sample

STATUS:
output looks largely good so far, but the bayesian unfolding could be better
there are better kReg values for the SVD unfolding to pick other than 20, figure it out!

CURRENTLY WORKING ON:
making more Bayesian Unfolding plots in both codes, want to examine the fitting across several different kIter values






**********
// DijetResponse --------------------------------------------------------

FILE(s):
DijetResponse.{C,h}
DijetResponse_printPlots.C
analyze_deriveResponseOutput.C

GOAL(s) OF CODE: 
Compute detector response to dijets using various methods: 
Aj/Relative response, missEt+jetpt/MPF response, missEt+matchedphoton/MPFAbsPho response, etc.

STATUS:
compiles + runs without issue on both ppData and ppMC events
Rel Response looks okay, except @ low pt, why?
MPF look has the flattest curves at 1 across pt bins for both ppData and ppMC
MPFAbsPho has flat curves with a tinier eta range, less than one

CURRENTLY WORKING ON:
making plots pretty + interpreting the output
understanding the low Rel Response curve for low jetpt on ppData and ppMC