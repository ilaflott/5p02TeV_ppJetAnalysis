//---------------------------------------------------------------------------------------

ian w laflotte
started on 6/15/2016
last updated: (7/12/16)

The primary analysis folder for 5p02TeV_ppJetAnalysis. This analysis is built on code inherited from Raghav Kunnawalkam's RAA_5p02_code repository. 
See below for notes regarding the use of scripts/code, goals, progress, etc. Some directories are split into data/MC versions (e.g. readFiles, filelists), some not (e.g. DijetResponse). 

WORKFLOW(s)
for QA/SanityChecks :  
input Forest files -> 
if(inputForestFiles are MC){jetWeights->weightsInTextFormat->hardCodeInto(readFiles_ppMC)} -> 
readFiles -> rootFiles -> readFiles_printPlots -> pdf plots

for Unfolding : 
hists in readFiles rootFiles -> unfoldSpectra/unfoldDataSpectra.C and/or unfoldSpectra/doMCClosureTests.C -> pdfPlots + unfolding RootFile output

for dijet response calculations : 
inputForestFiles -> DijetResponse(deriveResponse) -> rootfiles -> DijetResponse (sumResponse) -> one rootfile -> 
if(needDebug) analyze_deriveResponseOutput.C -> text output regarding quality of DijetResponse(sum+derive) output 


//----------------------------------------------------------------------------------------

//// unfoldSpectra

FILE(s):
unfoldDataSpectra.C
doMCClosureTests.C
unfoldSpectra.h
utilities.h

GOAL(s): 
remove detector distortions+effects w/ Bayesian+SVD unfolding techniques and a spectra "training" sample

PROGRESS: 
both Bayesian+SVD unfolding techniques create output without crashing quickly
both unfoldDataSpectra.C and doMCClosureTests.C run last on ppData=HighPtJetTriggered dataset and ppMC=Py6 QCDDJet samples w/ 15 < pthat < inf 

CURRENTLY WORKING ON:
deciding if output is sensible


//----------------------------------------------------------------------------------------

//// readFiles 

FILE(s):
ppData/readFiles_ppData.{C,h}
ppData/readFiles_ppData_printPlots.{C,h}
ppMC/readFiles_ppMC.{C,h}
ppMC/readFiles_ppMC_printPlots.{C,h}

GOAL(s):
Make QA/triggerCombination/vz plots for sanity check for ppData+ppMC, create input for unfoldSpectra steps, apply weights to MC events

PROGRESS:
ppData: JetQA/vz/trigger+combo plots created with HighPtJet Triggered datasets up on github in readFiles/ppData folder
ppMC: JetQA/vz/pthat/unfolding plots created with QCDDijet Py6 MC up on github in readFiles/ppMC folder

CURRENTLY WORKING ON:
Making the plots prettier, more readable (a constant work in progress)


//----------------------------------------------------------------------------------------

//// DijetResponse 

FILE(s):
DijetResponse.{C,h}
analyze_deriveResponseOutput.C

GOAL(s): 
Compute detector response to dijets using various methods: 
Aj/Relative response, missEt+jetpt/MPF response, missEt+matchedphoton/MPFAbsPho response, etc.

PROGRESS:


CURRENTLY WORKING ON:



//----------------------------------------------------------------------------------------

//// jetWeights

FILE(s):

GOAL(s): 
Compute Dijet weights for ppMC in different pthat bins

WORKFLOW(s):
jetWeights can be run on ppMC without input from any other script/code

PROGRESS:
all tests run well, dijetQCD MC weights calculated + hardcoded into readFiles_ppMC

CURRENTLY WORKING ON:
N/A



//----------------------------------------------------------------------------------------

//// "2ndTier"

TO WORK/RUN: normResponseMatrix, plot_JESClosure, pp_JER_plot

GOAL(s): not 100% certain yet, but i know this has to do with rooUnfold and the consistency of the corrections we derive for the distributions. More work here is needed.
 
WORKFLOW(s): 
this is the "2nd tier" of the analysis code; must run readFiles/DijetResponse/unfoldSpectra first

PROGRESS:

CURRENTLY WORKING ON:



//----------------------------------------------------------------------------------------

//// bonus round: jetObservables 

GOAL(s): 
some macros dedicated to computing the quanities in question (jetPull, jetThrust, njettiness) and exploring their usefuless as experimental observables. 

CURRENTLY WORKING ON:
not a priority: should look at later when there's some time to get a little distractesd



//----------------------------------------------------------------------------------------

//// RooUnfold Library Implementation

RooUnfold is a tad tricky to implement at cgate/MIT. So i've noted down the entire process here. I tested it in CMSSW_7_5_8/ROOT6 on hidsk

*** Step 1. 

setupenv and download trunk

In a SCRAM-based CMSSW-area's src folder; download the RooUnfold svn trunk with

     source /osg/app/cmssoft/cms/cmsset_default.sh
     cmsenv
     voms-proxy-init -voms cms
     svn co https://svnsrv.desy.de/public/unfolding/RooUnfold/trunk  #downloads the directory "trunk"
     mv trunk RooUnfold-svnTrunk  #rename for convenience

*** Step 2. 

build/compile RooUnfold library/binaries

self-explanatory; do
         
     cd RooUnfold-svnTrunk
     make
     make bin

*** Step 3. export environment variables 

bash won't know where to look for relevant RooUnfold binaries at run-time when code is run unless we alter LD_LIBRARY_PATH manually. its also useful to tell bash exactly where RooUnfold-svnTrunk is for compiling code later; do

     export ROOUNFOLD_PATH=<your RooUnfold directory in a SCRAM-based CMSSW are>
     export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:ROOUNFOLD_PATH

*** Step 4. Compiling+running code that uses RooUnfold

to use RooUnfold in code, make sure to #include the proper RooUnfold headers (e.g. #include "RooUnfoldBayes.h"). 
To compile code including these libraries/headers/functions, one must point g++ to not only the locations of the headers; but also the location of the RooUnfold binary, and to use the RooUnfold library/binary flag
the relvant g++ compilation therefore is

     # NAME is a bash variable that holds the name of the script, e.g. unfoldSpectra.C
     g++ -L"${ROOUNFOLD_PATH}" -I"${ROOUNFOLD_PATH}" -lRooUnfold -I"${ROOUNFOLD_PATH}/src" $NAME $(root-config --cflags --libs) -Werror -Wall -O2 -o "${NAME/%.C/}.exe"

here, the part that says 

-L"${ROOUNFOLD_PATH}" -I"${ROOUNFOLD_PATH}" -lRooUnfold

and 

-I"${ROOUNFOLD_PATH}/src"

points g++ to the lcoation of the RooUnfold libraries/binaries/headers for linking at compile time. I make this a bash function called "rooUnfoldCompile", and it's very similar to my rootcompile function; it differs only by the aforementioned two pieces

********************
*** MY OWN SETUP ***
********************

in my ${CMSSW_BASE}/src do 

   svn co https://svnsrv.desy.de/public/unfolding/RooUnfold/trunk  #downloads the directory "trunk"
   mv trunk RooUnfold-svnTrunk  #rename for convenience
   cd RooUnfold-svnTrunk
   make
   make bin

in .bashrc -> rooUnfoldCompile() function as described above, defined at login

in setupenv.sh -> LD_LIBRARY_PATH and ROOUNFOLD_PATH defined, do

   source setupenv.sh  # no need to enter grid pass

now I should be able to compile code that uses RooUnfold with
   
   rooUnfoldCompile myUnfoldCode.C  # should yield myUnfoldCode.exe
