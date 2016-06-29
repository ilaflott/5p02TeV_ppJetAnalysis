README_more.md = (README.md) - (bugs + details + github hosting)
this file is private notes for me only and the main one that gets updated. The main parts will regularly be copied/pasted into README.md that is worth pushing to the repo

started on 6/15/2016, ian laflotte

The primary analysis folder for 5p02TeV_ppJetAnalysis. This analysis is built on code inherited from Raghav Kunnawalkam's RAA_5p02_code repository. 
Some pieces of code are split into data/MC versions (e.g. readFiles, filelists), some not (e.g. DijetResponse). 
See below for notes regarding the use of scripts/code, goals, progress, etc.



//----------------------------------------------------------------------------------------

//// readFiles 
//// (these notes last updated 6/15/16)

GOAL(s): Make QA/triggerCombination/vz plots for ppData,MC. Compute weights for ppMC.

WORKFLOW(s)
ppData: can be run on data as is, no input files necessary.
ppMC: must use output from jetWeights in order to run properly. Otherwise the MC spectra are unweighted and nonsensical.

CURRENTLY WORKING ON
how to do offline "HLT_AK4CaloJet{40,60,80,100}_Eta5p1_v" on the ppMC since the hltanalysis/hltobject trees aren't there
- perhaps i should figure out how to run that part myself?
Making printPlots for readFiles_pp{MC,Data} print more informative, detailed, and eye-pleasing plots




//----------------------------------------------------------------------------------------

//// DijetResponse 
//// (these notes last updated 6/15/16)

GOAL(s): Compute detector response to dijets using various methods (use Aj/Relative response, missEt+jetpt/MPF response, missEt+matchedphoton/MPFAbsPho response).

WORKFLOW(s)
works on both ppData and ppMC

CURRENTLY WORKING ON
still *not quite sure* if ppMC needs weights to run properly; i'm leaning towards no but i should look a bit more into it 




//----------------------------------------------------------------------------------------

//// unfoldSpectra
//// (these notes last updated 6/15/16)

GOAL(s): remove detector distortions+effects due to nonlinear behaviors and reconstruction inefficiencies. 

WORKFLOW(s)
both readFiles ppMC AND ppData output files are fed into -> unfoldSpectra.C

CURRENTLY WORKING ON
getting output from unfoldDataSpectra
getting doMCClosureTest to compile




//----------------------------------------------------------------------------------------

//// other
//// (these notes last updated 6/15/16)

TO WORK/RUN: normResponseMatrix, plot_JESClosure, pp_JER_plot

GOAL(s): not 100% certain yet, but i know this has to do with rooUnfold and the consistency of the corrections we derive for the distributions. More work here is needed.
 
WORKFLOW(s): 
this is the "2nd tier" of the analysis code; must run readFiles/DijetResponse/unfoldSpectra first

CURRENTLY WORKING ON:





//----------------------------------------------------------------------------------------

//// jetWeights
//// (these notes last updated 6/15/16)

GOAL(s): Compute weights for ppMC.

WORKFLOW(s)
jetWeights can be run on ppMC without input from any other script/code

CURRENTLY WORKING ON
N/A, all tests ran well, dijetQCD MC weights calculated.




//----------------------------------------------------------------------------------------

//// bonus round: jetObservables 
//// (these notes last updated 6/15/16)

GOAL: some macros dedicated to computing the quanities in question (jetPull, jetThrust, njettiness) and exploring their usefuless as experimental observables. 

WORKFLOW(s)
should these macros read forest files, or readFiles output? i'm thinking ppMC will need the weights; i think they should probably read the forest

CURRENTLY WORKING ON
not until i hit my priorities (readFiles_pp{MC,Data},DijetResponse,unfoldSpectra,jetWeights,pieces of code currently in other category.)

//----------------------------------------------------------------------------------------




//// RooUnfold Implementation
//// (these notes last updated 6/15/16)

RooUnfold is a tad tricky to implement at cgate/MIT. So i've noted down the entire process here. I tested it in CMSSW_7_5_8/ROOT6 on hidsk

********************
*** FROM SCRATCH ***
********************

Step 1. 
setupenv and download trunk

In a SCRAM-based CMSSW-area's src folder; download the RooUnfold svn trunk with

     source /osg/app/cmssoft/cms/cmsset_default.sh
     cmsenv
     voms-proxy-init -voms cms
     svn co https://svnsrv.desy.de/public/unfolding/RooUnfold/trunk  #downloads the directory "trunk"
     mv trunk RooUnfold-svnTrunk  #rename for convenience

Step 2. 
build/compile RooUnfold library/binaries

self-explanatory; do
         
     cd RooUnfold-svnTrunk
     make
     make bin

Step 3. export environment variables 

bash won't know where to look for relevant RooUnfold binaries at run-time when code is run unless we alter LD_LIBRARY_PATH manually. 
its also useful to tell bash exactly where RooUnfold-svnTrunk is for compiling code later; do

     export ROOUNFOLD_PATH=<your RooUnfold directory in a SCRAM-based CMSSW are>
     export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:ROOUNFOLD_PATH

Step 4. Compiling+running code that uses RooUnfold

to use RooUnfold in code, make sure to #include the proper RooUnfold headers (e.g. #include "RooUnfoldBayes.h"). 
To compile code including these libraries/headers/functions, one must point g++ to not only the locations of the headers; but also the location of the RooUnfold binary, and to use the RooUnfold library/binary flag
the relvant g++ compilation therefore is

     # NAME is a bash variable that holds the name of the script, e.g. unfoldSpectra.C
     g++ -L"${ROOUNFOLD_PATH}" -I"${ROOUNFOLD_PATH}" -lRooUnfold -I"${ROOUNFOLD_PATH}/src" $NAME $(root-config --cflags --libs) -Werror -Wall -O2 -o "${NAME/%.C/}.exe"

here, the part that says 
-L"${ROOUNFOLD_PATH}" -I"${ROOUNFOLD_PATH}" -lRooUnfold
and 
-I"${ROOUNFOLD_PATH}/src"
points g++ to the lcoation of the RooUnfold libraries/binaries/headers for linking at compile time.
I make this a bash function called "rooUnfoldCompile", and it's very similar to my rootcompile function; it differs only by the aforementioned two pieces

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
