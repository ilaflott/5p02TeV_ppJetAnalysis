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
not a priority: should look at later when there's some time to get a little distracted