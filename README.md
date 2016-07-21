
ian w laflotte, started on 6/15/2016, last updated on 7/21/16

**********
The primary analysis folder for 5p02TeV_ppJetAnalysis. This analysis is built on code inherited from Raghav Kunnawalkam's RAA_5p02_code repository. 
See below for notes regarding the use of scripts/code, goals, progress, etc. Some directories are split into data/MC versions (e.g. readFiles, filelists), some not (e.g. DijetResponse). 

**********
// OVERALL STATUS --------------------------------------------------------

positive, most code compiles+runs without issue. 
Most output looks sensbile (DijetResponse_printPlots, readFiles_ppd{Data,MC}_printPlots)
Issue: Py6 QCDDijetMC for pp doesnt't have an emulated HLT/L1 sector, can't use/test readFiles_ppMC_normRespMatrix.C
putting MinBias samples into the spectra+fitter and then re-running the analysis shouldn't be an issue

**********
// readFiles --------------------------------------------------------

FILE(s):
pp{Data,MC}/readFiles_pp{Data,MC}.{C,h}
pp{Data,MC}/readFiles_pp{Data,MC}_printPlots.C
ppMC/readFiles_ppMC_deriveJERUnc.{C,h}
ppMC/readFiles_ppMC_normRespMatrix.{C,h}
ppMC/readFiles_ppMC_JERandJES_plotNfit.{C,h}
ppMC/readFiles_ppMC_JERandJES_fitsNfuncs.h
ppMC/readFiles_ppMC_JERandJES_binsNconsts.h

GOAL(s) OF CODE:
create pp{Data,MC} input for unfoldSpectra steps and apply pthat/vz/etc. 
apply weights computed from jetWeights.C to ppMC
make QA/triggerCombination/vz/pthat plots for sanity check for pp{Data,MC}
plotting, fitting, and deriving JER/JES + response matrix normalization for MC events

ppData STATUS:
JetQA/vz/trigger+combo plots created with HighPtJet Triggered datasets up on github in readFiles/ppData folder
output ntuples can be fed into unfolding code just fine

ppMC STATUS:
JetQA/vz/pthat/unfolding plots created with QCDDijet Py6 MC up on github in readFiles/ppMC folder
Dijet QCD MC appears to be correctly weighted
output ntuples can be fed into unfolding code just fine
normRespMatrix.C DOES NOT WORK! QCDDijet Py6 MC doesn't have req. HLT info
JERandJES_plotNfit/fitsNfuncs/binsNconsts.{C,h} scripts compile+work fine. unsure if output is sane
deriveJERUnc.{C,h} compiles+work fine. unsure if output is sane

CURRENTLY WORKING ON:
Making plots from printPlots more enlightening (more info on drawn PDF file, better binning, maybe!)
interpreting output from the JER/JES scripts, also making it prettier + more readable



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


**********
// jetWeights --------------------------------------------------------

FILE(s):
jetWeights.C

GOAL(s) OF CODE: 
Compute Dijet weights for ppMC in different pthat bins, output the event weights into a text file

STATUS:
all tests run well, dijetQCD MC weights calculated then hardcoded into readFiles_ppMC

CURRENTLY WORKING ON:



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
unfoldDataSpectra.C and do MCClosureTests

CURRENTLY WORKING ON:
interpreting + beautifying the output plots from both unfoldDataSpectra.C and doMCClosureTests.C


**********
// jetObservables --------------------------------------------------------

DIRECTORY(s):
jetPull/
jetThrust/
nSubJettiness/

GOAL(s) OF CODE: 
some macros dedicated to computing the quanities in question (jetPull, jetThrust, njettiness) and exploring their usefuless as experimental observables. 

CURRENTLY WORKING ON, STATUS, CURRENTLY WORKING ON:
N/A, later.

**********
// workflows --------------------------------------------------------

section describing the intended procedures the code + it's output and products should follow. Warning; this section contains some abuse of C++ syntax

**********
WORKFLOW(s) for QA/SanityChecks :  

if(inputForestFiles are MC){
ForestFiles into-> jetWeights.C outputs-> weightsInTextFormat hardCodeInto-> readFiles_ppMC.C
} 
inputForestFiles into-> readFiles.C outputs-> rootFiles into-> readFiles_printPlots.C outputs-> pdf plots


**********
WORKFLOW(s) for JER uncertainty :

readFiles_ppMC.C rootFiles into-> readFiles_ppMC_deriveJERUnc.C outputs-> pdf plots + rootFiles


**********
WORKFLOW(s) for for JER/JES plotting + fitting :

readFiles_ppMC rootFiles into-> readFiles_ppMC_JERandJES_plotNfit.C outputs-> pdf plots + rootFiles


**********
WORKFLOW(s) for data unfolding : 

hists from both data AND MC readFiles rootFiles into-> 
unfoldDataSpectra.C outputs->
data unfolding RootFile + plots output


**********
WORKFLOW(s) for MCClosure tests : 

hists from readFiles_ppMC rootFiles into-> 
doMCClosureTests.C outputs->
data unfolding RootFile + plots output


**********
WORKFLOW(s) for DijetResponse derive+sum : 

inputForestFiles into DijetResponse(deriveResponse) output-> deriveResponse rootFiles
if(suspect Bug/Error/other issue){
deriveResponse rootFiles into-> analyze_deriveResponseOutput.C outputs-> textToScreen about avg of/#entries in each bin
}
else { // deriveResponse stuff looks good
deriveResponse rootFiles into-> DijetResponse(sumResponse) outputs-> rootFiles and (optional) plots
sumResponse output rootFiles into-> DijetResponse_printPlots.C outputs-> sanity/QA plots 
