
ian w laflotte, started on 6/15/2016, last updated on 8/24/16






**********
The primary analysis folder for 5p02TeV_ppJetAnalysis. This analysis is built on code inherited from Raghav Kunnawalkam's RAA_5p02_code repository. 

See below for notes regarding the use of scripts/code, goals, progress, etc. Some directories are split into data/MC versions (e.g. readFiles, filelists), some not (e.g. DijetResponse). 

some of this may be more readable in github via the RAW format button.







**********
// OVERALL STATUS --------------------------------------------------------

Status: Good; beginning to look a bit closer and interpret analysis results to ask what we can be doing better

Issue: Py8 QCDDijetMC for pp doesnt't have an emulated HLT/L1 sector, can't use/test readFiles_ppMC_normRespMatrix.C
Issue: not vz-weighting the MC correctly






**********
// readFiles --------------------------------------------------------



FILE(s):
pp{Data,MC}/readFiles_pp{Data,MC}.{C,h}
ppData/readFiles_ppData_Run2MinBias.{C,h}

ppMC/readFiles_ppMC_deriveJERUnc.{C,h}
ppMC/normRespMatrix_ppMC.{C,h}
ppMC/jetWeights/jetWeights.{C,h}

printQAPlots/printQAPlots_readFiles_pp{Data,MC}.C



GOAL(s) OF CODE:
create pp{Data,MC} input files for unfoldSpectra, JERandJES, DijetResponse code(s)
make QA/triggerCombination/vz/pthat plots for sanity check for pp{Data,MC}
apply trigger weights to ppData
compute proper MC jet weights for application in readFiles_ppMC
apply vz weights from previous fits to ppMC
apply dijet weights computed from jetWeights.C to ppMC



ppData STATUS:
attempting to implement Run2MinBias version now
issue with Aj/xj computation; should do the first two *valid* jets


ppMC STATUS:
Dijet QCD MC appears to be correctly pthat-weighted, but not correctly vz-weighted

printQAPlots STATUS:
early stages, in progress...

jetWeights STATUS:
fine, no work to be done yet



CURRENTLY WORKING ON:
readFiles_ppData_Run2MinBias.{C,h} to begin to push our observations to lower pt
printQAPlots_Overlaid.C to better compare MC/Data plots







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
// WORKFLOW(s) --------------------------------------------------------

section describing the intended procedures the code + it's output and products should follow. Warning; this section contains some abuse of C++ syntax







**********
// for QA/SanityChecks -------------------------------------------------------- 

if(inputForestFiles are MC){
ForestFiles into-> jetWeights.C outputs-> weightsInTextFormat hardCodeInto-> readFiles_ppMC.C
} 
inputForestFiles into-> readFiles.C outputs-> rootFiles into-> readFiles_printPlots.C outputs-> pdf plots







**********
// for JER uncertainty :

readFiles_ppMC.C rootFiles into-> readFiles_ppMC_deriveJERUnc.C outputs-> pdf plots + rootFiles







**********
// for JER/JES plotting + fitting :

readFiles_ppMC rootFiles into-> readFiles_ppMC_JERandJES_plotNfit.C outputs-> pdf plots + rootFiles







**********
// for data unfolding --------------------------------------------------------

hists from both data AND MC readFiles rootFiles into-> 
unfoldDataSpectra.C outputs->
data unfolding RootFile + plots output







**********
// for MCClosure tests --------------------------------------------------------

hists from readFiles_ppMC rootFiles into-> 
doMCClosureTests.C outputs->
data unfolding RootFile + plots output







**********
// for DijetResponse derive+sum --------------------------------------------------------

inputForestFiles into DijetResponse(deriveResponse) output-> deriveResponse rootFiles
if(suspect Bug/Error/other issue){
deriveResponse rootFiles into-> analyze_deriveResponseOutput.C outputs-> textToScreen about avg of/#entries in each bin
}
else { // deriveResponse stuff looks good
deriveResponse rootFiles into-> DijetResponse(sumResponse) outputs-> rootFiles and (optional) plots
sumResponse output rootFiles into-> DijetResponse_printPlots.C outputs-> sanity/QA plots 
