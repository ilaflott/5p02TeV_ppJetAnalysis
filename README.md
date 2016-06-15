started on 6/15/2016, ian laflotte

The primary analysis folder for 5p02TeV_ppJetAnalysis. This analysis is built on code inherited from Raghav Kunnawalkam's RAA_5p02_code repository. Each folder represents a major arm of the analysis. Some are split into data/MC versions (e.g. readFiles, filelists), some not (e.g. DijetResponse). See below for notes regarding the use of scripts/code, goals, progress, etc.




//----------------------------------------------------------------------------------------

//// readFiles 
//// (these notes last updated 6/15/16)

GOAL(s): Make QA/triggerCombination/vz plots for ppData,MC. Compute weights for ppMC.

WORKFLOW(s)
ppData: can be run on data as is, no input files necessary.
ppMC: must use output from readFiles/ppMC/jetWeights in order to run properly. Otherwise the MC spectra are unweighted and nonsensical.

CURRENTLY WORKING ON
writing+testing+optimizing readFiles_ppMC and jetWeights
Making printPlots for readFiles_ppData print more informative, detailed, and eye-pleasing plots





//----------------------------------------------------------------------------------------

//// jetWeights
//// (these notes last updated 6/15/16)

GOAL(s): Compute weights for ppMC.

WORKFLOW(s)
jetWeights can be run on ppMC without input from any other script/code

CURRENTLY WORKING ON
writing+testing+optimizing jetWeights




//----------------------------------------------------------------------------------------

//// DijetResponse 
//// (these notes last updated 6/15/16)

GOAL(s): Compute detector response to dijets using various methods (use Aj/Relative response, missEt+jetpt/MPF response, missEt+matchedphoton/MPFAbsPho response).

WORKFLOW(s)
ppData: requires no outside input to run. Run deriveResponse on ppData via condor, then hadd the output, then run rumResponse on the hadd'd output.
ppMC: to be written

CURRENTLY WORKING ON
will add notes to ppMC workflow once i know explicitly what it needs. Should be very similar if not idential to ppData workflow.




//----------------------------------------------------------------------------------------

//// rooUnfold 
//// (these notes last updated 6/15/16)

GOAL(s): remove detector distortions+effects due to nonlinear behaviors and reconstruction inefficiencies. 

WORKFLOW(s)
not entirely certain, but i believe i need both MC and Data input to run this. Then i thing i need to apply the output to the samples

CURRENTLY WORKING ON
understanding how it works, how to use, etc.




//----------------------------------------------------------------------------------------

//// other
//// (these notes last updated 6/15/16)

TO WORK/RUN: normResponseMatrix, plot_JESClosure, pp_JER_plotunsure

GOAL(s): not 100% certain yet, but i know this has to do with rooUnfold and the consistency of the corrections we derive for the distributions. More work here is needed.
 
WORKFLOW(s): 
N/A

CURRENTLY WORKING ON:
N/A




//----------------------------------------------------------------------------------------

//// bonus round: jetObservables 
//// (these notes last updated 6/15/16)

GOAL: some macros dedicated to computing the quanities in question (jetPull, jetThrust, njettiness) and exploring their usefuless as experimental observables. 

WORKFLOW(s)
anticipated: ppData/MC versions of each. MC versions probably need the weights from jetWeights.

CURRENTLY WORKING ON
not until i hit my priorities (readFiles,DijetResponse,rooUnfold, and more, for both ppData and ppMC)