// Author; Ian Laflotte
// Inspired by code inhereted fro,, Raghav K., Kurt J.
// 
// reads and writes hists from pp data/MC readFiles output to pdf files
// produces jet, trigger, and evt QA spectra
// compile with...
// g++ printAPlots.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readForests_ppData.exe

////////// (initializa/declara)tions //////////
// C++, C, etc.
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
// ROOTSYS
#include <TSystem.h>
#include <TProfile.h>
// I/O
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TEventList.h>
#include <TCut.h>
// utilities
#include <TMath.h>
#include <TRandom3.h>
#include <TStopwatch.h>
// plotting
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TROOT.h>
// histos
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>

// custom
#include "printPlots_strings.h" //const strings, string methods, etc.
#include "printPlots_binsNconsts.h" //variable binning arrays, consts [e.g. int. luminosity]
#include "printPlots_draw.h"  //anything to do with graphics, titles, constant-binnings, axis ranges, styles, etc.
#include "printPlots_TH1.h" //methods that create/alter TH1s
#include "printPlots_evtCounts.h" //methods that have to do with printing out event counts, calculating integrated luminosity etc.

const int CANVX=1000;
const int CANVY=800;

#include "printPlots_jetPlots.h"    
#include "printPlots_jetTrigPlots.h"
#include "printPlots_jetIDPlots.h"
#include "printPlots_jetMultPlots.h"
//#include "printPlots_fitsNfuncs.h" //deprecated
//#include "printPlots_ppMC.h"  //deprecated
//#include "printPlots_ppData.h" // deprecated



