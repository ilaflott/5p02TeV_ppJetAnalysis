// C/C++
#include <cstdlib>
#include <cstring>
#include <iostream>

// file IO/settings
#include <TROOT.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TFile.h>

// fit + boundaries + plots
#include <TMath.h>
#include <TH1F.h>
#include <TF1.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>

// draw/print
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TLine.h>
#include <TLatex.h>

#include "JERandJES_binsNconsts.h"
#include "JERandJES_fitsNfuncs.h"
#include "JERandJES_IO.h"

// run settings
const bool printDebug=true;
const bool doBjets = false; //true;
const int iSave=1;
const int wJetID=1;  
//int id = wJetID;
