// Author; Ian Laflotte
// Inspired by code inhereted fro,, Raghav K., Kurt J.
// 
// reads and writes hists from pp data/MC readFiles output to pdf files
// produces jet, trigger, and evt QA spectra
// compile with...
// g++ printQAPlots.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readForests_ppData.exe

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
#include "ppMC_etc.h"
#include "ppData_etc.h"

//// combos
const std::string etaWidth = "20_eta_20";

// variable names, string array
const std::string var[] = {   
  "jtpt" ,  "rawpt",  
  "jteta", "jtphi", "jty", //normalize to data
  "trkMax", "trkSum", "trkHardSum", 
  "chMax",  "chSum",  "chHardSum", 
  "phMax",  "phSum",  "phHardSum", 
  "neMax",  "neSum", 
  "eMax",   "eSum", 
  "muMax",  "muSum", 
  "Aj",     "xj" , "dphi", 
  "leadJetPt", "subleadJetPt"
};
const int N_vars = sizeof(var)/sizeof(std::string);

// variable names, string array
const std::string var_xAx_Titles[] = {   
  "p_{t}^{reco} (GeV)" ,  "p_{t}^{raw} (GeV)",  
  "#eta_{jet}", "#phi_{jet}", "y_{jet}" ,
  "trkMax/p_{t}^{raw}", "trkMax/p_{t}^{raw}", "trkHardSum/p_{t}^{raw}", 
  "chMax/p_{t}^{raw}",  "chSum/p_{t}^{raw}",  "chHardSum/p_{t}^{raw}", 
  "phMax/p_{t}^{raw}",  "phSum/p_{t}^{raw}",  "phHardSum/p_{t}^{raw}", 
  "neMax/p_{t}^{raw}",  "neSum/p_{t}^{raw}", 
  "eMax/p_{t}^{raw}",   "eSum/p_{t}^{raw}", 
  "muMax/p_{t}^{raw}",  "muSum/p_{t}^{raw}", 
  "A_{j}",     "x_{j}" , "#Delta #phi", 
  "leadJet P_{t} (GeV)", "subleadJet P_{t} (GeV)"
};
//const int N_vars_xAx_Titles = sizeof(var_xAx_Titles)/sizeof(std::string);

const int var_xAx_reBin[]={                                   
  5, 5,     //reco/rawjtpt //500,0,500	   
  2,  2, 2,    //jteta, phi, jty   //60 ,-3,+3  and  //100,-4,+4  
  5, 5, 5,  //trks	   //200,0,2	   				   
  5, 5, 5,  //ch  
  5, 5, 5,  //ph					           
  5, 5,     //ne		
  5, 5,     //e					      
  5, 5,     //mu					      
  2, 2,     //xj,Aj	//100,0,1				     
  1,        //dphi	//50 ,0,+4		   				      
  5, 5      //sub/lead recojttpt, //500,0,500				      
};

const float jetQAxmax[]={
  (500. ), (500.),
  (3.   ), (4.  ), (3.  ),
  (2.00 ), (2.00), (2.00),
  (2.00 ), (2.00), (2.00),
  (2.00 ), (2.00), (2.00),
  (2.00 ), (2.00), 
  (2.00 ), (2.00), 
  (2.00 ), (2.00), 
  (1.   ), (1.  ), 
  (3.2  ),
  (500. ), (500.)
};


const float jetQAxmin[]={
  30., 30.,  //  0.  ,  0.,
  -3. , -4., -3.,
  0.  ,  0., 0.,
  0.  ,  0., 0.,
  0.  ,  0., 0.,
  0.  ,  0., 
  0.  ,  0., 
  0.  ,  0., 
  0.  ,  0., 
  2.0 ,
  0. , 0.
};

const int jetTrigQABinning=10;
const int jetSpectraRapBinning=10;

// gen variable names for Th2 profiles//MCEff Ploits
const std::string genVars[] = {   
  "pt","eta","phi" //,"drjt"
};
const int N_genVars = sizeof(genVars)/sizeof(std::string);

const std::string genVars_ptrat[] = {   
  "eta","phi" ,"drjt"
};
const int N_genVars_ptrat = sizeof(genVars_ptrat)/sizeof(std::string);

// trigs
const std::string HLTName[] = {
  "HLT40","HLT60","HLT80","HLT100", /*"HLTComb",*/ "TrgCombTest"
};

// trigs
const std::string HLTName_Leg[] = {
  "HLTak4Calo40","HLTak4Calo60","HLTak4Calo80","HLTak4Calo100", /*"HLTSum/4",*/ "HLTComb" 
};
const int N_trigs=sizeof(HLTName)/sizeof(std::string);



//// MC+DATA COMMON BINNING
// ---------------------------------------------------------------------------------------------------------------

// binning arrays
const int ptbins[] = { 15, 30, 50, 80, 120, 170, 220, 300, 500 };
const int nbins_pt = sizeof(ptbins)/sizeof(int)-1;//above values define edges of bins, not centers, so subtract one

const double JEC_ptbins[] = {
  17, 22, 27,    //15-30
  33, 39, 47,    //30-50
  55, 64, 74,    //50-80
  84, 97, 114,   //80-120
  133, 153,      //120-170
  174, 196,      //170-220
  220, 245, 272, //220-300
  300, 350, 400, //300-500
  550, 790, 1000 //500-inf
};
const int nbins_JEC_ptbins = sizeof(JEC_ptbins)/sizeof(double)-1;

const double etabins[] = {
  -5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, 
  -3.839, -3.664, -3.489, -3.314, -3.139, 
  -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, -2.043, 
  -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, 
  -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, 
  -0.087, +0.000, +0.087, 
  +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783, +0.879, +0.957, 
  +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566, +1.653, +1.740, +1.830, +1.930, 
  +2.043, +2.172, +2.322, +2.500, +2.650, +2.853, +2.964, 
  +3.139, +3.314, +3.489, +3.664, +3.839, 
  +4.013, +4.191, +4.363, +4.538, +4.716, +4.889, +5.191
};
const int nbins_eta = sizeof(etabins)/sizeof(double)-1;

const float rapbins[]={ 
  0.0, 0.5,
  1.0, 1.5,
  2.0, 2.5,
  3.0 //, 3.2,
  //3.7
};
const int nbins_rap=sizeof(rapbins)/sizeof(float)-1;

//// HELPER FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------

void divideBinWidth(TH1 *h){
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();++i){//binsX loop 
    Float_t val = h->GetBinContent(i);
    Float_t valErr = h->GetBinError(i);
    val/=h->GetBinWidth(i);
    valErr/=h->GetBinWidth(i);
    h->SetBinContent(i,val);
    h->SetBinError(i,valErr);
  }//end nbinsX loop
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  return;
}
