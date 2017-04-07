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
#include "printPlots_ppMC.h"
#include "printPlots_ppData.h"
#include "printPlots_binsNconsts.h"
#include "printPlots_fitsNfuncs.h"

const std::string CMSSW_BASE=
  "/net/hisrv0001/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/";
const std::string SCRATCH_BASE=
  "/export/d00/scratch/ilaflott/5p02TeV_ppJetAnalysis/";

const std::string inputDir=SCRATCH_BASE+"readForests/";
//const std::string inputDir=CMSSW_BASE+"readForests/outputCondor/";

const std::string outputDir=CMSSW_BASE+"doAnalysis/printPlots/output/";


// variable names for QA Plots
const std::string var[]={
  "jtpt"  , "rawpt" ,  //jets, 0-3=4
  "jteta" , "jtphi" , //"jty",
  "trkMax", "trkSum", "trkHardSum", //jet constituents, 4-21=18
  "chMax" , "chSum" , "chHardSum",
  "phMax" , "phSum" , "phHardSum",
  "neMax" , "neSum" ,
  "eMax"  , "eSum"  ,
  "muMax" , "muSum" ,
  "neN", "chN", "sumN",
  "Aj" , "xj" , "dphi", //dijet variables, 22-26=5
  "leadJetPt", "subleadJetPt"
};
const int N_vars=sizeof(var)/sizeof(std::string);
const int jetQA_varStart=0;
const int jetConstits_varStart=4;
const int dijet_varStart=22;

// variable names, string array
const std::string var_xAx_Titles[] = {   
  "p_{t}^{reco} (GeV)" ,  "p_{t}^{raw} (GeV)",  
  "#eta_{jet}", "#phi_{jet}", //"y_{jet}" ,
  "trkMax/p_{t}^{raw}", "trkMax/p_{t}^{raw}", "trkHardSum/p_{t}^{raw}", 
  "chMax/p_{t}^{raw}",  "chSum/p_{t}^{raw}",  "chHardSum/p_{t}^{raw}", 
  "phMax/p_{t}^{raw}",  "phSum/p_{t}^{raw}",  "phHardSum/p_{t}^{raw}", 
  "neMax/p_{t}^{raw}",  "neSum/p_{t}^{raw}", 
  "eMax/p_{t}^{raw}",   "eSum/p_{t}^{raw}", 
  "muMax/p_{t}^{raw}",  "muSum/p_{t}^{raw}", 
  "N_{neu}^{jet}",   "N_{chg}^{jet}",   "N_{sum}^{jet}",
  "A_{j}",     "x_{j}" , "#Delta #phi", 
  "leadJet P_{t} (GeV)", "subleadJet P_{t} (GeV)"
};
//const int N_vars_xAx_Titles = sizeof(var_xAx_Titles)/sizeof(std::string);

const int var_xAx_reBin[]={                                   
  5, 5,     //reco/rawjtpt //500,0,500	   
  2,  2, //2,    //jteta, phi, jty   //60 ,-3,+3  and  //100,-4,+4  
  5, 5, 5,  //trks	   //200,0,2	   				   
  5, 5, 5,  //ch  
  5, 5, 5,  //ph					           
  5, 5,     //ne		
  5, 5,     //e					      
  5, 5,     //mu					      
  2, 2, 2,     // N Const.
  2, 2,     //xj,Aj	//100,0,1				     
  1,        //dphi	//50 ,0,+4		   				      
  5, 5      //sub/lead recojttpt, //500,0,500				      
};

const float jetQAxmax[]={
  (500. ), (500.),
  (5.   ), (3.5  ), //(5.  ),
  (2.00 ), (2.00), (2.00),
  (2.00 ), (2.00), (2.00),
  (2.00 ), (2.00), (2.00),
  (2.00 ), (2.00), 
  (2.00 ), (2.00), 
  (2.00 ), (2.00), 
  (30),   (30),   (30),
  (1.   ), (1.  ), 
  (3.2  ),
  (500. ), (500.)
};


const float jetQAxmin[]={
  30., 30.,  //  0.  ,  0.,
  -5. , -3.5, //-5.,
  0.  ,  0., 0.,
  0.  ,  0., 0.,
  0.  ,  0., 0.,
  0.  ,  0., 
  0.  ,  0., 
  0.  ,  0., 
  0.  ,  0., 0.,
  0.  ,  0., 
  2.0 ,
  0. , 0.
};

//const int jetTrigQABinning=10;
const int jetSpectraRapBinning=20;


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
