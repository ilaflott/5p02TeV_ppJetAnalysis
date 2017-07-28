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

//const std::string inputDir=CMSSW_BASE+"readForests/outputCondor/";
//const std::string inputDir=SCRATCH_BASE+"readForests/6.07.17_outputCondor/";
const std::string inputDir=SCRATCH_BASE+"readForests/7.25.17_outputCondor/";
const std::string outputDir=CMSSW_BASE+"doAnalysis/printPlots/output/";


// variable names for QA Plots
const std::string var[]={
  //jets, 0-3=4 vars
  "jtpt"  , "rawpt" ,
  "jteta" , "jtphi" ,

  //jet constituents, 4-30=27 vars
  "trkN", "trkSum", "trkMax", //trk
  "trkHardN", "trkHardSum",

  "phN", "phSum", "phMax",    //ph
  "phHardN", "phHardSum",

  "chN", "chSum", "chMax",    //ch
  "chHardN", "chHardSum",

  "neN" ,  "neSum" , "neMax" , //ne
  "eN"  ,  "eSum"  , "eMax"  ,    //e
  "muN" ,  "muSum" , "muMax" , //mu

  "neuMult", "chMult", "numConst",

  // dijets, 31-35=5 vars
  "Aj" , "xj" , "dphi",    //JERS might not like me
  "leadJetPt", "subleadJetPt"  
};
const int N_vars=sizeof(var)/sizeof(std::string);
const int jetQA_varStart=0;
const int jetConstits_varStart=4;
const int dijet_varStart=31;

// variable names, string array
const std::string var_xAx_Titles[] = {   
  "p_{t}^{reco} (GeV)" ,  "p_{t}^{raw} (GeV)",  
  "#eta_{jet}", "#phi_{jet}", //"y_{jet}" ,

  "trkN", "trkSum/p_{t}^{raw}", "trkMax/p_{t}^{raw}", 
  "trkHardN", "trkHardSum/p_{t}^{raw}",

  "chN", "chSum/p_{t}^{raw}", "chMax/p_{t}^{raw}", 
  "chHardN", "chHardSum/p_{t}^{raw}",

  "phN", "phSum/p_{t}^{raw}", "phMax/p_{t}^{raw}", 
  "phHardN", "phHardSum/p_{t}^{raw}",

  "neN", "neSum/p_{t}^{raw}",  "neMax/p_{t}^{raw}", 
  "eN",   "eSum/p_{t}^{raw}",   "eMax/p_{t}^{raw}", 
  "muN", "muSum/p_{t}^{raw}",  "muMax/p_{t}^{raw}",

  "N_{neu}^{jet}",   "N_{chg}^{jet}",   "N_{sum}^{jet}",

  "A_{j}",     "x_{j}" , "#Delta #phi", 
  "leadJet P_{t} (GeV)", "subleadJet P_{t} (GeV)"
};
//const int N_vars_xAx_Titles = sizeof(var_xAx_Titles)/sizeof(std::string);

const int var_xAx_reBin[]={                                   
  10, 10,     //reco/rawjtpt  
  2,  2,      //jteta, phi    

  1, 5, 5,    //trks	      
  1, 5,                       

  1, 5, 5,    //ph	      				           
  1, 5,                       

  1, 5, 5,    //ch  	      
  1, 5,                       

  1, 5, 5,    //ne	      	
  1, 5, 5,    //e	      				      
  1, 5, 5,    //mu	      				      

  1, 1, 1,    // N Const.     

  2, 2, 1,    //xj,Aj,dphi    
  20, 20        //sub/lead jttpt
};

const float jetQAxmax[]={
  (1000. ), (1000.),             //reco/rawjtpt  
  (5.   ), (3.3  ),  		 //jteta, phi    

  (20 ), (2.00), (2.00), 	 //trks	      
  (20 ), (2.00), 		                 

  (20 ), (2.00), (2.00), 	 //ph	      
  (20 ), (2.00), 		                 

  (20 ), (2.00), (2.00),	 //ch  	      
  (20 ), (2.00), 		              
   
  (20 ), (2.00), (2.00),	 //ne	      
  (20 ), (2.00), (2.00),	 //e	      
  (20 ), (2.00), (2.00),	 //mu	      

  (40),   (40),   (50),   	 // N Const.     

  (1.   ), (1.  ),   (3.2  ), 	 //xj,Aj,dphi    
  (1000. ), (1000.) 		 //sub/lead jttpt
};


const float jetQAxmin[]={
  (50. ), (50.  ) ,           //reco/rawjtpt  
  (-5. ), (-3.3 ) , 	      //jteta, phi    

  (0.  ), (0.   ) , (0.),      //trks	      
  (0.  ), (0.   ) , 

  (0.  ), (0.   ) , (0.),      //ph	      
  (0.  ), (0.   ) , 

  (0.  ), (0.   ) , (0.),      //ch  	      
  (0.  ), (0.   ) , 

  (0.  ), (0.   ) , (0.),      //ne	      
  (0.  ), (0.   ) , (0.),      //e	      
  (0.  ), (0.   ) , (0.),      //mu	      

  (0.  ), (0.   ) , (0.),      // N Const.     

  (0.  ), (0.   ) , (0.),      //xj,Aj,dphi    
  (50. ), (50.  )          //sub/lead jttpt
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



