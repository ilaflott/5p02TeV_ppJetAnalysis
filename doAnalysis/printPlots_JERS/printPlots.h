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
#include <iomanip>
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
  "/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/";
const std::string SCRATCH_BASE=
  //"/cms/heavyion/obaron/5p02TeV_ppJetAnalysis/";
  "/cms/heavyion/ilaflott/T2_US_MIT_SCRATCH/5p02TeV_ppJetAnalysis/readForests/";
const std::string HOME_BASE=
  "/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots_JERS/";

  //  "/Users/ilaflott/Working/Dropbox/localAnalysis/src/5p02TeV_ppJetAnalysis-master/doAnalysis/printPlots/";

const std::string inputDir=SCRATCH_BASE;//+"8.16.17_outputCondor/";
//const std::string inputDir=HOME_BASE+"outputCondor/8.16.17_outputCondor/";

const std::string outputDir=HOME_BASE+"output/";

// misc text ------------------------
const std::string crossSectionAxTitle  ="#sigma (#mub)";
const std::string ddcrossSectionAxTitle="#frac{d^{2}#sigma}{dp_{T}d#||{y}} (#mub/GeV)";
const std::string dcrossSectionAxTitle="#frac{d#sigma^{}}{dp_{T}} (#mub/GeV)";
const std::string AUAxTitle  ="A.U.";
const std::string ratioTitle ="MC/Data";
const std::string jetIDratioTitle ="JetID/noJetID";
const std::string intLumiString="L_{int} = 27.4  ^{}pb^{ -1}";
const std::string sqrts2k15ppString="2015 pp  ^{}#sqrt{s} = 5.02 TeV";
const std::string Py8TuneString="Py8 Tune CUETP8M1";

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
  "jet p_{t} (GeV)" ,  "jet p_{t}^{raw} (GeV)",  
  "#eta_{jet}", "#phi_{jet}", //"y_{jet}" ,

  "trkN", "trkSum/p_{t}^{raw}", "trkMax/p_{t}^{raw}", 
  "trkHardN", "trkHardSum/p_{t}^{raw}",

  "phN", "phSum/p_{t}^{raw}", "phMax/p_{t}^{raw}", 
  "phHardN", "phHardSum/p_{t}^{raw}",

  "chN", "chSum/p_{t}^{raw}", "chMax/p_{t}^{raw}", 
  "chHardN", "chHardSum/p_{t}^{raw}",

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
  1,  2,      //jteta, phi    

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
  
  3, 3, 1,    //xj,Aj,dphi    
  20, 20        //sub/lead jttpt
};

const float jetQAxmax[]={
  (1000. ), (1000.),             //reco/rawjtpt  
  (5.   ), (3.2  ),  		 //jteta, phi    

  (60 ), (2.00), (2.00), 	 //trks	      
  (60 ), (2.00), 		                 

  (60 ), (2.00), (2.00), 	 //ph	      
  (60 ), (2.00), 		                 

  (60 ), (2.00), (2.00),	 //ch  	      
  (60 ), (2.00), 		              
   
  (20 ), (2.00), (2.00),	 //ne	      
  (10 ), (2.00), (2.00),	 //e	      
  (10 ), (2.00), (2.00),	 //mu	      

  (100),   (100),   (100),   	 // N Const.     

  (1.   ), (1.  ),   (3.2  ), 	 //xj,Aj,dphi    
  (1000. ), (1000.) 		 //sub/lead jttpt
};


const float jetQAxmin[]={
  (50. ), (50.  ) ,           //reco/rawjtpt  
  (-5. ), (-3.2 ) , 	      //jteta, phi    

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

  (0.  ), (0.   ) , (2.0),      //xj,Aj,dphi    
  (50. ), (50.  )          //sub/lead jttpt
};		



//const int jetTrigQABinning=10;
const int jetSpectraRapBinning=20;


//// HELPER FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------

//void divideBinWidth(TH1 *h){
//  h->Sumw2();
//  for (int i=0;i<=h->GetNbinsX();++i){//binsX loop 
//    Float_t val = h->GetBinContent(i);
//    Float_t valErr = h->GetBinError(i);
//    val/=h->GetBinWidth(i);
//    valErr/=h->GetBinWidth(i);
//    h->SetBinContent(i,val);
//    h->SetBinError(i,valErr);
//  }//end nbinsX loop
//  h->GetXaxis()->CenterTitle();
//  h->GetYaxis()->CenterTitle();
//  return;
//}



void divideBinWidth(TH1 *h){
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in divideBinWidth"<<std::endl<<std::endl;

  //h->Sumw2();
  if(funcDebug)h->Print("base");
  double nEntries=h->GetEntries();
  int numbins=h->GetNbinsX();

  //std::cout<<"#of bin should be = "<<nbins_pt<<std::endl;
  if(funcDebug)std::cout<<"#of bins = "<<numbins<<std::endl;
  if(funcDebug)std::cout<<"integral before = "<<h->Integral()<<std::endl;
  if(funcDebug)std::cout<<"entries before = "<<nEntries<<std::endl;

  for (int i=1;i<=numbins;i++)
    {

      if(funcDebug)std::cout<<std::endl<<"bin-i = "<< i <<std::endl;
      if(funcDebug)std::cout<<"bin-i width = "<<  h->GetBinWidth(i)  << std::endl;
      if(funcDebug)std::cout<<"bin-i pt range is "<<   h->GetBinLowEdge(i) << " to " << (h->GetBinLowEdge(i)+h->GetBinWidth(i)) << std::endl;
      if(funcDebug)std::cout<<std::endl<<"old val = "   <<      h->GetBinContent(i) << std::endl;
      if(funcDebug)std::cout<<"old valErr = "<<      h->GetBinError(i)   << std::endl;


      Float_t val    = h->GetBinContent(i);
      Float_t valErr = h->GetBinError(i);

      val    /= h->GetBinWidth(i);
      valErr /= h->GetBinWidth(i);

      h->SetBinContent(i,val);
      h->SetBinError(i,valErr);

      if(funcDebug)std::cout<<"new val = "   <<      h->GetBinContent(i) << std::endl;
      if(funcDebug)std::cout<<"new valErr = "<<      h->GetBinError(i)   << std::endl;

    }

  if(funcDebug)std::cout<<"integral after = "<<h->Integral()<<std::endl;
  if(funcDebug)std::cout<<"entries after = "<<h->GetEntries()<<std::endl;

  //h->SetEntries(nEntries);

  //std::cout<<"exiting function"<<std::endl;
  if(funcDebug)std::cout<<std::endl<<"divideBinWidth done"<<std::endl<<std::endl;
  return;
}

float getIntegralJetSpectra(TH1* h){
  bool funcDebug=false;

  float integral=0.;
  int numbins=h->GetNbinsX();
  std::cout<<"numbins="<<numbins<<std::endl;
  for(int i=1;i<=numbins;i++){
    
    Float_t val    = h->GetBinContent(i);
    if(funcDebug)std::cout<<"bin i="<<i<<std::endl;
    if(funcDebug)std::cout<<"bin width = "<<h->GetBinWidth(i)<<std::endl;
    if(funcDebug)std::cout<<"pt range of "<<h->GetBinLowEdge(i)<<" to "<<h->GetBinLowEdge(i)+h->GetBinWidth(i)<<std::endl;
    if(funcDebug)std::cout<<"val = "<<val<<std::endl;
    val    *= h->GetBinWidth(i);
    integral+=val;
    if(funcDebug)std::cout<<"integral now = "<<integral<<std::endl<<std::endl;
  }

  if(funcDebug)std::cout<<"final integral is:" << integral<<std::endl;
  if(funcDebug)std::cout<<"normal integral says:"<<h->Integral()<<std::endl;
  
  return integral;
}


float getIntegralErrJetSpectra(TH1* h){
  
  float integralErr=0.;
  
  int numbins=h->GetNbinsX();
  
  for(int i=1;i<=numbins;i++){
    
    Float_t valErr = h->GetBinError(i);
    
    valErr *= h->GetBinWidth(i);
    
    integralErr+=valErr*valErr;
  }
  
  integralErr=std::sqrt(integralErr);
  
  return integralErr;
}


//void makeJERSHists(TCanvas* pdfOutCanv=NULL, std::string thePDFFileName, 
void makeJERSHists(TCanvas* pdfOutCanv, //std::string thePDFFileName, 
		   TH1F* hrsp, std::string hrspTitle, std::string hrsp_XAxTitle,
		   TFile* fout=NULL, TDirectory* tdir=NULL){
  bool funcDebug=true;
  pdfOutCanv->cd();
  
  hrsp->SetTitle(hrspTitle.c_str());    
  //hrsp->SetTitleSize(0.02);    
  
  hrsp->SetMarkerStyle(8);    
  hrsp->SetMarkerSize(1.0);    
  
  hrsp->SetAxisRange(0.5, 1.5, "X");
  hrsp->GetXaxis()->SetTitle( hrsp_XAxTitle.c_str() );    
  
  hrsp->SetAxisRange(0.00001, hrsp->GetMaximum()*2., "Y");
  //hrsp->SetAxisRange(0.0, hrsp->GetMaximum()*1.1, "Y");
  hrsp->GetYaxis()->SetTitle("A.U.");    
  
  hrsp->Draw("E1");    

  //if((bool)pdfOutCanv)pdfOutCanv->Print(thePDFFileName.c_str());
  
  if((bool)fout){
    
    fout->cd();
    if((bool)tdir)tdir->cd();
    
    TCanvas* canv_fout=(TCanvas*)pdfOutCanv->DrawClone();
    if(funcDebug)std::cout<<"writing canvas to file"<<std::endl;
    canv_fout->SetTitle(hrspTitle.c_str());
    canv_fout->Write(  (((std::string)hrsp->GetName()) + "_canv").c_str());
  }
  
  return;
}

std::string float2string_prec(float num, int precision){
  std::stringstream stream;
  stream << std::fixed << std::setprecision(precision) << num;
  std::string floatstring = stream.str();
  return floatstring;
}

void makeFitProbStr(TF1* fit, std::string* prob_str){
  if(!fit)return;
  if(!prob_str)return;  
  //double prob=fit->GetProb();
  *(prob_str)=float2string_prec(fit->GetProb(), 2);    
  return;
}
void makeFitChi2Str(TF1* fit, std::string* chi2_str){
  if(!fit)return;
  if(!chi2_str)return;  
  *(chi2_str)=float2string_prec(fit->GetChisquare(), 2);    
  return;
}
