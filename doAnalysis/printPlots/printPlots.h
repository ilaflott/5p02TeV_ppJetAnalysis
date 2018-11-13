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
#include "printPlots_ppMC.h"
#include "printPlots_ppData.h"
#include "printPlots_binsNconsts.h"
#include "printPlots_strings.h"
#include "printPlots_draw.h"

// better way of organizing? not done yet; still kind of want to 4/20/18
//#include "printPlots_evtQAPlots.h"
//#include "printPlots_jetQAPlots.h"
//#include "printPlots_jetTrigPlots.h"
//#include "printPlots_jetIDPlots.h"
//#include "printPlots_jetMultPlots.h"
//#include "printPlots_fitsNfuncs.h"


//// HELPER FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------
bool str_replace(std::string& str, const std::string& from, const std::string& to) {
//bool str_replace(std::string str, const std::string from, const std::string to) {
  size_t start_pos = str.find(from);
  if(start_pos == std::string::npos)
    return false;
  str.replace(start_pos, from.length(), to);
  return true;
}

bool findsubstr(std::string& str, const std::string& substr) {
  //bool str_replace(std::string str, const std::string from, const std::string to) {
  bool substrfound=false;
  size_t substr_pos= str.find(substr);
  //size_t start_pos = str.find(from);
  if(substr_pos != std::string::npos)
    substrfound=true;  
  return substrfound;
}


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








///////////////////////////// EVT QA: COUNTS, EFF LUMI, VTX, NREF, MC PTHAT
void getEventCounts(TFile* fin, bool isData){
  bool funcDebug=false;
  std::cout<<" grabbing raw event counts "<<std::endl;   
  
  if(!fin){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  
  if(funcDebug)std::cout<<"opening hists"<<std::endl;
  TH1F *h_NEvents         = (TH1F*)fin->Get("NEvents");
  TH1F *h_NEvents_read    = (TH1F*)fin->Get("NEvents_read");      
  TH1F *h_NEvents_skimCut = (TH1F*)fin->Get("NEvents_skimCut");
  TH1F *h_NEvents_vzCut   = (TH1F*)fin->Get("NEvents_vzCut");
  
//  if(!h_NEvents        ||
//     !h_NEvents_read   ||
//     !h_NEvents_skimCut||
//     !h_NEvents_vzCut   ) assert(false);
  if(!h_NEvents        ) assert(false);
  if(!h_NEvents_read   ) assert(false);
  if(!h_NEvents_skimCut) assert(false);
  if(!h_NEvents_vzCut  ) assert(false);
  
  TH1F *h_NEvents_skipped = NULL;    
  if(isData){
    h_NEvents_skipped=(TH1F*)fin->Get("NEvents_skipped");      
    std::cout<<std::endl<<"/// Data Evt-Count Summary ///"<<std::endl<<std::endl;      
  }
  else std::cout<<std::endl<<"/// MC Evt-Count Summary ///"<<std::endl<<std::endl;      
  
  std::cout<<"Total Num of Events in file(s) opened       = " <<
    h_NEvents->GetEntries()<<std::endl;
  if(isData)  std::cout<<"Total Num of Events skipped from those file(s) = " <<
		h_NEvents_skipped->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read from those file(s) = " <<
    h_NEvents_read->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read passing skimCuts   = " <<
    h_NEvents_skimCut->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read passing vzCuts and skimCuts    = " <<
    h_NEvents_vzCut->GetEntries()<<std::endl<<std::endl;
  

//  TH1F *h_NEvents_trigd=NULL;
//  TH1F *h_NEvents_jet100=NULL,*h_NEvents_jet80=NULL,*h_NEvents_jet60=NULL,*h_NEvents_jet40=NULL;  
//  if(isData){
//    
//    h_NEvents_trigd=(TH1F*)fin->Get("NEvents_trigd");      
//
//    h_NEvents_jet40=(TH1F*)fin->Get("NEvents_jet40Trigd");      
//    h_NEvents_jet60=(TH1F*)fin->Get("NEvents_jet60Trigd");      
//    h_NEvents_jet80=(TH1F*)fin->Get("NEvents_jet80Trigd");      
//    h_NEvents_jet100=(TH1F*)fin->Get("NEvents_jet100Trigd");      
//    
//
//    std::cout<<"Total Num of Events trigd in those file(s) = " <<
//      h_NEvents_trigd->GetEntries()<<std::endl;
//   a ssert(false);      
//    std::cout<<"Total Num of Events jet40 in those file(s) = " <<
//      h_NEvents_jet40->GetEntries()<<std::endl;
//    std::cout<<"Total Num of Events jet60 in those file(s) = " <<
//      h_NEvents_jet60->GetEntries()<<std::endl;
//    std::cout<<"Total Num of Events jet80 in those file(s) = " <<
//      h_NEvents_jet80->GetEntries()<<std::endl;
//
//    if(h_NEvents_jet100)//not always made anymore
//      std::cout<<"Total Num of Events jet100 in those file(s) = " <<
//	h_NEvents_jet100->GetEntries()<<std::endl;
//  }
  
  return;

}



///////////////////////////// EVT QA: COUNTS, EFF LUMI, VTX, NREF, MC PTHAT
void getEventCountsTrkMult(TFile* fin, bool isData){
  bool funcDebug=false;
  std::cout<<" grabbing raw event counts "<<std::endl;   
  
  
  if(!fin){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  
  if(funcDebug)std::cout<<"opening hists"<<std::endl;
  TH1F *h_NEvents         = (TH1F*)fin->Get("NEvents");
  TH1F *h_NEvents_read    = (TH1F*)fin->Get("NEvents_read");      
  TH1F *h_NEvents_skimCut = (TH1F*)fin->Get("NEvents_skimCut");
  TH1F *h_NEvents_vzCut   = (TH1F*)fin->Get("NEvents_vzCut");
  TH1F *h_NEvents_trkCuts_1   = (TH1F*)fin->Get("NEvents_trkCuts_1");
  TH1F *h_NEvents_trkCuts_2   = (TH1F*)fin->Get("NEvents_trkCuts_2");
  TH1F *h_NEvents_trkCuts_3   = (TH1F*)fin->Get("NEvents_trkCuts_3");
  
  TH1F *h_NEvents_skipped = NULL;
  TH1F *h_NEvents_trigd=NULL;
  TH1F *h_NEvents_jet100=NULL,*h_NEvents_jet80=NULL,*h_NEvents_jet60=NULL,*h_NEvents_jet40=NULL;
  
  if(isData)
    { h_NEvents_skipped=(TH1F*)fin->Get("NEvents_skipped");      
      h_NEvents_trigd=(TH1F*)fin->Get("NEvents_trigd");      
      h_NEvents_jet40=(TH1F*)fin->Get("NEvents_jet40Trigd");      
      h_NEvents_jet60=(TH1F*)fin->Get("NEvents_jet60Trigd");      
      h_NEvents_jet80=(TH1F*)fin->Get("NEvents_jet80Trigd");      
      h_NEvents_jet100=(TH1F*)fin->Get("NEvents_jet100Trigd");      
    }
  
  if(isData)std::cout<<std::endl<<"/// Data Evt-Count Summary ///"<<std::endl<<std::endl;      
  else std::cout<<std::endl<<"/// MC Evt-Count Summary ///"<<std::endl<<std::endl;      
  
  std::cout<<"Total Num of Events in file(s) opened       = " <<
    h_NEvents->GetEntries()<<std::endl;
  if(isData)
    std::cout<<"Total Num of Events skipped from those file(s) = " <<
      h_NEvents_skipped->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read from those file(s) = " <<
	h_NEvents_read->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read passing skimCuts   = " <<
    h_NEvents_skimCut->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read passing vzCuts and skimCuts    = " <<
    h_NEvents_vzCut->GetEntries()<<std::endl<<std::endl;
  
  if(isData){
    std::cout<<"Total Num of Events trigd in those file(s) = " <<
      h_NEvents_trigd->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events jet40 in those file(s) = " <<
      h_NEvents_jet40->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events jet60 in those file(s) = " <<
      h_NEvents_jet60->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events jet80 in those file(s) = " <<
      h_NEvents_jet80->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events jet100 in those file(s) = " <<
      h_NEvents_jet100->GetEntries()<<std::endl;}
    
  std::cout<<"Total Num of Events read passing trkHiPur Cut (#1) = " <<
    h_NEvents_trkCuts_1->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read passing vtx spacing Cut (#2) = " <<
    h_NEvents_trkCuts_2->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read passing oneVtxInRange Cut (#3) = " <<
    h_NEvents_trkCuts_3->GetEntries()<<std::endl;
  
  return;

}



void printDataEventCountReport(TFile* fin){
  bool funcDebug=true;
  if(funcDebug)std::cout<<std::endl<<"printDataEventCountReport Called"<<std::endl<<std::endl;
  
  if(!fin){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  
   
  std::cout<<"--------------------------"<<std::endl;
  std::cout<<"--- Event Loop Summary ---"<<std::endl;
  std::cout<<"--------------------------"<<std::endl;
  std::cout<<"# evts (description) = / count / pscl weighted count" <<std::endl<<std::endl;

  TH1F *h_NEvents         = (TH1F*)fin->Get("NEvents");
  TH1F *h_NEvents_trigd_1=(TH1F*)fin->Get("NEvents_trigd_1");
  TH1F *h_NEvents_nopscl =(TH1F*)fin->Get("NEvents_nopscl");
  TH1F *h_NEvents_jet100=(TH1F*)fin->Get("NEvents_jet100Trigd"),    *h_NEvents_jet80=(TH1F*)fin->Get("NEvents_jet80Trigd"),    *h_NEvents_jet60=(TH1F*)fin->Get("NEvents_jet60Trigd"),    *h_NEvents_jet40=(TH1F*)fin->Get("NEvents_jet40Trigd");
  
  std::cout<<"# evts in file(s) opened    = / "    <<h_NEvents->GetBinContent(1)        << " / "<< h_NEvents->GetBinContent(2)<<std::endl; 
  std::cout<<"# evts trigd in those files = / "    <<h_NEvents_trigd_1->GetBinContent(1) << " / "<< h_NEvents_trigd_1->GetBinContent(2)<<std::endl;
  std::cout<<"# evts w pscl==1           = / " << h_NEvents_nopscl->GetBinContent(1) << " / " << h_NEvents_nopscl->GetBinContent(2) <<std::endl << std::endl;
  std::cout<<"# evts passing is40  = / " << h_NEvents_jet40->GetBinContent(1)  << " / " << h_NEvents_jet40->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing is60  = / " << h_NEvents_jet60->GetBinContent(1)  << " / " << h_NEvents_jet60->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing is80  = / " << h_NEvents_jet80->GetBinContent(1)  << " / " << h_NEvents_jet80->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing is100 = / " << h_NEvents_jet100->GetBinContent(1) << " / " << h_NEvents_jet100->GetBinContent(2) <<std::endl<<std::endl;  
  
  TH1F *h_NEvents_read    = (TH1F*)fin->Get("NEvents_read");        
  TH1F *h_NEvents_skipped = (TH1F*)fin->Get("NEvents_skipped");  
  TH1F *h_NEvents_trigd_2=(TH1F*)fin->Get("NEvents_trigd_2");
  
  std::cout<<"# evts in file(s) opened           = / "    <<h_NEvents->GetBinContent(1)        << " / "<< h_NEvents->GetBinContent(2)<<std::endl; //repeat on purpose
  std::cout<<"# evts read from those file(s)     = / " << h_NEvents_read->GetBinContent(1)    << " / " << h_NEvents_read->GetBinContent(2)<<std::endl;
  std::cout<<"# evts skipped from those file(s)  = / " << h_NEvents_skipped->GetBinContent(1) << " / " << h_NEvents_skipped->GetBinContent(2)<<std::endl;
  std::cout<<"# evts passing a trigger post-skip = / " << h_NEvents_trigd_2->GetBinContent(1) << " / " << h_NEvents_trigd_2->GetBinContent(2) <<std::endl << std::endl;
  
  TH1F *h_NEvents_skimCut = (TH1F*)fin->Get("NEvents_skimCut");
  TH1F *h_NEvents_vzCut   = (TH1F*)fin->Get("NEvents_vzCut");
  
  std::cout<<"# evts passing skimCuts  = / " << h_NEvents_skimCut->GetBinContent(1) << " / " << h_NEvents_skimCut->GetBinContent(2) << std::endl;
  std::cout<<"# evts passing vzCuts    = / " << h_NEvents_vzCut->GetBinContent(1)   << " / " << h_NEvents_vzCut->GetBinContent(2)   << std::endl<<std::endl;
  
  
  TH1F *hpp_HLT40trgEta=(TH1F*)fin->Get("isHLT40_trgEta"), *hpp_HLT60trgEta=(TH1F*)fin->Get("isHLT60_trgEta"), *hpp_HLT80trgEta=(TH1F*)fin->Get("isHLT80_trgEta"),   *hpp_HLT100trgEta=(TH1F*)fin->Get("isHLT100_trgEta"),*hpp_HLTCombtrgPt=(TH1F*)fin->Get("HLTComb_trgPt") ;
  TH1F *hpp_IncHLT40trgEta=(TH1F*)fin->Get("IncHLT40_trgEta"), *hpp_IncHLT60trgEta=(TH1F*)fin->Get("IncHLT60_trgEta"), *hpp_IncHLT80trgEta=(TH1F*)fin->Get("IncHLT80_trgEta"),   *hpp_IncHLT100trgEta=(TH1F*)fin->Get("IncHLT100_trgEta");
  
  std::cout<<"# evts passing trig  = / " << hpp_HLTCombtrgPt ->GetEntries() << " / " <<  hpp_HLTCombtrgPt->Integral()<<std::endl;
  std::cout<<"# evts passing is40  = / " << hpp_HLT40trgEta ->GetEntries()  << " / " <<  hpp_HLT40trgEta->Integral()<<std::endl;
  std::cout<<"# evts passing is60  = / " << hpp_HLT60trgEta ->GetEntries()  << " / " <<  hpp_HLT60trgEta->Integral()<<std::endl;
  std::cout<<"# evts passing is80  = / " << hpp_HLT80trgEta ->GetEntries()  << " / " <<  hpp_HLT80trgEta->Integral()<<std::endl;
  std::cout<<"# evts passing is100 = / " << hpp_HLT100trgEta->GetEntries()  << " / " <<  hpp_HLT100trgEta->Integral()<<std::endl<<std::endl;  
  
  std::cout<<"# evts passing each trigger with no exclusion" <<std::endl;
  std::cout<<"# evts passing HLT40  = / " << hpp_IncHLT40trgEta ->GetEntries()  << " / " <<  hpp_IncHLT40trgEta ->Integral()<<std::endl;
  std::cout<<"# evts passing HLT60  = / " << hpp_IncHLT60trgEta ->GetEntries()  << " / " <<  hpp_IncHLT60trgEta ->Integral()<<std::endl;
  std::cout<<"# evts passing HLT80  = / " << hpp_IncHLT80trgEta ->GetEntries()  << " / " <<  hpp_IncHLT80trgEta ->Integral()<<std::endl;
  std::cout<<"# evts passing HLT100 = / " << hpp_IncHLT100trgEta->GetEntries()  << " / " <<  hpp_IncHLT100trgEta->Integral()<<std::endl<<std::endl;  
  
  
  std::cout<<"------------------------------"<<std::endl;
  std::cout<<"--- Jet-Event Loop Summary ---"<<std::endl;
  std::cout<<"------------------------------"<<std::endl;
  std::cout<<"# evts (desc) = / count / pscl weighted count" <<std::endl<<std::endl;
  
  TH1F* h_NEvents_withJets           =(TH1F*)fin->Get("NEvents_withJets");	    
  TH1F* h_NEvents_withJets_jtptCut   =(TH1F*)fin->Get("NEvents_withJets_jtptCut");  
  TH1F* h_NEvents_withJets_jtetaCut1 =(TH1F*)fin->Get("NEvents_withJets_jtetaCut1");
  TH1F* h_NEvents_withJets_jtetaCut2 =(TH1F*)fin->Get("NEvents_withJets_jtetaCut2");
  std::cout<<"# evts w min. 1 jet                   = / " << h_NEvents_withJets->GetBinContent(1)           << " / " << h_NEvents_withJets->GetBinContent(2) <<std::endl; 
  std::cout<<"# evts w min. 1 jet w jtpt > jtPtCut = / " << h_NEvents_withJets_jtptCut->GetBinContent(1)   << " / " << h_NEvents_withJets_jtptCut->GetBinContent(2) << std::endl;
  std::cout<<"# evts w min. 1 jet w |jteta| < 4.7  = / " << h_NEvents_withJets_jtetaCut1->GetBinContent(1) << " / " << h_NEvents_withJets_jtetaCut1->GetBinContent(2) << std::endl; 
  std::cout<<"# evts w min. 1 jet w jtEtaCut(s)    = / " << h_NEvents_withJets_jtetaCut2->GetBinContent(1) << " / " << h_NEvents_withJets_jtetaCut2->GetBinContent(2) << std::endl;  
  
  TH1F* h_NEvents_withJets_JetIDCut = (TH1F*)fin->Get("NEvents_withJets_JetIDCut");
  if((bool)h_NEvents_withJets_JetIDCut)
    std::cout<<"# evts w min. 1 jet w JetIDCut       = / " << h_NEvents_withJets_JetIDCut->GetBinContent(1) << " / " << h_NEvents_withJets_JetIDCut->GetBinContent(2) << std::endl;          
  std::cout<<std::endl;  
  
  TH1F* h_NJets           =(TH1F*)fin->Get("NJets");	    
  TH1F* h_NJets_jtptCut   =(TH1F*)fin->Get("NJets_jtptCut");  
  TH1F* h_NJets_jtetaCut1 =(TH1F*)fin->Get("NJets_jtetaCut1");
  TH1F* h_NJets_jtetaCut2 =(TH1F*)fin->Get("NJets_jtetaCut2");
  std::cout<<"# jets read                   = / " << h_NJets->GetBinContent(1)           << " / " << h_NJets->GetBinContent(2) <<std::endl; 
  std::cout<<"# jets w jtpt > jtptCut      = / " << h_NJets_jtptCut->GetBinContent(1)   << " / " << h_NJets_jtptCut->GetBinContent(2) <<std::endl;   
  std::cout<<"# jets w |jteta|<4.7         = / " << h_NJets_jtetaCut1->GetBinContent(1) << " / " << h_NJets_jtetaCut1->GetBinContent(2) <<std::endl; 
  std::cout<<"# jets w passing jtEtaCut(s) = / " << h_NJets_jtetaCut2->GetBinContent(1) << " / " << h_NJets_jtetaCut2->GetBinContent(2) <<std::endl; 

  TH1F* h_NJets_JetIDCut =(TH1F*)fin->Get("NJets_JetIDCut");
  if((bool)h_NJets_JetIDCut)
    std::cout<<"# jets w jetIDCut            = / " << h_NJets_JetIDCut->GetBinContent(1) << " / " << h_NJets_JetIDCut->GetBinContent(2) <<std::endl; 
  

  /*
    if(isData)
    { h_NEvents_skipped=(TH1F*)fin->Get("NEvents_skipped");      
    h_NEvents_trigd=(TH1F*)fin->Get("NEvents_trigd");      
    h_NEvents_jet40=(TH1F*)fin->Get("NEvents_jet40Trigd");      
    h_NEvents_jet60=(TH1F*)fin->Get("NEvents_jet60Trigd");      
    h_NEvents_jet80=(TH1F*)fin->Get("NEvents_jet80Trigd");      
    h_NEvents_jet100=(TH1F*)fin->Get("NEvents_jet100Trigd");      
    }
    
    if(isData)std::cout<<std::endl<<"/// Data Evt-Count Summary ///"<<std::endl<<std::endl;      
    else std::cout<<std::endl<<"/// MC Evt-Count Summary ///"<<std::endl<<std::endl;      
    
    std::cout<<"Total Num of Events in file(s) opened       = " <<
    h_NEvents->GetEntries()<<std::endl;
    if(isData)
    std::cout<<"Total Num of Events skipped from those file(s) = " <<
    h_NEvents_skipped->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events read from those file(s) = " <<
    h_NEvents_read->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events read passing skimCuts   = " <<
    h_NEvents_skimCut->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events read passing vzCuts and skimCuts    = " <<
    h_NEvents_vzCut->GetEntries()<<std::endl<<std::endl;
    
    if(isData){
    std::cout<<"Total Num of Events trigd in those file(s) = " <<
    h_NEvents_trigd->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events jet40 in those file(s) = " <<
    h_NEvents_jet40->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events jet60 in those file(s) = " <<
    h_NEvents_jet60->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events jet80 in those file(s) = " <<
    h_NEvents_jet80->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events jet100 in those file(s) = " <<
    h_NEvents_jet100->GetEntries()<<std::endl;}
  */  

    
  return;

}

void printJetCountReport(TFile* fin, bool isData){
  //void printJetCountReport(TFile* fin){
  if(isData){
    std::cout<<std::endl;
    std::cout<<"------------- DATA -----------"<<std::endl;
    std::cout<<"--- Jet-Event Loop Summary ---"<<std::endl;
    std::cout<<"------------------------------"<<std::endl;
    std::cout<<std::endl;
  }
  else {
    std::cout<<std::endl;
    std::cout<<"-------------  MC  -----------"<<std::endl;
    std::cout<<"--- Jet-Event Loop Summary ---"<<std::endl;
    std::cout<<"------------------------------"<<std::endl;
    std::cout<<std::endl;
  }
  TH1F* h_NJets           =(TH1F*)fin->Get("NJets");	    
  std::cout<<"# jets read                   = / " << h_NJets->GetBinContent(1)           << " / " << h_NJets->GetBinContent(2) <<std::endl; 
  if(isData){
    TH1F* h_NJets_jtptCut   =(TH1F*)fin->Get("NJets_jtptCut");  
    TH1F* h_NJets_jtetaCut1 =(TH1F*)fin->Get("NJets_jtetaCut1");
    TH1F* h_NJets_jtetaCut2 =(TH1F*)fin->Get("NJets_jtetaCut2");
    
    std::cout<<"# jets w jtpt > jtptCut       = / " << h_NJets_jtptCut->GetBinContent(1)   << " / " << h_NJets_jtptCut->GetBinContent(2) <<std::endl;   
    std::cout<<"# jets w |jteta|<4.7          = / " << h_NJets_jtetaCut1->GetBinContent(1) << " / " << h_NJets_jtetaCut1->GetBinContent(2) <<std::endl; 
    std::cout<<"# jets w passing jtEtaCut(s)  = / " << h_NJets_jtetaCut2->GetBinContent(1) << " / " << h_NJets_jtetaCut2->GetBinContent(2) <<std::endl; 
  }
  else { //isMC
    TH1F* h_NJets_kmatCut=(TH1F*)fin->Get("NJets_kmatCut");
    std::cout<<"# jets passing kinematic cuts = / " << h_NJets_kmatCut->GetBinContent(1)   << " / " << h_NJets_kmatCut->GetBinContent(2) <<std::endl;   
  }
  
  TH1F* h_NJets_JetIDCut =(TH1F*)fin->Get("NJets_JetIDCut");
  if((bool)h_NJets_JetIDCut)
    std::cout<<"# jets w jetIDCut             = / " << h_NJets_JetIDCut->GetBinContent(1) << " / " << h_NJets_JetIDCut->GetBinContent(2) <<std::endl; 
  
  return;
}

long double computeEffLumi(TFile* finData){
  bool funcDebug=false;
  if(funcDebug)std::cout<<"computeEffLumi called."<<std::endl;
  long double effIntgrtdLumi_vz=1., LumiEff_vz=1.;// eff-> effective, Eff->efficiency
  
  TH1F *h_NEvents_vzCut   = (TH1F*)finData->Get("NEvents_vzCut");
  TH1F *h_NEvents_read    = (TH1F*)finData->Get("NEvents_read"); 

  std::cout<<std::endl<<"integrated luminosity     ="<<intgrtdLumi<<std::endl;               
  
  //if(funcDebug)std::cout<<"total # evts post vzCut = " << h_NEvents_vzCut->GetEntries() << std::endl;
  //if(funcDebug)std::cout<<"total # evts post vzCut = " << h_NEvents_vzCut->GetBinContent(1) << std::endl;
  if(funcDebug)std::cout<<"total # (weighted) evts post vzCut = " << h_NEvents_vzCut->GetBinContent(2) << std::endl;     
  //if(funcDebug)std::cout<<"total # evts read = " << h_NEvents_read->GetEntries() << std::endl;
  //if(funcDebug)std::cout<<"total # evts read = " << h_NEvents_read->GetBinContent(1) << std::endl;
  if(funcDebug)std::cout<<"total # (weighted) evts read = " << h_NEvents_read->GetBinContent(2) << std::endl;     
  
  LumiEff_vz = h_NEvents_vzCut->GetBinContent(2)/h_NEvents_read->GetBinContent(2);                                    
  //LumiEff_vz = h_NEvents_vzCut->GetBinContent(1)/h_NEvents_read->GetBinContent(1);                                    
  std::cout<<"lumi efficiency, vz cuts             ="<<LumiEff_vz<<std::endl;                                   
  effIntgrtdLumi_vz=intgrtdLumi*LumiEff_vz;                                                                   
  std::cout<<"effective ingrtdLumi for dataset+cuts="<<effIntgrtdLumi_vz<<std::endl;                            
  
  return effIntgrtdLumi_vz;
}

long double computeEffLumiTrkMult(TFile* finData){

  long double effIntgrtdLumi_vz=1., LumiEff_vz=1.;// eff-> effective, Eff->efficiency
  
  std::cout<<std::endl<<"dataset integrated Luminosity (microbarns) ="<<intgrtdLumi<<std::endl;               
  TH1F *h_NEvents_vzCut   = (TH1F*)finData->Get("NEvents_trkCuts_3");
  TH1F *h_NEvents_read    = (TH1F*)finData->Get("NEvents_read");      
  
  LumiEff_vz = h_NEvents_vzCut->GetEntries()/h_NEvents_read->GetEntries();                                    
  effIntgrtdLumi_vz=intgrtdLumi*LumiEff_vz;                                                                   
  
  std::cout<<std::endl<<"lumi efficiency, vz cuts="<<LumiEff_vz<<std::endl;                                   
  std::cout<<"effective ingrtdLumi for dataset+cuts="<<effIntgrtdLumi_vz<<std::endl;                          
  //LumiEff_vz=0.999;//temp
  return effIntgrtdLumi_vz;
}


const int CANVX=1000;
const int CANVY=800;

void printEvtVtxQAHist( TFile* finData , std::string inDataHistName, 
			TFile* finMC   , std::string inMCHistName  ,
			int rebin, std::string thePDFFileName, long double theLumi , TFile* fout=NULL) {
  
  if(!finData || !finMC ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(funcDebug)std::cout<<"creating temporary canvas for printing Evt plots..."<<std::endl;
  TCanvas *temp_canvEvt = new TCanvas(("vtxCanv"+inDataHistName).c_str(), "temp Evt Canv withLog", CANVX,CANVY);
  temp_canvEvt->cd();
  
  TPad *evtpad1 = new TPad("evtpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  TPad *evtpad2 = new TPad("evtpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
  
  evtpad1->SetBottomMargin(0); 
  
  evtpad1->Draw();
  evtpad1->cd();
  
  // open hist
  std::cout<<" getting hist(s): "<<inDataHistName<<std::endl;
  
  TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inDataHistName.c_str() );
  if(!theDataEvtQAHist){ std::cout<<"input data hist not found! skipping hist"<<std::endl;
    return;}
  theDataEvtQAHist=(TH1F*)theDataEvtQAHist->TH1::Rebin( (int)rebin, (inDataHistName+"_rebin").c_str());
  theDataEvtQAHist->Scale( 1./theDataEvtQAHist->GetBinWidth(1) );
  theDataEvtQAHist->Scale( 1./theLumi );
  dataHistStyle(theDataEvtQAHist);
  
  
  TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inMCHistName.c_str() );
  if(!theMCEvtQAHist){ std::cout<<"input MC hist not found! skipping hist"<<std::endl;
    return;}
  theMCEvtQAHist=(TH1F*)theMCEvtQAHist->TH1::Rebin( (int)rebin, (inMCHistName+"_rebin").c_str());
  theMCEvtQAHist->Scale( 1./theMCEvtQAHist->GetBinWidth(1) );
  theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );	
  MCHistStyle(theMCEvtQAHist);
  
  


  // title+axes
  
  std::string h_Title     =PDStatsString_2;
  std::string h_XAx_Title;// ="v_{z}^{evt} (cm)";
  float vxyr_ymin=1.4e-02, vxyr_ymax=1.e+06;
  
  if( inDataHistName.find("Vz") != std::string::npos ){
    h_XAx_Title ="v_{z}^{evt} (cm)";}
  else if( inDataHistName.find("Vx") != std::string::npos ){
    evtpad1->SetLogy(1);
    h_XAx_Title ="v_{x}^{evt} (cm)";
    theDataEvtQAHist->SetAxisRange(vxyr_ymin,vxyr_ymax,"Y");
    theMCEvtQAHist  ->SetAxisRange(vxyr_ymin,vxyr_ymax,"Y");
  }
  else if( inDataHistName.find("Vy") != std::string::npos ){
    evtpad1->SetLogy(1);
    h_XAx_Title ="v_{y}^{evt} (cm)";
    theDataEvtQAHist->SetAxisRange(vxyr_ymin,vxyr_ymax,"Y");
    theMCEvtQAHist  ->SetAxisRange(vxyr_ymin,vxyr_ymax,"Y");
  }
  else if( inDataHistName.find("Vr") != std::string::npos ){
    evtpad1->SetLogy(1);
    h_XAx_Title ="v_{r}^{evt} (cm)";
    theDataEvtQAHist->SetAxisRange(vxyr_ymin,vxyr_ymax,"Y");
    theMCEvtQAHist  ->SetAxisRange(vxyr_ymin,vxyr_ymax,"Y");
  }
  std::string h_YAx_Title ="A.U.";
  



  
  Double_t dataMax = theDataEvtQAHist->GetMaximum(), MCMax=theMCEvtQAHist->GetMaximum();
  std::cout<<"dataMax="<<dataMax<<std::endl;
  std::cout<<"MCMax="<<    MCMax<<std::endl;
  if(dataMax>MCMax){
    std::cout<<"dataMax > MC Max"<<std::endl;
    theDataEvtQAHist->SetTitle (    h_Title.c_str() );
    theDataEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
    theDataEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
    
    theDataEvtQAHist->Draw("E"); 
    theMCEvtQAHist->Draw("HIST E SAME"); 
    theDataEvtQAHist->Draw("E SAME"); //because i want data points written over MC's red line hist
  }
  else {
    theMCEvtQAHist->SetTitle (    h_Title.c_str() );
    theMCEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
    theMCEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
    //theMCEvtQAHist->SetAxisRange(.001,10.*MCMax,"Y");
    theMCEvtQAHist->Draw("HIST E"); 
    theDataEvtQAHist->Draw("E SAME"); 
  }
  
  
  // legend
  float legx1=0.65, legx2=legx1+0.08;
  float legy1=0.74, legy2=legy1+0.09;

  TLegend* theEvtQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");
  theEvtQALeg->SetFillStyle(0);
  theEvtQALeg->AddEntry(theDataEvtQAHist, " data" ,"lp");
  theEvtQALeg->AddEntry(theMCEvtQAHist,   " MC"       ,"lp");
  theEvtQALeg->Draw();
  
  // text
  float pp1x=0.14,pp1y=0.82;
  TLatex *pp1=makeTLatex(pp1x,pp1y,promptRecoPDString.c_str());  pp1->Draw();
  TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),PDStatsString.c_str());  pp2->Draw();  
  TLatex *pp4=makeTLatex(pp1x,(pp1y-0.08),Py8TuneString.c_str() );  pp4->Draw();
  


  //ratios go here
  temp_canvEvt->cd(); //change back to main canvas
  
  evtpad2->Draw();
  evtpad2->cd();
  
  
  
  //hist for ratio; use clone
  
  TH1F *theRatio=(TH1F*)theMCEvtQAHist->Clone("MCClone4Ratio");
  TH1F *theDenom=(TH1F*)theDataEvtQAHist->Clone("DataClone4Ratio");
  evtpad2->SetTopMargin(0);
  evtpad2->SetBottomMargin(0.3);
  
  theRatio->GetYaxis()->SetTitleSize(15);
  theRatio->GetYaxis()->SetTitleFont(43);
  theRatio->GetYaxis()->SetTitleOffset(2);
  theRatio->GetYaxis()->SetLabelFont(43); 
  theRatio->GetYaxis()->SetLabelSize(13);
  theRatio->GetXaxis()->SetTitleSize(20);
  theRatio->GetXaxis()->SetTitleFont(43);
  theRatio->GetXaxis()->SetTitleOffset(4.);
  theRatio->GetXaxis()->SetLabelFont(43); 
  theRatio->GetXaxis()->SetLabelSize(15);
  
  theRatio->SetTitle("");
  theRatio->SetXTitle( h_XAx_Title.c_str() );
  theRatio->SetYTitle( ratioTitle.c_str() );

  theRatio->Divide(theDenom);
  theRatio->SetAxisRange(0.0,2.0, "Y");  
  
  theRatio->Draw("HIST E");
  
  float min = theRatio->GetBinLowEdge( 1 );
  float max = theRatio->GetBinLowEdge( (int)theRatio->GetNbinsX() ) 
    + theRatio->GetBinWidth( (int)theRatio->GetNbinsX() ) ;
  
  TLine* lineAtOne          = new TLine(min,1.0,max,1.0); 
  lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
  TLine* lineAtOneHalf      = new TLine(min,0.5,max,0.5); 
  lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
  TLine* lineAtOneEtOneHalf = new TLine(min,1.5,max,1.5); 
  lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
  
  lineAtOne->Draw("same");
  lineAtOneHalf->Draw("same");
  lineAtOneEtOneHalf->Draw("same");	
  
  // print to PDF file
  temp_canvEvt->Print( thePDFFileName.c_str() ); //ratio plots moved above.
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvEvt->DrawClone();
    //outcanv->SetTitle(("Evt Vtx QA Plots Canvas "+inDataHistName).c_str());
    //outcanv->Write(("EvtVtxQAPlots_"+inDataHistName+"_canv").c_str());    
    outcanv->SetTitle(("Evt Vtx QA Plots Canvas "+inMCHistName).c_str());
    outcanv->Write(("EvtVtxQAPlots_"+inMCHistName+"_canv").c_str());    
    
  }
  
  temp_canvEvt->Close();
  return;
}

































void printEvtNrefQAHist( TFile* finData , std::string inDataHistName, 
			 TFile* finMC   , std::string inMCHistName  ,
			 std::string thePDFFileName, std::string fullJetType, 
			 long double theLumi ) {
  
  if(!finData || !finMC ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(funcDebug)std::cout<<"creating temporary canvas for printing Evt plots..."<<std::endl;
  TCanvas *temp_canvEvt = new TCanvas("tempEvt", "temp Evt Canv withLog", CANVX, CANVY); //This is for the Jet QA Evt plots without logs
  temp_canvEvt->cd();
  
  TPad *evtpad1 = new TPad("evtpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  TPad *evtpad2 = new TPad("evtpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
  
  evtpad1->SetBottomMargin(0); 
  evtpad1->SetLogy(1);
  evtpad1->Draw();
  evtpad1->cd();
  
  if(inDataHistName=="hLeadJetPt"){
    evtpad1->SetGridx(1);
    evtpad1->SetGridy(1);
    evtpad1->SetLogx(1);

    evtpad2->SetGridx(1);
    evtpad2->SetGridy(0);
    evtpad2->SetLogx(1);
  }


  

  // PUT JET CUTS SECTION IN HERE
  TH1F* genjtptCut_h= (TH1F*)finMC->Get( "hGenJetPtCut" );
  std::string genjtptCut_str = std::to_string( (int) genjtptCut_h->GetMean() );
  std::cout<<"genjtptCut_str = "<<genjtptCut_str<<std::endl;
  
  TH1F* jtptCut_h= (TH1F*)finData->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}^{gen}>"+genjtptCut_str+" GeV, "+"p_{T}>"+jtptCut_str+" GeV";
  std::cout<<"jetCutString="<<jetCutString<<std::endl;

  TH1F* jtetaLoCut_h= (TH1F*)finData->Get( "hJetEtaCutLo" );
  std::stringstream etaLo; etaLo.precision(1);
  float jtetaLoCut_F=jtetaLoCut_h->GetMean();
  etaLo << std::fixed << jtetaLoCut_F;
  std::string jtetaLoCut_str = etaLo.str();
  std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;
  
  TH1F* jtetaHiCut_h= (TH1F*)finData->Get( "hJetEtaCutHi" );
  std::stringstream etaHi; etaHi.precision(1);
  float jtetaHiCut_F=jtetaHiCut_h->GetMean();
  etaHi << std::fixed << jtetaHiCut_F;
  std::string jtetaHiCut_str = etaHi.str();
  std::cout<<"jtetaHiCut_str = "<<jtetaHiCut_str<<std::endl;
  
  std::string jetEtaCutString=jtetaLoCut_str+"<|y|<"+jtetaHiCut_str;
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  if(inMCHistName=="hWNref"||inMCHistName=="hNref" ){
    jetEtaCutString="|y|<5.1";
  }
  if(inMCHistName=="hWjetsPEvt"||inMCHistName=="hjetPEvt" ){
    jetEtaCutString="|y|<4.7";
  }

  
  
  // open hist
  std::cout<<" drawing "<<inDataHistName<<std::endl<<std::endl;
  TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inDataHistName.c_str() );
  if(!theDataEvtQAHist){ std::cout<<"input data hist not found! skipping hist"<<std::endl;
    return;}
  //theDataEvtQAHist->Scale( 1./theDataEvtQAHist->GetBinWidth(1) );
  dataHistStyle(theDataEvtQAHist);


  
  TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inMCHistName.c_str() );
  if(!theMCEvtQAHist){ std::cout<<"input MC hist not found! skipping hist"<<std::endl;
    return;}
  //theMCEvtQAHist->Scale( 1./theMCEvtQAHist->GetBinWidth(1) );
  //theMCEvtQAHist->SetMarkerStyle(kMultiply);
  //theMCEvtQAHist->SetMarkerSize(0.90);
  //theMCEvtQAHist->SetMarkerColor( theMCOverlayMarkerColor );
  theMCEvtQAHist->SetLineColor( theMCLineColor );  
  
  //hist-type specific scaling
  if(inDataHistName!=inMCHistName)std::cout<<"WARNING: inDataHistName!=inMCHistName!=\"hLeadJetPt\"!!!!!!"<<std::endl;
  
  if(inMCHistName!="hLeadJetPt"){
    //if(inMCHistName=="hWjetsPEvt"||inMCHistName=="hjetsPEvt") 
    theDataEvtQAHist->SetBinContent(1, 0.);
    theMCEvtQAHist->SetBinContent(1, 0.);


    theMCEvtQAHist->Scale( theDataEvtQAHist->Integral(1, (int)theMCEvtQAHist->GetNbinsX() )/
			   theMCEvtQAHist->Integral(   1, (int)theMCEvtQAHist->GetNbinsX() ) 
			   );	
      //else theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );	
  }
  else{// if(inMCHistName=="hLeadJetPt"){
    theDataEvtQAHist=(TH1F*)theDataEvtQAHist->
      TH1::Rebin(nbins_pt_debug, "ppData_analysisBins_evtLeadJetPt", ptbins_debug);
    divideBinWidth(theDataEvtQAHist);
    //theDataEvtQAHist->Scale( 1./theLumi );
    
    theMCEvtQAHist=(TH1F*)theMCEvtQAHist->
      TH1::Rebin(nbins_pt_debug, "ppMC_analysisBins_evtLeadJetPt", ptbins_debug);
    divideBinWidth(theMCEvtQAHist);
    //theMCEvtQAHist->Scale(1000.);
    theMCEvtQAHist->Scale( theDataEvtQAHist->Integral("width")/
			   theMCEvtQAHist->Integral("width") );	
  }
  

  

  // title+axes
  std::string h_Title=PDStatsString_2;
  std::string h_XAx_Title;
  std::string h_YAx_Title ="A.U.";
  
  if(inMCHistName=="hLeadJetPt") h_XAx_Title="Event Maximum Jet p_{T}^{RECO} (GeV)";
  else h_XAx_Title="N_{jets}/Evt";
  
  theDataEvtQAHist->SetTitle( h_Title.c_str() );
  theMCEvtQAHist->SetTitle( h_Title.c_str() );
  theMCEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
  theMCEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
  
  theDataEvtQAHist->Draw("E"); 
  theMCEvtQAHist->Draw("HIST E SAME"); 
  theDataEvtQAHist->Draw("E SAME"); 

  // legend
  float legx1=0.74, legx2=legx1+0.11;
  float legy1=0.74, legy2=legy1+0.09;

  TLegend* theEvtQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");  
  theEvtQALeg->SetFillStyle(0);
  theEvtQALeg->AddEntry(theDataEvtQAHist, "Data" ,"lp");
  theEvtQALeg->AddEntry(theMCEvtQAHist,   "P8 MC"       ,"lp");
  theEvtQALeg->Draw();
  
  // text
  float pp1x=0.48,pp1y=0.84;
  TLatex *pp1=makeTLatex(  pp1x,pp1y,"pp 2015 promptReco, Jet80+LowerJets" );  pp1->Draw();
  TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),"L_{int} = 27.4 pb^{-1}, #sqrt{s} = 5.02 TeV");  pp2->Draw();
  TLatex *pp3=makeTLatex(pp1x,(pp1y-0.08),"Py8 Tune CUETP8M1" );  pp3->Draw();
  TLatex *pp4=makeTLatex(pp1x,(pp1y-0.12), fullJetType+"Jets" );  pp4->Draw();
  TLatex *pp5=makeTLatex(pp1x,(pp1y-0.16), jetCutString );  pp5->Draw();
  TLatex *pp6=makeTLatex(pp1x,(pp1y-0.20), jetEtaCutString );  pp6->Draw();
  
  //ratios go here
  temp_canvEvt->cd(); //change back to main canvas
  evtpad2->Draw();
  evtpad2->cd();
  
  //hist for ratio; use clone

  //TH1F *theRatio=(TH1F*)theDataEvtQAHist->Clone("DataVzClone4Ratio");
  TH1F *theRatio=(TH1F*)theMCEvtQAHist->Clone("MCClone4Ratio");
  evtpad2->SetTopMargin(0);
  evtpad2->SetBottomMargin(0.3);
  
  theRatio->GetYaxis()->SetTitleSize(15);
  theRatio->GetYaxis()->SetTitleFont(43);
  theRatio->GetYaxis()->SetTitleOffset(2);
  theRatio->GetYaxis()->SetLabelFont(43); 
  theRatio->GetYaxis()->SetLabelSize(13);

  theRatio->GetXaxis()->SetTitleSize(20);
  theRatio->GetXaxis()->SetTitleFont(43);
  theRatio->GetXaxis()->SetTitleOffset(4.);
  theRatio->GetXaxis()->SetLabelFont(43); 
  theRatio->GetXaxis()->SetLabelSize(15);
  
  theRatio->SetTitle("");
  theRatio->SetXTitle( h_XAx_Title.c_str() );
  theRatio->SetYTitle( ratioTitle.c_str() );
  theRatio->Divide(theDataEvtQAHist);
  theRatio->SetAxisRange(0.0,2.0, "Y");  
  
  theRatio->Draw("HIST E");
  
  float min = theRatio->GetBinLowEdge( 1 );
  float max = theRatio->GetBinLowEdge( (int)theRatio->GetNbinsX() ) 
    + theRatio->GetBinWidth( (int)theRatio->GetNbinsX() ) ;
  
  TLine* lineAtOne          = new TLine(min,1.0,max,1.0); 
  lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
  TLine* lineAtOneHalf      = new TLine(min,0.5,max,0.5); 
  lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
  TLine* lineAtOneEtOneHalf = new TLine(min,1.5,max,1.5); 
  lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
  
  lineAtOne->Draw("same");
  lineAtOneHalf->Draw("same");
  lineAtOneEtOneHalf->Draw("same");	
  
  // print to PDF file
  temp_canvEvt->Print( thePDFFileName.c_str() ); //ratio plots moved above.
  temp_canvEvt->Close();
  return;
}






void printMCEvtQAHist( TFile* finMC   , std::string inMCHistName , std::string thePDFFileName  , TFile* fout=NULL){
  
  if(!finMC ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(funcDebug)std::cout<<"creating temporary canvas for printing Evt plots..."<<std::endl;
  TCanvas *temp_canvEvt = new TCanvas(("tempEvt_"+inMCHistName).c_str(), "temp Evt Canv withLog", CANVX, CANVY); //This is for the Jet QA Evt plots without logs
  temp_canvEvt->SetLogy(1);
  temp_canvEvt->SetLogx(1);
  temp_canvEvt->cd();  


  std::cout<<" drawing "<<inMCHistName<<std::endl<<std::endl;  
  TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inMCHistName.c_str() );
  if(!theMCEvtQAHist){ std::cout<<"input MC hist not found! skipping hist"<<std::endl;
    return;}
  MCHistStyle(theMCEvtQAHist);
  
  std::string h_Title=PDStatsString_2;
  std::string h_XAx_Title="MC evt #hat{p}_{t} (GeV)"       ;
  std::string h_YAx_Title="A.U."       ;

  theMCEvtQAHist->SetTitle( h_Title.c_str() );
  theMCEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
  theMCEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
  theMCEvtQAHist->SetAxisRange(15.,1000.,"X");
  
  
  if (inMCHistName=="hWeightedpthat") 
    theMCEvtQAHist->SetAxisRange(1.e-7,5.e+5,"Y");
  else//pthat, unweighted
    theMCEvtQAHist->SetAxisRange(1.e+1,5.e+5,"Y");
  
  theMCEvtQAHist->Draw("HIST E");
  
  // test
  float t1x=.65,t1y=.81;
  TLatex *t1;
  if (inMCHistName=="hWeightedpthat") {
    t1=makeTLatex(t1x,t1y,"QCD MC, weighted");
  }
  else{
    t1=makeTLatex(t1x,t1y,"QCD MC, unweighted");
  }
  t1->Draw();	
  
  TLatex *t2=makeTLatex(t1x,t1y-0.04,"Py8 Tune CUETP8M1");
  t2->Draw();	
  
  // print to PDF file
  temp_canvEvt->Print( thePDFFileName.c_str() ); 
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvEvt->DrawClone();
    outcanv->SetTitle(("MC Evt QA Plots Canvas "+inMCHistName).c_str());
    outcanv->Write(("MCEvtQAPlots_"+inMCHistName+"_canv").c_str());    
    //temp_canvEvt->Write();
  }
  
  temp_canvEvt->Close();  
  return;
}




void printJetQAHist( TFile* finData , TFile* finMC, int j, bool doJetIDPlots,
		     std::string inHistName , std::string thePDFFileName , std::string fullJetType, 
		     long double theLumi  ,TFile* fout=NULL) {

  
  std::string jetIDInt;
  if(doJetIDPlots)jetIDInt="1";
  else jetIDInt="0";
  
  
  if(!finData || !finMC ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  
  
  // STRINGS FOR HISTS
  TH1F* genjtptCut_h= (TH1F*)finMC->Get( "hGenJetPtCut" );
  std::string genjtptCut_str = std::to_string( (int) genjtptCut_h->GetMean() );
  std::cout<<"genjtptCut_str = "<<genjtptCut_str<<std::endl;
  
  std::string genJetCutString="p_{T}^{GEN} > "+genjtptCut_str+" GeV";
  std::cout<<"genJetCutString="<<genJetCutString<<std::endl;
  
  TH1F* jtptCut_h= (TH1F*)finData->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}^{RECO} > "+jtptCut_str+" GeV";
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
  TH1F* jtptQACut_h= (TH1F*)finData->Get( "hJetQAPtCut" );
  std::string jtptQACut_str = std::to_string( (int) jtptQACut_h->GetMean() );
  std::cout<<"jtptQACut_str = "<<jtptQACut_str<<std::endl;
  
  TH1F* jtetaLoCut_h= (TH1F*)finData->Get( "hJetEtaCutLo" );
  std::stringstream etaLo; etaLo.precision(1);
  float jtetaLoCut_F=jtetaLoCut_h->GetMean();
  etaLo << std::fixed << jtetaLoCut_F;
  std::string jtetaLoCut_str = etaLo.str();
  std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;
  
  TH1F* jtetaHiCut_h= (TH1F*)finData->Get( "hJetEtaCutHi" );
  std::stringstream etaHi; etaHi.precision(1);
  float jtetaHiCut_F=jtetaHiCut_h->GetMean();
  etaHi << std::fixed << jtetaHiCut_F;
  std::string jtetaHiCut_str = etaHi.str();
  std::cout<<"jtetaHiCut_str = "<<jtetaHiCut_str<<std::endl;
  
  std::string jetEtaCutString;
  if(jtetaLoCut_str=="0.0") jetEtaCutString="|y| < "+jtetaHiCut_str;
  else jetEtaCutString=jtetaLoCut_str+" < |y| < "+jtetaHiCut_str;
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  
  float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  
  
  
  
  
  // OPEN FILES
  std::string inDataHist=inHistName;
  std::cout<<" opening input Data "<<inDataHist<<std::endl<<std::endl;      
  
  TH1F* theDataJetQAHist=(TH1F*) ( (TH1*)finData->Get(inDataHist.c_str()) );
  theDataJetQAHist=(TH1F*) (theDataJetQAHist->Clone((inDataHist+"_"+std::to_string(j)+"clone").c_str()) );
  if(!theDataJetQAHist){ std::cout<<"input data hist not found! skipping hist"<<std::endl;
    return;}
  
  theDataJetQAHist->Scale( 1./theLumi);
  theDataJetQAHist->Scale( 1./etaBinWidth);//eta bin width scaling
  
  
  std::string inMCHist=inHistName;
  float MCscale4Units=1./etaBinWidth;
  if(var[j]=="jtpt_forRes" || var[j]=="jtpt_L2Res" || var[j]=="jtpt_L3Res"){
    inMCHist="hJetQA_"+std::to_string((int)doJetIDPlots)+"wJetID_jtpt";
    //MCscale4Units=1.;
  }
  
  std::cout<<" opening input MC "<<inMCHist<<std::endl<<std::endl;      
  TH1F* theMCJetQAHist= (TH1F*) ( (TH1*)finMC->Get(inMCHist.c_str()) );  
  theMCJetQAHist=(TH1F*) (theMCJetQAHist->Clone((inMCHist+"_"+std::to_string(j)+"clone").c_str()) );
  if(!theMCJetQAHist){ std::cout<<"input MC hist not found! skipping hist"<<std::endl;
    return;}
  theMCJetQAHist->Scale(MCscale4Units);
  
  
  
  // MORE SCALING
  if(var[j]=="jtpt" || var[j]=="rawpt" ||
     var[j]=="jtpt_forRes" ||     var[j]=="jtpt_L2Res" ||     var[j]=="jtpt_L3Res" ||
     var[j]=="leadJetPt" || var[j]=="subleadJetPt") {
    
    theDataJetQAHist=(TH1F*)theDataJetQAHist->TH1::Rebin(nbins_pt_debug, ("ppData_analysisBins_"+var[j]).c_str(), ptbins_debug);
    divideBinWidth(theDataJetQAHist);      
    
    theMCJetQAHist=(TH1F*)theMCJetQAHist->TH1::Rebin(nbins_pt_debug, ("ppMC_analysisBins_"+var[j]).c_str(), ptbins_debug);
    divideBinWidth(theMCJetQAHist); 
      
  }
  
  else {
    
    theDataJetQAHist->TH1::Rebin(var_xAx_reBin[j]);
    theDataJetQAHist->Scale( 1./theDataJetQAHist->GetBinWidth(1) );          
    theDataJetQAHist->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X"); 
    
    theMCJetQAHist->TH1::Rebin(var_xAx_reBin[j]);
    theMCJetQAHist->Scale( 1./theMCJetQAHist->GetBinWidth(1) );          
    theMCJetQAHist->Scale( theDataJetQAHist->Integral()/theMCJetQAHist->Integral() );	        
    theMCJetQAHist->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X"); 
    
  }
  
  if(funcDebug)theDataJetQAHist->Print("base");       std::cout<<std::endl;  
  if(funcDebug)theMCJetQAHist->Print("base");       std::cout<<std::endl;
  


  
  // STYLE  
  dataHistStyle(theDataJetQAHist);
  MCHistStyle(theMCJetQAHist);
  
  
  
  // TITLES
  std::string h_XAx_Title= var_xAx_Titles[j];
  //std::string h_YAx_Title= dcrossSectionAxTitle;
  //std::string h_YAx_Title= AUAxTitle;//for APS DNP
  std::string h_YAx_Title= ddcrossSectionAxTitle;//for APS DNP
  std::string h_Title = PDStatsString_2; //for APS DNP
  
  theDataJetQAHist->SetTitle (    h_Title.c_str() );
  theDataJetQAHist->SetXTitle( h_XAx_Title.c_str() );
  theDataJetQAHist->SetYTitle( h_YAx_Title.c_str() );
  
  theMCJetQAHist->SetTitle (    h_Title.c_str() );
  theMCJetQAHist->SetXTitle( h_XAx_Title.c_str() );
  theMCJetQAHist->SetYTitle( h_YAx_Title.c_str() );
  
  
  
  TH1F* theRatio=(TH1F*)theMCJetQAHist->Clone(  ("MCJetHistClone4Ratio_"+var[j]).c_str());  
  
  
  if(funcDebug) std::cout<<"creating temporary canvas for printing Jet plots..."<<std::endl;
  TCanvas *temp_canvJet = new TCanvas("tempEvt", "temp Jet Canv withLog", CANVX, CANVY); //This is for the Jet QA Evt plots without logs
  temp_canvJet->cd();
  
    
  TPad *jetpad1 = new TPad("jetpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  jetpad1->SetGridy(0);
  jetpad1->SetGridx(0);
  jetpad1->SetLogy(1);


  TPad *jetpad2 = new TPad("jetpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
  jetpad2->SetGridx(1);
  jetpad2->SetLogy(0);
  
  
  if(var[j]=="jtpt"      || var[j]=="rawpt"          ||
     var[j]=="jtpt_forRes"|| var[j]=="jtpt_L2Res"|| var[j]=="jtpt_L3Res"|| 
     var[j]=="leadJetPt" || var[j]=="subleadJetPt"     ){
    
    jetpad1->SetLogx(1);  
    jetpad2->SetLogx(1);  
  }
  else {
    jetpad1->SetLogx(0);
    jetpad2->SetLogx(0);
  }
  
  jetpad1->SetBottomMargin(0);
  
  jetpad2->SetTopMargin(0);
  jetpad2->SetBottomMargin(0.3);
  
  //if(var[j]=="jtpt"      || var[j]=="rawpt"          ||
  //   var[j]=="leadJetPt" || var[j]=="subleadJetPt"     )
  //  jetpad2->SetLogx(1);  
  //else jetpad2->SetLogx(0);
  
  
  
  
  jetpad1->Draw();
  jetpad2->Draw();  
  

  
  // variable specific draw orders, canv, etc.
  if(var[j]=="jtpt"||var[j]=="rawpt" ||
     var[j]=="jtpt_forRes" ||     var[j]=="jtpt_L2Res" ||     var[j]=="jtpt_L3Res"      )         {
    
    jetpad1->cd();
    if(var[j]=="jtpt"||var[j]=="rawpt"){
      theMCJetQAHist->Draw("HIST E"); 
      theDataJetQAHist->Draw("E SAME"); 
    } 
    
    else if(var[j]=="jtpt_forRes"||var[j]=="jtpt_L2Res"||var[j]=="jtpt_L3Res"){
      theDataJetQAHist->Draw("E"); 
      theMCJetQAHist->Draw("HIST E SAME"); 
    } 
    
    float t1Loc1=0.50, t1Loc2=0.84;
    TLatex* t1= makeTLatex(  t1Loc1  ,  t1Loc2      ,  fullJetType+"Jets" );    t1->Draw();
    TLatex *t2= makeTLatex( (t1Loc1) , (t1Loc2-.05) , (jetCutString)      );    t2->Draw();	
    TLatex *t3= makeTLatex( (t1Loc1) , (t1Loc2-.10) , (jetEtaCutString)   );    t3->Draw();	    
    if(doJetIDPlots){
      TLatex *t4= makeTLatex( (t1Loc1), (t1Loc2-.15), "Jet ID Applied" );      t4->Draw();	}
    
    float legx1=0.70, legx2=legx1+0.11;
    float legy1=0.76, legy2=legy1+0.09;
    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->AddEntry(theDataJetQAHist,"Data","lp");
    theJetQALeg->AddEntry(theMCJetQAHist,  "P8 MC"  ,"lp");
    theJetQALeg->Draw(); 
    
    
    jetpad2->cd();
    
    
    TH1F* theDenom=(TH1F*)theDataJetQAHist->Clone(("DataHistClone4Ratio_"+var[j]).c_str());
    
    theRatio->SetLineColor( theMCLineColor );
    
    theRatio->SetAxisRange(0.5,1.5,"Y");
    
    theRatio->GetYaxis()->SetTitleSize(15);
    theRatio->GetYaxis()->SetTitleFont(43);
    theRatio->GetYaxis()->SetTitleOffset(2);
    
    theRatio->GetYaxis()->SetLabelFont(43); 
    theRatio->GetYaxis()->SetLabelSize(13);
    
    theRatio->GetXaxis()->SetTitleSize(20);
    theRatio->GetXaxis()->SetTitleFont(43);
    theRatio->GetXaxis()->SetTitleOffset(4.);
    
    theRatio->GetXaxis()->SetLabelFont(43); 
    theRatio->GetXaxis()->SetLabelSize(18);
    
    theRatio->SetTitle("");
    theRatio->SetXTitle( h_XAx_Title.c_str() );
    theRatio->SetYTitle(ratioTitle.c_str() );
    
    theRatio->Divide(theDenom);
    
    theRatio->GetXaxis()->SetMoreLogLabels(true);
    theRatio->GetXaxis()->SetNoExponent(true);
    
    theRatio->Draw("HIST E");
    
    TLine* lineAtOne          = new TLine(ptbins_debug[0],1.0,ptbins_debug[nbins_pt_debug],1.0); 
    lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
    TLine* lineAtOneHalf      = new TLine(ptbins_debug[0],0.5,ptbins_debug[nbins_pt_debug],0.5); 
    lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
    TLine* lineAtOneEtOneHalf = new TLine(ptbins_debug[0],1.5,ptbins_debug[nbins_pt_debug],1.5); 
    lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
    
    lineAtOne->Draw("same");
    lineAtOneHalf->Draw("same");
    lineAtOneEtOneHalf->Draw("same");
    
    
  }
  
  else if(var[j]=="jteta" ||
	  var[j]=="jtphi"	||
	  ( j>=jetConstits_varStart && j<dijet_varStart )){ 
    
    jetpad1->cd();
    //if(var[j]=="jtphi"){
      //theMCJetQAHist->SetAxisRange(  .02,.8,"Y");
      //theDataJetQAHist->SetAxisRange(.02,.8,"Y");
    //};
    
    theMCJetQAHist->Draw("HIST E"); 
    theDataJetQAHist->Draw("E same");  
    
    float t1Loc1=0.54, t1Loc2=0.82;
    //float t1Loc1=0.50, t1Loc2=0.72; 	    
    
    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());    t1->Draw();
    TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.05),(genJetCutString).c_str());    t2->Draw();
    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.10),(jetCutString).c_str());    t3->Draw();
    TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.15),(jetEtaCutString).c_str());    t4->Draw();		
    if(doJetIDPlots){
      TLatex *t5=makeTLatex((t1Loc1),(t1Loc2-.20),"Jet ID Applied");      t5->Draw();	
    }
    
    float legx1=0.70, legx2=legx1+0.09;
    float legy1=0.76, legy2=legy1+0.09;

    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->AddEntry(theDataJetQAHist,"Data","lp");
    theJetQALeg->AddEntry(theMCJetQAHist,  "P8 MC"  ,"lp");
    theJetQALeg->Draw();
    
    //Ratio Plots
    
    jetpad2->cd();

    //TH1F* theRatio=(TH1F*)theMCJetQAHist->Clone("mcJetHistClone4Ratio");
    //    theRatio->SetMarkerStyle( kMultiply );
    //    theRatio->SetMarkerSize(0.99 );
    //    theRatio->SetMarkerColor( theMCColor );
    theRatio->SetLineColor( theMCLineColor );
    //theRatio->SetLineWidth( 2 );

    theRatio->SetAxisRange(0.,2.,"Y");
    theRatio->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
    
    theRatio->GetYaxis()->SetTitleSize(15);
    theRatio->GetYaxis()->SetTitleFont(43);
    theRatio->GetYaxis()->SetTitleOffset(2);
    theRatio->GetYaxis()->SetLabelFont(43); 
    theRatio->GetYaxis()->SetLabelSize(13);

    theRatio->GetXaxis()->SetTitleSize(20);
    theRatio->GetXaxis()->SetTitleFont(43);
    theRatio->GetXaxis()->SetTitleOffset(4.);
    theRatio->GetXaxis()->SetLabelFont(43); 
    theRatio->GetXaxis()->SetLabelSize(18);
    
    theRatio->SetTitle("");
    theRatio->SetXTitle( h_XAx_Title.c_str() );
    theRatio->SetYTitle(ratioTitle.c_str());
    theRatio->Divide(theDataJetQAHist);
    theRatio->Draw("HIST E");
	  
    TLine* lineAtOne          = new TLine(jetQAxmin[j],1.0,jetQAxmax[j],1.0); 
    lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
    TLine* lineAtOneHalf      = new TLine(jetQAxmin[j],0.5,jetQAxmax[j],0.5); 
    lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
    TLine* lineAtOneEtOneHalf = new TLine(jetQAxmin[j],1.5,jetQAxmax[j],1.5); 
    lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
    
    lineAtOne->Draw("same");
    lineAtOneHalf->Draw("same");
    lineAtOneEtOneHalf->Draw("same");
    
    
  }
  
  else if(var[j]=="Aj"||
	  var[j]=="xj"||
	  var[j]=="dphi"){ 
    
    jetpad1->cd();
    theMCJetQAHist->Draw("HIST E"); 
    theDataJetQAHist->Draw("E SAME"); 
    
    float t1Loc1=0.15, t1Loc2=0.7; 	    	  
    if(var[j]=="Aj") t1Loc1=.7;
    
    
    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());    t1->Draw();
    TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.04),"p_{T}^{ld} > 60 GeV");    t2->Draw();	
    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.08),"p_{T}^{subld} > 30GeV");    t3->Draw();	
    TLatex *t4=makeTLatex(t1Loc1, (t1Loc2-0.14), "#Delta#phi > #frac{2}{3}#pi");      t4->Draw();	
    if(doJetIDPlots){
      TLatex *t5=makeTLatex((t1Loc1),(t1Loc2-.18),"Jet ID Applied");      t5->Draw();	}
    
    float legx1=0.15, legx2=legx1+0.11;
    if(var[j]=="Aj") {legx1=.7; legx2=legx1+0.11;}
    float legy1=0.74, legy2=legy1+0.09;
    
    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->AddEntry(theDataJetQAHist,"Data","lp");
    theJetQALeg->AddEntry(theMCJetQAHist,  "P8 MC"  ,"lp");
    theJetQALeg->Draw(); 
    
    //Ratio Plots
    
    jetpad2->cd();

    //TH1F* theRatio=(TH1F*)theMCJetQAHist->Clone("MCJetHistClone4Ratio");
    //I need to redefine the axis titles here
    std::string h_XAx_ratTitle=var_xAx_Titles[j], h_YAx_ratTitle="Data/MC"    ;
    std::string h_ratTitle   ="JetQA, Data/MC";//, "+var[j];
    if(doJetIDPlots)h_ratTitle+=", w/ JetIDCut";      
    
    //theRatio->SetMarkerStyle( kMultiply );
    //theRatio->SetMarkerSize(0.99 );
    //theRatio->SetMarkerColor( theMCColor );
    theRatio->SetLineColor( theMCLineColor );
    //theRatio->SetLineWidth( 2 );
    
    theRatio->SetAxisRange(0.,2.,"Y");
    theRatio->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
    
    theRatio->GetYaxis()->SetTitleSize(15);
    theRatio->GetYaxis()->SetTitleFont(43);
    theRatio->GetYaxis()->SetTitleOffset(2);
    theRatio->GetYaxis()->SetLabelFont(43); 
    theRatio->GetYaxis()->SetLabelSize(13);
    theRatio->GetXaxis()->SetTitleSize(20);
    theRatio->GetXaxis()->SetTitleFont(43);
    theRatio->GetXaxis()->SetTitleOffset(4.);
    theRatio->GetXaxis()->SetLabelFont(43); 
    theRatio->GetXaxis()->SetLabelSize(15);
	  
    theRatio->SetTitle("");
    theRatio->SetXTitle( h_XAx_Title.c_str() );
    theRatio->SetYTitle(ratioTitle.c_str());
    theRatio->Divide(theDataJetQAHist);
    theRatio->Draw("E");
	  
    TLine* lineAtOne          = new TLine(jetQAxmin[j],1.0,jetQAxmax[j],1.0); 
    lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
    TLine* lineAtOneHalf      = new TLine(jetQAxmin[j],0.5,jetQAxmax[j],0.5); 
    lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
    TLine* lineAtOneEtOneHalf = new TLine(jetQAxmin[j],1.5,jetQAxmax[j],1.5); 
    lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
	  
    lineAtOne->Draw("same");
    lineAtOneHalf->Draw("same");
    lineAtOneEtOneHalf->Draw("same");
    
  }
      
  else if(var[j]=="leadJetPt"||var[j]=="subleadJetPt"){ 
    jetpad1->cd();
    theMCJetQAHist->Draw("HIST E"); 
    theDataJetQAHist->Draw("E SAME"); 
    
    float t1Loc1=0.7, t1Loc2=0.72; 	    	  
    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());    t1->Draw();
    TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.05),"p_{T}^{ld} > 60 GeV");    t2->Draw();	
    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.10),"p_{T}^{subld} > 30GeV");    t3->Draw();	
    TLatex *t4=makeTLatex(t1Loc1, (t1Loc2-0.17), "#Delta#phi > #frac{2}{3}#pi");	        t4->Draw();	
    if(doJetIDPlots){
      TLatex *t5=makeTLatex((t1Loc1),(t1Loc2-.22),"Jet ID Applied");      t5->Draw();	}
    
    float legx1=0.70, legx2=legx1+0.11;
    float legy1=0.76, legy2=legy1+0.09;
    
    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->AddEntry(theDataJetQAHist,"Data","lp");
    theJetQALeg->AddEntry(theMCJetQAHist,  "P8 MC"  ,"lp");
    theJetQALeg->Draw(); 
    
    //Ratio Plots
    
    jetpad2->cd();
    
    //TH1F* theRatio=(TH1F*)theMCJetQAHist->Clone("MCJetHistClone4Ratio");
    //theRatio->SetMarkerStyle( kMultiply );
    //theRatio->SetMarkerSize(0.99 );
    //theRatio->SetMarkerColor( theMCColor );
    theRatio->SetLineColor( theMCLineColor );
    //theRatio->SetLineWidth( 2 );
    
    theRatio->SetAxisRange(0.,2.,"Y");
    theRatio->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
    
    theRatio->GetYaxis()->SetTitleSize(15);
    theRatio->GetYaxis()->SetTitleFont(43);
    theRatio->GetYaxis()->SetTitleOffset(2);
    theRatio->GetYaxis()->SetLabelFont(43); 
    theRatio->GetYaxis()->SetLabelSize(13);
    theRatio->GetXaxis()->SetTitleSize(20);
    theRatio->GetXaxis()->SetTitleFont(43);
    theRatio->GetXaxis()->SetTitleOffset(4.);
    theRatio->GetXaxis()->SetLabelFont(43); 
    theRatio->GetXaxis()->SetLabelSize(15);
    
    theRatio->SetTitle ("");
    theRatio->SetXTitle( h_XAx_Title.c_str() );
    theRatio->SetYTitle( ratioTitle.c_str() );
    
    theRatio->Divide(theDataJetQAHist);
    theRatio->Draw("HIST E");
    
    TLine* lineAtOne          = new TLine(ptbins_debug[0],1.0,jetQAxmax[j],1.0); 
    lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
    TLine* lineAtOneHalf      = new TLine(ptbins_debug[0],0.5,jetQAxmax[j],0.5); 
    lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
    TLine* lineAtOneEtOneHalf = new TLine(ptbins_debug[0],1.5,jetQAxmax[j],1.5); 
    lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
    
    lineAtOne->Draw("same");
    lineAtOneHalf->Draw("same");
    lineAtOneEtOneHalf->Draw("same");
  }
  
  
  // print to PDF file
  temp_canvJet->Print( thePDFFileName.c_str() );
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJet->DrawClone();
    outcanv->SetTitle((inHistName+" Canvas").c_str());
    outcanv->Write((inHistName+"_canv").c_str());    
    
    theDataJetQAHist->Write(("DATATH1_"+inHistName).c_str());
    theMCJetQAHist  ->Write(("MCTH1_"+inHistName).c_str());    
    theRatio        ->Write(("RATIOTH1_"+inHistName).c_str());
  }
  temp_canvJet->Close();  
  return;
}





void printJetMultQATH1_vtx( TFile* finData ,  TFile* finMC, std::string inHistName, 
			    //			    TFile* finMC   , std::string inHistName  ,
			    std::string h_XAx_Title , int rebin, std::string thePDFFileName , long double theLumi) {
  
  if(!finData || !finMC ){    std::cout<<"input file(s) not found, cannot run printJetMultQATH1_vtx"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(funcDebug)std::cout<<"creating temporary canvas for printing jetMultQA plots..."<<std::endl;
  TCanvas *temp_canvEvt = new TCanvas(("jetMultCanv"+inHistName).c_str(), "temp jetMult Canv", 1000, 800); //This is for the Jet QA Evt plots without logs
  temp_canvEvt->cd();
  
  //if( ((std::size_t)inHistName.find("_N"))!=std::string::npos || 
  //    ((std::size_t)inHistName.find(""))!=std::string::npos 
  //    )
  //  temp_canvEvt->SetLogy(1);
  //if(findsubstr(inHistName, "_N") ||
  //findsubstr(inHistName, "vtx_r") )
  
  if( findsubstr(inHistName, "vtx_z") )
    temp_canvEvt->SetLogy(0);
  else 
    temp_canvEvt->SetLogy(1);  
  //else if( findsubstr(inHistName, "Nvtxtrks_vtx0") )
  //  temp_canvEvt->SetLogy(0);
  //else if( findsubstr(inHistName, "Nvtxtrks_vtxgt0") )
  //  temp_canvEvt->SetLogy(0);
  //else if( inHistName=="Nvtxtrks" )
  //  temp_canvEvt->SetLogy(0);
  
  
  //  TPad *evtpad1 = new TPad("evtpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  //  TPad *evtpad2 = new TPad("evtpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
  //  
  //  evtpad1->SetBottomMargin(0); 
  //  evtpad1->Draw();
  //  evtpad1->cd();
  
  
  // open hist
  std::cout<<" getting hist(s): "<<inHistName<<std::endl;
  
  TH1F* theDataJetMultTH1= (TH1F*)finData->Get( inHistName.c_str() );
  if(!theDataJetMultTH1) 
    std::cout<<"input data hist not found!"<<std::endl;
  
  TH1F* theMCJetMultTH1= (TH1F*)finMC->Get( inHistName.c_str() );
  if(!theMCJetMultTH1) 
    std::cout<<"input MC hist not found!"<<std::endl;
  
  if(!theDataJetMultTH1 || !theMCJetMultTH1) { 
    std::cout<<"skipping hist!"<<std::endl<<std::endl;
    return;   }
  
  theDataJetMultTH1=(TH1F*)theDataJetMultTH1->TH1::Rebin( (int)rebin, (inHistName+"_rebin").c_str());
  theDataJetMultTH1->Scale( 1./theDataJetMultTH1->GetBinWidth(1) );
  //theDataJetMultTH1->Scale( 1./theLumi );
  dataHistStyle(theDataJetMultTH1);
  
  
  theMCJetMultTH1=(TH1F*)theMCJetMultTH1->TH1::Rebin( (int)rebin, (inHistName+"_rebin").c_str());
  theMCJetMultTH1->Scale( 1./theMCJetMultTH1->GetBinWidth(1) );
  //theMCJetMultTH1->Scale( 1000. );
  theMCJetMultTH1->Scale( theDataJetMultTH1->Integral()/theMCJetMultTH1->Integral() );	
  MCHistStyle(theMCJetMultTH1);
  
  
  
  
  
  


  // title+axes
  
  std::string h_Title     =PDStatsString_2;
  std::string h_YAx_Title ="A.U.";  

  //std::string h_XAx_Title="";
  //if( inHistName.find("Vz") != std::string::npos ){
  //  h_XAx_Title ="v_{z}^{evt} (cm)";}
  //else if( inHistName.find("Vx") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{x}^{evt} (cm)";}
  //else if( inHistName.find("Vy") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{y}^{evt} (cm)";}
  //else if( inHistName.find("Vr") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{r}^{evt} (cm)";}
  
  //Double_t dataMax = theDataJetMultTH1->GetMaximum();//, MCMax=theMCJetMultTH1->GetMaximum();
  //std::cout<<"dataMax="<<dataMax<<std::endl;
  //std::cout<<"MCMax="<<    MCMax<<std::endl;
  //if(dataMax>MCMax){
  //  std::cout<<"dataMax > MC Max"<<std::endl;
  //  theDataJetMultTH1->SetTitle (    h_Title.c_str() );
  //  theDataJetMultTH1->SetXTitle( h_XAx_Title.c_str() );
  //  theDataJetMultTH1->SetYTitle( h_YAx_Title.c_str() );
  //  
  //  theDataJetMultTH1->Draw("E"); 
  //  theMCJetMultTH1->Draw("HIST E SAME"); 
  //  theDataJetMultTH1->Draw("E SAME"); //because i want data points written over MC's red line hist
  //}
  //else {
  //  theMCJetMultTH1->SetTitle (    h_Title.c_str() );
  //  theMCJetMultTH1->SetXTitle( h_XAx_Title.c_str() );
  //  theMCJetMultTH1->SetYTitle( h_YAx_Title.c_str() );
  //  theMCJetMultTH1->SetAxisRange(.001,10.*MCMax,"Y");
  //  theMCJetMultTH1->Draw("HIST E"); 
  //  theDataJetMultTH1->Draw("E SAME"); 
  //}
  
  theDataJetMultTH1->SetTitle (    h_Title.c_str() );
  theDataJetMultTH1->SetXTitle( h_XAx_Title.c_str() );
  theDataJetMultTH1->SetYTitle( h_YAx_Title.c_str() );
  
  theDataJetMultTH1->Draw("E"); 
  theMCJetMultTH1->Draw("HIST E SAME");

  // legend
  float legx1=0.80, legx2=legx1+0.10;
  float legy1=0.80, legy2=legy1+0.10;
  
  TLegend* theEvtQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");
  //theEvtQALeg->SetFillStyle(0);
  theEvtQALeg->AddEntry(theDataJetMultTH1, " data" ,"lp");
  theEvtQALeg->AddEntry(theMCJetMultTH1,   " MC"       ,"lp");
  theEvtQALeg->Draw();
  
  //// text
  //float pp1x=0.55,pp1y=0.85;
  //TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str());  pp1->Draw();
  //TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str());  pp2->Draw();  
  //TLatex *pp4=makeTLatex(pp1x,(pp1y-0.08),Py8TuneString.c_str() );  pp4->Draw();  
  
  
  //ratios go here
  //temp_canvEvt->cd(); //change back to main canvas
  
  //evtpad2->Draw();
  //evtpad2->cd();
  
  
  
  //hist for ratio; use clone
  
  //TH1F *theRatio=(TH1F*)theMCJetMultTH1->Clone("MCClone4Ratio");
  //TH1F *theDenom=(TH1F*)theDataJetMultTH1->Clone("DataClone4Ratio");
  //evtpad2->SetTopMargin(0);
  //evtpad2->SetBottomMargin(0.3);
  //
  //theRatio->GetYaxis()->SetTitleSize(15);
  //theRatio->GetYaxis()->SetTitleFont(43);
  //theRatio->GetYaxis()->SetTitleOffset(2);
  //theRatio->GetYaxis()->SetLabelFont(43); 
  //theRatio->GetYaxis()->SetLabelSize(13);
  //theRatio->GetXaxis()->SetTitleSize(20);
  //theRatio->GetXaxis()->SetTitleFont(43);
  //theRatio->GetXaxis()->SetTitleOffset(4.);
  //theRatio->GetXaxis()->SetLabelFont(43); 
  //theRatio->GetXaxis()->SetLabelSize(15);
  //
  //theRatio->SetTitle("");
  //theRatio->SetXTitle( h_XAx_Title.c_str() );
  //theRatio->SetYTitle( ratioTitle.c_str() );
  //
  //theRatio->Divide(theDenom);
  //theRatio->SetAxisRange(0.0,2.0, "Y");  
  //
  //theRatio->Draw("HIST E");
  //
  //float min = theRatio->GetBinLowEdge( 1 );
  //float max = theRatio->GetBinLowEdge( (int)theRatio->GetNbinsX() ) 
  //  + theRatio->GetBinWidth( (int)theRatio->GetNbinsX() ) ;
  //
  //TLine* lineAtOne          = new TLine(min,1.0,max,1.0); 
  //lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
  //TLine* lineAtOneHalf      = new TLine(min,0.5,max,0.5); 
  //lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
  //TLine* lineAtOneEtOneHalf = new TLine(min,1.5,max,1.5); 
  //lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
  //
  //lineAtOne->Draw("same");
  //lineAtOneHalf->Draw("same");
  //lineAtOneEtOneHalf->Draw("same");	
  
  // print to PDF file
  temp_canvEvt->Print( thePDFFileName.c_str() ); //ratio plots moved above.
  temp_canvEvt->Close();
  std::cout<<std::endl;
  return;
}





//void printJetMultQATH1_trk( TFile* finData , std::string inDataHistName, 
void printJetMultQATH1_trk( TFile* finData , TFile* finMC, std::string inHistName, 
			    std::string h_XAx_Title , int rebin, std::string thePDFFileName , long double theLumi) {
  
  if(!finData || !finMC ){    std::cout<<"input file(s) not found, cannot run printJetMultQATH1_vtx"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(funcDebug)std::cout<<"creating temporary canvas for printing jetMultQA plots..."<<std::endl;
  TCanvas *temp_canvEvt = new TCanvas(("jetMultCanv"+inHistName).c_str(), "temp jetMult Canv", 1000, 800); //This is for the Jet QA Evt plots without logs
  temp_canvEvt->cd();
  
  //if( ((std::size_t)inHistName.find("_N"))!=std::string::npos || 
  //    ((std::size_t)inHistName.find(""))!=std::string::npos 
  //    )
  //  temp_canvEvt->SetLogy(1);
  //if(findsubstr(inHistName, "_N") ||
  //findsubstr(inHistName, "vtx_r") )
  temp_canvEvt->SetLogy(1);
  if(findsubstr(inHistName, "trkPtPercError"))
    temp_canvEvt->SetLogx(0);
  else if(findsubstr(inHistName, "trkPtError"))
    temp_canvEvt->SetLogx(0);
  else if(findsubstr(inHistName, "trkPt"))
    temp_canvEvt->SetLogx(0);
  //  TPad *evtpad1 = new TPad("evtpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  //  TPad *evtpad2 = new TPad("evtpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
  //  
  //  evtpad1->SetBottomMargin(0); 
  //  evtpad1->Draw();
  //  evtpad1->cd();
  
  
  // open hist
  std::cout<<" getting hist(s): "<<inHistName<<std::endl;
  
  TH1F* theDataJetMultTH1= (TH1F*)finData->Get( inHistName.c_str() );
  if(!theDataJetMultTH1){ std::cout<<"input data hist not found! skipping hist"<<std::endl;
    return;}
  theDataJetMultTH1=(TH1F*)theDataJetMultTH1->TH1::Rebin( (int)rebin, (inHistName+"_rebin").c_str());
  theDataJetMultTH1->Scale( 1./theDataJetMultTH1->GetBinWidth(1) );
  //theDataJetMultTH1->Scale( 1./theLumi );
  dataHistStyle(theDataJetMultTH1);
  
  
  TH1F* theMCJetMultTH1= (TH1F*)finMC->Get( inHistName.c_str() );
  if(!theMCJetMultTH1){ std::cout<<"input MC hist not found! skipping hist"<<std::endl;
    return;}
  theMCJetMultTH1=(TH1F*)theMCJetMultTH1->TH1::Rebin( (int)rebin, (inHistName+"_rebin").c_str());
  theMCJetMultTH1->Scale( 1./theMCJetMultTH1->GetBinWidth(1) );
  //theMCJetMultTH1->Scale(1000.);
  theMCJetMultTH1->Scale( theDataJetMultTH1->Integral()/theMCJetMultTH1->Integral() );	
  MCHistStyle(theMCJetMultTH1);
  
  
  
  
  // title+axes
  
  std::string h_Title     =PDStatsString_2;
  std::string h_YAx_Title ="A.U.";  
  //std::string h_XAx_Title="";
  //if( inHistName.find("Vz") != std::string::npos ){
  //  h_XAx_Title ="v_{z}^{evt} (cm)";}
  //else if( inHistName.find("Vx") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{x}^{evt} (cm)";}
  //else if( inHistName.find("Vy") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{y}^{evt} (cm)";}
  //else if( inHistName.find("Vr") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{r}^{evt} (cm)";}
  
  //Double_t dataMax = theDataJetMultTH1->GetMaximum();//, MCMax=theMCJetMultTH1->GetMaximum();
  //std::cout<<"dataMax="<<dataMax<<std::endl;
  //std::cout<<"MCMax="<<    MCMax<<std::endl;
  //if(dataMax>MCMax){
  //  std::cout<<"dataMax > MC Max"<<std::endl;
  //  theDataJetMultTH1->SetTitle (    h_Title.c_str() );
  //  theDataJetMultTH1->SetXTitle( h_XAx_Title.c_str() );
  //  theDataJetMultTH1->SetYTitle( h_YAx_Title.c_str() );
  //  
  //  theDataJetMultTH1->Draw("E"); 
  //  theMCJetMultTH1->Draw("HIST E SAME"); 
  //  theDataJetMultTH1->Draw("E SAME"); //because i want data points written over MC's red line hist
  //}
  //else {
  //  theMCJetMultTH1->SetTitle (    h_Title.c_str() );
  //  theMCJetMultTH1->SetXTitle( h_XAx_Title.c_str() );
  //  theMCJetMultTH1->SetYTitle( h_YAx_Title.c_str() );
  //  theMCJetMultTH1->SetAxisRange(.001,10.*MCMax,"Y");
  //  theMCJetMultTH1->Draw("HIST E"); 
  //  theDataJetMultTH1->Draw("E SAME"); 
  //}
  
  if(findsubstr(inHistName, "trkPtPercError"))
    theDataJetMultTH1->SetAxisRange(0, 20., "X");
  
  theDataJetMultTH1->SetTitle (    h_Title.c_str() );
  theDataJetMultTH1->SetXTitle( h_XAx_Title.c_str() );
  theDataJetMultTH1->SetYTitle( h_YAx_Title.c_str() );
  
  theDataJetMultTH1->Draw("E"); 
  theMCJetMultTH1->Draw("HIST E SAME"); 
  
  // legend
  float legx1=0.80, legx2=legx1+0.10;
  float legy1=0.80, legy2=legy1+0.10;
  
  TLegend* theEvtQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");
  theEvtQALeg->SetFillStyle(0);
  theEvtQALeg->AddEntry(theDataJetMultTH1, " data" ,"lp");
  theEvtQALeg->AddEntry(theMCJetMultTH1,   " MC"       ,"lp");
  theEvtQALeg->Draw();
  
  // text
  float pp1x=0.55,pp1y=0.85;
  TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str());  pp1->Draw();
  TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str());  pp2->Draw();  
  //TLatex *pp4=makeTLatex(pp1x,(pp1y-0.08),Py8TuneString.c_str() );  pp4->Draw();
  


  //ratios go here
  //temp_canvEvt->cd(); //change back to main canvas
  
  //evtpad2->Draw();
  //evtpad2->cd();
  
  
  
  //hist for ratio; use clone
  
  //TH1F *theRatio=(TH1F*)theMCJetMultTH1->Clone("MCClone4Ratio");
  //TH1F *theDenom=(TH1F*)theDataJetMultTH1->Clone("DataClone4Ratio");
  //evtpad2->SetTopMargin(0);
  //evtpad2->SetBottomMargin(0.3);
  //
  //theRatio->GetYaxis()->SetTitleSize(15);
  //theRatio->GetYaxis()->SetTitleFont(43);
  //theRatio->GetYaxis()->SetTitleOffset(2);
  //theRatio->GetYaxis()->SetLabelFont(43); 
  //theRatio->GetYaxis()->SetLabelSize(13);
  //theRatio->GetXaxis()->SetTitleSize(20);
  //theRatio->GetXaxis()->SetTitleFont(43);
  //theRatio->GetXaxis()->SetTitleOffset(4.);
  //theRatio->GetXaxis()->SetLabelFont(43); 
  //theRatio->GetXaxis()->SetLabelSize(15);
  //
  //theRatio->SetTitle("");
  //theRatio->SetXTitle( h_XAx_Title.c_str() );
  //theRatio->SetYTitle( ratioTitle.c_str() );
  //
  //theRatio->Divide(theDenom);
  //theRatio->SetAxisRange(0.0,2.0, "Y");  
  //
  //theRatio->Draw("HIST E");
  //
  //float min = theRatio->GetBinLowEdge( 1 );
  //float max = theRatio->GetBinLowEdge( (int)theRatio->GetNbinsX() ) 
  //  + theRatio->GetBinWidth( (int)theRatio->GetNbinsX() ) ;
  //
  //TLine* lineAtOne          = new TLine(min,1.0,max,1.0); 
  //lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
  //TLine* lineAtOneHalf      = new TLine(min,0.5,max,0.5); 
  //lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
  //TLine* lineAtOneEtOneHalf = new TLine(min,1.5,max,1.5); 
  //lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
  //
  //lineAtOne->Draw("same");
  //lineAtOneHalf->Draw("same");
  //lineAtOneEtOneHalf->Draw("same");	
  
  // print to PDF file
  temp_canvEvt->Print( thePDFFileName.c_str() ); //ratio plots moved above.
  temp_canvEvt->Close();
  std::cout<<std::endl;
  return;
}




void printJetMultQATH2( TFile* finData , std::string inDataHistName, 
			//TFile* finMC   , std::string inMCHistName  ,
			std::string h_XAx_Title, int xrebin, 
			std::string h_YAx_Title , int yrebin,
			std::string thePDFFileName ) {
  
  //if(!finData || !finMC ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
  if(!finData ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(funcDebug)std::cout<<"creating temporary canvas for printing jetMultQA plots..."<<std::endl;
  TCanvas *temp_canvEvt = new TCanvas(("jetMultCanv"+inDataHistName).c_str(), "temp jetMult Canv", 1200, 800); //This is for the Jet QA Evt plots without logs
  temp_canvEvt->cd();
  
  //if( ((std::size_t)inDataHistName.find("_N"))!=std::string::npos || 
  //    ((std::size_t)inDataHistName.find(""))!=std::string::npos 
  //    )
  //  temp_canvEvt->SetLogy(1);
  //if(findsubstr(inDataHistName, "_N") ||
  //findsubstr(inDataHistName, "vtx_r") )
  //  temp_canvEvt->SetLogy(1);
  
  if(findsubstr(inDataHistName,"trkEtaVPhi"))
    temp_canvEvt->SetLogz(0);
  else
    temp_canvEvt->SetLogz(1);
  
  //  TPad *evtpad1 = new TPad("evtpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  //  TPad *evtpad2 = new TPad("evtpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
  //  
  //  evtpad1->SetBottomMargin(0); 
  //  evtpad1->Draw();
  //  evtpad1->cd();
  
  
  // open hist
  std::cout<<" getting hist(s): "<<inDataHistName<<std::endl;
  
  TH2F* theDataJetMultTH2= (TH2F*)finData->Get( inDataHistName.c_str() );
  if(!theDataJetMultTH2){ std::cout<<"input data hist not found! skipping hist"<<std::endl;
    return;}
  theDataJetMultTH2=(TH2F*)theDataJetMultTH2->TH2::Rebin2D( (int)xrebin, (int)yrebin, (inDataHistName+"_rebin").c_str());
  
  float scaleFactor = theDataJetMultTH2->GetXaxis()->GetBinWidth(1);
  scaleFactor *= theDataJetMultTH2->GetYaxis()->GetBinWidth(1);
  theDataJetMultTH2->Scale( 1./scaleFactor );
  //theDataJetMultTH2->Scale( 1./theLumi );
  //dataHistStyle(theDataJetMultTH2);
  
  
  //TH2F* theMCJetMultTH2= (TH2F*)finMC->Get( inMCHistName.c_str() );
  //if(!theMCJetMultTH2){ std::cout<<"input MC hist not found! skipping hist"<<std::endl;
  //  return;}
  //theMCJetMultTH2=(TH2F*)theMCJetMultTH2->TH1::Rebin( (int)rebin, (inMCHistName+"_rebin").c_str());
  //theMCJetMultTH2->Scale( 1./theMCJetMultTH2->GetBinWidth(1) );
  //theMCJetMultTH2->Scale( theDataJetMultTH2->Integral()/theMCJetMultTH2->Integral() );	
  //MCHistStyle(theMCJetMultTH2);
  
  


  // title+axes
  
  std::string h_Title     =PDStatsString_2;
  //std::string h_YAx_Title ="A.U.";  
  //std::string h_XAx_Title="";
  //if( inDataHistName.find("Vz") != std::string::npos ){
  //  h_XAx_Title ="v_{z}^{evt} (cm)";}
  //else if( inDataHistName.find("Vx") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{x}^{evt} (cm)";}
  //else if( inDataHistName.find("Vy") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{y}^{evt} (cm)";}
  //else if( inDataHistName.find("Vr") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{r}^{evt} (cm)";}
  
  //Double_t dataMax = theDataJetMultTH2->GetMaximum();//, MCMax=theMCJetMultTH2->GetMaximum();
  //std::cout<<"dataMax="<<dataMax<<std::endl;
  //std::cout<<"MCMax="<<    MCMax<<std::endl;
  //if(dataMax>MCMax){
  //  std::cout<<"dataMax > MC Max"<<std::endl;
  //  theDataJetMultTH2->SetTitle (    h_Title.c_str() );
  //  theDataJetMultTH2->SetXTitle( h_XAx_Title.c_str() );
  //  theDataJetMultTH2->SetYTitle( h_YAx_Title.c_str() );
  //  
  //  theDataJetMultTH2->Draw("E"); 
  //  theMCJetMultTH2->Draw("HIST E SAME"); 
  //  theDataJetMultTH2->Draw("E SAME"); //because i want data points written over MC's red line hist
  //}
  //else {
  //  theMCJetMultTH2->SetTitle (    h_Title.c_str() );
  //  theMCJetMultTH2->SetXTitle( h_XAx_Title.c_str() );
  //  theMCJetMultTH2->SetYTitle( h_YAx_Title.c_str() );
  //  theMCJetMultTH2->SetAxisRange(.001,10.*MCMax,"Y");
  //  theMCJetMultTH2->Draw("HIST E"); 
  //  theDataJetMultTH2->Draw("E SAME"); 
  //}
  
  theDataJetMultTH2->SetTitle (    h_Title.c_str() );
  theDataJetMultTH2->SetXTitle( h_XAx_Title.c_str() );
  theDataJetMultTH2->SetYTitle( h_YAx_Title.c_str() );
  
  theDataJetMultTH2->Draw("COLZ"); 
  
  //// legend
  //float legx1=0.80, legx2=legx1+0.10;
  //float legy1=0.80, legy2=legy1+0.10;
  //
  //TLegend* theEvtQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");
  //theEvtQALeg->SetFillStyle(0);
  //theEvtQALeg->AddEntry(theDataJetMultTH2, " data" ,"lp");
  //theEvtQALeg->AddEntry(theMCJetMultTH2,   " MC"       ,"lp");
  //theEvtQALeg->Draw();
  
  // text
  //float pp1x=0.55,pp1y=0.85;
  //TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str());  //pp1->Draw();
  //TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str()); // pp2->Draw();  
  //TLatex *pp4=makeTLatex(pp1x,(pp1y-0.08),Py8TuneString.c_str() );  pp4->Draw();
  
  
  
  //ratios go here
  //temp_canvEvt->cd(); //change back to main canvas
  
  //evtpad2->Draw();
  //evtpad2->cd();
  
  
  
  //hist for ratio; use clone
  
  //TH1F *theRatio=(TH1F*)theMCJetMultTH1->Clone("MCClone4Ratio");
  //TH1F *theDenom=(TH1F*)theDataJetMultTH1->Clone("DataClone4Ratio");
  //evtpad2->SetTopMargin(0);
  //evtpad2->SetBottomMargin(0.3);
  //
  //theRatio->GetYaxis()->SetTitleSize(15);
  //theRatio->GetYaxis()->SetTitleFont(43);
  //theRatio->GetYaxis()->SetTitleOffset(2);
  //theRatio->GetYaxis()->SetLabelFont(43); 
  //theRatio->GetYaxis()->SetLabelSize(13);
  //theRatio->GetXaxis()->SetTitleSize(20);
  //theRatio->GetXaxis()->SetTitleFont(43);
  //theRatio->GetXaxis()->SetTitleOffset(4.);
  //theRatio->GetXaxis()->SetLabelFont(43); 
  //theRatio->GetXaxis()->SetLabelSize(15);
  //
  //theRatio->SetTitle("");
  //theRatio->SetXTitle( h_XAx_Title.c_str() );
  //theRatio->SetYTitle( ratioTitle.c_str() );
  //
  //theRatio->Divide(theDenom);
  //theRatio->SetAxisRange(0.0,2.0, "Y");  
  //
  //theRatio->Draw("HIST E");
  //
  //float min = theRatio->GetBinLowEdge( 1 );
  //float max = theRatio->GetBinLowEdge( (int)theRatio->GetNbinsX() ) 
  //  + theRatio->GetBinWidth( (int)theRatio->GetNbinsX() ) ;
  //
  //TLine* lineAtOne          = new TLine(min,1.0,max,1.0); 
  //lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
  //TLine* lineAtOneHalf      = new TLine(min,0.5,max,0.5); 
  //lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
  //TLine* lineAtOneEtOneHalf = new TLine(min,1.5,max,1.5); 
  //lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
  //
  //lineAtOne->Draw("same");
  //lineAtOneHalf->Draw("same");
  //lineAtOneEtOneHalf->Draw("same");	
  
  // print to PDF file
  temp_canvEvt->Print( thePDFFileName.c_str() ); //ratio plots moved above.
  temp_canvEvt->Close();
  std::cout<<std::endl;
  return;
}




void printJetIDHist( TFile* fin_jetID , TFile* fin_nojetID, 
		     int j, bool isData,
		     std::string thePDFFileName , std::string fullJetType, 
		     long double theLumi,  
		     TFile* fout=NULL) {
  
  
  if(!fin_jetID || !fin_nojetID ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
    
  // PUT JET CUTS SECTION IN HERE
  TH1F* genjtptCut_h=NULL; 
  std::string genjtptCut_str;
  if(!isData){
    genjtptCut_h=(TH1F*)fin_jetID->Get( "hGenJetPtCut" );
    genjtptCut_str = std::to_string( (int) genjtptCut_h->GetMean() );
    std::cout<<"genjtptCut_str = "<<genjtptCut_str<<std::endl;}
  
  TH1F* jtptCut_h= (TH1F*)fin_jetID->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString;
  if(!isData)jetCutString="p_{T}^{gen}>"+genjtptCut_str+" GeV, "+"p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  else jetCutString="p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
  TH1F* jtptQACut_h= (TH1F*)fin_jetID->Get( "hJetQAPtCut" );
  std::string jtptQACut_str = std::to_string( (int) jtptQACut_h->GetMean() );
  std::cout<<"jtptQACut_str = "<<jtptQACut_str<<std::endl;
  
  TH1F* jtetaLoCut_h= (TH1F*)fin_jetID->Get( "hJetEtaCutLo" );
  std::stringstream etaLo; etaLo.precision(1);
  float jtetaLoCut_F=jtetaLoCut_h->GetMean();
  etaLo << std::fixed << jtetaLoCut_F;
  std::string jtetaLoCut_str = etaLo.str(); 
 std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;
 
 TH1F* jtetaHiCut_h= (TH1F*)fin_jetID->Get( "hJetEtaCutHi" );
  std::stringstream etaHi; etaHi.precision(1);
  float jtetaHiCut_F=jtetaHiCut_h->GetMean();
  etaHi << std::fixed << jtetaHiCut_F;
  std::string jtetaHiCut_str = etaHi.str();
  std::cout<<"jtetaHiCut_str = "<<jtetaHiCut_str<<std::endl;
  
  std::string jetEtaCutString=jtetaLoCut_str+"<|y|<"+jtetaHiCut_str;
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  
  float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  
  std::string sampleDescString;
  if(isData)sampleDescString=sqrts2k15ppString;
  else      sampleDescString=Py8TuneString;
  
  std::string inHistName="hJetQA_0wJetID_"+var[j];
  if(funcDebug)std::cout<<" opening input hist "<<inHistName<<std::endl<<std::endl;      
  
  TH1F* theNonJetIDHist=(TH1F*) ( (TH1*)fin_nojetID->Get(inHistName.c_str()) );
  if(!theNonJetIDHist){ std::cout<<"input hist not found! skipping hist"<<std::endl;
    return;}
  
  std::string inJetIDHistName="hJetQA_1wJetID_"+var[j];
  if(funcDebug)std::cout<<" opening input jetID hist "<<inJetIDHistName<<std::endl<<std::endl;      
  TH1F* theJetIDHist= (TH1F*) ( (TH1*)fin_jetID->Get(inJetIDHistName.c_str()) );  
  if(!theJetIDHist){ std::cout<<"input jetID hist not found! skipping hist"<<std::endl;
    return;}
  
  // MORE SCALING
  if(var[j]=="jtpt" || var[j]=="rawpt" ||
     var[j]=="leadJetPt" || var[j]=="subleadJetPt"      ||
     var[j]=="jtpt_forRes" ||     var[j]=="jtpt_L2Res" ||     var[j]=="jtpt_L3Res"      ) {
    
    theNonJetIDHist=(TH1F*)theNonJetIDHist->TH1::Rebin(nbins_pt_debug, "ppData_analysisBins_jtPt", ptbins_debug);
    divideBinWidth(theNonJetIDHist);      
    if(isData)theNonJetIDHist->Scale(1./theLumi);
    theNonJetIDHist->Scale(1./etaBinWidth);
    
    theJetIDHist=(TH1F*)theJetIDHist->TH1::Rebin(nbins_pt_debug, "ppMC_analysisBins_jtPt", ptbins_debug);
    divideBinWidth(theJetIDHist);      
    if(isData)theJetIDHist->Scale(1./theLumi);
    theJetIDHist->Scale(1./etaBinWidth);
  }
  
  else {
    
    theNonJetIDHist->TH1::Rebin(var_xAx_reBin[j]);
    theNonJetIDHist->Scale( 1./theNonJetIDHist->GetBinWidth(1) );          
    if(isData)theNonJetIDHist->Scale(1./theLumi);
    theNonJetIDHist->Scale(1./etaBinWidth);
    
    theJetIDHist->TH1::Rebin(var_xAx_reBin[j]);
    theJetIDHist->Scale( 1./theJetIDHist->GetBinWidth(1) );         
    if(isData)theJetIDHist->Scale(1./theLumi);
    theJetIDHist->Scale(1./etaBinWidth);
    //    theJetIDHist->Scale( theNonJetIDHist->Integral()/theJetIDHist->Integral() );	        
    
  }
  
  
  
  
  // TITLES
  std::string h_XAx_Title=var_xAx_Titles[j], h_YAx_Title=AUAxTitle;
  std::string h_Title="HIN Jet ID Efficiency";//, "+var[j];
  
  theNonJetIDHist->SetTitle (    h_Title.c_str() );
  theNonJetIDHist->SetXTitle( h_XAx_Title.c_str() );
  theNonJetIDHist->SetYTitle( h_YAx_Title.c_str() );
  
  theJetIDHist->SetTitle (    h_Title.c_str() );
  theJetIDHist->SetXTitle( h_XAx_Title.c_str() );
  theJetIDHist->SetYTitle( h_YAx_Title.c_str() );  
  
  
  // STYLE  
  TH1F* theRatio=(TH1F*)theJetIDHist->Clone("MCJetHistClone4Ratio");
  
  dataHistStyle(theNonJetIDHist);
  MCHistStyle(theJetIDHist);
  //ratioHistStyle(theRatio, h_XAx_Title, jetIDratioTitle);
  jetIDratioHistStyle(theRatio, h_XAx_Title, jetIDratioTitle);
  
  if(funcDebug)theNonJetIDHist->Print("base");       std::cout<<std::endl;  
  if(funcDebug)theJetIDHist->Print("base");       std::cout<<std::endl;
  
  if(funcDebug) std::cout<<"creating temporary canvas for printing Jet plots..."<<std::endl;
  TCanvas *temp_canvJet = new TCanvas("tempEvt", "temp Jet Canv withLog", CANVX, CANVY); //This is for the Jet QA Evt plots without logs
  temp_canvJet->cd();
  
  TPad *jetpad1 = new TPad("jetpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  jetpad1->SetLogy(1); 
  jetpad1->SetGridx(1);
  jetpad1->SetGridy(1);
  jetpad1->SetBottomMargin(0);
  
    
  TPad *jetpad2 = new TPad("jetpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
  jetpad2->SetTopMargin(0);
  jetpad2->SetBottomMargin(0.3);
  //jetpad2->SetLogy(0);
  jetpad2->SetLogy(1);
  
  
  if(var[j]=="jtpt"      || var[j]=="rawpt"          ||
     var[j]=="leadJetPt" || var[j]=="subleadJetPt"   ||
     var[j]=="jtpt_forRes" ||     var[j]=="jtpt_L2Res" ||     var[j]=="jtpt_L3Res"      ) {
    
    jetpad1->SetLogx(1);  
    jetpad2->SetLogx(1);  
    
    theJetIDHist->GetXaxis()->SetNoExponent(true);
    theJetIDHist->GetXaxis()->SetMoreLogLabels(true);
    
    theNonJetIDHist->GetXaxis()->SetNoExponent(true);
    theNonJetIDHist->GetXaxis()->SetMoreLogLabels(true);

    theRatio->GetXaxis()->SetNoExponent(true);
    theRatio->GetXaxis()->SetMoreLogLabels(true);
  }
 
  jetpad1->Draw();
  jetpad2->Draw();    


  // variable specific draw orders, canv, etc.
  if(var[j]=="jtpt"||var[j]=="rawpt" ||
     var[j]=="jtpt_forRes" ||     var[j]=="jtpt_L2Res" ||     var[j]=="jtpt_L3Res"      ) {
    
    jetpad1->cd();    
    theJetIDHist->Draw("HIST E"); 
    theNonJetIDHist->Draw("E SAME");  
    
    float t1Loc1=0.54, t1Loc2=0.82;
    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(sampleDescString).c_str());    t1->Draw();    
    TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.05),(fullJetType+"Jets").c_str());    t2->Draw();	    
    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.10),(jetCutString).c_str());    t3->Draw();	    
    TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.15),(jetEtaCutString).c_str());    t4->Draw();	
    
    float legx1=0.73, legx2=legx1+0.11;
    float legy1=0.76, legy2=legy1+0.09;
    
    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->AddEntry(theNonJetIDHist, "no JetID",   "lp");
    theJetQALeg->AddEntry(theJetIDHist,    "with JetID", "lp");
    theJetQALeg->Draw(); 
    
    jetpad2->cd();
    
    theRatio->Divide(theNonJetIDHist); //jetid/nojetid
    fracSubtracted(theRatio); //1-jetid/nojetid


    theRatio->Draw("HIST E");
    
    

    
  }
  
  else if(var[j]=="jteta" ||
	  var[j]=="jtphi"	||
	  ( j>=jetConstits_varStart && j<dijet_varStart )){ 
    
    
    jetpad1->cd();
    theJetIDHist->Draw("HIST E"); 
    theNonJetIDHist->Draw("E same");  
    
    float t1Loc1=0.54, t1Loc2=0.82;
    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(sampleDescString).c_str());    t1->Draw();
    TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.05),(fullJetType+"Jets").c_str());    t2->Draw();	
    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.10),(jetCutString).c_str());    t3->Draw();	
    TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.15),(jetEtaCutString).c_str());    t4->Draw();	
    
    float legx1=0.73, legx2=legx1+0.11;
    float legy1=0.76, legy2=legy1+0.09;
    
    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->AddEntry(theNonJetIDHist, "no JetID",   "lp");
    theJetQALeg->AddEntry(theJetIDHist,    "with JetID", "lp");
    
    theJetQALeg->Draw();
    
    jetpad2->cd();
    
    
    theRatio->Divide(theNonJetIDHist);
    fracSubtracted(theRatio); //1-jetid/nojetid
    theRatio->Draw("HIST E");
    
    
  }
  
  else if(var[j]=="Aj"||
	  var[j]=="xj"||
	  var[j]=="dphi"){ 
    
    jetpad1->cd();
    theJetIDHist->Draw("HIST E"); 
    theNonJetIDHist->Draw("E SAME"); 
    
    float t1Loc1=0.15, t1Loc2=0.7; 	    	  
    if(var[j]=="Aj") t1Loc1=.7;    
    
    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(sampleDescString).c_str());    t1->Draw();    
    TLatex *t2=makeTLatex(t1Loc1,t1Loc2-.04,(fullJetType+"Jets").c_str());    t2->Draw();    
    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.08),"p_{T}^{ld} > 60 GeV");    t3->Draw();		  
    TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.12),"p_{T}^{subld} > 30GeV");    t4->Draw();	    
    TLatex *t5=makeTLatex(t1Loc1, (t1Loc2-0.16), "#Delta#phi > #frac{2}{3}#pi");	        t5->Draw();	

    float legx1=0.15, legx2=legx1+0.11;
    if(var[j]=="Aj") {legx1=.7; legx2=legx1+0.11;}
    float legy1=0.74, legy2=legy1+0.09;
    
    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->AddEntry(theNonJetIDHist, "no JetID",   "lp");
    theJetQALeg->AddEntry(theJetIDHist,    "with JetID", "lp");
    theJetQALeg->Draw(); 
    
    //Ratio Plots
    
    jetpad2->cd();

    
    theRatio->Divide(theNonJetIDHist);
    fracSubtracted(theRatio); //1-jetid/nojetid
    theRatio->Draw("E");
    
    
  }
      
  else if(var[j]=="leadJetPt"||var[j]=="subleadJetPt"){ 
    
    jetpad1->cd();    
    theJetIDHist->Draw("HIST E"); 
    theNonJetIDHist->Draw("E SAME"); 
    
    float t1Loc1=0.7, t1Loc2=0.72; 	    	  

    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(sampleDescString).c_str());    t1->Draw();
    TLatex *t2=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());    t2->Draw();
    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.05),"p_{T}^{ld} > 60 GeV");    t3->Draw();	
    TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.10),"p_{T}^{subld} > 30GeV");    t4->Draw();	
    TLatex *t5=makeTLatex(t1Loc1, (t1Loc2-0.17), "#Delta#phi > #frac{2}{3}#pi");    t5->Draw();	

    float legx1=0.73, legx2=legx1+0.11;
    float legy1=0.76, legy2=legy1+0.09;
    
    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->AddEntry(theNonJetIDHist, "no JetID",   "lp");
    theJetQALeg->AddEntry(theJetIDHist,    "with JetID", "lp");
    theJetQALeg->Draw(); 
    
    //Ratio Plots
    
    jetpad2->cd();
    
    theRatio->Divide(theNonJetIDHist);
    fracSubtracted(theRatio); //1-jetid/nojetid
    theRatio->Draw("HIST E");
    
  }
  
  
  //TLines for ratio plots
  float lowLim, highLim;
  lowLim=theRatio->GetBinLowEdge(1);
  highLim=theRatio->GetBinLowEdge(theRatio->GetNbinsX())+theRatio->GetBinWidth(theRatio->GetNbinsX());
  
  TLine* lineAtOne          = new TLine(lowLim,1.0,highLim,1.0); 
  lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
  TLine* lineAtOneHalf      = new TLine(lowLim,0.5,highLim,0.5); 
  lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
  TLine* lineAtOneEtOneHalf = new TLine(lowLim,1.5,highLim,1.5); 
  lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
  
  //lineAtOne->Draw("same");
  //lineAtOneHalf->Draw("same");
  //lineAtOneEtOneHalf->Draw("same");

  
  // print to PDF file
  temp_canvJet->Print( thePDFFileName.c_str() ); 
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJet->DrawClone();
    
    std::string outcanvtitle="JetID Eff. "+var[j]+" Canvas";
    if(isData)outcanvtitle="Data "+outcanvtitle;
    else outcanvtitle="MC "+outcanvtitle;    
    outcanv->SetTitle(outcanvtitle.c_str());
    
    std::string outcanvname="JetIDEff_"+var[j];
    if(isData)outcanvname="Data_"+outcanvname;
    else outcanvname="MC_"+outcanvname;
    outcanv->Write(outcanvname.c_str());    
  }
  
  temp_canvJet->Close();  
  return;
}







const int rebinfactor=1;//for simpbins aka not-analysis-bins


void printJetTrigHist_wRatio( TFile* fin , bool usedHLT100, bool analysisRebin,
			      std::string thePDFFileName , std::string fullJetType , 
			      std::string trigType, std::string radius , bool usedHLTPF, bool didJetID , TFile* fout=NULL){
  
  if(!fin ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=true;
  //bool analysisRebin=true;
  
  if(trigType=="excl")std::cout<<"drawing exclusive trigger spectra + ratios"<<std::endl;
  else if(trigType=="incl")std::cout<<"drawing inclusive trigger spectra + ratios"<<std::endl;
  else {
    std::cout<<"trigger spectra + ratios not recognized, exit"<<std::endl; 
    return;}
  
  
  
  
  
  
  // PUT JET CUTS SECTION IN HERE
  TH1F* jtptCut_h= (TH1F*)fin->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
  TH1F* jtptQACut_h= (TH1F*)fin->Get( "hJetQAPtCut" );
  std::string jtptQACut_str = std::to_string( (int) jtptQACut_h->GetMean() );
  std::cout<<"jtptQACut_str = "<<jtptQACut_str<<std::endl;
  
  TH1F* jtetaLoCut_h= (TH1F*)fin->Get( "hJetEtaCutLo" );
  std::stringstream etaLo; etaLo.precision(1);
  float jtetaLoCut_F=jtetaLoCut_h->GetMean();
  etaLo << std::fixed << jtetaLoCut_F;
  std::string jtetaLoCut_str = etaLo.str();
  std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;
  
  TH1F* jtetaHiCut_h= (TH1F*)fin->Get( "hJetEtaCutHi" );
  std::stringstream etaHi; etaHi.precision(1);
  float jtetaHiCut_F=jtetaHiCut_h->GetMean();
  etaHi << std::fixed << jtetaHiCut_F;
  std::string jtetaHiCut_str = etaHi.str();
  std::cout<<"jtetaHiCut_str = "<<jtetaHiCut_str<<std::endl;
  
  std::string jetEtaCutString=jtetaLoCut_str+"<|y|<"+jtetaHiCut_str;
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  
  float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  
  std::string sampleDescString=sqrts2k15ppString;
  
  if(funcDebug)std::cout<<std::endl<<"creating temporary canvas for printing JetTrig plots..."<<std::endl;
  TCanvas *temp_canvJetTrig = new TCanvas(("tempJetTrigSpec"+trigType).c_str(), ("blahTrigPt"+trigType).c_str(), 1200, 1000);
  temp_canvJetTrig->cd();

  

  std::cout<<std::endl<<std::endl<<" ------------------------------------------------------------------- "<<std::endl;
  std::cout<<" --------------- DRAWING JET SPECTRA BY TRIGGER --------------- "<<std::endl;
  std::cout<<" ------------------------------------------------------------------- "<<std::endl<<std::endl<<std::endl;
  
  
  TPad *jetpad_excsp = new TPad("jetpad_excsp", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  jetpad_excsp->SetLogx(1);
  jetpad_excsp->SetLogy(1);
  jetpad_excsp->SetGridx(1);
  jetpad_excsp->SetGridy(1);
  
  jetpad_excsp->SetBottomMargin(0);
  jetpad_excsp->Draw();
  
  TPad *jetpad_excrat = new TPad("jetpad_excrat", "Ratio Pad", 0.0, 0.05, 1.0, 0.3);
  jetpad_excrat->SetLogx(1);
  jetpad_excrat->SetLogy(0);
  jetpad_excrat->SetGridx(1);
  jetpad_excrat->SetGridy(0);
  
  jetpad_excrat->SetTopMargin(0);
  jetpad_excrat->SetBottomMargin(0.3);
  
  jetpad_excrat->Draw();
  
  jetpad_excsp->cd();
  
  TLegend* JetTrigLegend=new TLegend(0.70,0.66,0.85,0.86, NULL,"brNDC");
  JetTrigLegend->SetFillStyle(0);
  
  for(int j=0; j<(N_trigs); j++){
    
    //open the hists + do scaling
    if(HLTName[j]=="HLT100"&&!usedHLT100)continue;
    
    //open the hists + do scaling
    std::string inHistName;
    if(trigType=="excl")inHistName="hpp_exc";
    else if(trigType=="incl")inHistName="hpp_";
    
    if(j==0)inHistName="hpp_";
    inHistName+=HLTName[j];
    if(didJetID)inHistName+="_wJetID";
    inHistName+="_"+radius+etaWidth;
    
    if(funcDebug)std::cout<<"inHistName="<<inHistName<<std::endl;
    
    TH1F* theJetTrigQAHist= (TH1F*) ( (TH1*)fin->Get(inHistName.c_str()) );
    if(funcDebug)theJetTrigQAHist->Print("base"); std::cout<<std::endl;
    theJetTrigQAHist->GetXaxis()->SetMoreLogLabels(true);
    theJetTrigQAHist->GetXaxis()->SetNoExponent(true);
    
    if(analysisRebin) {
      theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(nbins_pt_debug,(inHistName+"_spec_anabins").c_str(), ptbins_debug    );
      divideBinWidth(theJetTrigQAHist);
    }
    else {
      //theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(10, (inHistName+"_spec_rebin10").c_str() );
      //theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(1, (inHistName+"_spec_rebin1").c_str() );
      theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(rebinfactor, (inHistName+"_spec_rebin1").c_str() );
      theJetTrigQAHist->Scale(1./((float)rebinfactor));
    }
    theJetTrigQAHist->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");
	  
    
    if(funcDebug)std::cout<<"setting legend entry titles"<<std::endl;
    if(usedHLTPF){
      if(funcDebug)std::cout<<"using HLTPF Trig titles"<<std::endl;
      JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTPFName_Leg[j]).c_str(),"lp");
    }
    else {
      if(funcDebug)std::cout<<"using HLTCalo Trig titles"<<std::endl;
      JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTCaloName_Leg[j]).c_str(),"lp");
    }
    
    if( j==0 ){
      
      std::string h_Title;
      if(trigType=="excl")h_Title="Excl. Jet Spectra by Trigger";
      else if(trigType=="incl")h_Title="Incl. Jet Spectra by Trigger";      
      
      theJetTrigQAHist->SetTitle (    h_Title.c_str() );
      
      std::string h_YAx_Title= AUAxTitle;
      theJetTrigQAHist->SetYTitle( h_YAx_Title.c_str() );      
      
      trigSpectraHistStyle(theJetTrigQAHist, j);
      theJetTrigQAHist->Draw("E");
      
      
      float t1Loc1=0.54, t1Loc2=0.82;
      TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());      t1->Draw();
      TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.04),(jetCutString).c_str());      t2->Draw();
      TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.08),(jetEtaCutString).c_str());      t3->Draw();
      
      float pp1x=t1Loc1,pp1y=t1Loc2-0.12;
      TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str() );      pp1->Draw();
      TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str() );    pp2->Draw();      
    }
    else{
      trigSpectraHistStyle(theJetTrigQAHist, j);
      theJetTrigQAHist->Draw("E SAME");}
    
  }//end loop over trigs for spectra
  
  JetTrigLegend->Draw();
  
  
  
  std::cout<<std::endl<<std::endl<<" ------------------------------------------------------------------ "<<std::endl;
  std::cout<<" --------------- DRAWING JET RATIOS BY TRIGGER --------------- "<<std::endl;
  std::cout<<" ------------------------------------------------------------------ "<<std::endl<<std::endl<<std::endl;
  
  
  jetpad_excrat->cd();
  
  
  TH1F* theDenominator=NULL;
  for(int j=0; j<N_trigs; j++){
    
    //open the hists + do scaling
    if(HLTName[j]=="HLT100"&&!usedHLT100)continue;
    
    std::string inHistName;
    
    if( j==0 ){//combo hlt, the denominator
      
      //open the hists + do scaling
      inHistName="hpp_"+HLTName[j];
      if(didJetID)inHistName+="_wJetID";      
      inHistName+="_"+radius+etaWidth;
      
      if(funcDebug)std::cout<<"inHistName="<<inHistName<<std::endl<<std::endl;
      
      theDenominator = (TH1F*) ( (TH1*)fin->Get(inHistName.c_str()) );
      if(funcDebug)theDenominator->Print("base");
      
      if(analysisRebin) 
	theDenominator=(TH1F*)theDenominator->TH1::Rebin(nbins_pt_debug,(inHistName+"_rat_anabins").c_str(), ptbins_debug    );
      else {
	theDenominator=(TH1F*)theDenominator->TH1::Rebin(rebinfactor, (inHistName+"_rat_rebin10").c_str() );
	theDenominator->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");
      }
      
      trigRatioHistStyle(theDenominator,j);      
      
      
    }
    
    else if(j>0){//individial hlt paths, numerators
      
      if(trigType=="excl")      inHistName="hpp_exc"+HLTName[j];//+"_"+radius+etaWidth;
      else if(trigType=="incl") inHistName="hpp_"+HLTName[j];//+"_"+radius+etaWidth;
      if(didJetID)inHistName+="_wJetID";
      inHistName+="_"+radius+etaWidth;
      
      if(funcDebug)std::cout<<"inHistName="<<inHistName<<std::endl<<std::endl;
      
      TH1F* theJetTrigQAHist= (TH1F*) ( (TH1*)fin->Get(inHistName.c_str()) );
      if(funcDebug)theJetTrigQAHist->Print("base");
      
      if(analysisRebin) 
	theJetTrigQAHist= (TH1F*) theJetTrigQAHist->TH1::Rebin(nbins_pt_debug,(inHistName+"_rat_anabins").c_str(), ptbins_debug    );
      else {
	theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(rebinfactor, (inHistName+"_rat_rebin10").c_str() );
	theJetTrigQAHist->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");}
      
      trigRatioHistStyle(theJetTrigQAHist, j);
      
      theJetTrigQAHist->Divide(theDenominator);
      
      theJetTrigQAHist->SetAxisRange(0.0,2.0,"Y");

      std::string h_XAx_Title="jet p_{T} (GeV)     ";
      theJetTrigQAHist->SetXTitle( h_XAx_Title.c_str() );
      theJetTrigQAHist->GetXaxis()->SetMoreLogLabels(true);
      theJetTrigQAHist->GetXaxis()->SetNoExponent(true);
      
      std::string h_YAx_Title="Combo/Indiv.";
      theJetTrigQAHist->SetYTitle( h_YAx_Title.c_str() );
      
      theJetTrigQAHist->SetTitle("");	
      
      if(j==1)theJetTrigQAHist->Draw("E");
      else theJetTrigQAHist->Draw("E SAME");


      if( j==1 )  {
	
	float min = theJetTrigQAHist->GetBinLowEdge( 1 );
	float max = theJetTrigQAHist->GetBinLowEdge( (int)theJetTrigQAHist->GetNbinsX() ) 
	  + theJetTrigQAHist->GetBinWidth( (int)theJetTrigQAHist->GetNbinsX() ) ;
	
	TLine* lineAtOne          = new TLine(min,1.0,max,1.0);
	lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
	TLine* lineAtOneHalf      = new TLine(min,0.2,max,0.2);
	lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
	TLine* lineAtOneEtOneHalf = new TLine(min,1.8,max,1.8);
	lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
	
	lineAtOne->Draw("same");
	lineAtOneHalf->Draw("same");
	lineAtOneEtOneHalf->Draw("same");
	
      }



      


    }
    
  }//end loop over trigs for ratio
  
  temp_canvJetTrig->Print(thePDFFileName.c_str());
  
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJetTrig->DrawClone();
    
    std::string outcanvtitle = "HLTak4PF"; 
    if(usedHLT100)outcanvtitle+="40 to 100 ";
    else outcanvtitle+="40 to 80 ";
    outcanvtitle+=fullJetType+" Jet p_{T}";    
    if(analysisRebin)outcanvtitle+=" Ana. Bins ";
    else             outcanvtitle+=" Simp. Bins ";
    outcanvtitle +=trigType;

    outcanv->SetTitle((outcanvtitle+" Canvas").c_str());
    
    std::string outcanvname  = "HLTak4PF";
    if(usedHLT100)outcanvname+="40to100_";
    else          outcanvname+="40to80_";
    outcanvname+=fullJetType+"_jtpt_";
    if(analysisRebin) outcanvname+="anabins_";
    else              outcanvname+="simpbins_";
    outcanvname +=trigType;
  
    outcanv->Write((outcanvname+"_canv").c_str());    

  }
  
  //jetpad_excsp->Close();
  //jetpad_excrat->Close();
  temp_canvJetTrig->Close();
  
  
  
  
  return;
}






void printPFvCaloTrigHist_wRatio( TFile* fin , TFile* fin2, bool usedHLT100, bool usedHLTPF, bool analysisRebin, bool comboOnly,
				  std::string thePDFFileName , std::string fullJetType , 
				  std::string trigType, std::string radius , bool didJetID){
  
  if(!fin || !fin2){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  //bool analysisRebin=true;
  
  if(trigType=="excl")std::cout<<"drawing exclusive trigger spectra + ratios"<<std::endl;
  else if(trigType=="incl")std::cout<<"drawing inclusive trigger spectra + ratios"<<std::endl;
  else {
    std::cout<<"trigger spectra + ratios not recognized, exit"<<std::endl; 
    return;}
  
  
  
  
  
  
  // PUT JET CUTS SECTION IN HERE
  TH1F* jtptCut_h= (TH1F*)fin->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
  TH1F* jtptQACut_h= (TH1F*)fin->Get( "hJetQAPtCut" );
  std::string jtptQACut_str = std::to_string( (int) jtptQACut_h->GetMean() );
  std::cout<<"jtptQACut_str = "<<jtptQACut_str<<std::endl;
  
  TH1F* jtetaLoCut_h= (TH1F*)fin->Get( "hJetEtaCutLo" );
  std::stringstream etaLo; etaLo.precision(1);
  float jtetaLoCut_F=jtetaLoCut_h->GetMean();
  etaLo << std::fixed << jtetaLoCut_F;
  std::string jtetaLoCut_str = etaLo.str();
  std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;
  
  TH1F* jtetaHiCut_h= (TH1F*)fin->Get( "hJetEtaCutHi" );
  std::stringstream etaHi; etaHi.precision(1);
  float jtetaHiCut_F=jtetaHiCut_h->GetMean();
  etaHi << std::fixed << jtetaHiCut_F;
  std::string jtetaHiCut_str = etaHi.str();
  std::cout<<"jtetaHiCut_str = "<<jtetaHiCut_str<<std::endl;
  
  std::string jetEtaCutString=jtetaLoCut_str+"<|y|<"+jtetaHiCut_str;
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  
  float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  
  std::string sampleDescString=sqrts2k15ppString;
  
  if(funcDebug)std::cout<<std::endl<<"creating temporary canvas for printing JetTrig plots..."<<std::endl;
  TCanvas *temp_canvJetTrig = new TCanvas(("tempPFvCaloTrigSpec"+trigType).c_str(), ("blahTrigPt"+trigType).c_str(), 1200, 1000);
  temp_canvJetTrig->cd();
  
  
  
  std::cout<<std::endl<<std::endl<<" ------------------------------------------------------------------- "<<std::endl;
  std::cout<<" --------------- DRAWING PFvCalo SPECTRA BY TRIGGER --------------- "<<std::endl;
  std::cout<<" ------------------------------------------------------------------- "<<std::endl<<std::endl<<std::endl;
  
  
  TPad *jetpad_excsp = new TPad("PFvCalopad_excsp", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  jetpad_excsp->SetLogx(1);
  jetpad_excsp->SetLogy(1);
  jetpad_excsp->SetGridx(1);
  jetpad_excsp->SetGridy(1);
  
  jetpad_excsp->SetBottomMargin(0);
  jetpad_excsp->Draw();
  
  TPad *jetpad_excrat = new TPad("PFvCalopad_excrat", "Ratio Pad", 0.0, 0.05, 1.0, 0.3);
  jetpad_excrat->SetLogx(1);
  jetpad_excrat->SetLogy(0);
  jetpad_excrat->SetGridx(1);
  jetpad_excrat->SetGridy(0);
  
  jetpad_excrat->SetTopMargin(0);
  jetpad_excrat->SetBottomMargin(0.3);
  
  jetpad_excrat->Draw();
  
  jetpad_excsp->cd();
  
  TLegend* JetTrigLegend=new TLegend(0.70,0.66,0.85,0.86, NULL,"brNDC");
  JetTrigLegend->SetFillStyle(0);

  for(int j=0; j<(N_trigs); j++){
    
    //open the hists + do scaling
    if(HLTName[j]=="HLT100"&&!usedHLT100)continue;
    
    //open the hists + do scaling
    std::string inHistName;
    if(trigType=="excl")inHistName="hpp_exc";
    else if(trigType=="incl")inHistName="hpp_";
    
    if(j==0)inHistName="hpp_";
    else {if(comboOnly)continue; }//temporary
    inHistName+=HLTName[j];
    if(didJetID)inHistName+="_wJetID";
    inHistName+="_"+radius+etaWidth;
    
    if(funcDebug)std::cout<<"inHistName="<<inHistName<<std::endl;
    
    TH1F* theJetTrigQAHist= (TH1F*) ( (TH1*)fin->Get(inHistName.c_str()) );
    if(funcDebug)theJetTrigQAHist->Print("base"); std::cout<<std::endl;

    TH1F* theJetTrigQAHist2= (TH1F*) ( (TH1*)fin2->Get(inHistName.c_str()) );
    if(funcDebug)theJetTrigQAHist2->Print("base"); std::cout<<std::endl;
    
    if(analysisRebin) {
      theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(nbins_pt_debug,(inHistName+"_spec_anabins").c_str(), ptbins_debug    );
      theJetTrigQAHist2=(TH1F*)theJetTrigQAHist2->TH1::Rebin(nbins_pt_debug,(inHistName+"2_spec_anabins").c_str(), ptbins_debug    );}
    else {
      theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(10, (inHistName+"_spec_rebin10").c_str() );
      theJetTrigQAHist->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");      
      
      theJetTrigQAHist2=(TH1F*)theJetTrigQAHist2->TH1::Rebin(10, (inHistName+"2_spec_rebin10").c_str() );
      theJetTrigQAHist2->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");          }
    divideBinWidth(theJetTrigQAHist);
    divideBinWidth(theJetTrigQAHist2);
    
    trigSpectraHistStyle(theJetTrigQAHist, j);    
    trigSpectra2HistStyle(theJetTrigQAHist2, j);    
    
    if(usedHLTPF){
      JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTPFName_Leg[j]).c_str(),"lp");
      JetTrigLegend->AddEntry(theJetTrigQAHist2,(HLTCaloName_Leg[j]).c_str(),"lp");}
    
    else {
      JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTCaloName_Leg[j]).c_str(),"lp");
      JetTrigLegend->AddEntry(theJetTrigQAHist2,(HLTPFName_Leg[j]).c_str(),"lp");}
    
    if( j==0 ){
      
      std::string h_Title;
      if(trigType=="excl")h_Title="Excl. Jet Spectra by Trigger";
      else if(trigType=="incl")h_Title="Incl. Jet Spectra by Trigger";      
      
      theJetTrigQAHist->SetTitle (    h_Title.c_str() );
      
      std::string h_YAx_Title= AUAxTitle;
      theJetTrigQAHist->SetYTitle( h_YAx_Title.c_str() );
      
      theJetTrigQAHist->Draw("E");      
      theJetTrigQAHist2->Draw("E SAME");      
      
      float t1Loc1=0.54, t1Loc2=0.82;
      TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());      t1->Draw();      
      TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.04),(jetCutString).c_str());      t2->Draw();      
      TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.08),(jetEtaCutString).c_str());      t3->Draw();
      
      float pp1x=t1Loc1,pp1y=t1Loc2-0.12;
      TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str() );      pp1->Draw();      
      TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str() );      pp2->Draw();  
    }
    else{
      trigSpectraHistStyle(theJetTrigQAHist, j);
      theJetTrigQAHist->Draw("E SAME");
      trigSpectra2HistStyle(theJetTrigQAHist2, j);
      theJetTrigQAHist2->Draw("E SAME");
    }
    
  }//end loop over trigs for spectra
  
  JetTrigLegend->Draw();
  
  
  
  std::cout<<std::endl<<std::endl<<" ------------------------------------------------------------------ "<<std::endl;
  std::cout<<" --------------- DRAWING JET RATIOS BY TRIGGER --------------- "<<std::endl;
  std::cout<<" ------------------------------------------------------------------ "<<std::endl<<std::endl<<std::endl;
  
  
  jetpad_excrat->cd();
  
  
  TH1F* theDenominator=NULL;
  TH1F* theNumerator=NULL;
  for(int j=0; j<N_trigs; j++){
    
    //open the hists + do scaling
    if(HLTName[j]=="HLT100"&&!usedHLT100)continue;
    
    std::string inHistName;
    
    
    
    //open the hists + do scaling
    if(trigType=="excl")      inHistName="hpp_exc"+HLTName[j];//+"_"+radius+etaWidth;
    else if(trigType=="incl") inHistName="hpp_"+HLTName[j];//+"_"+radius+etaWidth;
    
    if(j==0)inHistName="hpp_"+HLTName[j];    
    else{ if(comboOnly)continue;} 
    
    if(didJetID)inHistName+="_wJetID";      
    inHistName+="_"+radius+etaWidth;
    
    if(funcDebug)std::cout<<"inHistName="<<inHistName<<std::endl<<std::endl;
    
    
    
    theDenominator = (TH1F*) ( (TH1*)fin->Get(inHistName.c_str()) );
    if(funcDebug)theDenominator->Print("base");
    
    if(analysisRebin) 
      theDenominator=(TH1F*)theDenominator->TH1::Rebin(nbins_pt_debug,(inHistName+"_rat_anabins").c_str(), ptbins_debug    );
    else {
      theDenominator=(TH1F*)theDenominator->TH1::Rebin(10, (inHistName+"_rat_rebin10").c_str() );    
      theDenominator->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");      }
    
    trigRatioHistStyle(theDenominator,j);      
    
    
    

    theNumerator = (TH1F*) ( (TH1*)fin2->Get(inHistName.c_str()) );
    if(funcDebug)theNumerator->Print("base");
    
    if(analysisRebin) 
      theNumerator=(TH1F*)theNumerator->TH1::Rebin(nbins_pt_debug,(inHistName+"2_rat_anabins").c_str(), ptbins_debug    );
    else {
      theNumerator=(TH1F*)theNumerator->TH1::Rebin(10, (inHistName+"2_rat_rebin10").c_str() );    
      theNumerator->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");          
    }

    trigRatioHistStyle2(theNumerator,j);      
    
    
    theNumerator->Divide(theDenominator);
    
    
    
    theNumerator->SetAxisRange(0.0,2.0,"Y");
    
    std::string h_XAx_Title="jet p_{T} (GeV)     ";
    theNumerator->SetXTitle( h_XAx_Title.c_str() );
    theNumerator->GetXaxis()->SetMoreLogLabels(true);
    theNumerator->GetXaxis()->SetNoExponent(true);
    
    std::string h_YAx_Title;
    if(usedHLTPF)h_YAx_Title="Calo/PF";
    else h_YAx_Title="PF/Calo";
    
    
    theNumerator->SetTitle("");	
    theNumerator->SetYTitle( (h_YAx_Title).c_str() );    
    theNumerator->GetYaxis()->SetTitleOffset( .40 );        
    theNumerator->GetYaxis()->SetTitleSize( 0.08 );        
    if(j==0){
      theNumerator->Draw("E");}
    else theNumerator->Draw("E SAME");    



    if( j==0 )  {     
      float min = theNumerator->GetBinLowEdge( 1 );
      float max = theNumerator->GetBinLowEdge( (int)theNumerator->GetNbinsX() ) 
	+ theNumerator->GetBinWidth( (int)theNumerator->GetNbinsX() ) ;
      
      TLine* lineAtOne          = new TLine(min,1.0,max,1.0);
      lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
      TLine* lineAtOneHalf      = new TLine(min,0.2,max,0.2);
      lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
      TLine* lineAtOneEtOneHalf = new TLine(min,1.8,max,1.8);
      lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
      
      lineAtOne->Draw("same");
      lineAtOneHalf->Draw("same");
      lineAtOneEtOneHalf->Draw("same");
      
    }
    
    
  }
  
  
  temp_canvJetTrig->Print(thePDFFileName.c_str());

  
  temp_canvJetTrig->Close();
  
  
  
  
  return;
}





void printTrigPtHist( TFile* fin , bool usedHLT100, bool analysisRebin,
		      std::string thePDFFileName , std::string fullJetType ,
		      std::string trigType, std::string radius , bool usedHLTPF, TFile* fout=NULL){
  
  if(!fin ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(trigType=="excl")std::cout<<"drawing exclusive trigger spectra + ratios"<<std::endl;
  else if(trigType=="incl")std::cout<<"drawing inclusive trigger spectra + ratios"<<std::endl;
  else {
    std::cout<<"trigger spectra + ratios not recognized, exit"<<std::endl; 
    return;}
  
  // PUT JET CUTS SECTION IN HERE
  TH1F* jtptCut_h= (TH1F*)fin->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
  TH1F* jtptQACut_h= (TH1F*)fin->Get( "hJetQAPtCut" );
  std::string jtptQACut_str = std::to_string( (int) jtptQACut_h->GetMean() );
  std::cout<<"jtptQACut_str = "<<jtptQACut_str<<std::endl;
  
  TH1F* jtetaLoCut_h= (TH1F*)fin->Get( "hJetEtaCutLo" );
  std::stringstream etaLo; etaLo.precision(1);
  float jtetaLoCut_F=jtetaLoCut_h->GetMean();
  etaLo << std::fixed << jtetaLoCut_F;
  std::string jtetaLoCut_str = etaLo.str();
  std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;
  
  TH1F* jtetaHiCut_h= (TH1F*)fin->Get( "hJetEtaCutHi" );
  std::stringstream etaHi; etaHi.precision(1);
  float jtetaHiCut_F=jtetaHiCut_h->GetMean();
  etaHi << std::fixed << jtetaHiCut_F;
  std::string jtetaHiCut_str = etaHi.str();
  std::cout<<"jtetaHiCut_str = "<<jtetaHiCut_str<<std::endl;
  
  std::string jetEtaCutString=jtetaLoCut_str+"<|y|<"+jtetaHiCut_str;
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  
  float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  
  std::string sampleDescString=sqrts2k15ppString;

  
  if(funcDebug)std::cout<<std::endl<<"creating temporary canvas for printing JetTrig plots..."<<std::endl;
  TCanvas *temp_canvJetTrig = new TCanvas(("tempJetTrigPt"+trigType).c_str(),("blahTrigPt"+trigType).c_str() , 1200, 1000);
  temp_canvJetTrig->SetLogx(1);
  temp_canvJetTrig->SetLogy(1);
  temp_canvJetTrig->SetGridx(1);
  temp_canvJetTrig->SetGridy(1);
  temp_canvJetTrig->cd();
  
  std::cout<<std::endl<<std::endl<<" ------------------------------------------------------------------- "<<std::endl;
  std::cout<<" --------------- DRAWING TRIGGER JET PT SPECTRA FOR TRIGGERS --------------- "<<std::endl;
  std::cout<<" ------------------------------------------------------------------- "<<std::endl<<std::endl<<std::endl;
  
  
  
  //TH1F* jettrigcomb_clone=NULL;
  
  TLegend* JetTrigLegend=new TLegend(0.70,0.70,0.85,0.86, NULL,"brNDC");
  JetTrigLegend->SetFillStyle(0);
  
  if(funcDebug)     std::cout<<"drawing evt trigPt Spectra, looping over trigs"<<std::endl<<std::endl;
  
  for(int j=0; j<N_trigs; j++){
    
    //open the hists + do scaling
    if(HLTName[j]=="HLT100"&&!usedHLT100)continue;
    //if(HLTName[j]=="HLT80")continue;//debug
    //if(HLTName[j]=="HLT60")continue;//debug
    
    std::string inHistName;
    if(trigType=="excl") inHistName="is"+HLTName[j]+"_trgPt";
    else if (trigType=="incl") inHistName="Inc"+HLTName[j]+"_trgPt";    
    
    if(HLTName[j]=="HLTComb")inHistName=HLTName[j]+"_trgPt";
    
    std::cout<<std::endl<<"inHistName="<<inHistName<<std::endl;
    
    TH1F* theJetTrigQAHist= (TH1F*) ( (TH1*)fin->Get(inHistName.c_str()) );
    theJetTrigQAHist->GetXaxis()->SetNoExponent(true);
    theJetTrigQAHist->GetXaxis()->SetMoreLogLabels(true);
    theJetTrigQAHist->Print("base");
    std::cout<<std::endl<<std::endl;
    
    if(analysisRebin){
      theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(nbins_pt_debug,(inHistName+"_spec_anabins").c_str(), ptbins_debug    );
      divideBinWidth(theJetTrigQAHist);
    }
    else{
      theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(rebinfactor, (inHistName+"_rebin").c_str());
      theJetTrigQAHist->Scale(1./((float)rebinfactor));
      //theJetTrigQAHist->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");          
      theJetTrigQAHist->SetAxisRange(40,1400,"X");          
    }
    
    
    trigPtHistStyle(theJetTrigQAHist,j);
    
    if(usedHLTPF) JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTPFName_Leg[j]).c_str(),"lp");
    else JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTCaloName_Leg[j]).c_str(),"lp");
    
    //if(j==0){
    if(HLTName[j]=="HLTComb"){
      
      std::string h_Title;
      if(trigType=="excl")h_Title="is[40,60,80] Event Trigger Jet Pt";
      else if(trigType=="incl")h_Title="Trigger Jet Pt for HLT[40,60,80]";
      
      
      theJetTrigQAHist->SetTitle (    h_Title.c_str() );
      theJetTrigQAHist->SetYTitle( (AUAxTitle).c_str() );
      theJetTrigQAHist->SetXTitle( "trigger jet p_{T} (GeV)   ");
      
      //theJetTrigQAHist->SetAxisRange(0.1,100000000, "Y");
      
      theJetTrigQAHist->Draw("E");
      //jettrigcomb_clone=(TH1F*)theJetTrigQAHist->Clone("HLTComb_cloneHist");
      
      
      float t1Loc1=0.54, t1Loc2=0.82;
      
      TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());
      t1->SetTextFont(63);
      t1->SetTextColor(kBlack);
      t1->SetTextSize(14);
      t1->SetLineWidth(1);
      t1->SetNDC();
      //t1->Draw();

      TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.04),("0.0< #||{y} <5.1"));
      t2->SetTextFont(63);
      t2->SetTextColor(kBlack);
      t2->SetTextSize(14);
      t2->SetLineWidth(1);
      t2->SetNDC();
      t2->Draw();

      float pp1x=t1Loc1,pp1y=t1Loc2-0.08;
      TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str() );
      pp1->SetTextFont(63);
      pp1->SetTextColor(kBlack);
      pp1->SetTextSize(15);
      pp1->SetLineWidth(1);
      pp1->SetNDC();
      pp1->Draw();
      
      TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str() );
      pp2->SetTextFont(63);
      pp2->SetTextColor(kBlack);
      pp2->SetTextSize(15);
      pp2->SetLineWidth(1);
      pp2->SetNDC();
      pp2->Draw();  
    }//if HLTComb    
    else {
      theJetTrigQAHist->Draw("E SAME");
    }//if not HLTComb
    
    
    
  }//end loop over trigs for spectra
  
  //jettrigcomb_clone->Draw("E SAME");
  JetTrigLegend->Draw();
  




  temp_canvJetTrig->Print(thePDFFileName.c_str());
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJetTrig->DrawClone();

    std::string outcanvtitle = "HLTak4PF"; 
    if(usedHLT100)outcanvtitle+="40 to 100 ";
    else outcanvtitle+="40 to 80 ";
    outcanvtitle+=fullJetType+" Trig. Jet p_{T}";    
    if(analysisRebin)outcanvtitle+=" Ana. Bins ";
    else             outcanvtitle+=" Simp. Bins ";
    outcanvtitle +=trigType;

    outcanv->SetTitle((outcanvtitle+" Canvas").c_str());
    
    std::string outcanvname  = "HLTak4PF";
    if(usedHLT100)outcanvname+="40to100_";
    else          outcanvname+="40to80_";
    outcanvname+=fullJetType+"_trgjtpt_";
    if(analysisRebin) outcanvname+="anabins_";
    else              outcanvname+="simpbins_";
    outcanvname +=trigType;

    outcanv->Write((outcanvname+"_canv").c_str());    

    //std::string outcanvtitle="HLTak4PFJet",outcanvname="HLTak4PFJet";//= "HLT"fullJetType+
    //if(usedHLT100)outcanvname+="40to100_";
    //else outcanvname+="40to80_";
    //if(usedHLT100)outcanvtitle+="40 to 100 ";
    //else outcanvtitle+="40 to 80 ";
    //outcanvname+=fullJetType+"Jets_";
    //outcanvtitle+=fullJetType+"Jets ";    
    //outcanv->SetTitle((outcanvtitle+" Canvas").c_str());
    //outcanv->Write((outcanvname+"canv").c_str());    
  }  
  temp_canvJetTrig->Close();

  return;
}







void printTrigEtaHist( TFile* fin , bool usedHLT100,
		       std::string thePDFFileName , std::string fullJetType ,
		       std::string trigType, std::string radius , bool usedHLTPF, TFile* fout=NULL){
  
  if(!fin ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(trigType=="excl")std::cout<<"drawing exclusive trigger spectra + ratios"<<std::endl;
  else if(trigType=="incl")std::cout<<"drawing inclusive trigger spectra + ratios"<<std::endl;
  else {
    std::cout<<"trigger spectra + ratios not recognized, exit"<<std::endl; 
    return;}
  
  // PUT JET CUTS SECTION IN HERE
  TH1F* jtptCut_h= (TH1F*)fin->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
  TH1F* jtptQACut_h= (TH1F*)fin->Get( "hJetQAPtCut" );
  std::string jtptQACut_str = std::to_string( (int) jtptQACut_h->GetMean() );
  std::cout<<"jtptQACut_str = "<<jtptQACut_str<<std::endl;
  
  TH1F* jtetaLoCut_h= (TH1F*)fin->Get( "hJetEtaCutLo" );
  std::stringstream etaLo; etaLo.precision(1);
  float jtetaLoCut_F=jtetaLoCut_h->GetMean();
  etaLo << std::fixed << jtetaLoCut_F;
  std::string jtetaLoCut_str = etaLo.str();
  std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;
  
  TH1F* jtetaHiCut_h= (TH1F*)fin->Get( "hJetEtaCutHi" );
  std::stringstream etaHi; etaHi.precision(1);
  float jtetaHiCut_F=jtetaHiCut_h->GetMean();
  etaHi << std::fixed << jtetaHiCut_F;
  std::string jtetaHiCut_str = etaHi.str();
  std::cout<<"jtetaHiCut_str = "<<jtetaHiCut_str<<std::endl;
  
  std::string jetEtaCutString=jtetaLoCut_str+"<|y|<"+jtetaHiCut_str;
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  
  float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  
  std::string sampleDescString=sqrts2k15ppString;
  
  if(funcDebug)std::cout<<std::endl<<"creating temporary canvas for printing JetTrig plots..."<<std::endl;
  
  
  
  
  
  TCanvas *temp_canvJetTrig = new TCanvas(("tempJetTrigEta"+trigType).c_str(),("blahTrigEta"+trigType).c_str() , 1200, 1000);
  temp_canvJetTrig->SetLogx(0);
  temp_canvJetTrig->SetLogy(1);
  temp_canvJetTrig->SetGridx(1);
  temp_canvJetTrig->SetGridy(1);
  temp_canvJetTrig->cd();

  

  std::cout<<std::endl<<std::endl<<" ------------------------------------------------------------------- "<<std::endl;
  std::cout<<" --------------- DRAWING TRIGGER JET ETA SPECTRA FOR TRIGGERS --------------- "<<std::endl;
  std::cout<<" ------------------------------------------------------------------- "<<std::endl<<std::endl<<std::endl;
  
  
  
  
  
  TLegend* JetTrigLegend=new TLegend(0.50,0.35,0.65,0.46, NULL,"brNDC");
  JetTrigLegend->SetFillStyle(0);
  
  if(funcDebug)     std::cout<<"drawing evt trigEta Spectra, looping over trigs"<<std::endl<<std::endl;
  
  for(int j=0; j<N_trigs; j++){
    
    //open the hists + do scaling
    if(HLTName[j]=="HLT100"&&!usedHLT100)continue;
    
    std::string inHistName;
    if(trigType=="excl") inHistName="is"+HLTName[j]+"_trgEta";
    else if (trigType=="incl") inHistName="Inc"+HLTName[j]+"_trgEta";
    if(j==0)inHistName=HLTName[j]+"_trgEta";
    
    if(funcDebug)std::cout<<std::endl<<"inHistName="<<inHistName<<std::endl;
    
    TH1F* theJetTrigQAHist= (TH1F*) ( (TH1*)fin->Get(inHistName.c_str()) );
    if(funcDebug)theJetTrigQAHist->Print("base");
    std::cout<<std::endl<<std::endl;
    
    
    theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(2, (inHistName+"_rebin").c_str());
    theJetTrigQAHist->Scale(1./2.);    
    


    trigEtaHistStyle(theJetTrigQAHist,j);
    
    theJetTrigQAHist->SetAxisRange(800,40000000, "Y");

    if(usedHLTPF) JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTPFName_Leg[j]).c_str(),"lp");
    else JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTCaloName_Leg[j]).c_str(),"lp");
    
    if(j==0){
      
      std::string h_Title;//   ="Exclusive Event Trigger Jet Eta";//, Data ONLY";
      if(trigType=="excl")h_Title="is[40,60,80] Event Trigger Jet Eta";
      else if(trigType=="incl")h_Title="Trigger Jet Eta for HLT[40,60,80]";
      
      
      theJetTrigQAHist->SetTitle (    h_Title.c_str() );
      theJetTrigQAHist->SetYTitle( (AUAxTitle).c_str() );
      theJetTrigQAHist->SetXTitle( "trigger jet #||{y}   ");

      theJetTrigQAHist->Draw("E");



      float t1Loc1=0.35, t1Loc2=0.50;
      //TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());      //t1->Draw();      
      TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.04),("0.0< #||{y} <5.1"));      t2->Draw();
      
      float pp1x=t1Loc1,pp1y=t1Loc2-0.08;
      TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str() );      pp1->Draw();
      TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str() );      pp2->Draw();  
    }
    
    else
      theJetTrigQAHist->Draw("E SAME");
    
    
    
  }//end loop over trigs for spectra
  
  JetTrigLegend->Draw();
  
  temp_canvJetTrig->Print(thePDFFileName.c_str());
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJetTrig->DrawClone();
    
    std::string outcanvtitle = "HLTak4PF"; 
    if(usedHLT100)outcanvtitle+="40 to 100 ";
    else outcanvtitle+="40 to 80 ";
    outcanvtitle+=fullJetType+" Trig. Jet #eta ";    
    //if(analysisRebin)outcanvtitle+=" Ana. Bins ";
    //else             outcanvtitle+=" Simp. Bins ";
    outcanvtitle +=trigType;

    outcanv->SetTitle((outcanvtitle+" Canvas").c_str());
    
    std::string outcanvname  = "HLTak4PF";
    if(usedHLT100)outcanvname+="40to100_";
    else          outcanvname+="40to80_";
    outcanvname+=fullJetType+"_trgjteta_";
    //if(analysisRebin) outcanvname+="anabins_";
    //else              outcanvname+="simpbins_";
    outcanvname +=trigType;
    
    outcanv->Write((outcanvname+"_canv").c_str());    

  }

  temp_canvJetTrig->Close();

  return;
}





TH1F* makeAsymmEtaHist(TH1F* theEtaHist,std::string trigType) {
  //  bool funcDebug=false;
  std::string hTitle=theEtaHist->GetName();
  int NBinsNewHist=theEtaHist->GetNbinsX()/2;
  float min = 0. ; // = theEtaHist->GetBinLowEdge( 1 );   
  float max = theEtaHist->GetBinLowEdge( (int)theEtaHist->GetNbinsX() )
    + theEtaHist->GetBinWidth( (int)theEtaHist->GetNbinsX() ) ;
  
  if(theEtaHist->GetNbinsX()%2 != 0) 
    std::cout<<std::endl<<"WARNING: AsymmEtaHist will have odd # of bins! Will not look good!"<<std::endl<<std::endl;
  
  TH1F* newAsymmEtaHist= new TH1F( (hTitle+"_asymmHist_"+trigType).c_str(), "",   NBinsNewHist, min,max);
  //newAsymmEtaHist->Print("base");
  
  for(int i = 1; i<=NBinsNewHist; i++){
    

    float negSideContent, negSideErr;
    float posiSideContent, posiSideErr;
    
    negSideContent=theEtaHist->GetBinContent(NBinsNewHist-i+1);
    negSideErr=theEtaHist      ->GetBinError(NBinsNewHist-i+1);
    
    posiSideContent=theEtaHist->GetBinContent( NBinsNewHist+i );
    posiSideErr=theEtaHist      ->GetBinError( NBinsNewHist+i );
    
    if(i==(NBinsNewHist*2-i+1) )continue;
    
    float asymmBinContent,asymmBinErr; 
    asymmBinContent=posiSideContent/negSideContent;
    asymmBinErr=
      (posiSideContent/negSideContent)* std::sqrt( 
						  (negSideErr/negSideContent)*(negSideErr/negSideContent) +
						  (posiSideErr/posiSideContent)*(posiSideErr/posiSideContent)    
						   );   
    
        
    if( asymmBinContent != asymmBinContent ||
	asymmBinErr != asymmBinErr
	){
      continue;
    }
    else{

      newAsymmEtaHist->SetBinContent( i, asymmBinContent );
      newAsymmEtaHist->SetBinError(   i, asymmBinErr   );
    }
    
    
    
    
  }
  
  
  return newAsymmEtaHist;
  
}




void printTrigEtaAsymmHist( TFile* fin , bool usedHLT100,
			    std::string thePDFFileName , std::string fullJetType ,
			    std::string trigType, std::string radius , bool usedHLTPF, TFile* fout=NULL){
  
  if(!fin ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(trigType=="excl")std::cout<<"drawing exclusive trigger spectra + ratios"<<std::endl;
  else if(trigType=="incl")std::cout<<"drawing inclusive trigger spectra + ratios"<<std::endl;
  else {
    std::cout<<"trigger spectra + ratios not recognized, exit"<<std::endl; 
    return;}
  
  // PUT JET CUTS SECTION IN HERE
  TH1F* jtptCut_h= (TH1F*)fin->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
  TH1F* jtptQACut_h= (TH1F*)fin->Get( "hJetQAPtCut" );
  std::string jtptQACut_str = std::to_string( (int) jtptQACut_h->GetMean() );
  std::cout<<"jtptQACut_str = "<<jtptQACut_str<<std::endl;
  
  TH1F* jtetaLoCut_h= (TH1F*)fin->Get( "hJetEtaCutLo" );
  std::stringstream etaLo; etaLo.precision(1);
  float jtetaLoCut_F=jtetaLoCut_h->GetMean();
  etaLo << std::fixed << jtetaLoCut_F;
  std::string jtetaLoCut_str = etaLo.str();
  std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;
  
  TH1F* jtetaHiCut_h= (TH1F*)fin->Get( "hJetEtaCutHi" );
  std::stringstream etaHi; etaHi.precision(1);
  float jtetaHiCut_F=jtetaHiCut_h->GetMean();
  etaHi << std::fixed << jtetaHiCut_F;
  std::string jtetaHiCut_str = etaHi.str();
  std::cout<<"jtetaHiCut_str = "<<jtetaHiCut_str<<std::endl;
  
  std::string jetEtaCutString=jtetaLoCut_str+"<|y|<"+jtetaHiCut_str;
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  
  float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  
  std::string sampleDescString=sqrts2k15ppString;
  
  if(funcDebug)std::cout<<std::endl<<"creating temporary canvas for printing JetTrig plots..."<<std::endl;
  
  
  
  
  
  TCanvas *temp_canvJetTrig = new TCanvas(("tempJetTrigAsymmEta"+trigType).c_str(),("blahTrigAsymmEta"+trigType).c_str() , 1200, 1000);
  temp_canvJetTrig->SetLogx(0);
  temp_canvJetTrig->SetLogy(0);
  temp_canvJetTrig->SetGridx(1);
  temp_canvJetTrig->SetGridy(1);
  temp_canvJetTrig->cd();

  
  
  std::cout<<std::endl<<std::endl<<" ------------------------------------------------------------------- "<<std::endl;
  std::cout<<" --------------- DRAWING TRIGGER JET ASYMM ETA HISTS FOR TRIGGERS --------------- "<<std::endl;
  std::cout<<" ------------------------------------------------------------------- "<<std::endl<<std::endl<<std::endl;
  
  
  
  
  
  TLegend* JetTrigLegend=new TLegend(0.70,0.70,0.85,0.86, NULL,"brNDC");
  JetTrigLegend->SetFillStyle(0);
  
  if(funcDebug)     std::cout<<"drawing evt trigEta AsymmSpectra, looping over trigs"<<std::endl<<std::endl;
  
  for(int j=0; j<N_trigs; j++){
    
    //open the hists + do scaling
    if(HLTName[j]=="HLT100"&&!usedHLT100)continue;
    
    std::string inHistName;
    if(trigType=="excl") inHistName="is"+HLTName[j]+"_trgEta";
    else if (trigType=="incl") inHistName="Inc"+HLTName[j]+"_trgEta";    
    if(j==0)inHistName=HLTName[j]+"_trgEta";
    if(funcDebug)std::cout<<std::endl<<"inHistName="<<inHistName<<std::endl;
    
    TH1F* theJetTrigQAHist= (TH1F*) ( (TH1*)fin->Get(inHistName.c_str()) );
    if(funcDebug)theJetTrigQAHist->Print("base");
    std::cout<<std::endl<<std::endl;
    
    
    theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(10, (inHistName+"_rebin").c_str());
    theJetTrigQAHist->Scale(1./10.);    
        
    
    
    
    TH1F* theJetTrigAsymmQAHist= makeAsymmEtaHist(theJetTrigQAHist,trigType) ; //(TH1*)fin->Get(inHistName.c_str()) );
    if(funcDebug)theJetTrigAsymmQAHist->Print("base");
    trigEtaHistStyle(theJetTrigAsymmQAHist,j);    
    
    if(usedHLTPF) JetTrigLegend->AddEntry(theJetTrigAsymmQAHist,(HLTPFName_Leg[j]).c_str(),"lp");
    else JetTrigLegend->AddEntry(theJetTrigAsymmQAHist,(HLTCaloName_Leg[j]).c_str(),"lp");
    
    theJetTrigAsymmQAHist->SetAxisRange(0.0,2.0,"Y");
    theJetTrigAsymmQAHist->SetAxisRange(0.0,5.0,"X");
    if(j==0){
      
      std::string h_Title;//   ="Exclusive Event Trigger Jet Eta";//, Data ONLY";
      if(trigType=="excl")h_Title="is[40,60,80] Event Trigger Jet Eta Asymm";
      else if(trigType=="incl")h_Title="Trigger Jet Eta Asymm for HLT[40,60,80]";
      
      
      theJetTrigAsymmQAHist->SetTitle (    h_Title.c_str() );
      theJetTrigAsymmQAHist->SetYTitle( "+y bin / -y bin" );
      theJetTrigAsymmQAHist->SetXTitle( "trigger jet #||{y}  ");
      
      //theJetTrigAsymmQAHist->SetAxisRange(0.1,100000000, "Y");
      theJetTrigAsymmQAHist->Draw("E");
      
      float t1Loc1=0.54, t1Loc2=0.82;
      //TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());      //t1->Draw();
      TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.04),("0.0< #||{y} <5.1"));      t2->Draw();
      
      float pp1x=t1Loc1,pp1y=t1Loc2-0.08;
      TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str() );       pp1->Draw();
      TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str() );   pp2->Draw();  
    }    
    else
      theJetTrigAsymmQAHist->Draw("E SAME");    

  }//end loop over trigs for spectra
  
  JetTrigLegend->Draw();  

  temp_canvJetTrig->Print(thePDFFileName.c_str());  
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJetTrig->DrawClone();
    
    std::string outcanvtitle = "HLTak4PF"; 
    if(usedHLT100)outcanvtitle+="40 to 100 ";
    else outcanvtitle+="40 to 80 ";
    outcanvtitle+=fullJetType+" Trig. Jet #eta Asymm. ";    
//    if(analysisRebin)outcanvtitle+=" Ana. Bins ";
//    else             outcanvtitle+=" Simp. Bins ";
    outcanvtitle +=trigType;

    outcanv->SetTitle((outcanvtitle+" Canvas").c_str());
    
    std::string outcanvname  = "HLTak4PF";
    if(usedHLT100)outcanvname+="40to100_";
    else          outcanvname+="40to80_";
    outcanvname+=fullJetType+"_trgjteta_asymm_";
//    if(analysisRebin) outcanvname+="anabins_";
//    else              outcanvname+="simpbins_";
    outcanvname +=trigType;
    
    outcanv->Write((outcanvname+"_canv").c_str());    


//    std::string outcanvtitle="HLTak4PFJet",outcanvname="HLTak4PFJet";//= "HLT"fullJetType+
//    if(usedHLT100)outcanvname+="40to100_";
//    else outcanvname+="40to80_";
//    if(usedHLT100)outcanvtitle+="40 to 100 ";
//    else outcanvtitle+="40 to 80 ";
//    outcanvname+=fullJetType+"Jets_";
//    outcanvtitle+=fullJetType+"Jets ";    
//    outcanv->SetTitle((outcanvtitle+" Canvas").c_str());
//    outcanv->Write((outcanvname+"canv").c_str());    
  }
  temp_canvJetTrig->Close();
  return;

}



void printJetSpectraRapHists( TFile* finData , TFile* finMC, bool doJetIDPlots,
			      std::string thePDFFileName , std::string fullJetType, 
			      long double theLumi  ) {
  bool funcDebug=true;
  bool drawMC=(bool)finMC;
  //if( ((bool)finData) && !((bool)finMC) )
  //  std::cout<<"drawing data jetspectraraphists only"<<std::endl;
  //else if( ((bool)finData) && ((bool)finMC) ) 
  //  std::cout<<"drawing data + MC RECO jetspectraraphists only"<<std::endl;
  //else if( !((bool)finData) && ((bool)finMC) )
  //  std::cout<<"drawing MC RECO + GEN jetspectraraphists only"<<std::endl;
  
  
  if(funcDebug)std::cout<<std::endl<<"creating temporary canvas for printing Jet rapBins plots..."<<std::endl;
  TCanvas *temp_canvJetRapBins = new TCanvas("tempJetRapBins", "tempJetRapBins", 2000, 1500);
  temp_canvJetRapBins->SetLogy(1);    
  temp_canvJetRapBins->SetLogx(1);    
  
  float TLegx1=0.77 ,TLegx2=TLegx1+0.13;
  float TLegy1=0.54, TLegy2=TLegy1+0.36;
  TLegend* theJetRapHistLeg=new TLegend(TLegx1, TLegy1, TLegx2, TLegy2, NULL,"brNDC");
  
  temp_canvJetRapBins->cd();
  
  std::cout<<"drawing dual-differential cross section plot"<<std::endl;	
  std::cout<<"looping over "<<nbins_abseta<<" rapidity bins"<<std::endl;      
  
  float power=std::pow(10., 7.);//std::pow(10., 8.);
  int bincounter=8;
  float totalIntXsec=0.0;
  
  for(int rapbin=0; rapbin<8; ++rapbin){	
    
    
    if(rapbin!=6){
      power/=10.;
      bincounter=bincounter-1;
      //if(rapbin>5)continue;
    }
    
    else{  std::cout<<"rapbin==6, continuing on..."<<std::endl;
      continue;  }
    
    //open the hists + rebin
    std::string inHistName="hJetSpectraRap";	
    if(doJetIDPlots) inHistName+="_wJetID"; 
    inHistName+="_bin"+std::to_string(rapbin);
    if(funcDebug)std::cout<<"inHistName="<<inHistName<<std::endl;	
    
    
    
    TH1F* theJetSpectraRapHist= (TH1F*) finData->Get( inHistName.c_str() );
    theJetSpectraRapHist=(TH1F*)theJetSpectraRapHist->
      TH1::Rebin(nbins_pt_debug, (inHistName+"_rebin").c_str(), ptbins_debug);
    divideBinWidth(theJetSpectraRapHist);
    if(rapbin==7) theJetSpectraRapHist->Scale(1/3.); //3.2 - 4.7
    theJetSpectraRapHist->Scale( 1./theLumi);
    theJetSpectraRapHist->Scale( power );
    
    float integral=theJetSpectraRapHist->Integral(); //correct this, integral doesnt account for bin width
    std::cout<<"integral="<<integral/power<<std::endl;
    totalIntXsec+=(integral/power);
    std::cout<<"totalIntXsec="<<totalIntXsec<<std::endl;
    
    
    theJetSpectraRapHist->SetMarkerStyle( kFullDotLarge ) ;
    theJetSpectraRapHist->SetMarkerColor( theRapOverlayMarkerColor[rapbin] );
    theJetSpectraRapHist->SetLineColor( kBlack );
    
    
    // MC Hist
    std::string inMCHistName="hJetSpectraRap";	
    if(doJetIDPlots) inMCHistName+="_wJetID"; 
    inMCHistName+="_bin"+std::to_string(rapbin);
    if(funcDebug)std::cout<<"inMCHistName="<<inMCHistName<<std::endl;	
    
    
    TH1F* theMCJetSpectraRapHist=NULL;
    if(drawMC){
      theMCJetSpectraRapHist= (TH1F*)finMC->Get(inMCHistName.c_str()) ;
      theMCJetSpectraRapHist=(TH1F*)theMCJetSpectraRapHist->
	TH1::Rebin(nbins_pt_debug, (inMCHistName+"_rebin").c_str(), ptbins_debug);
      divideBinWidth(theMCJetSpectraRapHist);
      if(rapbin==7) theMCJetSpectraRapHist->Scale(1./3.); //3.2 - 4.7
      theMCJetSpectraRapHist->Scale( 1000. ); 
      theMCJetSpectraRapHist->Scale( power ); //i think the MC pthat weights are in units of nanobarns? multiply by 1000 
      //theMCJetSpectraRapHist->Scale( theJetSpectraRapHist->Integral()/theMCJetSpectraRapHist->Integral() );
      
      theMCJetSpectraRapHist->SetLineColor( theRapOverlayMarkerColor[rapbin] );
    }
    
    std::stringstream stream1, stream2;
    stream1.precision(1); stream1 << std::fixed << absetabins[rapbin];
    stream2.precision(1); stream2 << std::fixed << absetabins[rapbin+1];
    
    if(rapbin==0){
      //      std::string rapHistLegDesc="(|y| < "+stream2.str()+") x 10^{"+std::to_string(bincounter)+"}";	
      std::string rapHistLegDesc="(|y| < "+stream2.str()+") x 10^{"+std::to_string(bincounter-1)+"}";	
      theJetRapHistLeg->AddEntry(theJetSpectraRapHist,(rapHistLegDesc).c_str(),"lp");
      if(drawMC)theJetRapHistLeg->AddEntry(theMCJetSpectraRapHist,"PY8 QCD MC","l");
    }
    else{
      //std::string rapHistLegDesc="("+stream1.str()+" < |y| <"+stream2.str()+") x10^{"+std::to_string(bincounter)+"}";	
      std::string rapHistLegDesc="("+stream1.str()+" < |y| <"+stream2.str()+") x10^{"+std::to_string(bincounter-1)+"}";	
      theJetRapHistLeg->AddEntry(theJetSpectraRapHist,(rapHistLegDesc).c_str(),"lp");
      if(drawMC)theJetRapHistLeg->AddEntry(theMCJetSpectraRapHist,"PY8 QCD MC","l");
    }
    
    if(rapbin==0){
      
      std::string h_Title   ="Jet Pt in |y| bins";
      if(doJetIDPlots) h_Title+=", w/ JetIDCut";
      theJetSpectraRapHist->SetTitle (    h_Title.c_str() );
      
      std::string h_XAx_Title="Jet p^{RECO}_{T} (GeV)", h_YAx_Title="d^{2}#sigma / d|y| dp_{T} (#mu b)";//crossSectionAxTitle;
      theJetSpectraRapHist->SetXTitle( h_XAx_Title.c_str() );
      theJetSpectraRapHist->SetYTitle( h_YAx_Title.c_str() );
      theJetSpectraRapHist->SetAxisRange(std::pow(10.,-10.),std::pow(10.,5.5),"Y");
      //theJetSpectraRapHist->SetAxisRange(std::pow(10.,-9.),std::pow(10.,6.5),"Y");
      //theJetSpectraRapHist->SetAxisRange(60.,840.,"X");
      
      theJetSpectraRapHist->GetXaxis()->SetMoreLogLabels(true);
      theJetSpectraRapHist->GetXaxis()->SetNoExponent(true);
      theJetSpectraRapHist->Draw("E");     	  
      if(drawMC)theMCJetSpectraRapHist->Draw("HIST E SAME");     	  
      theJetRapHistLeg->Draw(); 
    }
    else { 
      theJetSpectraRapHist->Draw("E SAME");
      if(drawMC)theMCJetSpectraRapHist->Draw("HIST E SAME");
      }
    
  }//end rapbin hist loop



  // STRINGS FOR HISTS

  TH1F* jtptCut_h= (TH1F*)finData->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}^{RECO} > "+jtptCut_str+" GeV";
  std::cout<<"jetCutString="<<jetCutString<<std::endl;

  
  std::cout<<"FINAL TOTAL INT JET XSEC="<<totalIntXsec<<std::endl;
  float t1Loc1=0.582, t1Loc2=0.86, deltaLoc=0.025;
  TLatex *pp1=makeTLatex(t1Loc1,t1Loc2,"pp 2015 promptReco, #sqrt{s} = 5.02 TeV");  pp1->Draw();  
  TLatex *pp2=makeTLatex(t1Loc1,(t1Loc2-deltaLoc),"Jet80+LowerJets, L_{int} = 27.4 pb^{-1}");  pp2->Draw();
  
  TLatex *t1=makeTLatex(t1Loc1,(t1Loc2-2*deltaLoc),(fullJetType+"Jets").c_str());    t1->Draw();
  TLatex *t2=makeTLatex(t1Loc1,(t1Loc2-3*deltaLoc),(jetCutString).c_str());    t2->Draw();
  if(drawMC)    { 
    TH1F* genjtptCut_h= (TH1F*)finMC->Get( "hGenJetPtCut" );
    std::string genjtptCut_str = std::to_string( (int) genjtptCut_h->GetMean() );
    std::cout<<"genjtptCut_str = "<<genjtptCut_str<<std::endl;
    
    std::string genJetCutString="p_{T}^{GEN} > "+genjtptCut_str+" GeV";
    std::cout<<"genJetCutString="<<genJetCutString<<std::endl;    
    
    TLatex *t3=makeTLatex(t1Loc1,(t1Loc2-4*deltaLoc),(genJetCutString).c_str());    
    t3->Draw(); 
  }
  if(doJetIDPlots)    { 
      TLatex *t4=makeTLatex(t1Loc1,(t1Loc2-5*deltaLoc),"Jet ID Applied");      t4->Draw();	  
  }

  temp_canvJetRapBins->Print( thePDFFileName.c_str() );
  temp_canvJetRapBins->Close();  
  return;
}







void printMCJetSpectraRapHists(  TFile* finMC, bool doJetIDPlots,
			      std::string thePDFFileName , std::string fullJetType  ) {
  bool funcDebug=true;
  
  
  if(funcDebug)std::cout<<std::endl<<"creating temporary canvas for printing Jet rapBins plots..."<<std::endl;
  TCanvas *temp_canvJetRapBins = new TCanvas("tempJetRapBins", "tempJetRapBins", 2000, 1500);
  temp_canvJetRapBins->SetLogy(1);    
  temp_canvJetRapBins->SetLogx(1);    
  
  TLegend* theJetRapHistLeg=new TLegend(0.72,0.60,0.85,0.86, NULL,"brNDC");
  
  temp_canvJetRapBins->cd();
  
  std::cout<<"drawing dual-differential cross section plot"<<std::endl;	
  std::cout<<"looping over "<<nbins_abseta<<" rapidity bins"<<std::endl;      
  
  float power=std::pow(10., 8.);
  int bincounter=8;
  float totalIntXsec=0.0;
  
  for(int rapbin=0; rapbin<8; ++rapbin){	
    
    
    if(rapbin!=6){
      power/=10;
      bincounter=bincounter-1;
      if(rapbin>5)continue;
    }
    
    else{  std::cout<<"rapbin==6, continuing on..."<<std::endl;
      continue;  }
    
    //open the hists + rebin
    std::string inHistName="hJetSpectraRap";	
    if(doJetIDPlots) inHistName+="_wJetID"; 
    inHistName+="_bin"+std::to_string(rapbin);
    if(funcDebug)std::cout<<"inHistName="<<inHistName<<std::endl;	
    
        
    TH1F* theJetSpectraRapHist= (TH1F*) finMC->Get( inHistName.c_str() );
    theJetSpectraRapHist=(TH1F*)theJetSpectraRapHist->
      TH1::Rebin(nbins_pt_debug, (inHistName+"_rebin").c_str(), ptbins_debug);
    divideBinWidth(theJetSpectraRapHist);
    if(rapbin==7) theJetSpectraRapHist->Scale(1/3.); //3.2 - 4.7
    theJetSpectraRapHist->Scale( power*1000 );
    
    float integral=theJetSpectraRapHist->Integral(); //correct this, integral doesnt account for bin width
    std::cout<<"integral="<<integral/power<<std::endl;
    totalIntXsec+=(integral/power);
    std::cout<<"totalIntXsec="<<totalIntXsec<<std::endl;
    
    
    theJetSpectraRapHist->SetMarkerStyle( kFullDotLarge ) ;
    theJetSpectraRapHist->SetMarkerColor( theRapOverlayMarkerColor[rapbin] );
    theJetSpectraRapHist->SetLineColor( kBlack );
    
    
    // MC Hist//hJetSpectraRap_1wJetID_bin5_genpt;1
    std::string inMCHistName="hJetSpectraRap";	
    if(doJetIDPlots) inMCHistName+="_1wJetID"; 
    else inMCHistName+="_0wJetID"; 
    inMCHistName+="_bin"+std::to_string(rapbin)+"_genpt";
    if(funcDebug)std::cout<<"inMCHistName="<<inMCHistName<<std::endl;	
    
    
    TH1F* theMCJetSpectraRapHist=NULL;
    theMCJetSpectraRapHist= (TH1F*)finMC->Get(inMCHistName.c_str()) ;
    theMCJetSpectraRapHist=(TH1F*)theMCJetSpectraRapHist->
      TH1::Rebin(nbins_pt_debug, (inMCHistName+"_rebin").c_str(), ptbins_debug);
    divideBinWidth(theMCJetSpectraRapHist);
    if(rapbin==7) theMCJetSpectraRapHist->Scale(1/3.); //3.2 - 4.7
    theMCJetSpectraRapHist->Scale( power * 1000 ); //i think the MC pthat weights are in units of nanobarns? multiply by 1000 
    //theMCJetSpectraRapHist->Scale( theJetSpectraRapHist->Integral()/theMCJetSpectraRapHist->Integral() );
    
    theMCJetSpectraRapHist->SetLineColor( theRapOverlayMarkerColor[rapbin] );
    
    
    std::stringstream stream1, stream2;
    stream1.precision(1); stream1 << std::fixed << absetabins[rapbin];
    stream2.precision(1); stream2 << std::fixed << absetabins[rapbin+1];
    
    if(rapbin==0){
      std::string rapHistLegDesc="(|y| < "+stream2.str()+") x 10^{"+std::to_string(bincounter)+"}";	
      theJetRapHistLeg->AddEntry(theJetSpectraRapHist,(rapHistLegDesc).c_str(),"lp");
      theJetRapHistLeg->AddEntry(theMCJetSpectraRapHist,"PY8 QCD MC GEN LEVEL","l");
    }
    else{
      std::string rapHistLegDesc="("+stream1.str()+" < |y| <"+stream2.str()+") x10^{"+std::to_string(bincounter)+"}";	
      theJetRapHistLeg->AddEntry(theJetSpectraRapHist,(rapHistLegDesc).c_str(),"lp");
      theJetRapHistLeg->AddEntry(theMCJetSpectraRapHist,"PY8 QCD MC GEN LEVEL","l");
    }
    
    if(rapbin==0){
      
      std::string h_Title   ="Jet Pt in |y| bins";
      if(doJetIDPlots) h_Title+=", w/ JetIDCut";
      theJetSpectraRapHist->SetTitle (    h_Title.c_str() );
      
      std::string h_XAx_Title="Jet p^{RECO}_{T} (GeV)", h_YAx_Title="d^{2}#sigma / d|y| dp_{T} (#mu b)";//crossSectionAxTitle;
      theJetSpectraRapHist->SetXTitle( h_XAx_Title.c_str() );
      theJetSpectraRapHist->SetYTitle( h_YAx_Title.c_str() );
      theJetSpectraRapHist->SetAxisRange(std::pow(10.,-9.),std::pow(10.,6.5),"Y");
      //theJetSpectraRapHist->SetAxisRange(60.,840.,"X");
      
      theJetSpectraRapHist->GetXaxis()->SetMoreLogLabels(true);
      theJetSpectraRapHist->GetXaxis()->SetNoExponent(true);
      theJetSpectraRapHist->Draw("E");     	  
      theMCJetSpectraRapHist->Draw("HIST E SAME");     	  
      theJetRapHistLeg->Draw(); 
    }
    else { 
      theJetSpectraRapHist->Draw("E SAME");
      theMCJetSpectraRapHist->Draw("HIST E SAME");
      }
    
  }//end rapbin hist loop
  
  std::cout<<"FINAL TOTAL INT JET XSEC="<<totalIntXsec<<std::endl;

  // STRINGS FOR HISTS
  TH1F* genjtptCut_h= (TH1F*)finMC->Get( "hGenJetPtCut" );
  std::string genjtptCut_str = std::to_string( (int) genjtptCut_h->GetMean() );
  std::cout<<"genjtptCut_str = "<<genjtptCut_str<<std::endl;
  
  std::string genJetCutString="p_{T}^{GEN} > "+genjtptCut_str+" GeV";
  std::cout<<"genJetCutString="<<genJetCutString<<std::endl;
  
  TH1F* jtptCut_h= (TH1F*)finMC->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}^{RECO} > "+jtptCut_str+" GeV";
  std::cout<<"jetCutString="<<jetCutString<<std::endl;

  float t1Loc1=0.62, t1Loc2=0.83, deltaLoc=.03; 	    
  TLatex *pp1=makeTLatex(t1Loc1,(t1Loc2),"PYTHIA8 ONLY, #sqrt{s} = 5.02 TeV"); pp1->Draw();
  TLatex *t1=makeTLatex(t1Loc1,(t1Loc2-deltaLoc),(fullJetType+"Jets").c_str());    t1->Draw();
  TLatex *t2=makeTLatex(t1Loc1,(t1Loc2-2*deltaLoc),(genJetCutString).c_str());    t2->Draw();
  TLatex *t3=makeTLatex(t1Loc1,(t1Loc2-3*deltaLoc),(jetCutString).c_str());    t3->Draw();
  if(doJetIDPlots){
    TLatex *t4=makeTLatex(t1Loc1,(t1Loc2-4*deltaLoc),"Jet ID Applied");      t4->Draw();	
  }
  
  temp_canvJetRapBins->Print( thePDFFileName.c_str() );
  temp_canvJetRapBins->Close();  
  return;
}

