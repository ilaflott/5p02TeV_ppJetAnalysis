// Raghav Kunnawalkam Elayavalli
// Nov 24th 2015
// CERN
// for questions or comments: raghav.k.e at CERN dot CH

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TEventList.h>
#include <TSystem.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"


// static const int nbins_pt = 39;
// static const double boundaries_pt[nbins_pt+1] = {
//   3, 4, 5, 7, 9, 12, 
//   15, 18, 21, 24, 28,
//   32, 37, 43, 49, 56,
//   64, 74, 84, 97, 114,
//   133, 153, 174, 196,
//   220, 245, 272, 300, 
//   330, 362, 395, 430,
//   468, 507, 548, 592,
//   638, 686, 1000 
// };

//these are the only radii we are interested for the RAA analysis: 2,3,4,5
//static const int no_radius = 3; 
//static const int list_radius[no_radius] = {2,3,4};
static const int trigValue = 4;
static const char trigName [trigValue][256] = {"HLT40","HLT60","HLT80","Combined"};
static const char * etaWidth = (char*)"20_eta_20";

// divide by bin width
void divideBinWidth(TH1 *h){
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();i++){
    Float_t val = h->GetBinContent(i);
    Float_t valErr = h->GetBinError(i);
    val/=h->GetBinWidth(i);
    valErr/=h->GetBinWidth(i);
    h->SetBinContent(i,val);
    h->SetBinError(i,valErr);
  }//binsX loop 
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

const int ptbins[] = {15, 30, 50, 80, 120, 170, 220, 300, 500};
const int nbins_pt = sizeof(ptbins)/sizeof(int) -1;


double trigComb(bool *trg, int *pscl, double pt){

  double weight=0;

  if(trg[3] && pt>=100) weight = pscl[3];
  if(trg[2] && pt>=80 && pt<100) weight = pscl[2];
  if(trg[1] && pt>=60 && pt<80) weight = pscl[1];
  if(trg[0] && pt>=40 && pt<60) weight = pscl[0];
  
  return weight;
  
}

using namespace std;

void RAA_read_data_pp(int startfile = 0,
		      int endfile = 1,
		      int radius = 4,
		      std::string jetType = "PF",
		      std::string kFoname="test_output.root"){
  
  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  bool printDebug = true;
  if(printDebug)cout<<"radius = "<<radius<<endl;
  
  TDatime date;

  std::string infile_Forest;

  infile_Forest = "pp_data_forests.txt";
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;
  
  if(printDebug)cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
  }

  const int N = 8; //6

  TChain * jetpp[N];

  string dir[N];
  dir[0] = "hltanalysis";
  dir[1] = "skimanalysis";
  dir[2] = Form("ak%d%sJetAnalyzer",radius, jetType.c_str());
  //dir[3] = "ppTrack";
  dir[3] = "hiEvtAnalyzer";
  dir[4] = "hltobject"; // hlt40
  dir[5] = "hltobject"; // hlt60
  dir[6] = "hltobject"; // hlt80
  dir[7] = "hltobject"; // hlt100

  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    //"trackTree"
    "HiTree",
    "HLT_AK4CaloJet40_Eta5p1_v",
    "HLT_AK4CaloJet60_Eta5p1_v",
    "HLT_AK4CaloJet80_Eta5p1_v",
    "HLT_AK4CaloJet100_Eta5p1_v"
  };

  for(int t = 0;t<N;t++){
    jetpp[t] = new TChain(string(dir[t]+"/"+trees[t]).data());
  }//tree loop ends

  for(int ifile = startfile; ifile<endfile; ++ifile){

    instr_Forest>>filename_Forest;

    jetpp[0]->Add(filename_Forest.c_str());
    jetpp[1]->Add(filename_Forest.c_str());
    jetpp[2]->Add(filename_Forest.c_str());
    jetpp[3]->Add(filename_Forest.c_str());
    jetpp[4]->Add(filename_Forest.c_str());
    jetpp[5]->Add(filename_Forest.c_str());
    jetpp[6]->Add(filename_Forest.c_str());
    jetpp[7]->Add(filename_Forest.c_str());
    
    if(printDebug)cout << "Tree loaded  " << string(dir[0]+"/"+trees[0]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpp[0]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[1]+"/"+trees[1]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpp[1]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[2]+"/"+trees[2]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpp[2]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[3]+"/"+trees[3]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpp[3]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[4]+"/"+trees[4]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpp[4]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[5]+"/"+trees[5]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpp[5]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[6]+"/"+trees[6]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpp[6]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[7]+"/"+trees[7]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpp[7]->GetEntries() << endl;
  }
  
  jetpp[2]->AddFriend(jetpp[0]);
  jetpp[2]->AddFriend(jetpp[1]);
  jetpp[2]->AddFriend(jetpp[3]);
  jetpp[2]->AddFriend(jetpp[4]);
  jetpp[2]->AddFriend(jetpp[5]);
  jetpp[2]->AddFriend(jetpp[6]);
  jetpp[2]->AddFriend(jetpp[7]);
  
  // Forest files 
  int nref_F;
  float pt_F[1000];
  float rawpt_F[1000];
  float eta_F[1000];
  float phi_F[1000];
  float chMax_F[1000];
  float trkMax_F[1000];
  float chSum_F[1000];
  int chN_F[1000];
  float phSum_F[1000];
  float neSum_F[1000];
  float trkSum_F[1000];
  float phMax_F[1000];
  float neMax_F[1000];
  float eMax_F[1000];
  float muMax_F[1000];
  float eSum_F[1000];
  float muSum_F[1000];
  float jtpu_F[1000];
  float chHardMax_F[1000];
  float trkHardMax_F[1000];
  float chHardSum_F[1000];
  float phHardSum_F[1000];
  float trkHardSum_F[1000];
  float phHardMax_F[1000];
  int jet40_F;
  int jet60_F;
  int jet80_F;
  int jet100_F;
  int jet40_p_F;
  int jet60_p_F;
  int jet80_p_F;
  int jet100_p_F;
  int jet40_l1seed_p_F;
  int jet60_l1seed_p_F;
  int jet80_l1seed_p_F;
  int jet100_l1seed_p_F;  
  ULong64_t evt_F;
  UInt_t run_F;
  UInt_t lumi_F;
  float vz_F;
  int pBeamScrapingFilter_F;
  int pHBHENoiseFilter_F;
  int pprimaryvertexFilter_F;

  jetpp[3]->SetBranchAddress("evt",&evt_F);
  jetpp[3]->SetBranchAddress("run",&run_F);
  jetpp[3]->SetBranchAddress("lumi",&lumi_F);
  jetpp[3]->SetBranchAddress("vz",&vz_F);

  jetpp[1]->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_F);
  jetpp[1]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilter_F);
  jetpp[1]->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryvertexFilter_F);

  jetpp[2]->SetBranchAddress("nref",&nref_F);
  jetpp[2]->SetBranchAddress("jtpt",&pt_F);
  jetpp[2]->SetBranchAddress("jteta",&eta_F);
  jetpp[2]->SetBranchAddress("jtphi",&phi_F);
  jetpp[2]->SetBranchAddress("rawpt",&rawpt_F);
  jetpp[2]->SetBranchAddress("jtpu",&jtpu_F);
  jetpp[2]->SetBranchAddress("chargedMax",&chMax_F);
  jetpp[2]->SetBranchAddress("chargedSum",&chSum_F);
  jetpp[2]->SetBranchAddress("chargedN",&chN_F);
  //jetpp[2]->SetBranchAddress("chargedHardMax",&chMax_F);
  jetpp[2]->SetBranchAddress("chargedHardSum",&chSum_F);
  jetpp[2]->SetBranchAddress("trackMax",&trkMax_F);
  jetpp[2]->SetBranchAddress("trackSum",&trkSum_F);
  //jetpp[2]->SetBranchAddress("trackHardMax",&trkMax_F);
  jetpp[2]->SetBranchAddress("trackHardSum",&trkSum_F);
  jetpp[2]->SetBranchAddress("photonMax",&phMax_F);
  jetpp[2]->SetBranchAddress("photonSum",&phSum_F);
  //jetpp[2]->SetBranchAddress("photonHardMax",&phMax_F);
  jetpp[2]->SetBranchAddress("photonHardSum",&phSum_F);
  jetpp[2]->SetBranchAddress("neutralMax",&neMax_F);
  jetpp[2]->SetBranchAddress("neutralSum",&neSum_F);

  jetpp[2]->SetBranchAddress("eSum",&eSum_F);
  jetpp[2]->SetBranchAddress("eMax",&eMax_F);
  jetpp[2]->SetBranchAddress("muSum",&muSum_F);
  jetpp[2]->SetBranchAddress("muMax",&muMax_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet40_Eta5p1_v1",&jet40_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet40_Eta5p1_v1_Prescl",&jet40_p_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet60_Eta5p1_v1",&jet60_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet60_Eta5p1_v1_Prescl",&jet60_p_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1",&jet80_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1_Prescl",&jet80_p_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet100_Eta5p1_v1",&jet100_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet100_Eta5p1_v1_Prescl",&jet100_p_F);
  jetpp[0]->SetBranchAddress("L1_SingleJet28_BptxAND_Prescl",&jet40_l1seed_p_F);
  jetpp[0]->SetBranchAddress("L1_SingleJet40_BptxAND_Prescl",&jet60_l1seed_p_F);
  jetpp[0]->SetBranchAddress("L1_SingleJet48_BptxAND_Prescl",&jet80_l1seed_p_F);
  jetpp[0]->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&jet100_l1seed_p_F);
    
  std::vector<float> *trgObjpt_40 = 0;
  jetpp[4]->SetBranchAddress("pt",&trgObjpt_40);
  std::vector<float> *trgObjpt_60 = 0;
  jetpp[5]->SetBranchAddress("pt",&trgObjpt_60);
  std::vector<float> *trgObjpt_80 = 0;
  jetpp[6]->SetBranchAddress("pt",&trgObjpt_80);
  std::vector<float> *trgObjpt_100 = 0;
  jetpp[7]->SetBranchAddress("pt",&trgObjpt_100);

  
  // Declare the output File and the necessary histograms after that:
  // std::string outdir="";
  // std::string outfile=outdir+kFoname;
  TFile *fout = new TFile(kFoname.c_str(),"RECREATE");
  fout->cd();

  // Add the Jet ID plots:
  // list of variables:
  std::string var[21] = {"jtpt" ,"rawpt", "jteta", "jtphi", "trkMax", "trkSum", "trkHardSum", "chMax", "chSum", "chHardSum","phMax", "phSum", "phHardSum", "neMax", "neSum", "eMax", "eSum", "muMax", "muSum","Aj","xj"};
  TH1F * hJetQA[2][21];

  for(int k = 0; k<2; ++k){
    for(int j = 0; j<21; ++j){
      if(j==2) hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),Form(";%s;",var[j].c_str()),100, -5, +5);
      else if(j==3) hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),Form(";%s;",var[j].c_str()),100, -4, +4);
      else if(j<=1)hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),Form(";%s;",var[j].c_str()),500, 0, 500);
      else if(j>=4)hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),Form(";%s;",var[j].c_str()),200, 0, 2);
    }
  }

  TH1F * hVz = new TH1F("hVz","",200, -20, 20);

  TH1F *hpp_TrgObj100[2];
  TH1F *hpp_TrgObj80[2];
  TH1F *hpp_TrgObj60[2];
  TH1F *hpp_TrgObj40[2];
  TH1F *hpp_TrgObjComb[2];
  TH1F *hpp_CombJetpT[2];

  // TH1F *hpp_JEC_TrgObj80;
  // TH1F *hpp_JEC_TrgObj60;
  // TH1F *hpp_JEC_TrgObj40;
  // TH1F *hpp_JEC_TrgObjComb;

  // TH1F *hpp_JEC_minus_TrgObj80;
  // TH1F *hpp_JEC_minus_TrgObj60;
  // TH1F *hpp_JEC_minus_TrgObj40;
  // TH1F *hpp_JEC_minus_TrgObjComb;

  // TH1F *hpp_JEC_gaus_TrgObj80;
  // TH1F *hpp_JEC_gaus_TrgObj60;
  // TH1F *hpp_JEC_gaus_TrgObj40;
  // TH1F *hpp_JEC_gaus_TrgObjComb;
  
  // TH1F *hpp_Smear_TrgObj80;
  // TH1F *hpp_Smear_TrgObj60;
  // TH1F *hpp_Smear_TrgObj40;
  // TH1F *hpp_Smear_TrgObjComb;
  
  // TH1F *hpp_anaBin_TrgObj80;
  // TH1F *hpp_anaBin_TrgObj60;
  // TH1F *hpp_anaBin_TrgObj40;
  // TH1F *hpp_anaBin_TrgObjComb;

  hpp_TrgObj100[0] = new TH1F(Form("hpp_HLT100_noJetID_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 100 R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_TrgObj80[0] = new TH1F(Form("hpp_HLT80_noJetID_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_TrgObj60[0] = new TH1F(Form("hpp_HLT60_noJetID_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_TrgObj40[0] = new TH1F(Form("hpp_HLT40_noJetID_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_TrgObjComb[0] = new TH1F(Form("hpp_HLTComb_noJetID_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_CombJetpT[0] = new TH1F(Form("hpp_TrgCombTest_noJetID_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra KurtMethod R%d %s ",radius,etaWidth),2000, 0, 2000);

  hpp_TrgObj100[1] = new TH1F(Form("hpp_HLT100_JetID_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 100 R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_TrgObj80[1] = new TH1F(Form("hpp_HLT80_JetID_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_TrgObj60[1] = new TH1F(Form("hpp_HLT60_JetID_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_TrgObj40[1] = new TH1F(Form("hpp_HLT40_JetID_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_TrgObjComb[1] = new TH1F(Form("hpp_HLTComb_JetID_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_CombJetpT[1] = new TH1F(Form("hpp_TrgCombTest_JetID_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra KurtMethod R%d %s ",radius,etaWidth),2000, 0, 2000);

  
  // hpp_JEC_TrgObj80 = new TH1F(Form("hpp_JEC_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_TrgObj60 = new TH1F(Form("hpp_JEC_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_TrgObj40 = new TH1F(Form("hpp_JEC_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_TrgObjComb = new TH1F(Form("hpp_JEC_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),1000, 0, 1000);

  // hpp_JEC_minus_TrgObj80 = new TH1F(Form("hpp_JEC_minus_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_minus_TrgObj60 = new TH1F(Form("hpp_JEC_minus_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_minus_TrgObj40 = new TH1F(Form("hpp_JEC_minus_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_minus_TrgObjComb = new TH1F(Form("hpp_JEC_minus_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),1000, 0, 1000);

  // hpp_JEC_gaus_TrgObj80 = new TH1F(Form("hpp_JEC_gaus_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_gaus_TrgObj60 = new TH1F(Form("hpp_JEC_gaus_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_gaus_TrgObj40 = new TH1F(Form("hpp_JEC_gaus_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_gaus_TrgObjComb = new TH1F(Form("hpp_JEC_gaus_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),1000, 0, 1000);

  
  // hpp_Smear_TrgObj80 = new TH1F(Form("hpp_Smear_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_Smear_TrgObj60 = new TH1F(Form("hpp_Smear_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_Smear_TrgObj40 = new TH1F(Form("hpp_Smear_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_Smear_TrgObjComb = new TH1F(Form("hpp_Smear_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),1000, 0, 1000);

  // hpp_anaBin_TrgObj80 = new TH1F(Form("hpp_anaBin_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  // hpp_anaBin_TrgObj60 = new TH1F(Form("hpp_anaBin_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  // hpp_anaBin_TrgObj40 = new TH1F(Form("hpp_anaBin_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  // hpp_anaBin_TrgObjComb = new TH1F(Form("hpp_anaBin_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);

  // now start the event loop for each file. 
  
  Long64_t nentries = jetpp[0]->GetEntries();
  //Long64_t nGoodEvt = 0;
  //if(printDebug) nentries = 100;
  
  cout<<"Running through all the "<<nentries<<" events now"<<endl;

  TRandom3 rnd; 
  //rnd.SetSeed(endfile);

  // variables for event summary
  Long64_t NEvents = 0;
  Long64_t NEvents_100 = 0;
  Long64_t NEvents_80 = 0;
  Long64_t NEvents_60 = 0;
  Long64_t NEvents_40 = 0;  
  
  for(int nEvt = 0; nEvt < nentries; ++ nEvt) {

    if(nEvt%10000 == 0)cout<<nEvt<<"/"<<nentries<<endl;
    //if(printDebug)cout<<"nEvt = "<<nEvt<<endl;

    for(int i = 0; i<N;++i)
      jetpp[i]->GetEntry(nEvt);
    
    if(pHBHENoiseFilter_F==0) continue; 
    if(pBeamScrapingFilter_F==0) continue; 
    if(pprimaryvertexFilter_F==0) continue; 
    if(fabs(vz_F)>15) continue;

    hVz->Fill(vz_F);

    bool is40 = false;
    bool is60 = false;
    bool is80 = false;
    bool is100 = false;

    // vector <float> arry_TrgObj;
    // if(trgObjpt_40->size()!=0)
    //   for(unsigned i = 0; i<trgObjpt_40->size(); ++i) arry_TrgObj.push_back(trgObjpt_40->at(i));
    // if(trgObjpt_60->size()!=0) 
    //   for(unsigned i = 0; i<trgObjpt_60->size(); ++i) arry_TrgObj.push_back(trgObjpt_60->at(i));
    // if(trgObjpt_80->size()!=0) 
    //   for(unsigned i = 0; i<trgObjpt_80->size(); ++i) arry_TrgObj.push_back(trgObjpt_80->at(i));
    // if(trgObjpt_100->size()!=0)
    //   for(unsigned i = 0; i<trgObjpt_100->size(); ++i) arry_TrgObj.push_back(trgObjpt_100->at(i));

    // // find the max in the array
    // float TrgObjMax = 0.0;
    // for(unsigned i = 0; i<arry_TrgObj.size();++i){
    //   if(arry_TrgObj[i] > TrgObjMax){
    // 	TrgObjMax = arry_TrgObj[i];
    //   }
    // }

    // if(jet40_F && TrgObjMax>=40 && TrgObjMax<60) is40 = true;
    // if(jet60_F && TrgObjMax>=60 && TrgObjMax<80) is60 = true;
    // if(jet80_F && TrgObjMax>=80 && TrgObjMax<100) is80 = true;
    // if(jet100_F && TrgObjMax>=100) is100 = true;
    // int weight_jet40  = jet40_p_F * jet40_l1seed_p_F; 
    // int weight_jet60  = jet60_p_F * jet60_l1seed_p_F; 
    // int weight_jet80  = jet80_p_F * jet80_l1seed_p_F; 
    // int weight_jet100 = jet100_p_F * jet100_l1seed_p_F;

    // double weight_eS = 1.0;
    // if(is40) weight_eS = weight_jet40;
    // if(is60) weight_eS = weight_jet60;
    // if(is80) weight_eS = weight_jet80;
    // if(is100) weight_eS = weight_jet100;

    // arry_TrgObj.clear();

    
    bool trgDec[4] = {(bool)jet40_F, (bool)jet60_F, (bool)jet80_F, (bool)jet100_F};
    int treePrescl[4] = {(jet40_p_F*jet40_l1seed_p_F), (jet60_p_F*jet60_l1seed_p_F), (jet80_p_F*jet80_l1seed_p_F), (jet100_p_F*jet100_l1seed_p_F)};
    
    int maxtrg= -1;
    for(int ii=4; ii>=0; ii--){
      if(trgDec[ii]==1){
	maxtrg=ii;
	break;
      }
    }
    
    double triggerPt=0.;
    if(jet40_F){
      for(unsigned int itt=0; itt<trgObjpt_40->size(); itt++){
	if(trgObjpt_40->at(itt)>triggerPt) triggerPt = trgObjpt_40->at(itt);
      }
    }
    if(jet60_F){
      for(unsigned int itt=0; itt<trgObjpt_60->size(); itt++){
	if(trgObjpt_60->at(itt)>triggerPt) triggerPt = trgObjpt_60->at(itt);
      }
    }
    if(jet80_F){
      for(unsigned int itt=0; itt<trgObjpt_80->size(); itt++){
	if(trgObjpt_80->at(itt)>triggerPt) triggerPt = trgObjpt_80->at(itt);
      }
    }
    if(jet100_F){
      for(unsigned int itt=0; itt<trgObjpt_100->size(); itt++){
	if(trgObjpt_100->at(itt)>triggerPt) triggerPt = trgObjpt_100->at(itt);
      }
    }
    
    double weight_eS = trigComb(trgDec, treePrescl, triggerPt);
    
    if(trgDec[3] && triggerPt>=100) is100 = true;
    if(trgDec[2] && triggerPt>=80 && triggerPt<100) is80 = true;
    if(trgDec[1] && triggerPt>=60 && triggerPt<80) is60 = true;
    if(trgDec[0] && triggerPt>=40 && triggerPt<60) is40 = true;

    NEvents++;
    if(is100) NEvents_100++;
    if(is80) NEvents_80++;
    if(is60) NEvents_60++;
    if(is40) NEvents_40++;
    

    for(int jet = 0; jet<nref_F; ++jet){

      if(fabs(eta_F[jet]) > 2) continue;
      
      float recpt = pt_F[jet];
      float rawpt = rawpt_F[jet];
      if(recpt<=15) continue;
      
      // float JEC_upShift   = (float)recpt;
      // float JEC_downShift = (float)recpt;
      // float JEC_gaussmear = (float)recpt;

      
      // if(jet40_F == 1 && jet60_F==0 && jet80_F==0) {
      // 	hpp_TrgObj40->Fill(recpt, jet40_p_F);
      // // 	//hpp_JEC_TrgObj40->Fill(recpt * (1 + rnd.Gaus(0, 0.01)), jet40_p_F);
      // // 	hpp_JEC_TrgObj40->Fill(JEC_upShift, jet40_p_F);
      // // 	hpp_JEC_minus_TrgObj40->Fill(JEC_downShift, jet40_p_F);
      // // 	hpp_JEC_gaus_TrgObj40->Fill(JEC_gaussmear, jet40_p_F);
      // // 	hpp_Smear_TrgObj40->Fill(recpt+ rnd.Gaus(0,0.01), jet40_p_F);
      // // 	hpp_anaBin_TrgObj40->Fill(recpt, jet40_p_F);
      // }
    
      // if(jet60_F == 1 && jet80_F==0) {
      // 	hpp_TrgObj60->Fill(recpt);
      // // 	//hpp_JEC_TrgObj60->Fill(recpt * (1 + rnd.Gaus(0, 0.01)));
      // // 	hpp_JEC_TrgObj60->Fill(JEC_upShift);
      // // 	hpp_JEC_minus_TrgObj60->Fill(JEC_downShift);
      // // 	hpp_JEC_gaus_TrgObj60->Fill(JEC_gaussmear);
      // // 	hpp_Smear_TrgObj60->Fill(recpt+rnd.Gaus(0,0.01));
      // // 	hpp_anaBin_TrgObj60->Fill(recpt);
      // }

      if(1){
	hJetQA[0][0]->Fill(recpt, weight_eS);
	hJetQA[0][1]->Fill(rawpt_F[jet], weight_eS);
	hJetQA[0][2]->Fill(eta_F[jet], weight_eS);
	hJetQA[0][3]->Fill(phi_F[jet], weight_eS);
	hJetQA[0][4]->Fill(trkMax_F[jet]/recpt, weight_eS);
	hJetQA[0][5]->Fill(trkSum_F[jet]/recpt, weight_eS);
	hJetQA[0][6]->Fill(trkHardSum_F[jet]/recpt, weight_eS);
	hJetQA[0][7]->Fill(chMax_F[jet]/recpt, weight_eS);
	hJetQA[0][8]->Fill(chSum_F[jet]/recpt, weight_eS);
	hJetQA[0][9]->Fill(chHardSum_F[jet]/recpt, weight_eS);
	hJetQA[0][10]->Fill(phMax_F[jet]/recpt, weight_eS);
	hJetQA[0][11]->Fill(phSum_F[jet]/recpt, weight_eS);
	hJetQA[0][12]->Fill(phHardSum_F[jet]/recpt, weight_eS);
	hJetQA[0][13]->Fill(neMax_F[jet]/recpt, weight_eS);
	hJetQA[0][14]->Fill(neSum_F[jet]/recpt, weight_eS);
	hJetQA[0][15]->Fill(eMax_F[jet]/recpt, weight_eS);
	hJetQA[0][16]->Fill(eSum_F[jet]/recpt, weight_eS);
	hJetQA[0][17]->Fill(muMax_F[jet]/recpt, weight_eS);
	hJetQA[0][18]->Fill(muSum_F[jet]/recpt, weight_eS);

	hpp_CombJetpT[0]->Fill(recpt, weight_eS);
	
	if(is40) hpp_TrgObj40[0]->Fill(recpt, treePrescl[0]);
	if(is60) hpp_TrgObj60[0]->Fill(recpt, treePrescl[1]);
	if(is80) hpp_TrgObj80[0]->Fill(recpt, treePrescl[2]);
	if(is100) hpp_TrgObj100[0]->Fill(recpt, treePrescl[3]);
	// apply JetID
	// charged hadron fraction > 0
	// charged hadron multiplicity > 0
	// charged EM fraction < 0.99
	
	if(jet == 0){
	  hJetQA[0][19]->Fill((float)(pt_F[0]-pt_F[1])/(pt_F[0]+pt_F[1]), weight_eS);
	  hJetQA[0][20]->Fill((float)(pt_F[1]/pt_F[0]), weight_eS);
	}
	
	if( chSum_F[jet]/recpt>0 && neSum_F[jet]/recpt<0.99 && chN_F[jet]>0 && phSum_F[jet]/recpt<0.99){
	  hJetQA[1][0]->Fill(recpt, weight_eS);
	  hJetQA[1][1]->Fill(rawpt_F[jet], weight_eS);
	  hJetQA[1][2]->Fill(eta_F[jet], weight_eS);
	  hJetQA[1][3]->Fill(phi_F[jet], weight_eS);
	  hJetQA[1][4]->Fill(trkMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][5]->Fill(trkSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][6]->Fill(trkHardSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][7]->Fill(chMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][8]->Fill(chSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][9]->Fill(chHardSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][10]->Fill(phMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][11]->Fill(phSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][12]->Fill(phHardSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][13]->Fill(neMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][14]->Fill(neSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][15]->Fill(eMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][16]->Fill(eSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][17]->Fill(muMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][18]->Fill(muSum_F[jet]/recpt, weight_eS);

	  hpp_CombJetpT[1]->Fill(recpt, weight_eS);
	  
	  if(is40) hpp_TrgObj40[1]->Fill(recpt, treePrescl[0]);
	  if(is60) hpp_TrgObj60[1]->Fill(recpt, treePrescl[1]);
	  if(is80) hpp_TrgObj80[1]->Fill(recpt, treePrescl[2]);
	  if(is100) hpp_TrgObj100[1]->Fill(recpt, treePrescl[3]);
	  
	  if(jet == 0){
	    hJetQA[1][19]->Fill((float)(pt_F[0]-pt_F[1])/(pt_F[0]+pt_F[1]), weight_eS);
	    hJetQA[1][20]->Fill((float)(pt_F[1]/pt_F[0]), weight_eS);
	  }  
	}
		
      }
      // if(jet80_F == 1) {      

      // 	hpp_TrgObj80->Fill(recpt);
      // // 	hpp_JEC_TrgObj80->Fill(JEC_upShift);
      // // 	hpp_JEC_minus_TrgObj80->Fill(JEC_downShift);
      // // 	hpp_JEC_gaus_TrgObj80->Fill(JEC_gaussmear);
      // // 	hpp_Smear_TrgObj80->Fill(recpt+rnd.Gaus(0,0.01));
      // // 	hpp_anaBin_TrgObj80->Fill(recpt);
      // }
      
    }// jet loop


	
    //if(printDebug)cout<<endl;

  }// event loop


  hpp_TrgObjComb[0]->Add(hpp_TrgObj100[0]);
  hpp_TrgObjComb[0]->Add(hpp_TrgObj80[0]);
  hpp_TrgObjComb[0]->Add(hpp_TrgObj60[0]);
  hpp_TrgObjComb[0]->Add(hpp_TrgObj40[0]);

  hpp_TrgObjComb[1]->Add(hpp_TrgObj100[1]);
  hpp_TrgObjComb[1]->Add(hpp_TrgObj80[1]);
  hpp_TrgObjComb[1]->Add(hpp_TrgObj60[1]);
  hpp_TrgObjComb[1]->Add(hpp_TrgObj40[1]);
  
  // divideBinWidth(hpp_TrgObjComb);
  // divideBinWidth(hpp_TrgObj80);
  // divideBinWidth(hpp_TrgObj60);
  // divideBinWidth(hpp_TrgObj40);
  
  // hpp_JEC_TrgObjComb->Add(hpp_JEC_TrgObj80);
  // hpp_JEC_TrgObjComb->Add(hpp_JEC_TrgObj60);
  // hpp_JEC_TrgObjComb->Add(hpp_JEC_TrgObj40);

  // divideBinWidth(hpp_JEC_TrgObjComb);
  // divideBinWidth(hpp_JEC_TrgObj80);
  // divideBinWidth(hpp_JEC_TrgObj60);
  // divideBinWidth(hpp_JEC_TrgObj40);

  // hpp_JEC_minus_TrgObjComb->Add(hpp_JEC_minus_TrgObj80);
  // hpp_JEC_minus_TrgObjComb->Add(hpp_JEC_minus_TrgObj60);
  // hpp_JEC_minus_TrgObjComb->Add(hpp_JEC_minus_TrgObj40);

  // divideBinWidth(hpp_JEC_minus_TrgObjComb);
  // divideBinWidth(hpp_JEC_minus_TrgObj80);
  // divideBinWidth(hpp_JEC_minus_TrgObj60);
  // divideBinWidth(hpp_JEC_minus_TrgObj40);

  // hpp_JEC_gaus_TrgObjComb->Add(hpp_JEC_gaus_TrgObj80);
  // hpp_JEC_gaus_TrgObjComb->Add(hpp_JEC_gaus_TrgObj60);
  // hpp_JEC_gaus_TrgObjComb->Add(hpp_JEC_gaus_TrgObj40);

  // divideBinWidth(hpp_JEC_gaus_TrgObjComb);
  // divideBinWidth(hpp_JEC_gaus_TrgObj80);
  // divideBinWidth(hpp_JEC_gaus_TrgObj60);
  // divideBinWidth(hpp_JEC_gaus_TrgObj40);
  
  // hpp_Smear_TrgObjComb->Add(hpp_Smear_TrgObj80);
  // hpp_Smear_TrgObjComb->Add(hpp_Smear_TrgObj60);
  // hpp_Smear_TrgObjComb->Add(hpp_Smear_TrgObj40);

  // divideBinWidth(hpp_Smear_TrgObjComb);
  // divideBinWidth(hpp_Smear_TrgObj80);
  // divideBinWidth(hpp_Smear_TrgObj60);
  // divideBinWidth(hpp_Smear_TrgObj40);

  // hpp_anaBin_TrgObjComb->Add(hpp_anaBin_TrgObj80);
  // hpp_anaBin_TrgObjComb->Add(hpp_anaBin_TrgObj60);
  // hpp_anaBin_TrgObjComb->Add(hpp_anaBin_TrgObj40);

  // divideBinWidth(hpp_anaBin_TrgObjComb);
  // divideBinWidth(hpp_anaBin_TrgObj80);
  // divideBinWidth(hpp_anaBin_TrgObj60);
  // divideBinWidth(hpp_anaBin_TrgObj40);

  
  hpp_TrgObj40[0]->Print("base");
  hpp_TrgObj60[0]->Print("base");
  hpp_TrgObj80[0]->Print("base");
  hpp_TrgObj100[0]->Print("base");
  hpp_TrgObjComb[0]->Print("base");
  hpp_CombJetpT[0]->Print("base");
  
  hpp_TrgObj40[1]->Print("base");
  hpp_TrgObj60[1]->Print("base");
  hpp_TrgObj80[1]->Print("base");
  hpp_TrgObj100[1]->Print("base");
  hpp_TrgObjComb[1]->Print("base");
  hpp_CombJetpT[1]->Print("base");
  

  cout<<"Event Summary"<<endl;
  cout<<"Total Number of Events read   = "<<NEvents<<endl;
  cout<<"Total Number of HLT100 events = "<<NEvents_100<<endl;
  cout<<"Total Number of HLT80 events = "<<NEvents_80<<endl;
  cout<<"Total Number of HLT60 events = "<<NEvents_60<<endl;
  cout<<"Total Number of HLT40 events = "<<NEvents_40<<endl;
  
  fout->Write();

  trgObjpt_40->clear();
  trgObjpt_60->clear();
  trgObjpt_80->clear();
  trgObjpt_100->clear();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}//macro end
