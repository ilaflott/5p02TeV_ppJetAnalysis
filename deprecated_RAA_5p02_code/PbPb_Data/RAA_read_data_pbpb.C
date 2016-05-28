// Raghav Kunnawalkam Elayavalli
// Nov 24 2015
// Rutgers, @CERN
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


static const int nbins_pt = 39;
static const double boundaries_pt[nbins_pt+1] = {
  3, 4, 5, 7, 9, 12, 
  15, 18, 21, 24, 28,
  32, 37, 43, 49, 56,
  64, 74, 84, 97, 114,
  133, 153, 174, 196,
  220, 245, 272, 300, 
  330, 362, 395, 430,
  468, 507, 548, 592,
  638, 686, 1000 
};

//these are the only radii we are interested for the RAA analysis: 2,3,4,5
//static const int no_radius = 3; 
//static const int list_radius[no_radius] = {2,3,4};

// static const int nbins_cent = 6;
// static const Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 2.5 to get your actual centrality % (old 2011 data) 

// 0-10, 10-30, 30-50, 50-100
static const int nbins_cent = 4;
static const Double_t boundaries_cent[nbins_cent+1] = {0, 20, 60, 100, 200};// multiply by 2.5 to get your actual centrality % (old 2011 data) 

//now we have to multiply by 5, since centrality goes from 0-200. 
// static const Double_t ncoll[nbins_cent] = { 1660, 1310, 745, 251, 62.8, 10.8 };
// static const int trigValue = 4;
// static const char trigName [trigValue][256] = {"HLT55","HLT65","HLT80","Combined"};
// static const Float_t effecPrescl = 2.047507;
// static const char * etaWidth = (char*)"20_eta_20";

int findBin(int bin)
{
  int ibin=-1;
  //! centrality is defined as 0.5% bins of cross section
  //! in 0-200 bins               
  // if(bin<10)ibin=0; //! 0-5%
  // else if(bin>=10  && bin<20 )ibin=1; //! 5-10%
  // else if(bin>=20  && bin<60 )ibin=2;  //! 10-30%
  // else if(bin>=60  && bin<100)ibin=3;  //! 30-50%
  // else if(bin>=100 && bin<140)ibin=4;  //! 50-70%
  // else if(bin>=140 && bin<180)ibin=5;  //! 70-90%
  // else if(bin>=180 && bin<200)ibin=6;  //! 90-100%

  if(bin < 20) ibin = 0;
  else if(bin >= 20 && bin < 60) ibin = 1;
  else if(bin >= 60 && bin < 100) ibin = 2;
  else if(bin >= 100 && bin < 200) ibin = 3;
  return ibin;
}

static const int ptBoundary[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 500};
static const int ptSelection = sizeof(ptBoundary)/sizeof(int) -1;

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

struct Jet{
  int id;
  float pt;
};
bool compare_pt(Jet jet1, Jet jet2);
bool compare_pt(Jet jet1, Jet jet2){
  return jet1.pt > jet2.pt;
}

static const Double_t ptbins_jec[] = {40.0, 46, 54.5, 60, 69, 79, 90.5, 105.5, 123.5, 143, 163.5, 185, 208, 232, 258, 286, 331};
static const int nbins_jec = sizeof(ptbins_jec)/sizeof(Double_t) -1;
static const Double_t jec_shift[] = {0.0205, 0.0179, 0.0160, 0.0146, 0.0133, 0.0122, 0.0113, 0.0104, 0.0097, 0.0091, 0.0086, 0.0082, 0.0082, 0.0077, 0.0075, 0.0074, 0.0073};

using namespace std;

void RAA_read_data_pbpb(int startfile = 0,
			int endfile = 1,
			int radius = 4,
			std::string algo = "Vs",
			std::string jet_type = "PF",
			std::string kFoname="test_output.root"){
  
  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  bool printDebug = false;
  if(printDebug)cout<<"radius = "<<radius<<endl;
  
  TDatime date;

  std::string infile_Forest;

  infile_Forest = "jetRAA_PbPb_data_forests.txt";
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;
  
  if(printDebug)cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
  }

  const int N = 4; //6

  TChain * jetpbpb[N];

  // Float_t ScaleFactor_55_NeqScale[nbins_cent] = {53.6 * 1e6, 53.6 * 1e6, 214.4 * 1e6, 214.4 * 1e6, 214.4 * 1e6, 214.4 * 1e6}; 
  // Float_t ScaleFactor_65_NeqScale[nbins_cent] = {53.6 * 1e6, 53.6 * 1e6, 214.4 * 1e6, 214.4 * 1e6, 214.4 * 1e6, 214.4 * 1e6}; 
  // Float_t ScaleFactor_80_NeqScale[nbins_cent] = {53.6 * 1e6, 53.6 * 1e6, 214.4 * 1e6, 214.4 * 1e6, 214.4 * 1e6, 214.4 * 1e6};   
  
  string dir[N];
  dir[0] = "hltanalysis";
  dir[1] = "skimanalysis";
  dir[2] = Form("ak%s%d%sJetAnalyzer",algo.c_str(),radius,jet_type.c_str());
  dir[3] = "hiEvtAnalyzer";
  // dir[4] = "hltobject";

  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    "HiTree"
    // , "jetObjTree"
  };

  for(int t = 0;t<N;t++){
    jetpbpb[t] = new TChain(string(dir[t]+"/"+trees[t]).data());
  }//tree loop ends

  for(int ifile = startfile; ifile<endfile; ++ifile){

    instr_Forest>>filename_Forest;

    if(printDebug)cout<<"HiForest filename = "<<filename_Forest.c_str()<<endl;

    jetpbpb[0]->Add(filename_Forest.c_str());
    jetpbpb[1]->Add(filename_Forest.c_str());
    jetpbpb[2]->Add(filename_Forest.c_str());
    jetpbpb[3]->Add(filename_Forest.c_str());
    
    if(printDebug)cout << "Tree loaded  " << string(dir[0]+"/"+trees[0]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[0]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[1]+"/"+trees[1]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[1]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[2]+"/"+trees[2]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[2]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[3]+"/"+trees[3]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[3]->GetEntries() << endl;
  }
  
  jetpbpb[2]->AddFriend(jetpbpb[0]);
  jetpbpb[2]->AddFriend(jetpbpb[1]);
  jetpbpb[2]->AddFriend(jetpbpb[3]);
  
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
  int L1_sj36_F;
  int L1_sj52_F;
  int L1_sj36_p_F;
  int L1_sj52_p_F;
  int jet40_p_F;
  int jet60_p_F;
  int jet80_p_F;
  float vz_F;
  int evt_F;
  int run_F;
  int lumi_F;
  int hiNpix_F;
  int hiNpixelTracks_F;
  int hiBin_F;
  float hiHF_F;
  float hiZDC_F;
  float hiZDCplus_F;
  float hiZDCminus_F;
  int hiNtracks_F;
  int hiNtracksPtCut_F;
  int hiNtracksEtaCut_F;
  int hiNtracksEtaPtCut_F;
  int pcollisionEventSelection_F;
  int pHBHENoiseFilter_F;
  
  jetpbpb[3]->SetBranchAddress("evt",&evt_F);
  jetpbpb[3]->SetBranchAddress("run",&run_F);
  jetpbpb[3]->SetBranchAddress("lumi",&lumi_F);
  jetpbpb[3]->SetBranchAddress("hiBin",&hiBin_F);
  jetpbpb[3]->SetBranchAddress("hiHF", &hiHF_F);
  jetpbpb[3]->SetBranchAddress("hiNpix",&hiNpix_F);
  jetpbpb[3]->SetBranchAddress("hiNpixelTracks",&hiNpixelTracks_F);
  jetpbpb[3]->SetBranchAddress("hiNtracks",&hiNtracks_F);
  jetpbpb[3]->SetBranchAddress("hiNtracksPtCut",&hiNtracksPtCut_F);
  jetpbpb[3]->SetBranchAddress("hiNtracksEtaCut",&hiNtracksEtaCut_F);
  jetpbpb[3]->SetBranchAddress("hiNtracksEtaPtCut",&hiNtracksEtaPtCut_F);
  jetpbpb[3]->SetBranchAddress("hiZDC", &hiZDC_F);
  jetpbpb[3]->SetBranchAddress("hiZDCplus", &hiZDCplus_F);
  jetpbpb[3]->SetBranchAddress("hiZDCminus", &hiZDCminus_F);
  jetpbpb[3]->SetBranchAddress("vz",&vz_F);
  jetpbpb[1]->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection_F);
  jetpbpb[1]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_F);
  //jetpbpb[0]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_F);
  //jetpbpb[0]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_F);
  jetpbpb[2]->SetBranchAddress("nref",&nref_F);
  jetpbpb[2]->SetBranchAddress("jtpt",pt_F);
  jetpbpb[2]->SetBranchAddress("jteta",eta_F);
  jetpbpb[2]->SetBranchAddress("jtphi",phi_F);
  jetpbpb[2]->SetBranchAddress("rawpt",rawpt_F);
  jetpbpb[2]->SetBranchAddress("jtpu",jtpu_F);
  
  jetpbpb[2]->SetBranchAddress("chargedMax",&chMax_F);
  jetpbpb[2]->SetBranchAddress("chargedSum",&chSum_F);
  jetpbpb[2]->SetBranchAddress("chargedN",&chN_F);
  //jetpbpb[2]->SetBranchAddress("chargedHardMax",&chMax_F);
  jetpbpb[2]->SetBranchAddress("chargedHardSum",&chSum_F);
  jetpbpb[2]->SetBranchAddress("trackMax",&trkMax_F);
  jetpbpb[2]->SetBranchAddress("trackSum",&trkSum_F);
  //jetpbpb[2]->SetBranchAddress("trackHardMax",&trkMax_F);
  jetpbpb[2]->SetBranchAddress("trackHardSum",&trkSum_F);
  jetpbpb[2]->SetBranchAddress("photonMax",&phMax_F);
  jetpbpb[2]->SetBranchAddress("photonSum",&phSum_F);
  //jetpbpb[2]->SetBranchAddress("photonHardMax",&phMax_F);
  jetpbpb[2]->SetBranchAddress("photonHardSum",&phSum_F);
  jetpbpb[2]->SetBranchAddress("neutralMax",&neMax_F);
  jetpbpb[2]->SetBranchAddress("neutralSum",&neSum_F);
  
  jetpbpb[2]->SetBranchAddress("eSum",eSum_F);
  jetpbpb[2]->SetBranchAddress("eMax",eMax_F);
  jetpbpb[2]->SetBranchAddress("muSum",muSum_F);
  jetpbpb[2]->SetBranchAddress("muMax",muMax_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet40_v1",&jet40_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet40_v1_Prescl",&jet40_p_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet60_v1",&jet60_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet60_v1_Prescl",&jet60_p_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet80_v1",&jet80_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet80_v1_Prescl",&jet80_p_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet36_BptxAND",&L1_sj36_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&L1_sj36_p_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet52_BptxAND",&L1_sj52_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&L1_sj52_p_F);

  // TFile * fJERJES = TFile::Open("JERJES_PF_PbPb_pp.root");

  // TH1F * hJES_mean[nbins_cent];
  // hJES_mean[0] = (TH1F*)fJERJES->Get(Form("ak%dJetAnalyzer/hMean_PbPb_R%d_PuPF_05",radius, radius));
  // hJES_mean[1] = (TH1F*)fJERJES->Get(Form("ak%dJetAnalyzer/hMean_PbPb_R%d_PuPF_510",radius, radius));
  // hJES_mean[2] = (TH1F*)fJERJES->Get(Form("ak%dJetAnalyzer/hMean_PbPb_R%d_PuPF_1030",radius, radius));
  // hJES_mean[3] = (TH1F*)fJERJES->Get(Form("ak%dJetAnalyzer/hMean_PbPb_R%d_PuPF_3050",radius, radius));
  // hJES_mean[4] = (TH1F*)fJERJES->Get(Form("ak%dJetAnalyzer/hMean_PbPb_R%d_PuPF_5070",radius, radius));
  // hJES_mean[5] = (TH1F*)fJERJES->Get(Form("ak%dJetAnalyzer/hMean_PbPb_R%d_PuPF_7090",radius, radius));  
  
  // Declare the output File and the necessary histograms after that:
  // std::string outdir="";
  // std::string outfile=outdir+kFoname;
  TFile *fout = new TFile(kFoname.c_str(),"RECREATE");
  fout->cd();

  // TH1F * hpbpb_Jet80[nbins_cent];
  // TH1F * hpbpb_Jet60[nbins_cent];
  // TH1F * hpbpb_Jet40[nbins_cent];
  // TH1F * hpbpb_JetComb[nbins_cent];

  // TH1F * hpbpb_JEC_Jet80[nbins_cent];
  // TH1F * hpbpb_JEC_Jet60[nbins_cent];
  // TH1F * hpbpb_JEC_Jet40[nbins_cent];
  // TH1F * hpbpb_JEC_JetComb[nbins_cent];

  // TH1F * hpbpb_JEC_minus_Jet80[nbins_cent];
  // TH1F * hpbpb_JEC_minus_Jet60[nbins_cent];
  // TH1F * hpbpb_JEC_minus_Jet40[nbins_cent];
  // TH1F * hpbpb_JEC_minus_JetComb[nbins_cent];


  // TH1F * hpbpb_JEC_gaus_Jet80[nbins_cent];
  // TH1F * hpbpb_JEC_gaus_Jet60[nbins_cent];
  // TH1F * hpbpb_JEC_gaus_Jet40[nbins_cent];
  // TH1F * hpbpb_JEC_gaus_JetComb[nbins_cent];
  
  // TH1F * hpbpb_Smear_Jet80[nbins_cent];
  // TH1F * hpbpb_Smear_Jet60[nbins_cent];
  // TH1F * hpbpb_Smear_Jet40[nbins_cent];
  // TH1F * hpbpb_Smear_JetComb[nbins_cent];

  std::string var[19] = {"jtpt" ,"rawpt", "jteta", "jtphi", "trkMax", "trkSum", "trkHardSum", "chMax", "chSum", "chHardSum","phMax", "phSum", "phHardSum", "neMax", "neSum", "eMax", "eSum", "muMax", "muSum" };
  TH1F * hJetQA[2][19][nbins_cent];
  
  for(int i = 0;i<nbins_cent;++i){

    /*
    hpbpb_Jet80[i] = new TH1F(Form("hpbpb_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_Jet60[i] = new TH1F(Form("hpbpb_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_Jet40[i] = new TH1F(Form("hpbpb_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet60 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_JetComb[i] = new TH1F(Form("hpbpb_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);

    hpbpb_JEC_Jet80[i] = new TH1F(Form("hpbpb_JEC_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_JEC_Jet60[i] = new TH1F(Form("hpbpb_JEC_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_JEC_Jet40[i] = new TH1F(Form("hpbpb_JEC_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet60 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_JEC_JetComb[i] = new TH1F(Form("hpbpb_JEC_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);

    hpbpb_JEC_minus_Jet80[i] = new TH1F(Form("hpbpb_JEC_minus_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_JEC_minus_Jet60[i] = new TH1F(Form("hpbpb_JEC_minus_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_JEC_minus_Jet40[i] = new TH1F(Form("hpbpb_JEC_minus_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet60 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_JEC_minus_JetComb[i] = new TH1F(Form("hpbpb_JEC_minus_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);

    hpbpb_JEC_gaus_Jet80[i] = new TH1F(Form("hpbpb_JEC_gaus_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_JEC_gaus_Jet60[i] = new TH1F(Form("hpbpb_JEC_gaus_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_JEC_gaus_Jet40[i] = new TH1F(Form("hpbpb_JEC_gaus_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet60 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_JEC_gaus_JetComb[i] = new TH1F(Form("hpbpb_JEC_gaus_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);

    
    hpbpb_Smear_Jet80[i] = new TH1F(Form("hpbpb_Smear_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_Smear_Jet60[i] = new TH1F(Form("hpbpb_Smear_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_Smear_Jet40[i] = new TH1F(Form("hpbpb_Smear_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet60 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_Smear_JetComb[i] = new TH1F(Form("hpbpb_Smear_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    */
    
    for(int k = 0; k<2; ++k){
      for(int j = 0; j<19; ++j){
	if(j==2) hJetQA[k][j][i] = new TH1F(Form("hJetQA_%dwJetID_%s_cent%d",k,var[j].c_str(),i),Form(";%s;",var[j].c_str()),100, -5, +5);
	else if(j==3) hJetQA[k][j][i] = new TH1F(Form("hJetQA_%dwJetID_%s_cent%d",k,var[j].c_str(),i),Form(";%s;",var[j].c_str()),100, -4, +4);
	else if(j<=1)hJetQA[k][j][i] = new TH1F(Form("hJetQA_%dwJetID_%s_cent%d",k,var[j].c_str(),i),Form(";%s;",var[j].c_str()),500, 0, 500);
	else if(j>=4)hJetQA[k][j][i] = new TH1F(Form("hJetQA_%dwJetID_%s_cent%d",k,var[j].c_str(),i),Form(";%s;",var[j].c_str()),200, 0, 2);
      }
    }
    
  }


  
  // float ue_fluctuation[nbins_cent] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  // if(radius == 4) {
  //   ue_fluctuation[0] = 5.0;
  //   ue_fluctuation[1] = 1.5;
  //   ue_fluctuation[2] = 0.3;
  //   ue_fluctuation[3] = 0.1;
  //   ue_fluctuation[4] = 0.1;
  //   ue_fluctuation[5] = 0.01;
  // }
  // if(radius == 3) {
  //   ue_fluctuation[0] = 3.17;
  //   ue_fluctuation[1] = 0.7;
  //   ue_fluctuation[2] = 0.1;
  //   ue_fluctuation[3] = 0.1;
  //   ue_fluctuation[4] = 0.1;
  //   ue_fluctuation[5] = 0.01;
  // }
  // if(radius == 2) {
  //   ue_fluctuation[0] = 1.5;
  //   ue_fluctuation[1] = 0.4;
  //   ue_fluctuation[2] = 0.2;
  //   ue_fluctuation[3] = 0.1;
  //   ue_fluctuation[4] = 0.1;
  //   ue_fluctuation[5] = 0.01;
  // }

  // float  ue_fluctuation = 6.0;
  // if(radius == 4) ue_fluctuation = 6.0;
  // if(radius == 3) ue_fluctuation = 4.0;
  // if(radius == 2) ue_fluctuation = 2.0;

  // the UE Smear corresponds to doing recpt = recpt * (1 + rnd.Gaus(0, (float)ue_fluctuation[cBin]/recpt);
  
  // now start the event loop for each file. 
  if(printDebug) cout<<"Running through all the events now"<<endl;
  Long64_t nentries = jetpbpb[0]->GetEntries();
  if(printDebug) nentries = 10;
  TRandom3 rnd;
  rnd.SetSeed(endfile);

  for(int nEvt = 0; nEvt < nentries; ++ nEvt) {

    if(nEvt%10000 == 0)cout<<nEvt<<"/"<<nentries<<endl;
    if(printDebug)cout<<"nEvt = "<<nEvt<<endl;
    
    jetpbpb[0]->GetEntry(nEvt);
    jetpbpb[1]->GetEntry(nEvt);
    jetpbpb[2]->GetEntry(nEvt);
    jetpbpb[3]->GetEntry(nEvt);

    if(printDebug) cout<<"forest values = "<<hiBin_F<<", "<<evt_F<<", "<<run_F<<", "<<lumi_F<<", "<<vz_F<<endl;
    
    if(pHBHENoiseFilter_F == 0) continue; 
    if(pcollisionEventSelection_F == 0) continue; 
    if(fabs(vz_F)>15) continue;
    // if(!isGoodEvent_eS) continue; 

    int cBin = findBin(hiBin_F);//tells us the centrality of the event.
    if(cBin == -1) continue;

    // int jetCounter = 0;
    
    // for(int g = 0;g<nref_F;g++){
      
    //   if(eta_F[g]>=-2 && eta_F[g]<2){ //to select inside 
	
    // 	if(pt_F[g]>=50) jetCounter++;
	
    //   }//eta selection cut
      
    // }// jet loop

    
    // jet_select->GetEntry(nGoodEvt);
    // ++nGoodEvt;
    Float_t weight_eS = 1.0; 

    //if(cBin == nbins_cent) continue;
    // // apply the correct supernova selection cut rejection here: 
    // if(hiNpix_F >= (38000 - 500*jetCounter)){
    //   if(printDebug) cout<<"removed this supernova event"<<endl;
    //   continue;
    // }

    // if((float)hiZDCminus_F/1350 >= (90 - 0.0204 * hiHF_F)) continue;

    // //! Sort the jetTree jets according to pT
    // std::vector < Jet > vJet;
    // for(int jet2 = 0; jet2<nref_jS; ++jet2){
    //   //cout <<"\t \t jetTree *** "<< jet2 <<  ", pT " << pt_jS[jet2] <<  ", chSum : "<< chSum_jS[jet2] << endl;
    //   Jet ijet;
    //   ijet.id = jet2;
    //   ijet.pt = pt_jS[jet2];
    //   vJet.push_back(ijet);
    // }
    // std::sort (vJet.begin(), vJet.end(), compare_pt);
    // std::vector < Jet >::const_iterator itr;

    // int jet=0;
    // for(itr=vJet.begin(); itr!=vJet.end(); ++itr, ++jet){

    //   int jetLoc = (*itr).id;
    //   if(isMultiMatch_jS[jetLoc]) {
    // 	++itr;
    // 	jetLoc = (*itr).id;
    // 	if(itr == vJet.end())  break;
    //   }
    //   if(fabs(eta_jS[jetLoc]) > 2) continue;
    //   //if(isPFElecCut_eS[jet] != 1) continue;
    //   // if(isMiMatch_eS[jet]) continue;
    //   if(pt_jS[jetLoc] <15) continue;

    //   bool PFElecCut = false;

    //   Float_t Sumcand = chSum_jS[jetLoc] + phSum_jS[jetLoc] + neSum_jS[jetLoc] + muSum_jS[jetLoc];
    //   if(isCaloMatch_jS[jetLoc] == 1){
    // 	if(calopt_jS[jetLoc]/pt_jS[jetLoc] > 0.5 && calopt_jS[jetLoc]/pt_jS[jetLoc] <= 0.85 && eMax_jS[jetLoc]/Sumcand < ((Float_t)18/7 *(Float_t)calopt_jS[jetLoc]/pt_jS[jetLoc] - (Float_t)9/7)) PFElecCut = true;
    // 	if(calopt_jS[jetLoc]/pt_jS[jetLoc] > 0.85) PFElecCut = true;
    // 	if(calopt_jS[jetLoc]/pt_jS[jetLoc] <= 0.5 && eMax_jS[jetLoc]/Sumcand < 0.05) PFElecCut = true;
    //   }
    //   if(isCaloMatch_jS[jetLoc] == 0)
    // 	if(eMax_jS[jetLoc]/Sumcand < 0.05 ) PFElecCut = true;

    //   // if(!PFElecCut) continue;
      
    //   // if(printDebug && (fabs(eta_jS[jet] > 2))) cout<<"jets with |eta| > 2 in jetTree"<<endl;
    //   // if(printDebug && (fabs(eta_F[jet] > 2)))  cout<<"jets with |eta| > 2 in Forest"<<endl;

      
    //   if(printDebug && index_eS[jet] >= 0 )cout<<jet<<", hiForest pT = "<<pt_F[jet]<<", jetTree pT = "<<pt_jS[jetLoc]<<", electronCut = "<<isPFElecCut_eS[jetLoc]<<", Calo pT = "<<calopt_F[index_eS[jet]]<<", onFly flag calculation = "<<PFElecCut<<", eMax from hiForest = "<<eMax_F[jet]<<", eMax from jet Tree = "<<eMax_jS[jetLoc]<<endl;      // if(printDebug)cout<<jet<<", hiForest pT = "<<pt_F[jet]<<", jetTree pT = "<<pt_jS[jetLoc]<<", electronCut = "<<isPFElecCut_eS[jetLoc]<<", eMax from hiForest = "<<eMax_F[jet]<<", eMax from jet Tree = "<<eMax_jS[jetLoc]<<endl;

    for( int jet = 0; jet<nref_F; jet++ ){

      if( fabs(eta_F[jet]) > 2.0 ) continue;
      //if( isPFElecCut_eS[jet] != 1 ) continue;

      //if( chMax_F[jet] < 7 && trkMax_F[jet] < 7 && neMax_F[jet] < 8 ) continue;
      // if( trkMax_F[jet] < 7 && neMax_F[jet] < 8 ) continue;

      // bool PFElecCut = false;
      // Float_t calopfpt = (float) calopt_eS[jet]/pt_F[jet];
      // Float_t Sumcand = chSum_F[jet] + phSum_F[jet] + neSum_F[jet] + muSum_F[jet];
      float recpt     = pt_F[jet];
      float rawpt = rawpt_F[jet];
      if(recpt<=30) continue;

// #if 0
//       int recpt_loc = 0;      
//       for(int j = 0; j<ptSelection; ++j)
// 	if(recpt > ptBoundary[j]) recpt_loc = j;      
//       if(isClMatch_eS[jet] == 1) heMaxSumcand_vs_calopfpt_RecopTSelection[recpt_loc][cBin]->Fill(calopfpt, (float)eMax_F[jet]/Sumcand);
//       else heMaxSumcand_vs_calopfpt_RecopTSelection[recpt_loc][cBin]->Fill(0.0, (float)eMax_F[jet]/Sumcand);
// #endif
      
//       if(isClMatch_eS[jet] == 1){
//     	if(calopfpt > 0.5 && calopfpt <= 0.85 && eMax_F[jet]/Sumcand < ((Float_t)18/7 *(Float_t)calopfpt - (Float_t)9/7)) PFElecCut = true;
//     	if(calopfpt > 0.85) PFElecCut = true;
//     	if(calopfpt <= 0.5 && eMax_F[jet]/Sumcand < 0.05) PFElecCut = true;
//       }
//       if(isClMatch_eS[jet] == 0)
//     	if(eMax_F[jet]/Sumcand < 0.05 ) PFElecCut = true;

//       if(isPFElecCut_eS[jet] != PFElecCut && printDebug) 
// 	cout<<" pf e cut not same, run = "<<run_F<<" "<<run_eS<<", event = "<<evt_F<<" "<<evt_eS<<" , lumi = "<<lumi_F<<" "<<lumi_eS<<endl;
	  
//       if(!PFElecCut) continue;

      // also need to cut on the high pT jets from the Jet 55 sample. unmatched PF jets with pfpt > 110 and calopt < 30 including -999, check on their high track content. at the moment dont worry about this.  

      // float JEC_Smear = recpt;
      // float JEC_minus_Smear = recpt;
      // float JEC_gaus_Smear = recpt;
      // float UE_Smear = recpt;
      // // this is some percentage (+ or -) from 0 to 1-JES
      // //float RanJES = (float)(rnd.Gaus(0.0,fabs(1.-hJES_mean[cBin]->GetBinContent(hJES_mean[cBin]->FindBin(recpt)))));
      // // this is some scale up and down (+ or -) from 0 to ue_fluct, so if ue_fluct is 10, the jet can change by 1sigma -10 or +10
      // float RanUE =  (float)(rnd.Gaus(0.0,(float)ue_fluctuation[cBin]));
      // UE_Smear = recpt+RanUE;
      // float JEC_5pcent_Smear = (float)recpt * (1. + (float)0.05/nbins_cent*(nbins_cent-cBin));
      // if(radius == 2 )JEC_Smear = (float)recpt * (1. + (float)0.01/nbins_cent*(nbins_cent-cBin));
      // if(radius == 3 )JEC_Smear = (float)recpt * (1. + (float)0.01/nbins_cent*(nbins_cent-cBin));
      // if(radius == 4 )JEC_Smear = (float)recpt * (1. + (float)0.01/nbins_cent*(nbins_cent-cBin));

      // float Shift = 0.0;
      // for(int bin = 0; bin<nbins_jec; ++bin)
      // 	if(recpt > ptbins_jec[bin])
      // 	  Shift = jec_shift[bin];
      
      // JEC_Smear = (float)recpt * (1.0 + Shift);
      // JEC_minus_Smear = (float)recpt * (1.0 - Shift);
      
      // JEC_gaus_Smear = (float)recpt * (rnd.Gaus(1, Shift));

      //JEC_Smear = recpt*(1 + RanJES);      
      // if(recpt >= 105 && recpt < 106){	
      // 	hJet_100GeV_input[cBin]->Fill(recpt);
      // 	hJet_100GeV_JECSmear[cBin]->Fill(JEC_Smear);
      // 	hJet_100GeV_JEC_5pcentSmear[cBin]->Fill(JEC_5pcent_Smear);
      // 	hJet_100GeV_UESmear[cBin]->Fill(UE_Smear);
      // }
 
      // float recpt     = pt_F[jet];
      // float JEC_Smear = recpt;
      // float UE_Smear = recpt;
      // float pm = (float)(rnd.Gaus(0.0, 1.0));
      // if(pm > 0) {
      // 	JEC_Smear = (float)(recpt * (1. + (float)(rnd.Gaus(0.0, (float)fabs(1. - hJES_mean[cBin]->GetBinContent(hJES_mean[cBin]->FindBin(recpt)))))));
      // 	UE_Smear  = (float)(recpt * (1. + (float)(rnd.Gaus(0.0, (float)ue_fluctuation[cBin]/recpt))));
      // }
      // if(pm < 0){
      // 	JEC_Smear = (float)(recpt * (1. - (float)(rnd.Gaus(0.0, (float)fabs(1. - hJES_mean[cBin]->GetBinContent(hJES_mean[cBin]->FindBin(recpt)))))));
      // 	UE_Smear  = (float)(recpt * (1. - (float)(rnd.Gaus(0.0, (float)ue_fluctuation[cBin]/recpt))));
      // }
      
      // if(jet60_F == 1) hpbpb_HLT65[cBin]->Fill(recpt);
      // if(jet40_F == 1) hpbpb_HLT55[cBin]->Fill(recpt);
      
      // if(jet40_F == 1 && jet60_F == 0 && jet80_F == 0){
      // 	//if(recpt > 140) continue;
      // 	hpbpb_Jet40[cBin]->Fill(recpt, jet40_p_F* wght);
      // 	hpbpb_Smear_Jet40[cBin]->Fill(UE_Smear, jet40_p_F* wght);
      // 	hpbpb_JEC_Jet40[cBin]->Fill(JEC_Smear, jet40_p_F* wght);
      // 	hpbpb_JEC_minus_Jet40[cBin]->Fill(JEC_minus_Smear, jet40_p_F* wght);
      // 	hpbpb_JEC_gaus_Jet40[cBin]->Fill(JEC_gaus_Smear, jet40_p_F* wght);

      // 	// hpbpb_TrgObj55[cBin]->Fill(recpt, jet40_p_F* wght);
      // 	// hpbpb_raw_TrgObj55[cBin]->Fill(rawpt_F[jet], jet40_p_F* wght);
      // 	// hpbpb_anaBin_TrgObj55[cBin]->Fill(recpt, jet40_p_F* wght);
      // 	// hpbpb_JEC_TrgObj55[cBin]->Fill(JEC_Smear, jet40_p_F* wght);
      // 	// hpbpb_JEC_5pcent_TrgObj55[cBin]->Fill(JEC_5pcent_Smear, jet40_p_F* wght);
      // 	// hpbpb_Smear_TrgObj55[cBin]->Fill(UE_Smear, jet40_p_F* wght);
      // }
      // if(jet60_F == 1 && jet80_F == 0){
      // 	//if(recpt > 140) continue;
      // 	hpbpb_Jet60[cBin]->Fill(recpt,wght);
      // 	hpbpb_Smear_Jet60[cBin]->Fill(UE_Smear, wght);
      // 	hpbpb_JEC_Jet60[cBin]->Fill(JEC_Smear, wght);
      // 	hpbpb_JEC_minus_Jet60[cBin]->Fill(JEC_minus_Smear, wght);
      // 	hpbpb_JEC_gaus_Jet60[cBin]->Fill(JEC_gaus_Smear, wght);

      // 	// hpbpb_TrgObj65[cBin]->Fill(recpt, wght);
      // 	// hpbpb_raw_TrgObj65[cBin]->Fill(rawpt_F[jet], wght);
      // 	// hpbpb_anaBin_TrgObj65[cBin]->Fill(recpt, wght);
      // 	// hpbpb_JEC_TrgObj65[cBin]->Fill(JEC_Smear, wght);
      // 	// hpbpb_JEC_5pcent_TrgObj65[cBin]->Fill(JEC_5pcent_Smear, wght);
      // 	// hpbpb_Smear_TrgObj65[cBin]->Fill(UE_Smear, wght);
      // }
      // if(jet80_F == 1){
      // 	hpbpb_Jet80[cBin]->Fill(recpt, wght);
      // 	hpbpb_Smear_Jet80[cBin]->Fill(UE_Smear,wght);
      // 	hpbpb_JEC_Jet80[cBin]->Fill(JEC_Smear, wght);
      // 	hpbpb_JEC_minus_Jet80[cBin]->Fill(JEC_minus_Smear, wght);
      // 	hpbpb_JEC_gaus_Jet80[cBin]->Fill(JEC_gaus_Smear, wght);

      // 	// hpbpb_TrgObj80[cBin]->Fill(recpt, wght);
      // 	// hpbpb_raw_TrgObj80[cBin]->Fill(rawpt_F[jet], wght);
      // 	// hpbpb_anaBin_TrgObj80[cBin]->Fill(recpt, wght);
      // 	// hpbpb_JEC_TrgObj80[cBin]->Fill(JEC_Smear, wght);
      // 	// hpbpb_JEC_5pcent_TrgObj80[cBin]->Fill(JEC_5pcent_Smear, wght);
      // 	// hpbpb_Smear_TrgObj80[cBin]->Fill(UE_Smear, wght);
      // }


      if(1){
	hJetQA[0][0][cBin]->Fill(recpt, weight_eS);
	hJetQA[0][1][cBin]->Fill(rawpt_F[jet], weight_eS);
	hJetQA[0][2][cBin]->Fill(eta_F[jet], weight_eS);
	hJetQA[0][3][cBin]->Fill(phi_F[jet], weight_eS);
	hJetQA[0][4][cBin]->Fill(trkMax_F[jet]/recpt, weight_eS);
	hJetQA[0][5][cBin]->Fill(trkSum_F[jet]/recpt, weight_eS);
	hJetQA[0][6][cBin]->Fill(trkHardSum_F[jet]/recpt, weight_eS);
	hJetQA[0][7][cBin]->Fill(chMax_F[jet]/recpt, weight_eS);
	hJetQA[0][8][cBin]->Fill(chSum_F[jet]/recpt, weight_eS);
	hJetQA[0][9][cBin]->Fill(chHardSum_F[jet]/recpt, weight_eS);
	hJetQA[0][10][cBin]->Fill(phMax_F[jet]/recpt, weight_eS);
	hJetQA[0][11][cBin]->Fill(phSum_F[jet]/recpt, weight_eS);
	hJetQA[0][12][cBin]->Fill(phHardSum_F[jet]/recpt, weight_eS);
	hJetQA[0][13][cBin]->Fill(neMax_F[jet]/recpt, weight_eS);
	hJetQA[0][14][cBin]->Fill(neSum_F[jet]/recpt, weight_eS);
	hJetQA[0][15][cBin]->Fill(eMax_F[jet]/recpt, weight_eS);
	hJetQA[0][16][cBin]->Fill(eSum_F[jet]/recpt, weight_eS);
	hJetQA[0][17][cBin]->Fill(muMax_F[jet]/recpt, weight_eS);
	hJetQA[0][18][cBin]->Fill(muSum_F[jet]/recpt, weight_eS);

	// apply JetID
	// charged hadron fraction > 0
	// charged hadron multiplicity > 0
	// charged EM fraction < 0.99 
	
	if(chSum_F[jet]/recpt<0.99 && chSum_F[jet]/recpt>0 && neSum_F[jet]/recpt<0.99 && chN_F[jet]>0 && eSum_F[jet]/recpt<0.99){
	  hJetQA[1][0][cBin]->Fill(recpt, weight_eS);
	  hJetQA[1][1][cBin]->Fill(rawpt_F[jet], weight_eS);
	  hJetQA[1][2][cBin]->Fill(eta_F[jet], weight_eS);
	  hJetQA[1][3][cBin]->Fill(phi_F[jet], weight_eS);
	  hJetQA[1][4][cBin]->Fill(trkMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][5][cBin]->Fill(trkSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][6][cBin]->Fill(trkHardSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][7][cBin]->Fill(chMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][8][cBin]->Fill(chSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][9][cBin]->Fill(chHardSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][10][cBin]->Fill(phMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][11][cBin]->Fill(phSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][12][cBin]->Fill(phHardSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][13][cBin]->Fill(neMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][14][cBin]->Fill(neSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][15][cBin]->Fill(eMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][16][cBin]->Fill(eSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][17][cBin]->Fill(muMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][18][cBin]->Fill(muSum_F[jet]/recpt, weight_eS);
	}
	
      }

      
    }// jet loop
    if(printDebug)cout<<endl;


  }// event loop


  /*
  for(int i = 0; i<nbins_cent; ++i){


    hpbpb_Jet80[i]->Scale(1./ScaleFactor_80_NeqScale[i]);
    hpbpb_Jet60[i]->Scale(1./ScaleFactor_65_NeqScale[i]);
    hpbpb_Jet40[i]->Scale(1./ScaleFactor_55_NeqScale[i]);
    hpbpb_JetComb[i]->Add(hpbpb_Jet80[i]);
    hpbpb_JetComb[i]->Add(hpbpb_Jet60[i]);
    hpbpb_JetComb[i]->Add(hpbpb_Jet40[i]);
    
    divideBinWidth(hpbpb_JetComb[i]);
    divideBinWidth(hpbpb_Jet80[i]);
    divideBinWidth(hpbpb_Jet60[i]);
    divideBinWidth(hpbpb_Jet40[i]);

    hpbpb_JEC_Jet80[i]->Scale(1./ScaleFactor_80_NeqScale[i]);
    hpbpb_JEC_Jet60[i]->Scale(1./ScaleFactor_65_NeqScale[i]);
    hpbpb_JEC_Jet40[i]->Scale(1./ScaleFactor_55_NeqScale[i]);
    hpbpb_JEC_JetComb[i]->Add(hpbpb_JEC_Jet80[i]);
    hpbpb_JEC_JetComb[i]->Add(hpbpb_JEC_Jet60[i]);
    hpbpb_JEC_JetComb[i]->Add(hpbpb_JEC_Jet40[i]);
    
    divideBinWidth(hpbpb_JEC_JetComb[i]);
    divideBinWidth(hpbpb_JEC_Jet80[i]);
    divideBinWidth(hpbpb_JEC_Jet60[i]);
    divideBinWidth(hpbpb_JEC_Jet40[i]);

    hpbpb_JEC_minus_Jet80[i]->Scale(1./ScaleFactor_80_NeqScale[i]);
    hpbpb_JEC_minus_Jet60[i]->Scale(1./ScaleFactor_65_NeqScale[i]);
    hpbpb_JEC_minus_Jet40[i]->Scale(1./ScaleFactor_55_NeqScale[i]);
    hpbpb_JEC_minus_JetComb[i]->Add(hpbpb_JEC_minus_Jet80[i]);
    hpbpb_JEC_minus_JetComb[i]->Add(hpbpb_JEC_minus_Jet60[i]);
    hpbpb_JEC_minus_JetComb[i]->Add(hpbpb_JEC_minus_Jet40[i]);
    
    divideBinWidth(hpbpb_JEC_minus_JetComb[i]);
    divideBinWidth(hpbpb_JEC_minus_Jet80[i]);
    divideBinWidth(hpbpb_JEC_minus_Jet60[i]);
    divideBinWidth(hpbpb_JEC_minus_Jet40[i]);

    hpbpb_JEC_gaus_Jet80[i]->Scale(1./ScaleFactor_80_NeqScale[i]);
    hpbpb_JEC_gaus_Jet60[i]->Scale(1./ScaleFactor_65_NeqScale[i]);
    hpbpb_JEC_gaus_Jet40[i]->Scale(1./ScaleFactor_55_NeqScale[i]);
    hpbpb_JEC_gaus_JetComb[i]->Add(hpbpb_JEC_gaus_Jet80[i]);
    hpbpb_JEC_gaus_JetComb[i]->Add(hpbpb_JEC_gaus_Jet60[i]);
    hpbpb_JEC_gaus_JetComb[i]->Add(hpbpb_JEC_gaus_Jet40[i]);
    
    divideBinWidth(hpbpb_JEC_gaus_JetComb[i]);
    divideBinWidth(hpbpb_JEC_gaus_Jet80[i]);
    divideBinWidth(hpbpb_JEC_gaus_Jet60[i]);
    divideBinWidth(hpbpb_JEC_gaus_Jet40[i]);
    
    hpbpb_Smear_Jet80[i]->Scale(1./ScaleFactor_80_NeqScale[i]);
    hpbpb_Smear_Jet60[i]->Scale(1./ScaleFactor_65_NeqScale[i]);
    hpbpb_Smear_Jet40[i]->Scale(1./ScaleFactor_55_NeqScale[i]);
    hpbpb_Smear_JetComb[i]->Add(hpbpb_Smear_Jet80[i]);
    hpbpb_Smear_JetComb[i]->Add(hpbpb_Smear_Jet60[i]);
    hpbpb_Smear_JetComb[i]->Add(hpbpb_Smear_Jet40[i]);
    
    divideBinWidth(hpbpb_Smear_JetComb[i]);
    divideBinWidth(hpbpb_Smear_Jet80[i]);
    divideBinWidth(hpbpb_Smear_Jet60[i]);
    divideBinWidth(hpbpb_Smear_Jet40[i]);

    
    // hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj80[i]);
    // hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj65[i]);
    // hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj55[i]);
    
    // divideBinWidth(hpbpb_TrgObjComb[i]);
    // divideBinWidth(hpbpb_TrgObj80[i]);
    // divideBinWidth(hpbpb_TrgObj65[i]);
    // divideBinWidth(hpbpb_TrgObj55[i]);
    
    
    // hpbpb_NoPileup_TrgObjComb[i]->Add(hpbpb_NoPileup_TrgObj80[i]);
    // hpbpb_NoPileup_TrgObjComb[i]->Add(hpbpb_NoPileup_TrgObj65[i]);
    // hpbpb_NoPileup_TrgObjComb[i]->Add(hpbpb_NoPileup_TrgObj55[i]);
    
    // divideBinWidth(hpbpb_NoPileup_TrgObjComb[i]);
    // divideBinWidth(hpbpb_NoPileup_TrgObj80[i]);
    // divideBinWidth(hpbpb_NoPileup_TrgObj65[i]);
    // divideBinWidth(hpbpb_NoPileup_TrgObj55[i]);
    
    // hpbpb_944Scale_TrgObj80[i] = (TH1F*)hpbpb_TrgObj80[i]->Clone(Form("hpbpb_HLT80_944Scale_R%d_20_eta_20_cent%d", radius, i));
    // hpbpb_944Scale_TrgObj65[i] = (TH1F*)hpbpb_TrgObj65[i]->Clone(Form("hpbpb_HLT65_944Scale_R%d_20_eta_20_cent%d", radius, i));
    // hpbpb_944Scale_TrgObj55[i] = (TH1F*)hpbpb_TrgObj55[i]->Clone(Form("hpbpb_HLT55_944Scale_R%d_20_eta_20_cent%d", radius, i));

    // hpbpb_944Scale_TrgObj80[i]->Scale(1./ScaleFactor_944Scale[i]);
    // hpbpb_944Scale_TrgObj65[i]->Scale(1./ScaleFactor_944Scale[i]);
    // hpbpb_944Scale_TrgObj55[i]->Scale(1./ScaleFactor_944Scale[i]);
    // hpbpb_944Scale_TrgObjComb[i] = (TH1F*)hpbpb_944Scale_TrgObj80[i]->Clone(Form("hpbpb_HLTComb_944Scale_R%d_20_eta_20_cent%d",radius, i));
    // hpbpb_944Scale_TrgObjComb[i]->Add(hpbpb_944Scale_TrgObj65[i]);
    // hpbpb_944Scale_TrgObjComb[i]->Add(hpbpb_944Scale_TrgObj55[i]);
    
    // hpbpb_NeqScalePerCent_TrgObj80[i] = (TH1F*)hpbpb_TrgObj80[i]->Clone(Form("hpbpb_HLT80_NeqScalePerCent_R%d_20_eta_20_cent%d", radius, i));
    // hpbpb_NeqScalePerCent_TrgObj65[i] = (TH1F*)hpbpb_TrgObj65[i]->Clone(Form("hpbpb_HLT65_NeqScalePerCent_R%d_20_eta_20_cent%d", radius, i));
    // hpbpb_NeqScalePerCent_TrgObj55[i] = (TH1F*)hpbpb_TrgObj55[i]->Clone(Form("hpbpb_HLT55_NeqScalePerCent_R%d_20_eta_20_cent%d", radius, i));
    // hpbpb_NeqScalePerCent_TrgObj80[i]->Scale(1./ScaleFactor_80_NeqScalePerCent[i]);
    // hpbpb_NeqScalePerCent_TrgObj65[i]->Scale(1./ScaleFactor_65_NeqScalePerCent[i]);
    // hpbpb_NeqScalePerCent_TrgObj55[i]->Scale(1./ScaleFactor_55_NeqScalePerCent[i]);
    // hpbpb_NeqScalePerCent_TrgObjComb[i] = (TH1F*)hpbpb_NeqScalePerCent_TrgObj80[i]->Clone(Form("hpbpb_HLTComb_NeqScalePerCent_R%d_20_eta_20_cent%d",radius, i));
    // hpbpb_NeqScalePerCent_TrgObjComb[i]->Add(hpbpb_NeqScalePerCent_TrgObj65[i]);
    // hpbpb_NeqScalePerCent_TrgObjComb[i]->Add(hpbpb_NeqScalePerCent_TrgObj55[i]);

    
    // hpbpb_NeqScale_TrgObj80[i] = (TH1F*)hpbpb_TrgObj80[i]->Clone(Form("hpbpb_HLT80_NeqScale_R%d_20_eta_20_cent%d", radius, i));
    // hpbpb_NeqScale_TrgObj65[i] = (TH1F*)hpbpb_TrgObj65[i]->Clone(Form("hpbpb_HLT65_NeqScale_R%d_20_eta_20_cent%d", radius, i));
    // hpbpb_NeqScale_TrgObj55[i] = (TH1F*)hpbpb_TrgObj55[i]->Clone(Form("hpbpb_HLT55_NeqScale_R%d_20_eta_20_cent%d", radius, i));
    // hpbpb_NeqScale_TrgObj80[i]->Scale(1./ScaleFactor_80_NeqScale[i]);
    // hpbpb_NeqScale_TrgObj65[i]->Scale(1./ScaleFactor_65_NeqScale[i]);
    // hpbpb_NeqScale_TrgObj55[i]->Scale(1./ScaleFactor_55_NeqScale[i]);
    // hpbpb_NeqScale_TrgObjComb[i] = (TH1F*)hpbpb_NeqScale_TrgObj80[i]->Clone(Form("hpbpb_HLTComb_NeqScale_R%d_20_eta_20_cent%d",radius, i));
    // hpbpb_NeqScale_TrgObjComb[i]->Add(hpbpb_NeqScale_TrgObj65[i]);
    // hpbpb_NeqScale_TrgObjComb[i]->Add(hpbpb_NeqScale_TrgObj55[i]);

    
    // hpbpb_raw_TrgObj80[i]->Scale(1./ScaleFactor_80_NeqScale[i]);
    // hpbpb_raw_TrgObj65[i]->Scale(1./ScaleFactor_65_NeqScale[i]);
    // hpbpb_raw_TrgObj55[i]->Scale(1./ScaleFactor_55_NeqScale[i]);
    
    // hpbpb_raw_TrgObjComb[i]->Add(hpbpb_raw_TrgObj80[i]);
    // hpbpb_raw_TrgObjComb[i]->Add(hpbpb_raw_TrgObj65[i]);
    // hpbpb_raw_TrgObjComb[i]->Add(hpbpb_raw_TrgObj55[i]);

    // divideBinWidth(hpbpb_raw_TrgObjComb[i]);
    // divideBinWidth(hpbpb_raw_TrgObj80[i]);
    // divideBinWidth(hpbpb_raw_TrgObj65[i]);
    // divideBinWidth(hpbpb_raw_TrgObj55[i]);

    // hpbpb_anaBin_TrgObj80[i]->Scale(1./ScaleFactor_80_NeqScale[i]);
    // hpbpb_anaBin_TrgObj65[i]->Scale(1./ScaleFactor_65_NeqScale[i]);
    // hpbpb_anaBin_TrgObj55[i]->Scale(1./ScaleFactor_55_NeqScale[i]);

    // hpbpb_anaBin_TrgObjComb[i]->Add(hpbpb_anaBin_TrgObj80[i]);
    // hpbpb_anaBin_TrgObjComb[i]->Add(hpbpb_anaBin_TrgObj65[i]);
    // hpbpb_anaBin_TrgObjComb[i]->Add(hpbpb_anaBin_TrgObj55[i]);

    // divideBinWidth(hpbpb_anaBin_TrgObjComb[i]);
    // divideBinWidth(hpbpb_anaBin_TrgObj80[i]);
    // divideBinWidth(hpbpb_anaBin_TrgObj65[i]);
    // divideBinWidth(hpbpb_anaBin_TrgObj55[i]);

    // hpbpb_Smear_TrgObj80[i]->Scale(1./ScaleFactor_80_NeqScale[i]);
    // hpbpb_Smear_TrgObj65[i]->Scale(1./ScaleFactor_65_NeqScale[i]);
    // hpbpb_Smear_TrgObj55[i]->Scale(1./ScaleFactor_55_NeqScale[i]);

    // hpbpb_Smear_TrgObjComb[i]->Add(hpbpb_Smear_TrgObj80[i]);
    // hpbpb_Smear_TrgObjComb[i]->Add(hpbpb_Smear_TrgObj65[i]);
    // hpbpb_Smear_TrgObjComb[i]->Add(hpbpb_Smear_TrgObj55[i]);

    // divideBinWidth(hpbpb_Smear_TrgObjComb[i]);

    // hpbpb_JEC_TrgObj80[i]->Scale(1./ScaleFactor_80_NeqScale[i]);
    // hpbpb_JEC_TrgObj65[i]->Scale(1./ScaleFactor_65_NeqScale[i]);
    // hpbpb_JEC_TrgObj55[i]->Scale(1./ScaleFactor_55_NeqScale[i]);

    // hpbpb_JEC_TrgObjComb[i]->Add(hpbpb_JEC_TrgObj80[i]);
    // hpbpb_JEC_TrgObjComb[i]->Add(hpbpb_JEC_TrgObj65[i]);
    // hpbpb_JEC_TrgObjComb[i]->Add(hpbpb_JEC_TrgObj55[i]);

    // divideBinWidth(hpbpb_JEC_TrgObjComb[i]);

    // hpbpb_JEC_5pcent_TrgObj80[i]->Scale(1./ScaleFactor_80_NeqScale[i]);
    // hpbpb_JEC_5pcent_TrgObj65[i]->Scale(1./ScaleFactor_65_NeqScale[i]);
    // hpbpb_JEC_5pcent_TrgObj55[i]->Scale(1./ScaleFactor_55_NeqScale[i]);

    // hpbpb_JEC_5pcent_TrgObjComb[i]->Add(hpbpb_JEC_5pcent_TrgObj80[i]);
    // hpbpb_JEC_5pcent_TrgObjComb[i]->Add(hpbpb_JEC_5pcent_TrgObj65[i]);
    // hpbpb_JEC_5pcent_TrgObjComb[i]->Add(hpbpb_JEC_5pcent_TrgObj55[i]);

    // divideBinWidth(hpbpb_JEC_5pcent_TrgObjComb[i]);

    
  }
  */
  
  
  fout->Write();

  // myfile1.close();
  // myfile2.close();
  
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}//macro end
