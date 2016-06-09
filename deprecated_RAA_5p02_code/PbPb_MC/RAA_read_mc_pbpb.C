// Raghav Kunnawalkam Elayavalli
// Nov 24th 2015
// Rutgers @CERN
// for questions or comments: raghav.k.e at CERN dot CH

// 
// read all the MC files for PbPb and pp and make the required histograms for the analysis. 
// need to follow the same cuts used in the data analysis here as well. 
//

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

#define pi 3.14159265

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

float deltaphi(float phi1, float phi2)
{
  //float pi=TMath::Pi();
  
  float dphi = TMath::Abs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;

  return TMath::Abs(dphi);
}

// static const int nbins_cent = 6;
// static const Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
// //now we have to multiply by 5, since centrality goes from 0-200. 
// static const Double_t ncoll[nbins_cent] = { 1660, 1310, 745, 251, 62.8, 10.8 };
// static const int trigValue = 4;
// static const char trigName [trigValue][256] = {"HLT55","HLT65","HLT80","Combined"};
// static const char * etaWidth = (char*)"20_eta_20";

//const int ptbins_ana[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 500};
//const int nbins_ana = sizeof(ptbins_ana)/sizeof(int) -1;

const int ptbins_ana[] = {15, 30, 50, 80, 120, 170, 220, 280, 370, 500};
const int nbins_ana = sizeof(ptbins_ana)/sizeof(int) -1;

const double etabins_ana[] = {-3.000,
			     -2.500, -2.043, -1.740, -1.392,  
			     -1.131, -0.879, -0.609, -0.435, -0.261, -0.087,
			     +0.000, 
			     +0.087, +0.261, +0.435, +0.609,  
			     +0.879, +1.131, +1.392,  
			     +1.740, +2.043, +2.500,  
			     +3.000};
const int netabins_ana = sizeof(etabins_ana)/sizeof(double) -1;

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

// int findBin(int bin)
// {
//   int ibin=-1;
//   //! centrality is defined as 0.5% bins of cross section
//   //! in 0-200 bins               
//   if(bin<10)ibin=0; //! 0-5%
//   else if(bin>=10  && bin<20 )ibin=1; //! 5-10%
//   else if(bin>=20  && bin<60 )ibin=2;  //! 10-30%
//   else if(bin>=60  && bin<100)ibin=3;  //! 30-50%
//   else if(bin>=100 && bin<140)ibin=4;  //! 50-70%
//   else if(bin>=140 && bin<180)ibin=5;  //! 70-90%
//   else if(bin>=180 && bin<200)ibin=6;  //! 90-100%
//   return ibin;
// }


float deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float deta = eta1 - eta2;
  float dphi = fabs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;
  float dr = sqrt(pow(deta,2) + pow(dphi,2));
  return dr;
}

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

using namespace std;

void RAA_read_mc_pbpb(int startfile = 0,
		      int endfile = 9,
		      int radius = 4,
		      std::string algo="Vs",
		      std::string jetType="PF",
		      std::string kFoname="test.root"){
  
  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  bool printDebug = false;
  if(printDebug)cout<<"radius = "<<radius<<endl;
  
  TDatime date;

  std::string infile_Forest;

  //infile_Forest = "jetRAA_PbPb_mc_forests.txt";
  //infile_Forest = "pthat30_Vstest.txt";
  infile_Forest = "pthat80_Vstest.txt";
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;
  
  if(printDebug)cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
  }

  const int N = 4; //6

  TChain * jetpbpb[N];

  string dir[N];
  dir[0] = "hltanalysis";
  dir[1] = "skimanalysis";
  dir[2] = Form("ak%s%d%sJetAnalyzer", algo.c_str(),radius, jetType.c_str());
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
  
  // TFile * fweight = new TFile(Form("weights_pbpb_5p02TeV_PYHYD_%d.root",pthatBins2015[startfile]), "r");
  // TTree * weight_Tree = (TTree*)fweight->Get("weights");

  for(int ifile = startfile; ifile<endfile; ++ifile){

    instr_Forest>>filename_Forest;

    jetpbpb[0]->Add(filename_Forest.c_str());
    jetpbpb[1]->Add(filename_Forest.c_str());
    jetpbpb[2]->Add(filename_Forest.c_str());
    jetpbpb[3]->Add(filename_Forest.c_str());

    cout<<"filename: "<<filename_Forest<<endl;
    //cout<<"pthat = "<<pthat[startfile]<<endl;
    
    if(printDebug)cout << "Tree loaded  " << string(dir[0]+"/"+trees[0]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[0]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[1]+"/"+trees[1]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[1]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[2]+"/"+trees[2]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[2]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[3]+"/"+trees[3]).data() << endl;

    cout<<"Total number of events loaded in HiForest = "<<jetpbpb[2]->GetEntries()<<endl;

  }
  
  jetpbpb[2]->AddFriend(jetpbpb[0]);
  jetpbpb[2]->AddFriend(jetpbpb[1]);
  jetpbpb[2]->AddFriend(jetpbpb[3]);
  
  // Forest files 
  int nref_F;
  float pt_F[1000];
  float refpt_F[1000];
  float geneta_F[100];
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
  float chHardMax_F[1000];
  float trkHardMax_F[1000];
  float chHardSum_F[1000];
  float phHardSum_F[1000];
  float trkHardSum_F[1000];
  float phHardMax_F[1000];
  
  float jtpu_F[1000];
  int   subid_F[1000];
  float refdrjt_F[1000];
  int refparton_F[1000];
  float pthat_F;
  int jet55_F;
  int jet65_F;
  int jet80_F;
  int L1_sj36_F;
  int L1_sj52_F;
  int L1_sj36_p_F;
  int L1_sj52_p_F;
  int jet55_p_F;
  int jet65_p_F;
  int jet80_p_F;
  float vz_F;
  int evt_F;
  int run_F;
  int lumi_F;
  int hiNpix_F;
  int hiNpixelTracks_F;
  int hiBin_F;
  float hiHF_F;
  int hiNtracks_F;
  int hiNtracksPtCut_F;
  int hiNtracksEtaCut_F;
  int hiNtracksEtaPtCut_F;
  int pcollisionEventSelection_F;

  
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
  jetpbpb[3]->SetBranchAddress("vz",&vz_F);
  jetpbpb[1]->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection_F);
  //jetpbpb[0]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_F);
  //jetpbpb[0]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_F);
  //jetpbpb[0]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_F);
  jetpbpb[2]->SetBranchAddress("pthat",&pthat_F);
  jetpbpb[2]->SetBranchAddress("nref",&nref_F);
  jetpbpb[2]->SetBranchAddress("subid",subid_F);
  jetpbpb[2]->SetBranchAddress("refdrjt",refdrjt_F);
  jetpbpb[2]->SetBranchAddress("refparton_flavor",refparton_F);
  jetpbpb[2]->SetBranchAddress("refpt",refpt_F);
  jetpbpb[2]->SetBranchAddress("refeta",geneta_F);
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
  jetpbpb[0]->SetBranchAddress("HLT_HIJet55_v7",&jet55_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet55_v7_Prescl",&jet55_p_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet65_v7",&jet65_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet65_v7_Prescl",&jet65_p_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet80_v7",&jet80_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet80_v7_Prescl",&jet80_p_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet36_BptxAND",&L1_sj36_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&L1_sj36_p_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet52_BptxAND",&L1_sj52_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&L1_sj52_p_F);


  // double pthatweight;
  // int evt_w;
  // int lumi_w;
  // float vz_w;

  // weight_Tree->SetBranchAddress("pthatweight",&pthatweight);
  // weight_Tree->SetBranchAddress("evt",&evt_w);
  // weight_Tree->SetBranchAddress("lumi",&lumi_w);
  // weight_Tree->SetBranchAddress("vz",&vz_w);

  // jetpbpb[2]->AddFriend(weight_Tree);
  
  // get the power law fit weights for the gen pT spectra:
  // TFile * fpL = TFile::Open(Form("response_fit_weight_R%d.root",radius));
  // TH1F * hHistPrior_fit[nbins_cent];
 
  // for(int i = 0; i<nbins_cent; ++i)
  //   hHistPrior_fit[i] = (TH1F*)fpL->Get(Form("hHistPrior_fit_cent%d",i));


  std::string  etaWidth = "20_eta_20";
  
  TFile *fout = new TFile(kFoname.c_str(),"RECREATE");
  fout->cd();


  //get the spectra with the specific trigger object from the different files.
  // TH1F *hpbpb_Jet80_gen[nbins_cent],
  //   *hpbpb_Jet80_reco[nbins_cent],
  //   *hpbpb_Jet80_raw[nbins_cent],
  //   *hpbpb_Jet80_GenSmear[nbins_cent],
  //   *hpbpb_Jet80_RecoSmear[nbins_cent];
  // TH1F *hpbpb_Jet65_gen[nbins_cent],
  //   *hpbpb_Jet65_reco[nbins_cent],
  //   *hpbpb_Jet65_raw[nbins_cent],
  //   *hpbpb_Jet65_GenSmear[nbins_cent],
  //   *hpbpb_Jet65_RecoSmear[nbins_cent];
  // TH1F *hpbpb_Jet55_gen[nbins_cent],
  //   *hpbpb_Jet55_reco[nbins_cent],
  //   *hpbpb_Jet55_raw[nbins_cent],
  //   *hpbpb_Jet55_GenSmear[nbins_cent],
  //   *hpbpb_Jet55_RecoSmear[nbins_cent],
  //   * hpbpb_Jet55_reco_5TrigIneff_Smear[nbins_cent],
  //   * hpbpb_Jet55_reco_10TrigIneff_Smear[nbins_cent],
  //   * hpbpb_Jet55_gen_5TrigIneff_Smear[nbins_cent],
  //   * hpbpb_Jet55_gen_10TrigIneff_Smear[nbins_cent];
  // TH1F *hpbpb_JetComb_gen[nbins_cent],
  //   *hpbpb_JetComb_reco[nbins_cent],
  //   *hpbpb_JetComb_raw[nbins_cent],
  //   *hpbpb_JetComb_GenSmear[nbins_cent],
  //   *hpbpb_JetComb_RecoSmear[nbins_cent],
  //   * hpbpb_JetComb_reco_5TrigIneff_Smear[nbins_cent],
  //   * hpbpb_JetComb_reco_10TrigIneff_Smear[nbins_cent],
  //   * hpbpb_JetComb_gen_5TrigIneff_Smear[nbins_cent],
  //   * hpbpb_JetComb_gen_10TrigIneff_Smear[nbins_cent];

  // TH1F * hpbpb_Jet55Selection[nbins_cent];
  // TH1F * hpbpb_Jet65Selection[nbins_cent];
  
  // TH1F * hpbpb_JetComb_gen2pSmear[nbins_cent],
  //   * hpbpb_Jet80_gen2pSmear[nbins_cent],
  //   * hpbpb_Jet65_gen2pSmear[nbins_cent],
  //   * hpbpb_Jet55_gen2pSmear[nbins_cent];

  TH1F * hpbpb_pthat[nbins_cent];
  
  TH1F *hpbpb_gen[nbins_cent],*hpbpb_reco[nbins_cent];
  TH2F *hpbpb_matrix[nbins_cent];
  
  // TH2F * hpbpb_matrix_HLT[nbins_cent];
  // TH2F * hpbpb_matrix_HLT_Jet80[nbins_cent];
  // TH2F * hpbpb_matrix_HLT_Jet65[nbins_cent];
  // TH2F * hpbpb_matrix_HLT_Jet55[nbins_cent];
  // TH2F * hpbpb_matrix_recoTrunc_HLT[nbins_cent];
  
  // TH2F * hpbpb_Trans_matrix_HLT[nbins_cent];
  // TH2F * hpbpb_Trans_matrix_Smooth_HLT[nbins_cent];
  // TH2F * hpbpb_Trans_matrix_HLT_Jet80[nbins_cent];
  // TH2F * hpbpb_Trans_matrix_HLT_Jet65[nbins_cent];
  // TH2F * hpbpb_Trans_matrix_HLT_Jet55[nbins_cent];

  // TH2F * hpbpb_anaBin_matrix_HLT[nbins_cent];
  // TH2F * hpbpb_anaBin_Trans_matrix_HLT[nbins_cent];
  // TH2F * hpbpb_anaBin_matrix_HLT_Jet80[nbins_cent];
  // TH2F * hpbpb_anaBin_Trans_matrix_HLT_Jet80[nbins_cent];
  // TH2F * hpbpb_anaBin_matrix_HLT_Jet65[nbins_cent];
  // TH2F * hpbpb_anaBin_Trans_matrix_HLT_Jet65[nbins_cent];
  // TH2F * hpbpb_anaBin_matrix_HLT_Jet55[nbins_cent];
  // TH2F * hpbpb_anaBin_Trans_matrix_HLT_Jet55[nbins_cent];
  
  // TH2F * hpbpb_matrix_HLT_GenSmear[nbins_cent];
  // TH2F * hpbpb_matrix_HLT_gen2pSmear[nbins_cent];
  // TH2F * hpbpb_matrix_HLT_RecoSmear[nbins_cent];
  // TH2F * hpbpb_matrix_HLT_BothSmear[nbins_cent];
  // TH2F * hpbpb_matrix_HLT_5TrigIneff_Smear[nbins_cent];
  // TH2F * hpbpb_matrix_HLT_10TrigIneff_Smear[nbins_cent];

  // TH2F * hpbpb_Trans_matrix_HLT_GenSmear[nbins_cent];
  // TH2F * hpbpb_Trans_matrix_HLT_gen2pSmear[nbins_cent];
  // TH2F * hpbpb_Trans_matrix_HLT_RecoSmear[nbins_cent];
  // TH2F * hpbpb_Trans_matrix_HLT_BothSmear[nbins_cent];
  
  TH2F * hpbpb_mcclosure_matrix[nbins_cent];
  // TH2F * hpbpb_mcclosure_matrix_HLT[nbins_cent];
  // TH2F * hpbpb_10GeVBins_mcclosure_matrix_HLT[nbins_cent];
  
  // TH2F * hpbpb_anaBin_mcclosure_matrix_HLT[nbins_cent];
  // TH2F * hpbpb_mcclosure_Trans_matrix_HLT[nbins_cent];
  // TH2F * hpbpb_anaBin_mcclosure_Trans_matrix_HLT[nbins_cent];

  // TH2F * hpbpb_mcclosure_matrix_HLT_Jet80[nbins_cent];
  // TH2F * hpbpb_anaBin_mcclosure_matrix_HLT_Jet80[nbins_cent];
  // TH2F * hpbpb_mcclosure_Trans_matrix_HLT_Jet80[nbins_cent];
  // TH2F * hpbpb_anaBin_mcclosure_Trans_matrix_HLT_Jet80[nbins_cent];

  // TH2F * hpbpb_mcclosure_matrix_HLT_Jet65[nbins_cent];
  // TH2F * hpbpb_anaBin_mcclosure_matrix_HLT_Jet65[nbins_cent];
  // TH2F * hpbpb_mcclosure_Trans_matrix_HLT_Jet65[nbins_cent];
  // TH2F * hpbpb_anaBin_mcclosure_Trans_matrix_HLT_Jet65[nbins_cent];

  // TH2F * hpbpb_mcclosure_matrix_HLT_Jet55[nbins_cent];
  // TH2F * hpbpb_anaBin_mcclosure_matrix_HLT_Jet55[nbins_cent];
  // TH2F * hpbpb_mcclosure_Trans_matrix_HLT_Jet55[nbins_cent];
  // TH2F * hpbpb_anaBin_mcclosure_Trans_matrix_HLT_Jet55[nbins_cent];


  // TH2F * hpbpb_mcclosure_matrix_HLT_5TrigIneff_Smear[nbins_cent];  
  // TH2F * hpbpb_mcclosure_matrix_HLT_10TrigIneff_Smear[nbins_cent];  
  TH1F * hpbpb_mcclosure_data[nbins_cent];
  TH1F * hpbpb_mcclosure_data_train[nbins_cent];
  // TH1F * hpbpb_mcclosure_JetComb_data[nbins_cent];
  // TH1F * hpbpb_mcclosure_JetComb_data_5TrigIneff_Smear[nbins_cent];
  // TH1F * hpbpb_mcclosure_JetComb_data_10TrigIneff_Smear[nbins_cent];
  // TH1F * hpbpb_mcclosure_Jet80_data[nbins_cent];
  // TH1F * hpbpb_mcclosure_Jet65_data[nbins_cent];
  // TH1F * hpbpb_mcclosure_Jet55_data[nbins_cent];
  // TH1F * hpbpb_mcclosure_Jet55_data_5TrigIneff_Smear[nbins_cent];
  // TH1F * hpbpb_mcclosure_Jet55_data_10TrigIneff_Smear[nbins_cent];
  // TH1F * hpbpb_mcclosure_JetComb_data_train[nbins_cent];
  // TH1F * hpbpb_mcclosure_JetComb_data_train_5TrigIneff_Smear[nbins_cent];
  // TH1F * hpbpb_mcclosure_JetComb_data_train_10TrigIneff_Smear[nbins_cent];
  // TH1F * hpbpb_mcclosure_Jet80_data_train[nbins_cent];
  // TH1F * hpbpb_mcclosure_Jet65_data_train[nbins_cent];
  // TH1F * hpbpb_mcclosure_Jet55_data_train[nbins_cent];
  // TH1F * hpbpb_mcclosure_Jet55_data_train_5TrigIneff_Smear[nbins_cent];
  // TH1F * hpbpb_mcclosure_Jet55_data_train_10TrigIneff_Smear[nbins_cent];
  TH1F * hpbpb_mcclosure_gen[nbins_cent];
  // TH1F * hpbpb_mcclosure_JetComb_gen[nbins_cent];
  // TH1F * hpbpb_mcclosure_JetComb_gen_5TrigIneff_Smear[nbins_cent];
  // TH1F * hpbpb_mcclosure_JetComb_gen_10TrigIneff_Smear[nbins_cent];
  // TH1F * hpbpb_mcclosure_Jet80_gen[nbins_cent];
  // TH1F * hpbpb_mcclosure_Jet65_gen[nbins_cent];
  // TH1F * hpbpb_mcclosure_Jet55_gen[nbins_cent];
  // TH1F * hpbpb_mcclosure_Jet55_gen_5TrigIneff_Smear[nbins_cent];
  // TH1F * hpbpb_mcclosure_Jet55_gen_10TrigIneff_Smear[nbins_cent];


  // TH1F * hpbpb_anaBin_Jet80_gen[nbins_cent],* hpbpb_anaBin_Jet80_reco[nbins_cent];
  // TH1F * hpbpb_anaBin_Jet65_gen[nbins_cent],* hpbpb_anaBin_Jet65_reco[nbins_cent];
  // TH1F * hpbpb_anaBin_Jet55_gen[nbins_cent],* hpbpb_anaBin_Jet55_reco[nbins_cent];
  // TH1F * hpbpb_anaBin_JetComb_gen[nbins_cent],* hpbpb_anaBin_JetComb_reco[nbins_cent];

  // TH2F * hSVD_Matrix_noWeight[nbins_cent];
  // TH2F * hSVD_Matrix_Weight[nbins_cent];
  // TH2F * hSVD_Trans_Matrix_noWeight[nbins_cent];
  // TH2F * hSVD_Trans_Matrix_Weight[nbins_cent];

  // Int_t nbins_jtpt[nbins_cent]       = {60, 60, 60, 60, 40, 24};
  // Double_t binLow_jtpt[nbins_cent]   = {50.0, 50.0, 50.0, 50.0, 50.0, 50.0};
  // Double_t binHigh_jtpt[nbins_cent]  = {350.0, 350.0, 350.0, 350.0, 250.0, 170.0};
  // // Int_t nbins_genpt      = 20;
  // // Double_t binLow_genpt  = 0.0;
  // // Double_t binHigh_genpt = 400.0;
  
  // Int_t nbins_genpt[nbins_cent]      = {40, 40, 40, 40, 30, 20};
  // Double_t binLow_genpt[nbins_cent]  = {20.0, 20.0, 20.0, 20.0, 20.0, 20.0};
  // Double_t binHigh_genpt[nbins_cent] = {420.0, 420.0, 420.0, 420.0, 320.0, 220.0};

  // TH2F * heMaxSumcand_vs_calopfpt_RecopTSelection[ptSelection][nbins_cent+1];
  // TH2F * heMaxSumcand_vs_calopfpt_RecoGenpTSelection[ptSelection][nbins_cent+1];

  TH1F * hJER[nbins_ana][nbins_cent];
  TH1F * hJER_eta_30pt50[netabins_ana][nbins_cent];
  TH1F * hJER_eta_80pt120[netabins_ana][nbins_cent];
  TH3F * hJEC[nbins_cent];

  TH1F * hdeltaphi[nbins_cent];

  std::string var[19] = {"jtpt" ,"rawpt", "jteta", "jtphi", "trkMax", "trkSum", "trkHardSum", "chMax", "chSum", "chHardSum","phMax", "phSum", "phHardSum", "neMax", "neSum", "eMax", "eSum", "muMax", "muSum" };
  TH1F * hJetQA[2][19][nbins_cent];
  
  // for(int i = 0; i<=nbins_cent; ++i){
  //   for(int j = 0; j < ptSelection; ++j){
  //     heMaxSumcand_vs_calopfpt_RecopTSelection[j][i] = new TH2F(Form("heMaxSumcand_vs_calopfpt_%d_RecopTSelection_%d_cent%d", ptBoundary[j], ptBoundary[j+1], i),"",100, 0, 2.5, 100, 0, 5);  
  //     heMaxSumcand_vs_calopfpt_RecoGenpTSelection[j][i] = new TH2F(Form("heMaxSumcand_vs_calopfpt_%d_RecoGenpTSelection_%d_cent%d", ptBoundary[j], ptBoundary[j+1], i),"",100, 0, 2.5, 100, 0, 5);              
  //   }
  // }

  // macros for quark vs gluon uncertainty from the MC.
  // TH1F * hpbpb_quark_Jet80_gen[nbins_cent];
  // TH1F * hpbpb_quark_Jet65_gen[nbins_cent];
  // TH1F * hpbpb_quark_Jet55_gen[nbins_cent];
  // TH1F * hpbpb_quark_JetComb_gen[nbins_cent];
  // TH1F * hpbpb_quark_Jet80_reco[nbins_cent];
  // TH1F * hpbpb_quark_Jet65_reco[nbins_cent];
  // TH1F * hpbpb_quark_Jet55_reco[nbins_cent];
  // TH1F * hpbpb_quark_JetComb_reco[nbins_cent];

  // TH1F * hpbpb_gluon_Jet80_gen[nbins_cent];
  // TH1F * hpbpb_gluon_Jet65_gen[nbins_cent];
  // TH1F * hpbpb_gluon_Jet55_gen[nbins_cent];
  // TH1F * hpbpb_gluon_JetComb_gen[nbins_cent];
  // TH1F * hpbpb_gluon_Jet80_reco[nbins_cent];
  // TH1F * hpbpb_gluon_Jet65_reco[nbins_cent];
  // TH1F * hpbpb_gluon_Jet55_reco[nbins_cent];
  // TH1F * hpbpb_gluon_JetComb_reco[nbins_cent];

  // TH1F * hpbpb_mcclosure_quark_Jet80_gen[nbins_cent];
  // TH1F * hpbpb_mcclosure_quark_Jet65_gen[nbins_cent];
  // TH1F * hpbpb_mcclosure_quark_Jet55_gen[nbins_cent];
  // TH1F * hpbpb_mcclosure_quark_JetComb_gen[nbins_cent];
  // TH1F * hpbpb_mcclosure_quark_Jet80_reco[nbins_cent];
  // TH1F * hpbpb_mcclosure_quark_Jet65_reco[nbins_cent];
  // TH1F * hpbpb_mcclosure_quark_Jet55_reco[nbins_cent];
  // TH1F * hpbpb_mcclosure_quark_JetComb_reco[nbins_cent];

  // TH1F * hpbpb_mcclosure_gluon_Jet80_gen[nbins_cent];
  // TH1F * hpbpb_mcclosure_gluon_Jet65_gen[nbins_cent];
  // TH1F * hpbpb_mcclosure_gluon_Jet55_gen[nbins_cent];
  // TH1F * hpbpb_mcclosure_gluon_JetComb_gen[nbins_cent];
  // TH1F * hpbpb_mcclosure_gluon_Jet80_reco[nbins_cent];
  // TH1F * hpbpb_mcclosure_gluon_Jet65_reco[nbins_cent];
  // TH1F * hpbpb_mcclosure_gluon_Jet55_reco[nbins_cent];
  // TH1F * hpbpb_mcclosure_gluon_JetComb_reco[nbins_cent];
  
  // TH2F * hpbpb_Trans_matrix_quark_HLT[nbins_cent];
  // TH2F * hpbpb_Trans_matrix_gluon_HLT[nbins_cent];

  // TH2F * hpbpb_mcclosure_Trans_matrix_quark_HLT[nbins_cent];
  // TH2F * hpbpb_mcclosure_Trans_matrix_gluon_HLT[nbins_cent];  
  
  for(int i = 0;i<nbins_cent;++i){

    for(int k = 0; k<2; ++k){
      for(int j = 0; j<19; ++j){
	if(j==2) hJetQA[k][j][i] = new TH1F(Form("hJetQA_%dwJetID_%s_cent%d",k,var[j].c_str(),i),Form(";%s;",var[j].c_str()),100, -5, +5);
	else if(j==3) hJetQA[k][j][i] = new TH1F(Form("hJetQA_%dwJetID_%s_cent%d",k,var[j].c_str(),i),Form(";%s;",var[j].c_str()),100, -4, +4);
	else if(j<=1)hJetQA[k][j][i] = new TH1F(Form("hJetQA_%dwJetID_%s_cent%d",k,var[j].c_str(),i),Form(";%s;",var[j].c_str()),500, 0, 500);
	else if(j>=4)hJetQA[k][j][i] = new TH1F(Form("hJetQA_%dwJetID_%s_cent%d",k,var[j].c_str(),i),Form(";%s;",var[j].c_str()),200, 0, 2);
      }
    }
    
    hdeltaphi[i] = new TH1F(Form("hdeltaphi_cent%d",i),";#delta #phi_{jt1, jt2};counts",100, 0, 3.5);
    
    hJEC[i] = new TH3F(Form("hJEC_cent%d",i),";raw p_{T};#eta;JEC",500, 0, 500, 200, -5, +5, 300, 0, 5);

    // hpbpb_quark_Jet80_gen[i] = new TH1F(Form("hpbpb_quark_Jet80_gen_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_quark_Jet65_gen[i] = new TH1F(Form("hpbpb_quark_Jet65_gen_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_quark_Jet55_gen[i] = new TH1F(Form("hpbpb_quark_Jet55_gen_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_quark_JetComb_gen[i] = new TH1F(Form("hpbpb_quark_JetComb_gen_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_quark_Jet80_reco[i] = new TH1F(Form("hpbpb_quark_Jet80_reco_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_quark_Jet65_reco[i] = new TH1F(Form("hpbpb_quark_Jet65_reco_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_quark_Jet55_reco[i] = new TH1F(Form("hpbpb_quark_Jet55_reco_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_quark_JetComb_reco[i] = new TH1F(Form("hpbpb_quark_JetComb_reco_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_gluon_Jet80_gen[i] = new TH1F(Form("hpbpb_gluon_Jet80_gen_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_gluon_Jet65_gen[i] = new TH1F(Form("hpbpb_gluon_Jet65_gen_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_gluon_Jet55_gen[i] = new TH1F(Form("hpbpb_gluon_Jet55_gen_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_gluon_JetComb_gen[i] = new TH1F(Form("hpbpb_gluon_JetComb_gen_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_gluon_Jet80_reco[i] = new TH1F(Form("hpbpb_gluon_Jet80_reco_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_gluon_Jet65_reco[i] = new TH1F(Form("hpbpb_gluon_Jet65_reco_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_gluon_Jet55_reco[i] = new TH1F(Form("hpbpb_gluon_Jet55_reco_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_gluon_JetComb_reco[i] = new TH1F(Form("hpbpb_gluon_JetComb_reco_R%d_cent%d", radius, i),"",1000, 0, 1000);

    // hpbpb_mcclosure_quark_Jet80_gen[i] = new TH1F(Form("hpbpb_mcclosure_quark_Jet80_gen_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_mcclosure_quark_Jet65_gen[i] = new TH1F(Form("hpbpb_mcclosure_quark_Jet65_gen_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_mcclosure_quark_Jet55_gen[i] = new TH1F(Form("hpbpb_mcclosure_quark_Jet55_gen_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_mcclosure_quark_JetComb_gen[i] = new TH1F(Form("hpbpb_mcclosure_quark_JetComb_gen_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_mcclosure_quark_Jet80_reco[i] = new TH1F(Form("hpbpb_mcclosure_quark_Jet80_reco_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_mcclosure_quark_Jet65_reco[i] = new TH1F(Form("hpbpb_mcclosure_quark_Jet65_reco_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_mcclosure_quark_Jet55_reco[i] = new TH1F(Form("hpbpb_mcclosure_quark_Jet55_reco_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_mcclosure_quark_JetComb_reco[i] = new TH1F(Form("hpbpb_mcclosure_quark_JetComb_reco_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_mcclosure_gluon_Jet80_gen[i] = new TH1F(Form("hpbpb_mcclosure_gluon_Jet80_gen_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_mcclosure_gluon_Jet65_gen[i] = new TH1F(Form("hpbpb_mcclosure_gluon_Jet65_gen_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_mcclosure_gluon_Jet55_gen[i] = new TH1F(Form("hpbpb_mcclosure_gluon_Jet55_gen_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_mcclosure_gluon_JetComb_gen[i] = new TH1F(Form("hpbpb_mcclosure_gluon_JetComb_gen_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_mcclosure_gluon_Jet80_reco[i] = new TH1F(Form("hpbpb_mcclosure_gluon_Jet80_reco_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_mcclosure_gluon_Jet65_reco[i] = new TH1F(Form("hpbpb_mcclosure_gluon_Jet65_reco_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_mcclosure_gluon_Jet55_reco[i] = new TH1F(Form("hpbpb_mcclosure_gluon_Jet55_reco_R%d_cent%d", radius, i),"",1000, 0, 1000);
    // hpbpb_mcclosure_gluon_JetComb_reco[i] = new TH1F(Form("hpbpb_mcclosure_gluon_JetComb_reco_R%d_cent%d", radius, i),"",1000, 0, 1000);

    // hpbpb_Trans_matrix_quark_HLT[i] = new TH2F(Form("hpbpb_Trans_matrix_quark_HLT_R%d_cent%d",radius,i),"",1000, 0, 1000, 1000, 0, 1000);
    // hpbpb_Trans_matrix_gluon_HLT[i] = new TH2F(Form("hpbpb_Trans_matrix_gluon_HLT_R%d_cent%d",radius,i),"",1000, 0, 1000, 1000, 0, 1000);
    // hpbpb_mcclosure_Trans_matrix_quark_HLT[i] = new TH2F(Form("hpbpb_mcclosure_Trans_matrix_quark_HLT_R%d_cent%d",radius,i),"",1000, 0, 1000, 1000, 0, 1000);
    // hpbpb_mcclosure_Trans_matrix_gluon_HLT[i] = new TH2F(Form("hpbpb_mcclosure_Trans_matrix_gluon_HLT_R%d_cent%d",radius,i),"",1000, 0, 1000, 1000, 0, 1000);

    for(int bin = 0; bin<netabins_ana; ++bin){
      hJER_eta_30pt50[bin][i] = new TH1F(Form("hJER_etabin_%d_30pt50_cent%d", bin, i),"",150, 0, 3);
      hJER_eta_80pt120[bin][i] = new TH1F(Form("hJER_etabin_%d_80pt120_cent%d", bin, i),"",150, 0, 3);
    }
    
    for(int bin = 0; bin<nbins_ana; ++bin){
      hJER[bin][i] = new TH1F(Form("hJER_%d_pt_%d_cent%d", ptbins_ana[bin], ptbins_ana[bin+1], i),"",150, 0, 3);
    }
    // hpbpb_Trans_matrix_Smooth_HLT[i] = new TH2F(Form("hpbpb_Trans_matrix_Smooth_HLT_R%d_%s_cent%d",radius,etaWidth.c_str(),i),
    // 						Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),
    // 						1000, 0, 1000,
    // 						1000, 0, 1000);
    // hpbpb_Trans_matrix_HLT_GenSmear[i] = new TH2F(Form("hpbpb_Trans_matrix_HLT_GenSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),
    // 						  Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),
    // 						  1000, 0, 1000,
    // 						  1000, 0, 1000);
    // hpbpb_Trans_matrix_HLT_gen2pSmear[i] = new TH2F(Form("hpbpb_Trans_matrix_HLT_gen2pSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),
    // 						    Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),
    // 						    1000, 0, 1000,
    // 						    1000, 0, 1000);
    // hpbpb_Trans_matrix_HLT_RecoSmear[i] = new TH2F(Form("hpbpb_Trans_matrix_HLT_RecoSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),
    // 						   Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),
    // 						   1000, 0, 1000,
    // 						   1000, 0, 1000);
    // hpbpb_Trans_matrix_HLT_BothSmear[i] = new TH2F(Form("hpbpb_Trans_matrix_HLT_BothSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),
    // 						   Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),
    // 						   1000, 0, 1000,
    // 						   1000, 0, 1000);
    // hpbpb_Trans_matrix_HLT[i] = new TH2F(Form("hpbpb_Trans_matrix_HLT_R%d_%s_cent%d",radius,etaWidth.c_str(),i),
    // 					 Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),
    // 					 1000, 0, 1000,
    // 					 1000, 0, 1000);
    // hpbpb_mcclosure_Trans_matrix_HLT[i] = new TH2F(Form("hpbpb_mcclosure_Trans_matrix_HLT_R%d_%s_cent%d",radius,etaWidth.c_str(),i),
    // 						   Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),
    // 						   1000, 0, 1000,
    // 						   1000, 0, 1000);
    // hpbpb_mcclosure_Trans_matrix_HLT_Jet55[i] = new TH2F(Form("hpbpb_mcclosure_Trans_matrix_HLT_Jet55_R%d_%s_cent%d",radius,etaWidth.c_str(),i),
    // 							 Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),
    // 							 1000, 0, 1000,
    // 							 1000, 0, 1000);
    // hpbpb_mcclosure_Trans_matrix_HLT_Jet80[i] = new TH2F(Form("hpbpb_mcclosure_Trans_matrix_HLT_Jet80_R%d_%s_cent%d",radius,etaWidth.c_str(),i),
    // 							 Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),
    // 							 1000, 0, 1000,
    // 							 1000, 0, 1000);
    // hpbpb_mcclosure_Trans_matrix_HLT_Jet65[i] = new TH2F(Form("hpbpb_mcclosure_Trans_matrix_HLT_Jet65_R%d_%s_cent%d",radius,etaWidth.c_str(),i),
    // 							 Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),
    // 							 1000, 0, 1000,
    // 							 1000, 0, 1000);

    // hpbpb_anaBin_mcclosure_Trans_matrix_HLT[i] = new TH2F(Form("hpbpb_anaBin_mcclosure_Trans_matrix_HLT_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000, 1000, 0, 1000);
    // hpbpb_anaBin_mcclosure_Trans_matrix_HLT_Jet80[i] = new TH2F(Form("hpbpb_anaBin_mcclosure_Trans_matrix_HLT_Jet80_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000, 1000, 0, 1000);
    // hpbpb_anaBin_mcclosure_Trans_matrix_HLT_Jet65[i] = new TH2F(Form("hpbpb_anaBin_mcclosure_Trans_matrix_HLT_Jet65_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000, 1000, 0, 1000);
    // hpbpb_anaBin_mcclosure_Trans_matrix_HLT_Jet55[i] = new TH2F(Form("hpbpb_anaBin_mcclosure_Trans_matrix_HLT_Jet55_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000, 1000, 0, 1000);
    // hpbpb_mcclosure_matrix_HLT_5TrigIneff_Smear[i] = new TH2F(Form("hpbpb_mcclosure_matrix_HLT_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000,1000, 0, 1000);
    // hpbpb_mcclosure_matrix_HLT_10TrigIneff_Smear[i] = new TH2F(Form("hpbpb_mcclosure_matrix_HLT_10TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000,1000, 0, 1000);
    
    
    hpbpb_gen[i] = new TH1F(Form("hpbpb_gen_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Gen refpt R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_JetComb_gen[i] = new TH1F(Form("hpbpb_JetComb_gen_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Gen refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet80_gen[i] = new TH1F(Form("hpbpb_Jet80_gen_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Gen refpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet65_gen[i] = new TH1F(Form("hpbpb_Jet65_gen_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Gen refpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet55_gen[i] = new TH1F(Form("hpbpb_Jet55_gen_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Gen refpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_JetComb_GenSmear[i] = new TH1F(Form("hpbpb_JetComb_GenSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_JetComb_gen2pSmear[i] = new TH1F(Form("hpbpb_JetComb_gen2pSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("gen2pSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet80_GenSmear[i] = new TH1F(Form("hpbpb_Jet80_GenSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("GenSmear refpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet80_gen2pSmear[i] = new TH1F(Form("hpbpb_Jet80_gen2pSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("gen2pSmear refpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet65_GenSmear[i] = new TH1F(Form("hpbpb_Jet65_GenSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("GenSmear refpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet65_gen2pSmear[i] = new TH1F(Form("hpbpb_Jet65_gen2pSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("gen2pSmear refpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet55_GenSmear[i] = new TH1F(Form("hpbpb_Jet55_GenSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("GenSmear refpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet55_gen2pSmear[i] = new TH1F(Form("hpbpb_Jet55_gen2pSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("gen2pSmear refpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_JetComb_gen_5TrigIneff_Smear[i] = new TH1F(Form("hpbpb_JetComb_gen_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet55_gen_5TrigIneff_Smear[i] = new TH1F(Form("hpbpb_Jet55_gen_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_JetComb_gen_10TrigIneff_Smear[i] = new TH1F(Form("hpbpb_JetComb_gen_10TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet55_gen_10TrigIneff_Smear[i] = new TH1F(Form("hpbpb_Jet55_gen_10TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_mcclosure_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("gen spectra for unfolding mc closure test R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_JetComb_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_JetComb_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("gen spectra for unfolding mc closure test trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_JetComb_gen_5TrigIneff_Smear[i] = new TH1F(Form("hpbpb_mcclosure_gen_JetComb_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("gen spectra for unfolding mc closure test trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_JetComb_gen_10TrigIneff_Smear[i] = new TH1F(Form("hpbpb_mcclosure_gen_JetComb_10TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("gen spectra for unfolding mc closure test trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_Jet80_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet80_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("gen spectra for unfolding mc closure test trigger 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_Jet65_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet65_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("gen spectra for unfolding mc closure test trigger 65 R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_Jet55_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet55_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("gen spectra for unfolding mc closure test trigger 55 R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_Jet55_gen_5TrigIneff_Smear[i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet55_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("gen spectra for unfolding mc closure test trigger 55 R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_Jet55_gen_10TrigIneff_Smear[i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet55_10TrigIneff_Smear[_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("gen spectra for unfolding mc closure test trigger 55 R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    
    hpbpb_reco[i] = new TH1F(Form("hpbpb_reco_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Reco jtpt R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_JetComb_reco[i] = new TH1F(Form("hpbpb_JetComb_reco_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("reco jtpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet80_reco[i] = new TH1F(Form("hpbpb_Jet80_reco_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("reco jtpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet65_reco[i] = new TH1F(Form("hpbpb_Jet65_reco_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("reco jtpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet55_reco[i] = new TH1F(Form("hpbpb_Jet55_reco_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("reco jtpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_JetComb_RecoSmear[i] = new TH1F(Form("hpbpb_JetComb_RecoSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("RecoSmear jtpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet80_RecoSmear[i] = new TH1F(Form("hpbpb_Jet80_RecoSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("RecoSmear jtpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet65_RecoSmear[i] = new TH1F(Form("hpbpb_Jet65_RecoSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("RecoSmear jtpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet55_RecoSmear[i] = new TH1F(Form("hpbpb_Jet55_RecoSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("RecoSmear jtpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_JetComb_reco_5TrigIneff_Smear[i] = new TH1F(Form("hpbpb_JetComb_reco_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet55_reco_5TrigIneff_Smear[i] = new TH1F(Form("hpbpb_Jet55_reco_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_JetComb_reco_10TrigIneff_Smear[i] = new TH1F(Form("hpbpb_JetComb_reco_10TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_Jet55_reco_10TrigIneff_Smear[i] = new TH1F(Form("hpbpb_Jet55_reco_10TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_mcclosure_data[i] =new TH1F(Form("hpbpb_mcclosure_data_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("data for unfolding mc closure test R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    hpbpb_mcclosure_data_train[i] =new TH1F(Form("hpbpb_mcclosure_data_train_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("data_train for unfolding mc closure test R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_JetComb_data[i] = new TH1F(Form("hpbpb_mcclosure_JetComb_data_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("data for unfolding mc closure test trigger combined  R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_JetComb_data_5TrigIneff_Smear[i] = new TH1F(Form("hpbpb_mcclosure_JetComb_data_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("data for unfolding mc closure test trigger combined  R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_JetComb_data_10TrigIneff_Smear[i] = new TH1F(Form("hpbpb_mcclosure_JetComb_data_10TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("data for unfolding mc closure test trigger combined  R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_Jet80_data[i] = new TH1F(Form("hpbpb_mcclosure_Jet80_data_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("data for unfolding mc closure test trigger 80  R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_Jet65_data[i] = new TH1F(Form("hpbpb_mcclosure_Jet65_data_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("data for unfolding mc closure test trigger 65  R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_Jet55_data[i] = new TH1F(Form("hpbpb_mcclosure_Jet55_data_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("data for unfolding mc closure test trigger 55  R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_Jet55_data_5TrigIneff_Smear[i] = new TH1F(Form("hpbpb_mcclosure_Jet55_data_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("data for unfolding mc closure test trigger 55  R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_Jet55_data_10TrigIneff_Smear[i] = new TH1F(Form("hpbpb_mcclosure_Jet55_data_10TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("data for unfolding mc closure test trigger 55  R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_JetComb_data_train[i] = new TH1F(Form("hpbpb_mcclosure_JetComb_data_train_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("data_train for unfolding mc closure test trigger combined  R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_Jet80_data_train[i] = new TH1F(Form("hpbpb_mcclosure_Jet80_data_train_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("data_train for unfolding mc closure test trigger 80  R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_Jet65_data_train[i] = new TH1F(Form("hpbpb_mcclosure_Jet65_data_train_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("data_train for unfolding mc closure test trigger 65  R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);
    // hpbpb_mcclosure_Jet55_data_train[i] = new TH1F(Form("hpbpb_mcclosure_Jet55_data_train_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("data_train for unfolding mc closure test trigger 55  R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000);

    
    // hSVD_Matrix_noWeight[i] = new TH2F(Form("hSVD_Matrix_noWeight_cent%d", i),"",30, 0, 300, 30, 0, 350);
    // hSVD_Matrix_Weight[i] = new TH2F(Form("hSVD_Matrix_Weight_cent%d", i),"",30, 0, 300, 25, 50, 300);
    // hSVD_Trans_Matrix_noWeight[i] = new TH2F(Form("hSVD_Trans_Matrix_noWeight_cent%d", i),"",25, 50, 300, 30, 0, 300);
    // hSVD_Trans_Matrix_Weight[i] = new TH2F(Form("hSVD_Trans_Matrix_Weight_cent%d", i),"",25, 50, 300, 30, 0, 300);

    // hpbpb_Jet55Selection[i] = new TH1F(Form("hpbpb_Jet55Selection_cent%d",i),"",1000,0,1000);
    // hpbpb_Jet65Selection[i] = new TH1F(Form("hpbpb_Jet65Selection_cent%d",i),"",1000,0,1000);

    hpbpb_pthat[i] = new TH1F(Form("pthat_cent%d",i),"",1000,0,1000);
    // hpbpb_pthat_newWeight[i] = new TH1F(Form("pthat_newWeight_cent%d",i),"",1000,0,1000);
    

    hpbpb_matrix[i] = new TH2F(Form("hpbpb_matrix_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    
    // hpbpb_matrix_HLT[i] = new TH2F(Form("hpbpb_matrix_HLT_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    // hpbpb_matrix_recoTrunc_HLT[i] = new TH2F(Form("hpbpb_matrix_recoTrunc_HLT_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,300,50,350);      
    // hpbpb_anaBin_matrix_HLT[i] = new TH2F(Form("hpbpb_anaBin_matrix_HLT_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
    // hpbpb_anaBin_Trans_matrix_HLT[i] = new TH2F(Form("hpbpb_anaBin_Trans_matrix_HLT_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);

    // hpbpb_matrix_HLT_Jet80[i] = new TH2F(Form("hpbpb_matrix_HLT_Jet80_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    // hpbpb_anaBin_matrix_HLT_Jet80[i] = new TH2F(Form("hpbpb_anaBin_matrix_HLT_Jet80_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
    // hpbpb_anaBin_Trans_matrix_HLT_Jet80[i] = new TH2F(Form("hpbpb_anaBin_Trans_matrix_HLT_Jet80_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);

    // hpbpb_Trans_matrix_HLT_Jet80[i] = new TH2F(Form("hpbpb_Trans_matrix_HLT_Jet80_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    // hpbpb_Trans_matrix_HLT_Jet65[i] = new TH2F(Form("hpbpb_Trans_matrix_HLT_Jet65_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    // hpbpb_Trans_matrix_HLT_Jet55[i] = new TH2F(Form("hpbpb_Trans_matrix_HLT_Jet55_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);

    // hpbpb_matrix_HLT_Jet65[i] = new TH2F(Form("hpbpb_matrix_HLT_Jet65_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    // hpbpb_anaBin_matrix_HLT_Jet65[i] = new TH2F(Form("hpbpb_anaBin_matrix_HLT_Jet65_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
    // hpbpb_anaBin_Trans_matrix_HLT_Jet65[i] = new TH2F(Form("hpbpb_anaBin_Trans_matrix_HLT_Jet65_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);

    // hpbpb_matrix_HLT_Jet55[i] = new TH2F(Form("hpbpb_matrix_HLT_Jet55_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    // hpbpb_anaBin_matrix_HLT_Jet55[i] = new TH2F(Form("hpbpb_anaBin_matrix_HLT_Jet55_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
    // hpbpb_anaBin_Trans_matrix_HLT_Jet55[i] = new TH2F(Form("hpbpb_anaBin_Trans_matrix_HLT_Jet55_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);

    
    
    // hpbpb_matrix_HLT_GenSmear[i] = new TH2F(Form("hpbpb_matrix_HLT_GenSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    // hpbpb_matrix_HLT_gen2pSmear[i] = new TH2F(Form("hpbpb_matrix_HLT_gen2pSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    // hpbpb_matrix_HLT_RecoSmear[i] = new TH2F(Form("hpbpb_matrix_HLT_RecoSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    // hpbpb_matrix_HLT_BothSmear[i] = new TH2F(Form("hpbpb_matrix_HLT_BothSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    // hpbpb_matrix_HLT_5TrigIneff_Smear[i] = new TH2F(Form("hpbpb_matrix_HLT_5TrigIneffSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    // hpbpb_matrix_HLT_10TrigIneff_Smear[i] = new TH2F(Form("hpbpb_matrix_HLT_10TrigIneffSmear_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);



    // hpbpb_anaBin_JetComb_gen[i] = new TH1F(Form("hpbpb_anaBin_JetComb_gen_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Gen refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    // hpbpb_anaBin_JetComb_reco[i] = new TH1F(Form("hpbpb_anaBin_JetComb_reco_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("reco jtpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    // hpbpb_anaBin_Jet80_gen[i] = new TH1F(Form("hpbpb_anaBin_Jet80_gen_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Gen refpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    // hpbpb_anaBin_Jet80_reco[i] = new TH1F(Form("hpbpb_anaBin_Jet80_reco_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("reco jtpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    // hpbpb_anaBin_Jet65_gen[i] = new TH1F(Form("hpbpb_anaBin_Jet65_gen_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Gen refpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    // hpbpb_anaBin_Jet65_reco[i] = new TH1F(Form("hpbpb_anaBin_Jet65_reco_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("reco jtpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    // hpbpb_anaBin_Jet55_gen[i] = new TH1F(Form("hpbpb_anaBin_Jet55_gen_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Gen refpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    // hpbpb_anaBin_Jet55_reco[i] = new TH1F(Form("hpbpb_anaBin_Jet55_reco_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("reco jtpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);

    // hpbpb_JetComb_raw[i] = new TH1F(Form("hpbpb_JetComb_raw_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("raw jtpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    // hpbpb_Jet80_raw[i] = new TH1F(Form("hpbpb_Jet80_raw_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("raw jtpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    // hpbpb_Jet65_raw[i] = new TH1F(Form("hpbpb_Jet65_raw_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("raw jtpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    // hpbpb_Jet55_raw[i] = new TH1F(Form("hpbpb_Jet55_raw_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("raw jtpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
  
    // ****** MC CLOSURE HISTOGRAMS *******

    hpbpb_mcclosure_matrix[i] = new TH2F(Form("hpbpb_mcclosure_matrix_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix for mcclosure refpt jtpt R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000,1000, 0, 1000);
    
    // hpbpb_mcclosure_matrix_HLT[i] = new TH2F(Form("hpbpb_mcclosure_matrix_HLT_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000,1000, 0, 1000);
    // hpbpb_10GeVBins_mcclosure_matrix_HLT[i] = new TH2F(Form("hpbpb_10GeVBins_mcclosure_matrix_HLT_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),100, 0, 1000, 100, 0, 1000);
    // hpbpb_anaBin_mcclosure_matrix_HLT[i] = new TH2F(Form("hpbpb_anaBin_mcclosure_matrix_HLT_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);

    // hpbpb_mcclosure_matrix_HLT_Jet80[i] = new TH2F(Form("hpbpb_mcclosure_matrix_HLT_Jet80_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000,1000, 0, 1000);
    // hpbpb_anaBin_mcclosure_matrix_HLT_Jet80[i] = new TH2F(Form("hpbpb_anaBin_mcclosure_matrix_HLT_Jet80_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);

    // hpbpb_mcclosure_matrix_HLT_Jet65[i] = new TH2F(Form("hpbpb_mcclosure_matrix_HLT_Jet65_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000,1000, 0, 1000);
    // hpbpb_anaBin_mcclosure_matrix_HLT_Jet65[i] = new TH2F(Form("hpbpb_anaBin_mcclosure_matrix_HLT_Jet65_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);

    // hpbpb_mcclosure_matrix_HLT_Jet55[i] = new TH2F(Form("hpbpb_mcclosure_matrix_HLT_Jet55_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),1000, 0, 1000,1000, 0, 1000);
    // hpbpb_anaBin_mcclosure_matrix_HLT_Jet55[i] = new TH2F(Form("hpbpb_anaBin_mcclosure_matrix_HLT_Jet55_R%d_%s_cent%d",radius,etaWidth.c_str(),i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth.c_str(),5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);


    
  }

  // float trgObjpt_PF_RMS[nbins_cent] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  // if(radius == 2){
  //   trgObjpt_PF_RMS[0] = 1.00;
  //   trgObjpt_PF_RMS[1] = 0.22;
  //   trgObjpt_PF_RMS[2] = 0.22;
  //   trgObjpt_PF_RMS[3] = 0.22;
  //   trgObjpt_PF_RMS[4] = 0.22;
  //   trgObjpt_PF_RMS[5] = 0.22;
  // }
  // if(radius == 3){
  //   trgObjpt_PF_RMS[0] = 8.10;
  //   trgObjpt_PF_RMS[1] = 6.89;
  //   trgObjpt_PF_RMS[2] = 0.26;
  //   trgObjpt_PF_RMS[3] = 0.22;
  //   trgObjpt_PF_RMS[4] = 0.22;
  //   trgObjpt_PF_RMS[5] = 0.22;
  // }
  // if(radius == 4){
  //   trgObjpt_PF_RMS[0] = 12.47;
  //   trgObjpt_PF_RMS[1] = 9.07;
  //   trgObjpt_PF_RMS[2] = 0.26;
  //   trgObjpt_PF_RMS[3] = 0.22;
  //   trgObjpt_PF_RMS[4] = 0.22;
  //   trgObjpt_PF_RMS[5] = 0.22;
  // }

  // TH1F * hpbpb_DijetSignalSelection_LeadPFJet_gt_A[nbins_cent];
  // TH2F * hdphiptcent_DijetSignalSelection_LeadPFJet_gt_A[nbins_cent];
  // TH1F * hpbpb_DijetBKGSelection_LeadPFJet_lt_A[nbins_cent];
  // TH2F * hdphiptcent_DijetBKGSelection_LeadPFJet_lt_A[nbins_cent];

  // for(int i = 0; i<nbins_cent; ++i){
  //   hdphiptcent_DijetSignalSelection_LeadPFJet_gt_A[i] = new TH2F(Form("hdphiptcent_DijetSignalSelection_LeadPFJet_gt_A_cent%d",i), "", 100,0,200,160,0,3.2);
  //   hdphiptcent_DijetBKGSelection_LeadPFJet_lt_A[i] = new TH2F(Form("hdphiptcent_DijetBKGSelection_LeadPFJet_lt_A_cent%d",i), "", 100,0,200,160,0,3.2);

  //   hpbpb_DijetSignalSelection_LeadPFJet_gt_A[i] = new TH1F(Form("hpbpb_DijetSignalSelection_LeadPFJet_gt_A_cent%d",i), "", 1000, 0, 1000);
  //   hpbpb_DijetBKGSelection_LeadPFJet_lt_A[i] = new TH1F(Form("hpbpb_DijetBKGSelection_LeadPFJet_lt_A_cent%d",i), "", 1000, 0, 1000);
  // }
  
  
  // TFile * cent_vz_weight = TFile::Open("data_mc_cent_vz_weight.root");
  // TH1F * hCent_weight = (TH1F*)cent_vz_weight->Get("hCent_weight");
  // TH1F * hVz_weight = (TH1F*)cent_vz_weight->Get("hVz_weight");
  
  // now start the event loop for each file. 
  
  if(printDebug) cout<<"Running through all the events now"<<endl;
  Long64_t nentries = jetpbpb[0]->GetEntries();
  //Long64_t nGoodEvt = 0;
  if(printDebug) nentries = 10;
  TRandom rnd;

  for(int nEvt = 0; nEvt < nentries; ++ nEvt) {

    if(nEvt%10000 == 0)
      cout<<nEvt<<"/"<<nentries<<endl;
    if(printDebug)cout<<"nEvt = "<<nEvt<<endl;
    jetpbpb[0]->GetEntry(nEvt);
    jetpbpb[1]->GetEntry(nEvt);
    jetpbpb[2]->GetEntry(nEvt);
    jetpbpb[3]->GetEntry(nEvt);
    // weight_Tree->GetEntry(nEvt);
    

    if(pcollisionEventSelection_F==0) continue; 
    if(fabs(vz_F)>15) continue;
    //if(!isGoodEvent_eS) continue; 
    
    // jet_select->GetEntry(nGoodEvt);
    // ++nGoodEvt;
 
    int cBin = findBin(hiBin_F);//tells us the centrality of the event.
    if(cBin==-1) continue;
    
// #if 0
//     hCentrality_fromHFBound->Fill(hfbin);


//     hNpix_vs_hiBin->Fill(hiBin_F, hiNpix_F);
//     hNpixelTracks_vs_hiBin->Fill(hiBin_F, hiNpixelTracks_F);
//     hNTracks_vs_hiBin->Fill(hiBin_F, hiNtracks_F);
//     hNTracksPtCut_vs_hiBin->Fill(hiBin_F, hiNtracksPtCut_F);
//     hNTracksEtaCut_vs_hiBin->Fill(hiBin_F, hiNtracksEtaCut_F);
//     hNTracksEtaPtCut_vs_hiBin->Fill(hiBin_F, hiNtracksEtaPtCut_F);

//     hNpix_vs_HF->Fill(hiHF_F, hiNpix_F);
//     hNpixelTracks_vs_HF->Fill(hiHF_F, hiNpixelTracks_F);
//     hNTracks_vs_HF->Fill(hiHF_F, hiNtracks_F);
//     hNTracksPtCut_vs_HF->Fill(hiHF_F, hiNtracksPtCut_F);
//     hNTracksEtaCut_vs_HF->Fill(hiHF_F, hiNtracksEtaCut_F);
//     hNTracksEtaPtCut_vs_HF->Fill(hiHF_F, hiNtracksEtaPtCut_F);
// #endif

//     // if(hiBin_eS != hiBin_F || evt_eS != evt_F || run_eS != run_F || lumi_eS != lumi_F || vz_eS != vz_F) cout<<"ERROR mismatch eS, F"<<endl;
//     // if(hiBin_eS != hiBin_jS || evt_eS != evt_jS || run_eS != run_jS || lumi_eS != lumi_jS || vz_eS != vz_jS) cout<<"ERROR mismatch eS, jS"<<endl;
//     // if(hiBin_F != hiBin_jS || evt_F != evt_jS || run_F != run_jS || lumi_F != lumi_jS || vz_F != vz_jS) cout<<"ERROR mismatch F, jS"<<endl;
//     // if(hiBin_F != hiBin_w || evt_F != evt_w || lumi_F != lumi_w || vz_F != vz_w) cout<<"ERROR mismatch forest and weight tree"<<endl;
    
//     // if(printDebug) cout<<"hibin hiForest = "<<hiBin_F<<", evtTree = "<<hiBin_eS<<", jetTree = "<<hiBin_jS<<endl;
//     // if(printDebug) cout<<"evt hiForest   = "<<evt_F<<", evtTree = "<<evt_eS<<", jetTree = "<<evt_jS<<endl;
//     // if(printDebug) cout<<"lumi hiForest  = "<<lumi_F<<", evtTree = "<<lumi_eS<<", jetTree = "<<lumi_jS<<endl;
//     // if(printDebug) cout<<"vz hiForest    = "<<vz_F<<", evtTree = "<<vz_eS<<", jetTree = "<<vz_jS<<endl;
    
//     // if(printDebug) cout<<"nref_F = "<<nref_F<<", nref_eS = "<<nref_eS<<", nref_jS = "<<nref_jS<<endl;

//     // if(nref_F != nref_eS) cout<<"ERROR mismatch in jet counts"<<endl;

//     double cent_weight = hCent_weight->GetBinContent(hCent_weight->FindBin(hiBin_F));
//     double vz_weight = hVz_weight->GetBinContent(hVz_weight->FindBin(vz_F));
  
//     double new_weight = pthatweight * cent_weight * vz_weight; 
//     //double new_weight = pthatweight; 
//     //double new_weight = weight_xsec[startfile] * cent_weight * vz_weight; 

    // double weight_eS = pthatweight;
    double weight_eS = 1.0;
    
    hpbpb_pthat[cBin]->Fill(pthat_F, weight_eS);
    
    // if(pthat_F < pthat[startfile] || pthat_F >= pthat[endfile]) continue;

    // double cweight = exp(- pow(varbin+1.11957e+01,2) / pow(1.34120e+01,2) / 2);

    // double weight = 1;
    // weight = cweight * fVz->Eval(vz_F);    
    // double scale = xsecs[]
    
    
    /*    //! Sort the jetTree jets according to pT, from the evtTree now. 
    std::vector < Jet > vJet;
    for(int jet2 = 0; jet2<nref_jS; ++jet2){
      if(printDebug) cout <<"\t \t jetTree *** "<< jet2 <<  ", pT " << pt_jS[jet2] <<", eta : "<< eta_jS[jet2] <<  ", chSum : "<< chSum_jS[jet2]  << endl;
      //if(printDebug) cout <<"\t \t evtTree *** "<< jet2 <<  ", pT " << pt_eS[jet2]<< endl;
      Jet ijet;
      ijet.id = jet2;
      ijet.pt = pt_jS[jet2];
      vJet.push_back(ijet);
    }
    std::sort (vJet.begin(), vJet.end(), compare_pt);
    std::vector < Jet >::const_iterator itr;

    // after sorting the jet array according to pT
    if(printDebug){
      cout<<"**** AFTER  SORTING ****"<<endl;
      for(itr=vJet.begin(); itr!=vJet.end(); ++itr){
	int jetLoc = (*itr).id;
	cout <<"\t \t jetTree *** "<< jetLoc <<  ", pT " << pt_jS[jetLoc]<< endl;
      }
    }
    int jet=0;
    if(printDebug) cout<<"Number of jets = "<< vJet.size()<<endl;
    for(itr=vJet.begin(); itr!=vJet.end(); ++itr, ++jet){

      bool PFElecCut = false;
      int jetLoc = (*itr).id;
      //if(isMiMatch_eS[jetLoc]) ++itr;
      if(isMultiMatch_jS[jetLoc]) {
	++itr;
	jetLoc = (*itr).id;
	if(itr == vJet.end())  break;
      }
      Float_t Sumcand = chSum_jS[jetLoc] + phSum_jS[jetLoc] + neSum_jS[jetLoc] + muSum_jS[jetLoc];
      
      if(fabs(eta_jS[jetLoc]) > 2) continue;
      // //if(isPFElecCut_eS[jet] != 1) continue;
      if(pt_jS[jetLoc] > 2 * pthat_F) continue;
      int refid = -1;
      refid = pfrefidx_jS[jetLoc];
      //if(refid < 0) continue;
      
      //if(printDebug) cout<<"pfrefidx_jS[jetLoc] = "<<pfrefidx_jS[jetLoc]<<", refdrjt_jS[refid] = "<<refdrjt_jS[refid]<<endl;

      Float_t delR = deltaR(eta_F[jet], phi_F[jet], refeta_F[jet], refphi_F[jet]);
      if(subid_F[jet] != 0 || refpt_F[jet] < 0 || delR > 0.2 || refpt_F[jet] > 2 * pthat_F) continue;

      if(pt_jS[jetLoc] != pt_F[jet]) {cout<<"Jets NOT Matched!!! ERROR!!!"<<endl; break;}
      
      //cout<<"refid = "<<refid<<", jetLoc = "<<jetLoc<<endl;

      if(printDebug) cout<<"jet = "<<jet<<", pt_F[jet] = "<<pt_F[jet]<<", jetLoc = "<<jetLoc<<", pt_jS[jetLoc] = "<<pt_jS[jetLoc]<<", calopt_jS[jetLoc] = "<<calopt_jS[jetLoc]<<", eMax/sumcand = "<<(Float_t)eMax_jS[jetLoc]/Sumcand<<", isCaloMatch_jS[jetLoc] = "<<isCaloMatch_jS[jetLoc]<<", PFEleccut = "<<PFElecCut<<", subid = "<<subid_jS[jetLoc]<<", refpt_F[jet] = "<< refpt_F[jet] <<", refid = "<< refid <<", " <<endl;

      // if(printDebug && (fabs(eta_jS[jet] > 2))) cout<<"jets with |eta| > 2 in jetTree"<<endl;
      // if(printDebug && (fabs(eta_F[jet] > 2)))  cout<<"jets with |eta| > 2 in Forest"<<endl;

      //if(printDebug && index_eS[jet] >= 0 )cout<<jet<<", hiForest pT = "<<pt_F[jet]<<", jetTree pT = "<<pt_jS[jetLoc]<<", electronCut = "<<isPFElecCut_eS[jetLoc]<<", Calo pT [indes_eS[jet]] = "<<calopt_F[index_eS[jet]]<<", Calo pT[jetLoc] "<<calopt_F[jetLoc]<<", onFly flag calculation = "<<PFElecCut<<", eMax from hiForest = "<<eMax_F[jet]<<", eMax from jet Tree = "<<eMax_jS[jetLoc]<<endl; 
      //if(printDebug)cout<<jet<<", hiForest pT = "<<pt_F[jet]<<", jetTree pT = "<<pt_jS[jetLoc]<<", electronCut = "<<isPFElecCut_eS[jetLoc]<<", Calo pT [index_eS[jet]] = "<<calopt_F[index_eS[jet]]<<", Calo pT[jetLoc] "<<calopt_F[jetLoc]<< ", calopt_jS[jetLoc] = "<<calopt_jS[jetLoc]<<", onFly flag calculation = "<<PFElecCut<<", eMax from hiForest = "<<eMax_F[jet]<<", eMax from jet Tree = "<<eMax_jS[jetLoc]<<endl; 

      // if(printDebug)cout<<jet<<", hiForest pT = "<<pt_F[jet]<<", jetTree pT = "<<pt_jS[jetLoc]<<", electronCut = "<<isPFElecCut_eS[jetLoc]<<", eMax from hiForest = "<<eMax_F[jet]<<", eMax from jet Tree = "<<eMax_jS[jetLoc]<<endl;
      */

    // weight_eS = new_weight;

    // std::vector <float> pt;
    // std::vector <float> phi;
    // std::vector <float> eta;
  
    for( int jet = 0; jet<nref_F; jet++ ){

      //cout  << "  subid : "  << subid_F[jet] << " refpt : "  << refpt_F[jet] << " pfpt :  " << pt_F[jet] << endl;
      if( fabs(eta_F[jet]) > 2.0 ) continue;
      if( subid_F[jet] != 0) continue;
      //if( isPFElecCut_eS[jet] != 1 ) continue; 
      if( refdrjt_F[jet] > 0.3 ) continue;
      if( pt_F[jet] > 3.0 * pthat_F ) continue;
      //if( pt_F[jet] > 3 * pthat[startfile] ) continue;


      // if( chMax_F[jet] < 7 && trkMax_F[jet] < 7 && neMax_F[jet] < 8 ) continue;
      // if( trkMax_F[jet] < 7 && neMax_F[jet] < 8 ) continue;

      // if( (fabs(refpt_F[jet]) > 3.0 * pthat_F) || (pt_F[jet]/refpt_F[jet] > 2.0) ) continue;
      // bool PFElecCut = false;
      // Float_t Sumcand = chSum_F[jet] + phSum_F[jet] + neSum_F[jet] + muSum_F[jet];
      // Float_t calopfpt = (float) calopt_eS[jet]/pt_F[jet];
      Float_t genpt = refpt_F[jet];
      Float_t recpt = pt_F[jet];
      Float_t rawpt = rawpt_F[jet];
      Float_t geneta = geneta_F[jet];
      Float_t delphi = deltaphi(phi_F[0], phi_F[1]);



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
      
      hJEC[cBin]->Fill(rawpt, eta_F[jet], (float)(recpt/rawpt));


      // int recpt_loc = 0;
      // int recgenpt_loc = 0;

      // for(int j = 0; j<ptSelection; ++j){
      // 	if(recpt > ptBoundary[j]) recpt_loc = j;
      // 	if(recpt > ptBoundary[j] && genpt > ptBoundary[j]) recgenpt_loc = j;	
      // }

      // if(isClMatch_eS[jet] == 1) {
      // 	heMaxSumcand_vs_calopfpt_RecopTSelection[recpt_loc][cBin]->Fill(calopfpt, (float)eMax_F[jet]/Sumcand);
      // 	heMaxSumcand_vs_calopfpt_RecoGenpTSelection[recgenpt_loc][cBin]->Fill(calopfpt, (float)eMax_F[jet]/Sumcand);
      // }else {
      // 	heMaxSumcand_vs_calopfpt_RecopTSelection[recpt_loc][cBin]->Fill(0.0, (float)eMax_F[jet]/Sumcand);
      // 	heMaxSumcand_vs_calopfpt_RecoGenpTSelection[recgenpt_loc][cBin]->Fill(0.0, (float)eMax_F[jet]/Sumcand);
      // }
      
      // if(isClMatch_eS[jet] == 1){
      // 	if(calopfpt > 0.5 && calopfpt <= 0.85 && eMax_F[jet]/Sumcand < ((Float_t)18/7 *(Float_t)calopfpt - (Float_t)9/7)) PFElecCut = true;
      // 	if(calopfpt > 0.85) PFElecCut = true;
      // 	if(calopfpt <= 0.5 && eMax_F[jet]/Sumcand < 0.05) PFElecCut = true;
      // }
      // if(isClMatch_eS[jet] == 0)
      // 	if(eMax_F[jet]/Sumcand < 0.05 ) PFElecCut = true;       

      // if(!PFElecCut) continue;

      // pt.push_back(recpt);
      // eta.push_back(eta_F[jet]);
      // phi.push_back(phi_F[jet]);


      hdeltaphi[cBin]->Fill(delphi);

      int etabin = -1;
      for(int bin = 0; bin<netabins_ana; ++bin){
	if(geneta > etabins_ana[bin]) etabin = bin;
      }
      if(etabin != -1) {
	if(genpt>=30.0 && genpt<50.0)hJER_eta_30pt50[etabin][cBin]->Fill((float)recpt/genpt);
	if(genpt>=80.0 && genpt<120.0)hJER_eta_80pt120[etabin][cBin]->Fill((float)recpt/genpt);
      }
      int ptbin = -1;
      for(int bin = 0; bin<nbins_ana; ++bin){
	if(genpt > ptbins_ana[bin]) ptbin = bin;
      }
      if(ptbin != -1){;
	hJER[ptbin][cBin]->Fill((float)recpt/genpt);
      }
      // if( isPFElecCut_eS[jet] != PFElecCut) cout<<"HUGE PROBLEM WITH ELECTRON CUT"<<endl;
      
      //if(genpt < 0) continue;
      
      hpbpb_gen[cBin]->Fill(genpt, weight_eS);
      hpbpb_reco[cBin]->Fill(recpt, weight_eS);
      hpbpb_matrix[cBin]->Fill(genpt, recpt, weight_eS);
      
      if(nEvt%2==0){
	hpbpb_mcclosure_data[cBin]->Fill(recpt, weight_eS);
      }else {
	hpbpb_mcclosure_matrix[cBin]->Fill(genpt, recpt, weight_eS);
	hpbpb_mcclosure_gen[cBin]->Fill(genpt, weight_eS);
	hpbpb_mcclosure_data_train[cBin]->Fill(recpt, weight_eS);
      }

      /*
      
      if(jet55_F == 1) hpbpb_Jet55Selection[cBin]->Fill(recpt, weight_eS * jet55_p_F);
      if(jet65_F == 1) hpbpb_Jet65Selection[cBin]->Fill(recpt, weight_eS * jet65_p_F);
      
      if(jet55_F == 1 && jet65_F==0 && jet80_F == 0){
	//weight_eS = weight_eS * jet55_p_F;

	hpbpb_Trans_matrix_Smooth_HLT[cBin]->Fill( recpt, genpt, weight_eS * jet55_p_F * hHistPrior_fit[cBin]->GetBinContent(hHistPrior_fit[cBin]->FindBin(genpt)));
	hpbpb_Trans_matrix_HLT_gen2pSmear[cBin]->Fill(recpt, genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), weight_eS * jet55_p_F * hHistPrior_fit[cBin]->GetBinContent(hHistPrior_fit[cBin]->FindBin(genpt)));
	hpbpb_Trans_matrix_HLT_GenSmear[cBin]->Fill( recpt ,genpt + rnd.Gaus(0,1), weight_eS * jet55_p_F * hHistPrior_fit[cBin]->GetBinContent(hHistPrior_fit[cBin]->FindBin(genpt)));
	hpbpb_Trans_matrix_HLT_RecoSmear[cBin]->Fill( recpt + rnd.Gaus(0,1),genpt, weight_eS * jet55_p_F * hHistPrior_fit[cBin]->GetBinContent(hHistPrior_fit[cBin]->FindBin(genpt)));
	hpbpb_Trans_matrix_HLT_BothSmear[cBin]->Fill( recpt + rnd.Gaus(0,1),genpt + rnd.Gaus(0,1), weight_eS * jet55_p_F * hHistPrior_fit[cBin]->GetBinContent(hHistPrior_fit[cBin]->FindBin(genpt)));

	
	hSVD_Matrix_noWeight[cBin]->Fill(genpt, recpt);
	hSVD_Matrix_Weight[cBin]->Fill(genpt, recpt, weight_eS * jet55_p_F);

	hSVD_Trans_Matrix_noWeight[cBin]->Fill(recpt, genpt);
	hSVD_Trans_Matrix_Weight[cBin]->Fill(recpt, genpt, weight_eS * jet55_p_F);
	
	hpbpb_Jet55_gen[cBin]->Fill(genpt, weight_eS * jet55_p_F);
	hpbpb_Jet55_GenSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), weight_eS * jet55_p_F);
	hpbpb_Jet55_gen2pSmear[cBin]->Fill(genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), weight_eS * jet55_p_F);
	hpbpb_Jet55_gen_5TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS * jet55_p_F);
	hpbpb_Jet55_gen_10TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.10/nbins_cent*(nbins_cent-cBin)), weight_eS * jet55_p_F);
	hpbpb_Jet55_reco_5TrigIneff_Smear[cBin]->Fill(recpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS * jet55_p_F);
	hpbpb_Jet55_reco_10TrigIneff_Smear[cBin]->Fill(recpt * (1. + 0.10/nbins_cent*(nbins_cent-cBin)), weight_eS * jet55_p_F);
	hpbpb_matrix_HLT_5TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS * jet55_p_F);
	hpbpb_matrix_HLT_10TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), weight_eS * jet55_p_F);


	hpbpb_Jet55_reco[cBin]->Fill(recpt, weight_eS * jet55_p_F);
	hpbpb_Jet55_RecoSmear[cBin]->Fill(recpt + rnd.Gaus(0,1), weight_eS * jet55_p_F);
	hpbpb_Jet55_raw[cBin]->Fill(rawpt, weight_eS * jet55_p_F);

	hpbpb_matrix_HLT_Jet55[cBin]->Fill(genpt, recpt, weight_eS * jet55_p_F);
	hpbpb_matrix_recoTrunc_HLT[cBin]->Fill(genpt, recpt, weight_eS * jet55_p_F);
	hpbpb_Trans_matrix_HLT_Jet55[cBin]->Fill( recpt, genpt, weight_eS * jet55_p_F);

	hpbpb_anaBin_matrix_HLT_Jet55[cBin]->Fill(genpt, recpt, weight_eS * jet55_p_F);
	hpbpb_anaBin_Trans_matrix_HLT_Jet55[cBin]->Fill( recpt,genpt, weight_eS * jet55_p_F);
	
	hpbpb_matrix_HLT_gen2pSmear[cBin]->Fill(genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), recpt, weight_eS * jet55_p_F);
	hpbpb_matrix_HLT_GenSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), recpt, weight_eS * jet55_p_F);
	hpbpb_matrix_HLT_RecoSmear[cBin]->Fill(genpt, recpt + rnd.Gaus(0,1), weight_eS * jet55_p_F);
	hpbpb_matrix_HLT_BothSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), recpt + rnd.Gaus(0,1), weight_eS * jet55_p_F);


	
	hpbpb_anaBin_Jet55_gen[cBin]->Fill(genpt, weight_eS * jet55_p_F);
	hpbpb_anaBin_Jet55_reco[cBin]->Fill(recpt, weight_eS * jet55_p_F);


	if(refparton_F[jet] < 9){
	  hpbpb_quark_Jet55_gen[cBin]->Fill(genpt, weight_eS * jet55_p_F);
	  hpbpb_quark_Jet55_reco[cBin]->Fill(recpt, weight_eS * jet55_p_F);
	  hpbpb_Trans_matrix_quark_HLT[cBin]->Fill(recpt, genpt, weight_eS * jet55_p_F);
	  if(nEvt%2==0){
	    hpbpb_mcclosure_Trans_matrix_quark_HLT[cBin]->Fill(recpt, genpt, weight_eS * jet55_p_F);
	  }else if(nEvt%2==1){
	    hpbpb_mcclosure_quark_Jet55_gen[cBin]->Fill(genpt, weight_eS * jet55_p_F);
	    hpbpb_mcclosure_quark_Jet55_reco[cBin]->Fill(recpt, weight_eS * jet55_p_F);
	  }
	}

	if(refparton_F[jet] == 21){
	  hpbpb_gluon_Jet55_gen[cBin]->Fill(genpt, weight_eS * jet55_p_F);
	  hpbpb_gluon_Jet55_reco[cBin]->Fill(recpt, weight_eS * jet55_p_F);
	  hpbpb_Trans_matrix_gluon_HLT[cBin]->Fill(recpt, genpt, weight_eS * jet55_p_F);
	  if(nEvt%2==0){
	    hpbpb_mcclosure_Trans_matrix_gluon_HLT[cBin]->Fill(recpt, genpt, weight_eS * jet55_p_F);
	  }else if(nEvt%2==1){
	    hpbpb_mcclosure_gluon_Jet55_gen[cBin]->Fill(genpt, weight_eS * jet55_p_F);
	    hpbpb_mcclosure_gluon_Jet55_reco[cBin]->Fill(recpt, weight_eS * jet55_p_F);
	  }
	}
	
	if(nEvt%2==0){
	  hpbpb_mcclosure_Jet55_data_5TrigIneff_Smear[cBin]->Fill(recpt* (1. + 0.5/nbins_cent*(nbins_cent-cBin)), weight_eS * jet55_p_F);
	  hpbpb_mcclosure_Jet55_data_10TrigIneff_Smear[cBin]->Fill(recpt* (1. + 0.1/nbins_cent*(nbins_cent-cBin)), weight_eS * jet55_p_F);
	  hpbpb_mcclosure_Jet55_data[cBin]->Fill(recpt, weight_eS * jet55_p_F);
	}else {
	  
	  hpbpb_mcclosure_matrix_HLT_5TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS * jet55_p_F);
	  hpbpb_mcclosure_matrix_HLT_10TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), weight_eS * jet55_p_F);

	  hpbpb_10GeVBins_mcclosure_matrix_HLT[cBin]->Fill(genpt, recpt, weight_eS * jet55_p_F);
	  
	  hpbpb_anaBin_mcclosure_matrix_HLT_Jet55[cBin]->Fill(genpt, recpt, weight_eS * jet55_p_F);
	  hpbpb_anaBin_mcclosure_Trans_matrix_HLT_Jet55[cBin]->Fill(recpt, genpt, weight_eS * jet55_p_F);
	  hpbpb_mcclosure_matrix_HLT_Jet55[cBin]->Fill(genpt, recpt, weight_eS * jet55_p_F);
	  hpbpb_mcclosure_Trans_matrix_HLT_Jet55[cBin]->Fill(recpt, genpt, weight_eS * jet55_p_F);

	  hpbpb_mcclosure_Jet55_gen_5TrigIneff_Smear[cBin]->Fill(genpt* (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS * jet55_p_F);
	  hpbpb_mcclosure_Jet55_gen[cBin]->Fill(genpt, weight_eS * jet55_p_F);
	  hpbpb_mcclosure_Jet55_gen_10TrigIneff_Smear[cBin]->Fill(genpt* (1. + 0.1/nbins_cent*(nbins_cent-cBin)), weight_eS * jet55_p_F);
	  hpbpb_mcclosure_Jet55_data_train[cBin]->Fill(recpt, weight_eS * jet55_p_F);
	}
	
	
      }// jet55 selection

      if(jet65_F == 1 && jet80_F == 0){

	hpbpb_Trans_matrix_Smooth_HLT[cBin]->Fill( recpt, genpt, weight_eS * hHistPrior_fit[cBin]->GetBinContent(hHistPrior_fit[cBin]->FindBin(genpt)));
	hpbpb_Trans_matrix_HLT_gen2pSmear[cBin]->Fill(recpt, genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), weight_eS * hHistPrior_fit[cBin]->GetBinContent(hHistPrior_fit[cBin]->FindBin(genpt)));
	hpbpb_Trans_matrix_HLT_GenSmear[cBin]->Fill( recpt ,genpt + rnd.Gaus(0,1), weight_eS * hHistPrior_fit[cBin]->GetBinContent(hHistPrior_fit[cBin]->FindBin(genpt)));
	hpbpb_Trans_matrix_HLT_RecoSmear[cBin]->Fill( recpt + rnd.Gaus(0,1),genpt, weight_eS * hHistPrior_fit[cBin]->GetBinContent(hHistPrior_fit[cBin]->FindBin(genpt)));
	hpbpb_Trans_matrix_HLT_BothSmear[cBin]->Fill( recpt + rnd.Gaus(0,1),genpt + rnd.Gaus(0,1), weight_eS * hHistPrior_fit[cBin]->GetBinContent(hHistPrior_fit[cBin]->FindBin(genpt)));

	
	hSVD_Matrix_noWeight[cBin]->Fill(genpt, recpt);
	hSVD_Matrix_Weight[cBin]->Fill(genpt, recpt, weight_eS);

	hSVD_Trans_Matrix_noWeight[cBin]->Fill(recpt, genpt);
	hSVD_Trans_Matrix_Weight[cBin]->Fill(recpt, genpt, weight_eS);
	
	hpbpb_Jet65_gen[cBin]->Fill(genpt, weight_eS);
	hpbpb_Jet65_GenSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), weight_eS);
	hpbpb_Jet65_gen2pSmear[cBin]->Fill(genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), weight_eS);
	hpbpb_matrix_HLT_5TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS);
	hpbpb_matrix_HLT_10TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), weight_eS);
	
	hpbpb_matrix_HLT_gen2pSmear[cBin]->Fill(genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), recpt, weight_eS);

	hpbpb_Jet65_reco[cBin]->Fill(recpt, weight_eS);
	hpbpb_Jet65_RecoSmear[cBin]->Fill(recpt + rnd.Gaus(0,1), weight_eS);
	hpbpb_Jet65_raw[cBin]->Fill(rawpt, weight_eS);
	
	hpbpb_matrix_HLT_Jet65[cBin]->Fill(genpt, recpt, weight_eS);
	hpbpb_matrix_recoTrunc_HLT[cBin]->Fill(genpt, recpt, weight_eS);
	hpbpb_Trans_matrix_HLT_Jet65[cBin]->Fill( recpt, genpt, weight_eS);
	//hpbpb_Trans_matrix_Smooth_HLT[cBin]->Fill( recpt, genpt, weight_eS * jet55_p_F * hHistPrior_fit[cBin]->GetBinContent(hHistPrior_fit[cBin]->FindBin(genpt)));
	hpbpb_anaBin_matrix_HLT_Jet65[cBin]->Fill(genpt, recpt, weight_eS);
	hpbpb_anaBin_Trans_matrix_HLT_Jet65[cBin]->Fill( recpt,genpt, weight_eS);

	// hpbpb_Trans_matrix_HLT_gen2pSmear[cBin]->Fill(recpt, genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), weight_eS);
	// hpbpb_Trans_matrix_HLT_GenSmear[cBin]->Fill( recpt ,genpt + rnd.Gaus(0,1), weight_eS);
	// hpbpb_Trans_matrix_HLT_RecoSmear[cBin]->Fill( recpt + rnd.Gaus(0,1),genpt, weight_eS);
	// hpbpb_Trans_matrix_HLT_BothSmear[cBin]->Fill( recpt + rnd.Gaus(0,1),genpt + rnd.Gaus(0,1), weight_eS);

	
	hpbpb_matrix_HLT_GenSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), recpt, weight_eS);
	hpbpb_matrix_HLT_RecoSmear[cBin]->Fill(genpt, recpt + rnd.Gaus(0,1), weight_eS);
	hpbpb_matrix_HLT_BothSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), recpt + rnd.Gaus(0,1), weight_eS);
	
	hpbpb_anaBin_Jet65_gen[cBin]->Fill(genpt, weight_eS);
	hpbpb_anaBin_Jet65_reco[cBin]->Fill(recpt, weight_eS);
	
	if(nEvt%2==0){
	  hpbpb_mcclosure_Jet65_data[cBin]->Fill(recpt, weight_eS);
	}else {
	  hpbpb_mcclosure_matrix_HLT_5TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS);
	  hpbpb_mcclosure_matrix_HLT_10TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), weight_eS);

	  hpbpb_10GeVBins_mcclosure_matrix_HLT[cBin]->Fill(genpt, recpt, weight_eS);
	  
	  hpbpb_anaBin_mcclosure_matrix_HLT_Jet65[cBin]->Fill(genpt, recpt, weight_eS);
	  hpbpb_anaBin_mcclosure_Trans_matrix_HLT_Jet65[cBin]->Fill(recpt, genpt, weight_eS);
	  hpbpb_mcclosure_matrix_HLT_Jet65[cBin]->Fill(genpt, recpt, weight_eS);
	  hpbpb_mcclosure_Trans_matrix_HLT_Jet65[cBin]->Fill(recpt, genpt, weight_eS);

	  hpbpb_mcclosure_Jet65_data_train[cBin]->Fill(recpt, weight_eS);
	  hpbpb_mcclosure_Jet65_gen[cBin]->Fill(genpt, weight_eS);

	}


	if(refparton_F[jet] < 9){
	  hpbpb_quark_Jet65_gen[cBin]->Fill(genpt, weight_eS);
	  hpbpb_quark_Jet65_reco[cBin]->Fill(recpt, weight_eS);
	  hpbpb_Trans_matrix_quark_HLT[cBin]->Fill(recpt, genpt, weight_eS);
	  if(nEvt%2==0){
	    hpbpb_mcclosure_Trans_matrix_quark_HLT[cBin]->Fill(recpt, genpt, weight_eS);
	  }else if(nEvt%2==1){
	    hpbpb_mcclosure_quark_Jet65_gen[cBin]->Fill(genpt, weight_eS);
	    hpbpb_mcclosure_quark_Jet65_reco[cBin]->Fill(recpt, weight_eS);
	  }
	}

	if(refparton_F[jet] == 21){
	  hpbpb_gluon_Jet65_gen[cBin]->Fill(genpt, weight_eS);
	  hpbpb_gluon_Jet65_reco[cBin]->Fill(recpt, weight_eS);
	  hpbpb_Trans_matrix_gluon_HLT[cBin]->Fill(recpt, genpt, weight_eS);
	  if(nEvt%2==0){
	    hpbpb_mcclosure_Trans_matrix_gluon_HLT[cBin]->Fill(recpt, genpt, weight_eS);
	  }else if(nEvt%2==1){
	    hpbpb_mcclosure_gluon_Jet65_gen[cBin]->Fill(genpt, weight_eS);
	    hpbpb_mcclosure_gluon_Jet65_reco[cBin]->Fill(recpt, weight_eS);
	  }
	}

	
      }// jet65 selection

      if(jet80_F == 1){
	
	hpbpb_Trans_matrix_Smooth_HLT[cBin]->Fill( recpt, genpt, weight_eS * hHistPrior_fit[cBin]->GetBinContent(hHistPrior_fit[cBin]->FindBin(genpt)));
	hpbpb_Trans_matrix_HLT_gen2pSmear[cBin]->Fill(recpt, genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), weight_eS * hHistPrior_fit[cBin]->GetBinContent(hHistPrior_fit[cBin]->FindBin(genpt)));
	hpbpb_Trans_matrix_HLT_GenSmear[cBin]->Fill( recpt ,genpt + rnd.Gaus(0,1), weight_eS * hHistPrior_fit[cBin]->GetBinContent(hHistPrior_fit[cBin]->FindBin(genpt)));
	hpbpb_Trans_matrix_HLT_RecoSmear[cBin]->Fill( recpt + rnd.Gaus(0,1),genpt, weight_eS * hHistPrior_fit[cBin]->GetBinContent(hHistPrior_fit[cBin]->FindBin(genpt)));
	hpbpb_Trans_matrix_HLT_BothSmear[cBin]->Fill( recpt + rnd.Gaus(0,1),genpt + rnd.Gaus(0,1), weight_eS * hHistPrior_fit[cBin]->GetBinContent(hHistPrior_fit[cBin]->FindBin(genpt)));


	
	hSVD_Matrix_noWeight[cBin]->Fill(genpt, recpt);
	hSVD_Matrix_Weight[cBin]->Fill(genpt, recpt, weight_eS);

	hSVD_Trans_Matrix_noWeight[cBin]->Fill(recpt, genpt);
	hSVD_Trans_Matrix_Weight[cBin]->Fill(recpt, genpt, weight_eS);
	
	hpbpb_Jet80_gen[cBin]->Fill(genpt, weight_eS);
	hpbpb_Jet80_GenSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), weight_eS);
	hpbpb_Jet80_gen2pSmear[cBin]->Fill(genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), weight_eS);
	hpbpb_matrix_HLT_5TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS);
	hpbpb_matrix_HLT_10TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), weight_eS);
	hpbpb_matrix_HLT_gen2pSmear[cBin]->Fill(genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), recpt, weight_eS);

	hpbpb_matrix_HLT_Jet80[cBin]->Fill(genpt, recpt, weight_eS);
	hpbpb_matrix_recoTrunc_HLT[cBin]->Fill(genpt, recpt, weight_eS);
	hpbpb_Trans_matrix_HLT_Jet80[cBin]->Fill( recpt, genpt, weight_eS);
	// hpbpb_Trans_matrix_Smooth_HLT[cBin]->Fill( recpt, genpt, weight_eS * jet55_p_F * hHistPrior_fit[cBin]->GetBinContent(hHistPrior_fit[cBin]->FindBin(genpt)));
	hpbpb_anaBin_matrix_HLT_Jet80[cBin]->Fill(genpt, recpt, weight_eS);
	hpbpb_anaBin_Trans_matrix_HLT_Jet80[cBin]->Fill( recpt,genpt, weight_eS);

	// hpbpb_Trans_matrix_HLT_gen2pSmear[cBin]->Fill(recpt, genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), weight_eS);
	// hpbpb_Trans_matrix_HLT_GenSmear[cBin]->Fill( recpt ,genpt + rnd.Gaus(0,1), weight_eS);
	// hpbpb_Trans_matrix_HLT_RecoSmear[cBin]->Fill( recpt + rnd.Gaus(0,1),genpt, weight_eS);
	// hpbpb_Trans_matrix_HLT_BothSmear[cBin]->Fill( recpt + rnd.Gaus(0,1),genpt + rnd.Gaus(0,1), weight_eS);

	
	hpbpb_Jet80_reco[cBin]->Fill(recpt, weight_eS);
	hpbpb_Jet80_RecoSmear[cBin]->Fill(recpt + rnd.Gaus(0,1), weight_eS);
	hpbpb_Jet80_raw[cBin]->Fill(rawpt, weight_eS);

	hpbpb_matrix_HLT_GenSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), recpt, weight_eS);
	hpbpb_matrix_HLT_RecoSmear[cBin]->Fill(genpt, recpt + rnd.Gaus(0,1), weight_eS);
	hpbpb_matrix_HLT_BothSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), recpt + rnd.Gaus(0,1), weight_eS);
	
	hpbpb_anaBin_Jet80_gen[cBin]->Fill(genpt, weight_eS);
	hpbpb_anaBin_Jet80_reco[cBin]->Fill(recpt, weight_eS);
        
	if(nEvt%2==0){
	  hpbpb_mcclosure_Jet80_data[cBin]->Fill(recpt, weight_eS);
	}else {
	  hpbpb_mcclosure_matrix_HLT_5TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS);
	  hpbpb_mcclosure_matrix_HLT_10TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), weight_eS);

	  hpbpb_10GeVBins_mcclosure_matrix_HLT[cBin]->Fill(genpt, recpt, weight_eS);
	  
	  hpbpb_anaBin_mcclosure_matrix_HLT_Jet80[cBin]->Fill(genpt, recpt, weight_eS);
	  hpbpb_anaBin_mcclosure_Trans_matrix_HLT_Jet80[cBin]->Fill(recpt, genpt, weight_eS);
	  hpbpb_mcclosure_matrix_HLT_Jet80[cBin]->Fill(genpt, recpt, weight_eS);
	  hpbpb_mcclosure_Trans_matrix_HLT_Jet80[cBin]->Fill(recpt, genpt, weight_eS);
	 
	  hpbpb_mcclosure_Jet80_data_train[cBin]->Fill(recpt, weight_eS);
	  hpbpb_mcclosure_Jet80_gen[cBin]->Fill(genpt, weight_eS);

	}

	if(refparton_F[jet] < 9){
	  hpbpb_quark_Jet80_gen[cBin]->Fill(genpt, weight_eS);
	  hpbpb_quark_Jet80_reco[cBin]->Fill(recpt, weight_eS);
	  hpbpb_Trans_matrix_quark_HLT[cBin]->Fill(recpt, genpt, weight_eS);
	  if(nEvt%2==0){
	    hpbpb_mcclosure_Trans_matrix_quark_HLT[cBin]->Fill(recpt, genpt, weight_eS);
	  }else if(nEvt%2==1){
	    hpbpb_mcclosure_quark_Jet80_gen[cBin]->Fill(genpt, weight_eS);
	    hpbpb_mcclosure_quark_Jet80_reco[cBin]->Fill(recpt, weight_eS);
	  }
	}

	if(refparton_F[jet] == 21){
	  hpbpb_gluon_Jet80_gen[cBin]->Fill(genpt, weight_eS);
	  hpbpb_gluon_Jet80_reco[cBin]->Fill(recpt, weight_eS);
	  hpbpb_Trans_matrix_gluon_HLT[cBin]->Fill(recpt, genpt, weight_eS);
	  if(nEvt%2==0){
	    hpbpb_mcclosure_Trans_matrix_gluon_HLT[cBin]->Fill(recpt, genpt, weight_eS);
	  }else if(nEvt%2==1){
	    hpbpb_mcclosure_gluon_Jet80_gen[cBin]->Fill(genpt, weight_eS);
	    hpbpb_mcclosure_gluon_Jet80_reco[cBin]->Fill(recpt, weight_eS);
	  }
	}

	
      }// jet80 selection

      */
      
    }// jet loop
    
    if(printDebug)cout<<endl;

    // if(pt.size() <= 1) continue;
    // double delphi_lead_sublead = deltaphi(phi[0], phi[1]);

    // bool isSignal = true;
    // if((pt[0] > 3 * trgObjpt_PF_RMS[cBin]) && (delphi_lead_sublead < 2 * (float)TMath::Pi()/3))
    //   isSignal = false;
    // if(pt[0] > 3 * trgObjpt_PF_RMS[cBin] && delphi_lead_sublead > 2 * (float)TMath::Pi()/3 && pt[1] < 1.6 * trgObjpt_PF_RMS[cBin]) isSignal = false;
    // if(pt[0] < 3 * trgObjpt_PF_RMS[cBin] ) isSignal = false;

    // for(unsigned  j = 1; j<pt.size(); ++j){
    //   float delphi = deltaphi (phi[0], phi[j]);
      
    //   if(isSignal) hdphiptcent_DijetSignalSelection_LeadPFJet_gt_A[cBin]->Fill(pt[0], delphi, weight_eS);
    //   if(!isSignal) hdphiptcent_DijetBKGSelection_LeadPFJet_lt_A[cBin]->Fill(pt[0], delphi, weight_eS);

    // }

    // for(unsigned j = 0; j<pt.size(); ++j){
    //   if(isSignal) hpbpb_DijetSignalSelection_LeadPFJet_gt_A[cBin]->Fill(pt[j], weight_eS);
    //   if(!isSignal) hpbpb_DijetBKGSelection_LeadPFJet_lt_A[cBin]->Fill(pt[j], weight_eS);
    // }
    
    // pt.clear();
    // eta.clear();
    // phi.clear();
    
  }// event loop

  // for(int i = 0; i<nbins_cent; ++i){

    // hpbpb_quark_JetComb_gen[i]->Add(hpbpb_quark_Jet80_gen[i]);
    // hpbpb_quark_JetComb_gen[i]->Add(hpbpb_quark_Jet65_gen[i]);
    // hpbpb_quark_JetComb_gen[i]->Add(hpbpb_quark_Jet55_gen[i]);

    // hpbpb_gluon_JetComb_gen[i]->Add(hpbpb_gluon_Jet80_gen[i]);
    // hpbpb_gluon_JetComb_gen[i]->Add(hpbpb_gluon_Jet65_gen[i]);
    // hpbpb_gluon_JetComb_gen[i]->Add(hpbpb_gluon_Jet55_gen[i]);

    // hpbpb_quark_JetComb_reco[i]->Add(hpbpb_quark_Jet80_reco[i]);
    // hpbpb_quark_JetComb_reco[i]->Add(hpbpb_quark_Jet65_reco[i]);
    // hpbpb_quark_JetComb_reco[i]->Add(hpbpb_quark_Jet55_reco[i]);

    // hpbpb_gluon_JetComb_reco[i]->Add(hpbpb_gluon_Jet80_reco[i]);
    // hpbpb_gluon_JetComb_reco[i]->Add(hpbpb_gluon_Jet65_reco[i]);
    // hpbpb_gluon_JetComb_reco[i]->Add(hpbpb_gluon_Jet55_reco[i]);

    // hpbpb_mcclosure_quark_JetComb_gen[i]->Add(hpbpb_mcclosure_quark_Jet80_gen[i]);
    // hpbpb_mcclosure_quark_JetComb_gen[i]->Add(hpbpb_mcclosure_quark_Jet65_gen[i]);
    // hpbpb_mcclosure_quark_JetComb_gen[i]->Add(hpbpb_mcclosure_quark_Jet55_gen[i]);

    // hpbpb_mcclosure_gluon_JetComb_gen[i]->Add(hpbpb_mcclosure_gluon_Jet80_gen[i]);
    // hpbpb_mcclosure_gluon_JetComb_gen[i]->Add(hpbpb_mcclosure_gluon_Jet65_gen[i]);
    // hpbpb_mcclosure_gluon_JetComb_gen[i]->Add(hpbpb_mcclosure_gluon_Jet55_gen[i]);

    // hpbpb_mcclosure_quark_JetComb_reco[i]->Add(hpbpb_mcclosure_quark_Jet80_reco[i]);
    // hpbpb_mcclosure_quark_JetComb_reco[i]->Add(hpbpb_mcclosure_quark_Jet65_reco[i]);
    // hpbpb_mcclosure_quark_JetComb_reco[i]->Add(hpbpb_mcclosure_quark_Jet55_reco[i]);

    // hpbpb_mcclosure_gluon_JetComb_reco[i]->Add(hpbpb_mcclosure_gluon_Jet80_reco[i]);
    // hpbpb_mcclosure_gluon_JetComb_reco[i]->Add(hpbpb_mcclosure_gluon_Jet65_reco[i]);
    // hpbpb_mcclosure_gluon_JetComb_reco[i]->Add(hpbpb_mcclosure_gluon_Jet55_reco[i]);
    
    // hpbpb_mcclosure_JetComb_data[i]->Add(hpbpb_mcclosure_Jet80_data[i]);
    // hpbpb_mcclosure_JetComb_data[i]->Add(hpbpb_mcclosure_Jet65_data[i]);
    // hpbpb_mcclosure_JetComb_data[i]->Add(hpbpb_mcclosure_Jet55_data[i]);
    // divideBinWidth(hpbpb_mcclosure_JetComb_data[i]);
        
    // hpbpb_mcclosure_JetComb_data_5TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet80_data[i]);
    // hpbpb_mcclosure_JetComb_data_5TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet65_data[i]);
    // hpbpb_mcclosure_JetComb_data_5TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet55_data_5TrigIneff_Smear[i]);
    // divideBinWidth(hpbpb_mcclosure_JetComb_data_5TrigIneff_Smear[i]);

    // hpbpb_mcclosure_JetComb_data_10TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet80_data[i]);
    // hpbpb_mcclosure_JetComb_data_10TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet65_data[i]);
    // hpbpb_mcclosure_JetComb_data_10TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet55_data_10TrigIneff_Smear[i]);
    // divideBinWidth(hpbpb_mcclosure_JetComb_data_10TrigIneff_Smear[i]);
    
    // hpbpb_mcclosure_JetComb_data_train[i]->Add(hpbpb_mcclosure_Jet80_data_train[i]);
    // hpbpb_mcclosure_JetComb_data_train[i]->Add(hpbpb_mcclosure_Jet65_data_train[i]);
    // hpbpb_mcclosure_JetComb_data_train[i]->Add(hpbpb_mcclosure_Jet55_data_train[i]);
    // divideBinWidth(hpbpb_mcclosure_JetComb_data_train[i]);
    
    // hpbpb_mcclosure_JetComb_gen[i]->Add(hpbpb_mcclosure_Jet80_gen[i]);    
    // hpbpb_mcclosure_JetComb_gen[i]->Add(hpbpb_mcclosure_Jet65_gen[i]);
    // hpbpb_mcclosure_JetComb_gen[i]->Add(hpbpb_mcclosure_Jet55_gen[i]);
    // divideBinWidth(hpbpb_mcclosure_JetComb_gen[i]);

    // hpbpb_mcclosure_JetComb_gen_10TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet80_gen[i]);
    // hpbpb_mcclosure_JetComb_gen_10TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet65_gen[i]);
    // hpbpb_mcclosure_JetComb_gen_10TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet55_gen_10TrigIneff_Smear[i]);
    // divideBinWidth(hpbpb_mcclosure_JetComb_gen_10TrigIneff_Smear[i]);

    // hpbpb_mcclosure_JetComb_gen_5TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet80_gen[i]);
    // hpbpb_mcclosure_JetComb_gen_5TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet65_gen[i]);
    // hpbpb_mcclosure_JetComb_gen_5TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet55_gen_5TrigIneff_Smear[i]);
    // divideBinWidth(hpbpb_mcclosure_JetComb_gen_5TrigIneff_Smear[i]);
    
    // hpbpb_JetComb_reco[i]->Add(hpbpb_Jet80_reco[i]);
    // hpbpb_JetComb_reco[i]->Add(hpbpb_Jet65_reco[i]);
    // hpbpb_JetComb_reco[i]->Add(hpbpb_Jet55_reco[i]);
    // divideBinWidth(hpbpb_JetComb_reco[i]);    

    // hpbpb_JetComb_reco_5TrigIneff_Smear[i]->Add(hpbpb_Jet80_reco[i]);
    // hpbpb_JetComb_reco_5TrigIneff_Smear[i]->Add(hpbpb_Jet65_reco[i]);
    // hpbpb_JetComb_reco_5TrigIneff_Smear[i]->Add(hpbpb_Jet55_reco_5TrigIneff_Smear[i]);
    
    // hpbpb_JetComb_gen_5TrigIneff_Smear[i]->Add(hpbpb_Jet80_gen[i]);
    // hpbpb_JetComb_gen_5TrigIneff_Smear[i]->Add(hpbpb_Jet65_gen[i]);
    // hpbpb_JetComb_gen_5TrigIneff_Smear[i]->Add(hpbpb_Jet55_gen_5TrigIneff_Smear[i]);
    
    // hpbpb_JetComb_reco_10TrigIneff_Smear[i]->Add(hpbpb_Jet80_reco[i]);
    // hpbpb_JetComb_reco_10TrigIneff_Smear[i]->Add(hpbpb_Jet65_reco[i]);
    // hpbpb_JetComb_reco_10TrigIneff_Smear[i]->Add(hpbpb_Jet55_reco_10TrigIneff_Smear[i]);

    // hpbpb_JetComb_gen_10TrigIneff_Smear[i]->Add(hpbpb_Jet80_gen[i]);
    // hpbpb_JetComb_gen_10TrigIneff_Smear[i]->Add(hpbpb_Jet65_gen[i]);
    // hpbpb_JetComb_gen_10TrigIneff_Smear[i]->Add(hpbpb_Jet55_gen_10TrigIneff_Smear[i]);

    // hpbpb_JetComb_RecoSmear[i]->Add(hpbpb_Jet80_RecoSmear[i]);
    // hpbpb_JetComb_RecoSmear[i]->Add(hpbpb_Jet65_RecoSmear[i]);
    // hpbpb_JetComb_RecoSmear[i]->Add(hpbpb_Jet55_RecoSmear[i]);
    // divideBinWidth(hpbpb_JetComb_RecoSmear[i]);

    // hpbpb_JetComb_GenSmear[i]->Add(hpbpb_Jet80_GenSmear[i]);
    // hpbpb_JetComb_GenSmear[i]->Add(hpbpb_Jet65_GenSmear[i]);    
    // hpbpb_JetComb_GenSmear[i]->Add(hpbpb_Jet55_GenSmear[i]);
    // divideBinWidth(hpbpb_JetComb_GenSmear[i]);

    // hpbpb_JetComb_gen2pSmear[i]->Add(hpbpb_Jet80_gen2pSmear[i]);
    // hpbpb_JetComb_gen2pSmear[i]->Add(hpbpb_Jet65_gen2pSmear[i]);
    // hpbpb_JetComb_gen2pSmear[i]->Add(hpbpb_Jet55_gen2pSmear[i]);
    // divideBinWidth(hpbpb_JetComb_gen2pSmear[i]);
    
    
    // hpbpb_JetComb_gen[i]->Add(hpbpb_Jet80_gen[i]);
    // hpbpb_JetComb_gen[i]->Add(hpbpb_Jet65_gen[i]);
    // hpbpb_JetComb_gen[i]->Add(hpbpb_Jet55_gen[i]);
    // divideBinWidth(hpbpb_JetComb_gen[i]);

    // hpbpb_JetComb_raw[i]->Add(hpbpb_Jet80_raw[i]);
    // hpbpb_JetComb_raw[i]->Add(hpbpb_Jet65_raw[i]);
    // hpbpb_JetComb_raw[i]->Add(hpbpb_Jet55_raw[i]);    
    // divideBinWidth(hpbpb_JetComb_raw[i]);
    
    // divideBinWidth(hpbpb_reco[i]);
    
    // divideBinWidth(hpbpb_gen[i]);

    // hpbpb_anaBin_JetComb_reco[i]->Add(hpbpb_anaBin_Jet80_reco[i]);
    // hpbpb_anaBin_JetComb_reco[i]->Add(hpbpb_anaBin_Jet65_reco[i]);    
    // hpbpb_anaBin_JetComb_reco[i]->Add(hpbpb_anaBin_Jet55_reco[i]);
    // divideBinWidth(hpbpb_anaBin_JetComb_gen[i]);
    
    // hpbpb_anaBin_JetComb_gen[i]->Add(hpbpb_anaBin_Jet80_gen[i]);
    // hpbpb_anaBin_JetComb_gen[i]->Add(hpbpb_anaBin_Jet65_gen[i]);
    // hpbpb_anaBin_JetComb_gen[i]->Add(hpbpb_anaBin_Jet55_gen[i]);    
    // divideBinWidth(hpbpb_anaBin_JetComb_reco[i]);

    
    // hpbpb_mcclosure_matrix_HLT[i]->Add(hpbpb_mcclosure_matrix_HLT_Jet80[i]);
    // hpbpb_mcclosure_matrix_HLT[i]->Add(hpbpb_mcclosure_matrix_HLT_Jet65[i]);
    // hpbpb_mcclosure_matrix_HLT[i]->Add(hpbpb_mcclosure_matrix_HLT_Jet55[i]);

    // hpbpb_anaBin_mcclosure_matrix_HLT[i]->Add(hpbpb_anaBin_mcclosure_matrix_HLT_Jet80[i]);
    // hpbpb_anaBin_mcclosure_matrix_HLT[i]->Add(hpbpb_anaBin_mcclosure_matrix_HLT_Jet65[i]);
    // hpbpb_anaBin_mcclosure_matrix_HLT[i]->Add(hpbpb_anaBin_mcclosure_matrix_HLT_Jet55[i]);
    
    // hpbpb_mcclosure_Trans_matrix_HLT[i]->Add(hpbpb_mcclosure_Trans_matrix_HLT_Jet80[i]);
    // hpbpb_mcclosure_Trans_matrix_HLT[i]->Add(hpbpb_mcclosure_Trans_matrix_HLT_Jet65[i]);
    // hpbpb_mcclosure_Trans_matrix_HLT[i]->Add(hpbpb_mcclosure_Trans_matrix_HLT_Jet55[i]);

    // hpbpb_anaBin_mcclosure_Trans_matrix_HLT[i]->Add(hpbpb_anaBin_mcclosure_Trans_matrix_HLT_Jet80[i]);
    // hpbpb_anaBin_mcclosure_Trans_matrix_HLT[i]->Add(hpbpb_anaBin_mcclosure_Trans_matrix_HLT_Jet65[i]);
    // hpbpb_anaBin_mcclosure_Trans_matrix_HLT[i]->Add(hpbpb_anaBin_mcclosure_Trans_matrix_HLT_Jet55[i]);

    // hpbpb_matrix_HLT[i]->Add(hpbpb_matrix_HLT_Jet80[i]);
    // hpbpb_matrix_HLT[i]->Add(hpbpb_matrix_HLT_Jet65[i]);
    // hpbpb_matrix_HLT[i]->Add(hpbpb_matrix_HLT_Jet55[i]);

    // hpbpb_anaBin_matrix_HLT[i]->Add(hpbpb_anaBin_matrix_HLT_Jet80[i]);
    // hpbpb_anaBin_matrix_HLT[i]->Add(hpbpb_anaBin_matrix_HLT_Jet65[i]);
    // hpbpb_anaBin_matrix_HLT[i]->Add(hpbpb_anaBin_matrix_HLT_Jet55[i]);

    // hpbpb_Trans_matrix_HLT[i]->Add(hpbpb_Trans_matrix_HLT_Jet80[i]);
    // hpbpb_Trans_matrix_HLT[i]->Add(hpbpb_Trans_matrix_HLT_Jet65[i]);
    // hpbpb_Trans_matrix_HLT[i]->Add(hpbpb_Trans_matrix_HLT_Jet55[i]);

    // hpbpb_anaBin_Trans_matrix_HLT[i]->Add(hpbpb_anaBin_Trans_matrix_HLT_Jet80[i]);
    // hpbpb_anaBin_Trans_matrix_HLT[i]->Add(hpbpb_anaBin_Trans_matrix_HLT_Jet65[i]);
    // hpbpb_anaBin_Trans_matrix_HLT[i]->Add(hpbpb_anaBin_Trans_matrix_HLT_Jet55[i]);


  // }

  fout->Write();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}//macro end
