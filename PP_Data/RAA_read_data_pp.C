// original author: Raghav Kunnawalkam Elayavalli
// Nov 24th 2015
// destroyed by Ian Laflotte
// Apr 29th 2016
// reads and writes jets from pp data forest files
// for producing quality assurance spectra, and for unfolding later

//C includes			    
#include <cstdlib>		    
#include <stdio.h>		    
#include <iostream>		    
#include <fstream>		    
#include <sstream>		    
#include <cmath>		    
				    
//load/read/write		    
#include <TSystem.h>		    
#include <TProfile.h>		    
#include <TFile.h>		    
#include <TTree.h>		    
#include <TNtuple.h>		    
				    
//utilities			    
#include <TChain.h>		    
#include <TCut.h>		    
#include <TMath.h>		    
#include <TRandom3.h>		    
#include <TStopwatch.h>		    
#include <TEventList.h>		    
				    
//plotting			    
#include <TStyle.h>		    
#include <TCanvas.h>		    
#include <TGraphErrors.h>	    
#include <TGraphAsymmErrors.h>	    
#include <TLegend.h>		    
#include <TLatex.h>		    
#include <TLine.h>		    

//histos				    
#include <TH1.h>		    
#include <TH1F.h>		    
#include <TF1.h>		    
#include <TH2.h>		    
#include <TH2F.h>                   
#include <TH3.h>		    

//declare helper functions
double trigComb(bool *trg, int *pscl, double pt);
void divideBinWidth(TH1 *h);

//switches
const bool printDebug = true; //int iii=0;    /*debug*/    //std::cout<<"here"<<iii<<std::endl;iii++;

//filelist to run over
const std::string infile_Forest = "testFile_2015_5p02TeVpp_data_forests.txt";

//default output file name
const std::string defaultName = "test_output.root"; 

//pt binning
const int ptbins[] = {15, 30, 50, 80, 120, 170, 220, 300, 500};
const int nbins_pt = sizeof(ptbins)/sizeof(int)-1;

//eta binning
static const char *etaWidth = (char*)"20_eta_20";

//root file directories + tree names
const std::string defDirNames[]={  "hltanalysis"    ,
				   "skimanalysis"   ,
				   "GARBAGE ENTRY"  , //later filled with Form("ak%d%sJetAnalyzer",radius, jetType.c_str()),
				   //"ppTrack"        ,
				   "hiEvtAnalyzer"  ,
				   "hltobject"      , // hlt40
				   "hltobject"      , // hlt60
				   "hltobject"      , // hlt80
				   "hltobject"      }; // hlt100

const std::string trees[]={ "HltTree"                    ,
			    "HltTree"                    ,
			    "t"                          , 
			    //"trackTree"                  ,
			    "HiTree"                     ,
			    "HLT_AK4CaloJet40_Eta5p1_v"  ,
			    "HLT_AK4CaloJet60_Eta5p1_v"  ,
			    "HLT_AK4CaloJet80_Eta5p1_v"  ,
			    "HLT_AK4CaloJet100_Eta5p1_v" };
const int N = sizeof(defDirNames)/sizeof(std::string);//note, for strings, don't take off 1?

// Add the Jet ID plots:
const std::string var[] = {"jtpt" , "rawpt" , "jteta"   , "jtphi", "trkMax", "trkSum"   , "trkHardSum", 
			   "chMax", "chSum", "chHardSum", "phMax", "phSum" , "phHardSum", "neMax"     , 
			   "neSum", "eMax" , "eSum"     , "muMax", "muSum" ,"Aj"        ,"xj"         };
const int N_vars = sizeof(var)/sizeof(std::string);

//RAA_read_data_pp
int RAA_read_data_pp(int startfile = 0 , int endfile = 1            , 
		     int radius = 4 , std::string jetType = "PF" , std::string kFoname = defaultName){ 
  
  //for timing the script
  //TDatime date;  
  TStopwatch timer;
  timer.Start();
  
  //plto settings for later
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  //basic inf othe screen
  std::cout<<"reading from file "<< startfile << " to file " << endfile<<std::endl;
  std::cout<<"radius = "    << radius<<std::endl;
  std::cout<<"N = "         << N<<std::endl;
    
  //initialize tree name array
  std::string dir[N];
  for(int i=0;i<N;i++){
    if(i!=2)dir[i]=defDirNames[i];
    else dir[i]=Form("ak%d%sJetAnalyzer",radius, jetType.c_str()); 
  }
  
  //the file list
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;  
  
  for(int ifile = 0;ifile<startfile;ifile++) instr_Forest >> filename_Forest;  //what's this for???
  
  //chain the trees together
  TChain* jetpp[N];
  for(int t = 0;t<N;t++)  jetpp[t] = new TChain(string(dir[t]+"/"+trees[t]).data());
  
  for(int ifile = startfile; ifile<endfile; ++ifile){
    instr_Forest>>filename_Forest;    
    for(int t = 0;t<N;t++)  jetpp[t]->Add(filename_Forest.c_str());
    
    if(printDebug){//check that the trees load properly
      for(int t = 0;t<N;t++){ 
	std::cout << "Tree loaded  " << string(dir[t]+"/"+trees[t]).data() << std::endl;
	std::cout << "Entries : "    << jetpp[t]->GetEntries()             << std::endl;
      }
      
//      std::cout << "Tree loaded  " << string(dir[1]+"/"+trees[1]).data() << std::endl;
//      std::cout << "Entries : " << jetpp[1]->GetEntries() << std::endl;
//      std::cout << "Tree loaded  " << string(dir[2]+"/"+trees[2]).data() << std::endl;
//      std::cout << "Entries : " << jetpp[2]->GetEntries() << std::endl;
//      std::cout << "Tree loaded  " << string(dir[3]+"/"+trees[3]).data() << std::endl;
//      std::cout << "Entries : " << jetpp[3]->GetEntries() << std::endl;
//      std::cout << "Tree loaded  " << string(dir[4]+"/"+trees[4]).data() << std::endl;
//      std::cout << "Entries : " << jetpp[4]->GetEntries() << std::endl;
//      std::cout << "Tree loaded  " << string(dir[5]+"/"+trees[5]).data() << std::endl;
//      std::cout << "Entries : " << jetpp[5]->GetEntries() << std::endl;
//      std::cout << "Tree loaded  " << string(dir[6]+"/"+trees[6]).data() << std::endl;
//      std::cout << "Entries : " << jetpp[6]->GetEntries() << std::endl;
//      std::cout << "Tree loaded  " << string(dir[7]+"/"+trees[7]).data() << std::endl;
//      std::cout << "Entries : " << jetpp[7]->GetEntries() << std::endl; 
    
    }
  }
  
  //friend jet pp trees to each other
  for(int t = 0;t<N;t++){ 
    if(t==2) continue;
    else jetpp[2]->AddFriend( jetpp[t] );
  }

  // HiEvtAnalyzer
  ULong64_t evt_F;
  UInt_t run_F;
  UInt_t lumi_F;
  float vz_F;
  
  // in skimanalysis
  int pBeamScrapingFilter_F;
  int pHBHENoiseFilter_F;
  int pprimaryvertexFilter_F;

  // in hltanalysis
  //L1
  int jet40_l1seed_p_F;
  int jet60_l1seed_p_F;
  int jet80_l1seed_p_F;
  int jet100_l1seed_p_F;  

  //HLT
  int jet40_F;
  int jet60_F;
  int jet80_F;
  int jet100_F;
  int jet40_p_F;
  int jet60_p_F;
  int jet80_p_F;
  int jet100_p_F;

  // in jet analyzer of choice
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

  //hiEvtAnalyzer
  jetpp[3]->SetBranchAddress("evt",&evt_F);
  jetpp[3]->SetBranchAddress("run",&run_F);
  jetpp[3]->SetBranchAddress("lumi",&lumi_F);
  jetpp[3]->SetBranchAddress("vz",&vz_F);

  //skimanalysis
  jetpp[1]->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_F);
  jetpp[1]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilter_F);
  jetpp[1]->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryvertexFilter_F);

  //hltanalysis
  jetpp[0]->SetBranchAddress("L1_SingleJet28_BptxAND_Prescl",&jet40_l1seed_p_F);
  jetpp[0]->SetBranchAddress("L1_SingleJet40_BptxAND_Prescl",&jet60_l1seed_p_F);
  jetpp[0]->SetBranchAddress("L1_SingleJet48_BptxAND_Prescl",&jet80_l1seed_p_F);
  jetpp[0]->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&jet100_l1seed_p_F);

  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet40_Eta5p1_v1",&jet40_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet40_Eta5p1_v1_Prescl",&jet40_p_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet60_Eta5p1_v1",&jet60_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet60_Eta5p1_v1_Prescl",&jet60_p_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1",&jet80_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1_Prescl",&jet80_p_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet100_Eta5p1_v1",&jet100_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet100_Eta5p1_v1_Prescl",&jet100_p_F);

  //ak%d%sJetAnalyzer
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
    
  //for specific HLT branches
  std::vector<float> *trgObjpt_40 = 0;   std::vector<float> *trgObjpt_60 = 0;   
  jetpp[4]->SetBranchAddress("pt",&trgObjpt_40);  jetpp[5]->SetBranchAddress("pt",&trgObjpt_60);  
  std::vector<float> *trgObjpt_80 = 0;   std::vector<float> *trgObjpt_100 = 0;
  jetpp[6]->SetBranchAddress("pt",&trgObjpt_80);  jetpp[7]->SetBranchAddress("pt",&trgObjpt_100);

  // Declare the output File and the necessary histograms after that:
  // std::string outdir="";
  // std::string outfile=outdir+kFoname;
  TFile *fout = new TFile(kFoname.c_str(),"RECREATE");
  fout->cd();

  TH1F * hJetQA[2][21];
  for(int k = 0; k<2; ++k){
    for(int j = 0; j<21; ++j){//woah else if statement, this looks dangerous, dont mess with until i know more
      if(j==2)      hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),Form(";%s;",var[j].c_str()),100, -5, +5  );  
      else if(j==3) hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),Form(";%s;",var[j].c_str()),100, -4, +4  );
      else if(j<=1) hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),Form(";%s;",var[j].c_str()),500,  0, 500 );
      else if(j>=4) hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),Form(";%s;",var[j].c_str()),200,  0, 2   );
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

  hpp_TrgObj100[0]  = new TH1F(Form("hpp_HLT100_noJetID_R%d_%s"     , radius, etaWidth), Form("Spectra from  Jet 100 R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_TrgObj80[0]   = new TH1F(Form("hpp_HLT80_noJetID_R%d_%s"      , radius, etaWidth), Form("Spectra from  Jet 80 R%d %s " ,radius,etaWidth),2000, 0, 2000);
  hpp_TrgObj60[0]   = new TH1F(Form("hpp_HLT60_noJetID_R%d_%s"      , radius, etaWidth), Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_TrgObj40[0]   = new TH1F(Form("hpp_HLT40_noJetID_R%d_%s"      , radius, etaWidth), Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_TrgObjComb[0] = new TH1F(Form("hpp_HLTComb_noJetID_R%d_%s"    , radius, etaWidth), Form("Trig Combined Spectra R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_CombJetpT[0]  = new TH1F(Form("hpp_TrgCombTest_noJetID_R%d_%s", radius, etaWidth), Form("Trig Combined Spectra KurtMethod R%d %s ",radius,etaWidth),2000, 0, 2000);

  hpp_TrgObj100[1]  = new TH1F(Form("hpp_HLT100_JetID_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 100 R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_TrgObj80[1]   = new TH1F(Form("hpp_HLT80_JetID_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_TrgObj60[1]   = new TH1F(Form("hpp_HLT60_JetID_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_TrgObj40[1]   = new TH1F(Form("hpp_HLT40_JetID_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_TrgObjComb[1] = new TH1F(Form("hpp_HLTComb_JetID_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),2000, 0, 2000);
  hpp_CombJetpT[1]  = new TH1F(Form("hpp_TrgCombTest_JetID_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra KurtMethod R%d %s ",radius,etaWidth),2000, 0, 2000);

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

  Long64_t nentries = jetpp[0]->GetEntries();
  //Long64_t nGoodEvt = 0;
  //if(printDebug) nentries = 100;
  Long64_t NEvents = 0;
  Long64_t NEvents_100 = 0;
  Long64_t NEvents_80 = 0;
  Long64_t NEvents_60 = 0;
  Long64_t NEvents_40 = 0;    
  
  TRandom3 rnd; //wonder what this is for...
  //rnd.SetSeed(endfile);

  std::cout<<"Running through all the "<<nentries<<" events now"<<std::endl;  
  for(int nEvt = 0; nEvt < nentries; ++ nEvt) {//event loop

    if(printDebug)std::cout<<"nEvt = "<<nEvt<<std::endl;
    if(printDebug && nEvt==10)break;

    //grab an entry from each ppjet tree
    for(int i = 0; i<N;++i) jetpp[i]->GetEntry(nEvt);
    
    //Evt/skim analysis cuts
    if( pHBHENoiseFilter_F==0     || 
        pBeamScrapingFilter_F==0  || 
        pprimaryvertexFilter_F==0 ||
        fabs(vz_F)>15              ) 
      continue;
    
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

    //trigger decision boolean array w/ prescale l1 seed 
    bool trgDec[4] = {(bool)jet40_F, (bool)jet60_F, (bool)jet80_F, (bool)jet100_F};
    int treePrescl[4] = {(jet40_p_F*jet40_l1seed_p_F), (jet60_p_F*jet60_l1seed_p_F), (jet80_p_F*jet80_l1seed_p_F), (jet100_p_F*jet100_l1seed_p_F)};    
    int maxtrg= -1;      
    for(int ii=4; ii>=0; ii--){
      if(trgDec[ii]==1){
	maxtrg=ii;
	break;
      }
    }
    
    //loop over objects that pass each trigger
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

    //come up with event weight for later use
    double weight_eS = trigComb(trgDec, treePrescl, triggerPt);
    
    // check trigger decisions for events + exclusivity between them
    if(trgDec[3] && triggerPt>=100                ) is100 = true;
    if(trgDec[2] && triggerPt>=80 && triggerPt<100) is80  = true;
    if(trgDec[1] && triggerPt>=60 && triggerPt<80 ) is60  = true;
    if(trgDec[0] && triggerPt>=40 && triggerPt<60 ) is40  = true;

    // count events and events per each trigger
    NEvents++;
    if(is100) NEvents_100++;
    if(is80)  NEvents_80++;
    if(is60)  NEvents_60++;
    if(is40)  NEvents_40++;
    
    //jet loop
    for(int jet = 0; jet<nref_F; ++jet){

      //pt too low or not in barrel, forget it.
      if( fabs(eta_F[jet]) > 2 ||
          pt_F[jet] <= 15       ) continue;     
      
      float recpt = pt_F[jet];
      // float rawpt = rawpt_F[jet];


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

      if(1){//why raghav why
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
	
	if(is40)  hpp_TrgObj40[0]->Fill(recpt, treePrescl[0]);
	if(is60)  hpp_TrgObj60[0]->Fill(recpt, treePrescl[1]);
	if(is80)  hpp_TrgObj80[0]->Fill(recpt, treePrescl[2]);
	if(is100) hpp_TrgObj100[0]->Fill(recpt, treePrescl[3]);

	//looking at lead+subleading jet
	if(jet == 0) {
	  hJetQA[0][19]->Fill( (float) (pt_F[0]-pt_F[1])/(pt_F[0]+pt_F[1]) , weight_eS ); // A_j
	  hJetQA[0][20]->Fill( (float) (pt_F[1]/pt_F[0])                   , weight_eS ); // x_j
	}

	// apply JetID
	if( chSum_F[jet]/recpt > 0    &&    // if  charged hadron pt frac > 0
	    neSum_F[jet]/recpt < 0.99 &&    // and neutral hadron pt frac < 99%
	    phSum_F[jet]/recpt < 0.99 &&    // and photon hadron pt frac < 99%
	    chN_F[jet] > 0             ) {  // and there's a charged particle   
	  
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
	  
	  if(is40)  hpp_TrgObj40[1]->Fill(recpt, treePrescl[0]);
	  if(is60)  hpp_TrgObj60[1]->Fill(recpt, treePrescl[1]);
	  if(is80)  hpp_TrgObj80[1]->Fill(recpt, treePrescl[2]);
	  if(is100) hpp_TrgObj100[1]->Fill(recpt, treePrescl[3]);
	  
	  if(jet == 0){
	    hJetQA[1][19]->Fill((float)(pt_F[0]-pt_F[1])/(pt_F[0]+pt_F[1]), weight_eS);
	    hJetQA[1][20]->Fill((float)(pt_F[1]/pt_F[0]), weight_eS);
	  }  

	}//end jet id cut
      }//end if(1) 

      // if(jet80_F == 1) {      
      // 	hpp_TrgObj80->Fill(recpt);
      // // 	hpp_JEC_TrgObj80->Fill(JEC_upShift);
      // // 	hpp_JEC_minus_TrgObj80->Fill(JEC_downShift);
      // // 	hpp_JEC_gaus_TrgObj80->Fill(JEC_gaussmear);
      // // 	hpp_Smear_TrgObj80->Fill(recpt+rnd.Gaus(0,0.01));
      // // 	hpp_anaBin_TrgObj80->Fill(recpt);
      // }      

    }// end jet loop
  }// end event loop
  
  hpp_TrgObjComb[0]->Add(hpp_TrgObj100[0]);
  hpp_TrgObjComb[0]->Add(hpp_TrgObj80[0] );
  hpp_TrgObjComb[0]->Add(hpp_TrgObj60[0] );
  hpp_TrgObjComb[0]->Add(hpp_TrgObj40[0] );

  hpp_TrgObjComb[1]->Add(hpp_TrgObj100[1]);
  hpp_TrgObjComb[1]->Add(hpp_TrgObj80[1] );
  hpp_TrgObjComb[1]->Add(hpp_TrgObj60[1] );
  hpp_TrgObjComb[1]->Add(hpp_TrgObj40[1] );
  
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
  
  std::cout<<"Event Summary"<<std::endl;
  std::cout<<"Total Number of Events read   = "<<NEvents<<std::endl;
  std::cout<<"Total Number of HLT100 events = "<<NEvents_100<<std::endl;
  std::cout<<"Total Number of HLT80 events = "<<NEvents_80<<std::endl;
  std::cout<<"Total Number of HLT60 events = "<<NEvents_60<<std::endl;
  std::cout<<"Total Number of HLT40 events = "<<NEvents_40<<std::endl;
  
  fout->Write();

  trgObjpt_40->clear();
  trgObjpt_60->clear();
  trgObjpt_80->clear();
  trgObjpt_100->clear();
  
  timer.Stop();

  std::cout<<"Macro finished: "<<std::endl;
  std::cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<std::endl;
  std::cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<std::endl;
  
  return 0;
}// end RAA_read_data_pp macro

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

/////////////////////////
// divide by bin width //
/////////////////////////
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

/////////////////////////
// trigger combination //
/////////////////////////
double trigComb(bool *trg, int *pscl, double pt){

  double weight=0;

  if(trg[3] && pt>=100 )          weight = pscl[3];
  if(trg[2] && pt>=80  && pt<100) weight = pscl[2];
  if(trg[1] && pt>=60  && pt<80 ) weight = pscl[1];
  if(trg[0] && pt>=40  && pt<60 ) weight = pscl[0];
  
  return weight;
}










