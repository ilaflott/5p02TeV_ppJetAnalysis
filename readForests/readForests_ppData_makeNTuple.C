// custom header
#include "readForests_makeNTuple.h"

// ppData switches
const bool useHLT100=false; // leave me off almost always
const bool useMBevts=false; // leave me off until Ian says to use me

//const std::string trgCombType="Calo"; // almost never used
const std::string trgCombType="PF"; // almost always used

const bool deepDebug=false; //put to true only if in trouble...

//lumi in Jet80/LowerJets PDs according to brilcalc, accounting for event selection efficiency [inv picobarns]
const long double effJet80Lumi  =(0.9974435 * 27.385867785) ;
const long double effLowJetsLumi=(0.9984391 * 27.108290680) ;

//// readForests_ppData_makeNTuple
// ---------------------------------------------------------------------------------------------------------------
int readForests_ppData_makeNTuple( std::string inFilelist , int startfile , int endfile , 
				   int radius , std::string jetType , bool debugMode ,
				   std::string outfile, float jtEtaCutLo, float jtEtaCutHi){
  
  // for monitoring performance + debugging
  TStopwatch timer;  timer.Start();
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl;
  
  // basic info the screen
  std::cout<<std::endl<<"///////////////////"<<std::endl;
  std::cout<<"reading filelist "<< inFilelist<<std::endl;
  std::cout<<"reading files #'s "<< startfile << " to " << endfile<<std::endl;
  std::cout<<"radius = " << radius;
  std::cout<<", jetType = " << jetType;
  std::cout<<", jtEtaCutLo = "<<jtEtaCutLo<<std::endl;
  std::cout<<", jtEtaCutHi = "<<jtEtaCutHi<<std::endl;
  std::cout<<", jetQAPtCut = "<<jetQAPtCut<<std::endl;
  std::cout<<", jtPtCut_Hi = "<<jtPtCut_Hi<<std::endl;  
  std::cout<<"///////////////////"<<std::endl<<std::endl;

  // plot settings for later
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  

  //check out what kind of filelist we're using
  bool filelistIsLowerJets=false; //to know if i should apply the duplicate check or not
  bool filelistIsJet80=false;  
  bool filelistIsMinBias=false;
  if( inFilelist.find("LowerJets")!=std::string::npos ) {
    filelistIsLowerJets=true; 
    std::cout<<"running LowerJets filelist!"<<std::endl; }
  else if( inFilelist.find("Jet80")    !=std::string::npos ) {
    filelistIsJet80=true;
    std::cout<<"running Jet80 filelist!"<<std::endl; }
  else if( inFilelist.find("MinBias")    !=std::string::npos ) {//minbias
    filelistIsMinBias=true;
    std::cout<<"running MinBias filelist!"<<std::endl; }
  else assert(false);
  assert( !(filelistIsJet80&&filelistIsLowerJets) );
  assert( !(filelistIsJet80&&filelistIsMinBias) );//minbias
  assert( !(filelistIsMinBias&&filelistIsLowerJets) );//minbias
  
  // jet tree name+directory formation
  std::string jetTreeName="ak"+std::to_string(radius)+jetType+"JetAnalyzer/t";
  if(debugMode)std::cout<<"looking at jetTree "<<jetTreeName<<std::endl;

  // initialize tree name array
  std::string trees[N_dataTrees];
  trees[0]=jetTreeName;
  for(int i=1;i<N_dataTrees;++i){
    if(i<6){
      trees[i]=dataTreeNames[i];
      if(debugMode){
	if(i==5)std::cout<<"trees[i="<<i<<"]="<<trees[i]<<" [NOT BEING USED THIS JOB]"<<std::endl;  
	else std::cout<<"trees[i="<<i<<"]="<<trees[i]<<std::endl;  
      }
    }
    else{
      if(trgCombType=="PF"){
	trees[i]=dataTreeNames[i]+PF_HLTBitStrings[i-6]+"_v";
	if(debugMode)std::cout<<"trees[i="<<i<<"]="<<trees[i]<<std::endl;  
      }
      else if(trgCombType=="Calo"){
	trees[i]=dataTreeNames[i]+Calo_HLTBitStrings[i-6]+"_v";
	if(debugMode)std::cout<<"trees[i="<<i<<"]="<<trees[i]<<std::endl;  
      }
      else assert(false);      
    }
  }
  //if(debugMode)for(int i=0;i<N_dataTrees;++i)
  //assert(false);

  // declare TChains for each tree + friend them
  TChain* jetpp[N_dataTrees];
  for(int t = 0;t<N_dataTrees;++t)  { 
    if(t==5)continue;
    jetpp[t] = new TChain( trees[t].data() );
    if(t>0){
      jetpp[0]->AddFriend( jetpp[t] );      
    }    
  }
  
  // open filelist, add files, including startfile+endfile, to chains
  std::cout<<"opening filelist: "<<inFilelist<<std::endl;
  std::ifstream instr_Forest(inFilelist.c_str(),std::ifstream::in);
  std::string filename_Forest;  
  std::string lastFileAdded=""; bool filesAdded=false; 
  for(int ifile = 0; ifile<=endfile; ++ifile){
    
    // grab a filename, check startfile and end of filelist condition
    instr_Forest>>filename_Forest; 
    if(ifile<startfile){ lastFileAdded=filename_Forest; 
      continue; }
    if(filename_Forest==lastFileAdded){ std::cout<<"end of filelist!"<<std::endl; 
      break; }
    
    std::cout<<"adding file #"<<ifile; 
    if(debugMode)std::cout<<", "<<filename_Forest; 
    
    std::cout<<", to each TChain in array"<<std::endl;
    for(int t = 0;t<N_dataTrees;++t) {
      if( t==5)continue;
      filesAdded=jetpp[t]->AddFile(filename_Forest.c_str()); 
      if(debugMode)
	jetpp[t]->TChain::Lookup(true);
    }
    lastFileAdded=filename_Forest;  
  }
  assert(filesAdded);//avoid segfault later
  
  
  //variables for writing output to TTree
  int run, lumi, evt;
  int njets;
  bool lowjetsPD, jet80PD;
  bool l1jet40, l1jet60, l1jet80;//L1 seed
  bool jet40, jet60, jet80;//L1seed * HLT decision
  int jet40_ps, jet60_ps, jet80_ps;//L1 prescale * HLT prescale
  bool is40, is60, is80;//final HLT decision after exclusion
  float evtweight;//final weight, including effective luminosity
  std::vector<float> rawpt, jetpt, jeteta, jetphi;
  std::vector<bool> jetIDtight, jetIDloose;
  
  // Declare the output file, declare TTree+branches after for easy write
  std::cout<<"opening output file "<<outfile<<std::endl;
  TFile *fout = new TFile(outfile.c_str(),"RECREATE");
  fout->cd();
  
  
  

  TTree* outputTree=new TTree("outTree","outTree");
  outputTree->Branch("run", &run, "run/I");
  outputTree->Branch("lumi", &lumi, "lumi/I");
  outputTree->Branch("evt", &evt, "evt/I");

  outputTree->Branch("njets", &njets, "njets/I");

  outputTree->Branch("lowjetsPD", &lowjetsPD, "lowjetsPD/O");
  outputTree->Branch("jet80PD", &jet80PD, "jet80PD/O");

  outputTree->Branch("l1jet40", &l1jet40, "l1jet40/O");
  outputTree->Branch("l1jet60", &l1jet60, "l1jet60/O");
  outputTree->Branch("l1jet80", &l1jet80, "l1jet80/O");

  outputTree->Branch("jet40", &jet40, "jet40/O");
  outputTree->Branch("jet60", &jet60, "jet60/O");
  outputTree->Branch("jet80", &jet80, "jet80/O");
  
  outputTree->Branch("jet40_ps", &jet40_ps, "jet40_ps/I");
  outputTree->Branch("jet60_ps", &jet60_ps, "jet60_ps/I");
  outputTree->Branch("jet80_ps", &jet80_ps, "jet80_ps/I");

  outputTree->Branch("is40", &is40, "is40/O");
  outputTree->Branch("is60", &is60, "is60/O");
  outputTree->Branch("is80", &is80, "is80/O");

  outputTree->Branch("evtweight",&evtweight,"evtweight/F");

  outputTree->Branch("rawpt" ,&rawpt);
  outputTree->Branch("jetpt" ,&jetpt);
  outputTree->Branch("jeteta",&jeteta);
  outputTree->Branch("jetphi",&jetphi);
  
  outputTree->Branch("jetIDtight",&jetIDtight);
  outputTree->Branch("jetIDloose",&jetIDloose);

  // declare cut hists
  TH1F *hJetPtCut        =new TH1F("hJetPtCut"      ,(std::to_string(jetQAPtCut)).c_str()   ,    100, 0.,100.); hJetPtCut->Fill(jetQAPtCut);           
  TH1F *hJetPtCut_Hi        =new TH1F("hJetPtCut_Hi"      ,(std::to_string(jtPtCut_Hi)).c_str()   ,    1000,1500.,2500.); hJetPtCut_Hi->Fill(jtPtCut_Hi);           
  TH1F *hJetEtaCutHi       =new TH1F("hJetEtaCutHi"     ,(std::to_string(jtEtaCutHi)).c_str()  ,    60,   0,6.  ); hJetEtaCutHi->Fill(jtEtaCutHi); 
  TH1F *hJetEtaCutLo       =new TH1F("hJetEtaCutLo"     ,(std::to_string(jtEtaCutLo)).c_str()  ,    60,   0,6.  ); hJetEtaCutLo->Fill(jtEtaCutLo);     
  
  TH1F *hHLT40ptThresh       =new TH1F("hHLT40ptThresh"      ,(std::to_string(HLTPFthresh[0])).c_str()   ,    200, 0.,200.); hHLT40ptThresh->Fill(HLTPFthresh[0]); 
  TH1F *hHLT60ptThresh       =new TH1F("hHLT60ptThresh"      ,(std::to_string(HLTPFthresh[1])).c_str()   ,    200, 0.,200.); hHLT60ptThresh->Fill(HLTPFthresh[1]); 
  TH1F *hHLT80ptThresh       =new TH1F("hHLT80ptThresh"      ,(std::to_string(HLTPFthresh[2])).c_str()   ,    200, 0.,200.); hHLT80ptThresh->Fill(HLTPFthresh[2]);   
  
  // EVENT LOOP + PREP
  // declare variables/arrays + set branch address for each input tree
  //JetAnalyzer, jets
  int nref_I;  
  float pt_F[1000];  float rawpt_F[1000];
  float eta_F[1000];   float phi_F[1000];  
  
  //tracks
  int trkN_I[1000], trkHardN_I[1000]; 
  float trkSum_F[1000],  trkMax_F[1000], trkHardSum_F[1000];

  //charged hadrons
  int chN_I[1000], chHardN_I[1000]; 
  float chSum_F[1000],  chMax_F[1000], chHardSum_F[1000];
  
  //photons  
  int phN_I[1000], phHardN_I[1000]; 
  float phSum_F[1000],  phMax_F[1000], phHardSum_F[1000];

  //leptons
  int eN_I[1000], muN_I[1000];
  float eSum_F[1000];   float eMax_F[1000];        
  float muSum_F[1000];  float muMax_F[1000];    

  //neutrals
  int neN_I[1000];
  float neSum_F[1000];  float neMax_F[1000];        

  //jets
  jetpp[0]->SetBranchAddress("nref",&nref_I);
  jetpp[0]->SetBranchAddress("jtpt",&pt_F);
  jetpp[0]->SetBranchAddress("rawpt",&rawpt_F);
  jetpp[0]->SetBranchAddress("jteta",&eta_F);
  jetpp[0]->SetBranchAddress("jtphi",&phi_F);
  
  //tracks
  jetpp[0]->SetBranchAddress("trackN",&trkN_I);  
  jetpp[0]->SetBranchAddress("trackSum",&trkSum_F);  
  jetpp[0]->SetBranchAddress("trackMax",&trkMax_F);  
  jetpp[0]->SetBranchAddress("trackHardN",&trkHardN_I);  
  jetpp[0]->SetBranchAddress("trackHardSum",&trkHardSum_F);
  
  //charged hadrons from PF (no HF-Calo ctontributions)
  jetpp[0]->SetBranchAddress("chargedN",&chN_I);  
  jetpp[0]->SetBranchAddress("chargedSum",&chSum_F);
  jetpp[0]->SetBranchAddress("chargedMax",&chMax_F);  
  jetpp[0]->SetBranchAddress("chargedHardN",&chHardN_I);  
  jetpp[0]->SetBranchAddress("chargedHardSum",&chHardSum_F);
  
  //photons from PF (no HF-Calo contributions)
  jetpp[0]->SetBranchAddress("photonN",&phN_I);  
  jetpp[0]->SetBranchAddress("photonSum",&phSum_F);  
  jetpp[0]->SetBranchAddress("photonMax",&phMax_F);  
  jetpp[0]->SetBranchAddress("photonHardN",&phHardN_I);  
  jetpp[0]->SetBranchAddress("photonHardSum",&phHardSum_F);
  
  //electrons from PF
  jetpp[0]->SetBranchAddress("eN",&eN_I);  
  jetpp[0]->SetBranchAddress("eSum",&eSum_F);  
  jetpp[0]->SetBranchAddress("eMax",&eMax_F);
  
  //muons from PF
  jetpp[0]->SetBranchAddress("muN",&muN_I);
  jetpp[0]->SetBranchAddress("muSum",&muSum_F);
  jetpp[0]->SetBranchAddress("muMax",&muMax_F);
  
  //neutral hadrons from PF (no HF e/gamma contributions)
  jetpp[0]->SetBranchAddress("neutralN",&neN_I);  
  jetpp[0]->SetBranchAddress("neutralSum",&neSum_F); 
  jetpp[0]->SetBranchAddress("neutralMax",&neMax_F);
  
  
  // HiEvtAnalyzer
  float vx_F,vy_F,vz_F;
  ULong64_t evt_I;   UInt_t run_I;   UInt_t lumi_I;   
  jetpp[1]->SetBranchAddress("vz",&vz_F);
  jetpp[1]->SetBranchAddress("vy",&vy_F);
  jetpp[1]->SetBranchAddress("vx",&vx_F);
  jetpp[1]->SetBranchAddress("evt",&evt_I);
  jetpp[1]->SetBranchAddress("run",&run_I);
  jetpp[1]->SetBranchAddress("lumi",&lumi_I);


  // skimanalysis
  int pBeamScrapingFilter_I, pHBHENoiseFilter_I, pprimaryvertexFilter_I,puvertexFilter_I;
  jetpp[2]->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_I);
  jetpp[2]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilter_I);
  jetpp[2]->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryvertexFilter_I);
  jetpp[2]->SetBranchAddress("pVertexFilterCutGtight",&puvertexFilter_I);
  

  // ppTrack/trackTree
  //vtxs
  int nVtx_I;
  int nTrkVtx_I[1000];
  float xVtx_F[1000];
  float yVtx_F[1000];
  float zVtx_F[1000];

  jetpp[3]->SetBranchAddress("nVtx",&nVtx_I);
  jetpp[3]->SetBranchAddress("nTrkVtx",&nTrkVtx_I);
  jetpp[3]->SetBranchAddress("xVtx",&xVtx_F);
  jetpp[3]->SetBranchAddress("yVtx",&yVtx_F);
  jetpp[3]->SetBranchAddress("zVtx",&zVtx_F);  
  
  //trks
  int     nTrk_I[1000];
  float  trkPt_F[10000];
  float trkEta_F[10000];
  float trkPhi_F[10000];

  jetpp[3]->SetBranchAddress("nTrk"   ,   &nTrk_I  );
  jetpp[3]->SetBranchAddress("trkPt"   ,   &trkPt_F  );
  jetpp[3]->SetBranchAddress("trkEta"   ,   &trkEta_F  );
  jetpp[3]->SetBranchAddress("trkPhi"   ,   &trkPhi_F  );
  //jetpp[3]->SetBranchAddress(""   ,   &  );  

  // hltanalysis
  // L1 Trigger Branches
  std::string L1Branches[N_L1Bits], L1PresclBranches[N_L1Bits];
  for(int i=0; i<N_L1Bits; ++i) {
    L1Branches[i] = L1BitStrings[i];
    L1PresclBranches[i] = L1BitStrings[i] + "_Prescl";  
    //if(debugMode){std::cout<<std::endl;
    //  std::cout<<"L1Branches      [i="<<i<<"]= "<<L1Branches[i]<<std::endl;
    //  std::cout<<"L1PresclBranches      [i="<<i<<"]= "<<L1PresclBranches[i]<<std::endl;    }    
  }




  //L1 seed used for HLT seeding HLTMB
  int mb_l1s_I;
  int mb_l1s_ps_I;
  jetpp[4]->SetBranchAddress( "L1_MinimumBiasHF1_OR", &mb_l1s_I);
  jetpp[4]->SetBranchAddress( "L1_MinimumBiasHF1_OR_Prescl", &mb_l1s_ps_I);
  
  //technical bit to perhaps be used for cleaning up min bias signal
  int L1_HF_coincidence_I;
  jetpp[4]->SetBranchAddress("L1Tech_HCAL_HF_coincidence_PM.v2",&L1_HF_coincidence_I);

  //other L1 Jets not used for seeding HLT40,60,80,100
  int l1jet24_l1s_I, l1jet20_l1s_I, l1jet16_l1s_I, l1jet12_l1s_I, l1jet8_l1s_I ;
  jetpp[4]->SetBranchAddress (  "L1_SingleJet8_BptxAND", & l1jet8_l1s_I  );
  jetpp[4]->SetBranchAddress (  "L1_SingleJet12_BptxAND", & l1jet12_l1s_I );
  jetpp[4]->SetBranchAddress (  "L1_SingleJet16_BptxAND", & l1jet16_l1s_I );
  jetpp[4]->SetBranchAddress (  "L1_SingleJet20_BptxAND", & l1jet20_l1s_I );
  jetpp[4]->SetBranchAddress (  "L1_SingleJet24_BptxAND", & l1jet24_l1s_I );

  //L1 Jets used for HLT seeding HLT40,60,80,100
  int jet40_l1s_I, jet60_l1s_I,  jet80_l1s_I,  jet100_l1s_I;   //L1 
  int jet40_l1s_ps_I, jet60_l1s_ps_I,  jet80_l1s_ps_I,  jet100_l1s_ps_I;     
  jetpp[4]->SetBranchAddress( L1Branches[0].c_str()  , &jet40_l1s_I);   //l1jet28
  jetpp[4]->SetBranchAddress( L1Branches[1].c_str()  , &jet60_l1s_I);   //l1jet40
  jetpp[4]->SetBranchAddress( L1Branches[2].c_str()  , &jet80_l1s_I);   //l1jet48
  jetpp[4]->SetBranchAddress( L1Branches[3].c_str()  , &jet100_l1s_I);  //l1jet52
  
  jetpp[4]->SetBranchAddress( L1PresclBranches[0].c_str()  , &jet40_l1s_ps_I);
  jetpp[4]->SetBranchAddress( L1PresclBranches[1].c_str()  , &jet60_l1s_ps_I);
  jetpp[4]->SetBranchAddress( L1PresclBranches[2].c_str()  , &jet80_l1s_ps_I);
  jetpp[4]->SetBranchAddress( L1PresclBranches[3].c_str()  , &jet100_l1s_ps_I);
  
  int l1jet60_l1s_I, l1jet68_l1s_I;
  jetpp[4]->SetBranchAddress( "L1_SingleJet60_BptxAND", &l1jet60_l1s_I);
  jetpp[4]->SetBranchAddress( "L1_SingleJet68_BptxAND", &l1jet68_l1s_I);




 
  
  
  // make sure # of L1/HLT bits is correct
  assert(N_HLTBits==N_L1Bits);  


  int MB_HF1ORp5_I, MB_HF1ORp5_p_I;
  jetpp[4]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part5_v1"        , &MB_HF1ORp5_I);
  jetpp[4]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part5_v1_Prescl" , &MB_HF1ORp5_p_I);
  
  
  std::string PF_HLTBranches[N_HLTBits],PF_HLTPresclBranches[N_HLTBits];
  std::string Calo_HLTBranches[N_HLTBits], Calo_HLTPresclBranches[N_HLTBits];
  double HLTthresh[N_HLTBits]={-1.};
  for(int i=0; i<N_HLTBits; ++i) {
    PF_HLTBranches[i] = PF_HLTBitStrings[i] + "_v1"; 
    PF_HLTPresclBranches[i] = PF_HLTBitStrings[i] + "_v1_Prescl";   
    Calo_HLTBranches[i] = Calo_HLTBitStrings[i] + "_v1"; 
    Calo_HLTPresclBranches[i] = Calo_HLTBitStrings[i] + "_v1_Prescl";     
    if(trgCombType=="PF") HLTthresh[i]=HLTPFthresh[i];	
    else if (trgCombType=="Calo")HLTthresh[i]=HLTCalothresh[i];
    else assert(false);
  }
  
  int Jet40_I, Jet60_I, Jet80_I, Jet100_I;   //primary HLT bits, either Calo or PF
  int Jet40_p_I, Jet60_p_I, Jet80_p_I, Jet100_p_I;    //prescales

  int Jet40_2_I, Jet60_2_I, Jet80_2_I, Jet100_2_I; //the other 40/60 bits from whichever one wasn't picked (for duplicate skip)
  if(trgCombType=="PF"){
    jetpp[4]->SetBranchAddress( PF_HLTBranches[0].c_str() , &Jet40_I);
    jetpp[4]->SetBranchAddress( PF_HLTBranches[1].c_str() , &Jet60_I);
    jetpp[4]->SetBranchAddress( PF_HLTBranches[2].c_str() , &Jet80_I);
    jetpp[4]->SetBranchAddress( PF_HLTBranches[3].c_str() , &Jet100_I);    
    jetpp[4]->SetBranchAddress( PF_HLTPresclBranches[0].c_str() , &Jet40_p_I);
    jetpp[4]->SetBranchAddress( PF_HLTPresclBranches[1].c_str() , &Jet60_p_I);
    jetpp[4]->SetBranchAddress( PF_HLTPresclBranches[2].c_str() , &Jet80_p_I);
    jetpp[4]->SetBranchAddress( PF_HLTPresclBranches[3].c_str() , &Jet100_p_I);
    jetpp[4]->SetBranchAddress( Calo_HLTBranches[0].c_str() , &Jet40_2_I);
    jetpp[4]->SetBranchAddress( Calo_HLTBranches[1].c_str() , &Jet60_2_I);  
    jetpp[4]->SetBranchAddress( Calo_HLTBranches[2].c_str() , &Jet80_2_I);
    jetpp[4]->SetBranchAddress( Calo_HLTBranches[3].c_str() , &Jet100_2_I);  }
  else if(trgCombType=="Calo"){
    jetpp[4]->SetBranchAddress( Calo_HLTBranches[0].c_str() , &Jet40_I);
    jetpp[4]->SetBranchAddress( Calo_HLTBranches[1].c_str() , &Jet60_I);
    jetpp[4]->SetBranchAddress( Calo_HLTBranches[2].c_str() , &Jet80_I);
    jetpp[4]->SetBranchAddress( Calo_HLTBranches[3].c_str() , &Jet100_I);    
    jetpp[4]->SetBranchAddress( Calo_HLTPresclBranches[0].c_str() , &Jet40_p_I);
    jetpp[4]->SetBranchAddress( Calo_HLTPresclBranches[1].c_str() , &Jet60_p_I);
    jetpp[4]->SetBranchAddress( Calo_HLTPresclBranches[2].c_str() , &Jet80_p_I);
    jetpp[4]->SetBranchAddress( Calo_HLTPresclBranches[3].c_str() , &Jet100_p_I);
    jetpp[4]->SetBranchAddress( PF_HLTBranches[0].c_str() , &Jet40_2_I);
    jetpp[4]->SetBranchAddress( PF_HLTBranches[1].c_str() , &Jet60_2_I);
    jetpp[4]->SetBranchAddress( PF_HLTBranches[2].c_str() , &Jet80_2_I);
    jetpp[4]->SetBranchAddress( PF_HLTBranches[3].c_str() , &Jet100_2_I);  }
  else assert(false);
  
  // other bits that i need to keep track of for full proper duplicate skip over jet80 events in the min bias dataset
  // pretty sure everyone of these bits are redundant w/ jet80 bits i.e. if they fire, so did the jet80 bits.
  // but i'm having issues... so i suspect no one and i suspect everyone.
  int CaloJet80_45_45_Eta2p1_I,  CaloJet80_Jet35_Eta0p7_I,  CaloJet80_Jet35_Eta1p1_I;
  int CaloJet100_Jet35_Eta0p7_I,  CaloJet100_Jet35_Eta1p1_I;
  int CaloJet110_I,  CaloJet120_I,  CaloJet150_I;
  jetpp[4]->SetBranchAddress("HLT_AK4CaloJet80_45_45_Eta2p1_v1", &CaloJet80_45_45_Eta2p1_I );
  jetpp[4]->SetBranchAddress("HLT_AK4CaloJet80_Jet35_Eta0p7_v1", &CaloJet80_Jet35_Eta0p7_I );
  jetpp[4]->SetBranchAddress("HLT_AK4CaloJet80_Jet35_Eta1p1_v1", &CaloJet80_Jet35_Eta1p1_I );
  jetpp[4]->SetBranchAddress("HLT_AK4CaloJet100_Jet35_Eta0p7_v1", &CaloJet100_Jet35_Eta0p7_I );
  jetpp[4]->SetBranchAddress("HLT_AK4CaloJet100_Jet35_Eta1p1_v1", &CaloJet100_Jet35_Eta1p1_I );
  jetpp[4]->SetBranchAddress("HLT_AK4CaloJet110_Eta5p1_v1", &CaloJet110_I );
  jetpp[4]->SetBranchAddress("HLT_AK4CaloJet120_Eta5p1_v1", &CaloJet120_I );
  jetpp[4]->SetBranchAddress("HLT_AK4CaloJet150_v1", &CaloJet150_I );
  
  int PFJet120_I, PFJet110_I;
  jetpp[4]->SetBranchAddress("HLT_AK4PFJet110_Eta5p1_v1", &PFJet110_I);
  jetpp[4]->SetBranchAddress("HLT_AK4PFJet120_Eta5p1_v1", &PFJet120_I);  





  std::vector<double>  *trgObjpt_40, *trgObjpt_60, *trgObjpt_80, *trgObjpt_100;
  jetpp[6]->SetBranchAddress("pt",&trgObjpt_40);
  jetpp[7]->SetBranchAddress("pt",&trgObjpt_60);  
  jetpp[8]->SetBranchAddress("pt",&trgObjpt_80);  
  jetpp[9]->SetBranchAddress("pt",&trgObjpt_100);
  
  std::vector<double>  *trgObjeta_40, *trgObjeta_60, *trgObjeta_80, *trgObjeta_100;
  jetpp[6]->SetBranchAddress("eta",&trgObjeta_40);
  jetpp[7]->SetBranchAddress("eta",&trgObjeta_60);  
  jetpp[8]->SetBranchAddress("eta",&trgObjeta_80);  
  jetpp[9]->SetBranchAddress("eta",&trgObjeta_100);
  

  
  //evtcounts check
  UInt_t NEvents_jetAnalyzr=jetpp[0]->GetEntries();   // preskim event count from files
  UInt_t NEvents_skimAnalyzr=jetpp[3]->GetEntries();   // preskim event count from files
  bool jetEventsExist= (bool)NEvents_jetAnalyzr, skimEventsExist=(bool)NEvents_skimAnalyzr;
  bool EventCountsEqual= (NEvents_jetAnalyzr==NEvents_skimAnalyzr);
  if( !EventCountsEqual || 
      !jetEventsExist || !skimEventsExist ){//if this fires, we're gonna segfault
    std::cout<<"jetAnalyzr # evts  = " << NEvents_jetAnalyzr  << std::endl;
    std::cout<<"skimAnalyzr # evts = " << NEvents_skimAnalyzr << std::endl; }
  assert( jetEventsExist && skimEventsExist && EventCountsEqual );
  
  UInt_t NEvents_allFiles=NEvents_jetAnalyzr;   // preskim event count from files  
  UInt_t NEvents_toRead= NEvents_allFiles;
  if(debugMode)NEvents_toRead/=5.;
  std::cout<<"reading "<<NEvents_toRead<<" out of "<< NEvents_allFiles<<" events"<<std::endl;  

  





  
  //begin event loop
  for(UInt_t nEvt = 0; nEvt < NEvents_toRead; ++nEvt) {
    

    run=-1, lumi=-1, evt=-1;
    njets=0;
    lowjetsPD=false,    jet80PD=false;
    l1jet40=false, l1jet60=false, l1jet80=false;
    jet40=false, jet60=false, jet80=false;
    jet40_ps=0, jet60_ps=0, jet80_ps=0;
    is40=false, is60=false, is80=false;
    evtweight=0.;
    
    rawpt.clear();
    jetpt.clear();
    jeteta.clear();
    jetphi.clear();
    jetIDtight.clear();
    jetIDloose.clear();

    float weight_eS=0.;
    if (nEvt%100000==0 || deepDebug) {
      //if (true) {
      float nEvt_fracdone=((float)nEvt)/((float)NEvents_toRead);
      float nEvt_percdone=nEvt_fracdone*100.;
      int   nEvt_intpercdone=(int)nEvt_percdone;                
      std::cout<<"from trees, grabbing Evt # = "<<nEvt<<", ~"<<nEvt_intpercdone<<"% of job completed"<<std::endl;
    }
    
    //grab an entry    
    jetpp[0]->GetEntry(nEvt);
    
    
    // skim/HiEvtAnalysis criteria
    if( (pHBHENoiseFilter_I    ==0    || 
	 pBeamScrapingFilter_I ==0 || 
	 pprimaryvertexFilter_I==0   ) )  continue;  
    
    if( (fabs(vz_F)>24. ) )     continue;
    
    
    bool trgDec[N_HLTBits]={ (bool) (Jet40_I * jet40_l1s_I) , (bool) ( Jet60_I * jet60_l1s_I ) ,
			     (bool) (Jet80_I * jet80_l1s_I) , (bool) (Jet100_I * jet100_l1s_I) };    
    int trgPscl[N_HLTBits]={ Jet40_p_I * jet40_l1s_ps_I ,  Jet60_p_I *  jet60_l1s_ps_I ,   
                             Jet80_p_I * jet80_l1s_ps_I , Jet100_p_I * jet100_l1s_ps_I };    
    
    run=run_I;
    lumi=lumi_I;
    evt=evt_I;
    
    lowjetsPD=filelistIsLowerJets;
    jet80PD=filelistIsJet80;
    
    l1jet40=(bool)jet40_l1s_I;
    l1jet60=(bool)jet60_l1s_I;
    l1jet80=(bool)jet80_l1s_I;
    
    jet40=trgDec[0];
    jet60=trgDec[1];
    jet80=trgDec[2];

    jet40_ps=trgPscl[0];
    jet60_ps=trgPscl[1];
    jet80_ps=trgPscl[2];
    

    
    double HLT40maxTrgPt=-1.,  HLT60maxTrgPt=-1.,  HLT80maxTrgPt=-1.;
    unsigned int trgObj40_size=trgObjpt_40->size(), trgObj60_size=trgObjpt_60->size();
    unsigned int trgObj80_size=trgObjpt_80->size();
    if(trgDec[2]){
      for(unsigned int itt=0; itt<trgObj80_size; ++itt){	
	double trgpt=trgObjpt_80->at(itt);
	if(trgpt > HLT80maxTrgPt)
	  HLT80maxTrgPt=trgpt;
      }
    }    
    if(trgDec[1]) {
      for(unsigned int itt=0; itt<trgObj60_size; ++itt){
	double trgpt=trgObjpt_60->at(itt);
	if(trgpt > HLT60maxTrgPt)
	  HLT60maxTrgPt=trgpt;
      }    
    }    
    if(trgDec[0])  {
      for(unsigned int itt=0; itt<trgObj40_size; ++itt)	{
	double trgpt=trgObjpt_40->at(itt);
	if(trgpt > HLT40maxTrgPt)
	  HLT40maxTrgPt=trgpt;
      }
    }
    
    
    
    //jet trigger exclusion
    if( trgDec[2] && !(HLT80maxTrgPt<HLTthresh[2])  ) 
      { is80  = true;  weight_eS=(float)trgPscl[2];}
    else if( trgDec[1] && !(HLT60maxTrgPt<HLTthresh[1])  && HLT60maxTrgPt<HLTthresh[2]  ) 
      { is60  = true;  weight_eS=(float)trgPscl[1];}
    else if( trgDec[0] && !(HLT40maxTrgPt<HLTthresh[0])  && HLT40maxTrgPt<HLTthresh[1]  ) 
      { is40  = true;  weight_eS=(float)trgPscl[0];}
    
    
    if(!is40 && !is60 && !is80)continue;      // if it doesnt fire any trigger exclusively
    if(lowjetsPD && !(is40 || is60))continue; // if lowjetsPD, only is40 and is60 decisions allowed through
    if(jet80PD && !is80)continue;             // if jet80 PD, only is80 decisions allowed through

    if(filelistIsJet80)
      weight_eS/=effJet80Lumi;
    else if(filelistIsLowerJets)
      weight_eS/=effLowJetsLumi;
    evtweight=effLowJetsLumi;

    int jetsPerEvent=0;
    for(int jet = 0; jet<nref_I; ++jet){
      
      
      
      
      // event+jet counting
      if( !(pt_F[jet] > jetQAPtCut) ){  continue;} //paranoia about float comparisons
      if( !(pt_F[jet] < jtPtCut_Hi) ){  continue;}      
      if( !(fabs(eta_F[jet]) < 4.7) ) {  continue;}      //this prett redundant; PDs dont seem to write jets past 4.7
      if( fabs(eta_F[jet]) < jtEtaCutLo ) { 	continue;}
      else if( !(fabs(eta_F[jet]) < jtEtaCutHi) ) { 	continue;}      
      
      int chMult  = chN_I[jet] + eN_I[jet] + muN_I[jet] ;
      int neuMult = neN_I[jet] + phN_I[jet] ;
      int numConst  = chMult + neuMult;
      
      // 13 TeV JetID criterion, loose or tight
      bool passesLooseJetID=false;
      if (!(fabs(eta_F[jet]) > 2.4)) 
	passesLooseJetID=(bool)jetID_00eta24( pt_F[jet], 
					      neSum_F[jet],  phSum_F[jet],  chSum_F[jet],  eSum_F[jet], muSum_F[jet],
					      numConst,  chMult);
      else if ( !(fabs(eta_F[jet])>2.7) && fabs(eta_F[jet])>2.4 ) 
	passesLooseJetID=(bool) jetID_24eta27( pt_F[jet],
					       neSum_F[jet],  phSum_F[jet], muSum_F[jet],
					       numConst);
      else if( !(fabs(eta_F[jet])>3.0) && fabs(eta_F[jet])>2.7 )
	passesLooseJetID=(bool) jetID_27eta30( pt_F[jet],
					       neSum_F[jet],  phSum_F[jet], 
					       neuMult);
      else  
	passesLooseJetID=(bool)jetID_32eta47( pt_F[jet], 
					      phSum_F[jet]);
      
      
      if(!passesLooseJetID){
	 continue;}
      
      bool passesTightJetID=false;
      if (!(fabs(eta_F[jet]) > 2.4)) 
	passesTightJetID=(bool)jetID_00eta24( pt_F[jet], 
					      neSum_F[jet],  phSum_F[jet],  chSum_F[jet],  eSum_F[jet], muSum_F[jet],
					      numConst,  chMult,true);
      else if ( !(fabs(eta_F[jet])>2.7) && fabs(eta_F[jet])>2.4 ) 
	passesTightJetID=(bool) jetID_24eta27( pt_F[jet],
					  neSum_F[jet],  phSum_F[jet], muSum_F[jet],
					       numConst,true);
      else if( !(fabs(eta_F[jet])>3.0) && fabs(eta_F[jet])>2.7 )
	passesTightJetID=(bool) jetID_27eta30( pt_F[jet],
					       neSum_F[jet],  phSum_F[jet], 
					       neuMult,true);
      else  
	passesTightJetID=(bool)jetID_32eta47( pt_F[jet], 
					      phSum_F[jet],true);
      
      rawpt.push_back(rawpt_F[jet]);
      jetpt.push_back(pt_F[jet]);
      jeteta.push_back(eta_F[jet]);
      jetphi.push_back(phi_F[jet]);
      jetIDtight.push_back(passesTightJetID);
      jetIDloose.push_back(passesLooseJetID);
      jetsPerEvent++;
      
    }//end incjetanalyzer jet loop

    njets=jetsPerEvent;    
    if(njets==0)continue;//if no jets in target region, don't fill the ntuple
    
    outputTree->Fill();
  }//end event loop
  
  std::cout<<std::endl<<"writing output file... "<<std::endl;
  fout->Write(); 
  //assert(false);

  
  std::cout<<"readForests_ppData_makeNTuple finished."<<std::endl;  timer.Stop();
  std::cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<std::endl;
  std::cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<std::endl;
  std::cout<<"output file written is  "<<outfile<<std::endl<<std::endl;
  return 0;
} 


////// main //////
int main(int argc, char *argv[]){
  
  // error, not enough arguments
  int rStatus = -1;
  if(argc!=8 && argc!=1 && argc!=10){
    std::cout<<std::endl<<std::endl;
    std::cout<<"for actually running, do..."<<std::endl;
    std::cout<<"./readForests_ppData_makeNTuple.exe ";
    std::cout<<"<inputFileList> <startFile> <endFile> ";
    std::cout<<"<jetRadius> <jetType> <debugMode> ";
    std::cout<<"<outputFilename> [<absEtaCutLo> <absEtaCutHi>]";
    std::cout<<std::endl<<std::endl;
    std::cout<<"where args in [] are optional."<<std::endl;
    std::cout<<std::endl<<std::endl;
    std::cout<<"rStatus="<<rStatus<<std::endl;
    std::cout<<std::endl<<std::endl;
    return rStatus;
  }
  
  // good input, run
  rStatus=1;
  if(argc==1) rStatus = readForests_ppData_makeNTuple();
  else if (argc==8){//read input argument vector
    std::string inputFileList=argv[1]; int startfile= atoi(argv[2]); int endfile= atoi(argv[3]);  
    int jetRadius= atoi(argv[4]); std::string jetType=argv[5];     bool debug=(bool)atoi(argv[6]);
    std::string outputFileName=argv[7];      
    rStatus = readForests_ppData_makeNTuple( inputFileList, startfile, endfile, 
				  jetRadius, jetType, debug,
				  outputFileName);
  }
  else if (argc==10){//read input argument vector
    std::string inputFileList=argv[1]; int startfile= atoi(argv[2]); int endfile= atoi(argv[3]);  
    int jetRadius= atoi(argv[4]); std::string jetType=argv[5];     bool debug=(bool)atoi(argv[6]);
    std::string outputFileName=argv[7];      float jtEtaCutLo= atof( argv[8] ) ; float jtEtaCutHi=atof( argv[9] ) ;
    rStatus = readForests_ppData_makeNTuple( inputFileList, startfile, endfile, 
				  jetRadius, jetType, debug,
					   outputFileName, jtEtaCutLo, jtEtaCutHi);
  }
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}
