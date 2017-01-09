// custom header
#include "readForests.h"

// ppData switches
const bool fillDataEvtQAHists=true, fillDataJetQAHists=true, fillDataJetIDHists=true; //most basic-level plots
const bool fillDataJetTrigQAHists=true; //data-specific
const bool fillDataJetSpectraRapHists=true; //other
const bool tightJetID=true;

//// readForests_ppData_jetPlots
// ---------------------------------------------------------------------------------------------------------------
int readForests_ppData_jetPlots( std::string inFilelist , int startfile , int endfile , 
			int radius , std::string jetType , bool debugMode ,
		        std::string outfile ){

  // for monitoring performance + debugging
  TStopwatch timer;  timer.Start();
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl;

  // basic info the screen
  std::cout<<std::endl<<"///////////////////"<<std::endl;
  std::cout<<"reading filelist "<< inFilelist<<std::endl;
  std::cout<<"reading files #'s "<< startfile << " to " << endfile<<std::endl;
  std::cout<<"radius = " << radius;
  std::cout<<", jetType = " << jetType;
  std::cout<<", debugMode = "<<debugMode<<std::endl;
  if(fastDebugMode&&debugMode)std::cout<<"fastDebugMode is ON"<<std::endl;
  std::cout<<"///////////////////"<<std::endl<<std::endl;

  // plot settings for later
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  

  //check out what kind of filelist we're using
  bool filelistIsLowerJets=false; //to know if i should apply the duplicate check or not
  if( inFilelist.find("HighPtLowerJets")!=std::string::npos ) {
    filelistIsLowerJets=true; 
    std::cout<<"running LowerJets filelist!"<<std::endl; }
  bool filelistIsJet80=false;
  if( inFilelist.find("HighPtJet80")    !=std::string::npos ) {
    filelistIsJet80=true;
    std::cout<<"running Jet80 filelist!"<<std::endl; }
  assert( !(filelistIsJet80&&filelistIsLowerJets) );
  

  // jet tree name+directory formation
  std::string jetTreeName="ak"+std::to_string(radius)+jetType+"JetAnalyzer/t";
  if(debugMode)std::cout<<"looking at jetTree "<<jetTreeName<<std::endl;

  // initialize tree name array
  std::string trees[N_dataTrees];
  trees[0]=jetTreeName;
  for(int i=1;i<N_dataTrees;++i)trees[i]=dataTreeNames[i];
  if(debugMode)for(int i=0;i<N_dataTrees;++i)std::cout<<"trees[i="<<i<<"]="<<trees[i]<<std::endl;  

  // declare TChains for each tree + friend them
  TChain* jetpp[N_dataTrees];
  for(int t = 0;t<N_dataTrees;++t)  { 
    jetpp[t] = new TChain( trees[t].data() );
    if(t>0)jetpp[0]->AddFriend( jetpp[t] );      }

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
    for(int t = 0;t<N_dataTrees;++t) filesAdded=jetpp[t]->AddFile(filename_Forest.c_str()); 
    lastFileAdded=filename_Forest;  }//end file loop
  assert(filesAdded);//avoid segfault later


  // Declare the output file, declare hists after for easy write
  std::cout<<"opening output file "<<outfile<<std::endl;
  TFile *fout = new TFile(outfile.c_str(),"RECREATE");
  fout->cd();

  // declare hists
  TH1F *hJetPtCut        =new TH1F("hJetPtCut"      ,(std::to_string(jtPtCut)).c_str()   ,    100, 0,100); hJetPtCut->Fill(jtPtCut);           
  TH1F *hJetEtaCut       =new TH1F("hJetEtaCut"     ,(std::to_string(jtEtaCut)).c_str()  ,    5,   0,5  ); hJetEtaCut->Fill(jtEtaCut);	        
  TH1F *hJetQAPtCut      =new TH1F("hJetQAPtCut"    ,(std::to_string(jetQAPtCut)).c_str(),    100, 0,100); hJetQAPtCut->Fill(jetQAPtCut);     
  TH1F *hLeadJetPtCut    =new TH1F("hLdJetPtCut"    ,(std::to_string(ldJetPtCut)).c_str(),    100, 0,100); hLeadJetPtCut->Fill(ldJetPtCut);     
  TH1F *hSubLeadJetPtCut =new TH1F("hSubldJetPtCut" ,(std::to_string(subldJetPtCut)).c_str(), 100, 0,100); hSubLeadJetPtCut->Fill(subldJetPtCut);
  TH1F *hPtAveCut        =new TH1F("hPtAveCut"      ,(std::to_string(ptAveCut)).c_str(),      100, 0,100); hPtAveCut->Fill(ptAveCut);  
 

  //evt count
  TH1F *h_NEvents         = new TH1F("NEvents"         , "NEvents", 1,0,2);
  TH1F *h_NEvents_read    = new TH1F("NEvents_read"    , "NEvents read", 1,0,2);
  TH1F *h_NEvents_skimCut = new TH1F("NEvents_skimCut" , "NEvents read post skimCut", 1,0,2);
  TH1F *h_NEvents_vzCut   = new TH1F("NEvents_vzCut"   , "NEvents read post vzCut AND skimCut", 1,0,2);

  //evt counts post all evt cuts. thrown out because no jets left in the event pre or post cuts
  TH1F *h_NEvents_withJets           = new TH1F("NEvents_withJets" , 
		  				"NEvents read post evt cuts, w/ jets", 1,0,2);
  TH1F *h_NEvents_withJets_kmatCut   = new TH1F("NEvents_withJets_kmatCut" , 
						"NEvents read post evt cuts, w/ jets post kmatCut", 1,0,2);
  TH1F *h_NEvents_withJets_JetIDCut = NULL;
  if(fillDataJetIDHists)
    h_NEvents_withJets_JetIDCut= new TH1F("NEvents_withJets_JetIDCut" , 
					  "NEvents read post evt cuts, w/ jets post kmatCut AND JetID Cut", 1,0,2);
  
  
  //jet counts
  TH1F *h_NJets          = new TH1F("NJets","NJets read", 1,0,2);
  TH1F *h_NJets_kmatCut  = new TH1F("NJets_kmatCut ","NJets read post kmatCut ", 1,0,2);
  TH1F *h_NJets_JetIDCut=NULL;
  if(fillDataJetIDHists)
    h_NJets_JetIDCut = new TH1F("NJets_JetIDCut","NJets read post JetIDCut AND kmatCut", 1,0,2);
  
  //vz
  TH1F *hVz=NULL, *hWVz=NULL;
  if(fillDataEvtQAHists){
    hVz = new TH1F("hVz","",  60,-15.,15.); 
    hWVz = new TH1F("hWeightedVz","",  60,-15.,15.);    }
  
  
  //jet QA
  TH1F *hJetQA[2][N_vars]={};
  if(fillDataJetQAHists)
    for(int k = 0; k<2; ++k){
      if(!fillDataJetIDHists && k==1)continue;	
      for(int j = 0; j<N_vars; ++j){	
	//binnings+hist declarations, etc.
	if(var[j]=="jtpt"||var[j]=="rawpt")
	  hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;", var[j].c_str()) , 500,0,500);       
	else if(var[j]=="jteta"||var[j]=="jty")  //rapidity + pseudorapidity
	  hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;",var[j].c_str()) , 100,-5,+5);
	else if(var[j]=="jtphi") 
	  hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;",var[j].c_str()) , 100,-4,+4);
	else if(var[j]=="dphi") 
	  hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;",var[j].c_str()) , 50,0,+4);
	else if(var[j]=="leadJetPt"||var[j]=="subleadJetPt") 
	  hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;", var[j].c_str()), 500,0,500);
	else if(var[j]=="neN"||var[j]=="sumN"||var[j]=="chN")
	  hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 30,0,30);         
	else //consituent + xj and Aj binnings
	  hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 200,0,2);         
     } }
  
  TH1F *hJetSpectraRap[2][nbins_rap]={};
  if(fillDataJetSpectraRapHists)
    for(int k = 0; k<2; ++k) {
      if(!fillDataJetIDHists && k==1)continue;	
      for(int j = 0; j<nbins_rap; ++j){
	std::string h_Title="hJetSpectraRap_";
	if(k==1)h_Title+="wJetID_";
	h_Title+="bin"+std::to_string(j);
	std::string h_Desc;
	std::stringstream stream1, stream2;	
	stream1.precision(1); stream1 << std::fixed << rapbins[j];
	stream2.precision(1); stream2 << std::fixed << rapbins[j+1];
	h_Desc="JetPt Spectra for "+stream1.str()+"<abs(y)<"+ stream2.str();	
	hJetSpectraRap[k][j]=new TH1F(h_Title.c_str(),h_Desc.c_str(), 2000,0,2000);  
      } }

  // trigger (also jet) level
  TH1F  *hpp_TrgObj40[2]={}, *hpp_TrgObj60[2]={}, *hpp_TrgObj80[2]={}, *hpp_TrgObj100[2]={};// separate trigs
  TH1F *hpp_TrgObjComb[2]={}, *hpp_CombJetpT[2]={};// combos of trigs
  //TH2F *hTrgQA_recpt_rawpt=NULL, *hTrgQA_recpt_trgpt=NULL, *hTrgQA_trgpt_rawpt=NULL;

  if(fillDataJetTrigQAHists){
    std::string hJetTrigQATitleArray[]={ "HLT40","HLT60",  "HLT80", "HLT100", "HLTComb","TrgCombTest" };
    const int nJetTrigQATitles=sizeof(hJetTrigQATitleArray)/sizeof(std::string);
    for(int j=0;j<2;j++){      
      if(!fillDataJetIDHists && j==1) continue;    
      std::cout<<"j=="<<j<<std::endl;
      for(int k=0; k<nJetTrigQATitles; k++){
	std::cout<<"k=="<<k<<std::endl;
	std::string hTitle="hpp_"+hJetTrigQATitleArray[k];
	if(j==1)hTitle+="_wJetID";
	hTitle+="_R"+std::to_string(radius)+"_"+etaWidth;
	std::cout<<"hTitle="<<hTitle<<std::endl;
	if(hJetTrigQATitleArray[k]=="HLT40")
	  hpp_TrgObj40[j]   = new TH1F(hTitle.c_str(), "HLT40Jet", 2000,0,2000);
	else if(hJetTrigQATitleArray[k]=="HLT60")
	  hpp_TrgObj60[j]   = new TH1F(hTitle.c_str(), "HLT60Jet", 2000,0,2000);
	else if(hJetTrigQATitleArray[k]=="HLT80")
	  hpp_TrgObj80[j]   = new TH1F(hTitle.c_str(), "HLT80Jet", 2000,0,2000);
	else if(hJetTrigQATitleArray[k]=="HLT100")
	  hpp_TrgObj100[j]  = new TH1F(hTitle.c_str(), "HLT100Jet", 2000,0,2000);	
	else if(hJetTrigQATitleArray[k]=="HLTComb")
	  hpp_TrgObjComb[j] = new TH1F(hTitle.c_str(), "HLTCombination", 2000,0,2000);
	else if(hJetTrigQATitleArray[k]=="TrgCombTest")
	  hpp_CombJetpT[j]  = new TH1F(hTitle.c_str(), "HLTAltCombo", 2000,0,2000);
      }}
  }

  // EVENT LOOP + PREP
  // declare variables/arrays + set branch address for each input tree
  //JetAnalyzer, jets
  int nref_I;  
  float pt_F[1000];  float rawpt_F[1000];
  float eta_F[1000];   float phi_F[1000];  
  float y_F[1000]; //float jtpu_F[1000]; 
  
  //tracks 
  float trkSum_F[1000]; float trkHardSum_F[1000];  
  float trkMax_F[1000]; 
  //charged particles          
  int chN_I[1000];  
  float chSum_F[1000]; float chHardSum_F[1000];  
  float chMax_F[1000]; 
  //photons  
  float phSum_F[1000]; float phHardSum_F[1000];   
  float phMax_F[1000]; 
  //leptons
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
  jetpp[0]->SetBranchAddress("jty",&y_F);
  //jetpp[0]->SetBranchAddress("jtpu",&jtpu_F);
  //tracks
  jetpp[0]->SetBranchAddress("trackSum",&trkSum_F);  
  jetpp[0]->SetBranchAddress("trackHardSum",&trkHardSum_F);
  jetpp[0]->SetBranchAddress("trackMax",&trkMax_F);  
  //charged hadrons
  jetpp[0]->SetBranchAddress("chargedN",&chN_I);  
  jetpp[0]->SetBranchAddress("chargedSum",&chSum_F);
  jetpp[0]->SetBranchAddress("chargedHardSum",&chHardSum_F);
  jetpp[0]->SetBranchAddress("chargedMax",&chMax_F);  
  //neutral hadrons
  jetpp[0]->SetBranchAddress("neutralN",&neN_I);  
  jetpp[0]->SetBranchAddress("neutralSum",&neSum_F); 
  jetpp[0]->SetBranchAddress("neutralMax",&neMax_F);
  //photons
  jetpp[0]->SetBranchAddress("photonSum",&phSum_F);  
  jetpp[0]->SetBranchAddress("photonHardSum",&phHardSum_F);
  jetpp[0]->SetBranchAddress("photonMax",&phMax_F);  
  //leptons
  jetpp[0]->SetBranchAddress("eSum",&eSum_F);  
  jetpp[0]->SetBranchAddress("eMax",&eMax_F);
  jetpp[0]->SetBranchAddress("muSum",&muSum_F);
  jetpp[0]->SetBranchAddress("muMax",&muMax_F);


  // HiEvtAnalyzer
  float vz_F;
  ULong64_t evt_I;   UInt_t run_I;   UInt_t lumi_I;   
  jetpp[1]->SetBranchAddress("vz",&vz_F);
  jetpp[1]->SetBranchAddress("evt",&evt_I);
  jetpp[1]->SetBranchAddress("run",&run_I);
  jetpp[1]->SetBranchAddress("lumi",&lumi_I);


  // skimanalysis
  int pBeamScrapingFilter_I, pHBHENoiseFilter_I, pprimaryvertexFilter_I;
  jetpp[2]->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_I);
  jetpp[2]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilter_I);
  jetpp[2]->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryvertexFilter_I);


  // hltanalysis
  assert(N_HLTBits==N_L1Bits);

  //string array formation for HLT branches
  std::string PF_HLTBranches[N_HLTBits], Calo_HLTBranches[N_HLTBits];
  std::string PF_HLTPresclBranches[N_HLTBits], Calo_HLTPresclBranches[N_HLTBits];
  std::string L1PresclBranches[N_L1Bits];

  for(int i=0; i<N_HLTBits; ++i) {
    L1PresclBranches[i] = L1BitStrings[i] + "_Prescl";  
    PF_HLTBranches[i] = PF_HLTBitStrings[i] + "_v1"; 
    PF_HLTPresclBranches[i] = PF_HLTBitStrings[i] + "_v1_Prescl";   
    Calo_HLTBranches[i] = Calo_HLTBitStrings[i] + "_v1"; 
    Calo_HLTPresclBranches[i] = Calo_HLTBitStrings[i] + "_v1_Prescl"; 
    
    if(debugMode){std::cout<<std::endl;
      std::cout<<"L1PresclBranches      [i="<<i<<"]= "<<L1PresclBranches[i]<<std::endl;
      std::cout<<"PF_HLTBranches        [i="<<i<<"]= "<<PF_HLTBranches[i]<<std::endl;
      std::cout<<"PF_HLTPresclBranches  [i="<<i<<"]= "<<PF_HLTPresclBranches[i]<<std::endl;
      std::cout<<"Calo_HLTBranches      [i="<<i<<"]= "<<Calo_HLTBranches[i]<<std::endl;
      std::cout<<"Calo_HLTPresclBranches[i="<<i<<"]= "<<Calo_HLTPresclBranches[i]<<std::endl<<std::endl; }  }
  
  int PFJet40_I, PFJet60_I, PFJet80_I, PFJet100_I;   //PF HLT
  int PFJet40_p_I, PFJet60_p_I, PFJet80_p_I, PFJet100_p_I;
  int CaloJet40_I, CaloJet60_I, CaloJet80_I, CaloJet100_I;   //Calo HLT
  int CaloJet40_p_I, CaloJet60_p_I, CaloJet80_p_I, CaloJet100_p_I;
  int jet40_l1s_ps_I, jet60_l1s_ps_I,  jet80_l1s_ps_I,  jet100_l1s_ps_I;   //L1 ps

  jetpp[3]->SetBranchAddress( PF_HLTBranches[0].c_str() , &PFJet40_I);
  jetpp[3]->SetBranchAddress( PF_HLTBranches[1].c_str() , &PFJet60_I);
  jetpp[3]->SetBranchAddress( PF_HLTBranches[2].c_str() , &PFJet80_I);
  jetpp[3]->SetBranchAddress( PF_HLTBranches[3].c_str() , &PFJet100_I);

  jetpp[3]->SetBranchAddress( PF_HLTPresclBranches[0].c_str() , &PFJet40_p_I);
  jetpp[3]->SetBranchAddress( PF_HLTPresclBranches[1].c_str() , &PFJet60_p_I);
  jetpp[3]->SetBranchAddress( PF_HLTPresclBranches[2].c_str() , &PFJet80_p_I);
  jetpp[3]->SetBranchAddress( PF_HLTPresclBranches[3].c_str() , &PFJet100_p_I);

  jetpp[3]->SetBranchAddress( Calo_HLTBranches[0].c_str() , &CaloJet40_I);
  jetpp[3]->SetBranchAddress( Calo_HLTBranches[1].c_str() , &CaloJet60_I);
  jetpp[3]->SetBranchAddress( Calo_HLTBranches[2].c_str() , &CaloJet80_I);
  jetpp[3]->SetBranchAddress( Calo_HLTBranches[3].c_str() , &CaloJet100_I);

  jetpp[3]->SetBranchAddress( Calo_HLTPresclBranches[0].c_str() , &CaloJet40_p_I);
  jetpp[3]->SetBranchAddress( Calo_HLTPresclBranches[1].c_str() , &CaloJet60_p_I);
  jetpp[3]->SetBranchAddress( Calo_HLTPresclBranches[2].c_str() , &CaloJet80_p_I);
  jetpp[3]->SetBranchAddress( Calo_HLTPresclBranches[3].c_str() , &CaloJet100_p_I);
  
  jetpp[3]->SetBranchAddress( L1PresclBranches[0].c_str()  , &jet40_l1s_ps_I);
  jetpp[3]->SetBranchAddress( L1PresclBranches[1].c_str()  , &jet60_l1s_ps_I);
  jetpp[3]->SetBranchAddress( L1PresclBranches[2].c_str()  , &jet80_l1s_ps_I);
  jetpp[3]->SetBranchAddress( L1PresclBranches[3].c_str()  , &jet100_l1s_ps_I);

  //ONE HLT path ONE tree ONE trig obj pt branch
  //e.g. trgObjpt_40 is filled with jet pt from the specific jet40 HLT tree/branch 
  std::vector<double>  *trgObjpt_40, *trgObjpt_60, *trgObjpt_80, *trgObjpt_100;
  jetpp[4]->SetBranchAddress("pt",&trgObjpt_40);
  jetpp[5]->SetBranchAddress("pt",&trgObjpt_60);  
  jetpp[6]->SetBranchAddress("pt",&trgObjpt_80);  
  jetpp[7]->SetBranchAddress("pt",&trgObjpt_100);

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
  for(UInt_t i=0;i < NEvents_allFiles; ++i) h_NEvents->Fill(1);
  
  UInt_t NEvents_read=0;
  if(fastDebugMode&&debugMode) NEvents_read = 10000*(endfile-startfile+1); 
  else NEvents_read = NEvents_allFiles;
  
  UInt_t NEvents_40=0, NEvents_60=0, NEvents_80=0, NEvents_100=0;  //trigger event counts
  std::cout<<"reading "<<NEvents_read<<" events"<<std::endl;  
  
  L2ResidualJES* L2JES = new L2ResidualJES(radius, 3, "pp5");
  L3ResidualJES* L3JES = new L3ResidualJES("pp5");

  float jetIDCut_neSum, jetIDCut_phSum;
  if(tightJetID){     jetIDCut_neSum=0.90;  jetIDCut_phSum=0.90;}
  else{     jetIDCut_neSum=0.99;  jetIDCut_phSum=0.99;}
  
  //begin event loop
  for(UInt_t nEvt = 0; nEvt < NEvents_read; ++nEvt) {
    
    //grab an entry
    if( (fastDebugMode || debugMode) && nEvt%100==0) std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;
    else if (nEvt%10000==0) std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;
    jetpp[0]->GetEntry(nEvt);
    h_NEvents_read->Fill(1);

    //duplicate skipping between LowerJets and Jet80
    if(filelistIsJet80)
      if( (bool)CaloJet40_I || (bool)CaloJet60_I || (bool)PFJet40_I || (bool)PFJet60_I ) {
	if(debugMode)std::cout<<"this event is in Jet80 AND LowerJets dataset.!"<<std::endl;
	if(debugMode)std::cout<<"Skipping event, will read it in LowerJets instead"<<std::endl;	continue; }
    
    
    // skim/HiEvtAnalysis criteria
    if( pHBHENoiseFilter_I==0     || 
        pBeamScrapingFilter_I==0  || 
        pprimaryvertexFilter_I==0  ) continue;    
    h_NEvents_skimCut->Fill(1);

    if( fabs(vz_F)>15. )     continue;
    h_NEvents_vzCut->Fill(1);
    
    //prescale/decision arrays, total prescale=L1_ps*HLT_ps
    bool PFtrgDec[N_HLTBits]   ={ (bool)PFJet40_I, (bool)PFJet60_I, (bool)PFJet80_I, (bool)PFJet100_I };
    int PFtrgPrescl[N_HLTBits] ={ PFJet40_p_I*jet40_l1s_ps_I , PFJet60_p_I*jet60_l1s_ps_I , 
				  PFJet80_p_I*jet80_l1s_ps_I , PFJet100_p_I*jet100_l1s_ps_I };
    
    bool CalotrgDec[N_HLTBits]   ={ (bool)CaloJet40_I, (bool)CaloJet60_I, (bool)CaloJet80_I, (bool)CaloJet100_I};
    int CalotrgPrescl[N_HLTBits] ={ CaloJet40_p_I*jet40_l1s_ps_I , CaloJet60_p_I*jet60_l1s_ps_I,
				    CaloJet80_p_I*jet80_l1s_ps_I , CaloJet100_p_I*jet100_l1s_ps_I }; 
    
    bool *trgDec=NULL; int *trgPscl=NULL;
    if(trgCombType=="Calo")   { trgDec=CalotrgDec ;   trgPscl=CalotrgPrescl ; }
    else if(trgCombType=="PF"){ trgDec=PFtrgDec   ;   trgPscl=PFtrgPrescl   ; }
    else assert(false);//should never fire if this is working right
    
    
    // grab largest triggerPt among HLT paths' trigger object collectiosn
    double maxTrgPt=0.; //, unsigned int maxTrgPtIndex=-1;
    unsigned int trgObj40_size=trgObjpt_40->size(), trgObj60_size=trgObjpt_60->size();
    unsigned int trgObj80_size=trgObjpt_80->size(), trgObj100_size=trgObjpt_100->size();
    
    if(trgDec[3])
      for(unsigned int itt=0; itt<trgObj100_size; ++itt)
    	if(trgObjpt_100->at(itt) > maxTrgPt) { maxTrgPt = trgObjpt_100->at(itt); } //maxTrgPtIndex=itt; }

    if(trgDec[2])
      for(unsigned int itt=0; itt<trgObj80_size; ++itt)
    	if(trgObjpt_80->at(itt) > maxTrgPt)  { maxTrgPt = trgObjpt_80->at(itt); } //maxTrgPtIndex=itt; }

    if(trgDec[1]) 
      for(unsigned int itt=0; itt<trgObj60_size; ++itt)
	if(trgObjpt_60->at(itt) > maxTrgPt)  { maxTrgPt = trgObjpt_60->at(itt); } //maxTrgPtIndex=itt; }

    if(trgDec[0])  
      for(unsigned int itt=0; itt<trgObj40_size; ++itt)	
	if(trgObjpt_40->at(itt) > maxTrgPt)  { maxTrgPt = trgObjpt_40->at(itt); } //maxTrgPtIndex=itt; }
    
    double trgPt=maxTrgPt; //unsigned int trgPtIndex=maxTrgPtIndex;

    // check trigger decisions for events + exclusivity between them, count events, assign prescale weight
    float weight_eS=0.;
    //float weight_eS = trigComb(trgDec, trgPscl, trgPt); // trig comb function replicates the procedure below

    bool is40  = false, is60  = false, is80  = false, is100 = false;
    if(      trgDec[3] && trgPt>=100.               ) 
      { is100 = true;  weight_eS=trgPscl[3]; }
    else if( trgDec[2] && trgPt>=80.  && trgPt<100. ) 
      { is80  = true;  weight_eS=trgPscl[2]; }
    else if( trgDec[1] && trgPt>=60.  && trgPt<80.  ) 
      { is60  = true;  weight_eS=trgPscl[1]; }
    else if( trgDec[0] && trgPt>=40.  && trgPt<60.  ) 
      { is40  = true;  weight_eS=trgPscl[0]; }
    
    if(is100) NEvents_100++;
    else if(is80)  NEvents_80++; 
    else if(is60)  NEvents_60++; 
    else if(is40)  NEvents_40++; 

    // fill evt vz histo
    if(fillDataEvtQAHists){
      hVz->Fill(vz_F, 1.0);    
      hWVz->Fill(vz_F, weight_eS);    }
    
    // JET LOOP 
    // for Aj/xj computation
    bool firstGoodJetFound=false, secondGoodJetFound=false; 
    float firstGoodJetPt=-1., secondGoodJetPt=-1.;
    float firstGoodJetPhi=-1., secondGoodJetPhi=-1.;
    
    bool firstGoodJetFound_wJetID=false, secondGoodJetFound_wJetID=false; 
    float firstGoodJetPt_wJetID=-1., secondGoodJetPt_wJetID=-1.;
    float firstGoodJetPhi_wJetID=-1., secondGoodJetPhi_wJetID=-1.;
    
    // for event counting + avoiding duplicate fills in dijet hists
    bool dijetHistsFilled=false, dijetHistsFilled_wJetID=false;
    bool hNEvts_withJets_Filled=false, hNEvts_withJets_kmatCut_Filled=false, hNEvts_withJets_JetIDCut_Filled=false; 

    // for trg/reco matching
    for(int jet = 0; jet<nref_I; ++jet){

      // event+jet counting
      h_NJets->Fill(1);
      if(!hNEvts_withJets_Filled){
	h_NEvents_withJets->Fill(1);
	hNEvts_withJets_Filled=true;      }      

      float recpt  = pt_F[jet];
      float receta = eta_F[jet];

      //if(debugMode)std::cout<<"recpt="<<recpt<<std::endl;
      recpt=L2JES->getCorrectedPt(recpt,receta);
      //if(debugMode)std::cout<<"after L2Residual, recpt="<<recpt<<std::endl;
      recpt=L3JES->getCorrectedPt(recpt);
      //if(debugMode)std::cout<<"after L3Residual, recpt="<<recpt<<std::endl;

      //// TEMP 10.12.16////
      // kmatCuts      
      if( recpt < jtPtCut ) continue;     

      float rawpt  = rawpt_F[jet];
      float recy   = y_F[jet];
      float recphi = phi_F[jet];
      
      // 13 TeV JetID criterion, loose or tight
      bool passesJetID=false;
      if(fillDataJetIDHists) {	
	if( fabs(recy)>3.0 ) 	  {
	  if( phSum_F[jet]/rawpt < 0.90 &&
	      neSum_F[jet]/rawpt<0.99 )//neN_I[jet] > 10) { 
	    passesJetID=true; 
	}
	else if ( 2.7<fabs(recy) && fabs(recy)<=3.0  )  {
	  if( phSum_F[jet]/rawpt < 0.90 && //neutral em
	      neN_I[jet] > 2  )
	    passesJetID=true;	
	} 	
	//if( fabs(recy)>2.7) {//jetIDv2
	//  if( phSum_F[jet]/rawpt<0.90 &&
	//      neSum_F[jet]/rawpt<0.99 &&
	//      ((phSum_F[jet]/rawpt>0.) || (neSum_F[jet]/rawpt>0.)) )
	//    passesJetID=true;
	//}
	else if ( 2.4<fabs(recy) && fabs(recy)<=2.7 ) {
	  if( neSum_F[jet]/rawpt    < jetIDCut_neSum &&  //neutral had 
	      phSum_F[jet]/rawpt    < jetIDCut_phSum && 
	      chN_I[jet]+neN_I[jet] > 1 ) 
	    passesJetID=true; 	  
	}
	else { //if(fabs(recy)<=2.4) //in the barrel, strictest
	  if( neSum_F[jet]/rawpt    < jetIDCut_neSum &&
	      phSum_F[jet]/rawpt    < jetIDCut_phSum &&
	      chN_I[jet]+neN_I[jet] > 1    &&
	      chSum_F[jet]/rawpt    > 0.   && //charged had 
	      eSum_F[jet]/rawpt     < 0.99 && 
	      chN_I[jet] > 0 )      	  
	    passesJetID=true;  
	} }
      
      //fill jetspectraRapHists w/ passing jetID criterion
      if( fillDataJetSpectraRapHists ) { 
	int theRapBin=-1;
	for(int rapbin=0;rapbin<nbins_rap;++rapbin)
	  if( rapbins[rapbin]<=fabs(recy)  && 		
	      fabs(recy)<rapbins[rapbin+1]    ) {
	    theRapBin=rapbin;
	    hJetSpectraRap[0][theRapBin]->Fill(recpt,weight_eS);  
	    if( passesJetID ) 
	      hJetSpectraRap[1][theRapBin]->Fill(recpt,weight_eS);    
	    break;
	  } }
      
      //second half of kmat cut      
      if( fabs(receta) > jtEtaCut ) continue;


      // jet/event counts
      h_NJets_kmatCut->Fill(1);
      if(!hNEvts_withJets_kmatCut_Filled){
	h_NEvents_withJets_kmatCut->Fill(1);
	hNEvts_withJets_kmatCut_Filled=true;  }      

      // trig plots
      if(fillDataJetTrigQAHists){
	if(true)  hpp_CombJetpT[0]->Fill(recpt, weight_eS); //kurts method
	if(trgDec[0])/*if(is40 )*/  hpp_TrgObj40[0]->Fill(recpt, trgPscl[0]);
	if(trgDec[1])/*if(is60 )*/  hpp_TrgObj60[0]->Fill(recpt, trgPscl[1]);
	if(trgDec[2])/*if(is80 )*/  hpp_TrgObj80[0]->Fill(recpt, trgPscl[2]);
	if(trgDec[3])/*if(is100)*/  hpp_TrgObj100[0]->Fill(recpt, trgPscl[3]); }	

      // jetQA noJetID
      if(fillDataJetQAHists){

	hJetQA[0][0]->Fill(recpt, weight_eS);
	hJetQA[0][1]->Fill(rawpt_F[jet], weight_eS);
	if(recpt>jetQAPtCut){
	  hJetQA[0][2]->Fill(eta_F[jet], weight_eS);
	  hJetQA[0][3]->Fill(phi_F[jet], weight_eS);
	  hJetQA[0][4]->Fill(y_F[jet], weight_eS);	  	  
	  hJetQA[0][5]->Fill(trkMax_F[jet]/rawpt, weight_eS);
	  hJetQA[0][6]->Fill(trkSum_F[jet]/rawpt, weight_eS);
	  hJetQA[0][7]->Fill(trkHardSum_F[jet]/rawpt, weight_eS);
	  hJetQA[0][8]->Fill(chMax_F[jet]/rawpt, weight_eS);
	  hJetQA[0][9]->Fill(chSum_F[jet]/rawpt, weight_eS);
	  hJetQA[0][10]->Fill(chHardSum_F[jet]/rawpt, weight_eS);
	  hJetQA[0][11]->Fill(phMax_F[jet]/rawpt, weight_eS);
	  hJetQA[0][12]->Fill(phSum_F[jet]/rawpt, weight_eS);
	  hJetQA[0][13]->Fill(phHardSum_F[jet]/rawpt, weight_eS);
	  hJetQA[0][14]->Fill(neMax_F[jet]/rawpt, weight_eS);
	  hJetQA[0][15]->Fill(neSum_F[jet]/rawpt, weight_eS);
	  hJetQA[0][16]->Fill(eMax_F[jet]/rawpt, weight_eS);
	  hJetQA[0][17]->Fill(eSum_F[jet]/rawpt, weight_eS);
	  hJetQA[0][18]->Fill(muMax_F[jet]/rawpt, weight_eS);
	  hJetQA[0][19]->Fill(muSum_F[jet]/rawpt, weight_eS);   
	  hJetQA[0][20]->Fill(neN_I[jet], weight_eS);   
	  hJetQA[0][21]->Fill(chN_I[jet], weight_eS);   
	  hJetQA[0][22]->Fill(chN_I[jet]+neN_I[jet], weight_eS);   	}
	
	  //looking for the first two good jets that meet the criteria specified
	  if ( !firstGoodJetFound ){
	    if(recpt>ldJetPtCut) { firstGoodJetFound=true;
	      firstGoodJetPt =recpt; 
	     firstGoodJetPhi=recphi; 
	    } 	  }
	  
	  else if ( !secondGoodJetFound && 
		    firstGoodJetFound      ) {
	    if (recpt>subldJetPtCut) { 
	      float checkdPhi=deltaphi(firstGoodJetPhi,recphi);
	      float ptAve=(firstGoodJetPt+recpt)/2.;
	      if( checkdPhi>dPhiCut && ptAve>ptAveCut) {
		if(recpt>firstGoodJetPt){secondGoodJetFound=true;
		  std::cout<<std::endl<<
		    "WARNING: picked wrong jet for lead jet! Swapping..."<<std::endl<<std::endl;
		  secondGoodJetPt  = firstGoodJetPt;          
		  secondGoodJetPhi = firstGoodJetPhi;
		  firstGoodJetPt   = recpt;          
		  firstGoodJetPhi  = recphi; }
		else { secondGoodJetFound=true;
		  secondGoodJetPt  = recpt;
		  secondGoodJetPhi = recphi; }
	      }	    }	  }
	  
	  if (!dijetHistsFilled &&
	      firstGoodJetFound && 
	      secondGoodJetFound ) { dijetHistsFilled=true;
	    float A_j=(firstGoodJetPt-secondGoodJetPt)
	      /(firstGoodJetPt+secondGoodJetPt);
	    float x_j=(secondGoodJetPt/firstGoodJetPt); 	
	    float dphi=deltaphi(firstGoodJetPhi,secondGoodJetPhi);
	    hJetQA[0][23]->Fill( A_j , weight_eS ); 
	    hJetQA[0][24]->Fill( x_j , weight_eS );       
	    hJetQA[0][25]->Fill( dphi , weight_eS );       
	    hJetQA[0][26]->Fill( firstGoodJetPt , weight_eS );       
	    hJetQA[0][27]->Fill( secondGoodJetPt , weight_eS );       } 
      }
	
            
      // apply JetID
      if(fillDataJetIDHists){
	if(!passesJetID) continue;

	// trig plots
	if(fillDataJetTrigQAHists){
	  if(true)  hpp_CombJetpT[1]->Fill(recpt, weight_eS); //kurts method
	  if(trgDec[0])/*(is40 )*/  hpp_TrgObj40[1]->Fill(recpt, trgPscl[0]);
	  if(trgDec[1])/*(is60 )*/  hpp_TrgObj60[1]->Fill(recpt, trgPscl[1]);
	  if(trgDec[2])/*(is80 )*/  hpp_TrgObj80[1]->Fill(recpt, trgPscl[2]);
	  if(trgDec[3])/*(is100)*/ hpp_TrgObj100[1]->Fill(recpt, trgPscl[3]); 	}

	// jetQA wJetID
	if(fillDataJetQAHists){
	  
	  hJetQA[1][0]->Fill(recpt, weight_eS);
	  hJetQA[1][1]->Fill(rawpt_F[jet], weight_eS);
	  if(recpt>jetQAPtCut){     
	    hJetQA[1][2]->Fill(eta_F[jet], weight_eS);
	    hJetQA[1][3]->Fill(phi_F[jet], weight_eS);
	    hJetQA[1][4]->Fill(y_F[jet], weight_eS);
	    hJetQA[1][5]->Fill(trkMax_F[jet]/rawpt, weight_eS);
	    hJetQA[1][6]->Fill(trkSum_F[jet]/rawpt, weight_eS);
	    hJetQA[1][7]->Fill(trkHardSum_F[jet]/rawpt, weight_eS);
	    hJetQA[1][8]->Fill(chMax_F[jet]/rawpt, weight_eS);
	    hJetQA[1][9]->Fill(chSum_F[jet]/rawpt, weight_eS);
	    hJetQA[1][10]->Fill(chHardSum_F[jet]/rawpt, weight_eS);
	    hJetQA[1][11]->Fill(phMax_F[jet]/rawpt, weight_eS);
	    hJetQA[1][12]->Fill(phSum_F[jet]/rawpt, weight_eS);
	    hJetQA[1][13]->Fill(phHardSum_F[jet]/rawpt, weight_eS);
	    hJetQA[1][14]->Fill(neMax_F[jet]/rawpt, weight_eS);
	    hJetQA[1][15]->Fill(neSum_F[jet]/rawpt, weight_eS);
	    hJetQA[1][16]->Fill(eMax_F[jet]/rawpt, weight_eS);
	    hJetQA[1][17]->Fill(eSum_F[jet]/rawpt, weight_eS);
	    hJetQA[1][18]->Fill(muMax_F[jet]/rawpt, weight_eS);
	    hJetQA[1][19]->Fill(muSum_F[jet]/rawpt, weight_eS);      	  
	    hJetQA[1][20]->Fill(neN_I[jet], weight_eS);   
	    hJetQA[1][21]->Fill(chN_I[jet], weight_eS);   
	    hJetQA[1][22]->Fill(chN_I[jet]+neN_I[jet], weight_eS);   	}

	  
	  
	  //looking for the first two good jets that meet the criteria specified
	  if ( !firstGoodJetFound_wJetID ){
	    if(recpt>ldJetPtCut) { firstGoodJetFound_wJetID=true;
	      firstGoodJetPt_wJetID =recpt; 
	      firstGoodJetPhi_wJetID=recphi; } }
	  
	  else if ( !secondGoodJetFound_wJetID && 
		    firstGoodJetFound_wJetID      ) {
	    if (recpt>subldJetPtCut) { 
	      float checkdPhi=deltaphi(firstGoodJetPhi_wJetID,recphi);
	      float ptAve=(firstGoodJetPt_wJetID+recpt)/2.;
	      if( checkdPhi>dPhiCut && ptAve>ptAveCut) {
		if(recpt>firstGoodJetPt_wJetID){secondGoodJetFound_wJetID=true;
		  std::cout<<std::endl<<
		    "WARNING: picked wrong jet for lead jet! Swapping..."<<std::endl<<std::endl;
		  secondGoodJetPt_wJetID  = firstGoodJetPt_wJetID;          
		  secondGoodJetPhi_wJetID = firstGoodJetPhi_wJetID;
		  firstGoodJetPt_wJetID   = recpt;          
		  firstGoodJetPhi_wJetID  = recphi; }
		else { secondGoodJetFound_wJetID=true;
		  secondGoodJetPt_wJetID  = recpt;
		  secondGoodJetPhi_wJetID = recphi; }}} }
	  
	  if (!dijetHistsFilled_wJetID &&
	      firstGoodJetFound_wJetID && 
	      secondGoodJetFound_wJetID ) { dijetHistsFilled_wJetID=true;
	    float A_j=(firstGoodJetPt_wJetID-secondGoodJetPt_wJetID)
	      /(firstGoodJetPt_wJetID+secondGoodJetPt_wJetID);
	    float x_j=(secondGoodJetPt_wJetID/firstGoodJetPt_wJetID); 	
	    float dphi=deltaphi(firstGoodJetPhi_wJetID,secondGoodJetPhi_wJetID);
	    hJetQA[1][23]->Fill( A_j , weight_eS ); 
	    hJetQA[1][24]->Fill( x_j , weight_eS );       
	    hJetQA[1][25]->Fill( dphi , weight_eS );       
	    hJetQA[1][26]->Fill( firstGoodJetPt_wJetID , weight_eS );       
	    hJetQA[1][27]->Fill( secondGoodJetPt_wJetID , weight_eS );       } 
	}//end fillDataJetQAhists
	
	// jet/event counts
	h_NJets_JetIDCut->Fill(1);
	if(!hNEvts_withJets_JetIDCut_Filled){
	  h_NEvents_withJets_JetIDCut->Fill(1);
	  hNEvts_withJets_JetIDCut_Filled=true;   }
	
      }//end fillDataAJetIDhists

    }//end jet loop
  }//end event loop
  std::cout<<std::endl;

  if(fillDataJetTrigQAHists){

    if(fillDataJetIDHists){
      hpp_TrgObj40[1]->Print("base");
      hpp_TrgObj60[1]->Print("base");
      hpp_TrgObj80[1]->Print("base");
      hpp_TrgObj100[1]->Print("base");
      
      hpp_TrgObjComb[1]->Add(hpp_TrgObj100[1]);
      hpp_TrgObjComb[1]->Add(hpp_TrgObj80[1] );
      hpp_TrgObjComb[1]->Add(hpp_TrgObj60[1] );
      hpp_TrgObjComb[1]->Add(hpp_TrgObj40[1] );
      hpp_TrgObjComb[0]->Scale(1./4.);

      hpp_TrgObjComb[1]->Print("base");      
      hpp_CombJetpT[1]->Print("base"); }

    hpp_TrgObj40[0]->Print("base");
    hpp_TrgObj60[0]->Print("base");
    hpp_TrgObj80[0]->Print("base");
    hpp_TrgObj100[0]->Print("base");
    
    hpp_TrgObjComb[0]->Add(hpp_TrgObj100[0]);
    hpp_TrgObjComb[0]->Add(hpp_TrgObj80[0] );
    hpp_TrgObjComb[0]->Add(hpp_TrgObj60[0] );
    hpp_TrgObjComb[0]->Add(hpp_TrgObj40[0] );
    hpp_TrgObjComb[0]->Scale(1./4.);

    hpp_TrgObjComb[0]->Print("base");
    hpp_CombJetpT[0]->Print("base");  }

  std::cout<<std::endl<<"/// Job Event-Loop Summary ///"<<std::endl<<std::endl;
  std::cout<<"Total Num of Events in file(s) opened       = " <<h_NEvents->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read from those file(s) = " <<h_NEvents_read->GetEntries()<<std::endl;

  if(!debugMode){
    std::cout<<"Total Num of Events read passing skimCuts                              = " <<
      h_NEvents_skimCut->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events read passing vzCuts and skimCuts (aka good events) = " <<
      h_NEvents_vzCut->GetEntries()<<std::endl<<std::endl;

    std::cout<<"Total Num of good Events w/ jets                             = " <<      
      h_NEvents_withJets->GetEntries()<<std::endl;  //note, approximate, some events make it through with nref=0 for some reason
    std::cout<<"Total Num of good Events, w/ jets, post kmatCuts             = " <<
      h_NEvents_withJets_kmatCut->GetEntries()<<std::endl; 
    if(fillDataJetIDHists) std::cout<<"Total Num of good Events, w/ jets, post kmatCut AND JetIDCut = " <<
      h_NEvents_withJets_JetIDCut->GetEntries()<<std::endl;                                                      }

  std::cout<<std::endl; 
  std::cout<<"Total Num of HLT40 events  = "<<NEvents_40<<std::endl;
  std::cout<<"Total Num of HLT60 events  = "<<NEvents_60<<std::endl;
  std::cout<<"Total Num of HLT80 events  = "<<NEvents_80<<std::endl;
  std::cout<<"Total Num of HLT100 events = "<<NEvents_100<<std::endl<<std::endl;
  

  std::cout<<std::endl<<"/// Job Jet-Loop Summary ///"<<std::endl<<std::endl;
  std::cout<<"Total Num of Jets read from good events                          = " <<
    h_NJets->GetEntries()<<std::endl;  
  std::cout<<"Total Num of Jets read from good events post kinCut              = " <<
    h_NJets_kmatCut->GetEntries()<<std::endl;
  if(fillDataJetIDHists) std::cout<<"Total Num of Jets read from good events post kinCut AND JetIDCut = " <<
    h_NJets_JetIDCut->GetEntries()<<std::endl<<std::endl;
  
  std::cout<<"writing output file..."<<std::endl;
  fout->Write(); 
  
  trgObjpt_40->clear();
  trgObjpt_60->clear();
  trgObjpt_80->clear();
  trgObjpt_100->clear();
  
  std::cout<<std::endl<<"readForests_ppData_jetPlots finished."<<std::endl;  timer.Stop();
  std::cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<std::endl;
  std::cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<std::endl;
  
  return 0;
}    // end readForests_ppData_jetPlots


////// main //////
int main(int argc, char *argv[]){
  
  // error, not enough arguments
  int rStatus = -1;
  if(argc!=8 && argc!=1){
    std::cout<<"for tests on default inputs, do..." <<std::endl;
    std::cout<<"./readForests_ppData_jetPlots.exe";
    std::cout<<std::endl<<std::endl;
    std::cout<<"for actually running, do..."<<std::endl;
    std::cout<<"./readForests_ppData_jetPlots.exe ";
    std::cout<<"<inputFileList> <startFile> <endFile> ";
    std::cout<<"<jetRadius> <jetType> <debugMode> ";
    std::cout<<"<outputFilename> ";
    std::cout<<std::endl<<std::endl;
    std::cout<<"rStatus="<<rStatus<<std::endl;
    return rStatus;
  }
  
  // good input, run
  rStatus=1;
  if(argc==1) rStatus = readForests_ppData_jetPlots();
  else{//read input argument vector
    std::string inputFileList=argv[1]; int startfile= atoi(argv[2]); int endfile= atoi(argv[3]);  
    int jetRadius= atoi(argv[4]); std::string jetType=argv[5];     bool debug=(bool)atoi(argv[6]);
    std::string outputFileName=argv[7];      
    rStatus = readForests_ppData_jetPlots( inputFileList, startfile, endfile, 
				  jetRadius, jetType, debug,
				  outputFileName);
  }
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}

