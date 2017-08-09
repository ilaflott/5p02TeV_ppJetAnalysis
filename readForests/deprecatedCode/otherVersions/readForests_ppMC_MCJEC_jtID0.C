// custom header
#include "readForests.h"

// ppMC switches
const bool fillMCEvtQAHists=true;
const bool fillMCUnfoldingHists=true;
const bool fillMCJECHists=true, fillMCEffHists=true;
const bool fillMCJetIDHists=false, doTightJetID=false;

//// readForests_ppMC_MCJEC
// ---------------------------------------------------------------------------------------------------------------
int readForests_ppMC_MCJEC(std::string inFilelist , int startfile , int endfile ,  
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
  std::cout<<"///////////////////"<<std::endl<<std::endl;
  
  // plot settings for later
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  
  // form jetTreeName
  std::string jetTreeName="ak"+std::to_string(radius)+jetType+"JetAnalyzer/t";
  if(debugMode)std::cout<<"looking at jetTree "<<jetTreeName<<std::endl;
  //std::string CaloJetTreeName="ak"+std::to_string(radius)+"CaloJetAnalyzer/t";
  //if(debugMode)std::cout<<"looking at CalojetTree "<<CaloJetTreeName<<std::endl;
  
  // initialize tree name array
  std::string trees[N_MCTrees];
  trees[0]=jetTreeName; //trees[1]=CaloJetTreeName;
  for(int i=1;i<N_MCTrees;++i)trees[i]=MCTreeNames[i];
  
  // declare TChains for each tree + friend  
  TChain* jetpp[N_MCTrees]={};
  for(int t = 0;t<N_MCTrees;++t)  { 
    jetpp[t] = new TChain( trees[t].data() );
    if(t>0)jetpp[0]->AddFriend( jetpp[t] );  }
  
  
  // open filelist
  std::cout<<"opening filelist: "<<inFilelist<<std::endl;
  std::ifstream instr_Forest(inFilelist.c_str(),std::ifstream::in);
  std::string filename_Forest;  
  std::string lastFileAdded=""; bool filesAdded=false; //Int_t filesCount=0;
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
    for(int t = 0;t<N_MCTrees;++t) filesAdded=jetpp[t]->AddFile(filename_Forest.c_str()); 
    lastFileAdded=filename_Forest;
  }//end file loop
  
  if(debugMode)std::cout<<"filesAdded="<<filesAdded<<std::endl;
  assert(filesAdded);//avoid segfault later
  
  
  // Declare the output file, declare hists after for easy write
  std::cout<<"opening output file "<<outfile<<std::endl;
  TFile *fout = new TFile(outfile.c_str(),"RECREATE");
  fout->cd();
  
  
  // declare hists
  /////   EVT+JET COUNTS AND METADATA   ///// 
  TH1F *hJetPtCut        =new TH1F("hJetPtCut"      ,(std::to_string(jtPtCut)).c_str()   ,    100, 0,100); hJetPtCut->Fill(jtPtCut);           
  TH1F *hJetEtaCut       =new TH1F("hJetEtaCut"     ,(std::to_string(jtEtaCut)).c_str()  ,    5,   0,5  ); hJetEtaCut->Fill(jtEtaCut);	        
  //TH1F *hJetQAPtCut      =new TH1F("hJetQAPtCut"    ,(std::to_string(jetQAPtCut)).c_str(),    100, 0,100); hJetQAPtCut->Fill(jetQAPtCut);     
  //TH1F *hLeadJetPtCut    =new TH1F("hLdJetPtCut"    ,(std::to_string(ldJetPtCut)).c_str(),    100, 0,100); hLeadJetPtCut->Fill(ldJetPtCut);     
  //TH1F *hSubLeadJetPtCut =new TH1F("hSubldJetPtCut" ,(std::to_string(subldJetPtCut)).c_str(), 100, 0,100); hSubLeadJetPtCut->Fill(subldJetPtCut);  
  //TH1F *hPtAveCut        =new TH1F("hPtAveCut"      ,(std::to_string(ptAveCut)).c_str(),      100, 0,100); hPtAveCut->Fill(ptAveCut);  

  TH1F *h_NEvents         = new TH1F("NEvents","NEvents", 1,0,2);
  TH1F *h_NEvents_read    = new TH1F("NEvents_read","NEvents read", 1,0,2);
  TH1F *h_NEvents_skimCut = new TH1F("NEvents_skimCut"      ,"NEvents read post skimCut", 1,0,2);
  TH1F *h_NEvents_vzCut   = new TH1F("NEvents_vzCut"        ,"NEvents read post vzCut AND skimCut", 1,0,2);
  
  //counts post all evt cuts, that are thrown out because of no good jets pre and post jet cuts
  TH1F *h_NEvents_withJets           = new TH1F("NEvents_withJets" , 
		  				"NEvents read post evt cuts, w/ jets", 1,0,2);
  TH1F *h_NEvents_withJets_kmatCut1   = new TH1F("NEvents_withJets_kmatCut1" , 
						"NEvents read post evt cuts, w/ jets post kmatCut1", 1,0,2);
  TH1F *h_NEvents_withJets_kmatCut2   = new TH1F("NEvents_withJets_kmatCut2" , 
						"NEvents read post evt cuts, w/ jets post kmatCut2", 1,0,2);
  
  //jet counts
  TH1F *h_NJets          = new TH1F("NJets","NJets read", 1,0,2);
  TH1F *h_NJets_kmatCut1  = new TH1F("NJets_kmatCut1 ","NJets read post kmatCut1 ", 1,0,2);
  TH1F *h_NJets_kmatCut2  = new TH1F("NJets_kmatCut2 ","NJets read post kmatCut2 ", 1,0,2);
  
  
  /////   EVENT QA   ///// 
  TH1F *hVz=NULL, *hpthatWVz=NULL, *hWVz=NULL ; //*hvzWVz=NULL, 
  TH1F *hpthat=NULL, *hWpthat=NULL;  
  if(fillMCEvtQAHists){
    hVz       = new TH1F("hVz","", 60,-15.,15.);//evtvz
    hpthatWVz = new TH1F("hpthatWeightedVz","", 60,-15.,15.);//pthat-weighted evtvz
    //hvzWVz    = new TH1F("hvzWeightedVz","", 60,-15.,15.);//vz-weighted evtvz
    hWVz      = new TH1F("hWeightedVz","", 60,-15.,15.);//pthat*vz-weighted evt vz
    hpthat    = new TH1F("hpthat","",1000,0,1000);//evt pthat, unweighted and weighted
    hWpthat   = new TH1F("hWeightedpthat","",1000,0,1000);  }
  
  
  /////  UNFOLDING   ///// 
  TH1F *hpp_gen[2]={}, *hpp_reco[2]={}; 
  TH2F* hpp_matrix[2]={};  //gen x reco = matrix
  TH1F *hpp_mcclosure_gen[2]={}, *hpp_mcclosure_data[2]={}, *hpp_mcclosure_data_train[2]={};
  TH2F *hpp_mcclosure_matrix[2]={};  
  if(fillMCUnfoldingHists){
    int jtID=0;
    if(fillMCJetIDHists)jtID=1;
    std::string hUnfTitleArray[]={ "gen", "reco", "matrix", "mcclosure_gen", "mcclosure_data", "mcclosure_data_train", "mcclosure_matrix"};
    const int nUnfTitles=sizeof(hUnfTitleArray)/sizeof(std::string);    
    for(int k=0; k<nUnfTitles; k++){     
      std::string hTitle="hpp_"+hUnfTitleArray[k];
      if(fillMCJetIDHists)hTitle+="_wJetID";
      hTitle+="_R"+std::to_string(radius)+"_"+etaWidth;      
      if(hUnfTitleArray[k]=="gen")	
	hpp_gen[jtID]    = new TH1F( hTitle.c_str(), "reco-matched genpt aka refpt", 1000,0,1000);
      else if(hUnfTitleArray[k]=="reco")	
	hpp_reco[jtID]    = new TH1F( hTitle.c_str(), "recopt", 1000,0,1000);
      else if(hUnfTitleArray[k]=="matrix")	
	hpp_matrix[jtID]    = new TH2F( hTitle.c_str(), "th2f reftpt jtpt", 1000, 0,1000, 1000, 0,1000);
      else if(hUnfTitleArray[k]=="mcclosure_gen")	
	hpp_mcclosure_gen[jtID]    = new TH1F( hTitle.c_str(), "reco-matched genpt aka refpt for mcclosure", 1000,0,1000);
      else if(hUnfTitleArray[k]=="mcclosure_data")	
	hpp_mcclosure_data[jtID]    = new TH1F( hTitle.c_str(), "reco-matched genpt aka refpt for mcclosure", 1000,0,1000);
      else if(hUnfTitleArray[k]=="mcclosure_data_train")	
	hpp_mcclosure_data_train[jtID]    = new TH1F( hTitle.c_str(), "reco-matched genpt aka refpt for mcclosure", 1000,0,1000);
      else if(hUnfTitleArray[k]=="mcclosure_matrix")	
	hpp_mcclosure_matrix[jtID]    = new TH2F( hTitle.c_str(), "th2f reftpt jtpt mcclosure", 1000, 0,1000, 1000, 0,1000);
    }  }
  
  
  /////   GEN/RECO MATCHING   ///// 
  TH2F *hpp_mceff_pt[2]={}, *hpp_mceff_eta[2]={}, *hpp_mceff_phi[2]={};//, *hpp_mceff_drjt[2]={};
  TH2F *hpp_mceff_ptrat_drjt[2]={}, *hpp_mceff_ptrat_eta[2]={}, *hpp_mceff_ptrat_phi[2]={};// all v. genpt/recpt  
  if(fillMCEffHists){
    int jtID=0;
    if(fillMCJetIDHists)jtID=1;    
    std::string hMCEffTitleArray[]={ "pt", "eta", "phi", "ptrat_drjt", "ptrat_eta", "ptrat_phi" };
    const int nMCEffTitles=sizeof(hMCEffTitleArray)/sizeof(std::string);    
    for(int k=0; k<nMCEffTitles; k++){      
      std::string hTitle="hpp_mceff_"+hMCEffTitleArray[k];
      if(fillMCJetIDHists)hTitle+="_wJetID";
      hTitle+="_R"+std::to_string(radius)+"_"+etaWidth;

      std::cout<<"k="<<k<<std::endl;
      std::cout<<"hMCEffTitleArray="<<hMCEffTitleArray[k]<<std::endl;
      
      if(hMCEffTitleArray[k]=="pt")	
	hpp_mceff_pt[jtID] =  new TH2F( hTitle.c_str(), "matched jets, genpt v recpt/genpt",    1000,0,1000, 50,0.,5.);
      else if(hMCEffTitleArray[k]=="eta")	
	hpp_mceff_eta[jtID] = new TH2F( hTitle.c_str(), "matched jets, geneta v receta/geneta", 40,-5.,5., 50,0.,5.);
      else if(hMCEffTitleArray[k]=="phi")	
	hpp_mceff_phi[jtID] = new TH2F( hTitle.c_str(), "matched jets, genphi v recphi/genphi", 40,-4.,4., 50,0.,5.);
      else if(hMCEffTitleArray[k]=="ptrat_drjt")	
	hpp_mceff_ptrat_drjt[jtID] = new TH2F( hTitle.c_str(), "matched jets, gendrjt v recpt/genpt", 10,0.,0.6, 50,0.,5.);   
      else if(hMCEffTitleArray[k]=="ptrat_eta")	
	hpp_mceff_ptrat_eta[jtID]  = new TH2F( hTitle.c_str(), "matched jets, geneta v recpt/genpt", 40,-5.,5., 50,0.,5.);
      else if(hMCEffTitleArray[k]=="ptrat_phi")	
	hpp_mceff_ptrat_phi[jtID]  = new TH2F( hTitle.c_str(), "matched jets, genphi v recpt/genpt", 40,-4.,4., 50,0.,5.);     
    }  }
  
  
  /////   JEC AND JER   /////
  //  TH3F * hJEC[2]={};
  TH1F * hJER_jtetaQA;
  //TH1F * hJEC_check[2][nbins_JEC_ptbins][nbins_eta]={};
  //TH1F * hJER[2][nbins_pt]={};  
  TH1F * hJER[2][nbins_pt2]={};  
  TH1F *hJER_eta_30pt50[2][nbins_eta]={}, *hJER_eta_150pt200[2][nbins_eta]={};
  TH1F *hJER_rapbins[2][nbins_rap]={};
  if(fillMCJECHists){        
    int jtID=0;
    if(fillMCJetIDHists)jtID=1;
    //hJEC[0]= new TH3F( ("hJEC_"+std::to_string(0)+"wJetID").c_str(), ";raw p_{T};#eta;JEC", 500,0.,500., 200,-5.,+5., 300,0.,5.);    
    for(int y = 0; y<nbins_eta; ++y){	               
      std::string hTitleJER_30to50="hJER_"+std::to_string(jtID)+"wJetID_"+
	"etabin"+std::to_string(y)+"_30_pt_50";//"etabin"+std::to_string(etabins[y])+"_30_pt_50";
      hJER_eta_30pt50[jtID][y] = new TH1F(hTitleJER_30to50.c_str(), "recpt/genpt, 30<genpt<50",100, 0, 2);      
      std::string hTitleJER_150to200="hJER_"+std::to_string(jtID)+"wJetID_"+
	"etabin"+std::to_string(y)+"_150_pt_200";
      hJER_eta_150pt200[jtID][y] = new TH1F(hTitleJER_150to200.c_str(), "recpt/genpt, 150<genpt<200",100, 0, 2);          
      /*for(int x = 0; x<nbins_JEC_ptbins; ++x) {
	std::string hTitleJEC_check="hJEC_"+std::to_string(jtID)+"wJetID_check_"+
	  "ptbin"+std::to_string(x)+"_etabin"+std::to_string(y);
	hJEC_check[jtID][x][y] = new TH1F(hTitleJEC_check.c_str(),
					  Form("rawpt/genpt, %2.0f<genpt<%2.0f, %2.4f<geneta<%2.4f", 
					       JEC_ptbins[x],JEC_ptbins[x+1], etabins[y],etabins[y+1]),
					  100,0.,3.);           }*/ }    
    //for(int bin = 0; bin<nbins_pt; ++bin) {      
    for(int bin = 0; bin<nbins_pt2; ++bin) {      
      std::string hTitleJER="hJER_"+std::to_string(jtID)+"wJetID_"+
	std::to_string((int)ptbins2[bin])+"_pt_"+std::to_string((int)ptbins2[bin+1]);      
      //	std::to_string((int)ptbins[bin])+"_pt_"+std::to_string((int)ptbins[bin+1]);      
      hJER[jtID][bin] = new TH1F(hTitleJER.c_str(),"matched jets, recpt/genpt, genptbins", 100,0.,2.); 
    } 
    for(int bin = 0; bin<nbins_rap; ++bin) {      
      std::string hTitleJER="hJER_"+std::to_string(jtID)+"wJetID_jty_bin"+std::to_string(bin);      
      hJER_rapbins[jtID][bin] = new TH1F(hTitleJER.c_str(),"matched jets, recpt/genpt, abs(y) bins", 100,0.,2.); 
    }
        
    hJER_jtetaQA=new TH1F("jtEtaSanityCheckJER", "jetEta Dist Check" , 100, -6., 6.);
    

  }
  
  
  // EVENT LOOP PREP
  // declare variables/arrays + set branch address for each input tree
  //JetAnalyzer, jets
  int nref_I;  float pthat_F;
  float pt_F[1000];  float eta_F[1000];   float phi_F[1000];  
  float y_F[1000]; float rawpt_F[1000];  //float jtpu_F[1000];
  //MC jet variables
  int subid_F[1000];    int refparton_F[1000];
  float refpt_F[1000];  float refeta_F[1000];  float refphi_F[1000];  
  float refdrjt_F[1000];
  //tracks 
  float trkSum_F[1000]; float trkHardSum_F[1000];  float trkMax_F[1000]; 
  //charged particles          
  int chN_I[1000];  
  float chSum_F[1000]; float chHardSum_F[1000];  float chMax_F[1000]; 
  //photons  
  float phSum_F[1000]; float phHardSum_F[1000];  float phMax_F[1000]; 
  //leptons
  float eSum_F[1000];   float eMax_F[1000];        
  float muSum_F[1000];  float muMax_F[1000];    
  //neutrals
  int neN_I[1000];
  float neSum_F[1000];  float neMax_F[1000];        
  
  //jets
  jetpp[0]->SetBranchAddress("nref",&nref_I);
  jetpp[0]->SetBranchAddress("pthat",&pthat_F);
  jetpp[0]->SetBranchAddress("jtpt",&pt_F);
  jetpp[0]->SetBranchAddress("jteta",&eta_F);
  jetpp[0]->SetBranchAddress("jty",&y_F);
  jetpp[0]->SetBranchAddress("jtphi",&phi_F);
  jetpp[0]->SetBranchAddress("rawpt",&rawpt_F);
  //jetpp[0]->SetBranchAddress("jtpu",&jtpu_F);
  //MC jet variables
  jetpp[0]->SetBranchAddress( "subid"	, &subid_F     );
  jetpp[0]->SetBranchAddress( "refparton_flavor" , &refparton_F );
  jetpp[0]->SetBranchAddress( "refpt"	, &refpt_F     );
  jetpp[0]->SetBranchAddress( "refeta"  , &refeta_F    );
  jetpp[0]->SetBranchAddress( "refphi"  , &refphi_F    );
  jetpp[0]->SetBranchAddress( "refdrjt"	, &refdrjt_F   );  
  //charged particles
  jetpp[0]->SetBranchAddress("chargedN",&chN_I);  
  jetpp[0]->SetBranchAddress("chargedSum",&chSum_F);
  jetpp[0]->SetBranchAddress("chargedHardSum",&chHardSum_F);
  jetpp[0]->SetBranchAddress("chargedMax",&chMax_F);  
  //tracks
  jetpp[0]->SetBranchAddress("trackSum",&trkSum_F);  
  jetpp[0]->SetBranchAddress("trackHardSum",&trkHardSum_F);
  jetpp[0]->SetBranchAddress("trackMax",&trkMax_F);  
  //photons
  jetpp[0]->SetBranchAddress("photonSum",&phSum_F);  
  jetpp[0]->SetBranchAddress("photonHardSum",&phHardSum_F);
  jetpp[0]->SetBranchAddress("photonMax",&phMax_F);  
  //leptons
  jetpp[0]->SetBranchAddress("eSum",&eSum_F);  
  jetpp[0]->SetBranchAddress("eMax",&eMax_F);
  jetpp[0]->SetBranchAddress("muSum",&muSum_F);
  jetpp[0]->SetBranchAddress("muMax",&muMax_F);
  //neutrals
  jetpp[0]->SetBranchAddress("neutralN",&neN_I); 
  jetpp[0]->SetBranchAddress("neutralSum",&neSum_F); 
  jetpp[0]->SetBranchAddress("neutralMax",&neMax_F);
  
  // HiEvtAnalyzer
  ULong64_t evt_I;   UInt_t run_I;   UInt_t lumi_I; float vz_F;
  jetpp[1]->SetBranchAddress("evt",&evt_I);
  jetpp[1]->SetBranchAddress("run",&run_I);
  jetpp[1]->SetBranchAddress("lumi",&lumi_I);
  jetpp[1]->SetBranchAddress("vz",&vz_F);
  
  // skimanalysis
  int pBeamScrapingFilter_I, pHBHENoiseFilter_I, pprimaryvertexFilter_I;
  jetpp[2]->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_I);
  jetpp[2]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilter_I);
  jetpp[2]->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryvertexFilter_I);
  
  
  // event count from files
  UInt_t NEvents_jetAnalyzr=jetpp[0]->GetEntries();   // preskim event count from files
  UInt_t NEvents_skimAnalyzr=jetpp[2]->GetEntries();   // preskim event count from files
  bool jetEventsExist= (bool)NEvents_jetAnalyzr, skimEventsExist=(bool)NEvents_skimAnalyzr;
  bool EventCountsEqual= (NEvents_jetAnalyzr==NEvents_skimAnalyzr);
  if( !EventCountsEqual ||
      !jetEventsExist  || !skimEventsExist ){//if this fires, we're gonna segfault
    std::cout<<"jetAnalyzr # evts  = " << NEvents_jetAnalyzr  << std::endl;
    std::cout<<"skimAnalyzr # evts = " << NEvents_skimAnalyzr << std::endl; }
  assert( jetEventsExist && skimEventsExist && EventCountsEqual );//exit if the counts don't make sense, help avoid garbage results
  
  UInt_t NEvents_allFiles=NEvents_jetAnalyzr;   // preskim event count from files
  for(UInt_t i=0;i < NEvents_allFiles; ++i) h_NEvents->Fill(1);
  
  UInt_t NEvents_read=0;
  if(debugMode)NEvents_read = (int)NEvents_allFiles/4 ;
  else NEvents_read = NEvents_allFiles;
  
  std::cout<<"reading "<<NEvents_read<<" events"<<std::endl;   

  float jetIDCut_neSum=0.99, jetIDCut_phSum=0.99;
  if(doTightJetID){     
    jetIDCut_neSum=0.90;  
    jetIDCut_phSum=0.90;  }

  for(UInt_t nEvt = 0; nEvt < NEvents_read; ++nEvt) {//event loop   
    
    if(nEvt%1000==0)std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;  
    jetpp[0]->GetEntry(nEvt);    
    h_NEvents_read->Fill(1);
    
    // skim/HiEvtAnalysis criteria
    if( pHBHENoiseFilter_I==0     || 
        pBeamScrapingFilter_I==0  || 
        pprimaryvertexFilter_I==0  ) continue;    
    h_NEvents_skimCut->Fill(1);    
    if( fabs(vz_F)>15.              ) continue;
    h_NEvents_vzCut->Fill(1);
        
    // compute weights, this way only works because of the previous vz>15. cut
    float vzWeight=1.;
    float vzStart=minbinValue_vzWeights, vzBinLeftSide=vzStart, vzBinRightSide=vzBinLeftSide+binsize_vzWeights;
    for( int i=0; i<nbins_vzWeights ; i++ ) { 
      if(vzBinLeftSide<vz_F && vz_F<=vzBinRightSide) {
	vzWeight=vzWeights[i];  
	break; }       
      vzBinLeftSide+=binsize_vzWeights;
      vzBinRightSide+=binsize_vzWeights;  }    
    float evtPthatWeight=0.;    
    for( int i=0; i<nbins_pthat && pthat_F>=pthatbins[i]; i++ ){ evtPthatWeight=pthatWeights[i]; }     
    float trigWeight=1.;
    float weight_eS=evtPthatWeight*trigWeight*vzWeight;              
    
    //vz hists
    if(fillMCEvtQAHists){
      hVz->Fill(vz_F, 1.);
      //hvzWVz->Fill(vz_F, vzWeight);
      hpthatWVz->Fill(vz_F, evtPthatWeight);
      hWVz->Fill(vz_F, weight_eS);      
      hpthat->Fill(pthat_F, 1.);
      hWpthat->Fill(pthat_F, weight_eS);}
        
    // for event counting + avoiding duplicate fills in dijet hists
    bool hNEvts_withJets_Filled=false;
    bool hNEvts_withJets_kmatCut1_Filled=false, hNEvts_withJets_kmatCut2_Filled=false; 
    for(int jet = 0; jet<nref_I; ++jet){

      // event+jet counting
      h_NJets->Fill(1);
      if(!hNEvts_withJets_Filled){
	h_NEvents_withJets->Fill(1);
	hNEvts_withJets_Filled=true;      }      
      
      //if reco jet w/o matched gen jet, skip.
      if( subid_F[jet]!=0 ) continue;
      
      float recpt  = pt_F[jet];
      float receta = eta_F[jet];
      float rawpt  = rawpt_F[jet];
      float recy = y_F[jet];
      float recphi = phi_F[jet];

      if( recpt<jtPtCut   ) continue;                 
      if( fabs(recy) > jtEtaCut) continue;

      // 13 TeV JetID criterion
      bool passesJetID=false; int jtID=0;
      if(fillMCJetIDHists) {	

	//if (fabs(recy)>3.0)passesJetID=false;
	if (fabs(recy)>3.0){
	  if(phSum_F[jet]/rawpt < 0.90 &&
	     neN_I[jet]>10) 
	    passesJetID=true;	}
	else if ( 2.7<fabs(recy) && fabs(recy)<=3.0  ){
	  if( phSum_F[jet]/rawpt < 0.90 && //neutral em
	      neN_I[jet] > 2  )
	    passesJetID=true;	} 	
	else if ( 2.4<fabs(recy) && fabs(recy)<=2.7){
	  if( neSum_F[jet]/rawpt    < jetIDCut_neSum &&  //neutral had 
	      phSum_F[jet]/rawpt    < jetIDCut_phSum && 
	      chN_I[jet]+neN_I[jet] > 1 ) 
	    passesJetID=true; 	  }	
	else { //if(fabs(recy)<=2.4) //in the barrel, strictest
	    if( neSum_F[jet]/rawpt    < jetIDCut_neSum &&
		phSum_F[jet]/rawpt    < jetIDCut_phSum &&
		chN_I[jet]+neN_I[jet] > 1    &&
		chSum_F[jet]/rawpt    > 0.   && //charged had 
		eSum_F[jet]/rawpt     < 0.99 && 
		chN_I[jet] > 0 )      	  
	      passesJetID=true;  }      

	if(passesJetID)jtID=1;
	else continue;      //continue because we do jtIDT, jtIDL, and jtID0 all in different job
      }      
      
      // jet/event counts
      h_NJets_kmatCut1->Fill(1);
      if(!hNEvts_withJets_kmatCut1_Filled){
	h_NEvents_withJets_kmatCut1->Fill(1);
	hNEvts_withJets_kmatCut1_Filled=true;      }

      float genpt= refpt_F[jet];
      float geneta = refeta_F[jet];
      float genphi = refphi_F[jet];
      float gendrjt = refdrjt_F[jet];     
      
      
      /////   JEC AND JER   ///// 
      if(fillMCJECHists ){
	
	// sanity check
	hJER_jtetaQA->Fill(receta,weight_eS);      
	//hJEC[jtID]->Fill( rawpt,eta_F[jet],(float)(recpt/rawpt), weight_eS );

	// dual-diff xsec bins
	int absetabin=-1;
	for(int bin=0; bin < nbins_rap; bin++){
	  if( fabs(geneta)>rapbins[bin] ) absetabin=bin; 	}
	if(fabs(geneta)>3.7)absetabin=-1;
	if(absetabin != -1 ){
	  hJER_rapbins[jtID][absetabin]->Fill( (float)(recpt/genpt), weight_eS);     }       

	// gen pt binnings
	int ptbin = -1; 
	//for(int bin = 0; bin<nbins_pt; ++bin) { // binning inherited from raghav
	//	  if(genpt > ptbins[bin]) ptbin = bin;}	
	for(int bin = 0; bin<nbins_pt2; ++bin)  { // binning recommended by raghav for 5TeVppJets
	  if(genpt > ptbins2[bin]) ptbin = bin; }	
	if(ptbin != -1) 
	  hJER[jtID][ptbin]->Fill( (float)(recpt/genpt), weight_eS);   	

	int etabin = -1;
	for(int bin = 0; bin<nbins_eta; ++bin)    {
	  if(geneta > etabins[bin]) etabin = bin; }
	if(geneta>5.191) etabin=-1;	
	if(etabin != -1) {
	  if(genpt >= 30  && genpt<50 ) 
	    hJER_eta_30pt50[jtID][etabin]->Fill((float)recpt/genpt, weight_eS);
	  else if(genpt >= 150 && genpt<200) 
	    hJER_eta_150pt200[jtID][etabin]->Fill((float)recpt/genpt, weight_eS);	

	  /*int binx = -1; //genpt JEC bin
	  for(int bin = 0; bin<nbins_JEC_ptbins; ++bin) 
	    if(genpt > JEC_ptbins[bin]) binx = bin;      	
	  //std::cout <<"binx="<<binx<<std::endl;	    	  
	  if(binx != -1)
	    hJEC_check[jtID][binx][etabin]->Fill( (float)(rawpt/genpt) , weight_eS); 	*/   
	} }
      
      /////   GEN/RECO MATCHING   /////
      if(fillMCEffHists){	
	hpp_mceff_pt[jtID]->Fill(genpt, recpt/genpt, weight_eS);
	hpp_mceff_eta[jtID]->Fill(geneta, receta/geneta, weight_eS);
	hpp_mceff_phi[jtID]->Fill(genphi, recphi/genphi, weight_eS);	
	hpp_mceff_ptrat_drjt[jtID]->Fill(gendrjt, recpt/genpt, weight_eS);	
	hpp_mceff_ptrat_eta[jtID]->Fill(geneta, recpt/genpt, weight_eS);
	hpp_mceff_ptrat_phi[jtID]->Fill(genphi, recpt/genpt, weight_eS);       }

      // jet/event counts
      h_NJets_kmatCut2->Fill(1);
      if(!hNEvts_withJets_kmatCut2_Filled){
	h_NEvents_withJets_kmatCut2->Fill(1);
	hNEvts_withJets_kmatCut2_Filled=true;      }
      
      /////   UNFOLDING   ///// 
      if(fillMCUnfoldingHists){
	hpp_gen[jtID]->Fill(genpt, weight_eS);
	hpp_reco[jtID]->Fill(recpt, weight_eS);
	hpp_matrix[jtID]->Fill(recpt, genpt, weight_eS);
	if(nEvt%2 == 0){
	  hpp_mcclosure_data[jtID]->Fill(recpt, weight_eS);      }
	else {
	  hpp_mcclosure_matrix[jtID]->Fill(recpt, genpt, weight_eS);
	  hpp_mcclosure_gen[jtID]->Fill(genpt, weight_eS);
	  hpp_mcclosure_data_train[jtID]->Fill(recpt, weight_eS); 	} 
      }      
    
      
    }//end jet loop
  }//end event loop
  
  std::cout<<std::endl;
  
  std::cout<<std::endl<<"/// Job Event-Loop Summary ///"<<std::endl<<std::endl;
  std::cout<<"Total Num of Events in file(s) opened       = " <<h_NEvents->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read from those file(s) = " <<h_NEvents_read->GetEntries()<<std::endl;
  
  
  std::cout<<"Total Num of Events read passing skimCuts                              = " <<
    h_NEvents_skimCut->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read passing vzCuts and skimCuts (aka good events) = " <<
    h_NEvents_vzCut->GetEntries()<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Total Num of good Events w/ jets                             = " <<
    h_NEvents_withJets->GetEntries()<<std::endl;	   
  std::cout<<"Total Num of good Events, w/ jets, post kmatCut1             = " <<
    h_NEvents_withJets_kmatCut1->GetEntries()<<std::endl; 
  std::cout<<"Total Num of good Events, w/ jets, post kmatCut2             = " <<
    h_NEvents_withJets_kmatCut2->GetEntries()<<std::endl ;
  std::cout<<std::endl;
  
  std::cout<<std::endl<<"/// Job Jet-Loop Summary ///"<<std::endl<<std::endl;
  std::cout<<"Total Num of Jets read from good events                          = " <<
    h_NJets->GetEntries()<<std::endl;  
  std::cout<<"Total Num of Jets read from good events post kmatCut1              = " <<
    h_NJets_kmatCut1->GetEntries()<<std::endl;
  std::cout<<"Total Num of Jets read from good events post kmatCut2              = " <<
    h_NJets_kmatCut2->GetEntries()<<std::endl;
  std::cout<<std::endl;
  

  std::cout<<"writing output file "<<outfile<<std::endl;
  fout->Write(); 


  std::cout<<std::endl<<"readForests_ppMC_MCJEC finished."<<std::endl;  timer.Stop();
  std::cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<std::endl;
  std::cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<std::endl;

return 0;
}// end readForests_ppMC_MCJEC
								 

////// main //////
// acts as the frontend control for .exe file
int main(int argc, char *argv[]){
  
  // error, not enough arguments
  int rStatus = -1;
  if(argc!=8 && argc!=1){
    std::cout<<"for tests on default inputs, do..." <<std::endl;
    std::cout<<"./readForests_ppMC_MCJEC.exe";
    std::cout<<std::endl<<std::endl;
    std::cout<<"for actually running, do..."<<std::endl;
    std::cout<<"./readForests_ppMC_MCJEC.exe ";
    std::cout<<"<inputFileList> <startFile> <endFile> ";
    std::cout<<"<jetRadius> <jetType> <debugMode> ";
    std::cout<<"<outputFilename> ";
    std::cout<<std::endl<<std::endl;
    std::cout<<"rStatus="<<rStatus<<std::endl;
    return rStatus;
  }
  
  // good input, run
  rStatus=1;
  if(argc==1) rStatus = readForests_ppMC_MCJEC();
  else{//read input argument vector
    std::string inputFileList=argv[1]; int startfile= atoi(argv[2]); int endfile= atoi(argv[3]);  
    int jetRadius= atoi(argv[4]); std::string jetType=argv[5];     bool debug=(bool)atoi(argv[6]);
    std::string outputFileName=argv[7];      
    rStatus = readForests_ppMC_MCJEC( inputFileList, startfile, endfile, 
				jetRadius, jetType, debug,
				outputFileName);  }
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}
