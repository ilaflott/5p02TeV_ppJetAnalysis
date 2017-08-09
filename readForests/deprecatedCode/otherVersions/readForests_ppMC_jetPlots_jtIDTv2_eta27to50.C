// custom header
#include "readForests.h"

// ppMC switches
const bool fillMCEvtQAHists=true, fillMCJetQAHists=true, fillMCJetIDHists=true; //most basic-level plots
const bool fillMCJetSpectraRapHists=true; //other
const bool tightJetID=true;

//// readForests_ppMC_jetPlots
// ---------------------------------------------------------------------------------------------------------------
int readForests_ppMC_jetPlots(std::string inFilelist , int startfile , int endfile ,  
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
  TH1F *hJetQAPtCut      =new TH1F("hJetQAPtCut"    ,(std::to_string(jetQAPtCut)).c_str(),    100, 0,100); hJetQAPtCut->Fill(jetQAPtCut);     
  TH1F *hLeadJetPtCut    =new TH1F("hLdJetPtCut"    ,(std::to_string(ldJetPtCut)).c_str(),    100, 0,100); hLeadJetPtCut->Fill(ldJetPtCut);     
  TH1F *hSubLeadJetPtCut =new TH1F("hSubldJetPtCut" ,(std::to_string(subldJetPtCut)).c_str(), 100, 0,100); hSubLeadJetPtCut->Fill(subldJetPtCut);  
  TH1F *hPtAveCut        =new TH1F("hPtAveCut"      ,(std::to_string(ptAveCut)).c_str(),      100, 0,100); hPtAveCut->Fill(ptAveCut);  
  
  TH1F *h_NEvents         = new TH1F("NEvents","NEvents", 1,0,2);
  TH1F *h_NEvents_read    = new TH1F("NEvents_read","NEvents read", 1,0,2);
  TH1F *h_NEvents_skimCut = new TH1F("NEvents_skimCut"      ,"NEvents read post skimCut", 1,0,2);
  TH1F *h_NEvents_vzCut   = new TH1F("NEvents_vzCut"        ,"NEvents read post vzCut AND skimCut", 1,0,2);
  
  //counts post all evt cuts, that are thrown out because of no good jets pre and post jet cuts
  TH1F *h_NEvents_withJets           = new TH1F("NEvents_withJets" , 
		  				"NEvents read post evt cuts, w/ jets", 1,0,2);
  TH1F *h_NEvents_withJets_kmatCut   = new TH1F("NEvents_withJets_kmatCut" , 
						"NEvents read post evt cuts, w/ jets post kmatCut", 1,0,2);
  TH1F *h_NEvents_withJets_JetIDCut = NULL;
  if(fillMCJetIDHists)
    h_NEvents_withJets_JetIDCut = new TH1F("NEvents_withJets_JetIDCut" , 
					   "NEvents read post evt cuts, w/ jets post kmatCut AND JetID Cut", 1,0,2);
  
  //jet counts
  TH1F *h_NJets          = new TH1F("NJets","NJets read", 1,0,2);
  TH1F *h_NJets_kmatCut  = new TH1F("NJets_kmatCut ","NJets read post kmatCut ", 1,0,2);
  TH1F *h_NJets_JetIDCut = NULL;
  if(fillMCJetIDHists)
    h_NJets_JetIDCut = new TH1F("NJets_JetIDCut","NJets read post JetIDCut AND kmatCut", 1,0,2);
  
  
  /////   EVENT QA HISTS   ///// 
  TH1F *hVz=NULL, *hpthatWVz=NULL, *hWVz=NULL ; //*hvzWVz=NULL, 
  TH1F *hpthat=NULL, *hWpthat=NULL;  
  if(fillMCEvtQAHists){
    hVz       = new TH1F("hVz","", 60,-15.,15.);//evtvz
    hpthatWVz = new TH1F("hpthatWeightedVz","", 60,-15.,15.);//pthat-weighted evtvz
    //hvzWVz    = new TH1F("hvzWeightedVz","", 60,-15.,15.);//vz-weighted evtvz
    hWVz      = new TH1F("hWeightedVz","", 60,-15.,15.);//pthat*vz-weighted evt vz
    hpthat    = new TH1F("hpthat","",1000,0,1000);//evt pthat, unweighted and weighted
    hWpthat   = new TH1F("hWeightedpthat","",1000,0,1000);  }
  
  
  /////   JET QA HISTS   ///// 
  //jet QA
  TH1F *hJetQA[2][N_vars]={};
  //MC JetQA
  bool fillgenJetQA=true&&fillMCJetQAHists;
  TH1F *hMCJetQA_genpt[2]={}, *hMCJetQA_geneta[2]={};
  TH1F *hMCJetQA_genrecpt[2]={}, *hMCJetQA_genreceta[2]={};
  if(fillMCJetQAHists){
    for(int k = 0; k<2; ++k){
      if(!fillMCJetIDHists && k==1)continue;
      //if( fillMCJetIDHists && k==0)continue;	
      for(int j = 0; j<N_vars; ++j){	
	//binnings+hist declarations, etc.
	if(var[j]=="jtpt"||var[j]=="rawpt") 
	  hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;", var[j].c_str()) , 500,0,500);       
	else if(var[j]=="jty"||var[j]=="jteta")  //rapidity + pseudorapidity
	  hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;",var[j].c_str()) , 100,-5,+5);
	else if(var[j]=="jtphi") 
	  hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;",var[j].c_str()) , 100,-4,+4);
	else if(var[j]=="dphi") 
	  hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;",var[j].c_str()) , 50,0,+4);
	else if(var[j]=="leadJetPt"||var[j]=="subleadJetPt") 
	  hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;", var[j].c_str()), 500,0,500);
	else if(var[j]=="neN"||var[j]=="sumN"||var[j]=="chN")
	  hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 30,0,30);         
	else //consituent binnings
	  hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 200,0,2); }
      
      if(fillgenJetQA){	
	hMCJetQA_genpt[k]    = new TH1F(Form("hJetQA_%dwJetID_genpt",k)     , "genpt",    1000, 0,1000 );
	hMCJetQA_geneta[k]   = new TH1F(Form("hJetQA_%dwJetID_geneta",k)    , "geneta", 200,  -6.,6. );
	hMCJetQA_genrecpt[k] = new TH1F(Form("hJetQA_%dwJetID_genrecpt",k)  , "genrecpt",    100, 0, 2. );
	hMCJetQA_genreceta[k]= new TH1F(Form("hJetQA_%dwJetID_genreceta",k) , "genreceta", 50,  -0.1, 2. ); }
      
    }}
  
  bool fillgenJetRapHists=true&&fillMCJetSpectraRapHists;  
  TH1F *hJetSpectraRap[2][nbins_rap]={};
  TH1F *hMCJetQA_rapBins_genpt[2][nbins_rap]={};
  TH1F *hMCJetQA_rapBins_geneta[2][nbins_rap]={};
  TH1F *hMCJetQA_rapBins_genrecpt[2][nbins_rap]={};
  TH1F *hMCJetQA_rapBins_genreceta[2][nbins_rap]={};
  if(fillMCJetSpectraRapHists){
    for(int k = 0; k<2; ++k){
      if(!fillMCJetIDHists && k==1)continue;	
      //      if( fillMCJetIDHists && k==0)continue;
      for(int j = 0; j<nbins_rap; ++j){
	std::string h_Title="hJetSpectraRap_";
	if(k==1)h_Title+="wJetID_";
	h_Title+="bin"+std::to_string(j);
	std::stringstream stream1, stream2;	
	stream1.precision(1); stream1 << std::fixed << rapbins[j];
	stream2.precision(1); stream2 << std::fixed << rapbins[j+1];
	std::string h_Desc="JetPt Spectra for "+stream1.str()+"<abs(y)<"+ stream2.str();	
	hJetSpectraRap[k][j]=new TH1F(h_Title.c_str(),h_Desc.c_str(), 2000,0,2000);  
	
	if(fillgenJetRapHists){
	  hMCJetQA_rapBins_genpt[k][j]     = new TH1F ( Form("hJetSpectraRap_%dwJetID_bin%d_genpt",k,j), "genpt", 1000,0,1000);
	  hMCJetQA_rapBins_geneta[k][j]    = new TH1F ( Form("hJetSpectraRap_%dwJetID_bin%d_geneta",k,j), "geneta", 200,-6.,6.);
	  hMCJetQA_rapBins_genrecpt[k][j]  = new TH1F ( Form("hJetSpectraRap_%dwJetID_bin%d_genrecpt",k,j),  "rec/gen,pt", 100, 0.,2.);
	  hMCJetQA_rapBins_genreceta[k][j] = new TH1F ( Form("hJetSpectraRap_%dwJetID_bin%d_genreceta",k,j), "rec/gen,eta", 105, -0.1,2.0);
	}
	
      }}}
  
  
  
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
  //if(fastDebugMode&&debugMode) NEvents_read = 1000*(endfile-startfile+1); 
  //else NEvents_read = NEvents_allFiles;
  NEvents_read = NEvents_allFiles;
  
  std::cout<<"reading "<<NEvents_read<<" events"<<std::endl;   
  
  float jetIDCut_neSum, jetIDCut_phSum;
  if(tightJetID){     jetIDCut_neSum=0.90;  jetIDCut_phSum=0.90;}
  else{     jetIDCut_neSum=0.99;  jetIDCut_phSum=0.99;}
  for(UInt_t nEvt = 0; nEvt < NEvents_read; ++nEvt) {//event loop   
    
    if( fastDebugMode &&
        debugMode     &&
        nEvt%1000==0     )std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;
    else if(nEvt%10000==0)std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;  
    
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
      else {
	vzBinLeftSide+=binsize_vzWeights;
	vzBinRightSide+=binsize_vzWeights; } }
    
    float evtPthatWeight=0.;    
    for( int i=0; i<nbins_pthat && pthat_F>=pthatbins[i]; i++ ){ evtPthatWeight=pthatWeights[i]; } 
    
    float trigWeight=1.;
    //trigWeight = trigComb(trgDec, treePrescl, triggerPt);    
    
    float weight_eS=evtPthatWeight*trigWeight*vzWeight;              
    
    //vz
    if(fillMCEvtQAHists){
      hVz->Fill(vz_F, 1.);
      //hvzWVz->Fill(vz_F, vzWeight);
      hpthatWVz->Fill(vz_F, evtPthatWeight);
      hWVz->Fill(vz_F, weight_eS);
      
      hpthat->Fill(pthat_F, 1.);
      hWpthat->Fill(pthat_F, weight_eS);}
    
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
    
    for(int jet = 0; jet<nref_I; ++jet)      {
      
      // event+jet counting
      h_NJets->Fill(1);
      if(!hNEvts_withJets_Filled){
	h_NEvents_withJets->Fill(1);
	hNEvts_withJets_Filled=true;      }      
      
      float recpt  = pt_F[jet];
      float receta = eta_F[jet];
      float genpt  = refpt_F[jet];
      float geneta = refeta_F[jet];
      
      //if reco jet w/o matched gen jet, skip.
      if( subid_F[jet]!=0 ) continue;
      if( recpt<jtPtCut   ) continue;           
      
      float rawpt  = rawpt_F[jet];      
      float recy = y_F[jet];
      float recphi = phi_F[jet];
      
      // 13 TeV JetID criterion, loose or tight
      // for loose; ne/phSum fractions = 0.99 for abs(recy)<2.7
      // for tight; ne/phSum fractions = 0.90 for abs(recy)<2.7
      bool passesJetID=false;
      
      if(fillMCJetIDHists) {		
	//	if( fabs(recy)>3.0 ) 	  {
	//	  if( phSum_F[jet]/rawpt < 0.90 &&
	//	      neN_I[jet] > 10) { 
	//	    passesJetID=true; } 
	//	}
	//	else if ( 2.7<fabs(recy) && fabs(recy)<=3.0  )  {
	//	  if( phSum_F[jet]/rawpt < 0.90 && //neutral em
	//	      neN_I[jet] > 2  )
	//	    passesJetID=true;	
	//	} 	
	if( fabs(recy)>2.7) {//jetIDv2
	  if( phSum_F[jet]/rawpt<0.90 &&
	      neSum_F[jet]/rawpt<0.99 &&
	      ((phSum_F[jet]/rawpt>0.) || (neSum_F[jet]/rawpt>0.)) )
	    passesJetID=true;
	}
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
      if(  fillMCJetSpectraRapHists ) { 
	int theRapBin=-1;
	for(int rapbin=0;rapbin<nbins_rap;++rapbin)
	  if( rapbins[rapbin]<=fabs(recy)  && 		
	      fabs(recy)<rapbins[rapbin+1]    ) {
	    
	    theRapBin=rapbin;
	    hJetSpectraRap[0][theRapBin]->Fill(recpt,weight_eS);  
	    
	    if(fillgenJetRapHists){
	      hMCJetQA_rapBins_genpt[0][theRapBin]     ->Fill(genpt,weight_eS);
	      hMCJetQA_rapBins_geneta[0][theRapBin]    ->Fill(geneta,weight_eS);
	      hMCJetQA_rapBins_genrecpt[0][theRapBin]  ->Fill(recpt/genpt,weight_eS);
	      hMCJetQA_rapBins_genreceta[0][theRapBin] ->Fill(receta/geneta,weight_eS);	    }
	    
	    if( passesJetID ){
	      hJetSpectraRap[1][theRapBin]->Fill(recpt,weight_eS); 
	      if(fillgenJetRapHists){
		hMCJetQA_rapBins_genpt[1][theRapBin]     ->Fill(genpt,weight_eS);
		hMCJetQA_rapBins_geneta[1][theRapBin]    ->Fill(geneta,weight_eS);
		hMCJetQA_rapBins_genrecpt[1][theRapBin]  ->Fill(recpt/genpt,weight_eS);
		hMCJetQA_rapBins_genreceta[1][theRapBin] ->Fill(receta/geneta,weight_eS);	    
	      }	    }	    
	    break; 
	  }      }
      
      //second half of kmat cut
      if( fabs(receta)<=2.7 ) continue;
      else if( fabs(receta)>(5.0) ) continue;
      //if( fabs(receta)>jtEtaCut ) continue;

      // jet/event counts
      h_NJets_kmatCut->Fill(1);
      if(!hNEvts_withJets_kmatCut_Filled){
	h_NEvents_withJets_kmatCut->Fill(1);
	hNEvts_withJets_kmatCut_Filled=true;      }

      /////   JETQA   ///// 
      if(fillMCJetQAHists){
	
	if(fillgenJetQA){
	  hMCJetQA_genpt[0]    ->Fill(genpt,weight_eS);
	  hMCJetQA_geneta[0]   ->Fill(geneta,weight_eS);
	  hMCJetQA_genrecpt[0] ->Fill(recpt/genpt,weight_eS);
	  hMCJetQA_genreceta[0]->Fill(receta/geneta,weight_eS);	}

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
	    firstGoodJetPhi=recphi; } }
	
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
		secondGoodJetPhi = recphi; }}} }
	
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
      }//end fillMCJetQAHists	      
      

      if(fillMCJetIDHists){

	if(!passesJetID) continue;
	
	/////   JETQA   /////
	if(fillMCJetQAHists){

	  if(fillgenJetQA){
	    hMCJetQA_genpt[1]    ->Fill(genpt,weight_eS);
	    hMCJetQA_geneta[1]   ->Fill(geneta,weight_eS);
	    hMCJetQA_genrecpt[1] ->Fill(recpt/genpt,weight_eS);
	    hMCJetQA_genreceta[1]->Fill(receta/geneta,weight_eS);	}	  

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
	    hJetQA[1][22]->Fill(chN_I[jet]+neN_I[jet], weight_eS); }
	  
	  
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
	}//end fillMCJetQAHists
	
	// jet/event counts
	h_NJets_JetIDCut->Fill(1);
	if(!hNEvts_withJets_JetIDCut_Filled){
	  h_NEvents_withJets_JetIDCut->Fill(1);
	  hNEvts_withJets_JetIDCut_Filled=true;      }
      }//end doJetIDHists
      
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
  std::cout<<"Total Num of good Events, w/ jets, post kmatCuts             = " <<
    h_NEvents_withJets_kmatCut->GetEntries()<<std::endl; 
  if(fillMCJetIDHists)std::cout<<"Total Num of good Events, w/ jets, post kmatCut AND JetIDCut = " <<
    h_NEvents_withJets_JetIDCut->GetEntries()<<std::endl;
  std::cout<<std::endl;
  
  std::cout<<std::endl<<"/// Job Jet-Loop Summary ///"<<std::endl<<std::endl;
  std::cout<<"Total Num of Jets read from good events                          = " <<
    h_NJets->GetEntries()<<std::endl;  
  std::cout<<"Total Num of Jets read from good events post kinCut              = " <<
    h_NJets_kmatCut->GetEntries()<<std::endl;
  if(fillMCJetIDHists)std::cout<<"Total Num of Jets read from good events post kinCut AND JetIDCut = " <<
    h_NJets_JetIDCut->GetEntries()<<std::endl;
  std::cout<<std::endl;


  std::cout<<"writing output file "<<outfile<<std::endl;
  fout->Write(); 


  std::cout<<std::endl<<"readForests_ppMC_jetPlots finished."<<std::endl;  timer.Stop();
  std::cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<std::endl;
  std::cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<std::endl;

return 0;
}// end readForests_ppMC_jetPlots
								 

////// main //////
// acts as the frontend control for .exe file
int main(int argc, char *argv[]){
  
  // error, not enough arguments
  int rStatus = -1;
  if(argc!=8 && argc!=1){
    std::cout<<"for tests on default inputs, do..." <<std::endl;
    std::cout<<"./readForests_ppMC_jetPlots.exe";
    std::cout<<std::endl<<std::endl;
    std::cout<<"for actually running, do..."<<std::endl;
    std::cout<<"./readForests_ppMC_jetPlots.exe ";
    std::cout<<"<inputFileList> <startFile> <endFile> ";
    std::cout<<"<jetRadius> <jetType> <debugMode> ";
    std::cout<<"<outputFilename> ";
    std::cout<<std::endl<<std::endl;
    std::cout<<"rStatus="<<rStatus<<std::endl;
    return rStatus;
  }
  
  // good input, run
  rStatus=1;
  if(argc==1) rStatus = readForests_ppMC_jetPlots();
  else{//read input argument vector
    std::string inputFileList=argv[1]; int startfile= atoi(argv[2]); int endfile= atoi(argv[3]);  
    int jetRadius= atoi(argv[4]); std::string jetType=argv[5];     bool debug=(bool)atoi(argv[6]);
    std::string outputFileName=argv[7];      
    rStatus = readForests_ppMC_jetPlots( inputFileList, startfile, endfile, 
				jetRadius, jetType, debug,
				outputFileName);  }
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}
