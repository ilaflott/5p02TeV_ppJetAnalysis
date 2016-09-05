// custom header
#include "readForests_ppData.h"

//// readForests_ppData
// ---------------------------------------------------------------------------------------------------------------
int readForests_ppData(int startfile , int endfile , std::string inFilelist , std::string outfile ,
		     int radius , std::string jetType , bool debugMode ){ 


  // for monitoring performance + debugging
  TStopwatch timer;  timer.Start();
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl;
  const bool fastDebugMode = (debugMode)&&true; //if debugMode is off, fastDebugMode shouldn't be on
  if(fastDebugMode)std::cout<<"fastDebugMode is ON"<<std::endl;

  // basic info the screen
  std::cout<<std::endl<<"///////////////////"<<std::endl;
  std::cout<<"reading files "<< startfile << " to " << endfile<<std::endl;
  std::cout<<"radius =  " << radius<<std::endl;
  std::cout<<"jetType =  " << jetType<<std::endl;
  if(debugMode)std::cout<<"N trees = " << N_trees<<std::endl;
  if(debugMode)std::cout<<"N vars  = " << N_vars<<std::endl;
  std::cout<<"///////////////////"<<std::endl<<std::endl;

  // plot settings for later
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  

  // form jetTreeName
  std::string jetTreeName="ak"+std::to_string(radius)+jetType+"JetAnalyzer/t";
  if(debugMode)std::cout<<"looking at jetTree "<<jetTreeName<<std::endl;

  // initialize tree name array
  std::string trees[N_trees];
  trees[0]=jetTreeName;
  for(int i=1;i<N_trees;++i)trees[i]=treeNames[i];

  // declare TChains for each tree + friend them
  TChain* jetpp[N_trees];
  for(int t = 0;t<N_trees;++t)  { 
    jetpp[t] = new TChain( trees[t].data() );
    if(t>0)jetpp[0]->AddFriend( jetpp[t] );  }

    
  // open filelist
  std::cout<<"opening filelist: "<<inFilelist<<std::endl;
  std::ifstream instr_Forest(inFilelist.c_str(),std::ifstream::in);
  std::string filename_Forest;  

  // add files, including startfile+endfile, to chains
  std::string lastFileAdded=""; bool treesAdded=false;
  for(int ifile = 0; ifile<=endfile; ++ifile){

    // grab a filename, check startfile and end of filelist condition
    instr_Forest>>filename_Forest; 
    if(ifile<startfile){ lastFileAdded=filename_Forest; 
      continue; }
    if(filename_Forest==lastFileAdded){ std::cout<<"end of filelist!"<<std::endl; 
      break; }

    std::cout<<"adding file #"<<ifile; 
    if(debugMode)std::cout<<", "<<filename_Forest; 

    //std::cout<<", to TChain jetpp[0]"<<std::endl;
    //jetpp[0]->Add(filename_Forest.c_str());
    //if(!treesAdded)treesAdded=true;     
    //lastFileAdded=filename_Forest;

    std::cout<<", to each TChain in array"<<std::endl;
    for(int t = 0;t<N_trees;++t) jetpp[t]->Add(filename_Forest.c_str()); 
    if(!treesAdded)treesAdded=true;     
    lastFileAdded=filename_Forest;
  }//end file loop
  assert(treesAdded);//avoid segfault later


  // declare variables/arrays + set branch address for each input tree
  //JetAnalyzer, jets
  int nref_I;  
  float pt_F[1000];  
  float eta_F[1000];   float phi_F[1000];  
  float rawpt_F[1000]; float jtpu_F[1000];
  //tracks 
  float trkSum_F[1000]; float trkHardSum_F[1000];  
  float trkMax_F[1000]; //float trkHardMax_F[1000];
  //charged particles          
  int chN_I[1000];  
  float chSum_F[1000]; float chHardSum_F[1000];  
  float chMax_F[1000]; //float chHardMax_F[1000]; 
  //photons  
  float phSum_F[1000]; float phHardSum_F[1000];   
  float phMax_F[1000]; //float phHardMax_F[1000];          
  //leptons
  float eSum_F[1000];   float eMax_F[1000];        
  float muSum_F[1000];  float muMax_F[1000];    
  //neutrals
  float neSum_F[1000];  float neMax_F[1000];        

  //jets
  jetpp[0]->SetBranchAddress("nref",&nref_I);
  jetpp[0]->SetBranchAddress("jtpt",&pt_F);
  jetpp[0]->SetBranchAddress("jteta",&eta_F);
  jetpp[0]->SetBranchAddress("jtphi",&phi_F);
  jetpp[0]->SetBranchAddress("rawpt",&rawpt_F);
  jetpp[0]->SetBranchAddress("jtpu",&jtpu_F);
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
  jetpp[0]->SetBranchAddress("neutralSum",&neSum_F); 
  jetpp[0]->SetBranchAddress("neutralMax",&neMax_F);


  // HiEvtAnalyzer
  ULong64_t evt_I;   UInt_t run_I;   UInt_t lumi_I;   
  float vz_F;

  jetpp[1]->SetBranchAddress("evt",&evt_I);
  jetpp[1]->SetBranchAddress("run",&run_I);
  jetpp[1]->SetBranchAddress("lumi",&lumi_I);
  jetpp[1]->SetBranchAddress("vz",&vz_F);


  // skimanalysis
  int pBeamScrapingFilter_I, pHBHENoiseFilter_I, pprimaryvertexFilter_I;
  jetpp[2]->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_I);
  jetpp[2]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilter_I);
  jetpp[2]->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryvertexFilter_I);


  // hltanalysis
  int jet40_l1s_ps_I, jet60_l1s_ps_I,  jet80_l1s_ps_I,  jet100_l1s_ps_I;   //L1
  int jet40_I, jet60_I, jet80_I, jet100_I;   //HLT
  int jet40_p_I, jet60_p_I, jet80_p_I, jet100_p_I;

  std::string HLTBranches[N_HLTBits]; 
  for(int i=0;i<N_L1Bits;i++) HLTBranches[i]=HLTBitStrings[i]+"_v1";
  jetpp[3]->SetBranchAddress( HLTBranches[0].c_str() , &jet40_I);
  jetpp[3]->SetBranchAddress( HLTBranches[1].c_str() , &jet60_I);
  jetpp[3]->SetBranchAddress( HLTBranches[2].c_str() , &jet80_I);
  jetpp[3]->SetBranchAddress( HLTBranches[3].c_str() , &jet100_I);
  
  std::string HLTPresclBranches[N_HLTBits]; 
  for(int i=0;i<N_HLTBits;i++) HLTPresclBranches[i]=HLTBitStrings[i]+"_v1_Prescl";
  jetpp[3]->SetBranchAddress( HLTPresclBranches[0].c_str() , &jet40_p_I);
  jetpp[3]->SetBranchAddress( HLTPresclBranches[1].c_str() , &jet60_p_I);
  jetpp[3]->SetBranchAddress( HLTPresclBranches[2].c_str() , &jet80_p_I);
  jetpp[3]->SetBranchAddress( HLTPresclBranches[3].c_str() , &jet100_p_I);

  std::string L1PresclBranches[N_L1Bits]  ; 
  for(int i=0;i<N_HLTBits;i++) L1PresclBranches[i]=L1BitStrings[i]+"_Prescl";
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


  // Declare the output file, declare hists after for easy write
  std::cout<<"opening output file "<<outfile<<std::endl;
  TFile *fout = new TFile(outfile.c_str(),"RECREATE");
  fout->cd();


  // declare hists
  //evt count
  TH1F *h_NEvents         = new TH1F("NEvents"         , "NEvents", 1,0,2);
  TH1F *h_NEvents_read    = new TH1F("NEvents_read"    , "NEvents read", 1,0,2);
  TH1F *h_NEvents_skimCut = new TH1F("NEvents_skimCut" , "NEvents read post skimCut", 1,0,2);
  TH1F *h_NEvents_vzCut   = new TH1F("NEvents_vzCut"   , "NEvents read post vzCut AND skimCut", 1,0,2);

  //TH1F *h_WNEvents         = new TH1F("WNEvents"        , "wghtdNEvents", 1,0,2);
  //TH1F *h_WNEvents         = new TH1F("WNEvents_read"   , "wghtdNEvents read", 1,0,2);
  //TH1F *h_WNEvents_skimCut = new TH1F("WNEvents_skimCut", "wghtdNEvents read post skimCut", 1,0,2);
  //TH1F *h_WNEvents_vzCut   = new TH1F("WNEvents_vzCut"  , "wghtdNEvents read post vzCut AND skimCut", 1,0,2);
  
  //evt counts post all evt cuts. thrown out because no jets left in the event pre or post cuts
  TH1F *h_NEvents_withJets           = new TH1F("NEvents_withJets" , 
		  				"NEvents read post evt cuts, w/ jets", 1,0,2);
  TH1F *h_NEvents_withJets_kmatCut   = new TH1F("NEvents_withJets_kmatCut" , 
						"NEvents read post evt cuts, w/ jets post kmatCut", 1,0,2);
  TH1F *h_NEvents_withJets_JetIDCut = new TH1F("NEvents_withJets_JetIDCut" , 
						"NEvents read post evt cuts, w/ jets post kmatCut AND JetID Cut", 1,0,2);

  //vz
  TH1F *hVz = new TH1F("hVz","",  60,-15.,15.); 
  TH1F *hWVz = new TH1F("hWeightedVz","",  60,-15.,15.);// W for weighted

  //jet counts
  TH1F *h_NJets          = new TH1F("NJets","NJets read", 1,0,2);
  TH1F *h_NJets_kmatCut  = new TH1F("NJets_kmatCut ","NJets read post kmatCut ", 1,0,2);
  TH1F *h_NJets_JetIDCut = new TH1F("NJets_JetIDCut","NJets read post JetIDCut AND kmatCut", 1,0,2);
  
  //TH1F *h_WNJets          = new TH1F("WNJets"         ,"wghtdNJets of selected events, no jet cuts", 1,0,2);
  //TH1F *h_WNJets_kmatCut  = new TH1F("WNJets_kmatCut ","wghtdNJets post kmatCut ", 1,0,2);// kmat for kinematics
  //TH1F *h_WNJets_JetIDCut = new TH1F("WNJets_JetIDCut","wghtdNJets post JetIDCut", 1,0,2);

  //jet QA
  TH1F *hJetQA[2][N_vars];
  for(int k = 0; k<2; ++k){ for(int j = 0; j<N_vars; ++j){

      //jtpt and rawpt and leadJetPt special binning
      if(j<=1||j==22) hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) ,
				        Form(";%s;", var[j].c_str()) , 500,0,500);
      //jteta special binning 
      else if(j==2) hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) ,
					     Form(";%s;",var[j].c_str()) , 100,-5,+5);
      //jtphi and dphi special binning 
      else if(j==3) hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) ,
					     Form(";%s;",var[j].c_str()) , 100,-4,+4);

      else if(j==21) hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) ,
					     Form(";%s;",var[j].c_str()) , 50,0,+4);
      //same binning for all others
      else /*if(j>=4)*/ hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) ,
						 Form(";%s;",var[j].c_str()) , 200,0,2);    }  }

  // trigger (also jet) level
  TH1F  *hpp_TrgObj40[2], *hpp_TrgObj60[2], *hpp_TrgObj80[2], *hpp_TrgObj100[2];// separate trigs
  TH1F *hpp_TrgObjComb[2], *hpp_CombJetpT[2];// combos of trigs

  hpp_TrgObj40[0]   = new TH1F(Form("hpp_HLT40_R%d_%s" , radius,etaWidth), 
			       Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s" , radius,etaWidth), 2000,0,2000);
  hpp_TrgObj60[0]   = new TH1F(Form("hpp_HLT60_R%d_%s" , radius,etaWidth), 
			       Form("Spectra from  Jet 60 && !jet80 R%d %s" , radius,etaWidth), 2000,0,2000);
  hpp_TrgObj80[0]   = new TH1F(Form("hpp_HLT80_R%d_%s" , radius,etaWidth), 
			       Form("Spectra from  Jet 80 R%d %s" , radius,etaWidth), 2000,0,2000);
  hpp_TrgObj100[0]  = new TH1F(Form("hpp_HLT100_R%d_%s" , radius,etaWidth), 
			       Form("Spectra from  Jet 100 R%d %s" , radius,etaWidth), 2000,0,2000);
  hpp_TrgObjComb[0] = new TH1F(Form("hpp_HLTComb_R%d_%s" , radius,etaWidth), 
			       Form("Trig Combined Spectra R%d %s" , radius,etaWidth), 2000,0,2000);
  hpp_CombJetpT[0]  = new TH1F(Form("hpp_TrgCombTest_R%d_%s" , radius,etaWidth), 
			       Form("Trig Comb Spectra kMethod R%d %s" , radius,etaWidth), 2000,0,2000);

  hpp_TrgObj40[1]   = new TH1F(Form("hpp_HLT40_wJetID_R%d_%s" , radius,etaWidth), 
			       Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s" , radius,etaWidth), 2000,0,2000);
  hpp_TrgObj60[1]   = new TH1F(Form("hpp_HLT60_wJetID_R%d_%s" , radius,etaWidth), 
			       Form("Spectra from  Jet 60 && !jet80 R%d %s" , radius,etaWidth), 2000,0,2000);
  hpp_TrgObj80[1]   = new TH1F(Form("hpp_HLT80_wJetID_R%d_%s" , radius,etaWidth), 
			       Form("Spectra from  Jet 80 R%d %s" , radius,etaWidth), 2000,0,2000);
  hpp_TrgObj100[1]  = new TH1F(Form("hpp_HLT100_wJetID_R%d_%s" , radius,etaWidth), 
			       Form("Spectra from  Jet 100 R%d %s" , radius,etaWidth), 2000,0,2000);
  hpp_TrgObjComb[1] = new TH1F(Form("hpp_HLTComb_wJetID_R%d_%s" , radius,etaWidth), 
			       Form("Trig Combined Spectra R%d %s" , radius,etaWidth), 2000,0,2000);
  hpp_CombJetpT[1]  = new TH1F(Form("hpp_TrgCombTest_wJetID_R%d_%s" , radius,etaWidth), 
			       Form("Trig Comb Spectra kMethod w/JetID R%d %s" , radius,etaWidth), 2000,0,2000);



  // EVENT LOOP
  UInt_t NEvents_jetAnalyzr=jetpp[0]->GetEntries();   // preskim event count from files
  UInt_t NEvents_skimAnalyzr=jetpp[3]->GetEntries();   // preskim event count from files
  if(!(NEvents_jetAnalyzr==NEvents_skimAnalyzr)){
    std::cout<<"jetAnalyzr # evts  = " << NEvents_jetAnalyzr  << std::endl;
    std::cout<<"skimAnalyzr # evts = " << NEvents_skimAnalyzr << std::endl;
    assert(NEvents_jetAnalyzr==NEvents_skimAnalyzr);  }

  UInt_t NEvents_allFiles=NEvents_jetAnalyzr;   // preskim event count from files
  for(UInt_t i=0;i < NEvents_allFiles; ++i) h_NEvents->Fill(1);
  
  UInt_t NEvents_read=0;
  if(fastDebugMode) NEvents_read = 1000*(endfile-startfile+1); 
  else NEvents_read = NEvents_allFiles;
  
  UInt_t NEvents_40=0, NEvents_60=0, NEvents_80=0, NEvents_100=0;  //trigger event counts

  std::cout<<"reading "<<NEvents_read<<" events"<<std::endl;  
  for(UInt_t nEvt = 0; nEvt < NEvents_read; ++nEvt) {//event loop   

    if(fastDebugMode&&nEvt%1000==0)std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;
    else if (nEvt%10000==0)std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;

    //for(int i=0;i<N_trees; ++i)jetpp[i]->GetEntry(nEvt);
    jetpp[0]->GetEntry(nEvt);
    h_NEvents_read->Fill(1);


    // skim/HiEvtAnalysis criteria
    if( pHBHENoiseFilter_I==0     || 
        pBeamScrapingFilter_I==0  || 
        pprimaryvertexFilter_I==0  ) continue;    
    h_NEvents_skimCut->Fill(1);

    if( fabs(vz_F)>15.              ) continue;
    h_NEvents_vzCut->Fill(1);
    

    // grab largest triggerPt among HLT paths' trigger object collectiosn
    double triggerPt=0.;
    unsigned int trgObj40_size =trgObjpt_40->size(), trgObj60_size =trgObjpt_60->size();
    unsigned int trgObj80_size =trgObjpt_80->size(), trgObj100_size=trgObjpt_100->size();
    
    int treePrescl[4]={ (jet40_p_I*jet40_l1s_ps_I), (jet60_p_I*jet60_l1s_ps_I), //total prescale=L1_ps*HLT_ps
			(jet80_p_I*jet80_l1s_ps_I), (jet100_p_I*jet100_l1s_ps_I) };    

    bool trgDec[4]={ (bool)jet40_I, (bool)jet60_I, //HLT decision=HLT&&L1=HLT&&true
		     (bool)jet80_I, (bool)jet100_I  }; 
    
    if(jet40_I){ for(unsigned int itt=0; itt<trgObj40_size; ++itt){
	if(trgObjpt_40->at(itt) > triggerPt) 
	  triggerPt = trgObjpt_40->at(itt);       } }
    if(jet60_I){ for(unsigned int itt=0; itt<trgObj60_size; ++itt){
	if(trgObjpt_60->at(itt) > triggerPt) 
	  triggerPt = trgObjpt_60->at(itt);       } }
    if(jet80_I){ for(unsigned int itt=0; itt<trgObj80_size; ++itt){
    	if(trgObjpt_80->at(itt) > triggerPt) 
	  triggerPt = trgObjpt_80->at(itt);       } }
    if(jet100_I){ for(unsigned int itt=0; itt<trgObj100_size; ++itt){
    	if(trgObjpt_100->at(itt) > triggerPt) 
	  triggerPt = trgObjpt_100->at(itt);      } }    

    
    // check trigger decisions for events + exclusivity between them, count events
    bool is40  = false, is60  = false, is80  = false, is100 = false;
    if(trgDec[0] && triggerPt>=40. && triggerPt<60. ) is40  = true;    
    if(trgDec[1] && triggerPt>=60. && triggerPt<80. ) is60  = true;    
    if(trgDec[2] && triggerPt>=80. && triggerPt<100.) is80  = true;    
    if(trgDec[3] && triggerPt>=100.                ) is100 = true;    

    if(is40)  NEvents_40++; 
    if(is60)  NEvents_60++; 
    if(is80)  NEvents_80++; 
    if(is100) NEvents_100++;
    //if(true)  NEvents_any++;

    
    //if(debugMode&&nEvt%1000==0)std::cout<<"cpu'ing trigger weight for evt# "<<nEvt<<std::endl;
    double weight_eS = trigComb(trgDec, treePrescl, triggerPt); //kurt method
    

    // fill evt vz histo
    hVz->Fill(vz_F, 1.0);    
    hWVz->Fill(vz_F, weight_eS);


    // JET LOOP 
    // for Aj/xj computation
    bool firstGoodJetFound=false, secondGoodJetFound=false; 
    float firstGoodJetPt=-1., secondGoodJetPt=-1.;
    float firstGoodJetPhi=-1., secondGoodJetPhi=-1.;
    bool firstGoodJetFound_wJetID=false, secondGoodJetFound_wJetID=false; 
    float firstGoodJetPt_wJetID=-1., secondGoodJetPt_wJetID=-1.;
    float firstGoodJetPhi_wJetID=-1., secondGoodJetPhi_wJetID=-1.;
    
    // for event counting
    bool hNEvts_withJets_Filled=false;
    bool hNEvts_withJets_kmatCut_Filled=false, hNEvts_withJets_JetIDCut_Filled=false;

    for(int jet = 0; jet<nref_I; ++jet){

      // event+jet counting
      h_NJets->Fill(1);
      if(!hNEvts_withJets_Filled){
	h_NEvents_withJets->Fill(1);
	hNEvts_withJets_Filled=true;      }      


      float recpt = pt_F[jet];
      float rawpt = rawpt_F[jet];
      float phi = phi_F[jet];

      // kmatCut
      if( fabs(eta_F[jet]) > 2. ||
          recpt <= 15.           ) continue;     


      // trig plots
      if(true)  hpp_CombJetpT[0]->Fill(recpt, weight_eS); //kurts method
      if(is40)  hpp_TrgObj40[0]->Fill(recpt, treePrescl[0]);
      if(is60)  hpp_TrgObj60[0]->Fill(recpt, treePrescl[1]);
      if(is80)  hpp_TrgObj80[0]->Fill(recpt, treePrescl[2]);
      if(is100) hpp_TrgObj100[0]->Fill(recpt, treePrescl[3]);


      // jetQA noJetID
      hJetQA[0][0]->Fill(recpt, weight_eS);
      hJetQA[0][1]->Fill(rawpt_F[jet], weight_eS);

      if(recpt>50.){
	hJetQA[0][2]->Fill(eta_F[jet], weight_eS);
	hJetQA[0][3]->Fill(phi_F[jet], weight_eS);
	hJetQA[0][4]->Fill(trkMax_F[jet]/rawpt, weight_eS);
	hJetQA[0][5]->Fill(trkSum_F[jet]/rawpt, weight_eS);
	hJetQA[0][6]->Fill(trkHardSum_F[jet]/rawpt, weight_eS);
	hJetQA[0][7]->Fill(chMax_F[jet]/rawpt, weight_eS);
	hJetQA[0][8]->Fill(chSum_F[jet]/rawpt, weight_eS);
	hJetQA[0][9]->Fill(chHardSum_F[jet]/rawpt, weight_eS);
	hJetQA[0][10]->Fill(phMax_F[jet]/rawpt, weight_eS);
	hJetQA[0][11]->Fill(phSum_F[jet]/rawpt, weight_eS);
	hJetQA[0][12]->Fill(phHardSum_F[jet]/rawpt, weight_eS);
	hJetQA[0][13]->Fill(neMax_F[jet]/rawpt, weight_eS);
	hJetQA[0][14]->Fill(neSum_F[jet]/rawpt, weight_eS);
	hJetQA[0][15]->Fill(eMax_F[jet]/rawpt, weight_eS);
	hJetQA[0][16]->Fill(eSum_F[jet]/rawpt, weight_eS);
	hJetQA[0][17]->Fill(muMax_F[jet]/rawpt, weight_eS);
	hJetQA[0][18]->Fill(muSum_F[jet]/rawpt, weight_eS);      }


      // want to compute Aj/xj with the first two good jets that meet all criteria found
      if (!firstGoodJetFound && !secondGoodJetFound){
        if(recpt>60.) { firstGoodJetFound=true;
          firstGoodJetPt=recpt; 
	  firstGoodJetPhi=phi; }
        else { firstGoodJetFound=false; }           }

      else if ( firstGoodJetFound && !secondGoodJetFound ) {
        secondGoodJetFound=true;
        if (recpt>firstGoodJetPt) {std::cerr<<std::endl<<
	    "WARNING: picked wrong jet for lead jet! Swapping..."<<std::endl<<std::endl;
          secondGoodJetPt=firstGoodJetPt;          secondGoodJetPhi=firstGoodJetPhi;
          firstGoodJetPt=recpt;          firstGoodJetPhi=phi; }
        else { secondGoodJetPt=recpt;
          secondGoodJetPhi=phi; }                          }

      if( firstGoodJetFound && secondGoodJetFound ) {
	float A_j=(firstGoodJetPt-secondGoodJetPt)/(firstGoodJetPt+secondGoodJetPt);
	float x_j=(secondGoodJetPt/firstGoodJetPt); 	
	float dphi=deltaphi(firstGoodJetPhi,secondGoodJetPhi);
	hJetQA[0][19]->Fill( A_j , weight_eS ); 
	hJetQA[0][20]->Fill( x_j , weight_eS );       
	hJetQA[0][21]->Fill( dphi , weight_eS );       
	hJetQA[0][22]->Fill( firstGoodJetPt , weight_eS );       }
      

      // jet/event counts
      h_NJets_kmatCut->Fill(1);
      if(!hNEvts_withJets_kmatCut_Filled){
	h_NEvents_withJets_kmatCut->Fill(1);
	hNEvts_withJets_kmatCut_Filled=true;      }

      
      // apply JetID
      if (chSum_F[jet]/rawpt <= 0    ||
          neSum_F[jet]/rawpt >= 0.99 || 
          phSum_F[jet]/rawpt >= 0.99 ||
          chN_I[jet] <= 0             ) continue;
      

      // trig plots
      if(true)  hpp_CombJetpT[1]->Fill(recpt, weight_eS); //kurts method
      if(is40 )  hpp_TrgObj40[1]->Fill(recpt, treePrescl[0]);
      if(is60 )  hpp_TrgObj60[1]->Fill(recpt, treePrescl[1]);
      if(is80 )  hpp_TrgObj80[1]->Fill(recpt, treePrescl[2]);
      if(is100) hpp_TrgObj100[1]->Fill(recpt, treePrescl[3]);      

      // jetQA wJetID
      hJetQA[1][0]->Fill(recpt, weight_eS);
      hJetQA[1][1]->Fill(rawpt_F[jet], weight_eS);

      if(recpt>50.){     
	hJetQA[1][2]->Fill(eta_F[jet], weight_eS);
	hJetQA[1][3]->Fill(phi_F[jet], weight_eS);
	hJetQA[1][4]->Fill(trkMax_F[jet]/rawpt, weight_eS);
	hJetQA[1][5]->Fill(trkSum_F[jet]/rawpt, weight_eS);
	hJetQA[1][6]->Fill(trkHardSum_F[jet]/rawpt, weight_eS);
	hJetQA[1][7]->Fill(chMax_F[jet]/rawpt, weight_eS);
	hJetQA[1][8]->Fill(chSum_F[jet]/rawpt, weight_eS);
	hJetQA[1][9]->Fill(chHardSum_F[jet]/rawpt, weight_eS);
	hJetQA[1][10]->Fill(phMax_F[jet]/rawpt, weight_eS);
	hJetQA[1][11]->Fill(phSum_F[jet]/rawpt, weight_eS);
	hJetQA[1][12]->Fill(phHardSum_F[jet]/rawpt, weight_eS);
	hJetQA[1][13]->Fill(neMax_F[jet]/rawpt, weight_eS);
	hJetQA[1][14]->Fill(neSum_F[jet]/rawpt, weight_eS);
	hJetQA[1][15]->Fill(eMax_F[jet]/rawpt, weight_eS);
	hJetQA[1][16]->Fill(eSum_F[jet]/rawpt, weight_eS);
	hJetQA[1][17]->Fill(muMax_F[jet]/rawpt, weight_eS);
	hJetQA[1][18]->Fill(muSum_F[jet]/rawpt, weight_eS);      }

      // want to compute Aj/xj with the first two good jets that meet all criteria found
      if (!firstGoodJetFound_wJetID && !secondGoodJetFound_wJetID){
        if(recpt>60.) {
          firstGoodJetFound_wJetID=true;
          firstGoodJetPt_wJetID=recpt;
          firstGoodJetPhi_wJetID=phi; }
        else {
          firstGoodJetFound_wJetID=false; }  }

      else if (firstGoodJetFound_wJetID && !secondGoodJetFound_wJetID){
        secondGoodJetFound_wJetID=true;
        if(recpt>firstGoodJetPt_wJetID) {std::cout<<std::endl<<"WARNING: picked wrong jet for lead jet! Swapping..."<<std::endl<<std::endl;
          secondGoodJetPt_wJetID=firstGoodJetPt_wJetID;
          secondGoodJetPhi_wJetID=firstGoodJetPhi_wJetID;
          firstGoodJetPt_wJetID=recpt;
          firstGoodJetPhi_wJetID=phi; }
        else {
          secondGoodJetPt_wJetID=recpt;
          secondGoodJetPhi_wJetID=phi; }  }


      if( firstGoodJetFound_wJetID && secondGoodJetFound_wJetID ) {
	float A_j=(firstGoodJetPt_wJetID-secondGoodJetPt_wJetID)/
	  (firstGoodJetPt_wJetID+secondGoodJetPt_wJetID);
	float x_j=(secondGoodJetPt_wJetID/firstGoodJetPt_wJetID); 	
	float dphi=deltaphi(firstGoodJetPhi_wJetID,secondGoodJetPhi_wJetID);
	hJetQA[1][19]->Fill( A_j , weight_eS ); 
	hJetQA[1][20]->Fill( x_j , weight_eS );       
      	hJetQA[1][21]->Fill( dphi , weight_eS );       
      	hJetQA[1][22]->Fill( firstGoodJetPt_wJetID , weight_eS );       }
      

      // jet/event counts
      h_NJets_JetIDCut->Fill(1);
      if(!hNEvts_withJets_JetIDCut_Filled){
	h_NEvents_withJets_JetIDCut->Fill(1);
	hNEvts_withJets_JetIDCut_Filled=true;      }

    }//end jet loop
  }//end event loop
  std::cout<<std::endl;


  //raghavs method
  hpp_TrgObj40[0]->Print("base");
  hpp_TrgObj60[0]->Print("base");
  hpp_TrgObj80[0]->Print("base");
  hpp_TrgObj100[0]->Print("base");

  hpp_TrgObjComb[0]->Add(hpp_TrgObj100[0]);
  hpp_TrgObjComb[0]->Add(hpp_TrgObj80[0] );
  hpp_TrgObjComb[0]->Add(hpp_TrgObj60[0] );
  hpp_TrgObjComb[0]->Add(hpp_TrgObj40[0] );
  hpp_TrgObjComb[0]->Print("base");
  
  hpp_TrgObj40[1]->Print("base");
  hpp_TrgObj60[1]->Print("base");
  hpp_TrgObj80[1]->Print("base");
  hpp_TrgObj100[1]->Print("base");

  hpp_TrgObjComb[1]->Add(hpp_TrgObj100[1]);
  hpp_TrgObjComb[1]->Add(hpp_TrgObj80[1] );
  hpp_TrgObjComb[1]->Add(hpp_TrgObj60[1] );
  hpp_TrgObjComb[1]->Add(hpp_TrgObj40[1] );
  hpp_TrgObjComb[1]->Print("base");


  //kurts method??
  hpp_CombJetpT[0]->Print("base");
  hpp_CombJetpT[1]->Print("base");
  

  std::cout<<std::endl<<"/// Job Event-Loop Summary ///"<<std::endl<<std::endl;
  std::cout<<"Total Num of Events in file(s) opened       = " <<h_NEvents->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read from those file(s) = " <<h_NEvents_read->GetEntries()<<std::endl;

  if(!debugMode){
    std::cout<<"Total Num of Events read passing skimCuts                              = " <<
      h_NEvents_skimCut->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events read passing vzCuts and skimCuts (aka good events) = " <<
      h_NEvents_vzCut->GetEntries()<<std::endl;
    
    std::cout<<std::endl;
    std::cout<<"Total Num of good Events w/ jets                             = " <<
      h_NEvents_withJets->GetEntries()<<std::endl;	   
    std::cout<<"Total Num of good Events, w/ jets, post kmatCuts             = " <<
      h_NEvents_withJets_kmatCut->GetEntries()<<std::endl; 
    std::cout<<"Total Num of good Events, w/ jets, post kmatCut AND JetIDCut = " <<
      h_NEvents_withJets_JetIDCut->GetEntries()<<std::endl;
    }

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
  std::cout<<"Total Num of Jets read from good events post kinCut AND JetIDCut = " <<
    h_NJets_JetIDCut->GetEntries()<<std::endl;

  std::cout<<std::endl;


  std::cout<<"writing output file..."<<std::endl;
  fout->Write(); //this writes duplicates
  //fout->Write("",TObject::kOverwrite);//this does not


  if(debugMode)std::cout<<"misc clean up.."<<std::endl;
  trgObjpt_40->clear();
  trgObjpt_60->clear();
  trgObjpt_80->clear();
  trgObjpt_100->clear();

  std::cout<<std::endl<<"readForests_ppData finished."<<std::endl;  timer.Stop();
  std::cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<std::endl;
  std::cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<std::endl;

  return 0;
}// end readForests_ppData


////// main //////
// acts as the frontend control for .exe file
int main(int argc, char *argv[]){
  
  // error, not enough arguments
  int rStatus = -1;
  if(argc!=8 && argc!=1){
    std::cout<<"for tests on default inputs, do..." <<std::endl;
    std::cout<<"./readForests_ppData.exe";
    std::cout<<std::endl<<std::endl;
    std::cout<<"for actually running, do..."<<std::endl;
    std::cout<<"./readForests_ppData.exe ";
    std::cout<<"<startFile> <endFile> <inputFileList> <outputFilename> ";
    std::cout<<"<jetRadius> <jetType> <debugMode>";
    std::cout<<std::endl<<std::endl;
    std::cout<<"rStatus="<<rStatus<<std::endl;
    return rStatus;
  }
  
  // good input, run
  rStatus=1;
  if(argc==1) rStatus = readForests_ppData();
  else{//read input argument vector
    int startfile= atoi(argv[1]); int endfile= atoi(argv[2]); std::string inputFileList=argv[3]; std::string outputFileName=argv[4];
    int jetRadius= atoi(argv[5]); std::string jetType=argv[6];     bool debug=(bool)atoi(argv[7]);      
    rStatus = readForests_ppData( startfile, endfile, inputFileList, outputFileName,
				jetRadius, jetType, debug);
  }
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}



    //if(debugMode&&nEvt%250==0)std::cout <<"triggerPt ="<<triggerPt<<std::endl;
    //if(debugMode&&nEvt%250==0)std::cout<<"jet 40, trgDec[0]  = "<< trgDec[0]  <<std::endl;
    //if(debugMode&&nEvt%250==0)std::cout<<"jet 60, trgDec[1]  = "<< trgDec[1]  <<std::endl;
    //if(debugMode&&nEvt%250==0)std::cout<<"jet 80, trgDec[2]  = "<< trgDec[2]  <<std::endl;
    //if(debugMode&&nEvt%250==0)std::cout<<"jet 100, trgDec[3] = "<< trgDec[3]  <<std::endl<<std::endl;

    //if(debugMode&&nEvt%250==0)std::cout<<"after checking trigger pt "<<std::endl;
    //if(debugMode&&nEvt%250==0)std::cout<<"jet 40,  is40  = "<< is40  <<std::endl;
    //if(debugMode&&nEvt%250==0)std::cout<<"jet 60,  is60  = "<< is60  <<std::endl;
    //if(debugMode&&nEvt%250==0)std::cout<<"jet 80,  is80  = "<< is80  <<std::endl;
    //if(debugMode&&nEvt%250==0)std::cout<<"jet 100, is100 = "<< is100  <<std::endl<<std::endl;

