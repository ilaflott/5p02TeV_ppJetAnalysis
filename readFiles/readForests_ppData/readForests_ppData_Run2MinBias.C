// custom header
#include "readFiles_ppData_Run2MinBias.h"

//// readFiles_ppData
// ---------------------------------------------------------------------------------------------------------------
int readFiles_ppData_Run2MinBias(int startfile , int endfile , std::string inFilelist , std::string outfile ,
		     int radius , std::string jetType , bool debugMode ){ 

  // for monitoring performance
  TStopwatch timer;  timer.Start();

  // debug mode settings+warnings
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
  
  // form jetTreeName, initialize tree name array
  std::string trees[N_trees];
  std::string jetTreeName="ak"+std::to_string(radius)+jetType+"JetAnalyzer/t";
  trees[0]=jetTreeName;
  if(debugMode)std::cout<<"looking at jetTree "<<trees[0]<<std::endl;
  for(int i=1;i<N_trees;++i)trees[i]=treeNames[i];
  
  // declare TChains for each tree we're interested in
  TChain* jetpp[N_trees];
  for(int t = 0;t<N_trees;++t)  jetpp[t] = new TChain( trees[t].data() );
  
  // open filelist
  std::cout<<"filelist used is "<<inFilelist<<std::endl;
  std::ifstream instr_Forest(inFilelist.c_str(),std::ifstream::in);
  std::string filename_Forest;  

  // add files, including startfile+endfile, to chains
  std::string lastFileAdded=""; bool treesAdded=false;
  for(int ifile = 0; ifile<=endfile; ++ifile){//file loop
    // loop to starting file, check for end of filelist
    instr_Forest>>filename_Forest; 
    if(ifile<startfile){lastFileAdded=filename_Forest;continue;}
    if(filename_Forest==lastFileAdded){std::cout<<"end of filelist!"<<std::endl;break;}//end of filelist condition
    std::cout<<"adding file #"<<ifile;  if(debugMode)std::cout<<" "<<filename_Forest;  std::cout<<" to TChain(s)"<<std::endl;
    for(int t = 0;t<N_trees;++t){ jetpp[t]->Add(filename_Forest.c_str()); treesAdded=true; }
    lastFileAdded=filename_Forest;
  }//end file loop
  assert(treesAdded);//avoid segfault later

  // friend all trees to jet analyzer tree in jetpp[0]
  for(int t = 1;t<N_trees;++t)jetpp[0]->AddFriend( jetpp[t] );
  
  // declare variables/arrays for input files
  // ak${radius}${jetType}JetAnalyzer
  //jets
  int nref_F;  
  float pt_F[1000];     float eta_F[1000];  float phi_F[1000]; 
  float rawpt_F[1000];  float jtpu_F[1000];
  //charged particles
  int chN_F[1000];  
  float chSum_F[1000]; float chHardSum_F[1000];  float chMax_F[1000];  //float chHardMax_F[1000];   
  //tracks
  float trkSum_F[1000];  float trkHardSum_F[1000];  float trkMax_F[1000];  //float trkHardMax_F[1000];
  //photons
  float phSum_F[1000];  float phHardSum_F[1000];  float phMax_F[1000];  //float phHardMax_F[1000];   
  //e, mu, neutral particles
  float eSum_F[1000];  float eMax_F[1000];  
  float muSum_F[1000];  float muMax_F[1000];  
  float neSum_F[1000];  float neMax_F[1000];   

  // HiEvtAnalyzer
  ULong64_t evt_F;   UInt_t run_F;   UInt_t lumi_F;   float vz_F;

  // skimanalysis
  int pBeamScrapingFilter_F;   int pHBHENoiseFilter_F;   int pprimaryvertexFilter_F;

  // hltanalysis
  //L1, one prescale branch per L1 path we're curious about
  int jet40_l1seed_p_F;  int jet60_l1seed_p_F;  int jet80_l1seed_p_F;  int jet100_l1seed_p_F;  //prescale
  //HLT, two branches, one decision one prescale, per HLT path we're curious about
  int jet40_F;           int jet60_F;           int jet80_F;           int jet100_F;           //decision
  int jet40_p_F;         int jet60_p_F;         int jet80_p_F;         int jet100_p_F;         //prescale

  // four trees, four specific HLT trigger pt objects
  //currently (5/26/16): HLT_AK${radius}CaloJet${trgObjpt}_Eta${}_v
  std::vector<double> *trgObjpt_40  ;   
  std::vector<double> *trgObjpt_60  ;   
  std::vector<double> *trgObjpt_80  ;   
  std::vector<double> *trgObjpt_100 ;

  // set branch addresses for the input file treesp, map them to above variables
  // ak${radius}${jetType}JetAnalyzer
  //jets
  jetpp[0]->SetBranchAddress("nref",&nref_F);
  jetpp[0]->SetBranchAddress("jtpt",&pt_F);
  jetpp[0]->SetBranchAddress("jteta",&eta_F);
  jetpp[0]->SetBranchAddress("jtphi",&phi_F);
  jetpp[0]->SetBranchAddress("rawpt",&rawpt_F);
  jetpp[0]->SetBranchAddress("jtpu",&jtpu_F);
  //charged particles
  jetpp[0]->SetBranchAddress("chargedN",&chN_F);
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
  //neutral partcles assoc. w/ jet
  jetpp[0]->SetBranchAddress("neutralSum",&neSum_F);
  jetpp[0]->SetBranchAddress("neutralMax",&neMax_F);

  // hiEvtAnalyzer
  jetpp[1]->SetBranchAddress("EvtInfo.EvtNo",&evt_F);
  jetpp[1]->SetBranchAddress("EvtInfo.RunNo",&run_F);
  jetpp[1]->SetBranchAddress("EvtInfo.LumiNo",&lumi_F);
  jetpp[1]->SetBranchAddress("EvtInfo.PVz",&vz_F);

  // skimanalysis
  jetpp[2]->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_F);
  jetpp[2]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilter_F);
  jetpp[2]->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryvertexFilter_F);

  // hltanalysis, I should try figuring out a better way to do this part
  std::string HLTBranches[N_HLTBits]; 
  for(int i=0;i<N_L1Bits;i++) HLTBranches[i]=HLTBitStrings[i]+"_v1";
  jetpp[3]->SetBranchAddress( HLTBranches[0].c_str() , &jet40_F);
  jetpp[3]->SetBranchAddress( HLTBranches[1].c_str() , &jet60_F);
  jetpp[3]->SetBranchAddress( HLTBranches[2].c_str() , &jet80_F);
  jetpp[3]->SetBranchAddress( HLTBranches[3].c_str() , &jet100_F);
  
  std::string HLTPresclBranches[N_HLTBits]; 
  for(int i=0;i<N_HLTBits;i++) HLTPresclBranches[i]=HLTBitStrings[i]+"_v1_Prescl";
  jetpp[3]->SetBranchAddress( HLTPresclBranches[0].c_str() , &jet40_p_F);
  jetpp[3]->SetBranchAddress( HLTPresclBranches[1].c_str() , &jet60_p_F);
  jetpp[3]->SetBranchAddress( HLTPresclBranches[2].c_str() , &jet80_p_F);
  jetpp[3]->SetBranchAddress( HLTPresclBranches[3].c_str() , &jet100_p_F);

  std::string L1PresclBranches[N_L1Bits]  ; 
  for(int i=0;i<N_HLTBits;i++) L1PresclBranches[i]=L1BitStrings[i]+"_Prescl";
  jetpp[3]->SetBranchAddress( L1PresclBranches[0].c_str()  , &jet40_l1seed_p_F);
  jetpp[3]->SetBranchAddress( L1PresclBranches[1].c_str()  , &jet60_l1seed_p_F);
  jetpp[3]->SetBranchAddress( L1PresclBranches[2].c_str()  , &jet80_l1seed_p_F);
  jetpp[3]->SetBranchAddress( L1PresclBranches[3].c_str()  , &jet100_l1seed_p_F);

  // specific HLT path trees
  //ONE path ONE tree ONE pt branch (see trees[] in header)
  //e.g. trgObjpt_40 is filled with jet pt from the specific jet40 HLT tree/branch 
  jetpp[4]->SetBranchAddress("pt",&trgObjpt_40);
  jetpp[5]->SetBranchAddress("pt",&trgObjpt_60);  
  jetpp[6]->SetBranchAddress("pt",&trgObjpt_80);  
  jetpp[7]->SetBranchAddress("pt",&trgObjpt_100);

  // Declare the output File and the necessary histograms after that:
  std::cout<<"opening output file "<<outfile<<std::endl;
  TFile *fout = new TFile(outfile.c_str(),"RECREATE");
  fout->cd();

  // book jet variable histograms
  TH1F *hJetQA[2][N_vars];
  for(int k = 0; k<2; ++k){
    for(int j = 0; j<N_vars; ++j){
      //jtpt and rawpt special binning for QA
      if(j<=1) hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),
				       Form(";%s;",var[j].c_str()),500,  0, 500 );
      //jteta special binning for QA
      else if(j==2) hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),
					    Form(";%s;",var[j].c_str()),100, -5, +5  );
      //jtphi special binning for QA
      else if(j==3) hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),
					    Form(";%s;",var[j].c_str()),100, -4, +4  );
      //same binning for all other QA, alter later if needed
      else if(j>=4) hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),
					    Form(";%s;",var[j].c_str()),200,  0, 2   );
    }
  }

  // for vertex position, can use for vz weighting
  TH1F *hVz = new TH1F("hVz","",200, -20, 20);
  TH1F *hWeightedVz = new TH1F("hWeightedVz","",200, -20, 20);

  // for triggerpt spectra before/after jetID
  //consider renaming h_ppTrgObj, 'hpp' confusing
  TH1F *hpp_TrgObj40[2];
  TH1F *hpp_TrgObj60[2];
  TH1F *hpp_TrgObj80[2];
  TH1F *hpp_TrgObj100[2];
  TH1F *hpp_TrgObjComb[2];
  TH1F *hpp_CombJetpT[2];
  
  //40, 60, 80, 100
  hpp_TrgObj40[0]   = new TH1F(Form("hpp_HLT40_noJetID_R%d_%s"                        , radius, etaWidth), 
			       Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_TrgObj60[0]   = new TH1F(Form("hpp_HLT60_noJetID_R%d_%s"               , radius, etaWidth), 
			       Form("Spectra from  Jet 60 && !jet80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_TrgObj80[0]   = new TH1F(Form("hpp_HLT80_noJetID_R%d_%s"     , radius, etaWidth), 
			       Form("Spectra from  Jet 80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_TrgObj100[0]  = new TH1F(Form("hpp_HLT100_noJetID_R%d_%s"     , radius, etaWidth), 
			       Form("Spectra from  Jet 100 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  //combined
  hpp_TrgObjComb[0] = new TH1F(Form("hpp_HLTComb_noJetID_R%d_%s"    , radius, etaWidth), 
			       Form("Trig Combined Spectra R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_CombJetpT[0]  = new TH1F(Form("hpp_TrgCombTest_noJetID_R%d_%s"           , radius, etaWidth), 
			       Form("Trig Combined Spectra KurtMethod R%d %s " , radius, etaWidth), 2000, 0, 2000);
  //40, 60, 80, 100
  hpp_TrgObj40[1]   = new TH1F(Form("hpp_HLT40_JetID_R%d_%s"                          , radius, etaWidth), 
			       Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_TrgObj60[1]   = new TH1F(Form("hpp_HLT60_JetID_R%d_%s"                 , radius, etaWidth), 
			       Form("Spectra from  Jet 60 && !jet80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_TrgObj80[1]   = new TH1F(Form("hpp_HLT80_JetID_R%d_%s"       , radius, etaWidth), 
			       Form("Spectra from  Jet 80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_TrgObj100[1]  = new TH1F(Form("hpp_HLT100_JetID_R%d_%s"       , radius, etaWidth), 
			       Form("Spectra from  Jet 100 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  //combined
  hpp_TrgObjComb[1] = new TH1F(Form("hpp_HLTComb_JetID_R%d_%s"      , radius, etaWidth), 
			       Form("Trig Combined Spectra R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_CombJetpT[1]  = new TH1F(Form("hpp_TrgCombTest_JetID_R%d_%s"             , radius, etaWidth), 
			       Form("Trig Combined Spectra KurtMethod R%d %s " , radius, etaWidth), 2000, 0, 2000);

  // for counting
  Long64_t NEvents_40=0, NEvents_60=0, NEvents_80=0, NEvents_100=0;//trigger counts
  Long64_t NEvents=0;//post-skim event counts
  Long64_t nentries = jetpp[0]->GetEntries();

  if(debugMode) nentries = 1000*(endfile-startfile+1);  
  if(fastDebugMode) nentries = 10*(endfile-startfile+1);

  int NEvents_w2Jets_woQACut=0,       NEvents_wo2Jets_woQACut=0;
  int NEvents_w2Jets_wQACut=0,        NEvents_wo2Jets_wQACut=0; 
  int NEvents_w2Jets_wQACut_wJetID=0, NEvents_wo2Jets_wQACut_wJetID=0;

  // EVENT LOOP
  std::cout<<"looping over "<<nentries<<" events"<<std::endl;  
  for(int nEvt = 0; nEvt < nentries; ++nEvt) {//event loop   

    // grab entry from each tree 
    for(int i=0;i<N_trees; ++i)jetpp[i]->GetEntry(nEvt);
    if(debugMode&&nEvt%1000==0)std::cout<<"nEvt = "<<nEvt<<std::endl;

    // skim/HiEvtAnalysis criteria
    if( pHBHENoiseFilter_F==0     || //skim 
        pBeamScrapingFilter_F==0  || //skim
        pprimaryvertexFilter_F==0 || //skim
        fabs(vz_F)>15              ) continue;
    
    // total prescale array
    int treePrescl[4]={ (jet40_p_F*jet40_l1seed_p_F), (jet60_p_F*jet60_l1seed_p_F), 
			(jet80_p_F*jet80_l1seed_p_F), (jet100_p_F*jet100_l1seed_p_F) };    
    // HLT decision array
    bool trgDec[4]={ (bool)jet40_F, (bool)jet60_F, (bool)jet80_F, (bool)jet100_F  };
    
    // prefer highPt l1 trigger for jet if it matches multiple l1jet pt criteria, then compute weight
    double triggerPt=0.;
    unsigned int trgObj40_size =trgObjpt_40->size();
    unsigned int trgObj60_size =trgObjpt_60->size();
    unsigned int trgObj80_size =trgObjpt_80->size();
    unsigned int trgObj100_size=trgObjpt_100->size();
    
    if(jet40_F){
      for(unsigned int itt=0; itt<trgObj40_size; ++itt){
	if(trgObjpt_40->at(itt) > triggerPt) triggerPt = trgObjpt_40->at(itt);
      }
    }
    if(jet60_F){
      for(unsigned int itt=0; itt<trgObj60_size; ++itt){
	if(trgObjpt_60->at(itt) > triggerPt) triggerPt = trgObjpt_60->at(itt);
      }
    }
    if(jet80_F){
      for(unsigned int itt=0; itt<trgObj80_size; ++itt){
    	if(trgObjpt_80->at(itt) > triggerPt) triggerPt = trgObjpt_80->at(itt);
      }
    }
    if(jet100_F){
      for(unsigned int itt=0; itt<trgObj100_size; ++itt){
    	if(trgObjpt_100->at(itt) > triggerPt) triggerPt = trgObjpt_100->at(itt);
      }
    }    

    double weight_eS = trigComb(trgDec, treePrescl, triggerPt);

    // fill vz histo
    hVz->Fill(vz_F, 1.0);
    hWeightedVz->Fill(vz_F, weight_eS);

    if(debugMode&&nEvt%250==0)std::cout <<"triggerPt ="<<triggerPt<<std::endl;
    if(debugMode&&nEvt%250==0)std::cout<<"jet 40, trgDec[0]  = "<< trgDec[0]  <<std::endl;
    if(debugMode&&nEvt%250==0)std::cout<<"jet 60, trgDec[1]  = "<< trgDec[1]  <<std::endl;
    if(debugMode&&nEvt%250==0)std::cout<<"jet 80, trgDec[2]  = "<< trgDec[2]  <<std::endl;
    if(debugMode&&nEvt%250==0)std::cout<<"jet 100, trgDec[3] = "<< trgDec[3]  <<std::endl<<std::endl;

    // check trigger decisions for events + exclusivity between them, count events
    bool is40  = false, is60  = false, is80  = false, is100 = false;
    if(trgDec[0] && triggerPt>=40 && triggerPt<60 ) is40  = true;    
    if(trgDec[1] && triggerPt>=60 && triggerPt<80 ) is60  = true;    
    if(trgDec[2] && triggerPt>=80 && triggerPt<100) is80  = true;    
    if(trgDec[3] && triggerPt>=100                ) is100 = true;    

    if(debugMode&&nEvt%250==0)std::cout<<"after checking trigger pt "<<std::endl;
    if(debugMode&&nEvt%250==0)std::cout<<"jet 40,  is40  = "<< is40  <<std::endl;
    if(debugMode&&nEvt%250==0)std::cout<<"jet 60,  is60  = "<< is60  <<std::endl;
    if(debugMode&&nEvt%250==0)std::cout<<"jet 80,  is80  = "<< is80  <<std::endl;
    if(debugMode&&nEvt%250==0)std::cout<<"jet 100, is100 = "<< is100  <<std::endl<<std::endl;

    if(is40)  NEvents_40++; 
    if(is60)  NEvents_60++; 
    if(is80)  NEvents_80++; 
    if(is100) NEvents_100++;
    if(true)  NEvents++;
    
    // JET LOOP 
    for(int jet = 0; jet<nref_F; ++jet){//jet loop

      if(debugMode)if(jet==0&&nref_F>1)  NEvents_w2Jets_woQACut++; 
      if(debugMode)if(jet==0&&nref_F<=1) NEvents_wo2Jets_woQACut++; 

      // pt is too low or not in barrel, toss jet
      if( fabs(eta_F[jet]) > 2 ||
          pt_F[jet] <= 15       ) continue;     
      
      if(debugMode)if(jet==0&&nref_F>1)  NEvents_w2Jets_wQACut++; 
      if(debugMode)if(jet==0&&nref_F<=1) NEvents_wo2Jets_wQACut++; 
      
      // jetQA noJetID
      float recpt = pt_F[jet];
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

      
      if(is40)  hpp_TrgObj40[0]->Fill(recpt, treePrescl[0]);//recpt here+below prescale weighted
      if(is60)  hpp_TrgObj60[0]->Fill(recpt, treePrescl[1]);
      if(is80)  hpp_TrgObj80[0]->Fill(recpt, treePrescl[2]);
      if(is100) hpp_TrgObj100[0]->Fill(recpt, treePrescl[3]);
      if(true)  hpp_CombJetpT[0]->Fill(recpt, weight_eS);//recpt here is weighted-weighted
      
      // dijet imbalance, pt ratio computation, lead+sublead jet
      if( jet==0 && nref_F>1 ) {
	float A_j=(pt_F[0]-pt_F[1])/(pt_F[0]+pt_F[1]);
	float x_j=(pt_F[1]/pt_F[0]); 	
	hJetQA[0][19]->Fill( A_j , weight_eS ); 
	hJetQA[0][20]->Fill( x_j , weight_eS ); 
      	//hJetQA[0][19]->Fill((float)(pt_F[0]-pt_F[1])/(pt_F[0]+pt_F[1]), weight_eS);
	//hJetQA[0][20]->Fill((float)(pt_F[1]/pt_F[0]), weight_eS);
      }

      // apply JetID
      if (chSum_F[jet]/recpt <= 0    ||
          neSum_F[jet]/recpt >= 0.99 || 
          phSum_F[jet]/recpt >= 0.99 ||
          chN_F[jet] <= 0             ) continue;
      
      if(debugMode)if(jet==0&&nref_F>1)  NEvents_w2Jets_wQACut_wJetID++; 
      if(debugMode)if(jet==0&&nref_F<=1) NEvents_wo2Jets_wQACut_wJetID++; 
      
      // jetQA wJetID
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
      
      if(is40 )  hpp_TrgObj40[1]->Fill(recpt, treePrescl[0]);
      if(is60 )  hpp_TrgObj60[1]->Fill(recpt, treePrescl[1]);
      if(is80 )  hpp_TrgObj80[1]->Fill(recpt, treePrescl[2]);
      if(is100) hpp_TrgObj100[1]->Fill(recpt, treePrescl[3]);
      if(true)  hpp_CombJetpT[1]->Fill(recpt, weight_eS);


      if( jet==0 && nref_F>1 ){
	float A_j=(pt_F[0]-pt_F[1])/(pt_F[0]+pt_F[1]);
	float x_j=(pt_F[1]/pt_F[0]); 	
	hJetQA[1][19]->Fill( A_j , weight_eS ); 
	hJetQA[1][20]->Fill( x_j , weight_eS ); 
      }
      
    }//end jet loop
  }//end event loop
  std::cout<<std::endl;

  // before JetID
  hpp_TrgObjComb[0]->Add(hpp_TrgObj100[0]);
  hpp_TrgObjComb[0]->Add(hpp_TrgObj80[0] );
  hpp_TrgObjComb[0]->Add(hpp_TrgObj60[0] );
  hpp_TrgObjComb[0]->Add(hpp_TrgObj40[0] );

  // after JetID
  hpp_TrgObjComb[1]->Add(hpp_TrgObj100[1]);
  hpp_TrgObjComb[1]->Add(hpp_TrgObj80[1] );
  hpp_TrgObjComb[1]->Add(hpp_TrgObj60[1] );
  hpp_TrgObjComb[1]->Add(hpp_TrgObj40[1] );
  
  // before JetID
  hpp_TrgObj40[0]->Print("base");
  hpp_TrgObj60[0]->Print("base");
  hpp_TrgObj80[0]->Print("base");
  hpp_TrgObj100[0]->Print("base");
  hpp_TrgObjComb[0]->Print("base");//diff between this
  hpp_CombJetpT[0]->Print("base");//and this?
  
  // after JetID
  hpp_TrgObj40[1]->Print("base");
  hpp_TrgObj60[1]->Print("base");
  hpp_TrgObj80[1]->Print("base");
  hpp_TrgObj100[1]->Print("base");
  hpp_TrgObjComb[1]->Print("base");
  hpp_CombJetpT[1]->Print("base");
  
  std::cout<<std::endl<<"Job Event-Loop Summary"<<std::endl;
  if(debugMode)std::cout<<"Total Number of Events read preskim  = "<<nentries<<std::endl;
  std::cout<<"Total Number of Events read          = "<<NEvents<<std::endl<<std::endl;
  
  if(debugMode)std::cout<<"NEvts w jets 0and1, before QA cut      = "<<NEvents_w2Jets_woQACut       <<std::endl;
  if(debugMode)std::cout<<"NEvts w jets 0and1, after QA cut     = "<<NEvents_w2Jets_wQACut        <<std::endl;
  if(debugMode)std::cout<<"NEvts w jets 0and1, after QA+ID cuts = "<<NEvents_w2Jets_wQACut_wJetID <<std::endl<<std::endl;
  if(debugMode)std::cout<<"NEvts w/o jets 0and1, before QA cut      = "<<NEvents_wo2Jets_woQACut       <<std::endl;
  if(debugMode)std::cout<<"NEvts w/o jets 0and1, after QA cut     = "<<NEvents_wo2Jets_wQACut        <<std::endl;
  if(debugMode)std::cout<<"NEvts w/o jets 0and1, after QA+ID cuts = "<<NEvents_wo2Jets_wQACut_wJetID <<std::endl<<std::endl;
  
  std::cout<<"Total Number of HLT40 events  = "<<NEvents_40<<std::endl;
  std::cout<<"Total Number of HLT60 events  = "<<NEvents_60<<std::endl;
  std::cout<<"Total Number of HLT80 events  = "<<NEvents_80<<std::endl;
  std::cout<<"Total Number of HLT100 events = "<<NEvents_100<<std::endl<<std::endl;
  
  std::cout<<"writing output file..."<<std::endl;
  fout->Write();

  if(debugMode)std::cout<<"misc clean up.."<<std::endl;
  trgObjpt_40->clear();
  trgObjpt_60->clear();
  trgObjpt_80->clear();
  trgObjpt_100->clear();

  std::cout<<std::endl<<"readFiles_ppData finished."<<std::endl;

  timer.Stop();
  std::cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<std::endl;
  std::cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<std::endl;

  return 0;
}// end readFiles_ppData


////// main //////
// acts as the frontend control for .exe file
int main(int argc, char *argv[]){
  
  // error, not enough arguments
  int rStatus = -1;
  if(argc!=8 && argc!=1){
    std::cout<<"for tests on default inputs, do..." <<std::endl;
    std::cout<<"./readFiles_ppData.exe";
    std::cout<<std::endl<<std::endl;
    std::cout<<"for actually running, do..."<<std::endl;
    std::cout<<"./readFiles_ppData.exe ";
    std::cout<<"<startFile> <endFile> <inputFileList> <outputFilename> ";
    std::cout<<"<jetRadius> <jetType> <debugMode>";
    std::cout<<std::endl<<std::endl;
    std::cout<<"rStatus="<<rStatus<<std::endl;
    return rStatus;
  }
  
  // good input, run
  rStatus=1;
  if(argc==1) rStatus = readFiles_ppData_Run2MinBias();
  else{//read input argument vector
    int startfile= atoi(argv[1]); int endfile= atoi(argv[2]); std::string inputFileList=argv[3]; 
    std::string outputFileName=argv[4]; int jetRadius= atoi(argv[5]); std::string jetType=argv[6];     
    bool debug=(bool)atoi(argv[7]);      
    rStatus = readFiles_ppData_Run2MinBias( startfile, endfile, inputFileList, outputFileName,
				jetRadius, jetType, debug);
  }
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}



