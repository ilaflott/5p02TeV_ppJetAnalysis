// custom header
#include "readFiles_ppMC.h"

//// readFiles_ppMC
// ---------------------------------------------------------------------------------------------------------------
int readFiles_ppMC(int startfile , int endfile , std::string inFilelist , std::string outfile ,
		   int radius , bool debugMode ){ 
  
  // for monitoring performance
  TStopwatch timer;  timer.Start();
  
  // debug mode settings+warnings
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl;
  const bool fastDebugMode = (debugMode)&&false; //if debugMode is off, fastDebugMode shouldn't be on
  if(fastDebugMode)std::cout<<"fastDebugMode is ON"<<std::endl;
  
  // basic info the screen
  std::cout<<std::endl<<"///////////////////"<<std::endl;
  std::cout<<"reading files "<< startfile << " to " << endfile<<std::endl;
  std::cout<<"radius =  " << radius<<std::endl;
  if(debugMode)std::cout<<"N trees = " << N_trees<<std::endl;
  if(debugMode)std::cout<<"N vars  = " << N_vars<<std::endl;
  std::cout<<"///////////////////"<<std::endl<<std::endl;
  
  // plot settings for later
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  // form PF/CalojetTreeNames, initialize tree name array
  std::string trees[N_trees];
  std::string PFJetTreeName="ak"+std::to_string(radius)+"PFJetAnalyzer/t";
  std::string CaloJetTreeName="ak"+std::to_string(radius)+"CaloJetAnalyzer/t";
  trees[0]=PFJetTreeName;
  trees[1]=CaloJetTreeName;
  if(debugMode)std::cout<<"looking at jetTrees "<<trees[0]<<" and "<<trees[1]<<std::endl;
  for(int i=2;i<N_trees;++i)trees[i]=treeNames[i];
  
  
  //declare tchain, ifstream, and variables for the loop over the filelist
  TChain* jetpp[N_trees];
  for(int t = 0;t<N_trees;++t)  jetpp[t] = new TChain( trees[t].data() );
  std::ifstream instr_Forest(inFilelist.c_str(),std::ifstream::in);
  std::cout<<"filelist used is "<<inFilelist<<std::endl;
  std::string filename_Forest;    std::string lastFileAdded="";   bool treesAdded=false;
  for(int ifile = 0; ifile<=endfile; ++ifile){//file loop

    // loop to starting file, check for end of filelist
    instr_Forest>>filename_Forest; 
    if(ifile<startfile){ 
      lastFileAdded=filename_Forest; 
      continue; } // file # under minimum of range desired
    if(filename_Forest==lastFileAdded){ 
      std::cout<<"end of filelist!"<<std::endl; 
      break; } // file # is greater than # files in list

    // debug info just in case
    std::cout<<"adding file #"<<ifile;     
    if(debugMode)std::cout<<" "<<filename_Forest;  
    std::cout<<" to TChain(s)"<<std::endl;

    // add file to each TChain in the jetpp array
    for(int t = 0;t<N_trees;++t) jetpp[t]->Add(filename_Forest.c_str()); 
    lastFileAdded=filename_Forest;
    treesAdded=true; 
  }// end file loop
  assert(treesAdded);// avoid segfault later, just in case
  
  // friend all trees to akPFJetAnalyzer
  for(int t = 1;t<N_trees;++t)jetpp[0]->AddFriend( jetpp[t] );
  
  // declare variables/arrays for input files
  // akPFJetAnalyzer
  //jet variables
  int nref_F;  float pthat_F;
  float pt_F[1000];     float eta_F[1000];  float phi_F[1000]; 
  float rawpt_F[1000];  float jtpu_F[1000];
  //MC jet variables
  int subid_F[1000];    int refparton_F[1000];
  float refpt_F[1000];  float refeta_F[1000];   float refdrjt_F[1000];
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
  
  //// akCaloJetAnalyzer
  ////jets
  //int nref_F;  
  //float pt_F[1000];     float eta_F[1000];  float phi_F[1000]; 
  //float rawpt_F[1000];  float jtpu_F[1000];
  
  // HiEvtAnalyzer
  ULong64_t evt_F;   UInt_t run_F;   UInt_t lumi_F;   float vz_F;
  
  // skimanalysis
  int pBeamScrapingFilter_F;   int pHBHENoiseFilter_F;   int pprimaryvertexFilter_F;
  
  // set branch addresses for the input forests
  // ak${radius}PFJetAnalyzer
  //jets
  jetpp[0]->SetBranchAddress("nref",&nref_F);
  jetpp[0]->SetBranchAddress("pthat",&pthat_F);
  jetpp[0]->SetBranchAddress("jtpt",&pt_F);
  jetpp[0]->SetBranchAddress("jteta",&eta_F);
  jetpp[0]->SetBranchAddress("jtphi",&phi_F);
  jetpp[0]->SetBranchAddress("rawpt",&rawpt_F);
  jetpp[0]->SetBranchAddress("jtpu",&jtpu_F);
  //MC jet variables
  jetpp[0]->SetBranchAddress( "subid"	, &subid_F     );
  jetpp[0]->SetBranchAddress( "refparton_flavor" , &refparton_F );
  jetpp[0]->SetBranchAddress( "refpt"	, &refpt_F     );
  jetpp[0]->SetBranchAddress( "refeta"  , &refeta_F    );
  jetpp[0]->SetBranchAddress( "refdrjt"	, &refdrjt_F   );
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
  
  //// ak${radius}CaloJetAnalyzer
  ////jets
  //jetpp[1]->SetBranchAddress("nref",&calonref_F);
  //jetpp[1]->SetBranchAddress("jtpt",&calopt_F);
  //jetpp[1]->SetBranchAddress("jteta",&eta_F);
  //jetpp[1]->SetBranchAddress("jtphi",&phi_F);
  //jetpp[1]->SetBranchAddress("rawpt",&rawpt_F);
  //jetpp[1]->SetBranchAddress("jtpu",&jtpu_F);
  
  // hiEvtAnalyzer
  jetpp[2]->SetBranchAddress("evt",&evt_F);
  jetpp[2]->SetBranchAddress("run",&run_F);
  jetpp[2]->SetBranchAddress("lumi",&lumi_F);
  jetpp[2]->SetBranchAddress("vz",&vz_F);
  
  // skimanalysis
  jetpp[3]->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_F);
  jetpp[3]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilter_F);
  jetpp[3]->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryvertexFilter_F);
  
  // Declare the output File and the necessary histograms after that:
  std::cout<<"opening output file "<<outfile<<std::endl;
  TFile *fout = new TFile(outfile.c_str(),"RECREATE");
  fout->cd();
  
  // histos for weights+weighting
  TH1F *hpthat = new TH1F("hpthat","",1000,0,1000);
  TH1F *hWeightedpthat = new TH1F("hWeightedpthat","",1000,0,1000);
  TH1F *hVz = new TH1F("hVz","",200, -20, 20);
  TH1F *hpthatWeightedVz = new TH1F("hpthatWeightedVz","",200, -20, 20);
  TH1F *hvzWeightedVz = new TH1F("hvzWeightedVz","",200, -20, 20);
  TH1F *hWeightedVz = new TH1F("hWeightedVz","",200, -20, 20);

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
  
  // for unfolding + closure tests ;ate
  //3d plot of genpt v rec pt -> matrix
  TH1F* hpp_gen    = new TH1F( Form("hpp_gen_R%d_%s",radius,etaWidth), Form("Gen refpt R%d %s ",radius,etaWidth),1000,0,1000);
  TH1F* hpp_reco   = new TH1F( Form("hpp_reco_R%d_%s",radius,etaWidth), Form("Reco jtpt R%d %s ",radius,etaWidth),1000,0,1000);
  TH2F* hpp_matrix = new TH2F( Form("hpp_matrix_R%d_%s",radius,etaWidth), 
			       Form("Matrix refpt jtpt R%d %s ; reco pT; gen pT",radius,etaWidth), 1000,0,1000,1000,0,1000);

  //3d plot of "indep." sample genpt rec pt -> matrx
  TH1F* hpp_mcclosure_gen = new TH1F(Form("hpp_mcclosure_gen_R%d_%s",radius,etaWidth),
				     Form("gen spectra for unfolding mc closure test R%d %s ",radius,etaWidth), 1000, 0, 1000);
  TH1F* hpp_mcclosure_data = new TH1F(Form("hpp_mcclosure_data_R%d_%s",radius,etaWidth),
				      Form("data for unfolding mc closure test R%d %s ",radius,etaWidth), 1000, 0, 1000);
  TH2F* hpp_mcclosure_matrix = new TH2F(Form("hpp_mcclosure_matrix_R%d_%s",radius,etaWidth),
					Form("Matrix for mcclosure refpt jtpt R%d %s;reco pT; gen pT ",radius,etaWidth),
					1000,0,1000, 1000,0,1000);

  //the datasample the unfoldig is "closure-tested" on
  TH1F* hpp_mcclosure_data_train = new TH1F(Form("hpp_mcclosure_data_train_R%d_%s",radius,etaWidth),
					    Form("data_train for unfolding mc closure test R%d %s ",radius,etaWidth), 1000, 0, 1000);  

  // JEC
  TH3F * hJEC= new TH3F("hJEC",";raw p_{T};#eta;JEC",500, 0, 500, 200, -5, +5, 300, 0, 5);

  TH1F * hJEC_check[nbins_JEC_ptbins][nbins_eta];
  for(int x = 0; x<nbins_JEC_ptbins; ++x){ 
    for(int y = 0; y<nbins_eta; ++y){
      hJEC_check[x][y] = new TH1F(Form("hJEC_check_ptbin%d_etabin%d",x,y),
				  Form("rawpt/genpt %2.0f < genpt < %2.0f, %2.4f < geneta < %2.4f",JEC_ptbins[x], JEC_ptbins[x+1], etabins[y], etabins[y+1]),
				  100, 0, 3);    }  }

  // JER
  TH1F * hJER[nbins_pt];
  for(int bin = 0; bin<nbins_pt; ++bin) hJER[bin] = new TH1F(Form("hJER_%d_pt_%d", ptbins[bin], ptbins[bin+1]),"",100, 0, 2);

  TH1F* hJER_eta_30pt50[nbins_eta]; TH1F* hJER_eta_150pt200[nbins_eta];
  for(int bin = 0; bin<nbins_eta; ++bin){
    hJER_eta_30pt50[bin] = new TH1F(Form("hJER_etabin%d_30_pt_50", bin),
				    Form("rawpt/genpt 30 < genpt < 50, %2.4f < geneta < %2.4f", etabins[bin], etabins[bin+1]),100, 0, 2);
    hJER_eta_150pt200[bin] = new TH1F(Form("hJER_etabin%d_150_pt_200", bin),
				      Form("rawpt/genpt 150 < genpt < 200, %2.4f < geneta < %2.4f", etabins[bin], etabins[bin+1]),100, 0, 2);  }


  // for counting trigger events
  Long64_t NEvents_40=0, NEvents_60=0, NEvents_80=0, NEvents_100=0;//trigger counts
  Long64_t NEvents=0;//post-skim event counts

  // for debugging
  int NEvents_w2Jets_woQACut=0,       NEvents_wo2Jets_woQACut=0;
  int NEvents_w2Jets_wQACut=0,        NEvents_wo2Jets_wQACut=0; 
  int NEvents_w2Jets_wQACut_wJetID=0, NEvents_wo2Jets_wQACut_wJetID=0;

  // evt  vz weights from Raghav
  TF1 *fVzPP = new TF1("fVzPP","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVzPP->SetParameters(1.05602e00,5.74688e-03,-3.37288e-03,-1.44764e-05,8.59060e-07);  // Raghav's weights
  //fVzPP->SetParameters(+0.941, -0.0173, +3.23e-3, +3.61e-6, -1.04e-5);  // Anna's weights

  // a little prep for pt/eta bin searching later
  // if extra speed is needed (i dont think it is), can implement at a later time
  const float minAbsJetEta=2.0, minJetPt=15.0;
  //int etabin = -1;  //int binx = -1;  //int ptbin = -1;
  
  // EVENT LOOP
  Long64_t nentries=jetpp[0]->GetEntries();  
  if(debugMode) nentries=1000*(endfile-startfile+1);  
  if(fastDebugMode) nentries=10*(endfile-startfile+1);
  
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

    // this was in Raghav's code, looked important to him, but is it important to me? can't tell.   
    //// aply proxy for collision eventselection : position.Rho <= 2 && tracksSize >= 2
    //double rho = TMath::Sqrt((vx_F[0]*vx_F[0]) + (vy_F[0]*vy_F[0]));
    //if(rho > 2 || nTrk < 2) continue;

    // compute weights
    double vzWeight=1;           
    vzWeight = fVzPP->Eval(vz_F);
    
    double evtPthatWeight=0.;    
    for( int i=0; i<nbins_pthat && pthat_F>=pthatbins[i]; i++ ){ evtPthatWeight=pthatWeights[i]; }
 
    double weight_eS=1;          
    //weight_eS = trigComb(trgDec, treePrescl, triggerPt);    
    double finalWeight=evtPthatWeight*vzWeight*weight_eS;

    // fill (un)weighted evt histos
    hVz->Fill(vz_F, 1.);
    hpthatWeightedVz->Fill(vz_F, evtPthatWeight);
    hvzWeightedVz->Fill(vz_F, vzWeight);
    hWeightedVz->Fill(vz_F, finalWeight);
    hpthat->Fill(pthat_F, 1.);
    hWeightedpthat->Fill(pthat_F, finalWeight);

    // check trigger decisions for events + exclusivity between them, count events
    bool is40  = false, is60  = false, is80  = false, is100 = false;
    if(is40)  NEvents_40++; 
    if(is60)  NEvents_60++; 
    if(is80)  NEvents_80++; 
    if(is100) NEvents_100++;
    if(true)  NEvents++;
    
    // JET LOOP 
    for(int jet = 0; jet<nref_F; ++jet){

      if(debugMode){
	if(jet==0&&nref_F>1)  NEvents_w2Jets_woQACut++; 
	if(jet==0&&nref_F<=1) NEvents_wo2Jets_woQACut++;}

      //// MCjet specific QA cuts
      //if(subid_F[jet] != 0) continue;
      //if(pt_F[jet] > 3 * pthat_F) continue;
      //if(refdrjt_F[jet] > 0.3) continue;

      // std jet cuts
      if( fabs(eta_F[jet]) > minAbsJetEta ||
          pt_F[jet] <= minJetPt              ) continue;     
      
      if(debugMode){
	if(jet==0&&nref_F>1)  NEvents_w2Jets_wQACut++; 
	if(jet==0&&nref_F<=1) NEvents_wo2Jets_wQACut++; }
            
      // jetQA noJetID
      float genpt = refpt_F[jet];
      float rawpt = rawpt_F[jet];
      float recpt = pt_F[jet];
      float geneta = refeta_F[jet];
      //float absgeneta = fabs(geneta);
      
      // JER/JEC
      hJEC->Fill(rawpt, eta_F[jet], (float)(recpt/rawpt));

      int etabin = -1;
      for(int bin = 0; bin<nbins_eta; ++bin){
        if(geneta > etabins[bin]) etabin = bin;
      }
      if(etabin == -1) continue;
      if(genpt >= 30 && genpt<50) hJER_eta_30pt50[etabin]->Fill((float)recpt/genpt, finalWeight);
      if(genpt >= 150 && genpt<200) hJER_eta_150pt200[etabin]->Fill((float)recpt/genpt, finalWeight);

      int binx = -1; //genpt JEC bin
      for(int bin = 0; bin<nbins_JEC_ptbins; ++bin){
        if(genpt > JEC_ptbins[bin]) binx = bin;
      }
      if(binx == -1) continue;
      hJEC_check[binx][etabin]->Fill((float)rawpt/genpt);

      int ptbin = -1; //another genpt binnning
      for(int bin = 0; bin<nbins_pt; ++bin){
        if(genpt > ptbins[bin]) ptbin = bin;
      }
      if(ptbin == -1) continue;
      hJER[ptbin]->Fill((float)recpt/genpt, finalWeight);

      // for unfolding
      hpp_gen->Fill(genpt, finalWeight);
      hpp_reco->Fill(recpt, finalWeight);
      hpp_matrix->Fill(recpt, genpt, finalWeight);

      // for MC closure tests
      if(nEvt%2 == 0){
        hpp_mcclosure_data->Fill(recpt, finalWeight);
      }else {
        hpp_mcclosure_matrix->Fill(recpt, genpt, finalWeight);
        hpp_mcclosure_gen->Fill(genpt, finalWeight);
        hpp_mcclosure_data_train->Fill(recpt, finalWeight);
      }

      // jet QA, before jetID cuts
      hJetQA[0][0]->Fill(recpt, finalWeight);
      hJetQA[0][1]->Fill(rawpt_F[jet], finalWeight);
      hJetQA[0][2]->Fill(eta_F[jet], finalWeight);
      hJetQA[0][3]->Fill(phi_F[jet], finalWeight);
      hJetQA[0][4]->Fill(trkMax_F[jet]/recpt, finalWeight);
      hJetQA[0][5]->Fill(trkSum_F[jet]/recpt, finalWeight);
      hJetQA[0][6]->Fill(trkHardSum_F[jet]/recpt, finalWeight);
      hJetQA[0][7]->Fill(chMax_F[jet]/recpt, finalWeight);
      hJetQA[0][8]->Fill(chSum_F[jet]/recpt, finalWeight);
      hJetQA[0][9]->Fill(chHardSum_F[jet]/recpt, finalWeight);
      hJetQA[0][10]->Fill(phMax_F[jet]/recpt, finalWeight);
      hJetQA[0][11]->Fill(phSum_F[jet]/recpt, finalWeight);
      hJetQA[0][12]->Fill(phHardSum_F[jet]/recpt, finalWeight);
      hJetQA[0][13]->Fill(neMax_F[jet]/recpt, finalWeight);
      hJetQA[0][14]->Fill(neSum_F[jet]/recpt, finalWeight);
      hJetQA[0][15]->Fill(eMax_F[jet]/recpt, finalWeight);
      hJetQA[0][16]->Fill(eSum_F[jet]/recpt, finalWeight);
      hJetQA[0][17]->Fill(muMax_F[jet]/recpt, finalWeight);
      hJetQA[0][18]->Fill(muSum_F[jet]/recpt, finalWeight);
      if( jet==0 && nref_F>1 ) {
	float A_j=(pt_F[0]-pt_F[1])/(pt_F[0]+pt_F[1]);
	float x_j=(pt_F[1]/pt_F[0]); 	
	hJetQA[0][19]->Fill( A_j , finalWeight ); 
	hJetQA[0][20]->Fill( x_j , finalWeight ); 
      }

      // apply JetID cut
      if (chSum_F[jet]/recpt <= 0    ||
          neSum_F[jet]/recpt >= 0.99 || 
          phSum_F[jet]/recpt >= 0.99 ||
          chN_F[jet] <= 0             ) continue;
      
      if(debugMode){ 
	if(jet==0&&nref_F>1)  NEvents_w2Jets_wQACut_wJetID++; 
	if(jet==0&&nref_F<=1) NEvents_wo2Jets_wQACut_wJetID++; }
      
      // jet QA, after jetID cuts
      hJetQA[1][0]->Fill(recpt, finalWeight);
      hJetQA[1][1]->Fill(rawpt_F[jet], finalWeight);
      hJetQA[1][2]->Fill(eta_F[jet], finalWeight);
      hJetQA[1][3]->Fill(phi_F[jet], finalWeight);
      hJetQA[1][4]->Fill(trkMax_F[jet]/recpt, finalWeight);
      hJetQA[1][5]->Fill(trkSum_F[jet]/recpt, finalWeight);
      hJetQA[1][6]->Fill(trkHardSum_F[jet]/recpt, finalWeight);
      hJetQA[1][7]->Fill(chMax_F[jet]/recpt, finalWeight);
      hJetQA[1][8]->Fill(chSum_F[jet]/recpt, finalWeight);
      hJetQA[1][9]->Fill(chHardSum_F[jet]/recpt, finalWeight);
      hJetQA[1][10]->Fill(phMax_F[jet]/recpt, finalWeight);
      hJetQA[1][11]->Fill(phSum_F[jet]/recpt, finalWeight);
      hJetQA[1][12]->Fill(phHardSum_F[jet]/recpt, finalWeight);
      hJetQA[1][13]->Fill(neMax_F[jet]/recpt, finalWeight);
      hJetQA[1][14]->Fill(neSum_F[jet]/recpt, finalWeight);
      hJetQA[1][15]->Fill(eMax_F[jet]/recpt, finalWeight);
      hJetQA[1][16]->Fill(eSum_F[jet]/recpt, finalWeight);
      hJetQA[1][17]->Fill(muMax_F[jet]/recpt, finalWeight);
      hJetQA[1][18]->Fill(muSum_F[jet]/recpt, finalWeight);
      if( jet==0 && nref_F>1 ){
	float A_j=(pt_F[0]-pt_F[1])/(pt_F[0]+pt_F[1]);
	float x_j=(pt_F[1]/pt_F[0]); 	
	hJetQA[1][19]->Fill( A_j , finalWeight ); 
	hJetQA[1][20]->Fill( x_j , finalWeight ); 
      }
    }//end jet loop
  }//end event loop


  std::cout<<std::endl;
  std::cout<<std::endl<<"Job Event-Loop Summary"<<std::endl;
  if(debugMode)std::cout<<"Total Number of Events read preskim  = "<<nentries<<std::endl;
  std::cout<<"Total Number of Events read          = "<<NEvents<<std::endl<<std::endl;
  
  if(debugMode){std::cout<<"NEvts w jets 0and1, before QA cut      = "<<NEvents_w2Jets_woQACut       <<std::endl;
    std::cout<<"NEvts w jets 0and1, after QA cut     = "<<NEvents_w2Jets_wQACut        <<std::endl;
    std::cout<<"NEvts w jets 0and1, after QA+ID cuts = "<<NEvents_w2Jets_wQACut_wJetID <<std::endl<<std::endl;
    std::cout<<"NEvts w/o jets 0and1, before QA cut      = "<<NEvents_wo2Jets_woQACut       <<std::endl;
    std::cout<<"NEvts w/o jets 0and1, after QA cut     = "<<NEvents_wo2Jets_wQACut        <<std::endl;
    std::cout<<"NEvts w/o jets 0and1, after QA+ID cuts = "<<NEvents_wo2Jets_wQACut_wJetID <<std::endl<<std::endl;}
  
  std::cout<<"Total Number of \"HLT40\" events  = "<<NEvents_40<<std::endl;
  std::cout<<"Total Number of \"HLT60\" events  = "<<NEvents_60<<std::endl;
  std::cout<<"Total Number of \"HLT80\" events  = "<<NEvents_80<<std::endl;
  std::cout<<"Total Number of \"HLT100\" events = "<<NEvents_100<<std::endl<<std::endl;
  
  std::cout<<"writing output file..."<<std::endl;
  fout->Write();

  std::cout<<std::endl<<"readFiles_ppMC finished."<<std::endl;

  timer.Stop();
  std::cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<std::endl;
  std::cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<std::endl;

  return 0;
}// end readFiles_ppMC


////// main //////
// acts as the frontend control for .exe file
int main(int argc, char *argv[]){
  
  // error, not enough arguments
  int rStatus = -1;
  if(argc!=7 && argc!=1){
    std::cout<<"for tests on default inputs, do..." <<std::endl;
    std::cout<<"./readFiles_ppMC.exe";
    std::cout<<std::endl<<std::endl;
    std::cout<<"for actually running, do..."<<std::endl;
    std::cout<<"./readFiles_ppMC.exe ";
    std::cout<<"<startFile> <endFile> <inputFileList> <outputFilename> ";
    std::cout<<"<jetRadius> <debugMode>";
    //std::cout<<"<weightsFile> <jetRadius> <debugMode>";
    std::cout<<std::endl<<std::endl;
    std::cout<<"rStatus="<<rStatus<<std::endl;
    return rStatus;
  }
  
  // good input, run
  rStatus=1;
  if(argc==1) rStatus = readFiles_ppMC();
  else{//read input argument vector
    int startfile= atoi(argv[1]); int endfile= atoi(argv[2]); std::string inputFileList=argv[3]; std::string outputFileName=argv[4];
    int jetRadius= atoi(argv[5]); bool debug=(bool)atoi(argv[6]);      
    rStatus = readFiles_ppMC( startfile, endfile, inputFileList, outputFileName,
			      jetRadius, debug);
  }
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}
