#include "DijetResponse.h"


//// deriveResponse
// ---------------------------------------------------------------------------------------------------------------
// writes nbins_pt*nbins_eta "Response Histos" in each pt,eta bin
int deriveResponse(int startfile, int endfile, std::string infile_Forest, std::string outfile, // I/O
		   int radius, bool isMC, bool debugMode){ // other

  // set random seed
  srand(time( NULL ));
  
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl;
  const bool fastDebugMode = (debugMode&&false); //if debugMode is off, fastDebugMode shouldn't be on
  if(fastDebugMode)std::cout<<"fastDebugMode is ON"<<std::endl;

  std::cout<<std::endl<<"///////////////////"<<std::endl;
  std::cout<<"reading files "<< startfile << " to " << endfile<<std::endl;
  if(debugMode)std::cout<<"radius =  " << radius<<std::endl;
  if(debugMode)std::cout<<"isMC   = " << isMC <<std::endl;
  std::cout<<"///////////////////"<<std::endl<<std::endl;  

  // for hRelResponse
  int nEntriesA[nbins_pt][nbins_eta];
  double avgA[nbins_pt][nbins_eta];
  TH1F *avgAHisto[nbins_pt][nbins_eta];
  TH1D *h_nEntriesA[nbins_pt][nbins_eta];
  TH1D *h_avgA[nbins_pt][nbins_eta];
  
  // for hMPFResponse
  int nEntriesB[nbins_pt][nbins_eta];  
  double avgB[nbins_pt][nbins_eta];
  TH1F *avgBHisto[nbins_pt][nbins_eta];  
  TH1D *h_nEntriesB[nbins_pt][nbins_eta];
  TH1D *h_avgB[nbins_pt][nbins_eta];
  
  // for hMPFAbsPhoResponse
  int nEntriesAbsPho[nbins_pt][nbins_eta];  
  double avgAbsPho[nbins_pt][nbins_eta];
  TH1F *avgAbsPhoHisto[nbins_pt][nbins_eta];  
  TH1D *h_nEntriesAbsPho[nbins_pt][nbins_eta];
  TH1D *h_avgAbsPho[nbins_pt][nbins_eta];
  
  // book histos and initialize arrays
  for(int i=0; i<nbins_pt; i++){//pt loop
    for(int j=0; j<nbins_eta; j++){//eta loop
      // avgA
      nEntriesA[i][j] = 0;
      avgA[i][j]      = 0;      
      avgAHisto[i][j]   = new TH1F( Form( "avgAHisto_pt%d_eta%d",i,j),"",50,-1,1);
      h_avgA[i][j]      = new TH1D( Form( "h_avgA_%d_%d",i,j)        ,"",2, 0, 2);
      h_nEntriesA[i][j] = new TH1D( Form( "h_nEntriesA_%d_%d",i,j)   ,"",2, 0, 2);
      // avgB
      nEntriesB[i][j] = 0;
      avgB[i][j]      = 0;      
      avgBHisto[i][j]   = new TH1F( Form( "avgBHisto_pt%d_eta%d",i,j),"",50,-1,1);
      h_avgB[i][j]      = new TH1D( Form( "h_avgB_%d_%d",i,j)        ,"",2, 0, 2);
      h_nEntriesB[i][j] = new TH1D( Form( "h_nEntriesB_%d_%d",i,j)   ,"",2, 0, 2);
      // avgAbsPho
      nEntriesAbsPho[i][j] = 0;
      avgAbsPho[i][j]      = 0;      
      avgAbsPhoHisto[i][j]   = new TH1F( Form( "avgAbsPhoHisto_pt%d_eta%d",i,j), "", 50, -1, 1);
      h_avgAbsPho[i][j]      = new TH1D( Form( "h_avgAbsPho_%d_%d",i,j),         "", 2,  0,  2);
      h_nEntriesAbsPho[i][j] = new TH1D( Form( "h_nEntriesAbsPho_%d_%d",i,j),    "", 2,  0,  2);
    }//end eta loop
  }//end pt loop
  
  // open filelist, loop to the starting file
  std::cout<<"opening filelist "<<infile_Forest<<std::endl;
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;  
  for(int ifile = 0;ifile<startfile;++ifile) instr_Forest >> filename_Forest;
  
  // loop over specified files 
  std::cout<<"looping over specified files..."<<std::endl;
  for(int ifile = startfile; ifile<=endfile; ++ifile){//file loop

    // grab+open file
    instr_Forest>>filename_Forest;
    if(debugMode)std::cout<<"opening trees from file "<<filename_Forest<<std::endl;    
    TFile *fin = TFile::Open(filename_Forest.c_str());

    // get trees
    std::string jetTreeName="ak"+std::to_string(radius)+"PFJetAnalyzer/t";
    TTree *jtTree  = (TTree*)fin->Get( jetTreeName.c_str()       );
    TTree *pfTree  = (TTree*)fin->Get( "pfcandAnalyzer/pfTree"   );
    TTree *phoTree = (TTree*)fin->Get( "ggHiNtuplizer/EventTree" );

    //jtTree
    int nref; 
    float pt_F[1000],  eta_F[1000],  phi_F[1000];
    float rawpt_F[1000], m_F[1000], eSum[1000], phoSum[1000];
    //pfTree
    int nPFpart_;
    int pfId_[10000];				       
    float pfPt_[10000], pfEta_[10000], pfPhi_[10000];
    //std::vector<int> *pfId_;                          //HighPtJet80 dataset only         
    //std::vector<float> *pfPt_=0, *pfEta_=0, *pfPhi_=0;	    

    //phoTree
    int nPho;
    std::vector<float> *phoEt=0, *phoEta=0, *phoPhi=0; 

    // branch addresses
    jtTree->SetBranchAddress("nref",      &nref   );
    jtTree->SetBranchAddress("jtpt",      &pt_F    );
    jtTree->SetBranchAddress("jteta",     &eta_F   );
    jtTree->SetBranchAddress("jtphi",     &phi_F   );
    jtTree->SetBranchAddress("jtm",       &m_F     );
    jtTree->SetBranchAddress("rawpt",     &rawpt_F );
    jtTree->SetBranchAddress("eSum",      &eSum    );  
    jtTree->SetBranchAddress("photonSum", &phoSum  );  

    pfTree->SetBranchAddress("nPFpart", &nPFpart_);
    pfTree->SetBranchAddress("pfId",    &pfId_);
    pfTree->SetBranchAddress("pfPt",    &pfPt_);
    pfTree->SetBranchAddress("pfEta",   &pfEta_);
    pfTree->SetBranchAddress("pfPhi",   &pfPhi_);

    phoTree->SetBranchAddress("nPho",   &nPho);
    phoTree->SetBranchAddress("phoEt",  &phoEt);
    phoTree->SetBranchAddress("phoEta", &phoEta);
    phoTree->SetBranchAddress("phoPhi", &phoPhi);

    // grab jet tree entries, loop over jets
    const double jtptMin = 40.;
    int jetEvents=jtTree->GetEntries();
    if(debugMode) std::cout<<"jetEvents = "<< jetEvents <<std::endl;
    if(debugMode) std::cout<<"beginning jet loop..."<< std::endl;    
    for(int ientry=0; ientry<jetEvents; ientry++){//jet loop

      // get a jet entry from the trees
      jtTree->GetEntry(ientry);
      phoTree->GetEntry(ientry);
      pfTree->GetEntry(ientry);
      if(nPFpart_>10000) std::cout<<"warning! nPF: "<<nPFpart_<<std::endl;
      if(nPFpart_>10000) nPFpart_=10000;

      // avgAbsPho histo stuff //
      // doesn't quite come out with sensible numbers at the moment
      // find response of detector to jet with particle-flow photon
      bool foundPhoton=false;
      for(int i=0; i<nPho; i++){// photon loop
	// check photon quality, if we found a good one, loop over jets!
	float photonEt=phoEt->at(i); 
	float photonEta=phoEta->at(i); 
	float photonPhi=phoPhi->at(i);
	bool foundPair=false;
	if( photonEt<=30 ||
	    abs(photonEta)>=1.4442 )continue;
	if(debugMode&&!foundPhoton&&ientry%100==0)std::cout << "at least one photon for this jet..." << std::endl; 
	if(!foundPhoton)foundPhoton=true;
	for(int ijet=0; ijet<nref-1; ijet++){// loop over jets in event
	  if( pt_F[ijet]<jtptMin ||         // jetpt cut
              abs(photonPhi-phi_F[ijet])<=2.95 || //angular separation between the +jet
              pt_F[ijet+1]/photonEt>=0.2 ) continue;
	  if(debugMode&&!foundPair&&ientry%100==0)std::cout << "at least one jet-photon pair found" << std::endl; 
	  if(!foundPair)foundPair=true;
	  // missEt+findBin for this photon/jet pair event
	  TLorentzVector missEt = findMissEt(nPFpart_, pfId_, pfPt_, pfEta_, pfPhi_,
					     nref, pt_F, eta_F, phi_F, 		  
					     m_F, rawpt_F, eSum, phoSum);           
	  int etaBin = findBin(photonEta,nbins_eta,xbins_eta);
	  int ptBin = findBin(photonEt,nbins_pt,xbins_pt);
	  // dot products
	  TLorentzVector phoVec(photonEt,photonEta,photonPhi, 0.);
	  double missEtDotPhoVec = missEt.Dot(phoVec);
	  double phoVecMagSq = phoVec.Dot(phoVec);
	  double absPho = 1+missEtDotPhoVec/phoVecMagSq;
	  // fill
	  nEntriesAbsPho[ptBin][etaBin]++;
	  avgAbsPho[ptBin][etaBin] += ( absPho );
	  h_nEntriesAbsPho[ptBin][etaBin]->Fill(1);
	  h_avgAbsPho[ptBin][etaBin]->Fill( absPho );
	  avgAbsPhoHisto[ptBin][etaBin]->Fill( absPho );
	  if(debugMode&&false)std::cout<<"for photon #"<<i<<" and jet #"<<ijet<<std::endl;
	  if(debugMode&&false)std::cout<<"in bins [ptBin][etaBin]=]["<<ptBin<<"]["<<etaBin<<"]...."<<std::endl;
	  if(debugMode&&false)std::cout<<"missEt dot phoVec = "<< missEtDotPhoVec<<std::endl;
	  if(debugMode&&false)std::cout<<"phoVec mag^2 = "     << phoVecMagSq<<std::endl;
	  if(debugMode&&false)std::cout<<"absPho = "           << absPho<<std::endl;
	  if(debugMode&&false)std::cout<<"nEntriesAbsPho = "<< nEntriesAbsPho[ptBin][etaBin]<<std::endl;
	}//end jet loop
      }//end photon loop
      // end avgAbsPho histo stuff //

      // avgA and avgB histo stuff //
      // Need lead+sublead jets
      if(nref<2)continue;

      // randomly assign lead/sublead jet depending on eta distribution
      int rJet, pJet;
      if( abs(eta_F[0])<1.3 && abs(eta_F[1])<1.3 ){ 
	int randInt=rand();
	//if(debugMode)std::cout<< "rand() gave "<<randInt <<std::endl;
	pJet = randInt%2;  rJet = (randInt+1)%2;
      } 
      else if( abs(eta_F[0])<1.3 && abs(eta_F[1])<5. ) { pJet = 1; rJet = 0; }
      else if( abs(eta_F[1])<1.3 && abs(eta_F[0])<5. ) { pJet = 0; rJet = 1; }
      else continue;// lead/sublead jets no good, next event!
      //if(debugMode)std::cout<< "rJet is "<<rJet<<" and pJet is "<<pJet <<std::endl;
      
      // minjetpt cut
      if(pt_F[rJet]<jtptMin) continue;
      
      // angular separation cut
      double dphi = abs(phi_F[rJet]-phi_F[pJet]);
      if(dphi>(2*3.14159)) dphi-=(2*3.14159);
      const double dphiMin=2.7;
      if(dphi < dphiMin) continue;
      
      // rejection criteria for events w/ third jets
      double avgPt = 0.5*(pt_F[rJet]+pt_F[pJet]);
      if( nref>2 && 
	  pt_F[2]/avgPt>0.2 )continue;
      
      // find the bin the jet belongs to
      int etaBin = findBin(eta_F[pJet],nbins_eta,xbins_eta);
      int ptBin = findBin(avgPt,nbins_pt,xbins_pt);
      
      // avgA histo stuff //
      // uses dijet asymmetry
      nEntriesA[ptBin][etaBin]++;
      avgA[ptBin][etaBin] += 0.5*(pt_F[pJet]-pt_F[rJet])/avgPt;
      h_nEntriesA[ptBin][etaBin]->Fill(1);
      h_avgA[ptBin][etaBin]->Fill(0.5*(pt_F[pJet]-pt_F[rJet])/avgPt);
      avgAHisto[ptBin][etaBin]->Fill(0.5*(pt_F[pJet]-pt_F[rJet])/avgPt);
      // end avgA histo stuff //


      // avgB histo stuff //
      // histB stuff, "MPF Method", uses missEt
      TLorentzVector missEt = findMissEt(nPFpart_, pfId_, pfPt_, pfEta_, pfPhi_, 
					 nref, pt_F, eta_F, phi_F, 		  
					 m_F, rawpt_F, eSum, phoSum);           
      for(int i=0; i<nref; i++){//second jet loop
	int etaBin = findBin(eta_F[pJet],nbins_eta,xbins_eta);
	int ptBin = findBin(avgPt,nbins_pt,xbins_pt);
	TLorentzVector jetVec(pt_F[rJet],eta_F[rJet],phi_F[rJet],m_F[rJet]);
	
	avgB[ptBin][etaBin] += (missEt.Dot(jetVec) / ( 2*avgPt * jetVec.Mag()));
	h_avgB[ptBin][etaBin]->Fill((missEt.Dot(jetVec) / ( 2*avgPt * jetVec.Mag())));
	avgBHisto[ptBin][etaBin]->Fill((missEt.Dot(jetVec) / ( 2*avgPt * jetVec.Mag())));
	nEntriesB[ptBin][etaBin]++;
	h_nEntriesB[ptBin][etaBin]->Fill(1);
      }//end second jet loop
      // end avgB histo stuff //
      // end avgA and avgB histo stuff //


    }//exit jet loop
    if(debugMode)std::cout<<std::endl;
    fin->Close();
  }//end file loop

  // open output file
  std::cout<< "opening output file "<<outfile<<std::endl;
  TFile *fout = new TFile( outfile.c_str() ,"recreate");
  fout->cd();

  // write histos in serarate pt/eta bins
  std::cout<< "writing histograms..."<<std::endl;
  for(int i=0; i<nbins_pt; i++){
    for(int j=0; j<nbins_eta; j++){
      // write the avgA stuff
      avgAHisto[i][j]->Write();
      h_avgA[i][j]->Write();
      h_nEntriesA[i][j]->Write();
      // write the avgB stuff
      avgBHisto[i][j]->Write();
      h_avgB[i][j]->Write();
      h_nEntriesB[i][j]->Write();
      // write the avgAbsPho stuff
      avgAbsPhoHisto[i][j]->Write();
      h_avgAbsPho[i][j]->Write();
      h_nEntriesAbsPho[i][j]->Write();
    }
  }
  std::cout<< "closing file "<< outfile <<std::endl;
  fout->Close();
  return 0;
} //end deriveResponse


//// sumResponse
// ---------------------------------------------------------------------------------------------------------------
// reads output file with "Response Histos" in each pt,eta bin
// derives total response in each pt bin by summing the histos over eta bins
// doDraw option allows one to print the histos directly from the routine
int sumResponse(std::string filename, std::string outFileName,
                bool isMC, bool doDraw, bool debugMode){
  
  // debug mode
  const bool loopDebugMode = (debugMode)&&false;
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl;
  if(loopDebugMode)std::cout<<"loopDebugMode is ON"<<std::endl;
  std::cout<<std::endl<<"///////////////////"<<std::endl;
  std::cout<<"reading file "<<filename<<std::endl;
  if(loopDebugMode)std::cout<<"nbins_pt is  "<<nbins_pt  <<std::endl;
  if(loopDebugMode)std::cout<<"nbins_eta is "<<nbins_eta <<std::endl;
  if(debugMode)std::cout<<    "isMC is      "<<isMC      <<std::endl;
  std::cout<<"///////////////////"<<std::endl<<std::endl;  
  
  // open the input file, and an output file to write histos to
  TFile *fin = new TFile(filename.c_str());
  std::cout<< "opening output file "<<outFileName << std::endl;
  TFile *fout = new TFile(outFileName.c_str(),"RECREATE");
  fout->cd();

  // sum deriveResponse histos across eta bins in each pt bin, yield response
  if(debugMode)std::cout<<"summing response histos across eta bins in each pt bin..." << std::endl;
  for(int i=0; i<nbins_pt; i++){//pt loop

    //avgA and hRel
    int totEntriesA_ptbin_i=0;
    std::string  hRelResponseTitle_i="hRelResponse_pt"+std::to_string(i);
    TH1F* hRelResponse_i = new TH1F(hRelResponseTitle_i.c_str(),"",nbins_eta,xbins_eta);
    //avgB and hMPF
    int totEntriesB_ptbin_i=0;
    std::string  hMPFResponseTitle_i="hMPFResponse_pt"+std::to_string(i);
    TH1F* hMPFResponse_i = new TH1F(hMPFResponseTitle_i.c_str(),"",nbins_eta,xbins_eta);
    //avgAbsPho and hMPFAbsPho
    int totEntriesAbsPho_ptbin_i=0;
    std::string  hMPFAbsPhoResponseTitle_i="hMPFAbsPhoResponse_pt"+std::to_string(i);
    TH1F* hMPFAbsPhoResponse_i = new TH1F(hMPFAbsPhoResponseTitle_i.c_str(),"",nbins_eta,xbins_eta);

    // for i'th repsponse histo, set content+error of bin j+1
    for(int j=0; j<nbins_eta; j++){//eta loop

      //// avgA and hRelResponse ////
      //avgA[i][j] = avgA[i][j]/(double)nEntries[i][j];
      std::string avgAHistTitle_ij="avgAHisto_pt"+std::to_string(i)+"_eta"+std::to_string(j);
      TH1F* avgAHisto_ij = (TH1F*)fin->Get( avgAHistTitle_ij.c_str() );
      float avgAMean_ij=avgAHisto_ij->GetMean();
      int avgAEntries_ij=avgAHisto_ij->GetEntries(); 
      if(loopDebugMode)std::cout<< "avgAHistTitle    = " << avgAHistTitle_ij<< std::endl;
      if(loopDebugMode)std::cout<< "(avgAMean_ij, avgAEntries_ij) = ("<<avgAMean_ij<<", "<<avgAEntries_ij<<")"<<std::endl;
      float hRelResponse_i_jadd1 = ( (1+avgAMean_ij) / (1-avgAMean_ij) );
      float hRelResponse_i_jadd1_err = ( hRelResponse_i_jadd1*(1./sqrt(avgAEntries_ij)) );
      if(loopDebugMode)std::cout<< "hRelResponse_i_jadd1 = " << hRelResponse_i_jadd1 << std::endl;
      if(loopDebugMode)std::cout<< "hRelResponse_i_jadd1_err = " << hRelResponse_i_jadd1_err << std::endl;
      // protect against NaN
      if(avgAEntries_ij!=0) hRelResponse_i->SetBinContent( j+1 , hRelResponse_i_jadd1 );     
      else hRelResponse_i->SetBinContent( j+1 , 0 );
      if(avgAEntries_ij!=0) hRelResponse_i->SetBinError( j+1 , hRelResponse_i_jadd1_err );
      else hRelResponse_i->SetBinError( j+1 , 0 );
      totEntriesA_ptbin_i+=avgAEntries_ij;    

      //// avgB and hMPFResponse ////
      //avgB[i][j] = avgB[i][j]/(double)nEntries[i][j];
      std::string avgBHistTitle_ij="avgBHisto_pt"+std::to_string(i)+"_eta"+std::to_string(j);
      TH1F* avgBHisto_ij = (TH1F*)fin->Get( avgBHistTitle_ij.c_str() );
      float avgBMean_ij=avgBHisto_ij->GetMean();
      int avgBEntries_ij=avgBHisto_ij->GetEntries(); 
      if(loopDebugMode)std::cout<< "avgBHistTitle    = " << avgBHistTitle_ij<< std::endl;
      if(loopDebugMode)std::cout<< "(avgBMean_ij, avgBEntries_ij) = ("<<avgBMean_ij<<", "<<avgBEntries_ij<<")"<<std::endl;
      float hMPFResponse_i_jadd1 = ( (1+avgBMean_ij) / (1-avgBMean_ij) );
      float hMPFResponse_i_jadd1_err = ( hMPFResponse_i_jadd1*(1./sqrt(avgBEntries_ij)) );
      if(loopDebugMode)std::cout<< "hMPFResponse_i_jadd1 = " << hMPFResponse_i_jadd1 << std::endl;
      if(loopDebugMode)std::cout<< "hMPFResponse_i_jadd1_err = " << hMPFResponse_i_jadd1_err << std::endl;
      // protect against NaN
      if(avgBEntries_ij!=0) hMPFResponse_i->SetBinContent( j+1 , hMPFResponse_i_jadd1 );     
      else hMPFResponse_i->SetBinContent( j+1 , 0 );
      if(avgBEntries_ij!=0) hMPFResponse_i->SetBinError( j+1 , hMPFResponse_i_jadd1_err );
      else hMPFResponse_i->SetBinError( j+1 , 0 );
      totEntriesB_ptbin_i+=avgBEntries_ij;    

      //// avgAbsPho and hMPFAbsPhoResponse ////
      //avgAbsPho[i][j] = avgAbsPho[i][j]/(double)nEntries[i][j];
      std::string avgAbsPhoHistTitle_ij="avgAbsPhoHisto_pt"+std::to_string(i)+"_eta"+std::to_string(j);
      TH1F* avgAbsPhoHisto_ij = (TH1F*)fin->Get( avgAbsPhoHistTitle_ij.c_str() );
      float avgAbsPhoMean_ij=avgAbsPhoHisto_ij->GetMean();
      int avgAbsPhoEntries_ij=avgAbsPhoHisto_ij->GetEntries(); 
      if(loopDebugMode)std::cout<< "avgAbsPhoHistTitle    = " << avgAbsPhoHistTitle_ij<< std::endl;
      if(loopDebugMode)std::cout<< "(avgAbsPhoMean_ij, avgAbsPhoEntries_ij) = ("<<avgAbsPhoMean_ij<<", "<<avgAbsPhoEntries_ij<<")"<<std::endl;
      //      float hMPFAbsPhoResponse_i_jadd1 = ( (1+avgAbsPhoMean_ij) / (1-avgAbsPhoMean_ij) );
      float hMPFAbsPhoResponse_i_jadd1 = ( avgAbsPhoMean_ij );
      float hMPFAbsPhoResponse_i_jadd1_err = ( hMPFAbsPhoResponse_i_jadd1*(1./sqrt(avgAbsPhoEntries_ij)) );
      if(loopDebugMode)std::cout<< "hMPFAbsPhoResponse_i_jadd1 = " << hMPFAbsPhoResponse_i_jadd1 << std::endl;
      if(loopDebugMode)std::cout<< "hMPFAbsPhoResponse_i_jadd1_err = " << hMPFAbsPhoResponse_i_jadd1_err << std::endl;
      // protect against NaN
      if(avgAbsPhoEntries_ij!=0) hMPFAbsPhoResponse_i->SetBinContent( j+1 , hMPFAbsPhoResponse_i_jadd1 );     
      else hMPFAbsPhoResponse_i->SetBinContent( j+1 , 0 );
      if(avgAbsPhoEntries_ij!=0) hMPFAbsPhoResponse_i->SetBinError( j+1 , hMPFAbsPhoResponse_i_jadd1_err );
      else hMPFAbsPhoResponse_i->SetBinError( j+1 , 0 );
      totEntriesAbsPho_ptbin_i+=avgAbsPhoEntries_ij;    
    }//end eta loop

    // hRelResponse
    if(loopDebugMode)std::cout<< "setting entries of hRelResponse_i to totEntriesA_ptbin_i = " << totEntriesA_ptbin_i << std::endl;
    hRelResponse_i->SetEntries(totEntriesA_ptbin_i);
    if(debugMode)std::cout<< "writing hRelResponse_i..." << std::endl;
    hRelResponse_i->SetMarkerColor(colors[i]);
    hRelResponse_i->SetLineColor(colors[i]);
    std::string hRelResponse_i_Title="Rel, "+std::to_string((int)xbins_pt[i])+" < pt < "+std::to_string((int)xbins_pt[i+1])+" GeV";
    hRelResponse_i->SetTitle(hRelResponse_i_Title.c_str());    
    hRelResponse_i->Write();

    // hMPFResponse
    if(loopDebugMode)std::cout<< "setting entries of hMPFResponse_i to totEntriesB_ptbin_i = " << totEntriesB_ptbin_i << std::endl;
    hMPFResponse_i->SetEntries(totEntriesB_ptbin_i);
    if(debugMode)std::cout<< "writing hMPFResponse_i..." << std::endl;
    hMPFResponse_i->SetMarkerColor(colors[i]);
    hMPFResponse_i->SetLineColor(colors[i]);
    std::string hMPFResponse_i_Title="MPF, "+std::to_string((int)xbins_pt[i])+" < pt < "+std::to_string((int)xbins_pt[i+1])+" GeV";
    hMPFResponse_i->SetTitle(hMPFResponse_i_Title.c_str());    
    hMPFResponse_i->Write();

    // hMPFAbsPhoResponse
    if(loopDebugMode)std::cout<< "setting entries of hMPFAbsPhoResponse_i to totEntriesAbsPho_ptbin_i = " << totEntriesAbsPho_ptbin_i << std::endl;
    hMPFAbsPhoResponse_i->SetEntries(totEntriesAbsPho_ptbin_i);
    if(debugMode)std::cout<< "writing hMPFAbsPhoResponse_i..." << std::endl;
    hMPFAbsPhoResponse_i->SetMarkerColor(colors[i]);
    hMPFAbsPhoResponse_i->SetLineColor(colors[i]);
    std::string hMPFAbsPhoResponse_i_Title="MPFAbsPho, "+std::to_string((int)xbins_pt[i])+" < pt < "+std::to_string((int)xbins_pt[i+1])+" GeV";
    hMPFAbsPhoResponse_i->SetTitle(hMPFAbsPhoResponse_i_Title.c_str());    
    hMPFAbsPhoResponse_i->Write();
  }//end pt loop

  if(doDraw){
    if(debugMode)std::cout<<"drawing all responses.."<<std::endl;
    //// avgA and hRelResponse ////
    TCanvas *c0 = new TCanvas("hRelResponse_all","Rel, all avail. pt",600,600);
    TLegend *l0 = new TLegend(0.4,0.8,0.7,0.9);
    c0->cd(); 
    std::cout<<"drawing hRelResponse_all..."<<std::endl;
    for(int i=0;i<nbins_pt;++i){ 
      std::string  hRelResponseTitle_i="hRelResponse_pt"+std::to_string(i);
      if(loopDebugMode)std::cout<<"drawing "<<hRelResponseTitle_i<<std::endl;
      TH1F* hRelResponse_i = (TH1F*)(fout->Get(hRelResponseTitle_i.c_str())->Clone()) ;
      if(i==0){//canvas prefers the first histrogram draw for things like title, range, etc.
	hRelResponse_i->SetTitle("Rel, all pt Bins");
	TAxis* theYAxis=hRelResponse_i->GetYaxis();
	theYAxis->SetRangeUser(0.,2.75);
      }
      std::string type = "";
      if(!isMC) type = "Data"; 
      else type = "PYTHIA_CUETP8M1";      
      hRelResponse_i->Draw("same");
      if(isMC) l0->AddEntry("", type.c_str() ,"");
      std::string legTitle= std::to_string((int)xbins_pt[i])+" < pt < "+std::to_string((int)xbins_pt[i+1]);
      l0->AddEntry(hRelResponse_i,legTitle.c_str());
    }
    l0->Draw("same");

    //// avgB and hMPFResponse ////    
    TCanvas *c1 = new TCanvas("hMPFResponse_all","MPF, all avail. pt",600,600);
    TLegend *l1 = new TLegend(0.4,0.8,0.7,0.9);
    c1->cd(); 
    std::cout<<"drawing hMPFResponse_all..."<<std::endl;
    for(int i=0;i<nbins_pt;++i){ 
      std::string  hMPFResponseTitle_i="hMPFResponse_pt"+std::to_string(i);
      if(loopDebugMode)std::cout<<"drawing "<<hMPFResponseTitle_i<<std::endl;
      TH1F* hMPFResponse_i = (TH1F*)(fout->Get(hMPFResponseTitle_i.c_str())->Clone()) ;
      if(i==0){//canvas prefers the first histrogram draw for things like title, range, etc.
	hMPFResponse_i->SetTitle("MPF, all pt Bins");
	TAxis* theYAxis=hMPFResponse_i->GetYaxis();
	theYAxis->SetRangeUser(0.,2.75);
      }
      std::string type = "";
      if(!isMC) type = "Data"; 
      else type = "PYTHIA_CUETP8M1";      
      hMPFResponse_i->Draw("same");
      if(isMC) l1->AddEntry("", type.c_str(), "");
      std::string legTitle= std::to_string((int)xbins_pt[i])+" < pt < "+std::to_string((int)xbins_pt[i+1]);
      l1->AddEntry(hMPFResponse_i,legTitle.c_str());
    }
    l1->Draw("same");

    //// avgAbsPho and hMPFAbsPhoResponse ////    
    TCanvas *c2 = new TCanvas("hMPFAbsPhoResponse_all","MPFAbsPho, all avail. pt",600,600);
    TLegend *l2 = new TLegend(0.4,0.8,0.7,0.9);
    c2->cd(); 
    std::cout<<"drawing hMPFAbsPhoResponse_all..."<<std::endl;
    for(int i=0;i<nbins_pt;++i){ 
      std::string  hMPFAbsPhoResponseTitle_i="hMPFAbsPhoResponse_pt"+std::to_string(i);
      if(loopDebugMode)std::cout<<"drawing "<<hMPFAbsPhoResponseTitle_i<<std::endl;
      TH1F* hMPFAbsPhoResponse_i = (TH1F*)(fout->Get(hMPFAbsPhoResponseTitle_i.c_str())->Clone()) ;
      if(i==0){//canvas prefers the first histrogram draw for things like title, range, etc.
	hMPFAbsPhoResponse_i->SetTitle("MPFAbsPho, all pt Bins");
	TAxis* theYAxis=hMPFAbsPhoResponse_i->GetYaxis();
	theYAxis->SetRangeUser(0.,2.75);
      }
      std::string type = "";
      if(!isMC) type = "Data"; 
      else type = "PYTHIA_CUETP8M1";      
      hMPFAbsPhoResponse_i->Draw("same");
      if(isMC) l2->AddEntry("", type.c_str(), "");
      std::string legTitle= std::to_string((int)xbins_pt[i])+" < pt < "+std::to_string((int)xbins_pt[i+1]);
      l2->AddEntry(hMPFAbsPhoResponse_i,legTitle.c_str());
    }
    l2->Draw("same");

    // write the drawn canvases to the output file
    fout->cd();
    if(debugMode)std::cout<<"saving c0.."<<std::endl;
    c0->Write();
    if(debugMode)std::cout<<"saving c1.."<<std::endl;
    c1->Write();
    if(debugMode)std::cout<<"saving c2.."<<std::endl;
    c2->Write();
  }//end doDraw section

  std::cout<< "output file "<<outFileName<<" written"<<std::endl;  
  std::cout<< "closing..."<<std::endl;  
  fout->Close();
  return 0;  
} //end sumResponse


//// main
// ---------------------------------------------------------------------------------------------------------------
// front-end/steering/control for .exe file
int main(int argc, char *argv[]){

  // rStatus=-1 for bad input, 1 for not finishing job of choice
  int rStatus = -1;
  if( argc==1 ) {//no input arguments, error
    std::cout<<"need at least 1 argument!"<<std::endl;
    std::cout<<"rStatus="<<rStatus<<std::endl<<std::endl;
    return rStatus;
  } 
  
  bool doDerive=(atoi(argv[1])==0);
  bool doSum=(atoi(argv[1])==1);    
  if(!doDerive && !doSum ){//bad first argument, error
    std::cout<<"bad value for <function>!"<<std::endl;
    std::cout<<"use 0 for deriveResponse"<<std::endl;
    std::cout<<"use 1 for sumResponse"<<std::endl;
    std::cout<<"rStatus="<<rStatus<<std::endl<<std::endl;
    return rStatus;
  }  
  
  // expected usage
  if( argc==2 ){//default case
    rStatus=1;
    if(doDerive) rStatus=deriveResponse();
    if(doSum) rStatus=sumResponse();
    std::cout<<"rStatus="<<rStatus<<std::endl<<std::endl;
    return rStatus;
  }
  else if( argc==deriveRespArgCount && doDerive ){//deriveResponse case
    rStatus=1;
    int startfile=atoi(argv[2]); int endfile=atoi(argv[3]); 
    std::string inputFileList=argv[4]; std::string outputFileName=argv[5]; 
    int jetRadius=atoi(argv[6]); bool isMC=(bool)atoi(argv[7]); bool debug=(bool)atoi(argv[8]);      
    rStatus = deriveResponse( startfile, endfile, 
				   inputFileList, outputFileName,
				   jetRadius, isMC, debug         );    
    std::cout<<"rStatus="<<rStatus<<std::endl<<std::endl;
    return rStatus;
  }
  else if( argc==sumRespArgCount && doSum ){//sumResponse case
    rStatus=1;
    std::string inputFileName=argv[2]; std::string outFileName=argv[3];
    bool isMC=(bool)atoi(argv[4]); bool doDraw=(bool)atoi(argv[5]);  bool debugMode=(bool)atoi(argv[6]);
    rStatus = sumResponse( inputFileName, outFileName, 
				isMC, doDraw, debugMode     );    
    std::cout<<"rStatus="<<rStatus<<std::endl<<std::endl;
    return rStatus;
  }
  else {//usage is incorrect, instruct of proper use
    std::cout<<"not enough arguments ";//<<std::endl;
    std::cout<<"to use default arguments as inputs, do ./DijetResponse.exe <function> "<<std::endl<<std::endl;    
    if (doDerive){
      std::cout<<"for deriveResponse, do..."<<std::endl;
      std::cout<<"./DijetResponse.exe 0 ";
      std::cout<<"<startfile> <endfile> <inputFileList> <outputFileName> <jetRadius> <isMC> <debug>"<<std::endl;
    }
    if (doSum){
      std::cout<<"for sumResponse, do..."<<std::endl;
      std::cout<<"./DijetResponse.exe 1 ";
      std::cout<<"<inputFile> <outputFile> <isMC> <doDraw> <debugMode>"<<std::endl;
    }
    std::cout<<"rStatus="<<rStatus<<std::endl<<std::endl;
    return rStatus;
  }
}// end main
