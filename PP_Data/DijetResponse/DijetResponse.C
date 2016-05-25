// Originally written by: Kurt Jung, December 2015
// Overhaul, Ian Laflotte, May 1st 2016
// Macro for deriving the dijet relative and absolute response from pp data/MC forests
// Uses JME-13-004 as a reference - "A" and "B" formulae are defined there
//first look for photon balance (AN2013-179 / JME-13-004)

///////////////////////////////////////////
////////// includes/preprocessor //////////
///////////////////////////////////////////

// C includes			    
#include <cstdlib>
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <vector>
// ROOT includes
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TH1F.h>
// custom
#include "DijetResponse.h"

// deriveDijetResponse
// default inputs
/////
const int defStartFile=0;
const int defEndFile=1; //exclusive boundary, range of files run over doesnt include endfile
const std::string defInFilelist = "../filelists/5p02TeV_HighPtJet80_2Files_debug_forests.txt";
/////
const int defRadius=4;
const bool defIsMCderive=false;
/////
const std::string defDeriveOutputName = "deriveResponseOut.root"; 
const bool defDebugMode = true; 
/////
//int iii=0;
//std::cout<<"here"<<iii<<std::endl;iii++;

// derive dijet response is making empty avgA histos
// the macro
int deriveDijetResponse(int startfile = defStartFile, int endfile = defEndFile, std::string infile_Forest=defInFilelist,
			int radius = 4, bool isMC = defIsMCderive, 
			std::string outfile = defDeriveOutputName, bool debugMode = defDebugMode ){

  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl;
  const bool fastDebugMode = (debugMode)&&false; //if debugMode is off, fastDebugMode shouldn't be on
  if(fastDebugMode)std::cout<<"fastDebugMode is ON"<<std::endl;

  std::cout<<std::endl<<"///////////////////"<<std::endl;
  std::cout<<"reading files "<< startfile << " to " << endfile<<std::endl;
  if(debugMode)std::cout<<"radius =  " << radius<<std::endl;
  if(debugMode)std::cout<<"isMC   = " << isMC <<std::endl;
  std::cout<<"///////////////////"<<std::endl<<std::endl;  

  TH1F *avgAHisto[nbins_pt][nbins_eta];
  //TH1F *avgBHisto[nbins_pt][nbins_eta];
  //TH1F *hAvgAbsPhoResponse[nbins_pt][nbins_eta];
  //TH1F *hRelResponse[nbins_pt];
  //TH1F *hMPFResponse[nbins_pt];  
  //TH1F *hAbsPhoResponse[nbins_pt];
  //TH1F *hMPFAbsPhoResponse[nbins_pt];
  
  for(int i=0; i<nbins_pt; i++){//pt loop
    //hRelResponse[i] = new TH1F(Form("hRelResponse_pt%d",i),"",nbins_eta,xbins_eta);
    //hMPFResponse[i] = new TH1F(Form("hMPFResponse_pt%d",i),"",nbins_eta,xbins_eta);
    //hAbsPhoResponse[i] = new TH1F(Form("hAbsPhoResponse_pt%d",i),"",nbins_eta,xbins_eta);
    //hMPFAbsPhoResponse[i] = new TH1F(Form("hMPFAbsPhoResponse_pt%d",i),"",nbins_eta,xbins_eta);
    for(int j=0; j<nbins_eta; j++){//eta loop
      avgAHisto[i][j] = new TH1F(Form("avgAHisto_pt%d_eta%d",i,j),"",50,-1,1);
      //avgBHisto[i][j] = new TH1F(Form("avgBHisto_pt%d_eta%d",i,j),"",50,-1,1);
      //hAvgAbsPhoResponse[i][j] = new TH1F(Form("hAvgAbsPhoResponse_pt%d_eta%d",i,j),"",50,-1,1);
    }//end eta loop
  }//end pt loop


  int nEntriesA[nbins_pt][nbins_eta];
  double avgA[nbins_pt][nbins_eta];
  TH1D *h_avgA[nbins_pt][nbins_eta];
  TH1D *h_nEntriesA[nbins_pt][nbins_eta];
  
  //int nEntriesB[nbins_pt][nbins_eta];  
  //double avgB[nbins_pt][nbins_eta];
  //TH1D *h_avgB[nbins_pt][nbins_eta];
  //TH1D *h_nEntriesB[nbins_pt][nbins_eta];

  //int nEntriesAbs[nbins_pt][nbins_eta];
  //double avgAbsPhoResponse[nbins_pt][nbins_eta];
  //TH1D *h_avgAbsPhoResponse[nbins_pt][nbins_eta];
  //TH1D *h_nEntriesAbs[nbins_pt][nbins_eta];
  
  for(int i=0; i<nbins_pt; i++){//pt loop
    for(int j=0; j<nbins_eta; j++){//eta loop

      avgA[i][j]=0;      nEntriesA[i][j]=0;
      h_avgA[i][j]=new TH1D(Form("h_avgA_%d_%d",i,j),"",2, 0, 2);;
      h_nEntriesA[i][j]=new TH1D(Form("h_nEntriesA_%d_%d",i,j),"",2, 0, 2);;

      //avgB[i][j]=0;      nEntriesB[i][j]=0;
      //h_avgB[i][j]=new TH1D(Form("h_avgB_%d_%d",i,j),"",2, 0, 2);;
      //h_nEntriesB[i][j]=new TH1D(Form("h_nEntriesB_%d_%d",i,j),"",2, 0, 2);;

      //avgAbsPhoResponse[i][j]=0;      nEntriesAbs[i][j]=0;
      //h_avgAbsPhoResponse[i][j]=new TH1D(Form("h_avgAbsPhoResponse_%d_%d",i,j),"",2, 0, 2);;
      //h_nEntriesAbs[i][j]=new TH1D(Form("h_nEntriesAbs_%d_%d",i,j),"",2, 0, 2);;      

    }//end etaloop
  }//end pt loop
  
  // open filelist, loop to the starting file
  if(debugMode)std::cout<<"reading filelist "<<infile_Forest<<std::endl;
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;  
  for(int ifile = 0;ifile<startfile;++ifile) instr_Forest >> filename_Forest;
  
  //loop over specified files 
  for(int ifile = startfile; ifile<endfile; ++ifile){//file loop
    
    // grab file
    instr_Forest>>filename_Forest;
    if(atMIT)filename_Forest=hadoopDir+filename_Forest;
    else filename_Forest=xrootdDirect+filename_Forest;

    if(debugMode)std::cout<<"opening trees from file "<<filename_Forest<<std::endl;    
    TFile *fin = TFile::Open(filename_Forest.c_str());
    TTree *jtTree  = (TTree*)fin->Get( Form("ak%dPFJetAnalyzer/t",radius) );
    TTree *pfTree  = (TTree*)fin->Get( "pfcandAnalyzer/pfTree"            );
    TTree *phoTree = (TTree*)fin->Get( "ggHiNtuplizer/EventTree"          );

    int nref; 
    float pt_F[1000],  eta_F[1000],  phi_F[1000];
    float rawpt_F[1000], m_F[1000], eSum[1000], phoSum[1000];

    int nPFpart_;
    //int pfId[10000];
    //float pfPt[10000], pfEta[10000], pfPhi[10000];
    std::vector<int> *pfId_;
    std::vector<float> *pfPt_=0, *pfEta_=0, *pfPhi_=0;
    
    int nPho;
    std::vector<float> *phoEt=0, *phoEta=0, *phoPhi=0; 

    jtTree->SetBranchAddress("jtpt",pt_F);
    jtTree->SetBranchAddress("jteta",eta_F);
    jtTree->SetBranchAddress("jtphi",phi_F);
    jtTree->SetBranchAddress("jtm",m_F);
    jtTree->SetBranchAddress("rawpt",rawpt_F);
    jtTree->SetBranchAddress("nref",&nref);
    jtTree->SetBranchAddress("eSum",&eSum);
    jtTree->SetBranchAddress("photonSum",&phoSum);

    pfTree->SetBranchAddress("nPFpart",&nPFpart_);
    pfTree->SetBranchAddress("pfId",&pfId_);
    pfTree->SetBranchAddress("pfPt",&pfPt_);//added ampersand to switch from float to vector<float>*
    pfTree->SetBranchAddress("pfEta",&pfEta_);
    pfTree->SetBranchAddress("pfPhi",&pfPhi_);

    phoTree->SetBranchAddress("nPho",&nPho);
    phoTree->SetBranchAddress("phoEt",&phoEt);
    phoTree->SetBranchAddress("phoEta",&phoEta);
    phoTree->SetBranchAddress("phoPhi",&phoPhi);

    // loop over jets
    int jetEntries=jtTree->GetEntries();
    if(fastDebugMode)jetEntries=100;
    if(debugMode)std::cout<<"jetEntries = "<< jetEntries <<std::endl;

    std::cout<<"beginning jet loop..."<< std::endl;    
    for(int ientry=0; ientry<jetEntries; ientry++){//jet loop
      //if(debugMode)std::cout<<"for [ptBin][etaBin]=]["<<ptBin<<"]["<<etaBin<<"]"<<std::endl;
      
      jtTree->GetEntry(ientry);
      pfTree->GetEntry(ientry);
      //if(debugMode)std::cout<<"successfully grabbed entry from pftree"<<std::endl;
      phoTree->GetEntry(ientry);

      if(nPFpart_>10000) std::cout<<"warning! nPF: "<<nPFpart_<<std::endl;

      //// avgPhoResponse histo stuff
      //// find response of detector to jet with particle-flow photon
      //for(int i=0; i<nPho; i++){
      //	if(phoEt->at(i)>30 && abs(phoEta->at(i))<1.4442){
      //	  // the jet that the photon actually corresponds to?
      //	  // why only loop up to nref-1?
      //	  for(int ijet=0; ijet<nref-1; ijet++){
      //	    if( pt_F[ijet]>=minJetPt && //jetpt cut
      //		abs(phoPhi->at(i)-phi_F[ijet])>2.95 && //angular separation between the +jet
      //		pt_F[ijet+1]/phoEt->at(i)<0.2 ){//jetpt/phoEt ratio
      //
      //	      if(debugMode)std::cout << "jet photon pair found" << std::endl;	    
      //	      TLorentzVector missEt = findMissEt(nPFpart_, pfId_, pfPt_, 
      //						 pfEta_, pfPhi_, 
      //						 eSum, phoSum, 
      //						 nref, pt_F, rawpt_F, 
      //						 eta_F, phi_F, m_F);
      //	      int etaBin = findBin(phoEta->at(i),nbins_eta,xbins_eta);
      //	      int ptBin = findBin(phoEt->at(i),nbins_pt,xbins_pt);
      //	      TLorentzVector phoVec(phoEt->at(i),phoEta->at(i),phoPhi->at(i),0.);
      //	      double num = missEt.Dot(phoVec);
      //	      
      //	      nEntriesAbs[ptBin][etaBin]++;
      //	      avgAbsPhoResponse[ptBin][etaBin] += (1+(num/phoVec.Dot(phoVec)));
      //	      if(debugMode)std::cout<<"avgAbsPhoResponse = "<< avgAbsPhoResponse[ptBin][etaBin]<<std::endl;
      //	      if(debugMode)std::cout<<"for [ptBin][etaBin]=]["<<ptBin<<"]["<<etaBin<<"]"<<std::endl;
      //	      h_nEntriesAbs[ptBin][etaBin]->Fill(1);
      //	      h_avgAbsPhoResponse[ptBin][etaBin]->Fill((1+(num/phoVec.Dot(phoVec))));
      //	      hAvgAbsPhoResponse[ptBin][etaBin]->Fill(1+(num/phoVec.Dot(phoVec)));
      //
      //	    }//end photon+jetcut
      //	  }//end jet loop
      //	}//end photon cut
      //}//end photon loop
      //// end avgPhoResponse histo stuff

      // standard Rrel method, need at least two jets
      // avgA and avgB histos
      if(nref>1){

	// assign two highest pt jets as rJet or pJet using eta/random# critera
	int rJet, pJet;
	if( abs(eta_F[0])<1.3 && abs(eta_F[1])<1.3 ){ 
	  int randInt=rand();// pseudorandom at best, fix with srand() in future, results could be slightly different
	  //if(debugMode)std::cout<< "rand() gave "<<randInt <<std::endl;
	  pJet = randInt%2;  rJet = (randInt+1)%2;
	}
	else if( abs(eta_F[0])<1.3 && abs(eta_F[1])<5. ){
	  pJet = 1; rJet = 0;
	}
	else if( abs(eta_F[1])<1.3 && abs(eta_F[0])<5. ){
	  pJet = 0; rJet = 1;
	}
	else continue; //does this skip to the end of the scope of the if statement, or the parent for statement? why do i find this confusing now?
	//if(debugMode)std::cout<< "rJet is "<<rJet<<" and pJet is "<<pJet <<std::endl;
	
	// minjetpt cut
	const double jtptMin = 40;
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

	// histA stuff, DIJET ASYMMETRY DIST
	nEntriesA[ptBin][etaBin]++;
	avgA[ptBin][etaBin] += 0.5*(pt_F[pJet]-pt_F[rJet])/avgPt;
	h_nEntriesA[ptBin][etaBin]->Fill(1);
	h_avgA[ptBin][etaBin]->Fill(0.5*(pt_F[pJet]-pt_F[rJet])/avgPt);
	avgAHisto[ptBin][etaBin]->Fill(0.5*(pt_F[pJet]-pt_F[rJet])/avgPt);

	// histB stuff, "MPF Method", jet loop and missEt
	TLorentzVector missEt = findMissEt(nPFpart_, pfId_, pfPt_, pfEta_, pfPhi_, eSum, phoSum,
					   nref, pt_F, rawpt_F, eta_F, phi_F, m_F);
	for(int i=0; i<nref; i++){//jet loop
	  int etaBin = findBin(eta_F[pJet],nbins_eta,xbins_eta);
	  int ptBin = findBin(avgPt,nbins_pt,xbins_pt);

	  TLorentzVector jetVec(pt_F[rJet],eta_F[rJet],phi_F[rJet],m_F[rJet]);
	  avgB[ptBin][etaBin] += (missEt.Dot(jetVec) / ( 2*avgPt * jetVec.Mag()));
	  h_avgB[ptBin][etaBin]->Fill((missEt.Dot(jetVec) / ( 2*avgPt * jetVec.Mag())));
	  avgBHisto[ptBin][etaBin]->Fill((missEt.Dot(jetVec) / ( 2*avgPt * jetVec.Mag())));
	  nEntriesB[ptBin][etaBin]++;
	  h_nEntriesB[ptBin][etaBin]->Fill(1);
	}//end jet loop

      }//end 2 or more jets req
    }//exit jet loop
    std::cout<<std::endl;
    fin->Close();
  }//end file loop

  // open output file
  if(debugMode)std::cout<< "opening file "<<outfile<<std::endl;
  TFile *fout = new TFile(Form("%s",outfile.c_str()),"recreate");
  fout->cd();

  // write histos in serarate pt/eta bins
  if(debugMode)std::cout<< "writing histograms..."<<std::endl;
  for(int i=0; i<nbins_pt; i++){
    for(int j=0; j<nbins_eta; j++){
      // check content of A histos
      //write the avgA stuff
      avgAHisto[i][j]->Write();
      h_avgA[i][j]->Write();
      h_nEntriesA[i][j]->Write();

      ////write the avgB stuff
      //avgBHisto[i][j]->Write();
      //h_avgB[i][j]->Write();
      //h_nEntriesB[i][j]->Write();

      ////write the AbsPhoResponse stuff
      //hAvgAbsPhoResponse[i][j]->Write();
      //h_avgAbsPhoResponse[i][j]->Write();
      //h_nEntriesAbs[i][j]->Write();
    }
  }
  if(debugMode)std::cout<< "closing file "<< outfile <<std::endl;
  fout->Close();
  return 0;
} //end deriveDijetResponse


// sumDijetResponse
// default inputs
/////
const std::string defSumInput="testInput_sum_DijetResponse_5p02TeV_HighPtJet80_ak4-0to40.root";
const bool defIsMCsum=false;
const bool defDoDraw=true;
const std::string defSumOutput="sumResponseOut.root";
const bool defDebugModeSum=true;
/////

// the macro
int sumDijetResponse(std::string filename=defSumInput , bool isMC=defIsMCsum , bool doDraw=defDoDraw , std::string outFileName = defSumOutput, bool debugMode=defDebugModeSum){
//int sumDijetResponse(std::string filename=defSumInput , bool isMC=defIsMCsum , std::string outFileName = defSumOutput, bool debugMode=defDebugModeSum){
  
  // debug mode
  const bool loopDebugMode = (debugMode)&&true;
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl;
  if(loopDebugMode)std::cout<<"loopDebugMode is ON"<<std::endl;
  
  if(debugMode)std::cout<<     std::endl<<"///////////////////" <<std::endl;
  if(loopDebugMode)std::cout<< "nbins_pt is  " << nbins_pt      <<std::endl;
  if(loopDebugMode)std::cout<< "nbins_eta is " << nbins_eta     <<std::endl;
  if(debugMode)std::cout<<     "isMC is      " << isMC          <<std::endl;
  if(debugMode)std::cout<<     "///////////////////"<<std::endl <<std::endl;  
  

  // open the file
  std::cout<< "input file is " << filename <<std::endl;
  TFile *fin = new TFile(filename.c_str());

  // open an output file for output plots
  if(debugMode)std::cout<< "opening output file "<<outFileName << std::endl;
  TFile *fout = new TFile(outFileName.c_str(),"RECREATE");
  fout->cd();

  if(debugMode)std::cout<<"setting Response content+error..." << std::endl;

  for(int i=0; i<nbins_pt; i++){//pt loop

    int totEntriesA_ptbin_i=0;
    std::string  hRelResponseTitle_i="hRelResponse_pt"+std::to_string(i);
    TH1F* hRelResponse_i = new TH1F(hRelResponseTitle_i.c_str(),"",nbins_eta,xbins_eta);

    //int totEntriesB_ptbin_i=0;
    //hMPFResponse_i = new TH1F(Form("hMPFResponse_pt%d",i),"",nbins_eta,xbins_eta);

    //int totEntriesAbs_ptbin_i=0;
    //hAbsPhoResponse_i = new TH1F(Form("hAbsPhoResponse_pt%d",i),"",nbins_eta,xbins_eta);
    //hMPFAbsPhoResponse_i = new TH1F(Form("hMPFAbsPhoResponse_pt%d",i),"",nbins_eta,xbins_eta);

    // for i'th repsponse histo, set content+error of bin j+1
    for(int j=0; j<nbins_eta; j++){//eta loop

      ///////////////////////////////
      //// avgA and hRelResponse ////
      ///////////////////////////////

      std::string avgAHistTitle_ij="avgAHisto_pt"+std::to_string(i)+"_eta"+std::to_string(j);
      TH1F* avgAHisto_ij = (TH1F*)fin->Get( avgAHistTitle_ij.c_str() );

      //avgA[i][j] = avgA[i][j]/(double)nEntries[i][j];

      float avgAMean_ij=avgAHisto_ij->GetMean();
      int avgAEntries_ij=avgAHisto_ij->GetEntries(); 

      if(loopDebugMode)std::cout<< "avgAHistTitle    = " << avgAHistTitle_ij<< std::endl;
      if(loopDebugMode)std::cout<< "avgAMean_ij , avgAEntries_ij  =  " << avgAMean_ij <<",   "<< avgAEntries_ij << std::endl;
      
      float hRelResponse_i_jadd1 = ( (1+avgAMean_ij) / (1-avgAMean_ij) );
      float hRelResponse_i_jadd1_err = ( hRelResponse_i_jadd1*(1./sqrt(avgAEntries_ij)) );

      if(loopDebugMode)std::cout<< "hRelResponse_i_jadd1 = " << hRelResponse_i_jadd1 << std::endl;
      if(loopDebugMode)std::cout<< "hRelResponse_i_jadd1_err = " << hRelResponse_i_jadd1_err << std::endl;

      // protect against NaN, all created equal? we'll see...
      if(avgAEntries_ij!=0) hRelResponse_i->SetBinContent( j+1 , hRelResponse_i_jadd1 );     
      else hRelResponse_i->SetBinContent( j+1 , 0 );
      if(avgAEntries_ij!=0) hRelResponse_i->SetBinError( j+1 , hRelResponse_i_jadd1_err );
      else hRelResponse_i->SetBinError( j+1 , 0 );

      totEntriesA_ptbin_i+=avgAEntries_ij;    

      ///////////////////////////////
      //// avgB and hMPFResponse ////
      ///////////////////////////////

      //avgBHisto_i_j = (TH1F*)fin->Get( Form("avgBHisto_pt%d_eta%d",i,j) )->Clone(Form("avgBHisto_pt%d_eta%d",i,j));
      //hMPFResponse_i->SetBinContent( j+1 , (1+avgBHisto_i_j->GetMean())/(1-avgBHisto_i_j->GetMean()) );
      //hMPFResponse_i->SetBinError( j+1 , (hMPFResponse_i->GetBinContent(j+1))*(1./sqrt(avgBHisto_i_j->GetEntries())) );
      //totEntriesB_ptbin_i+=avgBHisto_i_j->GetEntries();
      //     avgB[i][j] = avgB[i][j]/(double)nEntriesB[i][j];      
      //     if(nEntriesB[i][j]) hMPFResponse[i]->SetBinContent(j+1,(1+avgB[i][j])/(1-avgB[i][j]));
      //     //if(nEntriesB[i][j]) hMPFResponse[i]->SetBinError(j+1,hMPFResponse[i]->GetBinContent(j+1)*(1./sqrt(nEntriesB[i][j])));
      //     if(nEntriesB[i][j]) hMPFResponse[i]->SetBinError(j+1,avgBHisto[i][j]->GetRMS()*(1./sqrt(nEntriesB[i][j])));
      //     else hMPFResponse[i]->SetBinContent(j+1,0);
      //     avgB[i][j] = avgB[i][j]/(double)nEntriesB[i][j];
      //     if(nEntriesB[i][j]) hMPFResponse[i]->SetBinContent(j+1,(1+avgB[i][j])/(1-avgB[i][j]));
      //     //if(nEntriesB[i][j]) hMPFResponse[i]->SetBinError(j+1,hMPFResponse[i]->GetBinContent(j+1)*(1./sqrt(nEntriesB[i][j])));
      //     if(nEntriesB[i][j]) hMPFResponse[i]->SetBinError(j+1,avgBHisto[i][j]->GetRMS()*(1./sqrt(nEntriesB[i][j])));
      //     else hMPFResponse[i]->SetBinContent(j+1,0);


      //////////////////////////////////////////////////
      //// avgAbsPhoResponse and hMPFAbsPhoResponse ////
      //////////////////////////////////////////////////

      //hAvgAbsPhoResponse_i_j = (TH1F*)fin->Get(Form("hAvgAbsPhoResponse_pt%d_eta%d",i,j))->Clone(Form("hAvgAbsPhoResponse_pt%d_eta%d",i,j));
      //hMPFAbsPhoResponse_i->SetBinContent( j+1 ,hAvgAbsPhoResponse_i_j->GetMean());
      //hMPFAbsPhoResponse_i->SetBinError( j+1 ,hMPFAbsPhoResponse_i->GetBinContent(j+1)*(1./sqrt(hAvgAbsPhoResponse_i_j->GetEntries())));      
      //totEntriesAbs_ptbin_i+=hAvgAbsPhoResponse_i_j->GetEntries();
      //     if(nEntriesAbs[i][j]) hMPFAbsPhoResponse[i]->SetBinContent(j+1,avgAbsPhoResponse[i][j]/nEntriesAbs[i][j]);
      //     //if(nEntriesAbs[i][j]) hMPFAbsPhoResponse[i]->SetBinError(j+1,hMPFAbsPhoResponse[i]->GetBinContent(j+1)*(1./sqrt(nEntriesAbs[i][j])));
      //     if(nEntriesAbs[i][j]) hMPFAbsPhoResponse[i]->SetBinError(j+1,hAvgAbsPhoResponse[i][j]->GetRMS()*(1./sqrt(nEntriesAbs[i][j])));
      //     else hMPFAbsPhoResponse[i]->SetBinContent(j+1,0);  //     avgA[i][j] = avgA[i][j]/(double)nEntries[i][j];
      //     if(nEntriesAbs[i][j]) hMPFAbsPhoResponse[i]->SetBinContent(j+1,avgAbsPhoResponse[i][j]/nEntriesAbs[i][j]);
      //     //if(nEntriesAbs[i][j]) hMPFAbsPhoResponse[i]->SetBinError(j+1,hMPFAbsPhoResponse[i]->GetBinContent(j+1)*(1./sqrt(nEntriesAbs[i][j])));
      //     if(nEntriesAbs[i][j]) hMPFAbsPhoResponse[i]->SetBinError(j+1,hAvgAbsPhoResponse[i][j]->GetRMS()*(1./sqrt(nEntriesAbs[i][j])));
      //     else hMPFAbsPhoResponse[i]->SetBinContent(j+1,0);

    }//end eta loop

    //set entries in RelResponse histo
    if(loopDebugMode)std::cout<< "setting entries of hRelResponse_i to totEntriesA_ptbin_i = " << totEntriesA_ptbin_i << std::endl;
    hRelResponse_i->SetEntries(totEntriesA_ptbin_i);

    if(debugMode)std::cout<< "writing hRelResponse_i" << std::endl;
    hRelResponse_i->SetMarkerColor(colors[i]);
    hRelResponse_i->SetLineColor(colors[i]);
    hRelResponse_i->SetTitle(Form("Rel, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));    
    hRelResponse_i->Write();

    // hMPFResponse_i->SetMarkerColor(colors[i]);
    // hMPFResponse_i->SetLineColor(colors[i]);
    // hMPFResponse_i->SetMarkerStyle(21);
    // hMPFResponse_i->SetLineStyle(2);
    // hMPFResponse_i->SetTitle(Form("MPF, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    // hMPFResponse_i->Write();

    // hMPFAbsPhoResponse_i->SetMarkerColor(colors[i]);
    // hMPFAbsPhoResponse_i->SetLineColor(colors[i]);
    // hMPFAbsPhoResponse_i->SetTitle(Form("MPF Abs, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    // hMPFAbsPhoResponse_i->Write();
    //
    //  

  }//end pt loop

  if(doDraw){
    std::cout<<"drawing all plots on one figure..."<<std::endl;

    //// avgA and hRelResponse ////
    std::cout<<"drawing hRelResponse_all..."<<std::endl;
    TCanvas *c0 = new TCanvas("hRelResponse_all","Rel, all avail. p_{t}",600,600);
    TLegend *l0 = new TLegend(0.4,0.8,0.7,0.9);
    c0->cd(); 

    for(int i=0;i<nbins_pt;++i){
      std::string  hRelResponseTitle_i="hRelResponse_pt"+std::to_string(i);
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
      if(isMC) l0->AddEntry("",Form("%s",type.c_str()),"");
      l0->AddEntry(hRelResponse_i,Form("%g<p_{T}<%g",xbins_pt[i],xbins_pt[i+1]));
    }
    l0->Draw("same");

    //// avgB and hMPFResponse ////
    //  // if(debugMode)std::cout<<"drawing hMPFResponse Histos..."<<std::endl;
    //  // TCanvas *c1 = new TCanvas("c1","",600,600);
    //  // c1->cd();
    //  // hMPFResponse[0]->Draw();
    //  // TLegend *l1 = new TLegend(0.4,0.8,0.7,0.9);
    //  // l1->AddEntry(hMPFResponse[0],Form("%g<p_{T}<%g",xbins_pt[0],xbins_pt[1]));
    //  // for(int i=1; i<nbins_pt; i++){
    //  //   std::cout<<i<<std::endl;
    //  //   hMPFResponse[i]->Draw("same");
    //  //   if(isMC) l1->AddEntry("",Form("%s",type.c_str()),"");
    //  //   l1->AddEntry(hMPFResponse[i],Form("%g<p_{T}<%g",xbins_pt[i],xbins_pt[i+1]));
    //  // }
    //  // l1->Draw("Same");
    //  // std::cout<<"saving c1"<<std::endl;
    //  // c1->SaveAs(Form("MPFresponse_%s.pdf",type.c_str()),"RECREATE");


    //// avgAbsPhoResponse and hMPFAbsPhoResponse ////    
    //  // if(debugMode)std::cout<<"drawing hMPFAbsPhoResponse Histos..."<<std::endl;
    //  // TCanvas *c2 = new TCanvas("c2","",600,600);
    //  // c2->cd();
    //  // hMPFAbsPhoResponse[0]->Draw();
    //  // TLegend *l2 = new TLegend(0.4,0.8,0.7,0.9);
    //  // l2->AddEntry(hMPFAbsPhoResponse[0],Form("%g<p_{T}<%g",xbins_pt[0],xbins_pt[1]));    
    //  // for(int i=1; i<nbins_pt; i++){
    //  //   std::cout<<i<<std::endl;
    //  //   hMPFAbsPhoResponse[i]->Draw("same");
    //  //   if(isMC) l1->AddEntry("",Form("%s",type.c_str()),"");
    //  //   l2->AddEntry(hMPFAbsPhoResponse[i],Form("%g<p_{T}<%g",xbins_pt[i],xbins_pt[i+1]));
    //  // }    
    //  // l2->Draw("Same");
    //  // std::cout<<"saving c2"<<std::endl;
    //  // c2->SaveAs(Form("MPFABSPhoresponse_%s.pdf",type.c_str()),"RECREATE");
    //}
    
    fout->cd();
    std::cout<<"saving c0.."<<std::endl;
    c0->Write();
  }//end doDraw section
  
  if(debugMode)std::cout<< "output file "<<outFileName<<" written"<<std::endl;  
  if(debugMode)std::cout<< "closing..."<<std::endl;  
  fout->Close();
  
  return 0;
  
} //end sumDijetResponse


//////////////////////////
////////// main //////////
//////////////////////////
// acts as the frontend control for .exe file
int main(int argc, char *argv[]){

  int rStatus = 1;   

  // bad input errors
  if( argc==1 ) {   // no input arguments
    std::cout<<"need at least 1 argument, to test with defaults, do..."<<std::endl;
    std::cout<<"./DijetResponse.exe <function>"<<std::endl<<std::endl;
    return rStatus;
  } 
  
  bool doDerive=(atoi(argv[1])==0);
  bool doSum=(atoi(argv[1])==1);    
  if(!doDerive && !doSum ){  // non-existant function called
    std::cout<<"bad value for <function>"<<std::endl;
    std::cout<<"use 0 for deriveDijetResponse"<<std::endl;
    std::cout<<"use 1 for sumDijetResponse"<<std::endl<<std::endl;
    return rStatus;
  }  
  
  // expected usage cases
  if( argc==2 ){   // expected use for default-value tests
    rStatus=-1;
    if(doDerive) rStatus=deriveDijetResponse();
    if(doSum) rStatus=sumDijetResponse();
    return rStatus;
  }
  else if( argc==deriveRespArgCount && doDerive ){   // expected use for deriveDijetResponse
    rStatus=-1;
    int startfile=atoi(argv[2]); int endfile=atoi(argv[3]); std::string inputFileList=argv[4];
    int jetRadius=atoi(argv[5]); bool isMC=(bool)atoi(argv[6]); 
    std::string outputFileName=argv[7]; bool debug=(bool)atoi(argv[8]);      
    rStatus = deriveDijetResponse( startfile, endfile, inputFileList,
				   jetRadius, isMC,
				   outputFileName, debug);    
    return rStatus;
  }
  else if( argc==sumRespArgCount && doSum ){   // expected use for sumDijetResponse
    rStatus=-1;
    std::string inputFileName=argv[2]; bool isMC=(bool)atoi(argv[3]); bool doDraw=(bool)atoi(argv[4]); std::string outFileName=argv[5]; bool debugMode=(bool)atoi(argv[6]);
    rStatus = sumDijetResponse( inputFileName, isMC, doDraw, outFileName, debugMode);    
    //    std::string inputFileName=argv[2]; bool isMC=(bool)atoi(argv[3]); std::string outFileName=argv[4]; bool debugMode=(bool)atoi(argv[5]);
    //    rStatus = sumDijetResponse( inputFileName, isMC, outFileName, debugMode);    
    return rStatus;
  }
  else {  // inputs don't match any expected use cases
    std::cout<<"not enough arguments ";//<<std::endl;
    std::cout<<"to use default arguments as inputs, do ./DijetResponse.exe <function> "<<std::endl<<std::endl;    
    if (doDerive){
      std::cout<<"for deriveDijetResponse, do..."<<std::endl;
      std::cout<<"./DijetResponse.exe 0 ";
      std::cout<<"<startfile> <endfile> <inputFileList> <jetRadius> <isMC> <outputFileName> <debug>";
    }
    if (doSum){
      std::cout<<"for sumDijetResponse, do..."<<std::endl;
      std::cout<<"./DijetResponse.exe 1 ";
      //      std::cout<<"<inputFile> <isMC> <outputFile> <debugMode>";
      std::cout<<"<inputFile> <isMC> <doDraw> <outputFile> <debugMode>";
    }
    std::cout<<std::endl<<std::endl;
    return rStatus;
  }
}//end main
