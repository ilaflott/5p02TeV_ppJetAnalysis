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
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
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
const std::string defOutputName = "test_output.root"; 
const bool defDebugMode = true; 
/////
//const std::string pp_MC_filelist   ="";
//const std::string pp_data_filelist ="testFile_2015_5p02TeV_pp_data_forests.txt";
//int iii=0;
//std::cout<<"here"<<iii<<std::endl;iii++;

// the macro
int deriveDijetResponse(int startfile = defStartFile, int endfile = defEndFile, std::string infile_Forest=defInFilelist,
			int radius = 4, bool isMC = defIsMCderive, 
			std::string outfile = defOutputName, bool debugMode = defDebugMode ){

  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl;
  const bool fastDebugMode = (debugMode)&&false; //if debugMode is off, fastDebugMode shouldn't be on
  if(fastDebugMode)std::cout<<"fastDebugMode is ON"<<std::endl;

  std::cout<<std::endl<<"///////////////////"<<std::endl;
  std::cout<<"reading files "<< startfile << " to " << endfile<<std::endl;
  if(debugMode)std::cout<<"radius =  " << radius<<std::endl;
  if(debugMode)std::cout<<"isMC   = " << isMC <<std::endl;
  std::cout<<"///////////////////"<<std::endl<<std::endl;  

  // TH1F *hRelResponse[nbins_pt];
  // TH1F *hMPFResponse[nbins_pt];  
  // TH1F *hAbsPhoResponse[nbins_pt];
  // TH1F *hMPFAbsPhoResponse[nbins_pt];

  TH1F *avgAHisto[nbins_pt][nbins_eta];
  TH1F *avgBHisto[nbins_pt][nbins_eta];
  TH1F *hAvgAbsPhoResponse[nbins_pt][nbins_eta];
  
  for(int i=0; i<nbins_pt; i++){//pt loop

    // hRelResponse[i] = new TH1F(Form("hRelResponse_pt%d",i),"",nbins_eta,xbins_eta);
    // hMPFResponse[i] = new TH1F(Form("hMPFResponse_pt%d",i),"",nbins_eta,xbins_eta);
    // hAbsPhoResponse[i] = new TH1F(Form("hAbsPhoResponse_pt%d",i),"",nbins_eta,xbins_eta);
    // hMPFAbsPhoResponse[i] = new TH1F(Form("hMPFAbsPhoResponse_pt%d",i),"",nbins_eta,xbins_eta);

    for(int j=0; j<nbins_eta; j++){//eta loop
      avgAHisto[i][j] = new TH1F(Form("avgAHisto_pt%d_eta%d",i,j),"",50,-1,1);
      avgBHisto[i][j] = new TH1F(Form("avgBHisto_pt%d_eta%d",i,j),"",50,-1,1);
      hAvgAbsPhoResponse[i][j] = new TH1F(Form("hAvgAbsPhoResponse_pt%d_eta%d",i,j),"",50,-1,1);
    }//end eta loop
  }//end pt loop


  double avgA[nbins_pt][nbins_eta];
  TH1D *h_avgA[nbins_pt][nbins_eta];
  int nEntries[nbins_pt][nbins_eta];//nEntries imples nEntriesA, i think
  TH1D *h_nEntries[nbins_pt][nbins_eta];

  double avgB[nbins_pt][nbins_eta];
  TH1D *h_avgB[nbins_pt][nbins_eta];
  int nEntriesB[nbins_pt][nbins_eta];
  TH1D *h_nEntriesB[nbins_pt][nbins_eta];

  double avgAbsPhoResponse[nbins_pt][nbins_eta];
  TH1D *h_avgAbsPhoResponse[nbins_pt][nbins_eta];
  int nEntriesAbs[nbins_pt][nbins_eta];
  TH1D *h_nEntriesAbs[nbins_pt][nbins_eta];
  
  for(int i=0; i<nbins_pt; i++){//pt loop
    for(int j=0; j<nbins_eta; j++){//eta loop

      avgA[i][j]=0;      nEntries[i][j]=0;
      h_avgA[i][j]=new TH1D(Form("h_avgA_%d_%d",i,j),"",2, 0, 2);;
      h_nEntries[i][j]=new TH1D(Form("h_nEntries_%d_%d",i,j),"",2, 0, 2);;

      avgB[i][j]=0;      nEntriesB[i][j]=0;
      h_avgB[i][j]=new TH1D(Form("h_avgB_%d_%d",i,j),"",2, 0, 2);;
      h_nEntriesB[i][j]=new TH1D(Form("h_nEntriesB_%d_%d",i,j),"",2, 0, 2);;

      avgAbsPhoResponse[i][j]=0;      nEntriesAbs[i][j]=0;
      h_avgAbsPhoResponse[i][j]=new TH1D(Form("h_avgAbsPhoResponse_%d_%d",i,j),"",2, 0, 2);;
      h_nEntriesAbs[i][j]=new TH1D(Form("h_nEntriesAbs_%d_%d",i,j),"",2, 0, 2);;      

    }//end etaloop
  }//end pt loop
  
  // open filelist, loop to the starting file
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;  
  if(debugMode)std::cout<<"filelist used is "<<infile_Forest<<std::endl;
  for(int ifile = 0;ifile<startfile;++ifile) instr_Forest >> filename_Forest;

  
  //loop over files 
  for(int ifile = startfile; ifile<endfile; ++ifile){//file loop
    
    instr_Forest>>filename_Forest;//grab file
    if(atMIT)filename_Forest=hadoopDir+filename_Forest;
    else filename_Forest=xrootdDirect+filename_Forest;

    std::cout<<"opening trees from file "<<filename_Forest<<std::endl;
    
    TFile *fin = TFile::Open(filename_Forest.c_str());
    TTree *jtTree  = (TTree*)fin->Get( Form("ak%dPFJetAnalyzer/t",radius) );
    TTree *pfTree  = (TTree*)fin->Get( "pfcandAnalyzer/pfTree"            );
    TTree *phoTree = (TTree*)fin->Get( "ggHiNtuplizer/EventTree"          );

    int nref; 
    float pt_F[1000],  eta_F[1000],  phi_F[1000], rawpt_F[1000], m_F[1000];
    float eSum[1000], phoSum[1000];

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
    if(fastDebugMode)jetEntries=1000;
    if(debugMode)std::cout<<"jetEntries = "<< jetEntries <<std::endl;

    std::cout<<"beginning jet loop..."<< std::endl;    
    for(int ientry=0; ientry<jetEntries; ientry++){//jet loop

      jtTree->GetEntry(ientry);
      pfTree->GetEntry(ientry);
      //if(debugMode)std::cout<<"successfully grabbed entry from pftree"<<std::endl;
      phoTree->GetEntry(ientry);

      if(nPFpart_>10000) std::cout<<"warning! nPF: "<<nPFpart_<<std::endl;
      const double minJetPt = 40;
      
      for(int i=0; i<nPho; i++){//loop over photons
	if(phoEt->at(i)>30 && abs(phoEta->at(i))<1.4442){//photons cut
	  //if(debugMode)std::cout << "photon found, Et: " << phoEt->at(i) << std::endl;

	  //find a corresponding jet
	  for(int ijet=0; ijet<nref-1; ijet++){//loop over jets	   
	    if( pt_F[ijet]>=minJetPt && //jet cut
		abs(phoPhi->at(i)-phi_F[ijet])>2.95 && //angular separation photon+jet
		pt_F[ijet+1]/phoEt->at(i)<0.2 ){//jetpt/phoEt ratio
	      if(debugMode)std::cout << "jet photon pair found" << std::endl;	    

	      TLorentzVector missEt = findMissEt(nPFpart_, pfId_, pfPt_, 
						 pfEta_, pfPhi_, 
						 eSum, phoSum, 
						 nref, pt_F, rawpt_F, 
						 eta_F, phi_F, m_F);
	      int etaBin = findBin(phoEta->at(i),nbins_eta,xbins_eta);
	      int ptBin = findBin(phoEt->at(i),nbins_pt,xbins_pt);
	      TLorentzVector phoVec(phoEt->at(i),phoEta->at(i),phoPhi->at(i),0.);
	      double num = missEt.Dot(phoVec);
	      
	      avgAbsPhoResponse[ptBin][etaBin] += (1+(num/phoVec.Dot(phoVec)));
	      h_avgAbsPhoResponse[ptBin][etaBin]->Fill((1+(num/phoVec.Dot(phoVec))));
	      hAvgAbsPhoResponse[ptBin][etaBin]->Fill(1+(num/phoVec.Dot(phoVec)));
	      nEntriesAbs[ptBin][etaBin]++;
	      h_nEntriesAbs[ptBin][etaBin]->Fill(1);

	    }//end photon+jetcut
	  }//end jet loop
	}//end photon cut
      }//end photon loop

      // standard Rrel method
      if(nref>=2){//jet cut

	//if(debugMode)std::cout<<"two or more jets!"<<std::endl;
	int rJet, pJet;
	if( abs(eta_F[0])<1.3 && abs(eta_F[1])<1.3 ){ 
	  pJet = rand()%2;  rJet = (pJet+1)%2;
	}
	else if( abs(eta_F[0])<1.3 && abs(eta_F[1])<5. ){
	  pJet = 1; rJet = 0;
	}
	else if( abs(eta_F[1])<1.3 && abs(eta_F[0])<5. ){
	  pJet = 0; rJet = 1;
	}
	else continue;
	
	double avgPt = 0.5*(pt_F[rJet]+pt_F[pJet]);
	if(nref>2){//3 or more jets
	  if(pt_F[2]/avgPt > 0.2){
	    continue;
	  }
	}

	if(pt_F[rJet]<minJetPt) continue;

	double dphi = abs(phi_F[rJet]-phi_F[pJet]);
	if(dphi>(2*3.14159)) dphi-=(2*3.14159);
	if(dphi < 2.7) continue;

	// find the bin the jet belongs to
	int etaBin = findBin(eta_F[pJet],nbins_eta,xbins_eta);
	int ptBin = findBin(avgPt,nbins_pt,xbins_pt);
	//if(etaBin>xbins_eta[nbins_eta] || ptBin>xbins_pt[nbins_pt]) std::cout << "OH NO! Bin Mismatch!" << std::endl;
	avgA[ptBin][etaBin] += 0.5*(pt_F[pJet]-pt_F[rJet])/avgPt;
	h_avgA[ptBin][etaBin]->Fill(0.5*(pt_F[pJet]-pt_F[rJet])/avgPt);
	avgAHisto[ptBin][etaBin]->Fill(0.5*(pt_F[pJet]-pt_F[rJet])/avgPt);
	nEntries[ptBin][etaBin]++;
	h_nEntries[ptBin][etaBin]->Fill(1);

	// Starting MPF Method, missET and a jet loop
	TLorentzVector missEt = findMissEt(nPFpart_, pfId_, pfPt_, 
					   pfEta_, pfPhi_, 
					   eSum, phoSum,
					   nref, pt_F, rawpt_F, 
					   eta_F, phi_F, m_F);
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
  TFile *fout = new TFile(Form("%s",outfile.c_str()),"recreate");
  fout->cd();

  // write histos in serarate pt/eta bins
  if(debugMode)std::cout<< "writing histograms to file" <<std::endl;
  for(int i=0; i<nbins_pt; i++){
    for(int j=0; j<nbins_eta; j++){

      avgAHisto[i][j]->Write();
      avgBHisto[i][j]->Write();
      hAvgAbsPhoResponse[i][j]->Write();

      h_avgA[i][j]->Write();
      h_nEntries[i][j]->Write();
      h_avgB[i][j]->Write();
      h_nEntriesB[i][j]->Write();
      h_avgAbsPhoResponse[i][j]->Write();

      h_nEntriesAbs[i][j]->Write();
    }
  }

  fout->Close();
  return 0;
} //end deriveDijetResponse


// sumDijetResponse
// default inputs
/////
const std::string defSumInput="relDijetResponse_data_calo.root";
const bool defIsMCsum=false;
const bool defDoDraw=false;
/////

// the macro
int sumDijetResponse(std::string filename=defSumInput , bool isMC=defIsMCsum , bool doDraw=defDoDraw ){

  TH1F *hRelResponse[nbins_pt];// Rel
  TH1F *hMPFResponse[nbins_pt];// MPF

  //TH1F *hAbsPhoResponse[nbins_pt]; // AbsPho
  //TH1F *hMPFAbsPhoResponse[nbins_pt]; // MPFAbsPho
  //TH1F *hAvgAbsPhoResponse[nbins_pt][nbins_eta]; //AvgAbsPho

  TH1F *avgAHisto[nbins_pt][nbins_eta]; // avg A
  TH1F *avgBHisto[nbins_pt][nbins_eta]; // avg B

  //indicate the type of thing we're reading
  std::string type = "";
  if(!isMC) type = "Data"; 
  else type = "PYTHIA_CUETP8M1";

  //open the file
  TFile *fin = new TFile(filename.c_str());
  
  //array initializing
  for(int i=0; i<nbins_pt; i++){ //pt loop initializes response histo arrays
    std::cout<<i<<std::endl;

    hRelResponse[i] = new TH1F(Form("hRelResponse_pt%d",i),"",nbins_eta,xbins_eta);
    hMPFResponse[i] = new TH1F(Form("hMPFResponse_pt%d",i),"",nbins_eta,xbins_eta);
    //hAbsPhoResponse[i] = new TH1F(Form("hAbsPhoResponse_pt%d",i),"",nbins_eta,xbins_eta);
    //hMPFAbsPhoResponse[i] = new TH1F(Form("hMPFAbsPhoResponse_pt%d",i),"",nbins_eta,xbins_eta);

    for(int j=0; j<nbins_eta; j++){//eta loop gets+clones avgAB histos from the input file
      std::cout<<j<<std::endl;

      avgAHisto[i][j] = (TH1F*)fin->Get( Form("avgAHisto_pt%d_eta%d",i,j) )->Clone(Form("avgAHisto_pt%d_eta%d",i,j));
      avgBHisto[i][j] = (TH1F*)fin->Get( Form("avgBHisto_pt%d_eta%d",i,j) )->Clone(Form("avgBHisto_pt%d_eta%d",i,j));
      //hAvgAbsPhoResponse[i][j] = (TH1F*)fin->Get(Form("hAvgAbsPhoResponse_pt%d_eta%d",i,j))->Clone(Form("hAvgAbsPhoResponse_pt%d_eta%d",i,j));

    }//end eta loop
  }//end pt loop

  //
  for(int i=0; i<nbins_pt; i++){//pt loop
    for(int j=0; j<nbins_eta; j++){//eta loop
      std::cout<<i<<" "<<j<<std::endl;
    
      //Ahisto stuff
      hRelResponse[i]->SetBinContent( j+1, ( 1+avgAHisto[i][j]->GetMean() )/( 1-avgAHisto[i][j]->GetMean() ) );
      hRelResponse[i]->SetBinError( j+1 , (hRelResponse[i]->GetBinContent(j+1)) * (1./sqrt(avgAHisto[i][j]->GetEntries()) ) );

      //Bhisto stuff
      hMPFResponse[i]->SetBinContent( j+1 ,(1+avgBHisto[i][j]->GetMean())/(1-avgBHisto[i][j]->GetMean()));
      hMPFResponse[i]->SetBinError( j+1 ,hMPFResponse[i]->GetBinContent(j+1)*(1./sqrt(avgBHisto[i][j]->GetEntries())));

      //PhoResponse stuff
      //hMPFAbsPhoResponse[i]->SetBinContent( j+1 ,hAvgAbsPhoResponse[i][j]->GetMean());
      //hMPFAbsPhoResponse[i]->SetBinError( j+1 ,hMPFAbsPhoResponse[i]->GetBinContent(j+1)*(1./sqrt(hAvgAbsPhoResponse[i][j]->GetEntries())));

      //count entries in avg AB histos
      int totEntriesA=0, totEntriesB=0;
      //int totEntriesAbs=0;

      for(int j=0; j<nbins_eta; j++){//eta loop

	totEntriesA+=avgAHisto[i][j]->GetEntries(); 
	totEntriesB+=avgBHisto[i][j]->GetEntries();
	//totEntriesAbs+=hAvgAbsPhoResponse[i][j]->GetEntries();

      }//eta loop

      //set entries in Rel/MPF response histo arrays after filling because... why?
      hRelResponse[i]->SetEntries(totEntriesA);
      hMPFResponse[i]->SetEntries(totEntriesB);
      //hMPFAbsPhoResponse[i]->SetEntries(totEntriesAbs);

    }//end eta loop
  }//end pt loop
  
  //open an output file for output plots if desired
  TFile *fout = new TFile(Form("simplePlots_%disMC.root",isMC),"recreate");
  fout->cd();

  //color matrix for later
  int color[5] = {1,2,4,8,20};

  for(int i=0; i<nbins_pt; i++){//pt loop
    std::cout<<i<<std::endl;

    hRelResponse[i]->SetMarkerColor(color[i]);
    hRelResponse[i]->SetLineColor(color[i]);
    hRelResponse[i]->SetTitle(Form("Rel, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    hRelResponse[i]->Write();

    hMPFResponse[i]->SetMarkerColor(color[i]);
    hMPFResponse[i]->SetLineColor(color[i]);
    hMPFResponse[i]->SetMarkerStyle(21);
    hMPFResponse[i]->SetLineStyle(2);
    //    hRelResponse[i]->SetTitle(Form("MPF, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    hMPFResponse[i]->SetTitle(Form("MPF, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    hMPFResponse[i]->Write();

    // hMPFAbsPhoResponse[i]->SetMarkerColor(color[i]);
    // hMPFAbsPhoResponse[i]->SetLineColor(color[i]);
    // hMPFAbsPhoResponse[i]->SetTitle(Form("MPF Abs, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    // hMPFAbsPhoResponse[i]->Write();
  }//end pt loop


  if(!doDraw) std::cout<<"not going to draw histograms"<<std::endl;
  else{
    std::cout<<"going to draw histograms"<<std::endl;

    TCanvas *c0 = new TCanvas("c0","",600,600);
    TLegend *l0 = new TLegend(0.4,0.8,0.7,0.9);

    hRelResponse[0]->Draw();
    l0->AddEntry(hRelResponse[0],Form("%g<p_{T}<%g",xbins_pt[0],xbins_pt[1]));

    for(int i=1; i<nbins_pt; i++){
      std::cout<<i<<std::endl;
      hRelResponse[i]->Draw("same");
      if(isMC) l0->AddEntry("",Form("%s",type.c_str()),"");
      l0->AddEntry(hRelResponse[i],Form("%g<p_{T}<%g",xbins_pt[i],xbins_pt[i+1]));
    }

    l0->Draw("Same");
    std::cout<<"saving c0"<<std::endl;
    //c0->SaveAs("test.png");
    //sprintf("Relresponse_%s.pdf",type.c_str());
    //std::cout<<type.c_str()<<std::endl;
    c0->SaveAs("test.pdf");
    std::cout<<"after savng"<<std::endl;
    //Printf("c0 saved");

    TCanvas *c1 = new TCanvas("c1","",600,600);
    c1->cd();

    TLegend *l1 = new TLegend(0.4,0.8,0.7,0.9);
    hMPFResponse[0]->Draw();
    l1->AddEntry(hMPFResponse[0],Form("%g<p_{T}<%g",xbins_pt[0],xbins_pt[1]));

    for(int i=1; i<nbins_pt; i++){
      std::cout<<i<<std::endl;
      hMPFResponse[i]->Draw("same");
      if(isMC) l1->AddEntry("",Form("%s",type.c_str()),"");
      l1->AddEntry(hMPFResponse[i],Form("%g<p_{T}<%g",xbins_pt[i],xbins_pt[i+1]));
    }
    l1->Draw("Same");

    std::cout<<"saving c1"<<std::endl;
    c1->SaveAs(Form("MPFresponse_%s.pdf",type.c_str()),"RECREATE");

    // TCanvas *c2 = new TCanvas("c2","",600,600);
    // c2->cd();
    // hMPFAbsPhoResponse[0]->Draw();
    //
    // TLegend *l2 = new TLegend(0.4,0.8,0.7,0.9);
    // l2->AddEntry(hMPFAbsPhoResponse[0],Form("%g<p_{T}<%g",xbins_pt[0],xbins_pt[1]));
    //
    // for(int i=1; i<nbins_pt; i++){
    //   std::cout<<i<<std::endl;
    //   hMPFAbsPhoResponse[i]->Draw("same");
    //   if(isMC) l1->AddEntry("",Form("%s",type.c_str()),"");
    //   l2->AddEntry(hMPFAbsPhoResponse[i],Form("%g<p_{T}<%g",xbins_pt[i],xbins_pt[i+1]));
    // }
    //
    // l2->Draw("Same");
    // std::cout<<"saving c2"<<std::endl;
    // c2->SaveAs(Form("MPFABSPhoresponse_%s.pdf",type.c_str()),"RECREATE");
  }

  //close the output file
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
  else if( argc==9 && doDerive ){   // expected use for deriveDijetResponse
    rStatus=-1;
    int startfile=atoi(argv[2]); int endfile=atoi(argv[3]); std::string inputFileList=argv[4];
    int jetRadius=atoi(argv[5]); bool isMC=(bool)atoi(argv[6]); 
    std::string outputFileName=argv[7]; bool debug=(bool)atoi(argv[8]);      
    rStatus = deriveDijetResponse( startfile, endfile, inputFileList,
				   jetRadius, isMC,
				   outputFileName, debug);    
    return rStatus;
  }
  else if( argc==5 && doSum ){   // expected use for sumDijetResponse
    rStatus=-1;
    std::string inputFileName=argv[2]; bool isMC=(bool)atoi(argv[3]); bool doDraw=(bool)atoi(argv[4]);
    rStatus = sumDijetResponse( inputFileName, isMC, doDraw);    
    return rStatus;
  }
  else {  // inputs don't match any expected use cases
    std::cout<<"not enough arguments ";//<<std::endl;
    if (doDerive){
      std::cout<<"for deriveDijetResponse, do..."<<std::endl;
      std::cout<<"./DijetResponse.exe 0 ";
      std::cout<<"<startfile> <endfile> <inputFileList> <jetRadius> <isMC> <outputFileName> <debug>";
    }
    if (doSum){
      std::cout<<"for sumDijetResponse, do..."<<std::endl;
      std::cout<<"./DijetResponse.exe 1 ";
      std::cout<<"<inputFile> <isMC> <doDraw>";
    }
    std::cout<<std::endl<<std::endl;
    return rStatus;
  }
}//end main
