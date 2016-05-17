// Originally written by: Kurt Jung, December 2015
// Overhaul, Ian Laflotte, May 1st 2016
// Macro for deriving the dijet relative and absolute response from pp data/MC forests
// Uses JME-13-004 as a reference - "A" and "B" formulae are defined there
//first look for photon balance (AN2013-179 / JME-13-004)

///////////////////////////////////////////
////////// includes/preprocessor //////////
///////////////////////////////////////////

// C includes			    
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <unistd.h>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TH1F.h"



///////////////////////////////////////////////
////////// (initializa/declara)tions //////////
///////////////////////////////////////////////

// helper functions
int findEtaBin(float eta, int nbins_eta, const double* xbins_eta);
int findPtBin(float pt, int nbins_pt, const double* xbins_pt);
TLorentzVector findMissEt(int nPFpart_, std::vector<int>* pfId_, std::vector<float>* pfPt_, 
			  std::vector<float>* pfEta_, std::vector<float>* pfPhi_,
			  float* eSum, float *phoSum,            
			  int nref,    float* pt_F, float* rawpt_F, 
			  float* eta_F, float* phi_F, float* m_F);
// main macros
int deriveDijetResponse(int startfile, int endfile, std::string infile_Forest,
			int radius, bool isMC, 
			std::string outfile, bool debugMode);
int sumDijetResponse(std::string filename, bool isMC);

//for convenience during testing, coding, etc.
const bool atMIT  = true;
const std::string hadoopDir    = "/mnt/hadoop/cms";
const std::string xrootdDirect = "root://cmsxrootd.fnal.gov/";

//initializations
const double minJetPt = 40;

const double xbins_pt[] = {40, 60, 80, 110, 200, 1000};
const int nbins_pt = sizeof(xbins_pt)/sizeof(double)-1;

const double xbins_eta[] = {-5.191, -4.716, -4.363, -4.013,
			    -3.664, -3.314, -2.964,
			    -2.650, -2.322, -2.043, -1.830,
			    -1.653, -1.479, -1.305,
			    -1.131, -0.957, -0.783, -0.609,
			    -0.435, -0.261, -0.087,
			    0.000,
			    0.087,  0.261,  0.435, 0.609,
			    0.783,  0.957,  1.131, 
			    1.305,  1.479,  1.653,  1.830,
			    2.043,  2.322,  2.650, 
			    2.964,  3.314,   3.664, 4.013,
			    4.363,  4.716,  5.191};
const int nbins_eta = sizeof(xbins_eta)/sizeof(double)-1;

const int colors[] = {1,2,4,8,20};
const int N_colors= sizeof(colors)/sizeof(int);



///////////////////////////////////////////
////////// derive dijet response //////////
///////////////////////////////////////////
// default inputs to the routine
/////
const int defaultStartFile=0;
const int defaultEndFile=1; //exclusive boundary, range of files run over doesnt include endfile
const std::string defaultInFilelist = "debugFileList_2Files.txt";//"HighPtJet80_5p02TeV_pp_data_forests.txt";
/////
const int defaultRadius=4;
const bool defaultIsMC=false;
/////
const std::string defaultOutputName = "test_output.root"; 
const bool defaultDebugMode = true; 
/////
//const std::string pp_MC_filelist   ="";
//const std::string pp_data_filelist ="testFile_2015_5p02TeV_pp_data_forests.txt";
int iii=0;
//std::cout<<"here"<<iii<<std::endl;iii++;

// the macro
int deriveDijetResponse(int startfile = defaultStartFile, int endfile = defaultEndFile, std::string infile_Forest=defaultInFilelist,
			int radius = 4, bool isMC = defaultIsMC, 
			std::string outfile = defaultOutputName, bool debugMode = defaultDebugMode ){

  if(debugMode)std::cout<<"debugMode is ON"<<std::endl;

  std::cout<<std::endl<<"///////////////////"<<std::endl;
  std::cout<<"reading from file "<< startfile << " to file " << endfile<<std::endl;
  if(debugMode)std::cout<<"radius =  " << radius<<std::endl;
  if(debugMode)std::cout<<"n pt bins   = " << nbins_pt<<std::endl;
  if(debugMode)std::cout<<"n eta bins  = " << nbins_eta<<std::endl;
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
  std::cout << "opening file + trees" << std::endl;
  
  //loop over files 
  for(int ifile = startfile; ifile<endfile; ++ifile){//file loop
    
    instr_Forest>>filename_Forest;//grab file
    if(atMIT)filename_Forest=hadoopDir+filename_Forest;
    else filename_Forest=xrootdDirect+filename_Forest;
    if(debugMode)std::cout<<"opening file "<<filename_Forest<<std::endl;
    
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

    int jetEntries = jtTree->GetEntries();
    //std::cout<<"Num. Jet Entries = "<<jetEntries< std::endl;
    std::cout<<"jetEntries = "<<jetEntries<<std::endl;

    if(jetEntries > 1000 && debugMode) std::cout<< " warning! jetEntries= "<< jetEntries << std::endl; 
    if(debugMode)jetEntries=1000;
    if(debugMode)std::cout<<"beginning jet loop"<< std::endl;    

    // loop over jets
    for(int ientry=0; ientry<jetEntries; ientry++){//jet loop

      jtTree->GetEntry(ientry);
      pfTree->GetEntry(ientry);
      //if(debugMode)std::cout<<"successfully grabbed entry from pftree"<<std::endl;
      phoTree->GetEntry(ientry);

      if(nPFpart_>10000 && debugMode) std::cout << " warning! nPF: "<< nPFpart_ << std::endl;

      for(int i=0; i<nPho; i++){//loop over photons, look for photon
	if(phoEt->at(i)>30 && abs(phoEta->at(i))<1.4442){//photons cut
	  if(debugMode)std::cout << "photon found, Et: " << phoEt->at(i) << std::endl;
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
	      int etaBin = findEtaBin(phoEta->at(i),nbins_eta,xbins_eta);
	      int ptBin = findPtBin(phoEt->at(i),nbins_pt,xbins_pt);
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
	if(debugMode)std::cout<<"two or more jets!"<<std::endl;
 
	int rJet, pJet;
	if( abs(eta_F[0])<1.3 && abs(eta_F[1])<1.3 ){ 
	  pJet = rand()%2;
	  rJet = (pJet+1)%2;
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

	// phi matching
	double dphi = abs(phi_F[rJet]-phi_F[pJet]);
	if(dphi>(2*3.14159)) dphi-=(2*3.14159);
	if(dphi < 2.7) continue;

	// find the bin the jet belongs to
	int etaBin = findEtaBin(eta_F[pJet],nbins_eta,xbins_eta);
	int ptBin = findPtBin(avgPt,nbins_pt,xbins_pt);
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
	  int etaBin = findEtaBin(eta_F[pJet],nbins_eta,xbins_eta);
	  int ptBin = findPtBin(avgPt,nbins_pt,xbins_pt);
	  TLorentzVector jetVec(pt_F[rJet],eta_F[rJet],phi_F[rJet],m_F[rJet]);
	  avgB[ptBin][etaBin] += (missEt.Dot(jetVec) / ( 2*avgPt * jetVec.Mag()));
	  h_avgB[ptBin][etaBin]->Fill((missEt.Dot(jetVec) / ( 2*avgPt * jetVec.Mag())));
	  avgBHisto[ptBin][etaBin]->Fill((missEt.Dot(jetVec) / ( 2*avgPt * jetVec.Mag())));
	  nEntriesB[ptBin][etaBin]++;
	  h_nEntriesB[ptBin][etaBin]->Fill(1);
	}//end jet loop
      }//end 2 or more jets req
    }//exit jet loop

    fin->Close();
  }//end file loop

  // open output file
  TFile *fout = new TFile(Form("%s",outfile.c_str()),"recreate");
  fout->cd();

  // write histos in serarate pt/eta bins
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



////////////////////////////////////////
////////// sum dijet response //////////
////////////////////////////////////////

// default inputs to the routine
/////
const std::string defaultSumInput="relDijetResponse_data_calo.root";
/////

// the macro
int sumDijetResponse(std::string filename=defaultSumInput , bool isMC=defaultIsMC ){
  //garbage commands for test/debug
  filename+="1";
  bool haha=true&&isMC;
  std::cout<<haha;
  //garbage command for test/debug
  return 0;
}



///////////////
//// findBin //
///////////////
//int findBin(float value, const double* xbins){
//  int nbins = sizeof(xbins)/sizeof(double)-1;
//  for(int i=0; i<nbins; i++) if(value<xbins[i+1]) return i;
//  return nbins_eta-1;
//}

////////////////
// findEtaBin //
////////////////
int findEtaBin(float eta, 
	       int nbins_eta, const double* xbins_eta){
  for(int i=0; i<nbins_eta; i++) if(eta<xbins_eta[i+1]) return i;
  return nbins_eta-1;
}

///////////////
// findPtBin //
///////////////
int findPtBin(float pt, 
	      int nbins_pt, const double* xbins_pt){
  for(int i=0; i<nbins_pt; i++) if(pt<xbins_pt[i+1]) return i;  
  return nbins_pt-1;
}

////////////////
// findMissEt //
////////////////
TLorentzVector findMissEt(int nPFpart_, std::vector<int> * pfId_, std::vector<float> * pfPt_, 
			  std::vector<float> * pfEta_, std::vector<float> * pfPhi_,
			  float* eSum, float *phoSum,
                          int nref, float* pt_F, float* rawpt_F, 
			  float* eta_F, float* phi_F, float* m_F){
  
  // missEt = - Sum(Et you can see)
  TLorentzVector missEt(0,0,0,0);

  //missing Et from particles
  for(int i=0; i<nPFpart_; i++){//PF Loop
    TLorentzVector pfp;
    double pfMass=0;
    if(pfId_->at(i)==1 || pfId_->at(i)==4) pfMass = 0.1395702; 
    //pfp.SetPtEtaPhiM(pfPt[i],pfEta[i],pfPhi[i],pfMass);
    pfp.SetPtEtaPhiM(pfPt_->at(i),pfEta_->at(i),pfPhi_->at(i),pfMass);
    pfp.SetPz(0.); //2d projection in x-y plane
    missEt += pfp;
  }//end PF loop
  
  //missing Et from jets
  for(int i=0; i<nref; i++){//begin jet loop
    if(pt_F[i]<=15) continue;//quick cut before expensive computations
    TLorentzVector jtTmp(pt_F[i],eta_F[i],phi_F[i],m_F[i]);
    double jetEnergy=jtTmp.E();
    if( (eSum[i]+phoSum[i])/jetEnergy<0.9 ){
      TLorentzVector jt;
      jt.SetPtEtaPhiM(pt_F[i]-rawpt_F[i],eta_F[i],phi_F[i],m_F[i]);
      jt.SetPz(0.); //project to x-y plane
      missEt += jt;
    }
  }//end jet loop

//  //missing Et from jets
//  for(int i=0; i<nref; i++){//begin jet loop
//    if(pf_F[i]<=15) continue;//quick cut before expensive computations
//    TLorentzVector jt(pt_F[i],eta_F[i],phi_F[i],m_F[i]);
//    if( (eSum[i]+phoSum[i])/jt.E()<0.9 ){
//      jt.SetPtEtaPhiM(pt_F[i]-rawpt_F[i],eta_F[i],phi_F[i],m_F[i]);
//      jt.SetPz(0.); //project to x-y plane
//      missEt += jt;
//    }
//  }//end jet loop

  missEt*=-1; // rotate 180 deg, to get missing Et, 
  missEt.SetE(0.);  // only care about Et, not E
  //std::cout << "missEt:  Mag: " << missEt.Mag() << std::endl;  
  return missEt;
}


//////////////////////////
////////// main //////////
//////////////////////////
// acts as the frontend control for .exe file
// could use some improvement
// note, argc->argument count, argv->argument vector
int main(int argc, char *argv[]){

  std::cout<<std::endl;
  int rStatus = 1;  

  // no input, error!
  if( argc==1 ) {
    std::cout<<"no input, give at least one argument"<<std::endl<<std::endl;
    return rStatus;
  } 
  
  // non-existant function called, error!
  bool doDerive=(atoi(argv[1])==0);
  bool doSum=(atoi(argv[1])==1);    
  if(!doDerive && !doSum ){
    std::cout<<"bad input, first argument must be either... "<<std::endl;
    std::cout<<"0 for deriveDijetResponse"<<std::endl;
    std::cout<<"1 for sumDijetResponse"<<std::endl<<std::endl;
  }  
  
  // expected use for default-value tests
  if( argc==2 ){ 
    rStatus=-1;
    if(doDerive) rStatus=deriveDijetResponse();
    if(doSum) rStatus=sumDijetResponse();
    return rStatus;
  }
  // expected use for deriveDijetResponse
  else if( argc==9 && doDerive ){ 
    rStatus=-1;

    int startfile=atoi(argv[2]); int endfile=atoi(argv[3]); std::string inputFileList=argv[4];
    int jetRadius=atoi(argv[5]); bool isMC=(bool)atoi(argv[6]); 
    std::string outputFileName=argv[7]; bool debug=(bool)atoi(argv[7]);      
    rStatus = deriveDijetResponse( startfile, endfile, inputFileList,
				   jetRadius, isMC,
				   outputFileName, debug);    return rStatus;
  }
  // expected use for sumDijetResponse
  else if( argc==4 && doSum ){ 
    rStatus=-1;
    std::string inputFileName=argv[2]; bool isMC=(bool)atoi(argv[3]);
    rStatus = sumDijetResponse( inputFileName, isMC);    return rStatus;
  }
  // inputs don't match any expected use cases
  else {
    std::cout<<"to test the script with defaults, do..."<<std::endl;
    std::cout<<"./DijetResponse.exe <function>"<<std::endl<<std::endl;
    if (doDerive){
      std::cout<<"for deriveDijetResponse, do..."<<std::endl;
      std::cout<<"./DijetResponse.exe 0 ";
      std::cout<<"<startfile> <endfile> <inputFileList> <jetRadius> <isMC> <outputFileName> <debug>"<<std::endl<<std::endl;
    }
    if (doSum){
      std::cout<<"for sumDijetResponse, do..."<<std::endl;
      std::cout<<"./DijetResponse.exe 1 ";
      std::cout<<"<inputFile> <isMC>"<<std::endl<<std::endl;
    }
    return rStatus;
  }

}
