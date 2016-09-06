#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <string.h>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TStopwatch.h"


// definition and default values
const std::string ppMCfilelistFolder="../filelists/ppMC/";
//const std::string ppDatafilelistFolder="../filelists/ppData/";
const std::string defDatasetFilelist=ppMCfilelistFolder+"5p02TeV_Py8_CUETP8M1_QCDjet15_2Files_debug_forests.txt";
const std::string defOutputWeightFile="evtPthatWeights_defOutput.txt";
int evtPthatWeights(std::string infile_Forest=defDatasetFilelist , std::string outputWeightFile=defOutputWeightFile);


// constants
//source for bins+weights: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiForest#Dijet_Cross_Sections_for_reweigh
const bool doVzCuts=false, doNoiseFilterCuts=false, doEvtCuts=doVzCuts||doNoiseFilterCuts;
const double pthatBins[] = {15.,30.,50.,80.,120.,170.,220.,280.,370.,460.,540.,9999.};
const double xs2015[] = {5.269E-01, 3.455E-02, 4.068E-03, 4.959E-04, 7.096E-05, 1.223E-05, 3.031E-06, 7.746E-07, 1.410E-07, 3.216E-08, 1.001E-08, 0.0};
const int Npthats=sizeof(pthatBins)/sizeof(double);//should be 12
const int NpthatBins=Npthats-1;//should be 11


// ---------------------------------------------------------------------------------------------------------------
// evtPthatWeights
// compute the MC QCD dijet weights based on input filelist.
int evtPthatWeights( std::string infile_Forest , std::string outputWeightFile ){
  TStopwatch timer;  timer.Start();

  // basic info the screen
  std::cout<<std::endl<<"///////////////////"<<std::endl;
  std::cout<<"reading filelist "<<infile_Forest<<std::endl;
  std::cout<<"///////////////////"<<std::endl<<std::endl;

  // initialize variables
  bool endOfFilelist=false; 
  int n[NpthatBins]={0};
  int fileCount=0;

  // open filelist, loop over files and close them to save memory
  std::cout<<"looping over files..."<<std::endl;
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);

  while(!endOfFilelist){
    std::string filename_Forest;

    instr_Forest>>filename_Forest;
    if(filename_Forest==""){
      std::cout<<"end of filelist!"<<std::endl; 
      endOfFilelist=true; 
      break;    }//end of filelist condition

    if(fileCount%10==0)std::cout<<"opening file #"<<fileCount<<",  "<<filename_Forest<<std::endl;
    TFile *fin = TFile::Open(filename_Forest.c_str());     

    // grab tree(s), declare branches if needed
    TTree *t = (TTree*)fin->Get( "ak4PFJetAnalyzer/t"  );

    TTree *evt=NULL;
    if(doEvtCuts){evt=(TTree*)fin->Get( "HiEvtAnalyzer/HiEvt"  );
      t->AddFriend(evt);}
    

    // loop over pthat bins for this file
    for(int i = 0; i < NpthatBins; ++i){

      
      std::string theCut=formCutString(pthatBins[i],pthatBins[i+1],doNoiseFilterCuts,doVzCuts);
	//	"(pthat>="+std::to_string((int)pthatBins[i])+".)&&"+
	//	"(pthat<"+std::to_string((int)pthatBins[i+1])+".)";//&&"+
	//      if(doNoiseFilterCuts)theCut+=
	//      			   "&&(pBeamScrapingFilter&&HBHENoiseFilterResultRun2Loose&&pPAprimaryVertexFilter)";
	//if(doVzCuts)theCut+=
	//		      "&&(vz<15.&&vz>-15.)";
      
      
      
      TCut pthatCut(theCut.c_str());
      n[i] += t->GetEntries(pthatCut);

    }
    fin->Close();
    if(fileCount%10==0)std::cout<<"closing file #"<<fileCount<<",  "<<filename_Forest<<std::endl;
    fileCount++;
  }//end file loop
  assert(endOfFilelist);//otherwise why bother, make sure i'm hitting all the files



  std::cout<<"all files read, computing weights..."<<std::endl;
  std::ofstream outputFile; outputFile.open(outputWeightFile);
  for(int i=0 ; i<NpthatBins;i++){

    float xsDiff=xs2015[i]-xs2015[i+1], pthatweight;

    //protect against NaN and compute weight+output the number to the output text file
    if(n[i])pthatweight=xsDiff/n[i];
    else pthatweight=-1;

    std::string theCut=formCutString(pthatBins[i],pthatBins[i+1],doNoiseFilterCuts,doVzCuts);
    //	"(pthat>="+std::to_string((int)pthatBins[i])+".)&&"+
    //	"(pthat<"+std::to_string((int)pthatBins[i+1])+".)";//&&"+
    //      if(doNoiseFilterCuts)theCut+=
    //      			   "&&(pBeamScrapingFilter&&HBHENoiseFilterResultRun2Loose&&pPAprimaryVertexFilter)";
    //if(doVzCuts)theCut+=
    //		      "&&(vz<15.&&vz>-15.)";

    std::cout<<"for cut = "<<theCut<<std::endl;
    std::cout<<"there are "<<n[i]<<" events and ";
    std::cout<<"the weight is "<<pthatweight<<" mb"<<std::endl<<std::endl;
    
    //output weights to text file
    outputFile<<theCut<<", n["<<i<<"]="<<n[i]<<", ";  
    outputFile<<"xsDiff="<<xsDiff<<" mb^-3, ";
    outputFile<<"weight="<<pthatweight; 
    outputFile<<"\n";   }//end pthatbin loop for weight CPU


  std::cout<<"done computing weights, closing output file!"<<std::endl;
  outputFile.close();

  timer.Stop();
  std::cout<<"CPU time per file (sec)  = "<<(Float_t)timer.CpuTime()/fileCount<<std::endl;
  std::cout<<"Real time per file (sec) = "<<(Float_t)timer.RealTime()/fileCount<<std::endl;
  
  return 0;
} // end evtPthatWeights

std::string formCutString(double pthat_i,double pthat_iadd1,bool doNoiseFilterCuts,bool doVzCuts){
  
  std::string theCut;//formCutString(pthatBins[i],pthatBins[i+1],doNoiseFilterCuts,doVzCuts);
  theCut="(pthat>="+std::to_string((int)pthatBins[i])+".)&&" +
    "(pthat<"+std::to_string((int)pthatBins[i+1])+".)";
  
  if(doNoiseFilterCuts){theCut+="&&(pBeamScrapingFilter"+
      "&&HBHENoiseFilterResultRun2Loose"
      "&&pPAprimaryVertexFilter)";}
  
  if(doVzCuts)theCut+="&&(vz<15.&&vz>-15.)";

  return theCut; }


////// main //////
// acts as the frontend control for .exe file
int main(int argc, char* argv[]){

  // error, not enough arguments
  int rStatus = -1;  
  if(argc!=3&&argc!=1){
    std::cout<<"to use evtPthatWeights, you need a filelist. Then do..."<<std::endl;
    std::cout<<"./evtPthatWeights.exe <inputFileList> <outputWeightFile>"<<std::endl;
    std::cout<<"or leave arguments blank for default tests"<<std::endl;
    std::cout<<"rStatus="<<rStatus<<std::endl;
    return rStatus;  }
  
  if (argc==1) rStatus=evtPthatWeights();
  else{  std::string datasetFilelist =argv[1], outputWeightFile=argv[2];
    rStatus = evtPthatWeights( datasetFilelist, outputWeightFile);  }
  
  return rStatus;
}
