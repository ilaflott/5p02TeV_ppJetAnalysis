#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TStopwatch.h"


// definition and default values

//const std::string filelistFolder="/net/hisrv0001/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/readForests/filelists/";
const std::string filelistFolder="filelists/";
//const std::string defMCFilelist="test_readForests_ppMC_Py8_CUETP8M1_Official_forests_acrossBins.txt";
const std::string defMCFilelist="5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_Official_forests.txt";
const std::string defOutputWeightFile="evtPthatWeights_defOutput.txt";
const int defRadius=4;
const std::string defJetType="PF";
const float defpthatMin=0.;
//const float defpthatMin=220;
// constants
//const bool doVzCuts=false, doNoiseFilterCuts=false, doEvtCuts=doVzCuts||doNoiseFilterCuts;

//semi-private ppMC https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiForest#Dijet_Cross_Sections_for_reweigh
//const double pthatBins[] = {15.,30.,50.,80.,120.,170.,220.,280.,370.,460.,540.,9999.};
//const double xs2015[] = {5.269E-01, 3.455E-02, 4.068E-03, 4.959E-04, 7.096E-05, 1.223E-05, 3.031E-06, 7.746E-07, 1.410E-07, 3.216E-08, 1.001E-08, 0.0};

//official ppMC , xsec from test files seems like xsec units are different?
const double pthatBins[] = {15.,30.,50.,80.,120.,170.,220.,280.,370.,9999.};
const double xs2015[] = { 5.115E+08 , 3.734E+07 , 4.005E+06 , 5.543E+05 , 6.423E+04, 1.346E+04 , 3.091E+03 , 7.597E+02 , 1.416E+02, 0.0 };

int evtPthatWeights(std::string inputFilelist=defMCFilelist , std::string outputWeightFile=defOutputWeightFile,
		    const int radius=defRadius, std::string jetType=defJetType, float pthatMin = defpthatMin );

std::string formCutString(double pthat_i,double pthat_iadd1);//,bool doNoiseFilterCuts,bool doVzCuts);

const int Npthats=sizeof(pthatBins)/sizeof(double);//should be 12
const int NpthatBins=Npthats-1;//should be 11


// ---------------------------------------------------------------------------------------------------------------
// evtPthatWeights
// compute the MC QCD dijet weights based on input filelist.
int evtPthatWeights( std::string inputFilelist , std::string outputWeightFile,
		     const int radius, std::string jetType , float pthatMin ){
  TStopwatch timer;  timer.Start();
  
  // basic info the screen
  std::cout<<std::endl<<"///////////////////"<<std::endl;
  std::cout<<"reading filelist "<<inputFilelist<<std::endl;
  std::cout<<"///////////////////"<<std::endl<<std::endl;

  // initialize variables
  bool endOfFilelist=false; 
  int n[NpthatBins]={0};
  int fileCount=0;

  // open filelist, loop over files and close them to save memory
  std::cout<<"making TChain w filelist..."<<std::endl;
  inputFilelist=filelistFolder+inputFilelist;
  //inputFilelist=inputFilelist
  std::ifstream instr_Forest(inputFilelist.c_str(),std::ifstream::in);

  std::string jetTreeName="ak"+std::to_string(radius)+jetType+"JetAnalyzer/t";  
  TChain * finMC = new TChain( jetTreeName.data() );

  while(!endOfFilelist){

    std::string filename_Forest;
    instr_Forest>>filename_Forest;
    if(filename_Forest==""){
      std::cout<<"end of filelist!"<<std::endl<<std::endl; 
      endOfFilelist=true; 
      break;    }//end of filelist condition
    
    if(fileCount%100==0)std::cout<<std::endl<<"adding to TChain, file "<<filename_Forest<<std::endl<<std::endl;
    else if(fileCount%1==0)std::cout<<"adding file #"<<fileCount<<", to TChain "<<std::endl;

    finMC->AddFile( filename_Forest.c_str() );
    
    fileCount++;
    
    //if (fileCount%10==0){
    //  endOfFilelist=true;
    //  break;
    //}

  }//end file loop
  assert(endOfFilelist);
  
  // loop over pthat bins for this file
  std::cout<<"looping over pthat bins..."<<std::endl;
  for(int i = 0; i < NpthatBins; ++i){
    if(pthatBins[i] < pthatMin )continue;
    std::string theCut=formCutString(pthatBins[i],pthatBins[i+1]);
    n[i] += finMC->GetEntries(theCut.c_str());    
  }


  std::cout<<"all files read and done with loop, the pthat weights are..."<<std::endl<<std::endl;
  std::ofstream outputFile; outputFile.open(outputWeightFile);
  for(int i=0 ; i<NpthatBins;i++){

    float xsDiff=xs2015[i]-xs2015[i+1], pthatweight;

    //protect against NaN and compute weight+output the number to the output text file
    if(n[i])pthatweight=xsDiff/n[i];
    else pthatweight=-1;

    std::string theCut=formCutString(pthatBins[i],pthatBins[i+1]);

    std::cout<<"for cut = "<<theCut<<std::endl;
    std::cout<<"there are "<<n[i]<<" events and ";
    std::cout<<"the weight is "<<pthatweight<<" pb"<<std::endl<<std::endl;
    
    //output weights to text file
    // this mess is for easy copy/pasting for later when wanting to run parts of job at a time.
    outputFile<<theCut;
    if(i<=2)outputFile<<"    ";
    else if(i<=3)outputFile<<"   ";
    else if(i<=7)outputFile<<"  ";
    else if(i<=8)outputFile<<" ";
    outputFile<<", n["<<i<<"]="<<n[i]<<", ";
    outputFile<<"   ";
    //outputFile<<"xsDiff="<<xsDiff<<" mb^, ";
    outputFile<<"xsDiff="<<xsDiff<<" pb^, ";
    outputFile<<"weight="<<pthatweight<<" pb^, ";
    outputFile<<"\n";   }//end pthatbin loop for weight CPU


  std::cout<<"done computing weights, closing output file "<< outputWeightFile <<std::endl;
  outputFile.close();

  timer.Stop();
  std::cout<<"CPU time per file (sec)  = "<<(Float_t)timer.CpuTime()/fileCount<<std::endl;
  std::cout<<"Real time per file (sec) = "<<(Float_t)timer.RealTime()/fileCount<<std::endl;
  
  return 0;
} // end evtPthatWeights

std::string formCutString(double pthat_i,double pthat_iadd1){//,bool doNoiseFilterCuts,bool doVzCuts){
  
  std::string theCut;//formCutString(pthatBins[i],pthatBins[i+1],doNoiseFilterCuts,doVzCuts);
  theCut="(pthat>="+std::to_string((int)pthat_i)+".)&&" +
    "(pthat<"+std::to_string((int)pthat_iadd1)+".)";
  
  return theCut; 
}


////// main //////
// acts as the frontend control for .exe file
int main(int argc, char* argv[]){

  // error, not enough arguments
  int rStatus = -1;  
  if(argc!=6 && argc!=5 && argc!=1){
    std::cout<<"to use evtPthatWeights, you need a filelist. Then do..."<<std::endl;
    std::cout<<"./evtPthatWeights.exe <inputFileList> <outputWeightFile> <radius> <jetType>"<<std::endl;
    std::cout<<"or leave arguments blank for default tests"<<std::endl;
    std::cout<<"rStatus="<<rStatus<<std::endl;
    return rStatus;  }
  
  if (argc==1) rStatus=evtPthatWeights();
  else if (argc==5){  std::string datasetFilelist =argv[1], outputWeightFile=argv[2];
    int radius = std::atoi(argv[3]); std::string jetType= argv[4];
    rStatus = evtPthatWeights( datasetFilelist, outputWeightFile, radius, jetType);  }
  else {  std::string datasetFilelist =argv[1], outputWeightFile=argv[2];
    int radius = std::atoi(argv[3]); std::string jetType= argv[4]; float pthatMin= std::atof(argv[4]);
    rStatus = evtPthatWeights( datasetFilelist, outputWeightFile, radius, jetType, pthatMin);  }
  return rStatus;
}
