#include "printQAPlots_ppDataOnly.h"

// Input relevant
const std::string CMSSW_BASE =
  "/net/hisrv0001/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/";
const std::string inputDir =
  "src/readFiles/ppData/saved_outputCondor/";

//highPtJetTrig 5.02TeV 2015 prompt reco
const double integratedLuminosity=25.8*pow(10.,9.)*0.99;// 25.8 inv. picobarns to inv. millibarns
const std::string inputFilename=
  //"readFiles_ppData_5p02TeV_HighPtJetTrig_06-10-16__allFiles_vzHasNoTrigWeight/HighPtJetTrig_ak4PF-allFiles.root";
  //"readFiles_ppData_5p02TeV_HighPtJetTrig_08-23-16__noTrigWeights/HighPtJetTrig_ak4PF-allFiles.root";
  "readFiles_ppData_5p02TeV_HighPtJetTrig_08-23-16__vzWithAndWithoutWeights/HighPtJetTrig_ak4PF-allFiles.root";
const std::string fullFilename=CMSSW_BASE+inputDir+inputFilename;

// output relevant
const std::string thePDFFileName=
  "printQAPlots_ppData_HighPtJets_PromptReco2015_JetQA-vz-TrigComb_8.23.16.pdf";

// other useful things
const std::string radius="4";
const std::string HLTName[]={"HLT40_","HLT60_","HLT80_","HLT100_","HLTComb_","TrgCombTest_"};
const std::string HLThTitle[]={"HLT40 ","HLT60 ","HLT80 ","HLT100 ","HLTComb ","Kurt Meth. HLTComb "};
const int N_HLTNames=sizeof(HLTName)/sizeof(std::string);

const bool debugMode=true;


int main(int argc, char *argv[]){

  int rStatus = -1;
  if( argc!=1 ) {//no input arguments, error
    std::cout<<"settings hard coded, just do ./readFiles_ppMC_printPlots.exe"<<std::endl;
    return rStatus;
  }
  rStatus=1;

  //input file
  std::string filename=fullFilename;
  std::cout<< "filename is " << filename <<std::endl;
  TFile *fin = new TFile(filename.c_str());
  
  //get output file ready
  TCanvas *temp_canv = new TCanvas("temp", "temp", 1200, 600);
  std::string open_thePDFFileName=thePDFFileName+"[";
  std::string close_thePDFFileName=thePDFFileName+"]";
  temp_canv->Print( open_thePDFFileName.c_str() );      //open the .pdf file we want to print to
  temp_canv->cd();


  // evt plots
  //vz
  {
    std::string theHistName="hVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;

    TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
    theJetQAHist->SetTitle("data Event QA, vz, without TrigWeight");
    theJetQAHist->SetYTitle("millibarn/bin");
    theJetQAHist->SetXTitle("vz (cm)");

    theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));
    theJetQAHist->Scale(1/integratedLuminosity);

    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  {
    std::string theHistName="hWeightedVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;

    TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
    theJetQAHist->SetTitle("data Event QA, vz, with TrigWeight");
    theJetQAHist->SetYTitle("millibarn/bin");
    theJetQAHist->SetXTitle("vz (cm)");

    theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));
    theJetQAHist->Scale(1/integratedLuminosity);

    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  {
    std::string theHistName_num="hWeightedVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName_num<<std::endl;
    TH1F* theJetQAHist_num= (TH1F*)fin->Get( theHistName_num.c_str() );

    std::string theHistName_den="hVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName_den<<std::endl;
    TH1F* theJetQAHist_den= (TH1F*)fin->Get( theHistName_den.c_str() );

    TH1F* theRatio=theJetQAHist_num;
    theRatio->Divide(theJetQAHist_den);
    theRatio->SetTitle("data Event QA, vz w/ over vz w/o TrigWeight");
    theRatio->SetYTitle("weighted/unweighted");
    theRatio->SetXTitle("vz (cm)");
    theRatio->SetAxisRange(0., 18., "Y");

    temp_canv->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  //----------------------

  temp_canv->SetLogy(1);

  // jetQA plots ----------------------
  //{0,1}wJetID
  for(int j=0;j<N_vars;j++){      
    for(int i=0;i<2;i++){
      //e.g. hJetQA_{0,1}wJetID_{jtpt,rawpt,jteta...}
      std::string theHistName="hJetQA_"+std::to_string(i)+"wJetID_"+var[j];
      if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;

      TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
      //temp_canv->SetLogy(1);

      if(var[j]=="jteta") theJetQAHist->SetAxisRange(-2.5,2.5,"X");
      else if(var[j]=="xj"||var[j]=="Aj") {
	theJetQAHist->SetAxisRange(0.,1.05,"X");
	temp_canv->SetLogy(0);
      }
      

      std::string theHistTitle="Data JetQA "+var[j]+" ";
      if(i==0)theHistTitle+="w/o JetID";
      else theHistTitle+="w/ JetID";
      theJetQAHist->SetTitle(theHistTitle.c_str());
      
      theJetQAHist->SetYTitle("Cross-section (millibarn)/bin");
      
      theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));
      theJetQAHist->Scale(1/integratedLuminosity);
      
      theJetQAHist->Draw();
      temp_canv->Print( thePDFFileName.c_str() );
    }
  }
  //----------------------
  
  temp_canv->SetLogy(1);
  // HLT/L1 and Comb plots ----------------------
  //hpp_{HLT40,HLT60,HLT80,HLT100,HLTComb}_{no,}JetID_R4_20_eta_20
  for(int j=0;j<N_HLTNames;j++){
    for(int i=0;i<2;i++){
      
      std::string theHistName="hpp_"+HLTName[j];
      if(i==0)theHistName+="no";
      theHistName+="JetID_R"+radius+"_"+(std::string)etaWidth;
      if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
      
      TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
      
      std::string theHistTitle=HLThTitle[j];
      if(i==0)theHistTitle+="w/o JetID R";
      else theHistTitle+="w/ JetID R";      
      theHistTitle+=radius+" "+(std::string)etaWidth;
      theJetQAHist->SetTitle( theHistTitle.c_str() );
      theJetQAHist->SetYTitle("Cross-section (millibarn)/bin");

      theJetQAHist->SetAxisRange(0,1300.,"X");
      theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));
      theJetQAHist->Scale(1/integratedLuminosity);

      theJetQAHist->Draw();
      temp_canv->Print( thePDFFileName.c_str() );
    }
  }
  //----------------------

  fin->Close();
  temp_canv->Print( close_thePDFFileName.c_str() );//close the output file
  return 0;
}
