#include "readFiles_ppMC.h"

// Input relevant
const std::string CMSSW_BASE = 
  "/net/hisrv0001/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/";
const std::string inputMCDir = 
  "src/readFiles/ppMC/saved_outputCondor/";
const std::string inputDataDir = 
  "src/readFiles/ppData/saved_outputCondor/";

// pythia8 tune CUETP8M1, QCDDijet sample in pthat bins starting at 15 GeV,
const std::string inputMCFilename = 
  "readFiles_ppMC_5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_2016-07-09/QCDjetAllPtBins_ak4PF-allFiles.root";
const std::string fullMCFilename = CMSSW_BASE+inputMCDir+inputMCFilename;

// for vz plot only, for now
const double integratedLuminosity=25.8*pow(10.,9.)*0.99;// 25.8 inv. picobarns to inv. millibarns
const std::string inputDataFilename = 
  "readFiles_ppData_5p02TeV_HighPtJetTrig_2016-06-10_allFiles/HighPtJetTrig_ak4PF-allFiles.root";
const std::string fullDataFilename = CMSSW_BASE+inputDataDir+inputDataFilename;

// output relevant
const std::string thePDFFileName = 
  "readFiles_ppMC_printPlots_PY8_TuneCUETP8M1-QCDDijet_JetQA-vz-pthat_8.22.16.pdf";

// other useful things
const std::string radius="R4_";
const std::string pfRad_etaWidth=radius+"_"+(std::string)etaWidth;

const bool debugMode=true;


int main(int argc, char *argv[]){

  int rStatus = -1;
  if( argc!=1 ) {//no input arguments, error
    std::cout<<"settings hard coded, just do ./readFiles_ppMC_printPlots.exe"<<std::endl;
    return rStatus;
  }
  rStatus=1;


  //input file
  std::string MCfilename=fullMCFilename;
  std::cout<< "MCfilename is " << MCfilename <<std::endl;
  TFile *finMC = new TFile(MCfilename.c_str());
  

  //get output file ready
  TCanvas *temp_canv = new TCanvas("temp", "temp", 1200, 600);
  std::string open_thePDFFileName=thePDFFileName+"[";
  std::string close_thePDFFileName=thePDFFileName+"]";
  temp_canv->Print( open_thePDFFileName.c_str() );      //open the .pdf file we want to print to
  temp_canv->cd();


  // evt plots----------------------
  //vz
  {
    std::string theHistName="hVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;

    TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );

    if(debugMode)std::cout<<"bin width is "<< theJetQAHist->GetBinWidth(1) <<std::endl;
    theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));

    theJetQAHist->SetTitle("MC event QA, vz w/o weights");
    theJetQAHist->SetYTitle("Num. MC Events/bin");
    theJetQAHist->SetXTitle("vz (cm)");

    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  {
    std::string theHistName="hWeightedVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );

    if(debugMode)std::cout<<"bin width is "<< theJetQAHist->GetBinWidth(1) <<std::endl;
    theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));

    theJetQAHist->SetTitle("MC event QA, vz w/ weights");
    theJetQAHist->SetYTitle("millibarns/bin");
    theJetQAHist->SetXTitle("vz (cm)");

    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }

  //vz ratio, no weights
  {// since bin widths the same, no need to scale by the bin width
    //input file
    std::string Datafilename=fullDataFilename;
    std::cout<< "Datafilename is " << Datafilename <<std::endl;
    TFile *finData = new TFile(Datafilename.c_str());

    std::string data_Vz_Name="hVz";
    TH1F* data_Vz= (TH1F*)finData->Get( data_Vz_Name.c_str() );
    
    std::string MC_Vz_Name="hVz";
    TH1F* MC_Vz= (TH1F*)finMC->Get( MC_Vz_Name.c_str() );
    MC_Vz->Scale(data_Vz->Integral()/MC_Vz->Integral());
    
    TH1F* Vz_ratio=data_Vz;
    Vz_ratio->Divide(MC_Vz);
    Vz_ratio->SetTitle("Vz ratio, data/MC w/o weights, MC norm. to data");
    Vz_ratio->SetXTitle("vz (cm)");
    Vz_ratio->SetYTitle("data/MC");
 
    Vz_ratio->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
    finData->Close();
  }
  //vz ratio, with weights
  {
    //input file
    std::string Datafilename=fullDataFilename;
    std::cout<< "Datafilename is " << Datafilename <<std::endl;
    TFile *finData = new TFile(Datafilename.c_str());

    std::string data_Vz_Name="hVz";
    TH1F* data_Vz= (TH1F*)finData->Get( data_Vz_Name.c_str() );
    data_Vz->Scale(1/integratedLuminosity);

    std::string MC_Vz_Name="hWeightedVz";
    TH1F* MC_Vz= (TH1F*)finMC->Get( MC_Vz_Name.c_str() );
    MC_Vz->Scale(data_Vz->Integral()/MC_Vz->Integral());
    
    TH1F* Vz_ratio=data_Vz;
    Vz_ratio->Divide(MC_Vz);
    Vz_ratio->SetTitle("Vz ratio, data/MC w/ weights, MC norm. to data");
    Vz_ratio->SetXTitle("vz (cm)");
    Vz_ratio->SetYTitle("data/weightedMC");
 
    Vz_ratio->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
    finData->Close();
  }


  temp_canv->SetLogy(1);

  //pthat
  {
    std::string theHistName="hpthat";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
    theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));

    theJetQAHist->SetTitle("MC event QA, pthat, not weighted");
    theJetQAHist->SetYTitle("Num. Events/bin");
    theJetQAHist->SetXTitle("pthat (GeV)");

    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  {
    std::string theHistName="hWeightedpthat";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
    theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));

    theJetQAHist->SetTitle("MC event QA, pthat, weighted");
    theJetQAHist->SetYTitle("Cross-section (millibarn)/bin");
    theJetQAHist->SetXTitle("pthat (GeV)");

    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  //----------------------



  // JetQA plots ----------------------
  //{1,0}wJetID
  for(int j=0;j<N_vars;j++){      
    for(int i=0;i<2;i++){
      //e.g. hJetQA_{0,1}wJetID_{jtpt,rawpt,jteta...}
      std::string theHistName="hJetQA_"+std::to_string(i)+"wJetID_"+var[j];
      if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;

      TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );

      // variable specific
      if(var[j]=="jteta") theJetQAHist->SetAxisRange(-2.5,2.5,"X");
      if(var[j]=="xj"||var[j]=="Aj") {
	theJetQAHist->SetAxisRange(0.,1.05,"X");
	temp_canv->SetLogy(0);
      }
      else temp_canv->SetLogy(1);

      theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));
      theJetQAHist->SetYTitle("Cross-section (millibarn)/bin");

      std::string theHistTitle="weighted MC JetQA "+var[j]+" ";
      if(i==0)theHistTitle+="w/o JetID";
      else theHistTitle+="w/ JetID";      
      theJetQAHist->SetTitle(theHistTitle.c_str());

      theJetQAHist->Draw();
      temp_canv->Print( thePDFFileName.c_str() );
    }
  }
  //----------------------
















  //// unfolding plots----------------------
  //
  ////unfolding inputs
  //{  
  //  std::string theHistName="hpp_gen_"+pfRad_etaWidth;
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}
  //{  
  //  std::string theHistName="hpp_reco_"+pfRad_etaWidth;
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}
  //{  
  //  std::string theHistName="hpp_matrix_"+pfRad_etaWidth;
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH2F* theJetQAHist= (TH2F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}
  //
  ////MCClosure inputs
  //{  
  //  std::string theHistName="hpp_mcclosure_gen_"+pfRad_etaWidth;
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}
  //{  
  //  std::string theHistName="hpp_mcclosure_data_"+pfRad_etaWidth;
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}
  //{  
  //  std::string theHistName="hpp_mcclosure_data_train_"+pfRad_etaWidth;
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}
  //{  
  //  std::string theHistName="hpp_mcclosure_matrix_"+pfRad_etaWidth;
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH2F* theJetQAHist= (TH2F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}

  //// Jet Corrections----------------------
  ////JEC
  //{  
  //  std::string theHistName="";
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}
  ////JER
  //{  
  //  std::string theHistName="";
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}



  finMC->Close();
  temp_canv->Print( close_thePDFFileName.c_str() );//close the output file
  return 0;
}

