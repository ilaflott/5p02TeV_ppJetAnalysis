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
  //"readFiles_ppMC_5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_07-09-16__firstRun_allFiles/QCDjetAllPtBins_ak4PF-allFiles.root";
  //"readFiles_ppMC_5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_08-23-16__RaghavsVzWeights_noTrigWeights/QCDjetAllPtBins_ak4PF-allFiles.root";
  "readFiles_ppMC_5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_08-23-16__AnnasVzWeights_noTrigWeights/QCDjetAllPtBins_ak4PF-allFiles.root";

const std::string fullMCFilename = CMSSW_BASE+inputMCDir+inputMCFilename;

// for vz plot only, for now
const double integratedLuminosity=25.8*pow(10.,9.)*0.99;// 25.8 inv. picobarns to inv. millibarns
const std::string inputDataFilename = 
  //"readFiles_ppData_5p02TeV_HighPtJetTrig_06-10-16__allFiles_vzHasNoTrigWeight/HighPtJetTrig_ak4PF-allFiles.root";
  //"readFiles_ppData_5p02TeV_HighPtJetTrig_08-23-16__noTrigWeights/HighPtJetTrig_ak4PF-allFiles.root";
  "readFiles_ppData_5p02TeV_HighPtJetTrig_08-23-16__vzWithAndWithoutWeights/HighPtJetTrig_ak4PF-allFiles.root";
const std::string fullDataFilename = CMSSW_BASE+inputDataDir+inputDataFilename;

// output relevant
const std::string thePDFFileName = 
  //"printQAPlots_ppMC_PY8_TuneCUETP8M1-QCDDijet_JetQA-vz-pthat_8.23.16_RaghavsVzWeights.pdf";
  "printQAPlots_ppMC_PY8_TuneCUETP8M1-QCDDijet_JetQA-vz-pthat_8.23.16_AnnasVzWeights.pdf";

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

    theJetQAHist->SetTitle("MC evt QA, vz w/o weights");
    theJetQAHist->SetYTitle("Num. MC Evts/bin");
    theJetQAHist->SetXTitle("vz (cm)");

    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  {
    std::string theHistName="hpthatWeightedVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );

    if(debugMode)std::cout<<"bin width is "<< theJetQAHist->GetBinWidth(1) <<std::endl;
    theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));

    theJetQAHist->SetTitle("MC evt QA, vz w/ pthat weights only");
    theJetQAHist->SetYTitle("millibarns/bin");
    theJetQAHist->SetXTitle("vz (cm)");

    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  {
    std::string theHistName="hvzWeightedVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );

    if(debugMode)std::cout<<"bin width is "<< theJetQAHist->GetBinWidth(1) <<std::endl;
    theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));

    theJetQAHist->SetTitle("MC event QA, vz w/ vz weights only");
    theJetQAHist->SetYTitle("millibarns/bin");
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

    theJetQAHist->SetTitle("MC event QA, vz w/ all weights");
    theJetQAHist->SetYTitle("millibarns/bin");
    theJetQAHist->SetXTitle("vz (cm)");

    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  

  //vz, MC Only ratio plots
  {
    std::string theHistName_num="hvzWeightedVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName_num<<std::endl;
    TH1F* theJetQAHist_num= (TH1F*)finMC->Get( theHistName_num.c_str() );

    std::string theHistName_den="hVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName_den<<std::endl;
    TH1F* theJetQAHist_den= (TH1F*)finMC->Get( theHistName_den.c_str() );

    TH1F* theRatio=theJetQAHist_num;
    theRatio->Divide(theJetQAHist_den);
    theRatio->SetTitle("MC Evt QA, MC Only, vzWeights only");
    theRatio->SetYTitle("weighted/unweighted");
    theRatio->SetXTitle("vz (cm)");

    theRatio->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  {
    std::string theHistName_num="hpthatWeightedVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName_num<<std::endl;
    TH1F* theJetQAHist_num= (TH1F*)finMC->Get( theHistName_num.c_str() );

    std::string theHistName_den="hVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName_den<<std::endl;
    TH1F* theJetQAHist_den= (TH1F*)finMC->Get( theHistName_den.c_str() );

    TH1F* theRatio=theJetQAHist_num;
    theRatio->Divide(theJetQAHist_den);
    theRatio->SetTitle("MC Evt QA, MC Only, vz, pthatWeights only");
    theRatio->SetYTitle("weighted/unweighted");
    theRatio->SetXTitle("vz (cm)");
    
    theRatio->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  {
    std::string theHistName_num="hWeightedVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName_num<<std::endl;
    TH1F* theJetQAHist_num= (TH1F*)finMC->Get( theHistName_num.c_str() );

    std::string theHistName_den="hVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName_den<<std::endl;
    TH1F* theJetQAHist_den= (TH1F*)finMC->Get( theHistName_den.c_str() );

    TH1F* theRatio=theJetQAHist_num;
    theRatio->Divide(theJetQAHist_den);
    theRatio->SetTitle("MC Evt QA, MC Only, vz, all Weights");
    theRatio->SetYTitle("weighted/unweighted");
    theRatio->SetXTitle("vz (cm)");
    
    theRatio->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }


  temp_canv->SetLogy(1);



  //pthat
  {
    std::string theHistName="hpthat";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
    theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));

    theJetQAHist->SetTitle("MC event QA, pthat, not weighted");
    theJetQAHist->SetYTitle("Num. Evts/bin");
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

  //vz, MC v Data ratio plots
  { 
    std::string Datafilename=fullDataFilename; 
    TFile *finData = new TFile(Datafilename.c_str());

    std::string theHistName_den="hVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName_den<<std::endl;
    TH1F* theJetQAHist_den= (TH1F*)finData->Get( theHistName_den.c_str() );
    theJetQAHist_den->Scale(1/integratedLuminosity);
    theJetQAHist_den->Scale(1/theJetQAHist_den->GetBinWidth(1));

    std::string theHistName_den2="hWeightedVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName_den<<std::endl;
    TH1F* theJetQAHist_den2= (TH1F*)finData->Get( theHistName_den.c_str() );
    theJetQAHist_den2->Scale(1/integratedLuminosity);
    theJetQAHist_den2->Scale(1/theJetQAHist_den2->GetBinWidth(1));
    
    {//MC, norm, no weights
      
      std::string theHistName_num="hVz";
      if(debugMode)std::cout<<"theHistName="<<theHistName_num<<std::endl;
      TH1F* theJetQAHist_num = (TH1F*)finMC->Get( theHistName_num.c_str() );      
      theJetQAHist_num->Scale(1/theJetQAHist_num->GetBinWidth(1));
      theJetQAHist_num->Scale(theJetQAHist_den->Integral()/theJetQAHist_num->Integral());
      
      TH1F* theRatio=theJetQAHist_num;
      theRatio->Divide(theJetQAHist_den);
      theRatio->SetTitle("evtQA, ratio, MC norm., data");
      theRatio->SetYTitle("MC/Data");
      theRatio->SetXTitle("vz (cm)");
      theRatio->SetAxisRange(0.,2.5,"Y");
      theRatio->SetMarkerColor(kRed+2);
      theRatio->Draw();


      TH1F* theJetQAHist_num2= (TH1F*)theJetQAHist_num->Clone("cloneHist");
      theJetQAHist_num2->Scale(theJetQAHist_den2->Integral()/theJetQAHist_num2->Integral());

      TH1F* theRatio2=theJetQAHist_num2;
      theRatio2->Divide(theJetQAHist_den2);
      theRatio2->SetMarkerColor(kAzure);
      theRatio2->Draw("same");

      TLegend* theLeg = new TLegend(0.1,0.1,0.3,0.3,NULL,"NBNDC");
      theLeg->AddEntry(theRatio,"MC/data w/o trigWghts");
      theLeg->AddEntry(theRatio2,"MC/data w/ trig weights");
      theLeg->Draw();

      temp_canv->Print( thePDFFileName.c_str() );
    }

    {//MC, norm, vz weights
      std::string theHistName_num="hvzWeightedVz";
      if(debugMode)std::cout<<"theHistName="<<theHistName_num<<std::endl;
      TH1F* theJetQAHist_num = (TH1F*)finMC->Get( theHistName_num.c_str() );
      theJetQAHist_num->Scale(1/theJetQAHist_num->GetBinWidth(1));
      theJetQAHist_num->Scale(theJetQAHist_den->Integral()/theJetQAHist_num->Integral());
      
      TH1F* theRatio=theJetQAHist_num;
      theRatio->Divide(theJetQAHist_den);
      theRatio->SetTitle("evtQA, ratio, MC norm. w/ vz weights, data");
      theRatio->SetYTitle("MC/Data");
      theRatio->SetXTitle("vz (cm)");
      theRatio->SetAxisRange(0.,2.5,"Y");
      theRatio->SetMarkerColor(kRed+2);
      theRatio->Draw();


      TH1F* theJetQAHist_num2= (TH1F*)theJetQAHist_num->Clone("cloneHist");
      theJetQAHist_num2->Scale(theJetQAHist_den2->Integral()/theJetQAHist_num2->Integral());

      TH1F* theRatio2=theJetQAHist_num2;
      theRatio2->Divide(theJetQAHist_den2);
      theRatio2->SetMarkerColor(kAzure);
      theRatio2->Draw("same");

      TLegend* theLeg = new TLegend(0.1,0.1,0.3,0.3,NULL,"NBNDC");
      theLeg->AddEntry(theRatio,"MC/data w/o trigWghts");
      theLeg->AddEntry(theRatio2,"MC/data w/ trig weights");
      theLeg->Draw();

      temp_canv->Print( thePDFFileName.c_str() );
    }

    {//MC, norm, pthat weights
      std::string theHistName_num="hpthatWeightedVz";
      if(debugMode)std::cout<<"theHistName="<<theHistName_num<<std::endl;
      TH1F* theJetQAHist_num = (TH1F*)finMC->Get( theHistName_num.c_str() );
      theJetQAHist_num->Scale(1/theJetQAHist_num->GetBinWidth(1));
      theJetQAHist_num->Scale(theJetQAHist_den->Integral()/theJetQAHist_num->Integral());
      
      TH1F* theRatio=theJetQAHist_num;
      theRatio->Divide(theJetQAHist_den);
      theRatio->SetTitle("evtQA, ratio, MC norm. w/ pthat weights, data");
      theRatio->SetYTitle("MC/Data");
      theRatio->SetXTitle("vz (cm)");
      theRatio->SetAxisRange(0.,2.5,"Y");
      theRatio->SetMarkerColor(kRed+2);
      theRatio->Draw();

      TH1F* theJetQAHist_num2= (TH1F*)theJetQAHist_num->Clone("cloneHist");
      theJetQAHist_num2->Scale(theJetQAHist_den2->Integral()/theJetQAHist_num2->Integral());

      TH1F* theRatio2=theJetQAHist_num2;
      theRatio2->Divide(theJetQAHist_den2);
      theRatio2->SetMarkerColor(kAzure);
      theRatio2->Draw("same");

      TLegend* theLeg = new TLegend(0.1,0.1,0.25,0.25,NULL,"NBNDC");
      theLeg->AddEntry(theRatio,"MC/data w/o trigWghts");
      theLeg->AddEntry(theRatio2,"MC/data w/ trig weights");
      theLeg->Draw();

      temp_canv->Print( thePDFFileName.c_str() );
    }

    {//MC, norm, all weights
      std::string theHistName_num="hWeightedVz";
      if(debugMode)std::cout<<"theHistName="<<theHistName_num<<std::endl;
      TH1F* theJetQAHist_num = (TH1F*)finMC->Get( theHistName_num.c_str() );
      theJetQAHist_num->Scale(1/theJetQAHist_num->GetBinWidth(1));
      theJetQAHist_num->Scale(theJetQAHist_den->Integral()/theJetQAHist_num->Integral());
      
      TH1F* theRatio=theJetQAHist_num;
      theRatio->Divide(theJetQAHist_den);
      theRatio->SetTitle("evtQA, ratio, MC norm. w/ all weights, data");
      theRatio->SetYTitle("MC/Data");
      theRatio->SetXTitle("vz (cm)");
      theRatio->SetAxisRange(0.,2.5,"Y");
      theRatio->SetMarkerColor(kRed+2);
      theRatio->Draw();

      TH1F* theJetQAHist_num2= (TH1F*)theJetQAHist_num->Clone("cloneHist");
      theJetQAHist_num2->Scale(theJetQAHist_den2->Integral()/theJetQAHist_num2->Integral());

      TH1F* theRatio2=theJetQAHist_num2;
      theRatio2->Divide(theJetQAHist_den2);
      theRatio2->SetMarkerColor(kAzure);
      theRatio2->Draw("same");

      TLegend* theLeg = new TLegend(0.1,0.1,0.25,0.25,NULL,"NBNDC");
      theLeg->AddEntry(theRatio,"MC/data w/o trigWghts");
      theLeg->AddEntry(theRatio2,"MC/data w/ trig weights");
      theLeg->Draw();

      temp_canv->Print( thePDFFileName.c_str() );
    }

  }
  



  finMC->Close();
  temp_canv->Print( close_thePDFFileName.c_str() );//close the output file
  return 0;
}

