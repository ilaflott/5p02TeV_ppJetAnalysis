#include "printQAPlots.h"

//input files
const std::string input_ppData_Filename=
  //"/old/readFiles_ppData_5p02TeV_HighPtJetTrig_08-23-16__vzWithAndWithoutWeights/HighPtJetTrig_ak4PF-allFiles.root";
  "readForests_ppData_5p02TeV_HighPtJetTrig_08-29-16__wEvtJetCounts_AjxjBugFix/HighPtJetTrig_ak4PF-allFiles.root";
const std::string ppData_fullFilename=CMSSW_BASE+input_ppData_dir+input_ppData_Filename;

const std::string input_ppMC_Filename=
  //"readFiles_ppMC_5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_08-23-16__RaghavsVzWeights_noTrigWeights/QCDjetAllPtBins_ak4PF-allFiles.root";
  //"/old/readFiles_ppMC_5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_08-23-16__AnnasVzWeights_noTrigWeights/QCDjetAllPtBins_ak4PF-allFiles.root";
  //"readForests_ppMC_5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_08-29-16__noJECorJER_wmceff_wEvtJetCounts_AjxjBugFix/Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root";
  "readForests_ppMC_5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_08-29-16__divByVzWeightTry/Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root";
const std::string ppMC_fullFilename=CMSSW_BASE+input_ppMC_dir+input_ppMC_Filename;


// output pdfs
const std::string outputDir=CMSSW_BASE+"/src/printQAPlots/output/";
const std::string output_ppMC_PDFname =outputDir+
  "printQAPlots_ppMC_PY8_TuneCUETP8M1-QCDDijet_QA-jets-evts-weights_8.29.16.pdf";
const std::string output_ppData_PDFname=outputDir+
  "printQAPlots_ppData_HighPtJets_PromptReco2015_QA-jets-evts-trig_8.29.16.pdf";
const std::string output_ppDataMCOverlaid_PDFname=outputDir+
  "printQAPlots_pp_DataMCOverlaid_HighPtJets_PromptReco2015_PY8-QCD_QA-jets-evts-trig_8.29.16.pdf";

//debug switches
const bool debugMode=true, dataDebugMode=true, MCDebugMode=true;


//other switches
const bool drawEvtQAPlots=true;
//const bool drawJetQAPlots=true;
//const bool drawTrigPlots=true;
//const bool drawJECandJERPlots=true;
//const bool doEventCounts=true;

const bool drawDataMCOverlays=true;
const int theDataOverlayMarkerColor=1;//black
const int theMCOverlayMarkerColor=1;
const int theDataOverlayLineColor=2;//red
const int theMCOverlayLineColor=4;//blue
const int optThirdOverlayLineColor=3;//green

const bool drawDataMCRatios=true;

const bool openMCFile=drawDataMCOverlays;
const bool openDataFile=drawDataMCOverlays;

int printQAPlots(){
 
  std::cout<<std::endl<<"printing QA Plots, now opening input files!!"<<std::endl<<std::endl;

  //input file(s)
  TFile *finMC=NULL, *finData=NULL;
  assert( openMCFile||openDataFile );

  if(openDataFile&&debugMode)  std::cout<<" now opening Data File "<<std::endl<<ppData_fullFilename<<std::endl<<std::endl;
  if(openDataFile) finData = new TFile(ppData_fullFilename.c_str());

  if(openMCFile&&debugMode)  std::cout<<" now opening MC File "<<std::endl<<ppMC_fullFilename<<std::endl<<std::endl;
  if(openMCFile) finMC = new TFile(ppMC_fullFilename.c_str());
  /*
  // evt counts----------------------
  if(!doEventCounts)std::cout<<"skipping evt QA plots..."<<std::endl<<std::endl;
  else{ std::cout<<" grabbing raw event counts "<<std::endl;

    {
      std::cout<<std::endl<<"/// MC Evt-Count Summary ///"<<std::endl<<std::endl;
      std::cout<<"Total Num of Events in file(s) opened       = " <<h_NEvents->GetEntries()<<std::endl;
      std::cout<<"Total Num of Events read from those file(s) = " <<h_NEvents_read->GetEntries()<<std::endl;
      
      if(debugMode){
	std::cout<<"Total Num of Events read passing skimCuts                              = " <<
	  h_NEvents_skimCut->GetEntries()<<std::endl;
	std::cout<<"Total Num of Events read passing vzCuts and skimCuts (aka good events) = " <<
	  h_NEvents_vzCut->GetEntries()<<std::endl;
	std::cout<<std::endl;
	std::cout<<"Total Num of good Events w/ jets                             = " <<
	  h_NEvents_withJets->GetEntries()<<std::endl;
	std::cout<<"Total Num of good Events, w/ jets, post kmatCuts             = " <<
	  h_NEvents_withJets_kmatCut->GetEntries()<<std::endl;
	std::cout<<"Total Num of good Events, w/ jets, post kmatCut AND JetIDCut = " <<
	  h_NEvents_withJets_JetIDCut->GetEntries()<<std::endl;      }
      
      std::cout<<std::endl<<"/// Job Jet-Loop Summary ///"<<std::endl<<std::endl;
      std::cout<<"Total Num of Jets read from good events                          = " <<
	h_NJets->GetEntries()<<std::endl;
      std::cout<<"Total Num of Jets read from good events post kinCut              = " <<
	h_NJets_kmatCut->GetEntries()<<std::endl;
      std::cout<<"Total Num of Jets read from good events post kinCut AND JetIDCut = " <<
	h_NJets_JetIDCut->GetEntries()<<std::endl;
      std::cout<<std::endl;
    }

    {
      std::cout<<std::endl<<"/// MC Evt-Count Summary Summary ///"<<std::endl<<std::endl;
      std::cout<<"Total Num of Events in file(s) opened       = " <<h_NEvents->GetEntries()<<std::endl;
      std::cout<<"Total Num of Events read from those file(s) = " <<h_NEvents_read->GetEntries()<<std::endl;
      
      if(debugMode){
	std::cout<<"Total Num of Events read passing skimCuts                              = " <<
	  h_NEvents_skimCut->GetEntries()<<std::endl;
	std::cout<<"Total Num of Events read passing vzCuts and skimCuts (aka good events) = " <<
	  h_NEvents_vzCut->GetEntries()<<std::endl;
	std::cout<<std::endl;
	std::cout<<"Total Num of good Events w/ jets                             = " <<
	  h_NEvents_withJets->GetEntries()<<std::endl;
	std::cout<<"Total Num of good Events, w/ jets, post kmatCuts             = " <<
	  h_NEvents_withJets_kmatCut->GetEntries()<<std::endl;
	std::cout<<"Total Num of good Events, w/ jets, post kmatCut AND JetIDCut = " <<
	  h_NEvents_withJets_JetIDCut->GetEntries()<<std::endl;      }
      
      std::cout<<std::endl<<"/// MC Jet-Count Summary Summary ///"<<std::endl<<std::endl;
      std::cout<<"Total Num of Jets read from good events                          = " <<
	h_NJets->GetEntries()<<std::endl;
      std::cout<<"Total Num of Jets read from good events post kinCut              = " <<
	h_NJets_kmatCut->GetEntries()<<std::endl;
      std::cout<<"Total Num of Jets read from good events post kinCut AND JetIDCut = " <<
	h_NJets_JetIDCut->GetEntries()<<std::endl;
      
      std::cout<<std::endl;
    }

    
  }// end evt counts
  */



  //get output file ready
  std::string thePDFFileName="";
  if(drawDataMCOverlays)thePDFFileName+=output_ppDataMCOverlaid_PDFname;  
  std::string open_thePDFFileName=thePDFFileName+"[";
  std::string close_thePDFFileName=thePDFFileName+"]";

  if(debugMode)std::cout<<"creating temporary canvas for printing..."<<std::endl;
  TCanvas *temp_canv = new TCanvas("temp", "temp", 1200, 600);

  if(debugMode)std::cout<<"printing temp canvas to open PDF file..."<<std::endl;
  temp_canv->Print( open_thePDFFileName.c_str() );      
  temp_canv->cd();




  // evt plots ----------------------
  if(!drawEvtQAPlots)std::cout<<"skipping evt QA plots..."<<std::endl<<std::endl;
  else {    std::cout<<" drawing evt QA Plots now! "<<std::endl<<std::endl;

    {//MC and data total weighted vz
      std::string inHistName="hWeightedVz";

      if(drawDataMCOverlays){   if(debugMode)std::cout<<" drawing DataMC Overlay for hist "<<inHistName<<std::endl;

	std::string h_Title   ="Data+MC Evt QA, vz w/ all weights";
	std::string h_XAx_Title="vz (cm)"                          ;
	std::string h_YAx_Title="millibarns/bin"                   ;
	
	TH1F *theMCEvtQAHist, *theDataEvtQAHist;
	theDataEvtQAHist= (TH1F*)finData->Get( inHistName.c_str() );
	theDataEvtQAHist->Scale( 1/theDataEvtQAHist->GetBinWidth(0) );
	theDataEvtQAHist->Scale( 1/intgrtdLumi );

	theDataEvtQAHist->SetMarkerStyle(kDot);
	theDataEvtQAHist->SetMarkerSize(1.1);
	theDataEvtQAHist->SetMarkerColor( theDataOverlayMarkerColor );
	theDataEvtQAHist->SetLineColor( theDataOverlayLineColor );

	theDataEvtQAHist->SetTitle (    h_Title.c_str() );
	theDataEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
	theDataEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
	theDataEvtQAHist->Draw();
	
	theMCEvtQAHist= (TH1F*)finMC->Get( inHistName.c_str() );
	theMCEvtQAHist->Scale( 1/theMCEvtQAHist->GetBinWidth(0) );
	theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );
	
	theMCEvtQAHist->SetMarkerStyle(kDot);
	theMCEvtQAHist->SetMarkerSize(1.1);
	theMCEvtQAHist->SetMarkerColor( theMCOverlayMarkerColor );
	theMCEvtQAHist->SetLineColor( theMCOverlayLineColor );
	theMCEvtQAHist->Draw("SAME");

	TLegend* theEvtQALeg=new TLegend(0.10,0.70,0.30,0.85, NULL,"brNDC");
	theEvtQALeg->AddEntry(theDataEvtQAHist,"HighPt Jet Trig Data","lp");
	theEvtQALeg->AddEntry(theMCEvtQAHist,"Py8 QCD Dijet","lp");
	
	{//pthat weighted only MC
	  TH1F* theMCEvtQAHist2= (TH1F*)finMC->Get( "hpthatWeightedVz" );
	  theMCEvtQAHist2->Scale( 1/theMCEvtQAHist2->GetBinWidth(0) );
	  theMCEvtQAHist2->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist2->Integral() );
	  
	  theMCEvtQAHist2->SetMarkerStyle(kDot);
	  theMCEvtQAHist2->SetMarkerSize(1.1);
	  theMCEvtQAHist2->SetMarkerColor( theMCOverlayMarkerColor );
	  theMCEvtQAHist2->SetLineColor( 3 );
	  theMCEvtQAHist2->Draw("SAME");

	  theEvtQALeg->AddEntry(theMCEvtQAHist2,"MC w/o vz weights","lp");
	}

	theEvtQALeg->Draw();
	temp_canv->Print( thePDFFileName.c_str() );	}    

      if(drawDataMCRatios){	
	std::string h_Title   ="DataMC Ratio Evt QA, vz w/ all weights";
	std::string h_XAx_Title="vz (cm)"                          ;
	std::string h_YAx_Title="Data/MC"                   ;
	
	TH1F *theMCEvtQAHist, *theDataEvtQAHist;
	theDataEvtQAHist= (TH1F*)finData->Get( inHistName.c_str() );
	theDataEvtQAHist->Scale( 1/theDataEvtQAHist->GetBinWidth(0) );
	theDataEvtQAHist->Scale( 1/intgrtdLumi );

	theDataEvtQAHist->SetTitle (    h_Title.c_str() );
	theDataEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
	theDataEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
	
	theMCEvtQAHist= (TH1F*)finMC->Get( inHistName.c_str() );
	theMCEvtQAHist->Scale( 1/theMCEvtQAHist->GetBinWidth(0) );
	theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );
	
	TLegend* theEvtQALeg=new TLegend(0.10,0.70,0.30,0.85, NULL,"brNDC");
	
	//pthat weighted only Data/MC ratio
	{ TH1F *theMCEvtQAHist2, *theDataEvtQAHist2;
	  theDataEvtQAHist2=(TH1F*)theDataEvtQAHist->Clone("DataVzClone");
	  
	  theDataEvtQAHist2->SetTitle (    h_Title.c_str() );
	  theDataEvtQAHist2->SetXTitle( h_XAx_Title.c_str() );
	  theDataEvtQAHist2->SetYTitle( h_YAx_Title.c_str() );
	  
	  //pthat weighted only MC
	  theMCEvtQAHist2= (TH1F*)finMC->Get( "hpthatWeightedVz" );
	  theMCEvtQAHist2->Scale( 1/theMCEvtQAHist2->GetBinWidth(0) );
	  theMCEvtQAHist2->Scale( theDataEvtQAHist2->Integral()/theMCEvtQAHist2->Integral() );
	  
	  TH1F *theRatio2=theDataEvtQAHist2;
	  theRatio2->Divide(theMCEvtQAHist2);
	  
	  theRatio2->SetMarkerStyle(kDot);
	  theRatio2->SetMarkerSize(1.1);
	  theRatio2->SetMarkerColor( kBlack );
	  theRatio2->SetLineColor( 8 );
	  theRatio2->SetAxisRange(0,1.5,"Y");
	  
	  theRatio2->Draw("");
	  
	  theEvtQALeg->AddEntry(theRatio2,"Data/MC w/o vz weights","lp"); }
	
	TH1F *theRatio=theDataEvtQAHist;
	theRatio->Divide(theMCEvtQAHist);
	
	theRatio->SetMarkerStyle(kDot);
	theRatio->SetMarkerSize(1.1);
	theRatio->SetMarkerColor( kBlack );
	theRatio->SetLineColor( 9 );
	theRatio->SetAxisRange(0,1.5,"Y");
	
	theRatio->Draw("SAME");


	theEvtQALeg->AddEntry(theRatio,"Data/MC","lp");
	theEvtQALeg->Draw();	

	temp_canv->Print( thePDFFileName.c_str() );	}    
    }//end total weighted vz
    


      
    /*
    //else if(drawMCPlots)
    // MC Only evt hists
    {//unweighted Vz
    std::string theHistName="hVz";
    
    TH1F* theEvtQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
    
    theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));
    theJetQAHist->SetTitle("MC evt QA, vz w/o weights");
    theJetQAHist->SetYTitle("Num. MC Evts/bin");
    theJetQAHist->SetXTitle("vz (cm)");
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
    }
    {//vzweighted vz
    std::string theHistName="hvzWeightedVz";
    
    TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
    theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));
    theJetQAHist->SetTitle("MC event QA, vz w/ vz weights only");
    theJetQAHist->SetYTitle("millibarns/bin");
    theJetQAHist->SetXTitle("vz (cm)");
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
    }
    {//pthat
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
    {//pthat weighted pthat
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

    //Data Only Hists
    //else if(drawDataPlots)
    {//Data unweighted Vz
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
    //----------------------
    */  

  }  //end evt plots
  


  
  
  
  
  
  /*
  // jet plots----------------------

  // MC JetQA 
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


  // Data jetQA
  for(int j=0;j<N_vars;j++){ for(int i=0;i<2;i++){      
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
  */
  /*
  // ratio plots----------------------
  {   //vz, MC v Data
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

    }*/
  


  /*
  // trig plots ----------------------
  temp_canv->SetLogy(1);
  // Data Only, HLT/L1/TrigComb Plots
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
  */

  if(debugMode)std::cout<<std::endl<<"closing input files"<<std::endl;
  finMC->Close();
  finData->Close();

  if(debugMode)std::cout<<std::endl<<"closing the PDF file"<<std::endl<<std::endl;
  temp_canv->Print( close_thePDFFileName.c_str() );//close the output file

  return 0;
}













// steering ----------------------

int main(int argc, char *argv[]){

  int rStatus = -1;
  if( argc!=1 ) {//no input arguments, error
    std::cout<<"debugging, for now just do ./printQAPlots.exe"<<std::endl;
    return rStatus;
  }
  rStatus=1;
  rStatus=printQAPlots();
  std::cout<<"done, rStatus="<<rStatus<<std::endl;
  return rStatus;
}












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





  ////data only ratio, vz
  //{
  //  std::string theHistName_num="hWeightedVz";
  //  if(debugMode)std::cout<<"theHistName="<<theHistName_num<<std::endl;
  //  TH1F* theJetQAHist_num= (TH1F*)fin->Get( theHistName_num.c_str() );
  //
  //  std::string theHistName_den="hVz";
  //  if(debugMode)std::cout<<"theHistName="<<theHistName_den<<std::endl;
  //  TH1F* theJetQAHist_den= (TH1F*)fin->Get( theHistName_den.c_str() );
  //
  //  TH1F* theRatio=theJetQAHist_num;
  //  theRatio->Divide(theJetQAHist_den);
  //  theRatio->SetTitle("data Event QA, vz w/ over vz w/o TrigWeight");
  //  theRatio->SetYTitle("weighted/unweighted");
  //  theRatio->SetXTitle("vz (cm)");
  //  theRatio->SetAxisRange(0., 18., "Y");
  //
  //  temp_canv->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}
  ////----------------------
