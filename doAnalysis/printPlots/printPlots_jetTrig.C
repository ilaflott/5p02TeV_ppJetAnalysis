#include "printPlots.h"


// code/job switches ------------------------
const bool debugMode=true, doEventCounts=true, doJetIDPlots=true;
const bool drawJetTrigQAPlots=true, drawJetRapBinsPlot=false;

const bool comparePFandCalo=false;
const bool usedHLT100=false;
//const bool usedMinBias=true;
//const bool usedHLTPF=true;//if false, HLTCalo was used
const int N_etabins=4;

// the macro ------------------------
int printPlots_jetTrig ( const std::string input_ppData_condorDir ,   const std::string output_PDFname_base){
  
  globalTrigHistStyle();
  
  //  bool usedHLTPF=false;
  //  if( input_ppData_condorDir.find("HLTPF") != std::string::npos )
  //    usedHLTPF=true;
  
  
  std::string fullJetType, radius;
  makeFullJetTypeString(&fullJetType, &radius, input_ppData_condorDir);
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;







  // OPEN INPUT SECTION
  TFile *finData=NULL;    
  bool usedMinBias;
  if(output_PDFname_base.find("wMB")!=std::string::npos)
    usedMinBias=true;
  else if(output_PDFname_base.find("noMB")!=std::string::npos)
    usedMinBias=false;
  else assert(false);
  
  if(usedMinBias)std::cout<<"running plots for LowJets+Jet80+MinBias PDs"<<std::endl;
  else std::cout<<"running plots for LowJets+Jet80+MinBias PDs"<<std::endl;
  
  std::string input_ppData_Filename;//="HighPtJetTrig_" +fullJetType+ "-allFiles.root";
  if(usedMinBias)input_ppData_Filename="HighPtJetTrig_" +fullJetType+ "-allFiles.root";
  else if (!usedMinBias)input_ppData_Filename="HighPtJetTrig_noMB_" +fullJetType+ "-allFiles.root";
  //input_ppData_Filename="HighPtJet80_" +fullJetType+ "-allFiles.root";
  //input_ppData_Filename="HighPtLowerJets_" +fullJetType+ "-allFiles.root";
  std::cout<<"input_ppData_Filename="<<input_ppData_Filename<<std::endl;
  //assert(false);
  
  const std::string ppData_fullFilename=inputDir+input_ppData_condorDir+input_ppData_Filename;
  
  //open the input files
  std::cout<<" now opening ppData: "<<std::endl<<input_ppData_Filename<<std::endl;
  std::cout<<" in directory: "<<inputDir+input_ppData_condorDir<<std::endl<<std::endl;
  finData = new TFile(ppData_fullFilename.c_str());      
  std::cout<<" ppData_fullFilename = "<<ppData_fullFilename	   <<std::endl;
  if(!finData){
    std::cout << " DATA file not found" << std::endl;
    std::cout << "inputDir               =" << inputDir << std::endl;
    std::cout << "input_ppData_condorDir =" << input_ppData_condorDir << std::endl;
    std::cout << "input_ppData_filename  =" << input_ppData_Filename << std::endl;
    std::cout << "exiting." << std::endl;
    assert(false);
  }
  bool usedHLTPF=false;
  if(input_ppData_condorDir.find("HLTPF")!=std::string::npos || 
     input_ppData_condorDir.find("HLTCalo")!=std::string::npos      ){//older condor output didn't have the TrigComb hist to signify which HLT combo was used. so use the condor directory name instead. 
    if(input_ppData_condorDir.find("HLTPF")!=std::string::npos)
      usedHLTPF=true;
    else
      usedHLTPF=false;
  }
  else{
    TH1F* trigComboName_hist=(TH1F*)finData->Get("TrigComb");
    std::string trigComboName=(std::string)trigComboName_hist->GetTitle();
    if(trigComboName.find("PF")!=std::string::npos)
      usedHLTPF=true;
    else
      usedHLTPF=false;  
    std::cout<<"trgComboName="<<trigComboName<<std::endl;
  }
  
  std::cout<<"usedHLTPF="<<usedHLTPF<<std::endl;
  ///assert(false);
  
  

  // OPEN OPTIONAL INPUT SECTION
  TFile *finData2=NULL;
  std::string input_ppData2_condorDir;
  std::string ppData2_fullFilename;
  
  if(comparePFandCalo){     
    std::cout<<"comparing HLT PF and HLT Calo..."<<std::endl;
    input_ppData2_condorDir=input_ppData_condorDir;
    
    std::cout<<"first, input_ppData2_condorDir="<<input_ppData2_condorDir<<std::endl;
    
    std::string HLTstr, HLTstr_replacement;    
    if(usedHLTPF) { std::cout<<"used HLTPF!"<<std::endl;
      HLTstr="HLTPF";
      HLTstr_replacement="HLTCalo";}
    else { std::cout<<"used HLTCalo!"<<std::endl;
      HLTstr="HLTCalo";
      HLTstr_replacement="HLTPF";}
    
    bool replaceSuccess = str_replace(input_ppData2_condorDir, HLTstr, HLTstr_replacement);    
    std::cout<<"... now, input_ppData2_condorDir="<<input_ppData2_condorDir<<std::endl;
    
    if(replaceSuccess){      //open the input files
      input_ppData2_condorDir="08.28.17_outputCondor/ppData_PFHLT/ppData_HighPtJetTrig_ak4PFJets_08-28-17_jetPlots_0.0eta2.0/";

      std::cout<<"opening file in directory: "<<inputDir+input_ppData2_condorDir<<std::endl<<std::endl;
      ppData2_fullFilename=inputDir+input_ppData2_condorDir+input_ppData_Filename;
      finData2 = new TFile(ppData2_fullFilename.c_str());      
      if(!finData2){
	std::cout << " DATA file not found" << std::endl;
	std::cout << "inputDir               =" << inputDir << std::endl;
	std::cout << "input_ppData_condorDir =" << input_ppData2_condorDir << std::endl;
	std::cout << "input_ppData_filename  =" << input_ppData_Filename << std::endl;
	std::cout << "exiting." << std::endl;
	assert(false); }    
    }
  }
  
  
  
  
  
  
  
  
  
  
  // GET OUTPUT PDF FILE READY
  std::string thePDFFileName=outputDir+fullJetType+"_"+output_PDFname_base+"_jetTrig.pdf";
  std::string open_thePDFFileName=thePDFFileName+"[";    std::string close_thePDFFileName=thePDFFileName+"]";
  std::string theROOTFileName=outputDir+fullJetType+"_"+output_PDFname_base+"_jetTrig.root";
  TFile* fout = new TFile(theROOTFileName.c_str(), "RECREATE");
  std::cout<<"theROOTFileName="<<theROOTFileName<<std::endl;
  // just for opening the pdf
  TCanvas *temp_canvOpen = new TCanvas("temp", "temp", 1200, 1000);
  temp_canvOpen->Print( open_thePDFFileName.c_str() );  
  //temp_canvOpen->UseCurrentStyle();
  temp_canvOpen->Close();
  
  // evtcounts/effective integrated luminosity ----------------------
  long double theLumi;
  if(doEventCounts){
    getEventCounts( (TFile*)finData, true );
    theLumi=computeEffLumi( (TFile*) finData);   }
  else {
    std::cout<<"skipping evt/jet QA counts + plots..."<<std::endl<<std::endl;
    theLumi=intgrtdLumi;}
  
  
  



  
  
  // jet trig plots----------------------
  if(!drawJetTrigQAPlots)
    std::cout<<std::endl<<"skipping jet trig QA plots..."<<std::endl<<std::endl;
  else {     
    
//    // ----------------------------------    
//    printJetTrigHist_wRatio(finData, usedHLT100, true,
//			    thePDFFileName, fullJetType, "excl", radius, usedHLTPF, doJetIDPlots , (TFile*) fout );
//    printJetTrigHist_wRatio(finData, usedHLT100, true,
//			    thePDFFileName, fullJetType, "incl", radius, usedHLTPF, doJetIDPlots , (TFile*) fout);   
//    // ----------------------------------    
//    printTrigPtHist(finData, usedHLT100, true,
//		    thePDFFileName, fullJetType, "excl", radius, usedHLTPF , (TFile*) fout);    
//    printTrigPtHist(finData, usedHLT100, true,
//		    thePDFFileName, fullJetType, "incl", radius, usedHLTPF , (TFile*) fout);
    
    // ----------------------------------    
    printTrigPtHist(finData, usedHLT100, false, usedMinBias, 
    		    thePDFFileName, fullJetType, "excl", radius, usedHLTPF , (TFile*) fout);
    printTrigPtHist(finData, usedHLT100, false, usedMinBias, 
    		    thePDFFileName, fullJetType, "incl", radius, usedHLTPF , (TFile*) fout);    
    // ----------------------------------    
    for(int etabin=0;etabin<N_etabins;etabin++)
      printJetTrigHist_wRatio(finData, usedHLT100, usedMinBias, false, etabin,
			      thePDFFileName, fullJetType, "excl", radius, usedHLTPF, doJetIDPlots , (TFile*) fout);
    for(int etabin=0;etabin<N_etabins;etabin++)
      printJetTrigHist_wRatio(finData, usedHLT100, usedMinBias, false, etabin,
			      thePDFFileName, fullJetType, "incl", radius, usedHLTPF, doJetIDPlots , (TFile*) fout);
    // ----------------------------------    
    //    printTrigEtaHist(finData, usedHLT100,
    //		     thePDFFileName, fullJetType, "excl", radius, usedHLTPF , (TFile*) fout);    
    //    printTrigEtaAsymmHist(finData, usedHLT100,
    //			  thePDFFileName, fullJetType, "excl", radius, usedHLTPF , (TFile*) fout);            
    //    // ----------------------------------    
    //    printTrigEtaHist(finData, usedHLT100, 
    //		     thePDFFileName, fullJetType, "incl", radius, usedHLTPF , (TFile*) fout);    
    //    printTrigEtaAsymmHist(finData, usedHLT100,
    //			  thePDFFileName, fullJetType, "incl", radius, usedHLTPF , (TFile*) fout);    
    //    // ----------------------------------
    
    if(comparePFandCalo){      
      printPFvCaloTrigHist_wRatio(finData, finData2, usedHLT100, usedHLTPF, true, false,
				  thePDFFileName, fullJetType, "excl", radius,  doJetIDPlots);
      printPFvCaloTrigHist_wRatio(finData, finData2, usedHLT100, usedHLTPF, false, false,
				  thePDFFileName, fullJetType, "excl", radius,  doJetIDPlots);      
      printPFvCaloTrigHist_wRatio(finData, finData2, usedHLT100, usedHLTPF, true, false,
				  thePDFFileName, fullJetType, "incl", radius,  doJetIDPlots);
      printPFvCaloTrigHist_wRatio(finData, finData2, usedHLT100, usedHLTPF, false, false, 
				  thePDFFileName, fullJetType, "incl", radius,  doJetIDPlots);
      //combo only comparison
      printPFvCaloTrigHist_wRatio(finData, finData2, usedHLT100, usedHLTPF, true, true,
				  thePDFFileName, fullJetType, "incl", radius,  doJetIDPlots);
      printPFvCaloTrigHist_wRatio(finData, finData2, usedHLT100, usedHLTPF, false, true,
				  thePDFFileName, fullJetType, "incl", radius,  doJetIDPlots);            
    }
  }//end else statement for jet trig plots
  
  
  
  

  
  
  
  //std::string jetCutString, jetEtaCutString, jtptQACut_str;
  //makeJetCutStrings(&jetCutString, &jetEtaCutString, &jtptQACut_str, fin);



  TH1F* jtptCut_h= (TH1F*)finData->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;

  
  TH1F* jtetaLoCut_h= (TH1F*)finData->Get( "hJetEtaCutLo" );
  std::stringstream etaLo; etaLo.precision(1);
  float jtetaLoCut_F=jtetaLoCut_h->GetMean();
  etaLo << std::fixed << jtetaLoCut_F;
  std::string jtetaLoCut_str = etaLo.str();
  std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;
  
  TH1F* jtetaHiCut_h= (TH1F*)finData->Get( "hJetEtaCutHi" );
  std::stringstream etaHi; etaHi.precision(1);
  float jtetaHiCut_F=jtetaHiCut_h->GetMean();
  etaHi << std::fixed << jtetaHiCut_F;
  std::string jtetaHiCut_str = etaHi.str();
  std::cout<<"jtetaHiCut_str = "<<jtetaHiCut_str<<std::endl;
  
  std::string jetPtCutString="p_{T} ^{}>"+jtptCut_str+"  ^{}GeV";
  std::cout<<"jetCutString="<<jetPtCutString<<std::endl;
    
  std::string jetEtaCutString=jtetaLoCut_str+" <  ^{}#||{y} < "+jtetaHiCut_str;
  
  float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  








  
  
  // jets rap-bins plots---------------------- // this is set to "false" right now
  if( !drawJetRapBinsPlot ) 
    std::cout<<std::endl<<"skipping jet rap bin plot..."<<std::endl<<std::endl;
  
  else { std::cout<<std::endl<<std::endl;
    std::cout<<" drawing jet rap bin Plots..."<<std::endl;
    
    if(debugMode)std::cout<<std::endl<<"creating temporary canvas for printing Jet rapBins plots..."<<std::endl;
    TCanvas *temp_canvJetRapBins = new TCanvas("tempJetRapBins", "tempJetRapBins", 1200, 1000);
    temp_canvJetRapBins->SetLogy(1);    
    temp_canvJetRapBins->SetLogx(1);    

    TLegend* theJetRapHistLeg=new TLegend(0.72,0.60,0.85,0.86, NULL,"brNDC");
    
    float integralArray[nbins_abseta]={0.};
    float intCrossSec=0.;

    float integralErrArray[nbins_abseta]={0.};
    float intCrossSecErr=0.;

    temp_canvJetRapBins->cd();
    
    float power=std::pow(10., 6.);
    int powerint = 6;
    
    std::cout<<"drawing dual-differential cross section plot...."<<std::endl;	
    //std::cout<<"looping over "<< nbins_abseta <<" rapidity bins"<<std::endl;

    for(int rapbin=0; rapbin<nbins_abseta; ++rapbin){	
      
      std::stringstream stream1, stream2;
      stream1.precision(1); stream1 << std::fixed << absetabins[rapbin];
      stream2.precision(1); stream2 << std::fixed << absetabins[rapbin+1];
      
      std::string rapHistLegDesc="("+stream1.str()+"<|y|<"+stream2.str()+") x10^{"+std::to_string(powerint)+"}";	
      
      if(debugMode)
	std::cout<<"eta range is..."<< std::endl<< rapHistLegDesc << std::endl<<std::endl;
      
      
      //open the hists + rebin
      std::string inHistName="hJetSpectraRap";	
      if(doJetIDPlots) inHistName+="_wJetID"; 
      inHistName+="_bin"+std::to_string(rapbin);
      
      if(rapbin==6){
	std::cout<<"rapbin==6, continuing on..."<<std::endl;
	continue;
      }
      
      if(debugMode)std::cout<<std::endl<<"inHistName="<<inHistName<<std::endl<<std::endl;	
      
      TH1F* theJetSpectraRapHist= (TH1F*) ( (TH1*)finData->Get(inHistName.c_str()) );
      theJetSpectraRapHist=(TH1F*)theJetSpectraRapHist->TH1::Rebin(nbins_pt_debug, "ppData_raphist_jetTrig", ptbins_debug);
      theJetSpectraRapHist->Print( "base" );
      
      //TH1F* theMCJetSpectraRapHist= (TH1F*) ( (TH1*)finMC->Get(inHistName.c_str()) );
      //theMCJetSpectraRapHist=(TH1F*)theMCJetSpectraRapHist->TH1::Rebin(nbins_pt_debug, "ppMC_raphist_jetTrig", ptbins_debug);
      //theMCJetSpectraRapHist->Print( "base" );
      

      //// necessary scaling
      divideBinWidth(theJetSpectraRapHist);
      theJetSpectraRapHist->Scale( 1./theLumi);
      if(rapbin==7)theJetSpectraRapHist->Scale( 1./3. );

      //// integral 
      integralArray[rapbin]=getIntegralJetSpectra(theJetSpectraRapHist);
      if(rapbin==7)integralArray[rapbin]*=3.;
      integralArray[rapbin]/=power;
      
      intCrossSec+=integralArray[rapbin];      
      
      integralErrArray[rapbin]=getIntegralErrJetSpectra(theJetSpectraRapHist);
      if(rapbin==7)integralErrArray[rapbin]*=3.;
      integralErrArray[rapbin]/=power;

      intCrossSecErr+=integralErrArray[rapbin]*integralErrArray[rapbin];



      theJetSpectraRapHist->Scale( power );


      //// necessary scaling
      //divideBinWidth(theMCJetSpectraRapHist);
      //if(rapbin==7)theMCJetSpectraRapHist->Scale( 1./3. );
      //theMCJetSpectraRapHist->Scale( power );      
      //theMCJetSpectraRapHist->Scale(1000.);



      //// cosmetics
      theJetSpectraRapHist->SetMarkerStyle( kCircle ) ;
      theJetSpectraRapHist->SetMarkerSize(0.90);
      theJetSpectraRapHist->SetMarkerColor( theRapOverlayMarkerColor[rapbin] );
      theJetSpectraRapHist->SetLineColor( kBlack );

      //theMCJetSpectraRapHist->SetMarkerStyle( kMultiply ) ;
      //theMCJetSpectraRapHist->SetMarkerSize(0.90);
      //theMCJetSpectraRapHist->SetMarkerColor( theRapOverlayMarkerColor[rapbin] );
      //theMCJetSpectraRapHist->SetLineColor( kBlack );
      
      theJetRapHistLeg->AddEntry(theJetSpectraRapHist,(rapHistLegDesc).c_str(),"lp");
      
      if(rapbin==0){
	
	std::string h_Title   ="Jet Pt in |y| bins";
	if(doJetIDPlots) h_Title+=", w/ JetIDCut";
	theJetSpectraRapHist->SetTitle (    h_Title.c_str() );
	
	std::string h_XAx_Title="jet p_{T} (GeV)    ", h_YAx_Title=ddcrossSectionAxTitle;
	theJetSpectraRapHist->SetXTitle( h_XAx_Title.c_str() );
	theJetSpectraRapHist->SetYTitle( h_YAx_Title.c_str() );
	
	theJetSpectraRapHist->GetYaxis()->SetLabelSize(0.02);
	theJetSpectraRapHist->GetYaxis()->SetTitleSize(0.025);
	theJetSpectraRapHist->GetYaxis()->SetTitleOffset(1.2);
	
	theJetSpectraRapHist->GetXaxis()->SetLabelSize(0.02);
	theJetSpectraRapHist->GetXaxis()->SetTitleSize(0.025);
	theJetSpectraRapHist->GetXaxis()->SetTitleOffset(.8);
	
	theJetSpectraRapHist->SetAxisRange(std::pow(10.,-9.3),std::pow(10.,5.3),"Y");
	
	theJetSpectraRapHist->Draw("E");     	  
	//theMCJetSpectraRapHist->Draw("E SAME");
	theJetRapHistLeg->Draw(); }
      else {
	theJetSpectraRapHist->Draw("E SAME"); }
	//	theMCJetSpectraRapHist->Draw("E SAME");       }
      
      power/=10.;
      powerint--;
      
    }//end rapbin hist loop
    

    float t1Loc1=0.47, t1Loc2=0.85; 	    
    TLatex *t0=new TLatex(t1Loc1,t1Loc2+.02,(Py8TuneString + " (crosses)").c_str());
    t0->SetTextFont(63);
    t0->SetTextColor(kBlack);
    t0->SetTextSize(12);
    t0->SetLineWidth(1);
    t0->SetNDC();
    t0->Draw();

    //float t1Loc1=0.47, t1Loc2=0.85; 	    
    TLatex *t1=new TLatex(t1Loc1,t1Loc2,(sqrts2k15ppString + " (open circles)").c_str());
    t1->SetTextFont(63);
    t1->SetTextColor(kBlack);
    t1->SetTextSize(12);
    t1->SetLineWidth(1);
    t1->SetNDC();
    t1->Draw();
    
    TLatex *pp1=new TLatex(t1Loc1,(t1Loc2-0.02),intLumiString.c_str());
    pp1->SetTextFont(63);
    pp1->SetTextColor(kBlack);
    pp1->SetTextSize(12);
    pp1->SetLineWidth(1);
    pp1->SetNDC();
    pp1->Draw();
    
    TLatex *pp2=new TLatex(t1Loc1,(t1Loc2-0.04),(jetPtCutString + ", "+fullJetType+"jets").c_str());
    pp2->SetTextFont(63);
    pp2->SetTextColor(kBlack);
    pp2->SetTextSize(12);
    pp2->SetLineWidth(1);
    pp2->SetNDC();
    pp2->Draw();
    
    
    temp_canvJetRapBins->Print( thePDFFileName.c_str() );
    
    
    
    std::cout<<std::endl;    std::cout<<std::endl;
    
    for(int i =0; i<nbins_abseta;i++){
      continue;
      if(i==6) continue;
      std::cout<<"for eta bin i="<<i<<std::endl;
      std::cout<<"integral="<<integralArray[i]<<" +/- "<< integralErrArray[i]<<" microbarns"<<std::endl<<std::endl;      
      if(i==7){
	std::cout<<std::endl<<std::endl<<std::endl;
	intCrossSecErr=std::sqrt(intCrossSecErr);
	std::cout<<"integrated cross section ~= "<<std::endl<<
	  intCrossSec<<" +/- "<< intCrossSecErr <<" microbarns"<<std::endl;
	std::cout<<"%err = "<< (((intCrossSecErr/intCrossSec)*100.)) << std::endl;      }
      
    }
  
  }
  





  
  
  
  
  if(debugMode)std::cout<<std::endl<<"closing the PDF file"<<std::endl;
  TCanvas *temp_canvClose = new TCanvas("tempClose", "tempClose", 1200, 1000);
  temp_canvClose->Print( close_thePDFFileName.c_str() );  
  temp_canvClose->Close();    
  

  if(debugMode)std::cout<<std::endl<<"closing input files"<<std::endl;
  //finMC->Close();
  finData->Close();  
  fout->Close();
  //assert(false);
  
  return 0;
}// end printplots




// steering ----------------------

int main(int argc, char *argv[]){
  
  int rStatus = -1;
  if( argc!=3 ) {//no input arguments, error
    std::cout<<"do ./printPlots_jetTrig.exe <ppDataCondorDir> <outputNameBase>"<<std::endl;
    return rStatus;
  }
  rStatus=1;
  
  
  if(argc==3) rStatus=printPlots_jetTrig( (const std::string) argv[1], (const std::string) argv[2]);
  
  std::cout<<"done, rStatus="<<rStatus<<std::endl;
  return rStatus;
  
}



///////////////////////
//// PROBLEM SECTION, do NOT delete comments or other commented code
//if(j==N_trigs)std::cout<<"warning! j=="<< N_trigs<<"=N_trigs, about to seg fault..."<<std::endl;// NEVER FIRES
//
//if(debugMode)std::cout<<std::endl<<"N_trigs="<<N_trigs<<" and j="<<j<<std::endl;
//if(debugMode)std::cout<<"j<(N_trigs-1)="<<(bool)(j<(N_trigs-1))<<std::endl;
//if(debugMode)std::cout<<"j<N_trigs="<<(bool)(j<(N_trigs))<<std::endl;
//if(debugMode)std::cout<<"j<(N_trigs+1)="<<(bool)(j<(N_trigs+1))<<std::endl;	
//
////assert(j!=N_trigs);// NEVER FIRES
////assert(j<N_trigs); // NEVER FIRES
//
////open the hists + do scaling 
//if(j==N_trigs)std::cout<<"warning! j=="<< N_trigs<<"=N_trigs, about to seg fault..."<<std::endl;// NEVER FIRES
//
//if(j==N_trigs){  // uncommenting this if/else statement fixes the pathological behavior
//  std::cout<<"warning! j==N_trigs, about to seg fault..1."<<std::endl;
//  continue;}  
//else std::cout<<" HLTName ="<<HLTName[j]<<std::endl;
//
////if(j==N_trigs){  // uncommenting this if/else statement fixes the pathological behavior
////  std::cout<<"warning! j==N_trigs, about to seg fault..1."<<std::endl;
////  continue;}	
////else std::cout<<" HLTName ="<<HLTName[j]<<std::endl;
//  
////open the hists + do scaling 
//if(j==N_trigs)std::cout<<"warning! j=="<< N_trigs<<"=N_trigs, about to seg fault..."<<std::endl;// NEVER FIRES
//// PROBLEM SECTION
///////////////////////


// old -------------------------------------
//const bool drawEvtQAPlots=true;
//const bool drawJetQAPlots=true;
//const bool drawJetConstituentPlots=true, drawDijetPlots=true;
//const bool drawDataMCOverlaysInput=false; 		 //this should always be false, but needs to be cleaned up later.
//const bool drawDataMCOverlays   = drawDataMCOverlaysInput;
//const bool drawDataMCRatios     = true;   //!drawDataMCOverlays;
//const bool drawTrigataOverlays = false;  //drawDataMCOverlays && drawJetTrigQAPlots;
//const bool drawTrigDataRatios   = false; //drawDataMCRatios   && drawJetTrigQAPlots;
//const int theDataOverlayLineColor=1, theMCOverlayLineColor=1;//, altOverlayLineColor=3; 
//const int theRatioLineColor=1;//,altRatioLineColor1=8, altRatioLineColor2=7;
//const int theDataOverlayMarkerColor=2, theMCOverlayMarkerColor=4;//,theRatioMarkerColor=9;
//const int theTrigCombMarkerColor=1, altTrigCombMarkerColor=12;
  
 



     //else if(drawTrigDataRatios){ temp_canvJetRatio->cd();
     //  
     //  TLegend* JetTrigRatioLegend=new TLegend(0.12,0.59,0.30,0.87, NULL,"brNDC");    
     //  
     //  if(debugMode)std::cout<<"drawing dataTrig QA Plots"<<std::endl;
     //  
     //  std::string comboPlotName="hpp_"+HLTName[N_trigs-1];
     //  if(doJetIDPlots)comboPlotName+="_wJetID";
     //  comboPlotName+="_"+radius+etaWidth;
     //  TH1F* theComboHist=(TH1F*)
     //	 ( (TH1*)finData->Get(comboPlotName.c_str()) )
     //	 ->TH1::Rebin(jetTrigQABinning/2);
     //  
     //  for(int j=0; j<(N_trigs-1); j++){ 
     //	 
     //	 std::string inHistName="hpp_"+HLTName[j];
     //	 if(doJetIDPlots)inHistName+="_wJetID";
     //	 inHistName+="_"+radius+etaWidth;	
     //	 TH1F* theRatio=(TH1F*)
     //	   ( (TH1*)finData->Get(inHistName.c_str()) )
     //	   ->TH1::Rebin(jetTrigQABinning/2);
     //	 
     //	 theRatio->SetMarkerStyle(theTrigOverlayMarker[j]);
     //	 theRatio->SetMarkerSize(0.90);
     //	 theRatio->SetMarkerColor( theTrigOverlayMarkerColor[j] );
     //	 theRatio->SetLineColor( theTrigOverlayLineColor[j] );	
     //	 
     //	 JetTrigRatioLegend->AddEntry(theRatio,(HLTName_Leg[j]+"/TrigComb.").c_str(),"lp");
     //	 
     //	 TH1F* comboHistClone=(TH1F*)theComboHist->Clone("comboHistClone4Div");
     //	 theRatio->Divide(comboHistClone);
     //	 
     //	 if(j==0){
     //	   std::string h_Title   ="JetTrigQA Ratios";
     //	   if(doJetIDPlots)    h_Title+=", w/ JetIDCut";      
     //	   theRatio->SetTitle (    h_Title.c_str() );
     //	   
     //	   std::string h_XAx_Title="reco JetPt (GeV)", h_YAx_Title="Trig/TrigCombo"    ; 
     //	   theRatio->SetAxisRange(0.,1000.,"X");
     //	   theRatio->SetAxisRange(0.,2.,"Y");
     //	   theRatio->SetXTitle( h_XAx_Title.c_str() );
     //	   theRatio->SetYTitle( h_YAx_Title.c_str() );
     //	   theRatio->Draw("E"); }	
     //	 else 
     //	   theRatio->Draw("E SAME");
     //  }     
     //  
     //  JetTrigRatioLegend->Draw();
     //  temp_canvJetRatio->Print(thePDFFileName.c_str());
     //}//end DRAW TRIG Ratios       

      //TH1F* theMCJetSpectraRapHist= (TH1F*)
      //  ( (TH1*)finMC->Get(inMCHistName.c_str()) )
      //  ->TH1::Rebin(jetSpectraRapBinning);
      //
      //theMCJetSpectraRapHist->Scale( 1./theMCJetSpectraRapHist->GetBinWidth(1) );
      //theMCJetSpectraRapHist->Scale( power );
      //theMCJetSpectraRapHist->Scale( theJetSpectraRapHist->Integral()/theMCJetSpectraRapHist->Integral() );
      //
      //theMCJetSpectraRapHist->SetMarkerStyle( 27 ) ;//theRapOverlayMarker[j]);
      //theMCJetSpectraRapHist->SetMarkerSize(0.90);
      //theMCJetSpectraRapHist->SetMarkerColor( theRapOverlayMarkerColor[rapbin] );
      //theMCJetSpectraRapHist->SetLineColor( -1 );//theRapOverlayLineColor[j]   );
      ////theJetSpectraRapHist->SetAxisRange(0.,800.,"X");	 
