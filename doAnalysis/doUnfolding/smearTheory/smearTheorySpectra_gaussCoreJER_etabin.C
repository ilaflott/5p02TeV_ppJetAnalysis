#include "smearTheorySpectra.h"

//const std::string SPLINE_STR="spline";
//const std::string FIT_STR="fit";
const int CANVX=1200, CANVY=1000;
const int xsecorder=0;//0--> LO, 1--> NLO, 2--> NNLO... etc.

const bool printBaseDebug=true;

//const int nEvents=1e+09;  ///10x typical
const int nEvents=1e+08;  ///typical
//const int nEvents=1e+07;  /// debug nevents
//const int nEvents=1e+06;  /// debug nevents
//const int nEvents=1e+05;  /// debug nevents

std::string ddxsec_yax="#frac{d^{2}#sigma}{dp_{T}dy} [nb/GeV]";

const std::string NPCorrFile="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory/NPCorr5TeV/NLOpNP_InclusiveJets5TeV.root";
const std::string in_JERFile_default="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots_JERS/output/ak4PF_MCJEC_00eta05.root";
const std::string in_NLOFile_default="fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC";

const bool useHistSigmaFit=true;
int debugInt=0;
void debugcout(){
  std::cout<<"location #"<<debugInt<<std::endl;
  debugInt++;
  return; 
}
//, std::string weightMode="spline" //"fit"// 
//valid fitTypes; modLog, 7TeV, spl3wgts, spl3wgts + ext string
int smearTheorySpectra_gaussCoreJER_etabin( std::string in_NLOfileString=in_NLOFile_default, std::string fitType="modLog" , 
					    std::string in_JERFile=in_JERFile_default,  int etabin=0,  std::string outputFile="smearTheory_gaussCoreJER_etabin_defOut.root"
					    ){
  
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);
 
  TVirtualFitter::SetDefaultFitter("Minuit2");
  //ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-02);      //default i think
  //ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-03);      
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04);      
  
  // Output File to write to
  const bool useSplineWeights=((bool)(fitType.find("spl3wgts")!=std::string::npos));
  const bool useSplineExt=(useSplineWeights && ((bool)(fitType.find("ext")!=std::string::npos)));
  const bool useFitWeights=!(useSplineWeights);
  const bool useModLogFit=useFitWeights && ((bool)(fitType=="modLog"));
  const bool use7TeVFit=useFitWeights && ((bool)(fitType=="7TeV"));
  if(!useModLogFit && !use7TeVFit && !useSplineWeights){
    std::cout<<"fit type not recognized! Exit!"<<std::endl; assert(false);}
  
  std::string outputTag="";
  if(etabin==0) outputTag="00eta05";
  else if(etabin==1 ) outputTag="05eta10";
  else if(etabin==2 ) outputTag="10eta15";
  else if(etabin==3 ) outputTag="15eta20";
  else if(etabin==4 ) outputTag="20eta25";
  else if(etabin==5 ) outputTag="25eta30";
  
  std::cout<<"opening output:"<<outputFile<<std::endl<<std::endl;
  TFile *outf    = new TFile(outputFile.c_str(), "RECREATE" );
  
  // Input, Theory spectrum
  std::string inputFile = "fNLOJetsSpectra/R04/"+in_NLOfileString+".root";
  std::cout<<"opening input file:"<<inputFile<<std::endl<<std::endl;
  TFile* fin_NLO=TFile::Open(inputFile.c_str());
  
  
  /////////////// NPs to NLO jet spectra ////////////////////////////////////////////////////////  
  std::cout<<"opening NP corr file " << NPCorrFile<<std::endl;
  TFile* fin_NP=TFile::Open(NPCorrFile.c_str());  
  std::string NPCorrFits_str;
  if(NPCorrFile.find("nnlo")!=std::string::npos)
    NPCorrFits_str="fNPC_POWPY8_R4_etabin";
  else//raghavs suggestion to use POWPY8 for NNLO thy.
    NPCorrFits_str="fNPC_HerwigEE4C_R4_etabin";
  std::cout<<"NPCorrFits_str = "<<NPCorrFits_str<<std::endl;  
  
  /////////////// Get NP fit functions
  TF1 *fNP = (TF1*)fin_NP->Get((NPCorrFits_str+std::to_string(etabin)).c_str());
  if(!fin_NP || !fNP){
    std::cout<<"error, NP file and/or fits not found."<<std::endl; assert(false);}
  
  /////////////// Create/Get JER fit function(s)
  std::cout<<"opening JER file " << in_JERFile<<std::endl;
  TFile* fin_JER=TFile::Open(in_JERFile.c_str());
  if(!fin_JER){
    std::cout<<"error, JER file not found or not open."<<std::endl; assert(false); }

  std::cout<<"etabin="<<etabin<<std::endl;  
  std::string absetarange_str=absetabins_str[etabin]+" < #||{y} < "+absetabins_str[etabin+1];  
  std::string thyname="h"+std::to_string(xsecorder)+"100"+std::to_string(etabin+1)+"00";
  
  std::string orderstring="LO";
  if(xsecorder==1)orderstring="N"+orderstring;
  else if(xsecorder==2)orderstring=="NN"+orderstring;
  
  //////////// NLO Thy calculation xsecions get from file ////////////////////////////////////////////////////////
  std::cout<<"plotting NLO xsecs!"<<  std::endl;  

  TH1D *theory = (TH1D*)fin_NLO->Get(thyname.c_str());
  theory->Sumw2(true);
  theory->Scale(1.e-03);//picobarns to nanobarns
  if(printBaseDebug)  theory->Print("Base"); 

  theory->SetTitle((orderstring+" Inclusive #sigma_{jet}, "+absetarange_str+";Jet p_{T};"+ddxsec_yax).c_str());
  theory->GetXaxis()->SetNoExponent(true);
  theory->GetXaxis()->SetMoreLogLabels(true);
  theory->SetLineColor(kTeal);
  if(printBaseDebug)theory->Print("base");  

  //errors, PDF systematic
  std::string thyPDFNegerrname="h"+std::to_string(xsecorder)+"100"+std::to_string(etabin+1)+"01";//hOSPxRx2 == -1. * hOSPxRx1 (symm)  
  TH1D *theoryPDFNegerr = (TH1D*)fin_NLO->Get(thyPDFNegerrname.c_str());
  theoryPDFNegerr->Sumw2(true);
  theoryPDFNegerr->Scale(1.e-03);//picobarns to nanobarns
  if(printBaseDebug)  theoryPDFNegerr->Print("Base"); 
  
  std::string thyPDFPoserrname="h"+std::to_string(xsecorder)+"100"+std::to_string(etabin+1)+"02";//hOSPxRx2 == -1. * hOSPxRx1 (symm)  
  TH1D *theoryPDFPoserr = (TH1D*)fin_NLO->Get(thyPDFPoserrname.c_str());
  theoryPDFPoserr->Sumw2(true);
  theoryPDFPoserr->Scale(1.e-03);//picobarns to nanobarns
  if(printBaseDebug)  theoryPDFNegerr->Print("Base"); 
  
  //set theory calculation errors to that of the PDF errs
  for(int i=1; i<=theory->GetNbinsX(); i++) theory->SetBinError(i,theoryPDFPoserr->GetBinContent(i));
  
  std::string thy2PtScaleNegerrname="h"+std::to_string(xsecorder)+"100"+std::to_string(etabin+1)+"06";//2 PT SCALE UNCERTAINTY, ASYMMETRIC
  TH1D *theory2PtScaleNegerr = (TH1D*)fin_NLO->Get(thy2PtScaleNegerrname.c_str());
  theory2PtScaleNegerr->Sumw2(true);
  theory2PtScaleNegerr->Scale(1.e-03);//picobarns to nanobarns
  if(printBaseDebug)  theory2PtScaleNegerr->Print("Base"); 
  
  std::string thy2PtScalePoserrname="h"+std::to_string(xsecorder)+"100"+std::to_string(etabin+1)+"07";//07 IS BROKEN WARNING
  TH1D *theory2PtScalePoserr = (TH1D*)fin_NLO->Get(thy2PtScalePoserrname.c_str());
  theory2PtScalePoserr->Sumw2(true);
  theory2PtScalePoserr->Scale(1.e-03);//picobarns to nanobarns
  if(printBaseDebug)  theory2PtScalePoserr->Print("Base"); 
  
  std::string thy6PtScaleNegerrname="h"+std::to_string(xsecorder)+"100"+std::to_string(etabin+1)+"08";//2 PT SCALE UNCERTAINTY, ASYMMETRIC
  TH1D *theory6PtScaleNegerr = (TH1D*)fin_NLO->Get(thy6PtScaleNegerrname.c_str());
  theory6PtScaleNegerr->Sumw2(true);
  theory6PtScaleNegerr->Scale(1.e-03);//picobarns to nanobarns
  if(printBaseDebug)  theory6PtScaleNegerr->Print("Base"); 
  
  std::string thy6PtScalePoserrname="h"+std::to_string(xsecorder)+"100"+std::to_string(etabin+1)+"09";//07 IS BROKEN WARNING
  TH1D *theory6PtScalePoserr = (TH1D*)fin_NLO->Get(thy6PtScalePoserrname.c_str());
  theory6PtScalePoserr->Sumw2(true);
  theory6PtScalePoserr->Scale(1.e-03);//picobarns to nanobarns
  if(printBaseDebug)  theory6PtScalePoserr->Print("Base"); 
  
  std::cout<<"done making "<<orderstring<<" hist."<<std::endl;
  //assert(false);
  //-----------------------------------
  
  //////////// Apply NPs to theory + Plot////////////////////////////////////////////////////////
  std::cout<<"plotting NP+NLO xsecs!"<<  std::endl;  
  TH1D *theory_NP  = (TH1D*) applyNPtoxsec(theory, fNP);
  theory_NP->SetName(("h"+outputTag+"_wNP").c_str());
  
  theory_NP->SetTitle(("NP+"+orderstring+" Inclusive #sigma_{jet}, "+absetarange_str+";Jet p_{T} [GeV];"+ddxsec_yax).c_str());
  theory_NP->GetXaxis()->SetNoExponent(true);
  theory_NP->GetXaxis()->SetMoreLogLabels(true);
  theory_NP->SetLineColor(kTeal);
  if(printBaseDebug)theory_NP->Print("base");  
  
  std::cout<<"done making NLO+NPhist."<<std::endl;
  //-----------------------------------
  
  TH1D* theory_clone=NULL;
  TH1D* theory_NP_clone=NULL;
  if(etabin==0 || etabin==1){//get rid of last bin
    theory_clone=(TH1D*)theory->TH1::Rebin(n_thybins_incl, 
					   ( ( (std::string) theory->GetName() ) +"_specialrebin").c_str() ,   thyBins_incl );    
    theory_NP_clone=(TH1D*)theory_NP->TH1::Rebin(n_thybins_incl, 
						 ( ( (std::string) theory_NP->GetName() ) +"_specialrebin").c_str() ,     thyBins_incl ); }
  else{
    theory_clone=(TH1D*)theory->Clone(( ( (std::string) theory->GetName() ) +"_0").c_str());
    theory_NP_clone=(TH1D*)theory_NP->Clone(( ( (std::string) theory_NP->GetName() ) +"_0").c_str());}
  
  theory_clone->SetLineColor(kCyan+4);    
  theory_NP_clone->SetLineColor(kCyan+4);    
  
  int n_thybins=theory_clone->GetNbinsX();
  std::cout<<"n_thybins="<<n_thybins<<std::endl;
  std::cout<<"n_thybins+1="<<n_thybins+1<<std::endl;
  std::cout<<std::endl;
  double thybins[n_thybins+1]={0.};//+1 because this array represents the bin edges, not the bins themselves
  for(int i=0; i<(n_thybins+1); i++){
    std::cout<<"i="<<i<<std::endl;
    std::cout<<"i+1="<<i+1<<std::endl;
    if(i==n_thybins)
      thybins[n_thybins]=theory->GetBinLowEdge(i)+theory->GetBinWidth(i);
    else{
      thybins[i]=theory->GetBinLowEdge(i+1);      
    }
    if(printBaseDebug)
      std::cout<<"thybins["<<i<<"]="<<thybins[i]<<std::endl;
  }

  
  TF1 *fJER_default = new TF1("fJER_def",gJER_ynew_str.c_str(),thybins[0], thybins[n_thybins]);  //for comparison only
  fJER_default->SetLineColor(kBlue);
  fJER_default->SetMinimum(0.); //do i need to set min/max really?
  fJER_default->SetMaximum(0.25);
  
  TF1 *fJER = NULL;
  if(useHistSigmaFit)
    fJER=(TF1*)fin_JER->Get("SigmaFit_h");
  else
    fJER=(TF1*)fin_JER->Get("SigmaFit_f");
  
  if(!fJER){
    std::cout<<"error, JER fit not open or not found!."<<std::endl; assert(false);}
  fJER->SetLineColor(kRed);
  fJER->SetMinimum(0.); //do i need to set min/max really?
  fJER->SetMaximum(0.25);
  std::cout<<std::endl;
  
  
  /////////////// plots to to check NPs 
  // 2x2 canv of |y| bins
  TLine* lineatone=new TLine(thybins[0],1.,thybins[n_thybins],1.); lineatone->SetLineStyle(7);
  TCanvas *plot_NP = new TCanvas("plot_NP", ("plot y"+std::to_string(etabin)+" NPs").c_str(),1200,1000);
  //plot_y0_to_y3_NPs->Divide(2,2);
  //-----------------------------------
  plot_NP->cd()->SetLogx(1);
  plot_NP->cd();
  TH1F *hNP = (TH1F*)( (TH1F*)fNP->GetHistogram()
		       )->Clone( ("hNP_y"+std::to_string(etabin)).c_str());
  hNP->SetTitle( ("NP Corr Fit for "+absetarange_str+";Jet p_{T};NP Corr. Factor" ).c_str() );
  hNP->GetXaxis()->SetNoExponent(true);
  hNP->GetXaxis()->SetMoreLogLabels(true);
  hNP->SetLineColor(kRed);
  hNP->DrawClone();  
  lineatone->SetX2( hNP->GetBinLowEdge(hNP->GetNbinsX())+ hNP->GetBinWidth(hNP->GetNbinsX()) );
  lineatone->DrawClone();
  //-----------------------------------  
  
  
  /////////////// plot Core p_T Resolution ////////////////////////////////////////////////////////
  std::cout<<"plotting JER!"<<  std::endl;  
  //TF1 is irritating to set titles for+blah... so use this  
  TH1D *hJER = (TH1D*)( (TH1D*)fJER->GetHistogram() )->Clone(("hJER_"+std::to_string(etabin)).c_str());
  hJER->SetTitle(("Gauss Core JER Fit, "+absetarange_str+"; Jet p_{T} ; #sigma / #mu from Fit").c_str());
  hJER->GetXaxis()->SetNoExponent(true);
  hJER->GetXaxis()->SetMoreLogLabels(true);
  TH1D *hJER_default = (TH1D*)(
			       (TH1D*)fJER_default->GetHistogram()
			       )->Clone("hJER_default");//for comparison only
  hJER_default->SetTitle( (" Gauss Core JER Fit,"+absetarange_str+
			   " ; Jet p_{T} ; #sigma / #mu from Fit").c_str());//sets hist title+xaxis+yaxis title
  hJER_default->GetXaxis()->SetNoExponent(true);
  hJER_default->GetXaxis()->SetMoreLogLabels(true);
  //-----------------------------------  
  TCanvas *plotJER = new TCanvas("plotJER", "plot JER",900,600);
  plotJER->cd()->SetLogx(1);
  plotJER->cd();  
  //hJER_default->DrawClone("HIST E"); basically for a debug
  hJER->DrawClone("HIST E");    
  //-----------------------------------

  //TO DO REPLACE ME W/TEXT SAYING SAME THING
  //TLegend *leg=new TLegend(0.53,0.50,0.85,0.85);
  //leg->AddEntry(hJER_ynew, "PY8 MC CUETP8M1 JER Fits New");
  //leg->AddEntry(hJER_default, "Default JER Fits from APS DNP 10/3/17");
  //leg->AddEntry(fJER,"PYTHIA8 CUETP8M1 @ 5.02 TeV", "");
  //leg->AddEntry(fJER,"ak4PFJets", "");
  //leg->AddEntry(fJER,"0<|y|<2.0", "l"); 
  //leg->SetTextFont(42);
  //leg->SetFillColor(kWhite);
  //leg->Draw();
  //-----------------------------------

  
  
  
  
  
  
  /////////////// draw xsec hists from file
  //2X2 CANV
  TCanvas *plot_NLOxsec = new TCanvas("plot_NLOxsec", "plot NLOxsec",1200,1000);
  //--------------------------------------------
  plot_NLOxsec->cd()->SetLogx(1);
  plot_NLOxsec->cd()->SetLogy(1);
  plot_NLOxsec->cd();
  theory->SetTitle((orderstring+" #sigma_{jet}, "+absetarange_str+";Jet p_{T} [GeV];"+ddxsec_yax).c_str());
  theory->GetXaxis()->SetNoExponent(true);
  theory->GetXaxis()->SetMoreLogLabels(true);
  theory->SetLineColor(kRed);
  theory->DrawClone("HIST E");    
  //--------------------------------------------
  
  //2X2 CANV THEORY NPNLO xsec
  TCanvas *plot_NPNLOxsec = new TCanvas("plot_NPNLOxsec", "plot NPNLOxsec",1200,1000);
  //--------------------------------------------
  plot_NPNLOxsec->cd()->SetLogx(1);
  plot_NPNLOxsec->cd()->SetLogy(1);
  plot_NPNLOxsec->cd();
  theory_NP->DrawClone("HIST E");      
  //--------------------------------------------
  
  TCanvas *plot_xsecPosErrs = new TCanvas("plot_xsecPosErrs", "plot xsec +Errors",1200,1000);
  plot_xsecPosErrs->cd()->SetLogx(1);
  plot_xsecPosErrs->cd()->SetLogy(1);
  plot_xsecPosErrs->cd();
  
  theoryPDFPoserr->SetLineColor(kBlue);       
  theory2PtScalePoserr->SetLineColor(kRed);
  theory6PtScalePoserr->SetLineColor(kMagenta);
  
  theory6PtScalePoserr->SetTitle(("Theory "+orderstring+" #sigma_{Jet} (+) Errors;Jet p_{T} [GeV]; (+) Errors [nb/GeV]").c_str());
  
  theory6PtScalePoserr->SetMinimum(1.e-06);
  theory6PtScalePoserr->SetMaximum(1.e-02);
  theory6PtScalePoserr->GetXaxis()->SetNoExponent(true);
  theory6PtScalePoserr->GetXaxis()->SetMoreLogLabels(true);  

  theory6PtScalePoserr->DrawClone("HIST");
  theory2PtScalePoserr->DrawClone("HIST SAME");  
  theoryPDFPoserr     ->DrawClone("HIST SAME");
  
  TLegend* poserrsleg=new TLegend(0.1,0.7,0.4,0.9);
  poserrsleg->SetBorderSize(0.);
  poserrsleg->SetFillStyle(0);  
  poserrsleg->AddEntry( theory6PtScalePoserr, "6 Pt. Scale Error(Asymm.)" , "l");
  poserrsleg->AddEntry( theory2PtScalePoserr, "2 Pt. Scale Error(Asymm.)" , "l");
  poserrsleg->AddEntry( theoryPDFPoserr     , "PDF Error (Symm.)" , "l");
  
  poserrsleg->Draw();
  
  //--------------------------------------------
  TCanvas *plot_xsecNegErrs = new TCanvas("plot_xsecNegErrs", "plot xsec +Errors",1200,1000);
  plot_xsecNegErrs->cd()->SetLogx(1);
  plot_xsecNegErrs->cd()->SetLogy(1);
  plot_xsecNegErrs->cd();
  
  theoryPDFNegerr->SetLineColor(kBlue);       
  theory2PtScaleNegerr->SetLineColor(kRed);
  theory6PtScaleNegerr->SetLineColor(kMagenta);
  
  theory6PtScaleNegerr->SetTitle(("Theory "+orderstring+" #sigma_{Jet} (-) Errors;Jet p_{T}[GeV]; #||{(-) Error} [nb/GeV]").c_str());
  
  theory6PtScaleNegerr->Scale(-1.);
  theory2PtScaleNegerr->Scale(-1.);
  theoryPDFNegerr     ->Scale(-1.);
  
  theory6PtScaleNegerr->SetMinimum(1.e-06);
  theory6PtScaleNegerr->SetMaximum(1.e-02);
  theory6PtScaleNegerr->GetXaxis()->SetNoExponent(true);
  theory6PtScaleNegerr->GetXaxis()->SetMoreLogLabels(true);
  
  
  theory6PtScaleNegerr->DrawClone("HIST");
  theory2PtScaleNegerr->DrawClone("HIST SAME");
  theoryPDFNegerr     ->DrawClone("HIST SAME");
  
  TLegend* negerrsleg=new TLegend(0.1,0.7,0.4,0.9);
  negerrsleg->SetBorderSize(0.);
  negerrsleg->SetFillStyle(0);  
  negerrsleg->AddEntry( theory6PtScaleNegerr, "6 Pt. Scale Error (Asymm.)" , "l");
  negerrsleg->AddEntry( theory2PtScaleNegerr, "2 Pt. Scale Error (Asymm.)" , "l");
  negerrsleg->AddEntry( theoryPDFNegerr     , "PDF Error (Symm.)" , "l");
  
  negerrsleg->Draw();
  
  
  

  
  

  


            

  
  ///Cubic Spline Interpolations using Cross sections; has issues for p_T beyond first/last bin centers
  TSpline3 *spline3=NULL, *spline3_NP=NULL;
  TF1      *spline3_ext=NULL, *spline3_NP_ext=NULL;
  TCanvas *plot_splines=NULL;
  TLegend *leg_spline=NULL;
  if(useSplineWeights){
    std::cout<<"using spline weights!!!"<<std::endl;
    if(useSplineExt)std::cout<<"also adding spline weight extension!!!"<<std::endl;


    /// check spectra with splines
    std::cout<<"drawing Tspline and cross section on canvas"<<std::endl;
    plot_splines = new TCanvas("plot_splines", "plot splines",1200,800);
    plot_splines->cd()->SetLogx(1);
    plot_splines->cd()->SetLogy(1);
    plot_splines->cd();
    
    // NLO //
    std::cout<<"creating TSpline of cross section hist"<<std::endl;  
    spline3 = new TSpline3( theory_clone);
    spline3->SetName( ( (std::string)theory_clone->GetName() + "_spline3").c_str() );
    spline3->SetLineColor(kAzure);  

    double y2=calc_spline3ext_y2(theory_clone);//this is also used to help set range for the plot ; so it's outside of the 'usesplineext' flag
    double x1=-1., x2=-1.;
    //double y2=theory_clone->GetBinContent(theory_clone->GetNbinsX());
    if(useSplineExt){
      //figure out starting values for parabolic extension of spline
      //x1=theory_clone->GetBinCenter(theory_clone->GetNbinsX());// bad idea, spline3's val at x1 is approx the val of the last bin's height. 
      if(fitType.find("v1")!=std::string::npos)x1=theory_clone->GetBinLowEdge(theory_clone->GetNbinsX()); //v1
      else if(fitType.find("v2")!=std::string::npos)x1=theory_clone->GetBinCenter(theory_clone->GetNbinsX()-1);//v2
      else if(fitType.find("v3")!=std::string::npos)x1=theory_clone->GetBinLowEdge(theory_clone->GetNbinsX()-1); //v3
      else assert(false);
      x2=theory_clone->GetBinLowEdge(theory_clone->GetNbinsX()) + theory_clone->GetBinWidth(theory_clone->GetNbinsX());                    
      
      //double xprime=theory_clone->GetBinCenter(theory_clone->GetNbinsX());
      //double xprime=theory_clone->GetBinCenter(theory_clone->GetNbinsX()) - 0.25*theory_clone->GetBinWidth(theory_clone->GetNbinsX());
      //double yprime=theory_clone->GetBinContent(theory_clone->GetNbinsX());
      
      spline3_ext= new TF1( ((std::string)spline3->GetName()+"_tf1ext").c_str(), "[0]*pow(x,2)+[1]*x+[2]", x1, x2);    
      make_spline3_ext(  (TSpline3*)spline3, (TF1*)spline3_ext , (double)x1, (double)x2, (double)y2);
      //make_spline3_extv2(  (TSpline3*)spline3, (TF1*)spline3_ext , (double)x1, (double)x2, (double)y2, (double)xprime, (double)yprime);    
      spline3_ext->SetLineColor(kAzure-1);
    }
    
    // NLO+NP //      
    spline3_NP = new TSpline3( theory_NP_clone);
    spline3_NP->SetName( ( (std::string)theory_NP_clone->GetName() + "_spline3").c_str() );
    spline3_NP->SetLineColor(kAzure-8);  
    
    
    if(useSplineExt){
      //double NP_y2=theory_NP_clone->GetBinContent(theory_NP_clone->GetNbinsX());
      double NP_y2=calc_spline3ext_y2(theory_NP_clone);    
      
      spline3_NP_ext= new TF1( ((std::string)spline3_NP->GetName()+"_tf1ext").c_str(), "[0]*pow(x,2)+[1]*x+[2]", x1, x2);    
      make_spline3_ext(  (TSpline3*)spline3_NP, (TF1*)spline3_NP_ext , (double)x1, (double)x2, (double)NP_y2);    
      spline3_NP_ext->SetLineColor(kAzure-9);
    }
      
    //draw
    theory_clone->SetTitle(("Cubic Splines for "+absetarange_str+";Jet p_{T};Smear Weight").c_str());    
    theory_clone->SetAxisRange( y2/2., theory_clone->GetMaximum()*2., "Y");
    
    theory_clone->DrawClone("HIST E");
    spline3->Draw("SAME");    
    if(useSplineExt)spline3_ext->Draw("SAME");
    
    theory_NP_clone->DrawClone("HIST E SAME");  
    spline3_NP->Draw("SAME");
    if(useSplineExt)spline3_NP_ext->Draw("SAME");
    
    leg_spline=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
    leg_spline->AddEntry(theory_clone ,    ("Weighted "+orderstring+" Jet Counts,"+absetarange_str           ).c_str() , "lp");
    leg_spline->AddEntry(spline3   ,       ("Cubic Spline "+orderstring+","+absetarange_str                  ).c_str() , "l");
    //
    if(useSplineExt)leg_spline->AddEntry(spline3_ext  ,    ("Parabolic Spline Extension "+orderstring+","+absetarange_str    ).c_str() , "l");
    //
    leg_spline->AddEntry(theory_NP_clone , ("Weighted NP+"+orderstring+" Jet Counts,"+absetarange_str        ).c_str() , "lp");
    leg_spline->AddEntry(spline3_NP   ,    ("Cubic Spline NP+"+orderstring+","+absetarange_str               ).c_str() , "l");
    if(useSplineExt)leg_spline->AddEntry(spline3_NP_ext  , ("Parabolic Spline Extension NP+"+orderstring+","+absetarange_str ).c_str() , "l");
    leg_spline->Draw();    
  }

  // for fits
  const int startbin=1, endbin=theory_clone->GetNbinsX();//NP and nonNP hist should always have the same binning
  const float xlo=theory_clone->GetBinLowEdge(startbin);
  const float xhi=theory_clone->GetBinLowEdge(endbin) + theory_clone->GetBinWidth(endbin);  
  
  
  ///Log fit: A(B/pt)^[C+ D Log(pt/B) + E(Log(pt/B))^2] using Cross sections; currently no known issues
  TF1 *logFit  =NULL, *modLogFit  =NULL, *modLog2Fit  =NULL, *modLog3Fit  =NULL, *modLog4Fit  =NULL,*modLog5Fit  =NULL;
  TF1* logFit_4Toy=NULL;

  TF1 *logFit_NP=NULL, *modLogFit_NP=NULL, *modLog2Fit_NP=NULL, *modLog3Fit_NP=NULL, *modLog4Fit_NP=NULL,*modLog5Fit_NP=NULL;
  TF1* logFit_NP_4Toy=NULL;
  
  TCanvas *plot_logfits=NULL,*plot_logfits_NP=NULL;
  TLegend *leg_logfits=NULL,*leg_logfits_NP=NULL;  
  //  if(useFitWeights){
  if(useFitWeights && useModLogFit){
    std::cout<<"using modlogfit weights!!!"<<std::endl;

    // --- NLO --- //
    /// check spectra with fits
    plot_logfits = new TCanvas("plot_logfits", "plot logfits",1200,800);
    plot_logfits->cd()->SetLogx(1);
    plot_logfits->cd()->SetLogy(1);
    plot_logfits->cd();
    
    std::cout<<"fitting hist theory_clone"<<std::endl;
    logFit    =new TF1("logfit    ","[0] * ( TMath::Power( ( [1]/x ) ,   [2] ) )",xlo,xhi);
    modLogFit =new TF1("modlogfit ","[0] * ( TMath::Power( ( [1]/x ) ,   [2] + [3]*TMath::Log10( x/[1] ) ) )",xlo,xhi);
    modLog2Fit=new TF1("modlog2fit","[0] * ( TMath::Power( ( [1]/x ) , ( [2] + [3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )  ) ) )",xlo,xhi);
    modLog3Fit=new TF1("modlog3fit","[0] * ( TMath::Power( ( [1]/x ) , ( [2] + [3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )  + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 )  ) ) )",xlo,xhi);
    modLog4Fit=new TF1("modlog4fit","[0] * ( TMath::Power( ( [1]/x ) , ( [2] + [3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )  + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 ) + [6]*TMath::Power( TMath::Log10( x/[1] ), 4 ) )               ) )",xlo,xhi);
    modLog5Fit=new TF1("modlog5fit","[0] * ( TMath::Power( ( [1]/x ) , ( [2] + [3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )  + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 ) + [6]*TMath::Power( TMath::Log10( x/[1] ), 4 ) + [7]*TMath::Power( TMath::Log10( x/[1] ), 5 ) )               ) )",xlo,xhi);

    
    makeModLogFits(theory_clone, logFit, modLogFit, modLog2Fit, modLog3Fit, modLog4Fit, modLog5Fit);
    
    if(etabin==0)logFit_4Toy=modLog5Fit;
    else if(etabin==1)logFit_4Toy=modLog2Fit;
    else if(etabin==2)logFit_4Toy=modLog4Fit;
    else if(etabin==3)logFit_4Toy=modLog2Fit;
    else if(etabin==4)logFit_4Toy=modLog4Fit;
    else if(etabin==5)logFit_4Toy=modLog4Fit;
    


    //draw
    std::cout<<"drawing fits and cross section on canvas"<<std::endl;
    theory_clone->SetTitle(("Modified Logarithm Fits for "+orderstring+", "+absetarange_str+";Jet p_{T};Smear Weight").c_str());    
    theory_clone->DrawClone("HIST E");    
    
    logFit->Draw("SAME");
    modLogFit->Draw("SAME");
    modLog2Fit->Draw("SAME");
    modLog3Fit->Draw("SAME");    
    modLog4Fit->Draw("SAME");    
    modLog5Fit->Draw("SAME");    
    
    leg_logfits=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
    
    leg_logfits->AddEntry(theory_clone , ("Weighted "+orderstring+" Jet Counts for "+absetarange_str).c_str() , "lp");//change desc TO DO
    leg_logfits->AddEntry(logFit      , "Log Fit" , "l");//change desc TO DO
    leg_logfits->AddEntry(modLogFit   , "Mod. Log Fit" , "l");
    leg_logfits->AddEntry(modLog2Fit   , "2nd Order Mod. Log Fit" , "l");    
    leg_logfits->AddEntry(modLog3Fit   , "3rd Order Mod. Log Fit" , "l");    
    leg_logfits->AddEntry(modLog4Fit   , "4th Order Mod. Log Fit" , "l");    
    leg_logfits->AddEntry(modLog5Fit   , "5th Order Mod. Log Fit" , "l");    
    leg_logfits->Draw();    
    
    
    // --- NLO+NP --- //
    /// check spectra with fits    
    plot_logfits_NP = new TCanvas("plot_logfits_NP", "plot logfits NP",1200,800);
    plot_logfits_NP->cd()->SetLogx(1);
    plot_logfits_NP->cd()->SetLogy(1);
    plot_logfits_NP->cd();
    
    std::cout<<"fitting hist theory_NP_clone"<<std::endl;
    logFit_NP    =new TF1("logfit_NP    ","[0] * ( TMath::Power( ( [1]/x ) , ( [2]                                                                                                                         ) ) )",xlo,xhi);
    modLogFit_NP =new TF1("modlogfit_NP ","[0] * ( TMath::Power( ( [1]/x ) , ( [2]+[3]*TMath::Log10( x/[1] )                                                                                               ) ) )",xlo,xhi);
    modLog2Fit_NP=new TF1("modlog2fit_NP","[0] * ( TMath::Power( ( [1]/x ) , ( [2]+[3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )                                                ) ) )",xlo,xhi);
    modLog3Fit_NP=new TF1("modlog3fit_NP","[0] * ( TMath::Power( ( [1]/x ) , ( [2]+[3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 ) + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 ) ) ) )",xlo,xhi);
    modLog4Fit_NP=new TF1("modlog4fit_NP","[0] * ( TMath::Power( ( [1]/x ) , ( [2]+[3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 ) + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 ) + [6]*TMath::Power( TMath::Log10( x/[1] ), 4 ) ) ) )",xlo,xhi);
    modLog5Fit_NP=new TF1("modlog5fit_NP","[0] * ( TMath::Power( ( [1]/x ) , ( [2]+[3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 ) + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 ) + [6]*TMath::Power( TMath::Log10( x/[1] ), 4 ) + [7]*TMath::Power( TMath::Log10( x/[1] ), 5 ) ) ) )",xlo,xhi); 
    
    
    makeModLogFits(theory_NP_clone, logFit_NP, modLogFit_NP, modLog2Fit_NP, modLog3Fit_NP, modLog4Fit_NP, modLog5Fit_NP);

    if(etabin==0)logFit_NP_4Toy=modLog5Fit_NP;
    if(etabin==1)logFit_NP_4Toy=modLog5Fit_NP;
    if(etabin==2)logFit_NP_4Toy=modLog4Fit_NP;
    if(etabin==3)logFit_NP_4Toy=modLog2Fit_NP;
    if(etabin==4)logFit_NP_4Toy=modLog4Fit_NP;
    if(etabin==5)logFit_NP_4Toy=modLog4Fit_NP;

    
    
    //draw
    std::cout<<"drawing NP fits and cross section on canvas"<<std::endl;
    theory_NP_clone->SetTitle( ("Modified Logarithm Fits for "+orderstring+"+NP,"+absetarange_str+" ;Jet p_{T};Smear Weight").c_str());        
    theory_NP_clone->DrawClone("HIST E");        
    
    logFit_NP->Draw("SAME");
    modLogFit_NP->Draw("SAME");
    modLog2Fit_NP->Draw("SAME");
    modLog3Fit_NP->Draw("SAME");    
    modLog4Fit_NP->Draw("SAME");    
    modLog5Fit_NP->Draw("SAME");    
    
    leg_logfits_NP=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
    leg_logfits_NP->AddEntry(theory_NP_clone , ("Weighted "+orderstring+"+NP Jet Counts for "+absetarange_str).c_str() , "lp");//change desc TO DO
    //leg_logfits_NP->AddEntry(theory_NP , "Weighted NLO+NP Jet Counts for"+absetarange_str , "lp");//change desc TO DO
    leg_logfits_NP->AddEntry(logFit_NP      , "Log Fit" , "l");//change desc TO DO
    leg_logfits_NP->AddEntry(modLogFit_NP   , "Mod. Log Fit" , "l");
    leg_logfits_NP->AddEntry(modLog2Fit_NP   , "2nd Order Mod. Log Fit" , "l");    
    leg_logfits_NP->AddEntry(modLog3Fit_NP   , "3rd Order Mod. Log Fit" , "l");    
    leg_logfits_NP->AddEntry(modLog4Fit_NP   , "4th Order Mod. Log Fit" , "l");    
    leg_logfits_NP->AddEntry(modLog5Fit_NP   , "5th Order Mod. Log Fit" , "l");    
    leg_logfits_NP->Draw();    

  }


  ///7 TeV fit: two variations using Cross sections, each w/ 3 free parameters min. note; sqrt(s) should be in GeV
  ///cite: CERN-PH-EP/2011-053, 2018/10/22, CMS-QCD-10-011, 
  ///    : "Measurement of the Inclusive Jet Cross Section in pp Collisions at sqrt(s)=7 TeV"
  ///    : arixiv:1106.0208v1 [hep-ex] 1 Jun 2011
  ///    : Submitted to Physical Review Letters
  /// A) ddxsec[p_T] = D * [(p_T/GeV)^A] * [(1.-(2. *p_T/sqrt(s)]*cosh(ymin))^B] * [exp(-1*C/(p_T/GeV))]
  /// B) ddxsec[p_T] = [integral(ddxsec[p_T])] * [(p_T/GeV)^A] * [(1 - (2 * p_T / sqrt(s)]*cosh(ymin))^B] * [exp(-1*C/(p_T/GeV))]
  TF1 *_7TeVFitA  =NULL, *_7TeVFitB  =NULL;
  TF1 *_7TeVFitA_NP=NULL, *_7TeVFitB_NP=NULL;
  
  TCanvas *plot_7tevfits=NULL,*plot_7tevfits_NP=NULL;
  TLegend  *leg_7tevfits=NULL, *leg_7tevfits_NP=NULL;  
  
  //const float xlo_7TeV=theory_clone->GetBinLowEdge(startbin);//trying out for the 7TeV style fit
  //const float xhi_7TeV=theory_clone->GetBinCenter(endbin);// + theory_clone->GetBinWidth(endbin);  
  
  

  if(useFitWeights && use7TeVFit){
    std::cout<<"using 7TeVfit weights!!!"<<std::endl;
    
    std::cout<<"making 7tev fit strings"<<std::endl;
    float _coshymin=2.*TMath::CosH(absetabins[etabin]);//putting the 2. in here for simplicity
    std::string _coshymin_str=std::to_string(_coshymin);
    std::cout<<"_coshymin_str="<<_coshymin_str<<std::endl;
    
    float norm=theory_clone->Integral("width");
    std::cout<<"norm="<<norm<<std::endl;

    std::string norm_str=std::to_string(norm);
    std::cout<<"norm_str="<<norm_str<<std::endl;

    float norm_NP=theory_NP_clone->Integral("width");
    std::cout<<"norm_NP="<<norm_NP<<std::endl;
    
    std::string norm_NP_str=std::to_string(norm_NP);
    std::cout<<"norm_NP_str="<<norm_NP_str<<std::endl;    
    
    /// A) ddxsec[p_T] = D * [(p_T/GeV)^A] * [(1.-(2.*p_T/sqrt(s)]*cosh(ymin))^B] * [exp(-1*C/(p_T/GeV))]
    //    std::string _7tevfit_str="(TMath::Power(x,-1.*[0]))*(TMath::Power((1.-(2.*x*"+_coshymin_str+"*/(5.02*1000.) )),[1]))* (TMath::Exp(-1.*([2]/x)))";
    std::string _7tevfit_str="(TMath::Power(x,-1.*[0]))*(TMath::Power((1.-("+_coshymin_str+"*x/(5.02*1000.) )),[1]))* (TMath::Exp(-1.*([2]/x)))";
    std::cout<<"_7tevfit_str="<<_7tevfit_str<<std::endl;
    
    std::string _7tevfitA_str="([3])* "+_7tevfit_str;    //this string is same for nonNP and NP
    std::cout<<"_7tevfitA_str="<<_7tevfitA_str<<std::endl;
    
    // --- NLO --- //    
    /// check spectra with fits
    plot_7tevfits = new TCanvas("plot_7tevfits", "plot 7tevfits",1200,800);
    plot_7tevfits->cd()->SetLogx(1);
    plot_7tevfits->cd()->SetLogy(1);
    plot_7tevfits->cd();
    


  /// B) ddxsec[p_T] = [integral(ddxsec[p_T])] * [(p_T/GeV)^A] * [(1 - (2 * p_T / sqrt(s)]*cosh(ymin))^B] * [exp(-1*C/(p_T/GeV))]
    std::string _7tevfitB_str="("+norm_str+")* "+_7tevfit_str;
    std::cout<<"_7tevfitB_str="<<_7tevfitB_str<<std::endl;    
    
    std::cout<<std::endl<<"making _7TeVFitA TF1"<<std::endl;
    _7TeVFitA =new TF1("_7tevfitA", (_7tevfitA_str     ).c_str() ,xlo,xhi);
    
    std::cout<<std::endl<<"making _7TeVFitB TF1"<<std::endl;
    _7TeVFitB =new TF1("_7tevfitB" , (_7tevfitB_str).c_str(),xlo,xhi);
    
    make7TeVFits(theory_clone,  _7TeVFitA, _7TeVFitB, norm);

    //draw
    std::cout<<"drawing fits and cross section on canvas"<<std::endl;
    theory_clone->SetTitle(("7 TeV Fits for "+orderstring+", "+absetarange_str+";Jet p_{T};Smear Weight").c_str());    
    theory_clone->DrawClone("HIST E");        
    _7TeVFitA->Draw("SAME");
    _7TeVFitB->Draw("SAME");
    
    leg_7tevfits=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");    
    leg_7tevfits->AddEntry(theory_clone , ("Weighted "+orderstring+" Jet Counts for"+absetarange_str).c_str() , "lp");//change desc TO DO
    leg_7tevfits->AddEntry(_7TeVFitB   , ("4 Param. 7 TeV Fit for "+orderstring+"").c_str() , "l");
    leg_7tevfits->AddEntry(_7TeVFitA   , ("3 Param. 7 TeV Fit for "+orderstring+"").c_str() , "l");    
    leg_7tevfits->Draw();    
    
    
    // --- NP+NLO --- //
    /// check spectra with fits
    plot_7tevfits_NP = new TCanvas("plot_7tevfits_NP", "plot 7tevfits NP",1200,800);
    plot_7tevfits_NP->cd()->SetLogx(1);
    plot_7tevfits_NP->cd()->SetLogy(1);
    plot_7tevfits_NP->cd();
    

    std::string _7tevfitB_NP_str="("+norm_NP_str+")* "+_7tevfit_str;
    std::cout<<"_7tevfitB_NP_str="<<_7tevfitB_NP_str<<std::endl;        
    
    std::cout<<"fitting hist theory_NP_clone"<<std::endl;
    _7TeVFitA_NP =new TF1("_7tevfitA_NP",(_7tevfitA_str       ).c_str(),xlo,xhi);
    _7TeVFitB_NP=new TF1( "_7tevfitB_NP" ,(_7tevfitB_NP_str).c_str(),xlo,xhi);
    make7TeVFits(theory_NP_clone, _7TeVFitA_NP, _7TeVFitB_NP, norm_NP);
    
    //draw    
    std::cout<<"drawing NP fits and cross section on canvas"<<std::endl;
    theory_NP_clone->SetTitle( ("7 TeV Fits for NP+"+orderstring+","+absetarange_str+";Jet p_{T};Smear Weight").c_str());    
    theory_NP_clone->DrawClone("HIST E");    
    
    _7TeVFitA_NP->Draw("SAME");
    _7TeVFitB_NP->Draw("SAME");
    
    leg_7tevfits_NP=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
    
    leg_7tevfits_NP->AddEntry(theory_NP_clone , ("Weighted NP+"+orderstring+" Jet Counts for"+absetarange_str).c_str() , "lp");//change desc TO DO
    leg_7tevfits_NP->AddEntry(_7TeVFitA_NP   ,  ("3 Param. 7 TeV Fit for NP+"+orderstring).c_str() , "l");
    leg_7tevfits_NP->AddEntry(_7TeVFitB_NP   ,  ("4 Param. 7 TeV Fit for NP+"+orderstring).c_str() , "l");    
    leg_7tevfits_NP->Draw();    
  }
  
  
  //-----------------------------------------------------------------------------------------//
  //////////////////////  START production of Smeared NLO spectra  ////////////////////////////
  //-----------------------------------------------------------------------------------------//
  std::cout<<"creating TH1 for toy NLO spectra generation, RooUnfoldResponse class, etc."<<std::endl<<std::endl;
  TH1D *theory_rnd       = new TH1D("theory_rnd","theory_rnd",             n_thybins    , thybins);   
  TH1D *smeared_rnd      = new TH1D("smeared_rnd","smeared_rnd",           n_thybins, thybins);   
  TH1D *smeared_rnd_test = new TH1D("smeared_rnd_test","smeared_rnd_test", n_thybins, thybins);   
  
  TH2D* response_th2=new TH2D("response_th2","response_th2",
			      (Int_t)n_thybins, (Double_t*)thybins,
			      (Int_t)n_thybins, (Double_t*)thybins);  

  //--------------------------------------------------------------------------------------------------------//
  if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
    makeToySpectra( (TH1D*)theory_clone, 
		    (TSpline3*)spline3, 
		    (TF1*)fJER, 
		    nEvents, 
		    theory_rnd, 
		    smeared_rnd, 
		    smeared_rnd_test, 
		    (TH2D*)response_th2,
		    (TF1*) spline3_ext); }//		    NULL); }
  else if(useFitWeights&&useModLogFit){//looks like modlog4Fit wins!
    makeToySpectra( (TH1D*)theory_clone, 
		    (TF1*)logFit_4Toy, //		    (TF1*)modLog4Fit, // (TF1*)modLog3Fit, // 		    
		    (TF1*)fJER, 
		    nEvents, 
		    theory_rnd, 
		    smeared_rnd, 
		    smeared_rnd_test, 
		    (TH2D*)response_th2); }  
  else if(useFitWeights&&use7TeVFit){//must decide between A and B TO DO
    makeToySpectra( (TH1D*)theory_clone, 
		    (TF1*)_7TeVFitA,// (TF1*)_7TeVFitB, //
		    (TF1*)fJER, 
		    nEvents, 
		    theory_rnd, 
		    smeared_rnd, 
		    smeared_rnd_test, 
		    (TH2D*)response_th2); }  

  //--------------------------------------------------------------------------------------------------------//

  divideBinWidth(theory_rnd);  //for normalization only.
  divideBinWidth(smeared_rnd); 
  divideBinWidth(smeared_rnd_test); 
  
  double norm_rnd=theory_clone->Integral()/theory_rnd->Integral();
  response_th2->Scale(   norm_rnd);  
  
  theory_rnd->Scale( norm_rnd);
  smeared_rnd->Scale(norm_rnd);  
  smeared_rnd_test->Scale(norm_rnd);  
  
  //plotting smeared+theory toy mc xsection
  theory_rnd->SetLineColor(kRed); 
  smeared_rnd->SetLineColor(kBlue);    
  smeared_rnd_test->SetLineColor(kBlue-5);    
  
  Float_t maxy=theory_rnd->GetMaximum();    std::cout<<"maxy="<<maxy<<std::endl;
  maxy*=2.;  std::cout<<"maxy="<<maxy<<std::endl;
  Float_t miny=smeared_rnd->GetMinimum();   std::cout<<"miny="<<miny<<std::endl;
  miny/=2.; std::cout<<"miny="<<miny<<std::endl;
  theory_rnd->SetAxisRange(miny,maxy,"Y");  
  theory_rnd->SetTitle( ((std::string)("Toy MC Incl. "+orderstring+" #sigma;;"+ddxsec_yax)).c_str() );
  
  
  ///// check plots True/Smeared 
  TCanvas *plot_true_smeared_rat = new TCanvas("plot_true_smeared_rat", "plot true smeared rat",1400,1200);  
  plot_true_smeared_rat->cd();

  
  TPad* overlay_pad=new TPad("overlay_pad","pad for spectra",0.0,0.3,1.0,1.0);
  overlay_pad->SetLogx(1);
  overlay_pad->SetLogy(1);
  overlay_pad->SetBottomMargin(0);
  overlay_pad->Draw();//tells the canvas theres a pad at the coordinates specified in the new constructor
  overlay_pad->cd();
  
  theory_rnd->DrawClone("HIST E");  
  smeared_rnd->DrawClone("HIST E SAME");
  smeared_rnd_test->DrawClone("HIST E SAME");
  
  TLegend* leg_smear=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
  leg_smear->AddEntry(theory_rnd  ,      ("Truth Incl. "+orderstring+" #sigma"            ).c_str(), "lp" );
  leg_smear->AddEntry(smeared_rnd ,      ("Smeared Incl. "+orderstring+" #sigma"          ).c_str(), "lp" );
  leg_smear->AddEntry(smeared_rnd_test , ("Opp Side Smeared Incl. "+orderstring+" #sigma" ).c_str(), "lp" );
  leg_smear->Draw();
  
  plot_true_smeared_rat->cd();//go back to main canvas before doing new pad
  
  TPad* ratio_pad=new TPad("ratio_pad","pad for ratio",0.0,0.05,1.0,0.30);
  ratio_pad->SetLogx(1);
  ratio_pad->SetLogy(0);
  ratio_pad->SetTopMargin(0);
  //ratio_pad->SetBottomMargin(0.3);
  ratio_pad->Draw();//tell canvas there's a pad at the specified coordinates in the new constructor
  ratio_pad->cd();
  
  std::cout<<"creating ratio hist of true/smeared"<<std::endl;
  TH1D *true_smeared_rat=(TH1D*)theory_rnd->Clone("true_smeared_rat");   
  if(printBaseDebug)true_smeared_rat->Print("base");
  true_smeared_rat->SetTitle("; Jet p_T (GeV); True/Smeared");
  true_smeared_rat->Divide(true_smeared_rat,smeared_rnd,1.,1.,"B");   
  true_smeared_rat->GetXaxis()->SetMoreLogLabels(true);
  true_smeared_rat->GetXaxis()->SetNoExponent(true);
  true_smeared_rat->SetAxisRange(0.6,1.4,"Y");
  std::cout<<"creating ratio hist of true/opp side smeared"<<std::endl;

  TH1D *true_smeared_test_rat=(TH1D*)theory_rnd->Clone("true_smeared_test_rat");   
  if(printBaseDebug)true_smeared_test_rat->Print("base");  
  true_smeared_test_rat->SetTitle("; Jet p_T (GeV); True/Smeared");
  true_smeared_test_rat->Divide(true_smeared_test_rat,smeared_rnd_test,1.,1.,"B");   
  true_smeared_test_rat->GetXaxis()->SetMoreLogLabels(true);
  true_smeared_test_rat->GetXaxis()->SetNoExponent(true);
  true_smeared_test_rat->SetAxisRange(0.6,1.4,"Y");
  true_smeared_test_rat->SetLineColor(kRed-2);

  true_smeared_rat->DrawClone("HIST E");        
  true_smeared_test_rat->DrawClone("HIST E SAME");        
  
  lineatone->SetX2(thybins[n_thybins]);
  lineatone->DrawClone();

  TCanvas* plot_response_th2=new TCanvas("plot_response_th2","plot response th2",1200, 1000);
  plot_response_th2->cd()->SetLogx(1);
  plot_response_th2->cd()->SetLogy(1);
  plot_response_th2->cd()->SetLogz(1);
  plot_response_th2->cd();
  
  response_th2->SetTitle(("Toy MC Response Matrix from JER Smeared "+orderstring+" #sigma_{jet}").c_str());

  response_th2->GetXaxis()->SetTitle("RECO Jet p_{T} [GeV]");
  response_th2->GetXaxis()->SetNoExponent(true);
  response_th2->GetXaxis()->SetMoreLogLabels(true);
  response_th2->GetYaxis()->SetTitle("GEN Jet p_{T} [GeV]");
  response_th2->GetYaxis()->SetNoExponent(true);
  response_th2->GetYaxis()->SetMoreLogLabels(true);
  
  setTH2_ZAxisRange(response_th2);  
  response_th2->Draw("COLZ");
  //-----------------------------------------------------------------------------------------//
  //////////////////////  END production of Smeared NLO spectra  ////////////////////////////
  //-----------------------------------------------------------------------------------------//


  




  ///////////////////////////////////////////////////////////////////////////////////////////// 
  ///////////////////////////////////////////////////////////////////////////////////////////// 
  ///////////////////////////////////////////////////////////////////////////////////////////// 
  ///////////////////////////////////////////////////////////////////////////////////////////// 









  //-----------------------------------------------------------------------------------------//
  //////////////////////  START production of Smeared NLO+NP spectra  /////////////////////////
  //-----------------------------------------------------------------------------------------//
  std::cout<<"creating TH1 for toy NP+"+orderstring+" spectra generation"<<std::endl<<std::endl;
  TH1D *theory_rnd_NP    = new TH1D("theory_rnd_NP","theory_rnd_NP", n_thybins, thybins);   
  TH1D *smeared_rnd_NP = new TH1D("smeared_rnd_NP","smeared_rnd_NP", n_thybins, thybins);     
  TH1D *smeared_rnd_NP_test = new TH1D("smeared_rnd_NP_test","smeared_rnd_NP_test", n_thybins, thybins);     

  TH2D* response_NP_th2=new TH2D("response_NP_th2","response_NP_th2",
				     (Int_t)n_thybins, (Double_t*)thybins,
				     (Int_t)n_thybins, (Double_t*)thybins);
  
  //----------------------------------------------------------------------------------------------------------------//
  if(useSplineWeights){
    makeToySpectra(theory_NP_clone, 
		   (TSpline3*)spline3_NP, 
		   (TF1*)fJER, 
		   nEvents, 
		   theory_rnd_NP, 
		   smeared_rnd_NP, 
		   smeared_rnd_NP_test, 
		   (TH2D*) response_NP_th2, 	
		   (TF1*) spline3_NP_ext); }		    //
  else if(useFitWeights&&useModLogFit){
    makeToySpectra(theory_NP_clone, 
		   (TF1*) logFit_NP_4Toy,//(TF1*)modLog4Fit_NP,
		   (TF1*)fJER, 
		   nEvents, 
		   theory_rnd_NP, 
		   smeared_rnd_NP, 
		   smeared_rnd_NP_test, 
		   (TH2D*) response_NP_th2);}
  else if(useFitWeights&&use7TeVFit){//must decide between A and B TO DO
    makeToySpectra( (TH1D*)theory_NP_clone, 
		    (TF1*)_7TeVFitA_NP,// (TF1*)_7TeVFitB_NP, //
		    (TF1*)fJER, 
		    nEvents, 
		    theory_rnd_NP, 
		    smeared_rnd_NP, 
		    smeared_rnd_NP_test, 
		    (TH2D*)response_NP_th2); }  

  //----------------------------------------------------------------------------------------------------------------//
  
  divideBinWidth(theory_rnd_NP);  //for normalization only.
  divideBinWidth(smeared_rnd_NP); 
  divideBinWidth(smeared_rnd_NP_test); 
  
  double norm_rnd_NP=theory_NP_clone->Integral()/theory_rnd_NP->Integral();
  response_NP_th2->Scale(   norm_rnd_NP); 
  
  theory_rnd_NP->Scale( norm_rnd_NP);
  smeared_rnd_NP->Scale(norm_rnd_NP);
  smeared_rnd_NP_test->Scale(norm_rnd_NP);
  
  theory_rnd_NP->SetLineColor(kRed);
  smeared_rnd_NP->SetLineColor(kBlue);    
  smeared_rnd_NP_test->SetLineColor(kBlue-5);    
  
  if(theory_rnd_NP->GetBinContent(theory_rnd_NP->GetNbinsX() )<0.)
    theory_rnd_NP->SetBinContent(theory_rnd_NP->GetNbinsX(),1.e-10);
  if(smeared_rnd_NP->GetBinContent(smeared_rnd_NP->GetNbinsX() )<0.)
    smeared_rnd_NP->SetBinContent(smeared_rnd_NP->GetNbinsX(),1.e-08);
  
  
  Float_t NPmaxy=theory_rnd_NP->GetMaximum();    std::cout<<"NPmaxy="<<NPmaxy<<std::endl;
  NPmaxy*=2.;  std::cout<<"NPmaxy="<<NPmaxy<<std::endl;
  Float_t NPminy=smeared_rnd_NP->GetMinimum();   std::cout<<"NPminy="<<NPminy<<std::endl;
  NPminy/=2.; std::cout<<"NPminy="<<NPminy<<std::endl;
  theory_rnd_NP->SetAxisRange(NPminy,NPmaxy,"Y");
  theory_rnd_NP->SetTitle( ((std::string)("Toy MC Incl. NP+"+orderstring+" #sigma;;"+ddxsec_yax)).c_str() );
  
  ///// check plots True/Smeared 
  TCanvas *plot_NP_true_smeared_rat = new TCanvas("plot_NP_true_smeared_rat", "plot NP true smeared rat",1400,1200);  
  plot_NP_true_smeared_rat->cd();
  
  TPad* NPoverlay_pad=new TPad("NPoverlay_pad","pad for spectra",0.0,0.3,1.0,1.0);
  NPoverlay_pad->SetLogx(1);
  NPoverlay_pad->SetLogy(1);
  NPoverlay_pad->SetBottomMargin(0);
  NPoverlay_pad->Draw();//tells the canvas theres a pad at the coordinates specified in the new constructor
  NPoverlay_pad->cd();
  
  theory_rnd_NP->DrawClone("HIST E");  
  smeared_rnd_NP->DrawClone("HIST E SAME");
  
  TLegend* NPleg_smear=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
  NPleg_smear->AddEntry(theory_rnd_NP  , ("Truth Incl. NP+"+orderstring+" #sigma"   ).c_str(), "lp" );
  NPleg_smear->AddEntry(smeared_rnd_NP , ("Smeared Incl. NP+"+orderstring+" #sigma" ).c_str(), "lp" );
  NPleg_smear->Draw();
  
  plot_NP_true_smeared_rat->cd();//go back to main canvas before doing new pad
  
  TPad* NPratio_pad=new TPad("NPratio_pad","pad for ratio",0.0,0.05,1.0,0.30);
  NPratio_pad->SetLogx(1);
  NPratio_pad->SetLogy(0);
  NPratio_pad->SetTopMargin(0);
  //NPratio_pad->SetBottomMargin(0.3);
  NPratio_pad->Draw();//tell canvas there's a pad at the specified coordinates in the new constructor
  NPratio_pad->cd();
  
  std::cout<<"creating ratio hist of true/smeared for NP+"+orderstring+""<<std::endl;
  TH1D *NP_true_smeared_rat=(TH1D*)theory_rnd_NP->Clone("NP_true_smeared_rat"); 
  if(printBaseDebug)NP_true_smeared_rat->Print("base");
  NP_true_smeared_rat->SetTitle("; Jet p_T (GeV); True/Smeared");
  NP_true_smeared_rat->Divide(NP_true_smeared_rat,smeared_rnd_NP,1.,1.,"B");     
  NP_true_smeared_rat->GetXaxis()->SetMoreLogLabels(true);
  NP_true_smeared_rat->GetXaxis()->SetNoExponent(true);
  NP_true_smeared_rat->SetAxisRange(0.6,1.4,"Y");  
  
  std::cout<<"creating ratio hist of true/opp side smeared for NP+"+orderstring+""<<std::endl;
  TH1D *NP_true_smeared_test_rat=(TH1D*)theory_rnd_NP->Clone("NP_true_smeared_test_rat"); 
  if(printBaseDebug)NP_true_smeared_test_rat->Print("base");
  NP_true_smeared_test_rat->SetTitle("; Jet p_T (GeV); True/Smeared");
  NP_true_smeared_test_rat->Divide(NP_true_smeared_test_rat,smeared_rnd_NP_test,1.,1.,"B");     
  NP_true_smeared_test_rat->GetXaxis()->SetMoreLogLabels(true);
  NP_true_smeared_test_rat->GetXaxis()->SetNoExponent(true);
  NP_true_smeared_test_rat->SetAxisRange(0.6,1.4,"Y");
  NP_true_smeared_test_rat->SetLineColor(kRed-2);  

  NP_true_smeared_rat->DrawClone("HIST E");        
  NP_true_smeared_test_rat->DrawClone("HIST E SAME");        
  
  lineatone->DrawClone();
  
  
  TCanvas* plot_response_NP_th2=new TCanvas("plot_response_NP_th2","plot response NP th2",1200, 1000);
  plot_response_NP_th2->cd()->SetLogx(1);
  plot_response_NP_th2->cd()->SetLogy(1);
  plot_response_NP_th2->cd()->SetLogz(1);
  plot_response_NP_th2->cd();
  
  response_NP_th2->SetTitle(("Toy MC Response Matrix from JER Smeared NP+"+orderstring+" #sigma_{jet}").c_str());
  response_NP_th2->GetXaxis()->SetTitle("RECO Jet p_{T} [GeV]");
  response_NP_th2->GetXaxis()->SetNoExponent(true);
  response_NP_th2->GetXaxis()->SetMoreLogLabels(true);
  
  response_NP_th2->GetYaxis()->SetTitle("GEN Jet p_{T} [GeV]");
  response_NP_th2->GetYaxis()->SetNoExponent(true);
  response_NP_th2->GetYaxis()->SetMoreLogLabels(true);
  
  setTH2_ZAxisRange(response_NP_th2);
  
  response_NP_th2->Draw("COLZ");
  //-----------------------------------------------------------------------------------------//
  //////////////////////  END production of Smeared NLO+NP spectra  ///////////////////////////
  //-----------------------------------------------------------------------------------------//

  
  
  // compare the toy mc truth w/ the original thy hist
  TCanvas *plot_true_rnd_rat = new TCanvas("plot_true_rnd_rat", "plot true rnd rat",1400,1200);  
  plot_true_rnd_rat->cd()->SetLogx(1);
  plot_true_rnd_rat->cd()->SetLogy(0);
  plot_true_rnd_rat->cd();
  std::string thy_rnd_rat_title="Toy MC v. Orig. Spectra; Jet p_{T}; Orig./Toy Spectra";
  if(useSplineWeights)
    thy_rnd_rat_title="Spline Weighted "+thy_rnd_rat_title;
  if(useFitWeights)
    thy_rnd_rat_title="Fit Weighted "+thy_rnd_rat_title;
  
  TH1D* theory_ratio_rnd=(TH1D*)theory_clone->Clone ( 
							      ( (std::string)theory->GetName() +"_debughist") .c_str()
							       );
  theory_ratio_rnd->Divide(theory_rnd);
  theory_ratio_rnd->SetLineColor(kBlack);
  theory_ratio_rnd->SetTitle(thy_rnd_rat_title.c_str());
  
  
  TH1D* theory_NP_ratio_rnd_NP=(TH1D*)theory_NP_clone->Clone ( 
									  ( (std::string)theory_NP->GetName() +"_debughist") .c_str()
									   );
  theory_NP_ratio_rnd_NP->Divide(theory_rnd_NP);
  theory_NP_ratio_rnd_NP->SetLineColor(kGray+1);
  theory_NP_ratio_rnd_NP->SetTitle(thy_rnd_rat_title.c_str());
  
  theory_ratio_rnd->SetAxisRange    (0.5,1.5,"Y");
  theory_NP_ratio_rnd_NP->SetAxisRange(0.5,1.5,"Y");
  

  std::cout<<""+orderstring+" THY/RND THY LAST BIN CONTENT = "<<theory_ratio_rnd->GetBinContent(theory_ratio_rnd->GetNbinsX()) <<std::endl;
  std::cout<<""+orderstring+" THY/RND THY LAST BIN CONTENT = "<<theory_NP_ratio_rnd_NP->GetBinContent(theory_NP_ratio_rnd_NP->GetNbinsX()) <<std::endl;
  
  
  theory_ratio_rnd->DrawClone("HIST E");
  theory_NP_ratio_rnd_NP->DrawClone("HIST E SAME");
  
  lineatone->Draw();
  
  
  
  TLegend* leg_tru_rnd_rat=new TLegend(0.10, 0.70, 0.35, 0.9, NULL,"BRNDC");
  leg_tru_rnd_rat->AddEntry( theory_ratio_rnd,       (orderstring+" Spectra"   ).c_str(), "lp");
  leg_tru_rnd_rat->AddEntry( theory_NP_ratio_rnd_NP, (orderstring+"+NP Spectra").c_str(), "lp");
  leg_tru_rnd_rat->Draw();
  
  

  


  

  ////// Write out histos
  outf->cd();
  
  // NP FITS 
  fNP->Write();

  // NP FIT CANV
  plot_NP->Write();  
  
  //JET ENERGY RESOLUTION
  fJER->Write();//FIT
  plotJER->Write();//CANV
  
  // NLO HISTS
  theory->Write();
  plot_NLOxsec->Write();

  // ERR HISTS
  theoryPDFPoserr     ->Write();
  theory2PtScalePoserr->Write();
  theory6PtScalePoserr->Write();
  plot_xsecPosErrs->Write();
  
  theoryPDFNegerr     ->Write();
  theory2PtScaleNegerr->Write();
  theory6PtScaleNegerr->Write();
  plot_xsecNegErrs->Write();
  
  // NP+NLO HISTS
  theory_NP->Write();    
  plot_NPNLOxsec->Write();
  

  // CLONES 4 SMEARING+TOYMC
  theory_clone->Write();
  theory_NP_clone->Write();
  
  
  if(useSplineWeights){  // SPLINES
    spline3->Write();
    spline3_NP->Write();
    if(useSplineExt)spline3_ext->Write();
    if(useSplineExt)spline3_NP_ext->Write();
    plot_splines->Write();  
  }
  else if(useFitWeights && useModLogFit){//FITS
    logFit->Write();
    modLogFit->Write();
    modLog2Fit->Write();
    modLog3Fit->Write();
    modLog4Fit->Write();
    modLog5Fit->Write();
    plot_logfits->Write();  
    
    logFit_NP->Write();
    modLogFit_NP->Write();
    modLog2Fit_NP->Write();
    modLog3Fit_NP->Write();
    modLog4Fit_NP->Write();
    modLog5Fit_NP->Write();
    plot_logfits_NP->Write();      
  }
  else if(useFitWeights && use7TeVFit){//FITS
    _7TeVFitA->Write();
    _7TeVFitB->Write();
    plot_7tevfits->Write();  
    
    _7TeVFitA_NP->Write();
    _7TeVFitB_NP->Write();
    plot_7tevfits_NP->Write();      
  }
  
  
  //TOY MC, JER SMEARED NLO SPECTRA
  theory_rnd->Write();
  theory_ratio_rnd->Write();
  
  smeared_rnd->Write(); 
  smeared_rnd_test->Write(); 

  true_smeared_rat->Write();  //ratio
  true_smeared_test_rat->Write();  //ratio
  plot_true_smeared_rat->Write();
  
  response_th2->Write();
  plot_response_th2->Write();  
    


  ////TOY MC NP+NLO SPECTRA, SMEARED+TRUE   
  theory_rnd_NP->Write();
  theory_NP_ratio_rnd_NP->Write();
  
  smeared_rnd_NP->Write(); 
  smeared_rnd_NP_test->Write(); 

  NP_true_smeared_rat->Write();  //ratio
  NP_true_smeared_test_rat->Write();  //ratio
  plot_NP_true_smeared_rat->Write();
  
  response_NP_th2->Write();
  plot_response_NP_th2->Write();    


  // toy spectra comparison w/ orig calculation
  plot_true_rnd_rat->Write();
  
  
  outf->Write();
  std::cout<<"output file "<<outf->GetName()<< " written!"<<std::endl;
  return 1;
  //    assert(false);    
}


//  steering ---------------------------------------------------------------------------------
int main(int argc, char* argv[]){  
  int rStatus = -1;   
  
  if(argc==6){
    rStatus=smearTheorySpectra_gaussCoreJER_etabin( 
						   (std::string) argv[1] ,  
						   (std::string) argv[2], 
						   (std::string) argv[3] ,
						   (int)std::atoi(argv[4]) ,
						   (std::string)argv[5] 						    ) ;      
  }
  
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}
  
  
















