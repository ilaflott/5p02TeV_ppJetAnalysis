#include "smearTheorySpectra.h"

//const std::string SPLINE_STR="spline";
//const std::string FIT_STR="fit";
const int CANVX=1200, CANVY=1000;
const bool printBaseDebug=true;
//const bool useSplineWeights=true;
//const bool useFitWeights=!useSplineWeights;

//const int nEvents=1e+09;  ///10x typical
//const int nEvents=1e+08;  ///typical
const int nEvents=1e+07;  /// debug nevents
//const int nEvents=1e+06;  /// debug nevents
//const int nEvents=1e+05;  /// debug nevents

std::string ddxsec_yax="#frac{d^{2}#sigma}{dp_{T}dy} [nb/GeV]";

const std::string NPCorrFile="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory/NPCorr5TeV/NLOpNP_InclusiveJets5TeV.root";
//const std::string JERCorrFile="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots_JERS/output/ak4PF_MCJEC_00eta20_08.07.18.root";
const std::string in_JERFile_default="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots_JERS/output/ak4PF_MCJEC_00eta20_08.07.18.root";

int debugInt=0;
void debugcout(){
  std::cout<<"location #"<<debugInt<<std::endl;
  debugInt++;
  return; 
}
//, std::string weightMode="spline" //"fit"// 
int smearTheorySpectra_gaussCoreJER( std::string in_NLOfileString, const bool useSplineWeights=true, const std::string fitType="modLog" , const std::string in_JERFile=in_JERFile_default
				     ){
  
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
 
  TVirtualFitter::SetDefaultFitter("Minuit2");
  //ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-02);      //default i think
  //ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-03);      
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04);      
  
  // Output File to write to
  const bool useFitWeights=!(useSplineWeights);
  const bool useModLogFit=useFitWeights && (fitType=="modLog");
  const bool use7TeVFit=useFitWeights && (fitType=="7TeV");
  if(useModLogFit && use7TeVFit){
    std::cout<<"something wrong"<<std::endl; assert(false);}

  std::string outputFile= in_NLOfileString; 
  if(useSplineWeights)    outputFile +="_spl3wgts";
  else if (useFitWeights) {
    if(useModLogFit)
      outputFile += "_modlogfitwgts";
    else if(use7TeVFit)
      outputFile += "_7tevfitwgts";
  }
  outputFile+="_gaussSmear_00eta20.root";
  
  std::cout<<"opening output:"<<outputFile<<std::endl<<std::endl;
  TFile *outf    = new TFile(outputFile.c_str(), "RECREATE" );
  
  // Input, Theory spectrum
  std::string inputFile = "fNLOJetsSpectra/R04/"+in_NLOfileString+".root";
  std::cout<<"opening input file:"<<inputFile<<std::endl<<std::endl;
  TFile* fin_NLO=TFile::Open(inputFile.c_str());
  

  /////////////// NPs to NLO jet spectra ////////////////////////////////////////////////////////  
  std::cout<<"opening NP corr file + fits! " << std::endl;
  TFile* fin_NP=TFile::Open(NPCorrFile.c_str());  
  std::string NPCorrFits_str;
  if(NPCorrFile.find("nnlo")!=std::string::npos)
    NPCorrFits_str="fNPC_POWPY8_R4_etabin";
  else//raghavs suggestion to use POWPY8 for NNLO thy.
    NPCorrFits_str="fNPC_HerwigEE4C_R4_etabin";
  std::cout<<"NPCorrFits_str = "<<NPCorrFits_str<<std::endl;  
  
  /////////////// Get NP fit functions
  TF1 *fNP_y0 = (TF1*)fin_NP->Get((NPCorrFits_str+"0").c_str());
  TF1 *fNP_y1 = (TF1*)fin_NP->Get((NPCorrFits_str+"1").c_str());
  TF1 *fNP_y2 = (TF1*)fin_NP->Get((NPCorrFits_str+"2").c_str());
  TF1 *fNP_y3 = (TF1*)fin_NP->Get((NPCorrFits_str+"3").c_str());
  if(!fin_NP ||
     !fNP_y0 || !fNP_y1 || !fNP_y2 || !fNP_y3 )
    {std::cout<<"error, NP file and/or fits not found."<<std::endl; assert(false);}
  std::string fNP_ynew_str= "(" + 
    ((std::string)fNP_y0->GetName()) + "+" + 
    ((std::string)fNP_y1->GetName()) + "+" +
    ((std::string)fNP_y2->GetName()) + "+" +
    ((std::string)fNP_y3->GetName()) +")/4.";  
  TF1 *fNP_ynew=new TF1("fNP_ynew",fNP_ynew_str.c_str() ,thyBins_incl[0],thyBins_incl[n_thybins_incl]);
  
  
  /////////////// Create/Get JER fit function(s)
  std::cout<<"opening JER file + fits! " << std::endl;
  TFile* fin_JER=TFile::Open(in_JERFile.c_str());
  if(!fin_JER){
    std::cout<<"error, JER file not found or not open."<<std::endl; assert(false); }
  
  TF1 *fJER_default = new TF1("fJER_ynew",gJER_ynew_str.c_str(),thyBins_incl[0], thyBins_incl[n_thybins_incl]);  //for comparison only
  fJER_default->SetLineColor(kBlue);
  fJER_default->SetMinimum(0.); //do i need to set min/max really?
  fJER_default->SetMaximum(0.25);
  
  TF1 *fJER_ynew = (TF1*)fin_JER->Get("SigmaFit_h");
  //  TF1 *fJER_ynew = (TF1*)fin_JER->Get("SigmaFit_f");
  if(!fJER_ynew){
    std::cout<<"error, JER fit not open or not found!."<<std::endl; assert(false);}
  fJER_ynew->SetLineColor(kRed);
  fJER_ynew->SetMinimum(0.); //do i need to set min/max really?
  fJER_ynew->SetMaximum(0.25);
  std::cout<<std::endl;
  
  



  /////////////// plots to to check NPs 
  // 2x2 canv of |y| bins
  TLine* lineatone=new TLine(thyBins_incl[0],1.,1000.,1.); lineatone->SetLineStyle(7);
  TCanvas *plot_y0_to_y3_NPs = new TCanvas("plot_y0_to_y3_NPs", "plot y0 to y3 NPs",1200,1000);
  plot_y0_to_y3_NPs->Divide(2,2);
  //-----------------------------------
  plot_y0_to_y3_NPs->cd(1)->SetLogx(1);
  plot_y0_to_y3_NPs->cd(1);
  TH1F *hNP_y0 = (TH1F*)( (TH1F*)fNP_y0->GetHistogram()
			  )->Clone("hNP_y0");
  hNP_y0->SetTitle("NP Corr Fit for 0.0 < #||{y} < 0.5;Jet p_{T};NP Corr. Factor");
  hNP_y0->GetXaxis()->SetNoExponent(true);
  hNP_y0->GetXaxis()->SetMoreLogLabels(true);
  hNP_y0->SetLineColor(kRed);
  hNP_y0->DrawClone();  
  lineatone->SetX2( hNP_y0->GetBinLowEdge(hNP_y0->GetNbinsX())+ hNP_y0->GetBinWidth(hNP_y0->GetNbinsX()) );
  lineatone->DrawClone();
  //-----------------------------------  
  plot_y0_to_y3_NPs->cd(2)->SetLogx(1);
  plot_y0_to_y3_NPs->cd(2);
  TH1F *hNP_y1 = (TH1F*)( (TH1F*)fNP_y1->GetHistogram()
			  )->Clone("hNP_y1");  
  hNP_y1->SetTitle("NP Corr Fit for 0.5 < #||{y} < 1.0;Jet p_{T};NP Corr. Factor");
  hNP_y1->GetXaxis()->SetNoExponent(true);
  hNP_y1->GetXaxis()->SetMoreLogLabels(true);
  hNP_y1->SetLineColor(kBlue);
  hNP_y1->DrawClone();  
  lineatone->SetX2( hNP_y1->GetBinLowEdge(hNP_y1->GetNbinsX())+ hNP_y1->GetBinWidth(hNP_y1->GetNbinsX()) );
  lineatone->DrawClone();
  //-----------------------------------
  plot_y0_to_y3_NPs->cd(3)->SetLogx(1);
  plot_y0_to_y3_NPs->cd(3);
  TH1F *hNP_y2 = (TH1F*)( (TH1F*)fNP_y2->GetHistogram()
			  )->Clone("hNP_y2");  
  hNP_y2->SetTitle("NP Corr Fit for 1.0 < #||{y} < 1.5;Jet p_{T};NP Corr. Factor");
  hNP_y2->GetXaxis()->SetNoExponent(true);
  hNP_y2->GetXaxis()->SetMoreLogLabels(true);
  hNP_y2->SetLineColor(kGreen);
  hNP_y2->DrawClone();  
  lineatone->SetX2( hNP_y2->GetBinLowEdge(hNP_y2->GetNbinsX())+ hNP_y2->GetBinWidth(hNP_y2->GetNbinsX()) );
  lineatone->DrawClone();
  
  //-----------------------------------
  plot_y0_to_y3_NPs->cd(4)->SetLogx(1);
  plot_y0_to_y3_NPs->cd(4);
  TH1F *hNP_y3 = (TH1F*)( (TH1F*)fNP_y3->GetHistogram()
			  )->Clone("hNP_y3");  
  hNP_y3->SetTitle("NP Corr Fit for 1.5 < #||{y} < 2.0;Jet p_{T};NP Corr. Factor");
  hNP_y3->GetXaxis()->SetNoExponent(true);
  hNP_y3->GetXaxis()->SetMoreLogLabels(true);
  hNP_y3->SetLineColor(kMagenta);
  hNP_y3->DrawClone();    
  lineatone->SetX2( hNP_y3->GetBinLowEdge(hNP_y3->GetNbinsX())+ hNP_y3->GetBinWidth(hNP_y3->GetNbinsX()) );
  lineatone->DrawClone();  
  //-----------------------------------
  

  // |y|<2.0
  TCanvas *plot_totNPs = new TCanvas("plot_totNPs", "plot totNPs",600,500);//1600,1200);
  plot_totNPs->SetLogx(1);
  plot_totNPs->SetLogy(0);
  plot_totNPs->cd();
  //-----------------------------------
  TH1F *hNP_ynew = (TH1F*)( (TH1F*)fNP_ynew->GetHistogram()
			    )->Clone("hNP_ynew");  
  hNP_ynew->SetTitle("Avg. NP Corr Fit for 0.0 < #||{y} < 2.0;Jet p_{T};NP Corr. Factor");
  hNP_ynew->GetXaxis()->SetNoExponent(true);
  hNP_ynew->GetXaxis()->SetMoreLogLabels(true);
  hNP_ynew->SetLineColor(kTeal);
  hNP_ynew->DrawClone();  //plot_totNPs->Update();  //why did this need to happen?
  lineatone->SetX2( hNP_ynew->GetBinLowEdge(hNP_ynew->GetNbinsX())+ hNP_ynew->GetBinWidth(hNP_ynew->GetNbinsX()) );
  lineatone->DrawClone();  
  //-----------------------------------  
  
  
  // |y|<2.0 + 0.5 size |y| bins on same canv
  TCanvas *plot_allNPs = new TCanvas("plot_allNPs", "plot allNPs",600,500);
  plot_allNPs->SetLogx(1);
  plot_allNPs->cd();
  //-----------------------------------  
  float NPfit_max=hNP_y1->GetMaximum();
  float NPfit_min=hNP_y3->GetMinimum();
  hNP_y0->SetAxisRange(NPfit_min/1.02,NPfit_max*1.02,"Y");
  //hNP_ynew->SetTitle("NP fits all #||{y} bins");
  hNP_y0->SetTitle("NP Fits, all #||{y} bins;Jet p_{T};NP Corr. Factor");
  hNP_y0  ->DrawClone("");
  hNP_y1  ->DrawClone("SAME");
  hNP_y2  ->DrawClone("SAME");
  hNP_y3  ->DrawClone("SAME");
  //hNP_ynew->DrawClone("SAME");
  //-----------------------------------    
  lineatone->SetX2( 1000.);//hNP_y0->GetBinLowEdge(hNP_y0->GetNbinsX())+ hNP_y0->GetBinWidth(hNP_y0->GetNbinsX()) );
  lineatone->DrawClone();
  //-----------------------------------    
  //TLegend* leg_allNPs=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
  TLegend* leg_allNPs=new TLegend(0.10, 0.70, 0.35, 0.9, NULL,"BRNDC");
  //leg_allNPs->AddEntry(hNP_ynew , "#||{y} < 2.0 (avg of others)" , "l");
  leg_allNPs->AddEntry(hNP_y0   , "0.0 < #||{y} < 0.5" , "l");
  leg_allNPs->AddEntry(hNP_y1   , "0.5 < #||{y} < 1.0" , "l");
  leg_allNPs->AddEntry(hNP_y2   , "1.0 < #||{y} < 1.5" , "l");
  leg_allNPs->AddEntry(hNP_y3   , "1.5 < #||{y} < 2.0" , "l");
  leg_allNPs->Draw();
  //-----------------------------------  
  hNP_y0->SetTitle("NP Corr Fit for 0.0 < #||{y} < 0.5;Jet p_{T};NP Corr. Factor");//undo the titling because this is for this plot only
  //-----------------------------------    
  
  
  
  /////////////// plot Core p_T Resolution ////////////////////////////////////////////////////////
  std::cout<<"plotting JER!"<<  std::endl;  
  TH1D *hJER_ynew = (TH1D*)(
			    (TH1D*)fJER_ynew->GetHistogram()
			    )->Clone("hJER_ynew");//TF1 is irritating to set titles for+blah... so use this
  hJER_ynew->SetTitle(" Gauss Core JER Fit #||{y} < 2.0 ; Jet p_{T} ; #sigma / #mu from Fit");//sets hist title+xaxis+yaxis title
  hJER_ynew->GetXaxis()->SetNoExponent(true);
  hJER_ynew->GetXaxis()->SetMoreLogLabels(true);
  TH1D *hJER_default = (TH1D*)(
			       (TH1D*)fJER_default->GetHistogram()
			       )->Clone("hJER_default");//for comparison only
  hJER_default->SetTitle(" Gauss Core JER Fit #||{y} < 2.0 ; Jet p_{T} ; #sigma / #mu from Fit");//sets hist title+xaxis+yaxis title
  hJER_default->GetXaxis()->SetNoExponent(true);
  hJER_default->GetXaxis()->SetMoreLogLabels(true);
  //-----------------------------------  
  TCanvas *plotJER = new TCanvas("plotJER", "plot JER",900,600);
  plotJER->cd()->SetLogx(1);
  plotJER->cd();  
  //hJER_default->DrawClone("HIST E");    
  hJER_ynew->DrawClone("HIST E");    
  //-----------------------------------
  //TO DO REPLACE ME W/TEXT SAYING SAME THING
  //TLegend *leg=new TLegend(0.53,0.50,0.85,0.85);
  //leg->AddEntry(hJER_ynew, "PY8 MC CUETP8M1 JER Fits New");
  //leg->AddEntry(hJER_default, "Default JER Fits from APS DNP 10/3/17");
  //leg->AddEntry(fJER_ynew,"PYTHIA8 CUETP8M1 @ 5.02 TeV", "");
  //leg->AddEntry(fJER_ynew,"ak4PFJets", "");
  //leg->AddEntry(fJER_ynew,"0<|y|<2.0", "l"); 
  //leg->SetTextFont(42);
  //leg->SetFillColor(kWhite);
  //leg->Draw();
  //-----------------------------------

  
  
  //////////// NLO Thy calculation xsecions get from file ////////////////////////////////////////////////////////
  std::cout<<"plotting NLO xsecs!"<<  std::endl;  
  TH1D *theory_y0 = (TH1D*)fin_NLO->Get("h0100100");
  theory_y0->Sumw2(true);
  if(printBaseDebug)  theory_y0->Print("Base");
  theory_y0->Scale(1.e-03);//picobarns to nanobarns
  TH1D *theory_y1 = (TH1D*)fin_NLO->Get("h0100200");
  theory_y1->Sumw2(true);
  if(printBaseDebug)  theory_y1->Print("Base");
  theory_y1->Scale(1.e-03);
  TH1D *theory_y2 = (TH1D*)fin_NLO->Get("h0100300");
  theory_y2->Sumw2(true);
  if(printBaseDebug)  theory_y2->Print("Base");
  theory_y2->Scale(1.e-03);
  TH1D *theory_y3 = (TH1D*)fin_NLO->Get("h0100400");
  theory_y3->Sumw2(true);
  if(printBaseDebug)  theory_y3->Print("Base");
  theory_y3->Scale(1.e-03);
  std::cout<<std::endl;
  //-----------------------------------
  std::cout <<"making 00eta20 hist from theory curves..."<<std::endl<<std::endl;
  TH1D *theory_ynew =    make00eta20Hist( (TH1D*) theory_y0, 
					  (TH1D*) theory_y1, 
					  (TH1D*) theory_y2, 
					  (TH1D*) theory_y3,
					  "h00eta20");  
  theory_ynew->SetTitle(("NLO Inclusive #sigma_{jet}, #||{y} < 2.0;Jet p_{T};"+ddxsec_yax).c_str());
  theory_ynew->GetXaxis()->SetNoExponent(true);
  theory_ynew->GetXaxis()->SetMoreLogLabels(true);
  theory_ynew->SetLineColor(kTeal);
  std::cout<<"done making hist."<<std::endl;
  if(printBaseDebug)theory_ynew->Print("base");  
  //-----------------------------------
  

  /////////////// draw xsec hists from file
  //2X2 CANV
  TCanvas *plot_y0_to_y3_NLOxsec = new TCanvas("plot_y0_to_y3_NLOxsec", "plot y0 to y3 NLOxsec",1200,1000);
  plot_y0_to_y3_NLOxsec->Divide(2,2);
  //--------------------------------------------
  plot_y0_to_y3_NLOxsec->cd(1)->SetLogx(1);
  plot_y0_to_y3_NLOxsec->cd(1)->SetLogy(1);
  plot_y0_to_y3_NLOxsec->cd(1);
  theory_y0->SetTitle(("NLO #sigma_{jet}, 0.0 < #||{y} < 0.5;Jet p_{T};"+ddxsec_yax).c_str());
  theory_y0->GetXaxis()->SetNoExponent(true);
  theory_y0->GetXaxis()->SetMoreLogLabels(true);
  theory_y0->SetLineColor(kRed);
  theory_y0->DrawClone("HIST E");    
  //--------------------------------------------
  plot_y0_to_y3_NLOxsec->cd(2)->SetLogx(1);
  plot_y0_to_y3_NLOxsec->cd(2)->SetLogy(1);
  plot_y0_to_y3_NLOxsec->cd(2);
  theory_y1->SetTitle(("NLO #sigma_{jet}, 0.5 < #||{y} < 1.0;Jet p_{T};"+ddxsec_yax).c_str());
  theory_y1->GetXaxis()->SetNoExponent(true);
  theory_y1->GetXaxis()->SetMoreLogLabels(true);
  theory_y1->SetLineColor(kBlue);
  theory_y1->DrawClone("HIST E");    
  //--------------------------------------------
  plot_y0_to_y3_NLOxsec->cd(3)->SetLogx(1);
  plot_y0_to_y3_NLOxsec->cd(3)->SetLogy(1);
  plot_y0_to_y3_NLOxsec->cd(3);
  theory_y2->SetTitle(("NLO #sigma_{jet}, 1.0 < #||{y} < 1.5;Jet p_{T};"+ddxsec_yax).c_str());
  theory_y2->GetXaxis()->SetNoExponent(true);
  theory_y2->GetXaxis()->SetMoreLogLabels(true);
  theory_y2->SetLineColor(kGreen);
  theory_y2->DrawClone("HIST E");    
  //--------------------------------------------
  plot_y0_to_y3_NLOxsec->cd(4)->SetLogx(1);
  plot_y0_to_y3_NLOxsec->cd(4)->SetLogy(1);
  plot_y0_to_y3_NLOxsec->cd(4);
  theory_y3->SetTitle(("NLO #sigma_{jet}, 1.5 < #||{y} < 2.0;Jet p_{T};"+ddxsec_yax).c_str());
  theory_y3->GetXaxis()->SetNoExponent(true);
  theory_y3->GetXaxis()->SetMoreLogLabels(true);
  theory_y3->SetLineColor(kMagenta);
  theory_y3->DrawClone("HIST E");    
  //--------------------------------------------
  

  // total NLO xsec (avg of other 4 [kinda])
  TCanvas *plot_totNLOxsec = new TCanvas("plot_totNLOxsec", "plot totNLOxsec",600,500);//1600,1200);
  plot_totNLOxsec->SetLogx(1);
  plot_totNLOxsec->SetLogy(1);
  plot_totNLOxsec->cd();
  theory_ynew->SetTitle("NLO Inclusive #sigma_{jet} #||{y} < 2.0");
  theory_ynew->DrawClone("HIST E");
  //--------------------------------------------
  
  // plot of all xsecs on one hist
  TCanvas *plot_allNLOxsec= new TCanvas("plot_allNLOxsec", "plot allNLOxsec",900,750);
  plot_allNLOxsec->SetLogx(1);
  plot_allNLOxsec->SetLogy(1);
  plot_allNLOxsec->cd();
  //--------------------------------------------
  theory_y0->Scale(1.e+03);  
  theory_y1->Scale(1.e+02);
  theory_y2->Scale(1.e+01);
  theory_y3->Scale(1.e+00);                                                                                      
  //--------------------------------------------  
  theory_y0->SetTitle("NLO Inclusive #sigma_{jet} all #||{y} bins"); 
  Double_t max=(Double_t)theory_y0->GetMaximum();
  Double_t min=(Double_t)theory_y3->GetMinimum();
  theory_y0->SetAxisRange((min/2.),(max*2.),"Y");
  theory_y0->DrawClone("HIST E");
  theory_y0->Print("base");
  theory_y1->DrawClone("HIST E SAME");
  theory_y2->DrawClone("HIST E SAME");
  theory_y3->DrawClone("HIST E SAME");  
  //-------------------------------------------- 
  TLegend* leg_allNLOxsec=new TLegend(0.60, 0.70, 0.9, 0.9, NULL,"BRNDC");
  leg_allNLOxsec->AddEntry(theory_y0   , "(0.0 < #||{y} < 0.5) x 10^{3}" , "lp");
  leg_allNLOxsec->AddEntry(theory_y1   , "(0.5 < #||{y} < 1.0) x 10^{2}" , "lp");
  leg_allNLOxsec->AddEntry(theory_y2   , "(1.0 < #||{y} < 1.5) x 10^{1}" , "lp");
  leg_allNLOxsec->AddEntry(theory_y3   , "(1.5 < #||{y} < 2.0) x 10^{0}" , "lp");
  leg_allNLOxsec->Draw();
  //--------------------------------------------  
  //unscale the NLO xsec histograms for later use...
  theory_y0->Scale(  1.e-03);  
  theory_y1->Scale(  1.e-02);
  theory_y2->Scale(  1.e-01);
  theory_y3->Scale(  1.e-00);                                                                                      
  theory_y0->SetTitle(("NLO #sigma_{jet}, 0.0 < #||{y} < 0.5;Jet p_{T};"+ddxsec_yax).c_str());
  //--------------------------------------------
  
  
  
  
  
  //////////// Apply NPs to theory + Plot////////////////////////////////////////////////////////
  std::cout<<"plotting NP+NLO xsecs!"<<  std::endl;  
  TH1D *theory_NPy0  = (TH1D*) applyNPtoxsec(theory_y0, fNP_y0);
  TH1D *theory_NPy1  = (TH1D*) applyNPtoxsec(theory_y1, fNP_y1);
  TH1D *theory_NPy2  = (TH1D*) applyNPtoxsec(theory_y2, fNP_y2);
  TH1D *theory_NPy3  = (TH1D*) applyNPtoxsec(theory_y3, fNP_y3);
  TH1D *theory_NPynew= (TH1D*) make00eta20Hist( (TH1D*) theory_NPy0, 
						(TH1D*) theory_NPy1, 
						(TH1D*) theory_NPy2, 
						(TH1D*) theory_NPy3,
						"h00eta20_wNP");  
  theory_NPynew->SetTitle(("NP+NLO Inclusive #sigma_{jet}, #||{y} < 2.0;Jet p_{T};"+ddxsec_yax).c_str());
  theory_NPynew->GetXaxis()->SetNoExponent(true);
  theory_NPynew->GetXaxis()->SetMoreLogLabels(true);
  theory_NPynew->SetLineColor(kTeal);
  

  //2X2 CANV THEORY NPNLO xsec
  TCanvas *plot_y0_to_y3_NPNLOxsec = new TCanvas("plot_y0_to_y3_NPNLOxsec", "plot y0 to y3 NPNLOxsec",1200,1000);
  plot_y0_to_y3_NPNLOxsec->Divide(2,2);
  //--------------------------------------------
  plot_y0_to_y3_NPNLOxsec->cd(1)->SetLogx(1);
  plot_y0_to_y3_NPNLOxsec->cd(1)->SetLogy(1);
  plot_y0_to_y3_NPNLOxsec->cd(1);
  theory_NPy0->DrawClone("HIST E");    
  //--------------------------------------------
  plot_y0_to_y3_NPNLOxsec->cd(2)->SetLogx(1);
  plot_y0_to_y3_NPNLOxsec->cd(2)->SetLogy(1);
  plot_y0_to_y3_NPNLOxsec->cd(2);
  theory_NPy1->DrawClone("HIST E");    
  //--------------------------------------------
  plot_y0_to_y3_NPNLOxsec->cd(3)->SetLogx(1);
  plot_y0_to_y3_NPNLOxsec->cd(3)->SetLogy(1);
  plot_y0_to_y3_NPNLOxsec->cd(3);
  theory_NPy2->DrawClone("HIST E");    
  //--------------------------------------------
  plot_y0_to_y3_NPNLOxsec->cd(4)->SetLogx(1);
  plot_y0_to_y3_NPNLOxsec->cd(4)->SetLogy(1);
  plot_y0_to_y3_NPNLOxsec->cd(4);
  theory_NPy3->DrawClone("HIST E");    
  //--------------------------------------------

  // total NLO xsec (avg of other 4 [kinda])
  TCanvas *plot_totNPNLOxsec = new TCanvas("plot_totNPNLOxsec", "plot totNPNLOxsec",600,500);//1600,1200);
  plot_totNPNLOxsec->SetLogx(1);
  plot_totNPNLOxsec->SetLogy(1);
  plot_totNPNLOxsec->cd();
  theory_NPynew->SetTitle("NP+NLO Inclusive #sigma_{jet} #||{y} < 2.0");
  theory_NPynew->DrawClone("HIST E");
  //--------------------------------------------
  
  // plot of all xsecs on one hist
  TCanvas *plot_allNPNLOxsec= new TCanvas("plot_allNPNLOxsec", "plot allNPNLOxsec",900,750);
  plot_allNPNLOxsec->SetLogx(1);
  plot_allNPNLOxsec->SetLogy(1);
  plot_allNPNLOxsec->cd();
  //--------------------------------------------
  theory_NPy0->Scale(1.e+03);  
  theory_NPy1->Scale(1.e+02);
  theory_NPy2->Scale(1.e+01);
  theory_NPy3->Scale(1.e+00);                                                                                      
  //--------------------------------------------  
  theory_NPy0->SetTitle("NP+NLO Inclusive #sigma_{jet}, all #||{y} bins");
  Double_t NPmax=(Double_t)theory_NPy0->GetMaximum();
  Double_t NPmin=(Double_t)theory_NPy3->GetMinimum();
  theory_NPy0->SetAxisRange((NPmin/2.),(NPmax*2.),"Y");
  theory_NPy0->DrawClone("HIST E");
  theory_NPy0->Print("base");
  theory_NPy1->DrawClone("HIST E SAME");
  theory_NPy2->DrawClone("HIST E SAME");
  theory_NPy3->DrawClone("HIST E SAME");  
  //-------------------------------------------- 
  TLegend* leg_allNPNLOxsec=new TLegend(0.60, 0.70, 0.9, 0.9, NULL,"BRNDC");
  leg_allNPNLOxsec->AddEntry(theory_NPy0   , "(0.0 < #||{y} < 0.5) x 10^{3}" , "lp");
  leg_allNPNLOxsec->AddEntry(theory_NPy1   , "(0.5 < #||{y} < 1.0) x 10^{2}" , "lp");
  leg_allNPNLOxsec->AddEntry(theory_NPy2   , "(1.0 < #||{y} < 1.5) x 10^{1}" , "lp");
  leg_allNPNLOxsec->AddEntry(theory_NPy3   , "(1.5 < #||{y} < 2.0) x 10^{0}" , "lp");
  leg_allNPNLOxsec->Draw();
  //--------------------------------------------  
  //unscale the NLO xsec histograms for later use...
  theory_NPy0->Scale(  1.e-03);  
  theory_NPy1->Scale(  1.e-02);
  theory_NPy2->Scale(  1.e-01);
  theory_NPy3->Scale(  1.e-00);                            

  theory_NPy0->SetTitle(("NP+NLO #sigma_{jet}, 0.0 < #||{y} < 0.5;Jet p_{T};"+ddxsec_yax).c_str());
  //--------------------------------------------





  

  //////////////// Use clones for smearing weights (fit or spline)
  //TH1D* theory_ynew_clone=(TH1D*)theory_ynew->TH1::Rebin(n_thybins_incl, 
  //							 ( ( (std::string) theory_ynew->GetName() ) +"_specialrebin").c_str() ,	 
  //							 thyBins_incl );  theory_ynew_clone->SetLineColor(kCyan+4);    
  //
  //TH1D* theory_NPynew_clone=(TH1D*)theory_NPynew->TH1::Rebin(n_thybins_incl, 
  //							     ( ( (std::string) theory_NPynew->GetName() ) +"_specialrebin").c_str() ,	 
  //							     thyBins_incl );    theory_NPynew_clone->SetLineColor(kCyan-6);  
  TH1D* theory_ynew_clone=(TH1D*)theory_ynew->Clone(( ( (std::string) theory_ynew->GetName() ) +"_0").c_str());
  theory_ynew_clone=(TH1D*)theory_ynew_clone->TH1::Rebin(n_thybins_incl, ( ( (std::string) theory_ynew_clone->GetName() ) +"_rebin").c_str() , thyBins_incl );  
  theory_ynew_clone->SetLineColor(kCyan+4);    

  TH1D* theory_NPynew_clone=(TH1D*)theory_NPynew->Clone(( ( (std::string) theory_NPynew->GetName() ) +"_0").c_str());
  theory_NPynew_clone=(TH1D*)theory_NPynew_clone->TH1::Rebin(n_thybins_incl, ( ( (std::string) theory_NPynew_clone->GetName() ) +"_rebin").c_str() , thyBins_incl );  
  theory_NPynew_clone->SetLineColor(kCyan+4);    
  
  
  ///Cubic Spline Interpolations using Cross sections; has issues for p_T beyond first/last bin centers
  TSpline3 *spline3_ynew=NULL, *spline3_NPynew=NULL;
  TF1      *spline3_ynew_ext=NULL, *spline3_NPynew_ext=NULL;
  TCanvas *plot_splines=NULL;
  TLegend *leg_spline=NULL;
  if(useSplineWeights){
    
    /// check spectra with splines
    std::cout<<"drawing Tspline and cross section on canvas"<<std::endl;
    plot_splines = new TCanvas("plot_splines", "plot splines",1200,800);
    plot_splines->cd()->SetLogx(1);
    plot_splines->cd()->SetLogy(1);
    plot_splines->cd();
    
    // NLO //
    std::cout<<"creating TSpline of cross section hist"<<std::endl;  
    spline3_ynew = new TSpline3( theory_ynew_clone);
    spline3_ynew->SetName( ( (std::string)theory_ynew_clone->GetName() + "_spline3").c_str() );
    spline3_ynew->SetLineColor(kAzure);  
    
    //figure out starting values for parabolic extension of spline
    //double x1=theory_ynew_clone->GetBinCenter(theory_ynew_clone->GetNbinsX());// bad idea, spline3's val at x1 is approx the val of the last bin's height. 
    //double x1=theory_ynew_clone->GetBinLowEdge(theory_ynew_clone->GetNbinsX()); //v1
    double x1=theory_ynew_clone->GetBinCenter(theory_ynew_clone->GetNbinsX()-1);//v2
    //double x1=theory_ynew_clone->GetBinLowEdge(theory_ynew_clone->GetNbinsX()-1); //v3
    double x2=theory_ynew_clone->GetBinLowEdge(theory_ynew_clone->GetNbinsX()) + theory_ynew_clone->GetBinWidth(theory_ynew_clone->GetNbinsX());
    
    
    
    //double ynew_y2=theory_ynew_clone->GetBinContent(theory_ynew_clone->GetNbinsX());
    double ynew_y2=calc_spline3ext_y2(theory_ynew_clone);
    //double xprime=theory_ynew_clone->GetBinCenter(theory_ynew_clone->GetNbinsX());
    //double xprime=theory_ynew_clone->GetBinCenter(theory_ynew_clone->GetNbinsX()) - 0.25*theory_ynew_clone->GetBinWidth(theory_ynew_clone->GetNbinsX());
    //double ynew_yprime=theory_ynew_clone->GetBinContent(theory_ynew_clone->GetNbinsX());
    
    spline3_ynew_ext= new TF1( ((std::string)spline3_ynew->GetName()+"_tf1ext").c_str(), "[0]*pow(x,2)+[1]*x+[2]", x1, x2);    
    make_spline3_ext(  (TSpline3*)spline3_ynew, (TF1*)spline3_ynew_ext , (double)x1, (double)x2, (double)ynew_y2);        //make_spline3_extv2(  (TSpline3*)spline3_ynew, (TF1*)spline3_ynew_ext , (double)x1, (double)x2, (double)ynew_y2, (double)xprime, (double)ynew_yprime);    
    spline3_ynew_ext->SetLineColor(kAzure-1);
    
    
    // NLO+NP //      
    spline3_NPynew = new TSpline3( theory_NPynew_clone);
    spline3_NPynew->SetName( ( (std::string)theory_NPynew_clone->GetName() + "_spline3").c_str() );
    spline3_NPynew->SetLineColor(kAzure-8);  
    
    //double NPynew_y2=theory_NPynew_clone->GetBinContent(theory_NPynew_clone->GetNbinsX());
    double NPynew_y2=calc_spline3ext_y2(theory_NPynew_clone);
    
    spline3_NPynew_ext= new TF1( ((std::string)spline3_NPynew->GetName()+"_tf1ext").c_str(), "[0]*pow(x,2)+[1]*x+[2]", x1, x2);    
    make_spline3_ext(  (TSpline3*)spline3_NPynew, (TF1*)spline3_NPynew_ext , (double)x1, (double)x2, (double)NPynew_y2);    
    spline3_NPynew_ext->SetLineColor(kAzure-9);
    
    //draw
    theory_ynew_clone->SetTitle("Cubic Splines for #||{y}<2.0;Jet p_{T};Smear Weight");    
    theory_ynew_clone->SetAxisRange( ynew_y2/2., theory_ynew_clone->GetMaximum()*2., "Y");
    
    theory_ynew_clone->DrawClone("HIST E");
    spline3_ynew->Draw("SAME");    
    spline3_ynew_ext->Draw("SAME");
    
    theory_NPynew_clone->DrawClone("HIST E SAME");  
    spline3_NPynew->Draw("SAME");
    spline3_NPynew_ext->Draw("SAME");
    
    leg_spline=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
    leg_spline->AddEntry(theory_ynew_clone , "Weighted NLO Jet Counts for #||{y} < 2.0" , "lp");
    leg_spline->AddEntry(spline3_ynew   , "Cubic Spline for NLO #||{y} < 2.0" , "l");
    leg_spline->AddEntry(spline3_ynew_ext  , "Parabolic Spline Extension for NLO #||{y} < 2.0" , "l");
    leg_spline->AddEntry(theory_NPynew_clone , "Weighted NP+NLO Jet Counts for #||{y} < 2.0" , "lp");
    leg_spline->AddEntry(spline3_NPynew   , "Cubic Spline for NP+NLO #||{y} < 2.0" , "l");
    leg_spline->AddEntry(spline3_NPynew_ext  , "Parabolic Spline Extension for NP+NLO #||{y} < 2.0" , "l");
    leg_spline->Draw();    
  }

  // for fits
  const int startbin=1, endbin=theory_ynew_clone->GetNbinsX();//ynew and NP ynew should always have the same binning
  //const float xlo=theory_ynew_clone->GetBinLowEdge(startbin);
  //const float xhi=theory_ynew_clone->GetBinLowEdge(endbin) + theory_ynew_clone->GetBinWidth(endbin);  
  const float xlo=theory_ynew_clone->GetBinLowEdge(startbin);//trying out for the 7TeV style fit
  const float xhi=theory_ynew_clone->GetBinCenter(endbin);// + theory_ynew_clone->GetBinWidth(endbin);  
  
  
    
  ///Log fit: A(B/pt)^[C+ D Log(pt/B) + E(Log(pt/B))^2] using Cross sections; currently no known issues
  TF1 *logFit_ynew  =NULL, *modLogFit_ynew  =NULL, *modLog2Fit_ynew  =NULL, *modLog3Fit_ynew  =NULL, *modLog4Fit_ynew  =NULL,*modLog5Fit_ynew  =NULL;
  TF1 *logFit_NPynew=NULL, *modLogFit_NPynew=NULL, *modLog2Fit_NPynew=NULL, *modLog3Fit_NPynew=NULL, *modLog4Fit_NPynew=NULL,*modLog5Fit_NPynew=NULL;
  
  TCanvas *plot_logfits_ynew=NULL,*plot_logfits_NPynew=NULL;
  TLegend *leg_logfits_ynew=NULL,*leg_logfits_NPynew=NULL;  
  //  if(useFitWeights){
  if(useFitWeights && useModLogFit){
    
    // --- NLO --- //
    /// check spectra with fits
    plot_logfits_ynew = new TCanvas("plot_logfits_ynew", "plot logfits ynew",1200,800);
    plot_logfits_ynew->cd()->SetLogx(1);
    plot_logfits_ynew->cd()->SetLogy(1);
    plot_logfits_ynew->cd();
    
    std::cout<<"fitting hist theory_ynew_clone"<<std::endl;
    logFit_ynew    =new TF1("logfit_ynew    ","[0] * ( TMath::Power( ( [1]/x ) ,   [2] ) )",xlo,xhi);
    modLogFit_ynew =new TF1("modlogfit_ynew ","[0] * ( TMath::Power( ( [1]/x ) ,   [2] + [3]*TMath::Log10( x/[1] ) ) )",xlo,xhi);
    modLog2Fit_ynew=new TF1("modlog2fit_ynew","[0] * ( TMath::Power( ( [1]/x ) , ( [2] + [3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )  ) ) )",xlo,xhi);
    modLog3Fit_ynew=new TF1("modlog3fit_ynew","[0] * ( TMath::Power( ( [1]/x ) , ( [2] + [3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )  + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 )  ) ) )",xlo,xhi);
    modLog4Fit_ynew=new TF1("modlog4fit_ynew","[0] * ( TMath::Power( ( [1]/x ) , ( [2] + [3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )  + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 ) + [6]*TMath::Power( TMath::Log10( x/[1] ), 4 ) )               ) )",xlo,xhi);
    modLog5Fit_ynew=new TF1("modlog5fit_ynew","[0] * ( TMath::Power( ( [1]/x ) , ( [2] + [3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )  + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 ) + [6]*TMath::Power( TMath::Log10( x/[1] ), 4 ) + [7]*TMath::Power( TMath::Log10( x/[1] ), 5 ) )               ) )",xlo,xhi);

    
    makeModLogFits(theory_ynew_clone, logFit_ynew, modLogFit_ynew, modLog2Fit_ynew, modLog3Fit_ynew, modLog4Fit_ynew, modLog5Fit_ynew);
    
    



    //draw
    std::cout<<"drawing ynew fits and cross section on canvas"<<std::endl;
    theory_ynew_clone->SetTitle("Modified Logarithm Fits for NLO #||{y}<2.0;Jet p_{T};Smear Weight");    
    theory_ynew_clone->DrawClone("HIST E");    
    //theory_ynew->SetTitle("Modified Logarithm Fits for NLO #||{y}<2.0;Jet p_{T};Smear Weight");    
    //theory_ynew->DrawClone("HIST E");    
    
    logFit_ynew->Draw("SAME");
    modLogFit_ynew->Draw("SAME");
    modLog2Fit_ynew->Draw("SAME");
    modLog3Fit_ynew->Draw("SAME");    
    modLog4Fit_ynew->Draw("SAME");    
    modLog5Fit_ynew->Draw("SAME");    
    
    leg_logfits_ynew=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
    
    leg_logfits_ynew->AddEntry(theory_ynew_clone , "Weighted NLO Jet Counts for #||{y} < 2.0" , "lp");//change desc TO DO
    //    leg_logfits_ynew->AddEntry(theory_ynew , "Weighted NLO Jet Counts for #||{y} < 2.0" , "lp");//change desc TO DO
    leg_logfits_ynew->AddEntry(logFit_ynew      , "Log Fit" , "l");//change desc TO DO
    leg_logfits_ynew->AddEntry(modLogFit_ynew   , "Mod. Log Fit" , "l");
    leg_logfits_ynew->AddEntry(modLog2Fit_ynew   , "2nd Order Mod. Log Fit" , "l");    
    leg_logfits_ynew->AddEntry(modLog3Fit_ynew   , "3rd Order Mod. Log Fit" , "l");    
    leg_logfits_ynew->AddEntry(modLog4Fit_ynew   , "4th Order Mod. Log Fit" , "l");    
    leg_logfits_ynew->AddEntry(modLog5Fit_ynew   , "5th Order Mod. Log Fit" , "l");    
    leg_logfits_ynew->Draw();    
    
    
    // --- NLO+NP --- //
    /// check spectra with fits    
    plot_logfits_NPynew = new TCanvas("plot_logfits_NPynew", "plot logfits NPynew",1200,800);
    plot_logfits_NPynew->cd()->SetLogx(1);
    plot_logfits_NPynew->cd()->SetLogy(1);
    plot_logfits_NPynew->cd();
    
    std::cout<<"fitting hist theory_NPynew_clone"<<std::endl;
    logFit_NPynew    =new TF1("logfit_NPynew    ","[0] * ( TMath::Power( ( [1]/x ) , ( [2]                                                                                                                         ) ) )",xlo,xhi);
    modLogFit_NPynew =new TF1("modlogfit_NPynew ","[0] * ( TMath::Power( ( [1]/x ) , ( [2]+[3]*TMath::Log10( x/[1] )                                                                                               ) ) )",xlo,xhi);
    modLog2Fit_NPynew=new TF1("modlog2fit_NPynew","[0] * ( TMath::Power( ( [1]/x ) , ( [2]+[3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )                                                ) ) )",xlo,xhi);
    modLog3Fit_NPynew=new TF1("modlog3fit_NPynew","[0] * ( TMath::Power( ( [1]/x ) , ( [2]+[3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 ) + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 ) ) ) )",xlo,xhi);
    modLog4Fit_NPynew=new TF1("modlog4fit_NPynew","[0] * ( TMath::Power( ( [1]/x ) , ( [2]+[3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 ) + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 ) + [6]*TMath::Power( TMath::Log10( x/[1] ), 4 ) ) ) )",xlo,xhi);
    modLog5Fit_NPynew=new TF1("modlog5fit_NPynew","[0] * ( TMath::Power( ( [1]/x ) , ( [2]+[3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 ) + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 ) + [6]*TMath::Power( TMath::Log10( x/[1] ), 4 ) + [7]*TMath::Power( TMath::Log10( x/[1] ), 5 ) ) ) )",xlo,xhi); 
    
    
    makeModLogFits(theory_NPynew_clone, logFit_NPynew, modLogFit_NPynew, modLog2Fit_NPynew, modLog3Fit_NPynew, modLog4Fit_NPynew, modLog5Fit_NPynew);
    
    
    //draw
    std::cout<<"drawing NPynew fits and cross section on canvas"<<std::endl;
    theory_NPynew_clone->SetTitle("Modified Logarithm Fits for NLO+NP #||{y}<2.0;Jet p_{T};Smear Weight");        
    theory_NPynew_clone->DrawClone("HIST E");        
    //theory_NPynew->SetTitle("Modified Logarithm Fits for NLO+NP #||{y}<2.0;Jet p_{T};Smear Weight");        
    //theory_NPynew->DrawClone("HIST E");        
    
    logFit_NPynew->Draw("SAME");
    modLogFit_NPynew->Draw("SAME");
    modLog2Fit_NPynew->Draw("SAME");
    modLog3Fit_NPynew->Draw("SAME");    
    modLog4Fit_NPynew->Draw("SAME");    
    modLog5Fit_NPynew->Draw("SAME");    
    
    leg_logfits_NPynew=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
    leg_logfits_NPynew->AddEntry(theory_NPynew_clone , "Weighted NLO+NP Jet Counts for #||{y} < 2.0" , "lp");//change desc TO DO
    //leg_logfits_NPynew->AddEntry(theory_NPynew , "Weighted NLO+NP Jet Counts for #||{y} < 2.0" , "lp");//change desc TO DO
    leg_logfits_NPynew->AddEntry(logFit_NPynew      , "Log Fit" , "l");//change desc TO DO
    leg_logfits_NPynew->AddEntry(modLogFit_NPynew   , "Mod. Log Fit" , "l");
    leg_logfits_NPynew->AddEntry(modLog2Fit_NPynew   , "2nd Order Mod. Log Fit" , "l");    
    leg_logfits_NPynew->AddEntry(modLog3Fit_NPynew   , "3rd Order Mod. Log Fit" , "l");    
    leg_logfits_NPynew->AddEntry(modLog4Fit_NPynew   , "4th Order Mod. Log Fit" , "l");    
    leg_logfits_NPynew->AddEntry(modLog5Fit_NPynew   , "5th Order Mod. Log Fit" , "l");    
    leg_logfits_NPynew->Draw();    

  }


  ///7 TeV fit: two variations using Cross sections, each w/ 3 free parameters min. note; sqrt(s) should be in GeV
  ///cite: CERN-PH-EP/2011-053, 2018/10/22, CMS-QCD-10-011, 
  ///    : "Measurement of the Inclusive Jet Cross Section in pp Collisions at sqrt(s)=7 TeV"
  ///    : arixiv:1106.0208v1 [hep-ex] 1 Jun 2011
  ///    : Submitted to Physical Review Letters
  /// A) ddxsec[p_T] = D * [(p_T/GeV)^A] * [(1.-(2.*p_T/sqrt(s)]*cosh(ymin))^B] * [exp(-1*C/(p_T/GeV))]
  /// B) ddxsec[p_T] = [integral(ddxsec[p_T])] * [(p_T/GeV)^A] * [(1 - (2 * p_T / sqrt(s)]*cosh(ymin))^B] * [exp(-1*C/(p_T/GeV))]
  TF1 *_7TeVFitA_ynew  =NULL, *_7TeVFitB_ynew  =NULL;
  TF1 *_7TeVFitA_NPynew=NULL, *_7TeVFitB_NPynew=NULL;
  
  TCanvas *plot_7tevfits_ynew=NULL,*plot_7tevfits_NPynew=NULL;
  TLegend  *leg_7tevfits_ynew=NULL, *leg_7tevfits_NPynew=NULL;  
  if(useFitWeights && use7TeVFit){

    std::cout<<"making 7tev fit strings"<<std::endl;
    std::string _7tevfit_str="(TMath::Power(x,-1.*[0]))* (TMath::Power( (1.-(2.*x/(5.02*1000.) )),[1]))* (TMath::Exp(-1.*([2]/x)))";
    std::cout<<"_7tevfit_str="<<_7tevfit_str<<std::endl;
    
    std::string _7tevfitA_str="([3])* "+_7tevfit_str;    //this string is same for ynew and NPynew
    std::cout<<"_7tevfitA_str="<<_7tevfitA_str<<std::endl;
    
    // --- NLO --- //    
    /// check spectra with fits
    plot_7tevfits_ynew = new TCanvas("plot_7tevfits_ynew", "plot 7tevfits ynew",1200,800);
    plot_7tevfits_ynew->cd()->SetLogx(1);
    plot_7tevfits_ynew->cd()->SetLogy(1);
    plot_7tevfits_ynew->cd();
    
    float norm_ynew=theory_ynew_clone->Integral("width")*4.;
    std::cout<<"norm_ynew="<<norm_ynew<<std::endl;
    std::string norm_ynew_str=std::to_string(norm_ynew);
    std::cout<<"norm_ynew_str="<<norm_ynew_str<<std::endl;
    std::string _7tevfitB_ynew_str="("+norm_ynew_str+")* "+_7tevfit_str;
    std::cout<<"_7tevfitB_ynew_str="<<_7tevfitB_ynew_str<<std::endl;    
    
    _7TeVFitA_ynew =new TF1("_7tevfitA_ynew", (_7tevfitA_str     ).c_str() ,xlo,xhi);
    _7TeVFitB_ynew =new TF1("_7tevfitB_ynew" , (_7tevfitB_ynew_str).c_str(),xlo,xhi);
    make7TeVFits(theory_ynew_clone,  _7TeVFitA_ynew, _7TeVFitB_ynew, norm_ynew);
    
    //draw
    std::cout<<"drawing ynew fits and cross section on canvas"<<std::endl;
    theory_ynew_clone->SetTitle("7 TeV Fits for NLO #||{y}<2.0;Jet p_{T};Smear Weight");    
    theory_ynew_clone->DrawClone("HIST E");        
    _7TeVFitA_ynew->Draw("SAME");
    _7TeVFitB_ynew->Draw("SAME");
    
    leg_7tevfits_ynew=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");    
    leg_7tevfits_ynew->AddEntry(theory_ynew_clone , "Weighted NLO Jet Counts for #||{y} < 2.0" , "lp");//change desc TO DO
    leg_7tevfits_ynew->AddEntry(_7TeVFitB_ynew   , "4 Param. 7 TeV Fit for NLO" , "l");
    leg_7tevfits_ynew->AddEntry(_7TeVFitA_ynew   , "3 Param. 7 TeV Fit for NLO" , "l");    
    leg_7tevfits_ynew->Draw();    
    
    
    // --- NP+NLO --- //
    /// check spectra with fits
    plot_7tevfits_NPynew = new TCanvas("plot_7tevfits_NPynew", "plot 7tevfits NPynew",1200,800);
    plot_7tevfits_NPynew->cd()->SetLogx(1);
    plot_7tevfits_NPynew->cd()->SetLogy(1);
    plot_7tevfits_NPynew->cd();
    
    float norm_NPynew=theory_NPynew_clone->Integral("width")*4.;
    std::cout<<"norm_NPynew="<<norm_NPynew<<std::endl;
    std::string norm_NPynew_str=std::to_string(norm_NPynew);
    std::cout<<"norm_NPynew_str="<<norm_NPynew_str<<std::endl;
    std::string _7tevfitB_NPynew_str="("+norm_NPynew_str+")* "+_7tevfit_str;
    std::cout<<"_7tevfitB_NPynew_str="<<_7tevfitB_NPynew_str<<std::endl;        
    
    std::cout<<"fitting hist theory_NPynew_clone"<<std::endl;
    _7TeVFitA_NPynew =new TF1("_7tevfitA_NPynew",(_7tevfitA_str       ).c_str(),xlo,xhi);
    _7TeVFitB_NPynew=new TF1( "_7tevfitB_NPynew" ,(_7tevfitB_NPynew_str).c_str(),xlo,xhi);
    make7TeVFits(theory_NPynew_clone, _7TeVFitA_NPynew, _7TeVFitB_NPynew, norm_NPynew);
    
    //draw


    std::cout<<"drawing NPynew fits and cross section on canvas"<<std::endl;
    theory_NPynew_clone->SetTitle("7 TeV Fits for NP+NLO #||{y}<2.0;Jet p_{T};Smear Weight");    
    theory_NPynew_clone->DrawClone("HIST E");    
    
    _7TeVFitA_NPynew->Draw("SAME");
    _7TeVFitB_NPynew->Draw("SAME");
    
    leg_7tevfits_NPynew=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
    
    leg_7tevfits_NPynew->AddEntry(theory_NPynew_clone , "Weighted NP+NLO Jet Counts for #||{y} < 2.0" , "lp");//change desc TO DO
    leg_7tevfits_NPynew->AddEntry(_7TeVFitA_NPynew   ,  "3 Param. 7 TeV Fit for NP+NLO" , "l");
    leg_7tevfits_NPynew->AddEntry(_7TeVFitB_NPynew   ,  "4 Param. 7 TeV Fit for NP+NLO" , "l");    
    leg_7tevfits_NPynew->Draw();    
    //assert(false);
    
  }
  
  
  //-----------------------------------------------------------------------------------------//
  //////////////////////  START production of Smeared NLO spectra  ////////////////////////////
  //-----------------------------------------------------------------------------------------//
  std::cout<<"creating TH1 for toy NLO spectra generation, RooUnfoldResponse class, etc."<<std::endl<<std::endl;
  TH1D *theory_rnd_ynew       = new TH1D("theory_rnd_ynew","theory_rnd_ynew",             n_thybins_incl    , thyBins_incl);   
  TH1D *smeared_rnd_ynew      = new TH1D("smeared_rnd_ynew","smeared_rnd_ynew",           n_smearedbins_incl, smearedBins_incl);   
  TH1D *smeared_rnd_ynew_test = new TH1D("smeared_rnd_ynew_test","smeared_rnd_ynew_test", n_smearedbins_incl, smearedBins_incl);   
  
  TH2D* response_ynew_th2=new TH2D("response_ynew_th2","response_ynew_th2",
				   (Int_t)n_smearedbins_incl, (Double_t*)smearedBins_incl,
				   (Int_t)n_thybins_incl, (Double_t*)thyBins_incl);  

  //--------------------------------------------------------------------------------------------------------//
  if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
    makeToySpectra( (TH1D*)theory_ynew_clone, 
		    (TSpline3*)spline3_ynew, 
		    (TF1*)fJER_ynew, 
		    nEvents, 
		    theory_rnd_ynew, 
		    smeared_rnd_ynew, 
		    smeared_rnd_ynew_test, 
		    (TH2D*)response_ynew_th2,
		    (TF1*) spline3_ynew_ext); }
  else if(useFitWeights&&useModLogFit){//looks like modlog4Fit wins!
    makeToySpectra( (TH1D*)theory_ynew_clone, 
		    (TF1*)modLog4Fit_ynew, // (TF1*)modLog3Fit_ynew, // 
		    (TF1*)fJER_ynew, 
		    nEvents, 
		    theory_rnd_ynew, 
		    smeared_rnd_ynew, 
		    smeared_rnd_ynew_test, 
		    (TH2D*)response_ynew_th2); }  
  else if(useFitWeights&&use7TeVFit){//must decide between A and B TO DO
    makeToySpectra( (TH1D*)theory_ynew_clone, 
		    (TF1*)_7TeVFitA_ynew,// (TF1*)_7TeVFitB_ynew, //
		    (TF1*)fJER_ynew, 
		    nEvents, 
		    theory_rnd_ynew, 
		    smeared_rnd_ynew, 
		    smeared_rnd_ynew_test, 
		    (TH2D*)response_ynew_th2); }  

  //--------------------------------------------------------------------------------------------------------//

  divideBinWidth(theory_rnd_ynew);  //for normalization only.
  theory_rnd_ynew->Scale(1./4.);//etabinwidth
  
  divideBinWidth(smeared_rnd_ynew); 
  smeared_rnd_ynew->Scale(1./4.);

  divideBinWidth(smeared_rnd_ynew_test); 
  smeared_rnd_ynew_test->Scale(1./4.);
  
  double norm_rnd_ynew=theory_ynew_clone->Integral()/theory_rnd_ynew->Integral();
  response_ynew_th2->Scale(   norm_rnd_ynew);
  
  theory_rnd_ynew->Scale( norm_rnd_ynew);
  multiplyBinWidth(theory_rnd_ynew);  //for creating resp matrix
  theory_rnd_ynew->Scale(4.);//etabinwidth
  
  smeared_rnd_ynew->Scale(norm_rnd_ynew);  
  multiplyBinWidth(smeared_rnd_ynew); 
  smeared_rnd_ynew->Scale(4.);  

  smeared_rnd_ynew_test->Scale(norm_rnd_ynew);  
  multiplyBinWidth(smeared_rnd_ynew_test); 
  smeared_rnd_ynew_test->Scale(4.);
  
  divideBinWidth(theory_rnd_ynew);  //for normalization only.
  theory_rnd_ynew->Scale(1./4.);//etabinwidth
  
  divideBinWidth(smeared_rnd_ynew); 
  smeared_rnd_ynew->Scale(1./4.);
  
  divideBinWidth(smeared_rnd_ynew_test); 
  smeared_rnd_ynew_test->Scale(1./4.);
  
  //plotting smeared+theory toy mc xsection
  theory_rnd_ynew->SetLineColor(kRed); 
  smeared_rnd_ynew->SetLineColor(kBlue);    
  smeared_rnd_ynew_test->SetLineColor(kBlue-5);    
  
  Float_t maxy=theory_rnd_ynew->GetMaximum();    std::cout<<"maxy="<<maxy<<std::endl;
  maxy*=2.;  std::cout<<"maxy="<<maxy<<std::endl;
  Float_t miny=smeared_rnd_ynew->GetMinimum();   std::cout<<"miny="<<miny<<std::endl;
  miny/=2.; std::cout<<"miny="<<miny<<std::endl;
  theory_rnd_ynew->SetAxisRange(miny,maxy,"Y");  
  theory_rnd_ynew->SetTitle( ((std::string)("Toy MC Incl. NLO #sigma;;"+ddxsec_yax)).c_str() );
  
  

  ///// check plots True/Smeared 
  TCanvas *plot_ynew_true_smeared_rat = new TCanvas("plot_ynew_true_smeared_rat", "plot ynew true smeared rat",1400,1200);  
  plot_ynew_true_smeared_rat->cd();

  
  TPad* overlay_pad=new TPad("overlay_pad","pad for spectra",0.0,0.3,1.0,1.0);
  overlay_pad->SetLogx(1);
  overlay_pad->SetLogy(1);
  overlay_pad->SetBottomMargin(0);
  overlay_pad->Draw();//tells the canvas theres a pad at the coordinates specified in the new constructor
  overlay_pad->cd();
  
  theory_rnd_ynew->DrawClone("HIST E");  
  smeared_rnd_ynew->DrawClone("HIST E SAME");
  smeared_rnd_ynew_test->DrawClone("HIST E SAME");
  
  TLegend* leg_smear=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
  leg_smear->AddEntry(theory_rnd_ynew  , "Truth Incl. NLO #sigma" , "lp" );
  leg_smear->AddEntry(smeared_rnd_ynew , "Smeared Incl. NLO #sigma" , "lp" );
  leg_smear->AddEntry(smeared_rnd_ynew_test , "Opp Side Smeared Incl. NLO #sigma" , "lp" );
  leg_smear->Draw();
  
  plot_ynew_true_smeared_rat->cd();//go back to main canvas before doing new pad
  
  TPad* ratio_pad=new TPad("ratio_pad","pad for ratio",0.0,0.05,1.0,0.30);
  ratio_pad->SetLogx(1);
  ratio_pad->SetLogy(0);
  ratio_pad->SetTopMargin(0);
  //ratio_pad->SetBottomMargin(0.3);
  ratio_pad->Draw();//tell canvas there's a pad at the specified coordinates in the new constructor
  ratio_pad->cd();
  
  std::cout<<"creating ratio hist of true/smeared"<<std::endl;
  TH1D *ynew_true_smeared_rat=(TH1D*)theory_rnd_ynew->Clone("ynew_true_smeared_rat");   
  if(printBaseDebug)ynew_true_smeared_rat->Print("base");
  ynew_true_smeared_rat->SetTitle("; Jet p_T (GeV); True/Smeared");
  ynew_true_smeared_rat->Divide(ynew_true_smeared_rat,smeared_rnd_ynew,1.,1.,"B");   
  ynew_true_smeared_rat->GetXaxis()->SetMoreLogLabels(true);
  ynew_true_smeared_rat->GetXaxis()->SetNoExponent(true);
  ynew_true_smeared_rat->SetAxisRange(0.6,1.4,"Y");
  std::cout<<"creating ratio hist of true/opp side smeared"<<std::endl;

  TH1D *ynew_true_smeared_test_rat=(TH1D*)theory_rnd_ynew->Clone("ynew_true_smeared_test_rat");   
  if(printBaseDebug)ynew_true_smeared_test_rat->Print("base");  
  ynew_true_smeared_test_rat->SetTitle("; Jet p_T (GeV); True/Smeared");
  ynew_true_smeared_test_rat->Divide(ynew_true_smeared_test_rat,smeared_rnd_ynew_test,1.,1.,"B");   
  ynew_true_smeared_test_rat->GetXaxis()->SetMoreLogLabels(true);
  ynew_true_smeared_test_rat->GetXaxis()->SetNoExponent(true);
  ynew_true_smeared_test_rat->SetAxisRange(0.6,1.4,"Y");
  ynew_true_smeared_test_rat->SetLineColor(kRed-2);

  ynew_true_smeared_rat->DrawClone("HIST E");        
  ynew_true_smeared_test_rat->DrawClone("HIST E SAME");        
  
  lineatone->SetX2(thyBins_incl[n_thybins_incl]);
  lineatone->DrawClone();

  TCanvas* plot_response_ynew_th2=new TCanvas("plot_response_ynew_th2","plot response ynew th2",1200, 1000);
  plot_response_ynew_th2->cd()->SetLogx(1);
  plot_response_ynew_th2->cd()->SetLogy(1);
  plot_response_ynew_th2->cd()->SetLogz(1);
  plot_response_ynew_th2->cd();
  
  response_ynew_th2->SetTitle("Toy MC Response Matrix from JER Smeared NLO #sigma_{jet}");

  response_ynew_th2->GetXaxis()->SetTitle("RECO Jet p_{T} [GeV]");
  response_ynew_th2->GetXaxis()->SetNoExponent(true);
  response_ynew_th2->GetXaxis()->SetMoreLogLabels(true);
  response_ynew_th2->GetYaxis()->SetTitle("GEN Jet p_{T} [GeV]");
  response_ynew_th2->GetYaxis()->SetNoExponent(true);
  response_ynew_th2->GetYaxis()->SetMoreLogLabels(true);
  
  setTH2_ZAxisRange(response_ynew_th2);  
  response_ynew_th2->Draw("COLZ");
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
  std::cout<<"creating TH1 for toy NP+NLO spectra generation"<<std::endl<<std::endl;
  TH1D *theory_rnd_NPynew    = new TH1D("theory_rnd_NPynew","theory_rnd_NPynew", n_thybins_incl, thyBins_incl);   
  TH1D *smeared_rnd_NPynew = new TH1D("smeared_rnd_NPynew","smeared_rnd_NPynew", n_smearedbins_incl, smearedBins_incl);     
  TH1D *smeared_rnd_NPynew_test = new TH1D("smeared_rnd_NPynew_test","smeared_rnd_NPynew_test", n_smearedbins_incl, smearedBins_incl);     

  TH2D* response_NPynew_th2=new TH2D("response_NPynew_th2","response_NPynew_th2",
				     (Int_t)n_smearedbins_incl, (Double_t*)smearedBins_incl,
				     (Int_t)n_thybins_incl, (Double_t*)thyBins_incl);
  
  //----------------------------------------------------------------------------------------------------------------//
  if(useSplineWeights){
    makeToySpectra(theory_NPynew_clone, 
		   (TSpline3*)spline3_NPynew, 
		   (TF1*)fJER_ynew, 
		   nEvents, 
		   theory_rnd_NPynew, 
		   smeared_rnd_NPynew, 
		   smeared_rnd_NPynew_test, 
		   (TH2D*) response_NPynew_th2, 
		   (TF1*) spline3_NPynew_ext); }		    //NULL ) ;}//
  else if(useFitWeights&&useModLogFit){
    makeToySpectra(theory_NPynew_clone, 
		   (TF1*)modLog4Fit_NPynew,
		   (TF1*)fJER_ynew, 
		   nEvents, 
		   theory_rnd_NPynew, 
		   smeared_rnd_NPynew, 
		   smeared_rnd_NPynew_test, 
		   (TH2D*) response_NPynew_th2);}
  else if(useFitWeights&&use7TeVFit){//must decide between A and B TO DO
    makeToySpectra( (TH1D*)theory_NPynew_clone, 
		    (TF1*)_7TeVFitA_NPynew,// (TF1*)_7TeVFitB_NPynew, //
		    (TF1*)fJER_ynew, 
		    nEvents, 
		    theory_rnd_NPynew, 
		    smeared_rnd_NPynew, 
		    smeared_rnd_NPynew_test, 
		    (TH2D*)response_NPynew_th2); }  

  //----------------------------------------------------------------------------------------------------------------//
  
  divideBinWidth(theory_rnd_NPynew);  //for normalization only.
  theory_rnd_NPynew->Scale(1./4.);//etabinwidth
  
  divideBinWidth(smeared_rnd_NPynew); 
  smeared_rnd_NPynew->Scale(1./4.);

  divideBinWidth(smeared_rnd_NPynew_test); 
  smeared_rnd_NPynew_test->Scale(1./4.);

  double norm_rnd_NPynew=theory_NPynew_clone->Integral()/theory_rnd_NPynew->Integral();
  response_NPynew_th2->Scale(   norm_rnd_NPynew); 

  theory_rnd_NPynew->Scale( norm_rnd_NPynew);
  smeared_rnd_NPynew->Scale(norm_rnd_NPynew);
  smeared_rnd_NPynew_test->Scale(norm_rnd_NPynew);
  
  multiplyBinWidth(theory_rnd_NPynew);  //for creating resp matrix
  theory_rnd_NPynew->Scale(4.);//etabinwidth

  multiplyBinWidth(smeared_rnd_NPynew); 
  smeared_rnd_NPynew->Scale(4.);

  multiplyBinWidth(smeared_rnd_NPynew_test); 
  smeared_rnd_NPynew_test->Scale(4.);
  
  divideBinWidth(theory_rnd_NPynew);  //for normalization only.
  theory_rnd_NPynew->Scale(1./4.);//etabinwidth

  divideBinWidth(smeared_rnd_NPynew); 
  smeared_rnd_NPynew->Scale(1./4.);

  divideBinWidth(smeared_rnd_NPynew_test); 
  smeared_rnd_NPynew_test->Scale(1./4.);

  theory_rnd_NPynew->SetLineColor(kRed);
  smeared_rnd_NPynew->SetLineColor(kBlue);    
  smeared_rnd_NPynew_test->SetLineColor(kBlue-5);    
  
  if(theory_rnd_NPynew->GetBinContent(theory_rnd_NPynew->GetNbinsX() )<0.)
    theory_rnd_NPynew->SetBinContent(theory_rnd_NPynew->GetNbinsX(),1.e-10);
  if(smeared_rnd_NPynew->GetBinContent(smeared_rnd_NPynew->GetNbinsX() )<0.)
    smeared_rnd_NPynew->SetBinContent(smeared_rnd_NPynew->GetNbinsX(),1.e-08);
  
  
  Float_t NPmaxy=theory_rnd_NPynew->GetMaximum();    std::cout<<"NPmaxy="<<NPmaxy<<std::endl;
  NPmaxy*=2.;  std::cout<<"NPmaxy="<<NPmaxy<<std::endl;
  Float_t NPminy=smeared_rnd_NPynew->GetMinimum();   std::cout<<"NPminy="<<NPminy<<std::endl;
  NPminy/=2.; std::cout<<"NPminy="<<NPminy<<std::endl;
  theory_rnd_NPynew->SetAxisRange(NPminy,NPmaxy,"Y");
  theory_rnd_NPynew->SetTitle( ((std::string)("Toy MC Incl. NP+NLO #sigma;;"+ddxsec_yax)).c_str() );
  
  ///// check plots True/Smeared 
  TCanvas *plot_NPynew_true_smeared_rat = new TCanvas("plot_NPynew_true_smeared_rat", "plot NPynew true smeared rat",1400,1200);  
  plot_NPynew_true_smeared_rat->cd();
  
  TPad* NPoverlay_pad=new TPad("NPoverlay_pad","pad for spectra",0.0,0.3,1.0,1.0);
  NPoverlay_pad->SetLogx(1);
  NPoverlay_pad->SetLogy(1);
  NPoverlay_pad->SetBottomMargin(0);
  NPoverlay_pad->Draw();//tells the canvas theres a pad at the coordinates specified in the new constructor
  NPoverlay_pad->cd();
  
  theory_rnd_NPynew->DrawClone("HIST E");  
  smeared_rnd_NPynew->DrawClone("HIST E SAME");
  
  TLegend* NPleg_smear=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
  NPleg_smear->AddEntry(theory_rnd_NPynew  , "Truth Incl. NP+NLO #sigma" , "lp" );
  NPleg_smear->AddEntry(smeared_rnd_NPynew , "Smeared Incl. NP+NLO #sigma" , "lp" );
  NPleg_smear->Draw();
  
  plot_NPynew_true_smeared_rat->cd();//go back to main canvas before doing new pad
  
  TPad* NPratio_pad=new TPad("NPratio_pad","pad for ratio",0.0,0.05,1.0,0.30);
  NPratio_pad->SetLogx(1);
  NPratio_pad->SetLogy(0);
  NPratio_pad->SetTopMargin(0);
  //NPratio_pad->SetBottomMargin(0.3);
  NPratio_pad->Draw();//tell canvas there's a pad at the specified coordinates in the new constructor
  NPratio_pad->cd();
  
  std::cout<<"creating ratio hist of true/smeared for NP+NLO"<<std::endl;
  TH1D *NPynew_true_smeared_rat=(TH1D*)theory_rnd_NPynew->Clone("NPynew_true_smeared_rat"); 
  if(printBaseDebug)NPynew_true_smeared_rat->Print("base");
  NPynew_true_smeared_rat->SetTitle("; Jet p_T (GeV); True/Smeared");
  NPynew_true_smeared_rat->Divide(NPynew_true_smeared_rat,smeared_rnd_NPynew,1.,1.,"B");     
  NPynew_true_smeared_rat->GetXaxis()->SetMoreLogLabels(true);
  NPynew_true_smeared_rat->GetXaxis()->SetNoExponent(true);
  NPynew_true_smeared_rat->SetAxisRange(0.6,1.4,"Y");  
  
  std::cout<<"creating ratio hist of true/opp side smeared for NP+NLO"<<std::endl;
  TH1D *NPynew_true_smeared_test_rat=(TH1D*)theory_rnd_NPynew->Clone("NPynew_true_smeared_test_rat"); 
  if(printBaseDebug)NPynew_true_smeared_test_rat->Print("base");
  NPynew_true_smeared_test_rat->SetTitle("; Jet p_T (GeV); True/Smeared");
  NPynew_true_smeared_test_rat->Divide(NPynew_true_smeared_test_rat,smeared_rnd_NPynew_test,1.,1.,"B");     
  NPynew_true_smeared_test_rat->GetXaxis()->SetMoreLogLabels(true);
  NPynew_true_smeared_test_rat->GetXaxis()->SetNoExponent(true);
  NPynew_true_smeared_test_rat->SetAxisRange(0.6,1.4,"Y");
  NPynew_true_smeared_test_rat->SetLineColor(kRed-2);  

  NPynew_true_smeared_rat->DrawClone("HIST E");        
  NPynew_true_smeared_test_rat->DrawClone("HIST E SAME");        
  
  lineatone->DrawClone();
  
  
  TCanvas* plot_response_NPynew_th2=new TCanvas("plot_response_NPynew_th2","plot response NPynew th2",1200, 1000);
  plot_response_NPynew_th2->cd()->SetLogx(1);
  plot_response_NPynew_th2->cd()->SetLogy(1);
  plot_response_NPynew_th2->cd()->SetLogz(1);
  plot_response_NPynew_th2->cd();
  
  response_NPynew_th2->SetTitle("Toy MC Response Matrix from JER Smeared NP+NLO #sigma_{jet}");
  response_NPynew_th2->GetXaxis()->SetTitle("RECO Jet p_{T} [GeV]");
  response_NPynew_th2->GetXaxis()->SetNoExponent(true);
  response_NPynew_th2->GetXaxis()->SetMoreLogLabels(true);
  
  response_NPynew_th2->GetYaxis()->SetTitle("GEN Jet p_{T} [GeV]");
  response_NPynew_th2->GetYaxis()->SetNoExponent(true);
  response_NPynew_th2->GetYaxis()->SetMoreLogLabels(true);
  
  setTH2_ZAxisRange(response_NPynew_th2);
  
  response_NPynew_th2->Draw("COLZ");
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
  
  TH1D* theory_ynew_ratio_rnd_ynew=(TH1D*)theory_ynew_clone->Clone ( 
							      ( (std::string)theory_ynew->GetName() +"_debughist") .c_str()
							       );
  theory_ynew_ratio_rnd_ynew->Divide(theory_rnd_ynew);
  theory_ynew_ratio_rnd_ynew->SetLineColor(kBlack);
  theory_ynew_ratio_rnd_ynew->SetTitle(thy_rnd_rat_title.c_str());

  
  TH1D* theory_NPynew_ratio_rnd_NPynew=(TH1D*)theory_NPynew_clone->Clone ( 
									  ( (std::string)theory_NPynew->GetName() +"_debughist") .c_str()
									   );
  theory_NPynew_ratio_rnd_NPynew->Divide(theory_rnd_NPynew);
  theory_NPynew_ratio_rnd_NPynew->SetLineColor(kGray+1);
  theory_NPynew_ratio_rnd_NPynew->SetTitle(thy_rnd_rat_title.c_str());
  
  theory_ynew_ratio_rnd_ynew->SetAxisRange    (0.1,1.5,"Y");
  theory_NPynew_ratio_rnd_NPynew->SetAxisRange(0.1,1.5,"Y");
  

  std::cout<<"NLO THY/RND THY LAST BIN CONTENT = "<<theory_ynew_ratio_rnd_ynew->GetBinContent(theory_ynew_ratio_rnd_ynew->GetNbinsX()) <<std::endl;
  std::cout<<"NLO THY/RND THY LAST BIN CONTENT = "<<theory_NPynew_ratio_rnd_NPynew->GetBinContent(theory_NPynew_ratio_rnd_NPynew->GetNbinsX()) <<std::endl;
  
  
  theory_ynew_ratio_rnd_ynew->DrawClone("HIST E");
  theory_NPynew_ratio_rnd_NPynew->DrawClone("HIST E SAME");
  
  TLegend* leg_tru_rnd_rat=new TLegend(0.10, 0.70, 0.35, 0.9, NULL,"BRNDC");
  leg_tru_rnd_rat->AddEntry( theory_ynew_ratio_rnd_ynew, "NLO Spectra", "lp");
  leg_tru_rnd_rat->AddEntry( theory_NPynew_ratio_rnd_NPynew, "NLO+NP Spectra", "lp");
  leg_tru_rnd_rat->Draw();
  
  

  


  

  ////// Write out histos
  outf->cd();
  
  // NP FITS 
  fNP_y0->Write();
  fNP_y1->Write();
  fNP_y2->Write();
  fNP_y3->Write();
  fNP_ynew->Write();

  // NP FIT CANV
  plot_y0_to_y3_NPs->Write();  
  plot_totNPs->Write();
  plot_allNPs->Write();
  
  //JET ENERGY RESOLUTION
  fJER_ynew->Write();//FIT
  plotJER->Write();//CANV
  
  // NLO HISTS
  theory_y0->Write();
  theory_y1->Write();
  theory_y2->Write();
  theory_y3->Write();  
  theory_ynew->Write();//made from the 4 above so it goes in the "opened" section
  
  plot_y0_to_y3_NLOxsec->Write();
  plot_totNLOxsec->Write();
  plot_allNLOxsec->Write();
  
  // NP+NLO HISTS
  theory_NPynew->Write();  
  theory_NPy0->Write();  
  theory_NPy1->Write();  
  theory_NPy2->Write();  
  theory_NPy3->Write();  
  
  plot_y0_to_y3_NPNLOxsec->Write();
  plot_totNPNLOxsec->Write();
  plot_allNPNLOxsec->Write();

  // CLONES 4 SMEARING+TOYMC
  theory_ynew_clone->Write();
  theory_NPynew_clone->Write();
  

  if(useSplineWeights){  // SPLINES
    spline3_ynew->Write();
    spline3_NPynew->Write();
    if((bool)spline3_ynew_ext)spline3_ynew_ext->Write();
    if((bool)spline3_NPynew_ext)spline3_NPynew_ext->Write();
    plot_splines->Write();  
  }
  else if(useFitWeights && useModLogFit){//FITS
    logFit_ynew->Write();
    modLogFit_ynew->Write();
    modLog2Fit_ynew->Write();
    modLog3Fit_ynew->Write();
    modLog4Fit_ynew->Write();
    modLog5Fit_ynew->Write();
    plot_logfits_ynew->Write();  
    
    logFit_NPynew->Write();
    modLogFit_NPynew->Write();
    modLog2Fit_NPynew->Write();
    modLog3Fit_NPynew->Write();
    modLog4Fit_NPynew->Write();
    modLog5Fit_NPynew->Write();
    plot_logfits_NPynew->Write();      
  }
  else if(useFitWeights && use7TeVFit){//FITS
    _7TeVFitA_ynew->Write();
    _7TeVFitB_ynew->Write();
    plot_7tevfits_ynew->Write();  
    
    _7TeVFitA_NPynew->Write();
    _7TeVFitB_NPynew->Write();
    plot_7tevfits_NPynew->Write();      
  }
  
  
  //TOY MC, JER SMEARED NLO SPECTRA
  theory_rnd_ynew->Write();
  theory_ynew_ratio_rnd_ynew->Write();
  
  smeared_rnd_ynew->Write(); 
  smeared_rnd_ynew_test->Write(); 

  ynew_true_smeared_rat->Write();  //ratio
  ynew_true_smeared_test_rat->Write();  //ratio
  plot_ynew_true_smeared_rat->Write();
  
  response_ynew_th2->Write();
  plot_response_ynew_th2->Write();  
    


  ////TOY MC NP+NLO SPECTRA, SMEARED+TRUE   
  theory_rnd_NPynew->Write();
  theory_NPynew_ratio_rnd_NPynew->Write();
  
  smeared_rnd_NPynew->Write(); 
  smeared_rnd_NPynew_test->Write(); 

  NPynew_true_smeared_rat->Write();  //ratio
  NPynew_true_smeared_test_rat->Write();  //ratio
  plot_NPynew_true_smeared_rat->Write();
  
  response_NPynew_th2->Write();
  plot_response_NPynew_th2->Write();    


  // toy spectra comparison w/ orig calculation
  plot_true_rnd_rat->Write();
  

  outf->Write();
  return 1;
  //    assert(false);    
}


//  steering ---------------------------------------------------------------------------------
int main(int argc, char* argv[]){  
  int rStatus = -1;   
  
  if(argc==2) {
    rStatus=smearTheorySpectra_gaussCoreJER( (std::string) argv[1] );      
  }
  else if(argc==3) {
    rStatus=smearTheorySpectra_gaussCoreJER( (std::string) argv[1] ,  ( (bool)std::atoi( argv[2]) ) )  ;      
  }
  else if(argc==4){
    rStatus=smearTheorySpectra_gaussCoreJER( (std::string) argv[1] ,  ( (bool)std::atoi( argv[2]) ) , (std::string) argv[3])  ;      
  }
  
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}
  
  
















