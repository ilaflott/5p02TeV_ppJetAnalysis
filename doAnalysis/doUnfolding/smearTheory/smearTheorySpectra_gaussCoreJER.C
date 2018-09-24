#include "smearTheorySpectra.h"

const bool printBaseDebug=true;
const int nEvents=1e+09;  /// Increase your stats here, typically 100M events are enough 
//const int nEvents=1e+07;  /// debug nevents
//const int nEvents=1e+06;  /// debug nevents
//const int nEvents=1e+05;  /// debug nevents
std::string ddxsec_yax="#frac{d^{2}#sigma}{dp_{T}dy} [nb/GeV]";

const std::string NPCorrFile="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory/NPCorr5TeV/NLOpNP_InclusiveJets5TeV.root";
const std::string JERCorrFile="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots_JERS/output/ak4PF_MCJEC_00eta20_08.07.18.root";


int debugInt=0;
void debugcout(){
  std::cout<<"location #"<<debugInt<<endl;
  debugInt++;
  return; }

void smearTheorySpectra_gaussCoreJER( string inputString ){
  
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  
  // Output File to write to
  string outputFile= inputString+"_gaussSmear_00eta20.root";
  cout<<"opening output:"<<outputFile<<endl<<endl;
  TFile *outf    = new TFile(outputFile.c_str(), "RECREATE" );
  
  // Input, Theory spectrum
  string inputFile = "fNLOJetsSpectra/R04/"+inputString+".root";
  cout<<"opening input file:"<<inputFile<<endl<<endl;
  TFile* fin_NLO=TFile::Open(inputFile.c_str());
  

  /////////////// NPs to NLO jet spectra ////////////////////////////////////////////////////////  
  cout<<"opening NP corr file + fits! " << endl;
  TFile* fin_NP=TFile::Open(NPCorrFile.c_str());  
  std::string NPCorrFits_str;
  if(NPCorrFile.find("nnlo")!=std::string::npos)
    NPCorrFits_str="fNPC_POWPY8_R4_etabin";//raghavs suggestion to use POWPY8 for NNLO thy.
  else
    NPCorrFits_str="fNPC_HerwigEE4C_R4_etabin";
  /////////////// Get NP fit functions
  //TF1 *fNP_y0 = new TF1("fNP_y0",NP_y0_str.c_str(), thyBins_incl[0], thyBins_incl[n_thybins_incl]);
  //TF1 *fNP_y1 = new TF1("fNP_y1",NP_y1_str.c_str(), thyBins_incl[0], thyBins_incl[n_thybins_incl]);
  //TF1 *fNP_y2 = new TF1("fNP_y2",NP_y2_str.c_str(), thyBins_incl[0], thyBins_incl[n_thybins_incl]);
  //TF1 *fNP_y3 = new TF1("fNP_y3",NP_y3_str.c_str(), thyBins_incl[0], thyBins_incl[n_thybins_incl]);
  //TF1 *fNP_ynew = new TF1("fNP_ynew",NP_ynew_str.c_str(), thyBins_incl[0], thyBins_incl[n_thybins_incl]);  //avg of other four
  //TF1 *fNP_y0 = (TF1*)fin_NP->Get("fNPC_POWPY8_R4_etabin0");
  //TF1 *fNP_y1 = (TF1*)fin_NP->Get("fNPC_POWPY8_R4_etabin1");
  //TF1 *fNP_y2 = (TF1*)fin_NP->Get("fNPC_POWPY8_R4_etabin2");
  //TF1 *fNP_y3 = (TF1*)fin_NP->Get("fNPC_POWPY8_R4_etabin3");
  TF1 *fNP_y0 = (TF1*)fin_NP->Get((NPCorrFits_str+"0").c_str());
  TF1 *fNP_y1 = (TF1*)fin_NP->Get((NPCorrFits_str+"1").c_str());
  TF1 *fNP_y2 = (TF1*)fin_NP->Get((NPCorrFits_str+"2").c_str());
  TF1 *fNP_y3 = (TF1*)fin_NP->Get((NPCorrFits_str+"3").c_str());
  if(!fin_NP ||
     !fNP_y0 || !fNP_y1 || !fNP_y2 || !fNP_y3 )
    {cout<<"error, NP file and/or fits not found."<<endl; assert(false);}
  std::string fNP_ynew_str= "(" + 
    ((std::string)fNP_y0->GetName()) + "+" + 
    ((std::string)fNP_y1->GetName()) + "+" +
    ((std::string)fNP_y2->GetName()) + "+" +
    ((std::string)fNP_y3->GetName()) +")/4.";  
  TF1 *fNP_ynew=new TF1("fNP_ynew",fNP_ynew_str.c_str() ,thyBins_incl[0],thyBins_incl[n_thybins_incl]);
  //TF1 *fNP_ynew = new TF1("fNP_ynew",NP_ynew_str.c_str(), thyBins_incl[0], thyBins_incl[n_thybins_incl]);  //avg of other four
  //-----------------------------------


  cout<<"opening JER file + fits! " << endl;
  TFile* fin_JER=TFile::Open(JERCorrFile.c_str());
  /////////////// Create/Get JER fit function(s)
  TF1 *fJER_default = new TF1("fJER_ynew",gJER_ynew_str.c_str(),thyBins_incl[0], thyBins_incl[n_thybins_incl]);  //for comparison only
  fJER_default->SetLineColor(kBlue);
  fJER_default->SetMinimum(0.); //do i need to set min/max really?
  fJER_default->SetMaximum(0.25);
  TF1 *fJER_ynew = (TF1*)fin_JER->Get("hSigmaFit");
  if(!fin_JER || !fJER_ynew)
    {cout<<"error, JER file and/or fits not found."<<endl; assert(false);}
  fJER_ynew->SetLineColor(kRed);
  fJER_ynew->SetMinimum(0.); //do i need to set min/max really?
  fJER_ynew->SetMaximum(0.25);
  
  
  //TFile *f1 = new TFile(inputFile.c_str());  
  cout<<endl;
  
  



  /////////////// plots to to check NPs 
  // 2x2 canv of |y| bins
  TLine* lineatone=new TLine(thyBins_incl[0],1.,1000.,1.); lineatone->SetLineStyle(7);
  TCanvas *plot_y0_to_y3_NPs = new TCanvas("plot_y0_to_y3_NPs", "plot_y0_to_y3_NPs",1200,1000);
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
  TCanvas *plot_totNPs = new TCanvas("plot_totNPs", "plot_totNPs",600,500);//1600,1200);
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
  TCanvas *plot_allNPs = new TCanvas("plot_allNPs", "plot_allNPs",600,500);
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
  hNP_y0->SetTitle("NP Corr Fit for 0.0 < #||{y} < 0.5;Jet p_{T};NP Corr. Factor");//undo the titling because this is for this plot only
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
  
  
  
  
  /////////////// Core p_T Resolution ////////////////////////////////////////////////////////
  cout<<"plotting JER!"<<  endl;  //cout<<""<<  << endl;
  
  //-----------------------------------  
  
  /////////////// Create plot to check
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
  TCanvas *plotJER = new TCanvas("plotJER", "plotJER",900,600);
  plotJER->cd()->SetLogx(1);
  plotJER->cd();  
  //hJER_default->DrawClone("HIST E");    
  hJER_ynew->DrawClone("HIST E");    
  //REPLACE ME W/TEXT SAYING SAME THING
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

  
  
  //////////// NLO Thy calculation xsecions ////////////////////////////////////////////////////////
  cout<<"plotting NLO xsecs!"<<  endl;  //cout<<""<<  << endl;
  
  /////////////// get NLO xsec hists from file
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
  cout<<endl;
  //-----------------------------------
  cout <<"making 00eta20 hist from theory curves..."<<endl<<endl;
  TH1D *theory_ynew =    make00eta20Hist( (TH1D*) theory_y0, 
					  (TH1D*) theory_y1, 
					  (TH1D*) theory_y2, 
					  (TH1D*) theory_y3 );  
  theory_ynew->SetTitle(("NLO Inclusive #sigma_{jet}, #||{y} < 2.0;Jet p_{T};"+ddxsec_yax).c_str());
  theory_ynew->GetXaxis()->SetNoExponent(true);
  theory_ynew->GetXaxis()->SetMoreLogLabels(true);
  theory_ynew->SetLineColor(kTeal);
  cout<<"done making hist."<<endl;
  if(printBaseDebug)theory_ynew->Print("base");  
  //-----------------------------------
  

  /////////////// draw xsec hists from file

  //2X2 CANV
  TCanvas *plot_y0_to_y3_NLOxsec = new TCanvas("plot_y0_to_y3_NLOxsec", "plot_y0_to_y3_NLOxsec",1200,1000);
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
  TCanvas *plot_totNLOxsec = new TCanvas("plot_totNLOxsec", "plot_totNLOxsec",600,500);//1600,1200);
  plot_totNLOxsec->SetLogx(1);
  plot_totNLOxsec->SetLogy(1);
  plot_totNLOxsec->cd();
  theory_ynew->DrawClone("HIST E");
  //--------------------------------------------
  
  // plot of all xsecs on one hist
  TCanvas *plot_allNLOxsec= new TCanvas("plot_allNLOxsec", "plot_allNLOxsec",900,750);
  plot_allNLOxsec->SetLogx(1);
  plot_allNLOxsec->SetLogy(1);
  plot_allNLOxsec->cd();
  //--------------------------------------------
  theory_ynew->Scale(1.e+04);  
  theory_y0->Scale(1.e+03);  
  theory_y1->Scale(1.e+02);
  theory_y2->Scale(1.e+01);
  theory_y3->Scale(1.e+00);                                                                                      
  //--------------------------------------------  
  //theory_ynew->SetTitle("NLO Inclusive #sigma_{jet} all #||{y} bins");
  theory_y0->SetTitle("NLO Inclusive #sigma_{jet} all #||{y} bins");
  Double_t max=(Double_t)theory_y0->GetMaximum();
  Double_t min=(Double_t)theory_y3->GetMinimum();
  //theory_ynew->SetAxisRange((min/2.),(max*2.),"Y");
  theory_y0->SetAxisRange((min/2.),(max*2.),"Y");
  //theory_ynew->DrawClone("HIST E");
  theory_y0->DrawClone("HIST E");
  //theory_y0->DrawClone("HIST E SAME");
  theory_y1->DrawClone("HIST E SAME");
  theory_y2->DrawClone("HIST E SAME");
  theory_y3->DrawClone("HIST E SAME");  
  //-------------------------------------------- 
  TLegend* leg_allNLOxsec=new TLegend(0.60, 0.70, 0.9, 0.9, NULL,"BRNDC");
  //leg_allNPs->AddEntry(hNP_ynew , "#||{y} < 2.0 (sum of others)" , "l");
  //leg_allNLOxsec->AddEntry(theory_ynew   , "(0.0 < #||{y} < 2.0) x 10^{4}" , "lp");
  leg_allNLOxsec->AddEntry(theory_y0   , "(0.0 < #||{y} < 0.5) x 10^{3}" , "lp");
  leg_allNLOxsec->AddEntry(theory_y1   , "(0.5 < #||{y} < 1.0) x 10^{2}" , "lp");
  leg_allNLOxsec->AddEntry(theory_y2   , "(1.0 < #||{y} < 1.5) x 10^{1}" , "lp");
  leg_allNLOxsec->AddEntry(theory_y3   , "(1.5 < #||{y} < 2.0) x 10^{0}" , "lp");
  leg_allNLOxsec->Draw();
  //--------------------------------------------  
  //unscale the NLO xsec histograms for later use...
  theory_ynew->Scale(1.e-04);  
  theory_y0->Scale(  1.e-03);  
  theory_y1->Scale(  1.e-02);
  theory_y2->Scale(  1.e-01);
  theory_y3->Scale(  1.e-00);                                                                                      
  theory_y0->SetTitle(("NLO #sigma_{jet}, 0.0 < #||{y} < 0.5;Jet p_{T};"+ddxsec_yax).c_str());
  //--------------------------------------------
  
  
  
  
  
  //////////// Apply NPs to theory ////////////////////////////////////////////////////////
  cout<<"plotting NP+NLO xsecs!"<<  endl;  //cout<<""<<  << endl;

  /// Cross section ynew
  TH1D *theory_NPy0=(TH1D*)applyNPtoxsec(theory_y0, fNP_y0);
  TH1D *theory_NPy1=(TH1D*)applyNPtoxsec(theory_y1, fNP_y1);
  TH1D *theory_NPy2=(TH1D*)applyNPtoxsec(theory_y2, fNP_y2);
  TH1D *theory_NPy3=(TH1D*)applyNPtoxsec(theory_y3, fNP_y3);
  
  //TH1D *theory_NPynew=(TH1D*)applyNPtoxsec(theory_ynew, fNP_ynew);
  TH1D *theory_NPynew=    make00eta20Hist( (TH1D*) theory_NPy0, 
					   (TH1D*) theory_NPy1, 
					   (TH1D*) theory_NPy2, 
					   (TH1D*) theory_NPy3 );  
  theory_NPynew->SetTitle(("NP+NLO Inclusive #sigma_{jet}, #||{y} < 2.0;Jet p_{T};"+ddxsec_yax).c_str());
  theory_NPynew->GetXaxis()->SetNoExponent(true);
  theory_NPynew->GetXaxis()->SetMoreLogLabels(true);
  theory_NPynew->SetLineColor(kTeal);
  

  //2X2 CANV THEORY NPNLO xsec
  TCanvas *plot_y0_to_y3_NPNLOxsec = new TCanvas("plot_y0_to_y3_NPNLOxsec", "plot_y0_to_y3_NPNLOxsec",1200,1000);
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
  TCanvas *plot_totNPNLOxsec = new TCanvas("plot_totNPNLOxsec", "plot_totNPNLOxsec",600,500);//1600,1200);
  plot_totNPNLOxsec->SetLogx(1);
  plot_totNPNLOxsec->SetLogy(1);
  plot_totNPNLOxsec->cd();
  theory_NPynew->SetTitle("NP+NLO Inclusive #sigma_{jet} #||{y} < 2.0");
  theory_NPynew->DrawClone("HIST E");
  //--------------------------------------------
  

  // plot of all xsecs on one hist
  TCanvas *plot_allNPNLOxsec= new TCanvas("plot_allNPNLOxsec", "plot_allNPNLOxsec",900,750);
  plot_allNPNLOxsec->SetLogx(1);
  plot_allNPNLOxsec->SetLogy(1);
  plot_allNPNLOxsec->cd();
  //--------------------------------------------
  theory_NPynew->Scale(1.e+04);  
  theory_NPy0->Scale(1.e+03);  
  theory_NPy1->Scale(1.e+02);
  theory_NPy2->Scale(1.e+01);
  theory_NPy3->Scale(1.e+00);                                                                                      
  //--------------------------------------------  
  //theory_NPynew->SetTitle("NP+NLO Inclusive #sigma_{jet} all #||{y} bins");
  theory_NPy0->SetTitle("NP+NLO Inclusive #sigma_{jet}, all #||{y} bins");
  Double_t NPmax=(Double_t)theory_NPy0->GetMaximum();
  //Double_t NPmax=(Double_t)theory_NPynew->GetMaximum();
  Double_t NPmin=(Double_t)theory_NPy3->GetMinimum();
  theory_NPy0->SetAxisRange((NPmin/2.),(NPmax*2.),"Y");
  //theory_NPynew->DrawClone("HIST E");
  //theory_NPy0->DrawClone("HIST E SAME");
  theory_NPy0->DrawClone("HIST E");
  theory_NPy1->DrawClone("HIST E SAME");
  theory_NPy2->DrawClone("HIST E SAME");
  theory_NPy3->DrawClone("HIST E SAME");  
  //-------------------------------------------- 
  TLegend* leg_allNPNLOxsec=new TLegend(0.60, 0.70, 0.9, 0.9, NULL,"BRNDC");
  //leg_allNPs->AddEntry(hNP_ynew , "#||{y} < 2.0 (sum of others)" , "l");
  //leg_allNPNLOxsec->AddEntry(theory_NPynew   , "(0.0 < #||{y} < 2.0) x 10^{4}" , "lp");
  leg_allNPNLOxsec->AddEntry(theory_NPy0   , "(0.0 < #||{y} < 0.5) x 10^{3}" , "lp");
  leg_allNPNLOxsec->AddEntry(theory_NPy1   , "(0.5 < #||{y} < 1.0) x 10^{2}" , "lp");
  leg_allNPNLOxsec->AddEntry(theory_NPy2   , "(1.0 < #||{y} < 1.5) x 10^{1}" , "lp");
  leg_allNPNLOxsec->AddEntry(theory_NPy3   , "(1.5 < #||{y} < 2.0) x 10^{0}" , "lp");
  leg_allNPNLOxsec->Draw();
  //--------------------------------------------  
  //unscale the NLO xsec histograms for later use...
  theory_NPynew->Scale(1.e-04);  
  theory_NPy0->Scale(  1.e-03);  
  theory_NPy1->Scale(  1.e-02);
  theory_NPy2->Scale(  1.e-01);
  theory_NPy3->Scale(  1.e-00);                            
  theory_NPy0->SetTitle(("NP+NLO #sigma_{jet}, 0.0 < #||{y} < 0.5;Jet p_{T};"+ddxsec_yax).c_str());
  //--------------------------------------------





  

  /////////////// Fitting Spectra with spline3
  /////////////// Not worning well for y6 bin because Klaus's spectra ends < 500 GeV
  /////////////// We must fit spectra by hand for y6 bin
  
  
  ///Create Cubic Splines using Cross sections
  cout<<"creating TSpline of cross section hist"<<endl;
  TH1D* theory_ynew_spl3clone=(TH1D*)theory_ynew->Clone(  
							( ( (std::string) theory_ynew->GetName() ) +"_specialrebin").c_str() 
							  );
  theory_ynew_spl3clone=(TH1D*)theory_ynew_spl3clone->TH1::Rebin(n_thybins_incl, 
								 theory_ynew_spl3clone->GetName(), 
								 thyBins_incl );
  theory_ynew_spl3clone->SetLineColor(kCyan+4);    
  
  TSpline3 *spline3_ynew = new TSpline3( theory_ynew_spl3clone);
  spline3_ynew->SetName( ( (std::string)theory_ynew_spl3clone->GetName() + "_spline3").c_str() );
  spline3_ynew->SetLineColor(kAzure);  
  
  
  
  TH1D* theory_NPynew_spl3clone=(TH1D*)theory_NPynew->Clone(  
							    ( ( (std::string) theory_NPynew->GetName() ) +"_specialrebin").c_str() 
							      );  
  theory_NPynew_spl3clone=(TH1D*)theory_NPynew_spl3clone->TH1::Rebin(n_thybins_incl, 
								     theory_NPynew_spl3clone->GetName(), 
								     thyBins_incl );  
  theory_NPynew_spl3clone->SetLineColor(kCyan-6);
  
  TSpline3 *spline3_NPynew = new TSpline3( theory_NPynew_spl3clone);
  spline3_NPynew->SetName( ( (std::string)theory_NPynew_spl3clone->GetName() + "_spline3").c_str() );
  spline3_NPynew->SetLineColor(kAzure-8);  
  
  /// check spectra with splines
  cout<<"drawing Tspline and cross section on canvas"<<endl;
  TCanvas *plot_splines = new TCanvas("plot_splines", "plot_splines",1200,800);
  plot_splines->cd()->SetLogx(1);
  plot_splines->cd()->SetLogy(1);
  plot_splines->cd();
  
  //draw
  theory_ynew_spl3clone->SetTitle("Cubic Spline Fits #||{y}<2.0;Jet p_{T};Smear Weight");
  theory_ynew_spl3clone->DrawClone("HIST E");
  spline3_ynew->Draw("SAME");
  theory_NPynew_spl3clone->DrawClone("HIST E SAME");  
  spline3_NPynew->Draw("SAME");

  
  TLegend* leg_spline=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
  leg_spline->AddEntry(theory_ynew_spl3clone , "Weighted NLO Jet Counts for #||{y} < 2.0" , "lp");
  leg_spline->AddEntry(spline3_ynew   , "Cubic Spline Fit for NLO #||{y} < 2.0" , "l");
  leg_spline->AddEntry(theory_NPynew_spl3clone , "Weighted NP+NLO Jet Counts for #||{y} < 2.0" , "lp");
  leg_spline->AddEntry(spline3_NPynew   , "Cubic Spline Fit for NP+NLO #||{y} < 2.0" , "l");
  leg_spline->Draw();
  
  
  
  
  //////////////////////  START production of Smeared NLO spectra     //////////////////////  
  cout<<"creating TH1 for toy NLO spectra generation, RooUnfoldResponse class, etc."<<endl<<endl;
  //TH1D *theory_rnd_ynew    = new TH1D("theory_rnd_ynew","theory_rnd_ynew", n_thybins_incl, thyBins_incl);   
  //TH1D *smeared_rnd_ynew = new TH1D("smeared_rnd_ynew","smeared_rnd_ynew", n_smearedbins_incl, smearedBins_incl);   
  //TH1D *smeared_rnd_ynew_test = new TH1D("smeared_rnd_ynew_test","smeared_rnd_ynew_test", n_smearedbins_incl, smearedBins_incl);   
  TH1D *theory_rnd_ynew       = new TH1D("theory_rnd_ynew","theory_rnd_ynew",             n_thybins_incl    , thyBins_incl);   
  TH1D *smeared_rnd_ynew      = new TH1D("smeared_rnd_ynew","smeared_rnd_ynew",           n_smearedbins_incl, smearedBins_incl);   
  TH1D *smeared_rnd_ynew_test = new TH1D("smeared_rnd_ynew_test","smeared_rnd_ynew_test", n_smearedbins_incl, smearedBins_incl);   
  
  //RooUnfoldResponse response_ynew(smeared_rnd_ynew,theory_rnd_ynew);   
  //TH2* response_ynew_th2=(response_ynew.Hresponse());
  TH2D* response_ynew_th2=new TH2D("response_ynew_th2","response_ynew_th2",
				   (Int_t)n_smearedbins_incl, (Double_t*)smearedBins_incl,
				   (Int_t)n_thybins_incl, (Double_t*)thyBins_incl);
  
  
  ///////////////////////////////////////////////////////////////////////////////////////////// 
  //-----------------------------------------------------------------------------------------//
  //-----------------------------------------------------------------------------------------//
  //makeToySpectra(theory_ynew_spl3clone, spline3_ynew, fJER_ynew,              //-------------//
  //nEvents, theory_rnd_ynew, smeared_rnd_ynew, (RooUnfoldResponse*)&response_ynew); //-------------//
  //makeToySpectra(theory_ynew_spl3clone, spline3_ynew, fJER_ynew,              //-------------//
  //		 nEvents, theory_rnd_ynew, smeared_rnd_ynew, (TH2D*)response_ynew_th2); //-------------//
  makeToySpectra(theory_ynew_spl3clone, spline3_ynew, fJER_ynew,              //-------------//
		 nEvents, theory_rnd_ynew, smeared_rnd_ynew, smeared_rnd_ynew_test, (TH2D*)response_ynew_th2); //-------------//
  //-----------------------------------------------------------------------------------------//
  //-----------------------------------------------------------------------------------------//
  ///////////////////////////////////////////////////////////////////////////////////////////// 
  divideBinWidth(theory_rnd_ynew);  //for normalization only.
  theory_rnd_ynew->Scale(1./4.);//etabinwidth

  divideBinWidth(smeared_rnd_ynew); 
  smeared_rnd_ynew->Scale(1./4.);

  divideBinWidth(smeared_rnd_ynew_test); 
  smeared_rnd_ynew_test->Scale(1./4.);
  
  double norm_rnd_ynew=theory_ynew_spl3clone->Integral()/theory_rnd_ynew->Integral();
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
  
  std::string resptitle="response_ynew";//might do this facier later.
  RooUnfoldResponse response_ynew(smeared_rnd_ynew,theory_rnd_ynew,response_ynew_th2,resptitle.c_str()); 
  
  //now that roounfoldresponse is made, can divide hists by bin widths for the final smeared toy MC cross section
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
  
  Float_t maxy=theory_rnd_ynew->GetMaximum();    cout<<"maxy="<<maxy<<endl;
  maxy*=2.;  cout<<"maxy="<<maxy<<endl;
  Float_t miny=smeared_rnd_ynew->GetMinimum();   cout<<"miny="<<miny<<endl;
  miny/=2.; cout<<"miny="<<miny<<endl;
  if(miny<0.){
    cout<<"WARNING! Bin contents of smeared hist is NEGATIVE! Not possible!"<<endl;
    //miny=fabs(miny);
    miny=10.;  }
  theory_rnd_ynew->SetAxisRange(miny,maxy,"Y");  
  theory_rnd_ynew->SetTitle( ((std::string)("Toy MC Incl. NLO #sigma;;"+ddxsec_yax)).c_str() );
  
  
  ///// check plots True/Smeared 
  TCanvas *plot_ynew_true_smeared_rat = new TCanvas("plot_ynew_true_smeared_rat", "plot_ynew_true_smeared_rat",1400,1200);  
  plot_ynew_true_smeared_rat->cd();

  
  TPad* overlay_pad=new TPad("overlay_pad","pad for spectra",0.0,0.3,1.0,1.0);
  overlay_pad->SetLogx(1);
  overlay_pad->SetLogy(1);
  overlay_pad->SetBottomMargin(0);
  overlay_pad->Draw();//tells the canvas theres a pad at the coordinates specified in the new constructor
  overlay_pad->cd();
  
  theory_rnd_ynew->Draw("HIST E");  
  smeared_rnd_ynew->Draw("HIST E SAME");
  smeared_rnd_ynew_test->Draw("HIST E SAME");
  
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
  
  cout<<"creating ratio hist of true/smeared"<<endl;
  TH1D *ynew_true_smeared_rat=(TH1D*)theory_rnd_ynew->Clone("ynew_true_smeared_rat");   
  if(printBaseDebug)ynew_true_smeared_rat->Print("base");
  ynew_true_smeared_rat->SetTitle("; Jet p_T (GeV); True/Smeared");
  ynew_true_smeared_rat->Divide(ynew_true_smeared_rat,smeared_rnd_ynew,1.,1.,"B");   
  ynew_true_smeared_rat->GetXaxis()->SetMoreLogLabels(true);
  ynew_true_smeared_rat->GetXaxis()->SetNoExponent(true);
  ynew_true_smeared_rat->SetAxisRange(0.6,1.4,"Y");
  
  
  cout<<"creating ratio hist of true/opp side smeared"<<endl;
  TH1D *ynew_true_smeared_test_rat=(TH1D*)theory_rnd_ynew->Clone("ynew_true_smeared_test_rat");   
  if(printBaseDebug)ynew_true_smeared_test_rat->Print("base");  
  ynew_true_smeared_test_rat->SetTitle("; Jet p_T (GeV); True/Smeared");
  ynew_true_smeared_test_rat->Divide(ynew_true_smeared_test_rat,smeared_rnd_ynew_test,1.,1.,"B");   
  ynew_true_smeared_test_rat->GetXaxis()->SetMoreLogLabels(true);
  ynew_true_smeared_test_rat->GetXaxis()->SetNoExponent(true);
  ynew_true_smeared_test_rat->SetAxisRange(0.6,1.4,"Y");
  ynew_true_smeared_test_rat->SetLineColor(kRed-2);

  
  ynew_true_smeared_rat->Draw("HIST E");        
  ynew_true_smeared_test_rat->Draw("HIST E SAME");        
  
  //TLine* lineatone=new TLine(thyBins_incl[0],1.,thyBins_incl[n_thybins_incl],1.);
  lineatone->SetX2(thyBins_incl[n_thybins_incl]);
  lineatone->DrawClone();
  
  
  TH1D* theory_ynew_ratio_rnd_ynew=(TH1D*)theory_rnd_ynew->Clone ( 
								  ( (std::string)theory_rnd_ynew->GetName() +"_rathist") .c_str()
								   );
  theory_ynew_ratio_rnd_ynew->Divide(theory_ynew_spl3clone);
  theory_ynew_ratio_rnd_ynew->SetAxisRange(0.95,1.05,"Y");
  theory_ynew_ratio_rnd_ynew->GetYaxis()->SetTitle("Toy MC / orig. NLO");
    
    
  
    
    

    
    
    
  TCanvas* plot_response_ynew_th2=new TCanvas("plot_response_ynew_th2","plot_response_ynew_th2",1200, 1000);
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
  


  //double respmin=response_ynew_th2->GetMinimum();
  //std::cout<<"respmin="<<respmin<<std::endl;
  //double respmax=response_ynew_th2->GetMaximum();
  //std::cout<<"respmax="<<respmin<<std::endl;
  //response_ynew_th2->SetAxisRange(respmin/2., respmax*2., "Z");

  //debugcout();//0
  
  setTH2_ZAxisRange(response_ynew_th2);
  
  //debugcout();
  
  response_ynew_th2->Draw("COLZ");
  
  //////////////////////  END production of Smeared NLO spectra     //////////////////////    
  






  //////////////////////  START production of Smeared NP+NLO spectra     //////////////////////  
  cout<<"creating TH1 for toy NP+NLO spectra generation"<<endl<<endl;
  TH1D *theory_rnd_NPynew    = new TH1D("theory_rnd_NPynew","theory_rnd_NPynew", n_thybins_incl, thyBins_incl);   
  TH1D *smeared_rnd_NPynew = new TH1D("smeared_rnd_NPynew","smeared_rnd_NPynew", n_smearedbins_incl, smearedBins_incl);     
  TH1D *smeared_rnd_NPynew_test = new TH1D("smeared_rnd_NPynew_test","smeared_rnd_NPynew_test", n_smearedbins_incl, smearedBins_incl);     
  cout<<"creating instance of RooUnfoldResponse class (response_NPynew)"<<endl;
  
  //RooUnfoldResponse response_NPynew(smeared_rnd_NPynew,theory_rnd_NPynew); 
  //TH2* response_NPynew_th2=(response_NPynew.Hresponse());
  TH2D* response_NPynew_th2=new TH2D("response_NPynew_th2","response_NPynew_th2",
				     (Int_t)n_smearedbins_incl, (Double_t*)smearedBins_incl,
				     (Int_t)n_thybins_incl, (Double_t*)thyBins_incl);
  
  ///////////////////////////////////////////////////////////////////////////////////////////// 
  //-----------------------------------------------------------------------------------------//
  //-----------------------------------------------------------------------------------------//
  //makeToySpectra(theory_NPynew_spl3clone, spline3_NPynew, fJER_ynew,              //-------------//
  //nEvents, theory_rnd_NPynew, smeared_rnd_NPynew, (RooUnfoldResponse*)&response_NPynew); //-------------//
  //makeToySpectra(theory_NPynew_spl3clone, spline3_NPynew, fJER_ynew,              //-------------//
  //nEvents, theory_rnd_NPynew, smeared_rnd_NPynew, (TH2D*)response_NPynew_th2); //-------------//
  makeToySpectra(theory_NPynew_spl3clone, spline3_NPynew, fJER_ynew,              //-------------//
		 nEvents, theory_rnd_NPynew, smeared_rnd_NPynew, smeared_rnd_NPynew_test, (TH2D*)response_NPynew_th2); //-------------//
  //-----------------------------------------------------------------------------------------//
  //-----------------------------------------------------------------------------------------//
  ///////////////////////////////////////////////////////////////////////////////////////////// 
  divideBinWidth(theory_rnd_NPynew);  //for normalization only.
  theory_rnd_NPynew->Scale(1./4.);//etabinwidth

  divideBinWidth(smeared_rnd_NPynew); 
  smeared_rnd_NPynew->Scale(1./4.);

  divideBinWidth(smeared_rnd_NPynew_test); 
  smeared_rnd_NPynew_test->Scale(1./4.);

  double norm_rnd_NPynew=theory_NPynew_spl3clone->Integral()/theory_rnd_NPynew->Integral();
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
  
  std::string NPresptitle="response_NPynew";//might do this facier later.
  RooUnfoldResponse response_NPynew(smeared_rnd_NPynew,theory_rnd_NPynew,response_NPynew_th2,NPresptitle.c_str()); 

  divideBinWidth(theory_rnd_NPynew);  //for normalization only.
  theory_rnd_NPynew->Scale(1./4.);//etabinwidth

  divideBinWidth(smeared_rnd_NPynew); 
  smeared_rnd_NPynew->Scale(1./4.);

  divideBinWidth(smeared_rnd_NPynew_test); 
  smeared_rnd_NPynew_test->Scale(1./4.);

  theory_rnd_NPynew->SetLineColor(kRed);
  smeared_rnd_NPynew->SetLineColor(kBlue);    
  smeared_rnd_NPynew_test->SetLineColor(kBlue-5);    
  
  Float_t NPmaxy=theory_rnd_NPynew->GetMaximum();    cout<<"NPmaxy="<<NPmaxy<<endl;
  NPmaxy*=2.;  cout<<"NPmaxy="<<NPmaxy<<endl;
  Float_t NPminy=smeared_rnd_NPynew->GetMinimum();   cout<<"NPminy="<<NPminy<<endl;
  NPminy/=2.; cout<<"NPminy="<<NPminy<<endl;
  if(NPminy<0.){
    cout<<"WARNING! Bin contents of smeared hist is NEGATIVE! Not possible!"<<endl;
    //NPminy=fabs(NPminy);
    NPminy=10.;
  }
  theory_rnd_NPynew->SetAxisRange(NPminy,NPmaxy,"Y");
  theory_rnd_NPynew->SetTitle( ((std::string)("Toy MC Incl. NP+NLO #sigma;;"+ddxsec_yax)).c_str() );
  
  ///// check plots True/Smeared 
  TCanvas *plot_NPynew_true_smeared_rat = new TCanvas("plot_NPynew_true_smeared_rat", "plot_NPynew_true_smeared_rat",1400,1200);  
  plot_NPynew_true_smeared_rat->cd();
  
  TPad* NPoverlay_pad=new TPad("NPoverlay_pad","pad for spectra",0.0,0.3,1.0,1.0);
  NPoverlay_pad->SetLogx(1);
  NPoverlay_pad->SetLogy(1);
  NPoverlay_pad->SetBottomMargin(0);
  NPoverlay_pad->Draw();//tells the canvas theres a pad at the coordinates specified in the new constructor
  NPoverlay_pad->cd();
  
  theory_rnd_NPynew->Draw("HIST E");  
  smeared_rnd_NPynew->Draw("HIST E SAME");

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
  
  cout<<"creating ratio hist of true/smeared for NP+NLO"<<endl;
  TH1D *NPynew_true_smeared_rat=(TH1D*)theory_rnd_NPynew->Clone("NPynew_true_smeared_rat"); 
  if(printBaseDebug)NPynew_true_smeared_rat->Print("base");
  NPynew_true_smeared_rat->SetTitle("; Jet p_T (GeV); True/Smeared");
  NPynew_true_smeared_rat->Divide(NPynew_true_smeared_rat,smeared_rnd_NPynew,1.,1.,"B");     
  NPynew_true_smeared_rat->GetXaxis()->SetMoreLogLabels(true);
  NPynew_true_smeared_rat->GetXaxis()->SetNoExponent(true);
  NPynew_true_smeared_rat->SetAxisRange(0.6,1.4,"Y");
  
  
  cout<<"creating ratio hist of true/opp side smeared for NP+NLO"<<endl;
  TH1D *NPynew_true_smeared_test_rat=(TH1D*)theory_rnd_NPynew->Clone("NPynew_true_smeared_test_rat"); 
  if(printBaseDebug)NPynew_true_smeared_test_rat->Print("base");
  NPynew_true_smeared_test_rat->SetTitle("; Jet p_T (GeV); True/Smeared");
  NPynew_true_smeared_test_rat->Divide(NPynew_true_smeared_test_rat,smeared_rnd_NPynew_test,1.,1.,"B");     
  NPynew_true_smeared_test_rat->GetXaxis()->SetMoreLogLabels(true);
  NPynew_true_smeared_test_rat->GetXaxis()->SetNoExponent(true);
  NPynew_true_smeared_test_rat->SetAxisRange(0.6,1.4,"Y");
  NPynew_true_smeared_test_rat->SetLineColor(kRed-2);  



  NPynew_true_smeared_rat->Draw("HIST E");        
  NPynew_true_smeared_test_rat->Draw("HIST E SAME");        
  
  lineatone->DrawClone();
  
  
  TH1D* theory_NPynew_ratio_rnd_NPynew=(TH1D*)theory_rnd_NPynew->Clone ( 
									( (std::string)theory_rnd_NPynew->GetName() +"_rathist") .c_str()
									 );
  theory_NPynew_ratio_rnd_NPynew->Divide(theory_NPynew_spl3clone);
  theory_NPynew_ratio_rnd_NPynew->GetYaxis()->SetTitle("Toy MC / orig. NLO+NP");
  theory_NPynew_ratio_rnd_NPynew->SetAxisRange(0.95,1.05,"Y");





  

  TCanvas* plot_response_NPynew_th2=new TCanvas("plot_response_NPynew_th2","plot_response_NPynew_th2",1200, 1000);
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
  
  //double NPrespmin=response_NPynew_th2->GetMinimum();
  //std::cout<<"NPrespmin="<<NPrespmin<<std::endl;
  //double NPrespmax=response_NPynew_th2->GetMaximum();
  //response_NPynew_th2->SetAxisRange(NPrespmin/2., NPrespmax*2., "Z");
  setTH2_ZAxisRange(response_NPynew_th2);
  
  response_NPynew_th2->Draw("COLZ");
  //////////////////////  END production of Smeared NLO spectra     //////////////////////    

  




  


  
  
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
  
  // SPLINES
  theory_ynew_spl3clone->Write();
  spline3_ynew->Write();
  
  theory_NPynew_spl3clone->Write();
  spline3_NPynew->Write();
  
  plot_splines->Write();
  
  
  
  //TOY MC, JER SMEARED NLO SPECTRA
  theory_rnd_ynew->Write();
  theory_ynew_ratio_rnd_ynew->Write();
  
  smeared_rnd_ynew->Write(); 
  smeared_rnd_ynew_test->Write(); 
  ynew_true_smeared_rat->Write();  //ratio
  ynew_true_smeared_test_rat->Write();  //ratio
  plot_ynew_true_smeared_rat->Write();
  
  response_ynew.Write();//response matrix
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
  
  response_NPynew.Write();
  response_NPynew_th2->Write();
  plot_response_NPynew_th2->Write();  
  
  return;

  
  outf->Write();
  return;
}


//  steering ---------------------------------------------------------------------------------
int main(int argc, char* argv[]){  
  int rStatus = -1;   
  
  if(argc==2) {
    smearTheorySpectra_gaussCoreJER( (string) argv[1] );  
    rStatus=0;  }
  else 
    rStatus=1;  
  
  cout<<"rStatus="<<rStatus<<endl;
  return rStatus;
}



















/*
  TH1D *True_incl_jet_NPynew    = new TH1D("True_incl_jet_NPynew","True_incl_jet_NPynew", n_thybins_incl, thyBins_incl);   
  True_incl_jet_NPynew->Sumw2();
  if(printBaseDebug)True_incl_jet_NPynew->Print("base");
  
  TH1D *NormBinWidth_True_incl_jet_NPynew    = new TH1D("NormBinWidth_True_incl_jet_NPynew","NormBinWidth_True_incl_jet_NPynew", n_thybins_incl, thyBins_incl); 
  NormBinWidth_True_incl_jet_NPynew->Sumw2();
  if(printBaseDebug)NormBinWidth_True_incl_jet_NPynew->Print("base");
  
  //TH1D *Smeared_incl_jet_NPynew = new TH1D("Smeared_incl_jet_NPynew","Smeared_incl_jet_NPynew", n_thybins_incl, thyBins_incl);   
  TH1D *Smeared_incl_jet_NPynew = new TH1D("Smeared_incl_jet_NPynew","Smeared_incl_jet_NPynew", n_smearedbins_incl, smearedBins_incl);   
  Smeared_incl_jet_NPynew->Sumw2();
  if(printBaseDebug)Smeared_incl_jet_NPynew->Print("base");
  
  TH1D *NormBinWidth_Smeared_incl_jet_NPynew = new TH1D("NormBinWidth_Smeared_incl_jet_NPynew","NormBinWidth_Smeared_incl_jet_NPynew", n_smearedbins_incl, smearedBins_incl);   
  NormBinWidth_Smeared_incl_jet_NPynew->Sumw2();
  if(printBaseDebug)NormBinWidth_Smeared_incl_jet_NPynew->Print("base");
  
  cout<<"constructing RooUnfoldResponse for NPynew"<<endl;
  RooUnfoldResponse response_NPynew(Smeared_incl_jet_NPynew,True_incl_jet_NPynew);
  
  
  

 */
