#include "smearTheorySpectra.h"

const bool printBaseDebug=true;
const int nEvents=1e+08;  /// Increase your stats here, typically 100M events are enough 
//const int nEvents=1e+07;  /// debug nevents
std::string ddxsec_yax="#frac{d^{2}#sigma}{dp_{T}dy} [unit]";

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
  TFile *f1 = new TFile(inputFile.c_str());  
  cout<<endl;
  
  /////////////// NPs to NLO jet spectra ////////////////////////////////////////////////////////
  cout<<"plotting NP fits! " << endl;
  
  /////////////// Get NP fit functions
  TF1 *fNP_y0 = new TF1("fNP_y0",NP_y0_str.c_str(), thyBins_incl[0], thyBins_incl[n_thybins_incl]);
  TF1 *fNP_y1 = new TF1("fNP_y1",NP_y1_str.c_str(), thyBins_incl[0], thyBins_incl[n_thybins_incl]);
  TF1 *fNP_y2 = new TF1("fNP_y2",NP_y2_str.c_str(), thyBins_incl[0], thyBins_incl[n_thybins_incl]);
  TF1 *fNP_y3 = new TF1("fNP_y3",NP_y3_str.c_str(), thyBins_incl[0], thyBins_incl[n_thybins_incl]);
  TF1 *fNP_ynew = new TF1("fNP_ynew",NP_ynew_str.c_str(), thyBins_incl[0], thyBins_incl[n_thybins_incl]);  //avg of other four
  //-----------------------------------

  /////////////// plots to to check NPs 
  // 2x2 canv of |y| bins
  TCanvas *plot_y0_to_y3_NPs = new TCanvas("plot_y0_to_y3_NPs", "plot_y0_to_y3_NPs",1200,1000);
  plot_y0_to_y3_NPs->Divide(2,2);
  //-----------------------------------
  plot_y0_to_y3_NPs->cd(1)->SetLogx(1);
  plot_y0_to_y3_NPs->cd(1);
  TH1F *hNP_y0 = (TH1F*)( (TH1F*)fNP_y0->GetHistogram()
			  )->Clone("hNP_y0");
  hNP_y0->SetTitle("NP Corr Fit for 0.0 < #||{y} < 0.5;Jet p_{T};Fit Val");
  hNP_y0->GetXaxis()->SetNoExponent(true);
  hNP_y0->GetXaxis()->SetMoreLogLabels(true);
  hNP_y0->SetLineColor(kRed);
  hNP_y0->DrawClone();  
  //-----------------------------------  
  plot_y0_to_y3_NPs->cd(2)->SetLogx(1);
  plot_y0_to_y3_NPs->cd(2);
  TH1F *hNP_y1 = (TH1F*)( (TH1F*)fNP_y1->GetHistogram()
			  )->Clone("hNP_y1");  
  hNP_y1->SetTitle("NP Corr Fit for 0.5 < #||{y} < 1.0;Jet p_{T};Fit Val");
  hNP_y1->GetXaxis()->SetNoExponent(true);
  hNP_y1->GetXaxis()->SetMoreLogLabels(true);
  hNP_y1->SetLineColor(kBlue);
  hNP_y1->DrawClone();  
  //-----------------------------------
  plot_y0_to_y3_NPs->cd(3)->SetLogx(1);
  plot_y0_to_y3_NPs->cd(3);
  TH1F *hNP_y2 = (TH1F*)( (TH1F*)fNP_y2->GetHistogram()
			  )->Clone("hNP_y2");  
  hNP_y2->SetTitle("NP Corr Fit for 1.0 < #||{y} < 1.5;Jet p_{T};Fit Val");
  hNP_y2->GetXaxis()->SetNoExponent(true);
  hNP_y2->GetXaxis()->SetMoreLogLabels(true);
  hNP_y2->SetLineColor(kGreen);
  hNP_y2->DrawClone();  
  //-----------------------------------
  plot_y0_to_y3_NPs->cd(4)->SetLogx(1);
  plot_y0_to_y3_NPs->cd(4);
  TH1F *hNP_y3 = (TH1F*)( (TH1F*)fNP_y3->GetHistogram()
			  )->Clone("hNP_y3");  
  hNP_y3->SetTitle("NP Corr Fit for 1.5 < #||{y} < 2.0;Jet p_{T};Fit Val");
  hNP_y3->GetXaxis()->SetNoExponent(true);
  hNP_y3->GetXaxis()->SetMoreLogLabels(true);
  hNP_y3->SetLineColor(kMagenta);
  hNP_y3->DrawClone();    
  //-----------------------------------
  

  // |y|<2.0
  TCanvas *plot_totNPs = new TCanvas("plot_totNPs", "plot_totNPs",600,500);//1600,1200);
  plot_totNPs->SetLogx(1);
  plot_totNPs->cd();
  //-----------------------------------
  TH1F *hNP_ynew = (TH1F*)( (TH1F*)fNP_ynew->GetHistogram()
			    )->Clone("hNP_ynew");  
  hNP_ynew->SetTitle("Avg. NP Corr Fit for 0.0 < #||{y} < 2.0;Jet p_{T};Fit Val");
  hNP_ynew->GetXaxis()->SetNoExponent(true);
  hNP_ynew->GetXaxis()->SetMoreLogLabels(true);
  hNP_ynew->SetLineColor(kTeal);
  hNP_ynew->DrawClone();  //plot_totNPs->Update();  //why did this need to happen?
  //-----------------------------------  
  
  
  // |y|<2.0 + 0.5 size |y| bins on same canv
  TCanvas *plot_allNPs = new TCanvas("plot_allNPs", "plot_allNPs",600,500);
  plot_allNPs->SetLogx(1);
  plot_allNPs->cd();
  //-----------------------------------  
  hNP_ynew->SetTitle("NP fits all #||{y} bins");
  hNP_ynew->DrawClone();
  hNP_y0  ->DrawClone("SAME");
  hNP_y1  ->DrawClone("SAME");
  hNP_y2  ->DrawClone("SAME");
  hNP_y3  ->DrawClone("SAME");
  //-----------------------------------    
  TLegend* leg_allNPs=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
  leg_allNPs->AddEntry(hNP_ynew , "#||{y} < 2.0 (avg of others)" , "l");
  leg_allNPs->AddEntry(hNP_y0   , "0.0 < #||{y} < 0.5" , "l");
  leg_allNPs->AddEntry(hNP_y1   , "0.5 < #||{y} < 1.0" , "l");
  leg_allNPs->AddEntry(hNP_y2   , "1.0 < #||{y} < 1.5" , "l");
  leg_allNPs->AddEntry(hNP_y3   , "1.5 < #||{y} < 2.0" , "l");
  leg_allNPs->Draw();
  //-----------------------------------  
  

  
  
  /////////////// Core p_T Resolution ////////////////////////////////////////////////////////
  cout<<"plotting JER!"<<  endl;  //cout<<""<<  << endl;
  
  /////////////// Create/Get JER fit function(s)
  TF1 *fJER_ynew = new TF1("fJER_ynew",gJER_ynew_str.c_str(),thyBins_incl[0], thyBins_incl[n_thybins_incl]);  
  fJER_ynew->SetLineColor(2);
  fJER_ynew->SetMinimum(0); 
  fJER_ynew->SetMaximum(0.25);
  //-----------------------------------  
  
  /////////////// Create plot to check
  TH1D *hJER_ynew = (TH1D*)(
			    (TH1D*)fJER_ynew->GetHistogram()
			    )->Clone("hJER_ynew");//TF1 is irritating to set titles for+blah... so use this
  hJER_ynew->SetTitle(" Gauss Core JER Fit #||{y} < 2.0 ; Jet p_{T} ; #sigma / #mu from Fit");//sets hist title+xaxis+yaxis title
  hJER_ynew->GetXaxis()->SetNoExponent(true);
  hJER_ynew->GetXaxis()->SetMoreLogLabels(true);
  //-----------------------------------  
  TCanvas *plotJER = new TCanvas("plotJER", "plotJER",700,550);
  plotJER->cd()->SetLogx(1);
  plotJER->cd();  
  hJER_ynew->DrawClone();    
  //REPLACE ME W/TEXT SAYING SAME THING
  //TLegend *leg=new TLegend(0.53,0.50,0.85,0.85);
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
  TH1D *theory_y0 = (TH1D*)f1->Get("h0100100");
  theory_y0->Sumw2(true);
  if(printBaseDebug)  theory_y0->Print("Base");
  //theory_y0->Scale(1.e-03);//picobarns to nanobarns
  TH1D *theory_y1 = (TH1D*)f1->Get("h0100200");
  theory_y1->Sumw2(true);
  if(printBaseDebug)  theory_y1->Print("Base");
  //theory_y1->Scale(1.e-03);
  TH1D *theory_y2 = (TH1D*)f1->Get("h0100300");
  theory_y2->Sumw2(true);
  if(printBaseDebug)  theory_y2->Print("Base");
  //theory_y2->Scale(1.e-03);
  TH1D *theory_y3 = (TH1D*)f1->Get("h0100400");
  theory_y3->Sumw2(true);
  if(printBaseDebug)  theory_y3->Print("Base");
  //theory_y3->Scale(1.e-03);
  cout<<endl;
  //-----------------------------------
  cout <<"making 00eta20 hist from theory curves..."<<endl<<endl;
  TH1D *theory_ynew = make00eta20Hist( (TH1D*) theory_y0, 
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
  theory_ynew->SetTitle("NLO Inclusive #sigma_{jet} all #||{y} bins");
  Double_t max=(Double_t)theory_ynew->GetMaximum();
  Double_t min=(Double_t)theory_y3->GetMinimum();
  theory_ynew->SetAxisRange((min/2.),(max*2.),"Y");
  theory_ynew->DrawClone("HIST E");
  theory_y0->DrawClone("HIST E SAME");
  theory_y1->DrawClone("HIST E SAME");
  theory_y2->DrawClone("HIST E SAME");
  theory_y3->DrawClone("HIST E SAME");  
  //-------------------------------------------- 
  TLegend* leg_allNLOxsec=new TLegend(0.60, 0.70, 0.9, 0.9, NULL,"BRNDC");
  //leg_allNPs->AddEntry(hNP_ynew , "#||{y} < 2.0 (sum of others)" , "l");
  leg_allNLOxsec->AddEntry(theory_ynew   , "(0.0 < #||{y} < 2.0) x 10^{4}" , "lp");
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
  //--------------------------------------------
  
  
  
  
  
  //////////// Apply NPs to theory ////////////////////////////////////////////////////////
  cout<<"plotting NP+NLO xsecs!"<<  endl;  //cout<<""<<  << endl;

  /// Cross section ynew
  TH1D *theory_NPy0=(TH1D*)applyNPtoxsec(theory_y0, fNP_y0);
  TH1D *theory_NPy1=(TH1D*)applyNPtoxsec(theory_y1, fNP_y1);
  TH1D *theory_NPy2=(TH1D*)applyNPtoxsec(theory_y2, fNP_y2);
  TH1D *theory_NPy3=(TH1D*)applyNPtoxsec(theory_y3, fNP_y3);

  TH1D *theory_NPynew=(TH1D*)applyNPtoxsec(theory_ynew, fNP_ynew);
  

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
  theory_NPynew->SetTitle("NP+NLO Inclusive #sigma_{jet} all #||{y} bins");
  //Double_t NPmax=(Double_t)theory_NPy0->GetMaximum();
  Double_t NPmax=(Double_t)theory_NPynew->GetMaximum();
  Double_t NPmin=(Double_t)theory_NPy3->GetMinimum();
  theory_NPynew->SetAxisRange((NPmin/2.),(NPmax*2.),"Y");
  theory_NPynew->DrawClone("HIST E");
  theory_NPy0->DrawClone("HIST E SAME");
  theory_NPy1->DrawClone("HIST E SAME");
  theory_NPy2->DrawClone("HIST E SAME");
  theory_NPy3->DrawClone("HIST E SAME");  
  //-------------------------------------------- 
  TLegend* leg_allNPNLOxsec=new TLegend(0.60, 0.70, 0.9, 0.9, NULL,"BRNDC");
  //leg_allNPs->AddEntry(hNP_ynew , "#||{y} < 2.0 (sum of others)" , "l");
  leg_allNPNLOxsec->AddEntry(theory_NPynew   , "(0.0 < #||{y} < 2.0) x 10^{4}" , "lp");
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
  //--------------------------------------------





  

  /////////////// Fitting Spectra with spline3
  /////////////// Not worning well for y6 bin because Klaus's spectra ends < 500 GeV
  /////////////// We must fit spectra by hand for y6 bin
  

  ///Create Cubic Splines using Cross sections
  cout<<"creating TSpline of cross section hist"<<endl;
  TH1D* theory_ynew_spl3clone=(TH1D*)theory_ynew->Clone(  
							( ( (std::string) theory_ynew->GetName() ) +"_spline3").c_str() 
							  );
  theory_ynew_spl3clone->SetLineColor(kCyan+4);
  multiplyBinWidth(theory_ynew_spl3clone);//pt bin width(s)
  
  //TSpline3 *spline3_ynew = new TSpline3( theory_ynew_spl3clone);
  TSpline3 *spline3_ynew = new TSpline3( theory_ynew_spl3clone,"",thyBins_incl[0],thyBins_incl[n_thybins_incl]);
  spline3_ynew->SetLineColor(kAzure);  
  
  
  TH1D* theory_NPynew_spl3clone=(TH1D*)theory_NPynew->Clone(  
							    ( ( (std::string) theory_NPynew->GetName() ) +"_spline3").c_str() 
							      );  
  theory_NPynew_spl3clone->SetLineColor(kCyan-6);
  multiplyBinWidth(theory_NPynew_spl3clone);
  
  TSpline3 *spline3_NPynew = new TSpline3( theory_NPynew_spl3clone);
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
  leg_spline->AddEntry(theory_ynew , "Weighted NLO Jet Counts for #||{y} < 2.0" , "lp");
  leg_spline->AddEntry(spline3_ynew   , "Cubic Spline Fit for NLO #||{y} < 2.0" , "l");
  leg_spline->AddEntry(theory_NPynew , "Weighted NP+NLO Jet Counts for #||{y} < 2.0" , "lp");
  leg_spline->AddEntry(spline3_NPynew   , "Cubic Spline Fit for NP+NLO #||{y} < 2.0" , "l");
  leg_spline->Draw();
  
  
  
  
  //////////////////////  Start production of Smeared spectra
  int tenth_nEvents=nEvents/10;//for debug text output
  
  ////intuitive option, problem; neg weight from TSpline for ptTrue>the center of last bin
  double ptmin_thy = thyBins_incl[0];
  double ptmax_thy = thyBins_incl[n_thybins_incl];
  //double ptmax_thy_splinevalid=thyBins_incl[n_thybins_incl] - 0.5*(thyBins_incl[n_thybins_incl] - thyBins_incl[n_thybins_incl-1]);//i'm sorry scientific integrity....
  //double ptmax_thy_splinevalid=877.3827;//found painfully by hand... 
  // double ptmax_thy_splinevalid= 877.382702298;//found painfully by hand...
  double ptmax_thy_splinevalid=  877.3827022981749337304790;
  ////over-correction option, problem; avoids neg weight, but last bin for smeared NLO spectra is basically lost. might be palatable because we're stat-limited there anyways.
  //double ptmin_thy = thyBins_incl[1];
  //double ptmax_thy = thyBins_incl[n_thybins_incl-1];
  //slightly-less-over-correction problem; avoids neg weight, maximises pt range possible w/ current NLO spectra, but last bin probably still not going to be usable. slightly prefer
  //double ptmin_thy = thyBins_incl[0]+0.5*(thyBins_incl[1]-thyBins_incl[0]);
  //double ptmax_thy = thyBins_incl[n_thybins_incl]-0.5*(thyBins_incl[n_thybins_incl]-thyBins_incl[n_thybins_incl-1]);  
  
  double ptmin_smeared = smearedBins_incl[0];
  double ptmax_smeared = smearedBins_incl[n_smearedbins_incl];
  
  TRandom3 *rnd = new TRandom3();
  
  //NLO smearing first
  cout<<"creating TH1 for toy NLO spectra generation, RooUnfoldResponse class, etc."<<endl<<endl;
  TH1D *theory_rnd_ynew    = new TH1D("theory_rnd_ynew","theory_rnd_ynew", n_thybins_incl, thyBins_incl);   
  TH1D *smeared_rnd_ynew = new TH1D("smeared_rnd_ynew","smeared_rnd_ynew", n_smearedbins_incl, smearedBins_incl);   
  RooUnfoldResponse response_ynew(smeared_rnd_ynew,theory_rnd_ynew); 
  
  int response_count=0, fake_count=0, miss_count=0;
  int debug_ptw_lt0_count=0, debug_ptTrue_gt843_count=0;//, debug_ptTrue_lt46_count=0;
  for(int i=0;i<nEvents;++i){      
    
    if(i%tenth_nEvents==0)
      cout<<"throwing random #'s for event # "<<i<<endl;
    
    double ptTrue  = rnd->Uniform(ptmin_thy,ptmax_thy);
    double pt_w      =  spline3_ynew->Eval(ptTrue);
    double sigma   = gJER_ynew(ptTrue); 
    //sigma*=1.079; //JER Scaling factor 8 TeV  Acoording to Mikko no scaling for the moment
    double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
    
    
    if(ptTrue>ptmax_thy_splinevalid){
      debug_ptTrue_gt843_count++;
      if(pt_w<0.){
	//cout<<"WARNING!!!! pt_w is NEGATIVE!!!!"<<endl;//happens because we ask spline for a value beyond last bin center.
	//cout<<"replacing w/ weight of largest possible pt value!! "<<ptmax_thy_splinevalid<<endl
	debug_ptw_lt0_count++;
	///cout<<"pt_w was = "<<pt_w<<endl;
	pt_w=spline3_ynew->Eval(ptmax_thy_splinevalid);
	//cout<<"pt_w now = "<<pt_w<<endl;      
      }    }      
    
    
    bool in_smearpt_range=(  ( ptSmeared>ptmin_smeared ) && ( ptSmeared<ptmax_smeared )  );
    bool in_trupt_range=(  ( ptTrue>ptmin_thy )        && ( ptTrue<ptmax_thy )  )  ;
    if (  in_smearpt_range && in_trupt_range){
      if(response_count%tenth_nEvents==0){
	cout<<"filling response w/ jet, response_count="<<response_count<<endl;
	cout<< "ptTrue / ptSmeared = "<< ptTrue << " / " << ptSmeared << endl;      }
      response_ynew.Fill(ptSmeared,ptTrue,pt_w); 
      response_count++;          }
    else if (  in_trupt_range ){
      if(miss_count%100000==0){
	cout<<"filling response w/ missed jet, miss_count="<<miss_count<<endl;
	cout<< "ptTrue / ptSmeared = "<< ptTrue << " / " << ptSmeared << endl;      }
      response_ynew.Miss(ptTrue,pt_w); 
      miss_count++;    }
    else if ( in_smearpt_range) {
      if(fake_count%100000==0){
	cout<<"filling response w/ 'fake' jet, fake_count="<<fake_count<<endl;
	cout<< "ptTrue / ptSmeared = "<< ptTrue << " / " << ptSmeared << endl;      }
      response_ynew.Fake(ptSmeared,pt_w); 
      fake_count++;          }
    
    if(in_trupt_range)
      theory_rnd_ynew    -> Fill(ptTrue,pt_w);
    if(in_smearpt_range)
      smeared_rnd_ynew -> Fill(ptSmeared,pt_w);
    
    
  }
  cout<<"done smearing ynew."<<endl;  
  cout<<"smear summary;"<<endl;
  cout<<"nEvents="<<nEvents<<endl;
  cout<<"response_count="<<response_count<<endl;
  cout<<"miss_count="<<miss_count<<endl;
  cout<<"fake_count="<<fake_count<<endl;
  
  cout<<"debug_ptTrue_gt843_count="<<debug_ptTrue_gt843_count<<endl;
  cout<<"debug_ptw_lt0_count="<<debug_ptw_lt0_count<<endl;
  
  
  
  
  divideBinWidth(theory_rnd_ynew);
  theory_rnd_ynew->SetLineColor(kRed);
  divideBinWidth(smeared_rnd_ynew);
  smeared_rnd_ynew->SetLineColor(kBlue);
  
  Float_t maxy=theory_rnd_ynew->GetMaximum();    cout<<"maxy="<<maxy<<endl;
  maxy*=2.;  cout<<"maxy="<<maxy<<endl;
  Float_t miny=smeared_rnd_ynew->GetMinimum();   cout<<"miny="<<miny<<endl;
  miny/=2.; cout<<"miny="<<miny<<endl;
  theory_rnd_ynew->SetAxisRange(miny,maxy,"Y");
  //theory_rnd_ynew->SetTitle( ((std::string)("Toy MC Incl. NLO #sigma;;"+ddxsec_yax)).c_str() );
  
  
  
  cout<<"creating ratio hist of true/smeared"<<endl;
  TH1D *ynew_true_smeared_rat=(TH1D*)theory_rnd_ynew->Clone("ynew_true_smeared_rat"); 
  if(printBaseDebug)ynew_true_smeared_rat->Print("base");
  ynew_true_smeared_rat->SetTitle("; Jet p_T (GeV); True/Smeared          ");
  ynew_true_smeared_rat->Divide(ynew_true_smeared_rat,smeared_rnd_ynew,1.,1.,"B");   
  




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

  TLegend* leg_smear=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
  leg_smear->AddEntry(theory_rnd_ynew  , "Truth Incl. NLO #sigma" , "lp" );
  leg_smear->AddEntry(smeared_rnd_ynew , "Smeared Incl. NLO #sigma" , "lp" );
  leg_smear->Draw();
  
  plot_ynew_true_smeared_rat->cd();//go back to main canvas before doing new pad
  
  TPad* ratio_pad=new TPad("ratio_pad","pad for ratio",0.0,0.05,1.0,0.30);
  ratio_pad->SetLogx(1);
  ratio_pad->SetLogy(0);
  ratio_pad->SetTopMargin(0);
  //ratio_pad->SetBottomMargin(0.3);
  ratio_pad->Draw();//tell canvas there's a pad at the specified coordinates in the new constructor
  ratio_pad->cd();
  ynew_true_smeared_rat->GetXaxis()->SetMoreLogLabels(true);
  ynew_true_smeared_rat->GetXaxis()->SetNoExponent(true);
  ynew_true_smeared_rat->SetAxisRange(0.6,1.4,"Y");
  ynew_true_smeared_rat->Draw("HIST E");        
  
  



  TH2* response_ynew_th2=(response_ynew.Hresponse());
  TCanvas* plot_response_ynew_th2=new TCanvas("plot_response_ynew_th2","plot_response_ynew_th2",1200, 1000);
  plot_response_ynew_th2->cd()->SetLogx(1);
  plot_response_ynew_th2->cd()->SetLogy(1);
  plot_response_ynew_th2->cd()->SetLogz(1);
  plot_response_ynew_th2->cd();
  
  response_ynew_th2->SetTitle("Toy MC Response Matrix from JER Smeared NLO #sigma_{jet}");
  response_ynew_th2->GetXaxis()->SetTitle("RECO Jet p_{T} [GeV]");
  response_ynew_th2->GetYaxis()->SetTitle("GEN Jet p_{T} [GeV]");
  response_ynew_th2->Draw("COLZ");
  
  
  
  
  
  
  
  TH1D *ynew_true_smeared_rat_test=(TH1D*)theory_rnd_ynew->Clone("ynew_true_smeared_rat_test"); 
  if(printBaseDebug)ynew_true_smeared_rat_test->Print("base");
  ynew_true_smeared_rat_test->SetTitle("Inclusive NLO #sigma_{jet}, True/Smeared TEST HIST ratio; Jet p_T (GeV); True/Smeared          ");
  ynew_true_smeared_rat_test->Divide(smeared_rnd_ynew);  
  
  ///// alt ratio plot test; want to see how the lack of the "B" option in divide changes the answer
  TCanvas *plot_ynew_true_smeared_rat_test = new TCanvas("plot_ynew_true_smeared_rat_test", "plot_ynew_true_smeared_rat_test",1600,800);  
  plot_ynew_true_smeared_rat_test->SetLogx(1);
  plot_ynew_true_smeared_rat_test->SetLogy(0);
  plot_ynew_true_smeared_rat_test->cd();
  ynew_true_smeared_rat_test->SetAxisRange(0.6,1.4,"Y");
  ynew_true_smeared_rat_test->Draw("HIST E");        
  
  
  
  //rebin 4 ratio only
  //cout<<"rebinning..."<<endl;
  //ynew_true_smeared_rat=(TH1D*)ynew_true_smeared_rat->TH1::Rebin(n_smearedbins_incl, "theory_rnd_ynew_clone_ratio_rebin", smearedBins_incl)    ;
  //if(printBaseDebug)ynew_true_smeared_rat->Print("base");  
  //cout<<"dividing..."<<endl;
  //ynew_true_smeared_rat->Divide(theory_rnd_ynew,smeared_rnd_ynew,1.,1.,"B");  
  //ynew_true_smeared_rat->Divide(smeared_rnd_ynew);  
  //ynew_true_smeared_rat->SetAxisRange(smearedBins_incl[0],smearedBins_incl[n_smearedbins_incl],"X"); 
  
  
      
  //NP+NLO smearing next
  cout<<"creating TH1 for toy NP+NLO spectra generation"<<endl<<endl;
  TH1D *theory_rnd_NPynew    = new TH1D("theory_rnd_NPynew","theory_rnd_NPynew", n_thybins_incl, thyBins_incl);   
  TH1D *smeared_rnd_NPynew = new TH1D("smeared_rnd_NPynew","smeared_rnd_NPynew", n_smearedbins_incl, smearedBins_incl);   
  
  cout<<"creating instance of RooUnfoldResponse class (response_NPynew)"<<endl;
  RooUnfoldResponse response_NPynew(smeared_rnd_NPynew,theory_rnd_NPynew); 
  
  
  




  


  
  
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
  spline3_ynew->Write();
  spline3_NPynew->Write();
  plot_splines->Write();

  
  // RESPONSE MATRICES GENERATED
  //NLO smeared mat
  response_ynew.Write("response_ynew");
  response_ynew_th2->Write();
  plot_response_ynew_th2->Write();

  //NLO smeared/truth response spectra
  theory_rnd_ynew->Write();
  smeared_rnd_ynew->Write(); 
  ynew_true_smeared_rat->Write();  

  plot_ynew_true_smeared_rat->Write();

  ynew_true_smeared_rat_test->Write();  
  plot_ynew_true_smeared_rat_test->Write();



  //NP+NLO smeared mat
  //NP+NLO smeared/truth response spectra
  
  
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
