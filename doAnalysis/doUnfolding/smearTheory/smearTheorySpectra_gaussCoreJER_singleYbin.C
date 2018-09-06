#include "smearTheorySpectra.h"

const bool printBaseDebug=true;
//const int nEvents=1e+08;  /// Increase your stats here, typically 100M events are enough 
const int nEvents=1e+07;  /// debug nevents
std::string ddxsec_yax="#frac{d^{2}#sigma}{dp_{T}dy} [unit]";

int debugInt=0;
void debugcout(){
  std::cout<<"location #"<<debugInt<<endl;
  debugInt++;
  return; }

void smearTheorySpectra_gaussCoreJER_singleYbin( string inputString, int etaBin ){
  
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  std::string etaBin_str=std::to_string(etaBin);
  
  // Output File to write to
  //string outputFile= inputString+"_gaussSmear_00eta20.root";
  string outputFile= inputString+"_gaussSmear_etaBin"+etaBin_str+".root";
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
  std::string NPstr;
  if(etaBin==1) NPstr=NP_y0_str;
  else if(etaBin==2)NPstr=NP_y1_str;
  else if(etaBin==3)NPstr=NP_y2_str;
  else if(etaBin==4)NPstr=NP_y3_str;
  else assert(false);
  std::string NPname="fNP_y"+std::to_string(etaBin-1);

  /////////////// get NLO xsec hists from file
  //doing this earlier in the singleYbin version because i need the pt range explicitly. 
  std::string histname="h0100"+etaBin_str+"00";
  TH1D *theory = (TH1D*)f1->Get(histname.c_str());
  if(!theory)assert(false);
  theory->Sumw2(true);
  if(printBaseDebug)  theory->Print("Base");
  //theory_y0->Scale(1.e-03);//picobarns to nanobarns
  
  if(etaBin==1||etaBin==2){//fix the last-huge-bin issue by rebinning to same binning minus that one bin.
    theory=(TH1D*)theory->TH1::Rebin(n_thybins_incl2, theory->GetName(), thyBins_incl2);    
  }  
  double ptmin_thy=theory->GetBinLowEdge(1);
  double ptmax_thy=theory->GetBinLowEdge(theory->GetNbinsX()) + theory->GetBinWidth(theory->GetNbinsX()) ;
  
  TF1 *fNP = new TF1(NPname.c_str(),NPstr.c_str(), ptmin_thy, ptmax_thy);
  
  /////////////// plots to to check NP
  TCanvas *plot_NP = new TCanvas( ("plot_NP_y"+etaBin_str).c_str(), 
				  ("plot_NP_y"+etaBin_str).c_str(), 800, 600);
  plot_NP->cd()->SetLogx(1);
  plot_NP->cd();
  TH1F *hNP = (TH1F*)( (TH1F*)fNP->GetHistogram()
		       )->Clone( ("hNP_y"+etaBin_str).c_str());
  std::string hNP_title="NP Corr Fit #||{y} bin #" + etaBin_str+";Jet p_{T};Fit Val";
  hNP->SetTitle(hNP_title.c_str());
  hNP->GetXaxis()->SetNoExponent(true);
  hNP->GetXaxis()->SetMoreLogLabels(true);
  hNP->SetLineColor(kRed);
  hNP->DrawClone();       
  
  
  /////////////// Core p_T Resolution ////////////////////////////////////////////////////////
  cout<<"plotting JER!"<<  endl;  //cout<<""<<  << endl;
  
  /////////////// Create/Get JER fit function(s)
  //TF1 *fJER_ynew = new TF1("fJER_ynew",gJER_ynew_str.c_str(),thyBins_incl[0], thyBins_incl[n_thybins_incl]);  
  TF1 *fJER_ynew = new TF1("fJER_ynew",gJER_ynew_str.c_str(),thyBins_incl[0], 1000); //temp fix until i come up with something smarter 
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
  
  
  
  
  std::string plotname="plot_NLOxsec_y"+etaBin_str;
  TCanvas *plot_NLOxsec = new TCanvas(plotname.c_str(), plotname.c_str(),800,600);
  //--------------------------------------------
  plot_NLOxsec->cd()->SetLogx(1);
  plot_NLOxsec->cd()->SetLogy(1);
  plot_NLOxsec->cd();
  theory->SetTitle(("NLO #sigma_{jet}, #||{y} bin #"+etaBin_str+";Jet p_{T};"+ddxsec_yax).c_str());
  theory->GetXaxis()->SetNoExponent(true);
  theory->GetXaxis()->SetMoreLogLabels(true);
  theory->SetLineColor(kRed);
  theory->DrawClone("HIST E");    
  //--------------------------------------------
  
  //////////// Apply NPs to theory ////////////////////////////////////////////////////////
  cout<<"plotting NP+NLO xsecs!"<<  endl;  //cout<<""<<  << endl;
  
  /// Cross section ynew
  TH1D *theory_NP=(TH1D*)applyNPtoxsec(theory, fNP);
  
  //2X2 CANV THEORY NPNLO xsec
  std::string NPplotname="plot_NPNLOxsec_y"+etaBin_str;
  TCanvas *plot_NPNLOxsec = new TCanvas(NPplotname.c_str(),NPplotname.c_str(),800,600);
  //--------------------------------------------
  plot_NPNLOxsec->cd()->SetLogx(1);
  plot_NPNLOxsec->cd()->SetLogy(1);
  plot_NPNLOxsec->cd();
  theory_NP->DrawClone("HIST E");    
  //--------------------------------------------
  
  
  ///Create Cubic Splines using Cross sections
  cout<<"creating TSpline of cross section hist"<<endl;
  TH1D* theory_spl3clone=(TH1D*)theory->Clone(  
					      ( ( (std::string) theory->GetName() ) +"_spline3").c_str() 
						);
  //theory_ynew_spl3clone=(TH1D*)theory_ynew_spl3clone->TH1::Rebin(n_thybins_incl2, 
  //(((std::string)theory_ynew_spl3clone->GetName()) + "_specialrebin").c_str(), 
  //thyBins_incl2 );
  theory_spl3clone->SetLineColor(kCyan+4);    
  
  TSpline3 *spline3 = new TSpline3( theory_spl3clone);
  spline3->SetLineColor(kAzure);
  
  
  
  TH1D* theory_NP_spl3clone=(TH1D*)theory_NP->Clone(  
						    ( ( (std::string) theory_NP->GetName() ) +"_spline3").c_str() 
						      );  
  theory_NP_spl3clone->SetLineColor(kCyan-6);
  
  TSpline3 *spline3_NP = new TSpline3( theory_NP_spl3clone);
  spline3_NP->SetLineColor(kAzure-8);  
  
  /// check spectra with splines
  cout<<"drawing Tspline and cross section on canvas"<<endl;
  std::string splinecanvtitle="plot_splines_y"+etaBin_str;
  TCanvas *plot_splines = new TCanvas(splinecanvtitle.c_str(),splinecanvtitle.c_str(),1200,800);
  plot_splines->cd()->SetLogx(1);
  plot_splines->cd()->SetLogy(1);
  plot_splines->cd();
  
  //draw
  std::string splinetitle="Cubic Spline Fits #||{y} bin #"+etaBin_str+";Jet p_{T};Smear Weight";
  theory_spl3clone->SetTitle(splinetitle.c_str());
  theory_spl3clone->DrawClone("HIST E");
  spline3->Draw("SAME");
  theory_NP_spl3clone->DrawClone("HIST E SAME");  
  spline3_NP->Draw("SAME");
  
  
  TLegend* leg_spline=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
  leg_spline->AddEntry(theory ,    (    "Weighted NLO Jet Counts for #||{y} bin #" + etaBin_str ).c_str(), "lp");
  leg_spline->AddEntry(spline3   , (       "Cubic Spline Fit for NLO #||{y} bin #" + etaBin_str ).c_str(), "l");
  leg_spline->AddEntry(theory_NP , ( "Weighted NP+NLO Jet Counts for #||{y} bin #" + etaBin_str ).c_str(), "lp");
  leg_spline->AddEntry(spline3_NP, (    "Cubic Spline Fit for NP+NLO #||{y} bin #" + etaBin_str ).c_str(), "l");
  leg_spline->Draw();
  
  
  






























  
  

  
  
  
  
  
  
  
  
  //////////////////////  START production of Smeared NLO spectra     //////////////////////  
  cout<<"creating TH1 for toy NLO spectra generation, RooUnfoldResponse class, etc."<<endl<<endl;
  //TH1D *theory_rnd    = new TH1D("theory_rnd","theory_rnd", n_thybins_incl, thyBins_incl);   
  //TH1D *smeared_rnd = new TH1D("smeared_rnd","smeared_rnd", n_smearedbins_incl, smearedBins_incl);   
  TH1D *theory_rnd    = (TH1D*)theory->Clone("theory_rnd");
  theory_rnd->Reset("M ICES");
  TH1D *smeared_rnd    = (TH1D*)theory->Clone("smeared_rnd");
  smeared_rnd->Reset("M ICES");
  
  RooUnfoldResponse response(smeared_rnd,theory_rnd); 
  
  makeToySpectra(theory_spl3clone, spline3, fJER_ynew, 
  		 nEvents, theory_rnd, smeared_rnd, &response);  
  
  divideBinWidth(theory_rnd);
  theory_rnd->SetLineColor(kRed);
  
  divideBinWidth(smeared_rnd);
  smeared_rnd->SetLineColor(kBlue);    
  
  Float_t maxy=theory_rnd->GetMaximum();    cout<<"maxy="<<maxy<<endl;
  maxy*=2.;  cout<<"maxy="<<maxy<<endl;
  Float_t miny=smeared_rnd->GetMinimum();   cout<<"miny="<<miny<<endl;
  miny/=2.; cout<<"miny="<<miny<<endl;
  if(miny<0.){
    cout<<"WARNING! Bin contents of smeared hist is NEGATIVE! Not possible!"<<endl;
    //miny=fabs(miny);
    miny=10.;
  }
  theory_rnd->SetAxisRange(miny,maxy,"Y");
  theory_rnd->SetTitle( ((std::string)("Toy MC Incl. NLO #sigma;;"+ddxsec_yax)).c_str() );
  
  cout<<"creating ratio hist of true/smeared"<<endl;
  TH1D *true_smeared_rat=(TH1D*)theory_rnd->Clone("true_smeared_rat"); 
  if(printBaseDebug)true_smeared_rat->Print("base");
  true_smeared_rat->SetTitle("; Jet p_T (GeV); True/Smeared          ");
  true_smeared_rat->Divide(true_smeared_rat,smeared_rnd,1.,1.,"B");   
  
  ///// check plots True/Smeared 
  std::string ratcanvtitle="plot_true_smeared_rat_y"+etaBin_str;
  TCanvas *plot_true_smeared_rat = new TCanvas(ratcanvtitle.c_str(),ratcanvtitle.c_str(),1400,1200);  
  plot_true_smeared_rat->cd();
  
  TPad* overlay_pad=new TPad("overlay_pad","pad for spectra",0.0,0.3,1.0,1.0);
  overlay_pad->SetLogx(1);
  overlay_pad->SetLogy(1);
  overlay_pad->SetBottomMargin(0);
  overlay_pad->Draw();//tells the canvas theres a pad at the coordinates specified in the new constructor
  overlay_pad->cd();
  
  theory_rnd->Draw("HIST E");  
  smeared_rnd->Draw("HIST E SAME");
  
  TLegend* leg_smear=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
  leg_smear->AddEntry(theory_rnd  , "Truth Incl. NLO #sigma" , "lp" );
  leg_smear->AddEntry(smeared_rnd , "Smeared Incl. NLO #sigma" , "lp" );
  leg_smear->Draw();
  
  plot_true_smeared_rat->cd();//go back to main canvas before doing new pad
  
  TPad* ratio_pad=new TPad("ratio_pad","pad for ratio",0.0,0.05,1.0,0.30);
  ratio_pad->SetLogx(1);
  ratio_pad->SetLogy(0);
  ratio_pad->SetTopMargin(0);
  //ratio_pad->SetBottomMargin(0.3);
  ratio_pad->Draw();//tell canvas there's a pad at the specified coordinates in the new constructor
  ratio_pad->cd();
  
  true_smeared_rat->GetXaxis()->SetMoreLogLabels(true);
  true_smeared_rat->GetXaxis()->SetNoExponent(true);
  true_smeared_rat->SetAxisRange(0.6,1.4,"Y");
  true_smeared_rat->Draw("HIST E");        
  
  TLine* lineatone=new TLine(ptmin_thy, 1., ptmax_thy, 1.);
  lineatone->Draw();
  
  TH2* response_th2=(response.Hresponse());
  std::string respcanvtitle="plot_response_th2_y"+etaBin_str;
  TCanvas* plot_response_th2=new TCanvas(respcanvtitle.c_str(),respcanvtitle.c_str(),1200, 1000);
  plot_response_th2->cd()->SetLogx(1);
  plot_response_th2->cd()->SetLogy(1);
  plot_response_th2->cd()->SetLogz(1);
  plot_response_th2->cd();
  
  response_th2->SetTitle(("Toy MC Response Matrix from JER Smeared NLO #sigma_{jet} #||{y} bin # "+etaBin_str).c_str());
  response_th2->GetXaxis()->SetTitle("RECO Jet p_{T} [GeV]");
  response_th2->GetYaxis()->SetTitle("GEN Jet p_{T} [GeV]");
  response_th2->Draw("COLZ");
  //////////////////////  END production of Smeared NLO spectra     //////////////////////    
  
  
  
  
  
  
  //////////////////////  START production of Smeared NP+NLO spectra     //////////////////////  
  cout<<"creating TH1 for toy NP+NLO spectra generation"<<endl<<endl;
  TH1D *theory_rnd_NP    = (TH1D*)theory->Clone(("theory_rnd_NP_y"+etaBin_str).c_str());
  theory_rnd_NP->Reset("M ICES");
  TH1D *smeared_rnd_NP    = (TH1D*)theory->Clone(("smeared_rnd_NP"+etaBin_str).c_str());
  smeared_rnd_NP->Reset("M ICES");
  
  cout<<"creating instance of RooUnfoldResponse class (response_NP)"<<endl;
  RooUnfoldResponse response_NP(smeared_rnd_NP,theory_rnd_NP); 
  
  makeToySpectra(theory_NP_spl3clone, spline3_NP, fJER_ynew, 
  		 nEvents, theory_rnd_NP, smeared_rnd_NP, &response_NP);  
  
  divideBinWidth(theory_rnd_NP);
  theory_rnd_NP->SetLineColor(kRed);
  
  divideBinWidth(smeared_rnd_NP);
  smeared_rnd_NP->SetLineColor(kBlue);    
  
  Float_t NPmaxy=theory_rnd_NP->GetMaximum();    cout<<"NPmaxy="<<NPmaxy<<endl;
  NPmaxy*=2.;  cout<<"NPmaxy="<<NPmaxy<<endl;
  Float_t NPminy=smeared_rnd_NP->GetMinimum();   cout<<"NPminy="<<NPminy<<endl;
  NPminy/=2.; cout<<"NPminy="<<NPminy<<endl;
  if(NPminy<0.){
    cout<<"WARNING! Bin contents of smeared hist is NEGATIVE! Not possible!"<<endl;
    //NPminy=fabs(NPminy);
    NPminy=10.;
  }
  theory_rnd_NP->SetAxisRange(NPminy,NPmaxy,"Y");
  theory_rnd_NP->SetTitle( ((std::string)("Toy MC Incl. NP+NLO #sigma;;"+ddxsec_yax)).c_str() );
  
  cout<<"creating ratio hist of true/smeared"<<endl;
  TH1D *NP_true_smeared_rat=(TH1D*)theory_rnd_NP->Clone("NP_true_smeared_rat"); 
  if(printBaseDebug)NP_true_smeared_rat->Print("base");
  NP_true_smeared_rat->SetTitle("; Jet p_T (GeV); True/Smeared          ");
  NP_true_smeared_rat->Divide(NP_true_smeared_rat,smeared_rnd_NP,1.,1.,"B");   
  
  ///// check plots True/Smeared 
  std::string NPratcanvtitle="plot_NP_true_smeared_rat_y"+etaBin_str;
  TCanvas *plot_NP_true_smeared_rat = new TCanvas(NPratcanvtitle.c_str(),NPratcanvtitle.c_str(),1400,1200);  
  plot_NP_true_smeared_rat->cd();
  
  TPad* NPoverlay_pad=new TPad("NPoverlay_pad","pad for spectra",0.0,0.3,1.0,1.0);
  NPoverlay_pad->SetLogx(1);
  NPoverlay_pad->SetLogy(1);
  NPoverlay_pad->SetBottomMargin(0);
  NPoverlay_pad->Draw();//tells the canvas theres a pad at the coordinates specified in the new constructor
  NPoverlay_pad->cd();
  
  theory_rnd_NP->Draw("HIST E");  
  smeared_rnd_NP->Draw("HIST E SAME");

  TLegend* NPleg_smear=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
  NPleg_smear->AddEntry(theory_rnd_NP  , "Truth Incl. NP+NLO #sigma" , "lp" );
  NPleg_smear->AddEntry(smeared_rnd_NP , "Smeared Incl. NP+NLO #sigma" , "lp" );
  NPleg_smear->Draw();
  
  plot_NP_true_smeared_rat->cd();//go back to main canvas before doing new pad
  
  TPad* NPratio_pad=new TPad("NPratio_pad","pad for ratio",0.0,0.05,1.0,0.30);
  NPratio_pad->SetLogx(1);
  NPratio_pad->SetLogy(0);
  NPratio_pad->SetTopMargin(0);
  //NPratio_pad->SetBottomMargin(0.3);
  NPratio_pad->Draw();//tell canvas there's a pad at the specified coordinates in the new constructor
  NPratio_pad->cd();
  
  NP_true_smeared_rat->GetXaxis()->SetMoreLogLabels(true);
  NP_true_smeared_rat->GetXaxis()->SetNoExponent(true);
  NP_true_smeared_rat->SetAxisRange(0.6,1.4,"Y");
  NP_true_smeared_rat->Draw("HIST E");        
  
  //TLine* lineatone=new TLine(thyBins_incl[0],1.,thyBins_incl[n_thybins_incl],1.);
  lineatone->Draw();
  
  
  TH2* response_NP_th2=(response_NP.Hresponse());
  std::string NPrespcanvtitle="plot_response_NP_th2_y"+etaBin_str;
  TCanvas* plot_response_NP_th2=new TCanvas(NPrespcanvtitle.c_str(),NPrespcanvtitle.c_str(),1200, 1000);
  plot_response_NP_th2->cd()->SetLogx(1);
  plot_response_NP_th2->cd()->SetLogy(1);
  plot_response_NP_th2->cd()->SetLogz(1);
  plot_response_NP_th2->cd();
  
  response_NP_th2->SetTitle("Toy MC Response Matrix from JER Smeared NP+NLO #sigma_{jet}");
  response_NP_th2->GetXaxis()->SetTitle("RECO Jet p_{T} [GeV]");
  response_NP_th2->GetYaxis()->SetTitle("GEN Jet p_{T} [GeV]");
  response_NP_th2->Draw("COLZ");
  //////////////////////  END production of Smeared NLO spectra     //////////////////////    

  




  


  
  
  ////// Write out histos
  outf->cd();
  
  // NP FITS 
  fNP->Write();

  // NP FIT CANV
  plot_NP->Write();  
  
  //JET ENERGY RESOLUTION
  fJER_ynew->Write();//FIT
  plotJER->Write();//CANV
  
  // NLO HISTS
  theory->Write();
  plot_NLOxsec->Write();
  
  // NP+NLO HISTS
  theory_NP->Write();    
  plot_NPNLOxsec->Write();
  
  // SPLINES
  spline3->Write();
  spline3_NP->Write();
  plot_splines->Write();
  
  
  
  //TOY MC, JER SMEARED NLO SPECTRA
  theory_rnd->Write();
  smeared_rnd->Write(); 
  true_smeared_rat->Write();  //ratio
  plot_true_smeared_rat->Write();

  response.Write("response");//response matrix
  response_th2->Write();
  plot_response_th2->Write();  

  
  ////TOY MC NP+NLO SPECTRA, SMEARED+TRUE   
  theory_rnd_NP->Write();
  smeared_rnd_NP->Write(); 
  NP_true_smeared_rat->Write();  //ratio
  plot_NP_true_smeared_rat->Write();

  response_NP.Write("response_NP");
  response_NP_th2->Write();
  plot_response_NP_th2->Write();  
  

  
  outf->Write();
  return;
}


//  steering ---------------------------------------------------------------------------------
int main(int argc, char* argv[]){  
  int rStatus = -1;   
  
  if(argc==3) {
    smearTheorySpectra_gaussCoreJER_singleYbin( (string) argv[1] , (int) std::atoi(argv[2]) );  
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
