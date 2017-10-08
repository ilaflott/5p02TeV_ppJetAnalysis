#include "smearTheorySpectra.h"



void smearTheorySpectra_gaussCoreJER( string inputString ){
  
  
  //gStyle->SetOptStat(0);
  //gROOT->ForceStyle();
  bool printBaseDebug=true;
  
  // Output
  string outputFile= inputString+"_gaussSmear_00eta20.root";
  cout<<"opening output:"<<outputFile<<endl<<endl;
  TFile *outf    = new TFile(outputFile.c_str(), "RECREATE" );
  
  // Input, Theory spectrum
  string inputFile = "fNLOJetsSpectra/"+inputString+".root";
  // TFile *f1 = new TFile("PredictionsFastNLO-ak7-CT14nlo.root");
  cout<<"opening input:"<<inputFile<<endl<<endl;
  TFile *f1 = new TFile(inputFile.c_str());
  
  
  TH1D *theory_y0 = (TH1D*)f1->Get("h0100100");
  if(printBaseDebug)  theory_y0->Print("Base");
  TH1D *theory_y1 = (TH1D*)f1->Get("h0100200");
  if(printBaseDebug)  theory_y1->Print("Base");
  TH1D *theory_y2 = (TH1D*)f1->Get("h0100300");
  if(printBaseDebug)  theory_y2->Print("Base");
  TH1D *theory_y3 = (TH1D*)f1->Get("h0100400");
  if(printBaseDebug)  theory_y3->Print("Base");
  cout<<endl;
  
  //  TH1D *theory_y4 = (TH1D*)f1->Get("InclusiveJet_5bin;1");
  //  TH1D *theory_y5 = (TH1D*)f1->Get("InclusiveJet_6bin;1");
  //  TH1D *theory_y6 = (TH1D*)f1->Get("InclusiveJet_7bin;1");
  
  cout <<"making 00eta20 hist from theory curves..."<<endl<<endl;
  TH1D *theory_ynew = make00eta20Hist( (TH1D*) theory_y0, 
				       (TH1D*) theory_y1, 
				       (TH1D*) theory_y2, 
				       (TH1D*) theory_y3 );  
  cout<<"done making hist."<<endl;
  if(printBaseDebug)theory_ynew->Print("base");
  
  cout<<endl;
  //<assert(false);
  
  /////////////// NPs to theory spectra ////////////////////////////////////////////////////////
  /////////////// Create functions just to check NPs
  
  TF1 *fNP_ynew = new TF1("fNP_ynew",NP_ynew_str.c_str(), thyBins_incl[0], thyBins_incl[n_thybins_incl]);
  
  /// check NPs
  TCanvas *plotNPs = new TCanvas("plotNPs", "plotNPs",1200,800);
  
  //  plotNPs->Divide(4,2);  
  plotNPs->cd();
  fNP_ynew->Draw();
  plotNPs->Update();
  
  
  /////////////// Core p_T Resolution ////////////////////////////////////////////////////////
  /////////////// Create functions for checking and using them during the smearing
  TF1 *fres_ynew = new TF1("fres_ynew",gres_ynew_str.c_str(),thyBins_incl[0], thyBins_incl[n_thybins_incl]);
  
  fres_ynew->SetLineColor(2);
  fres_ynew->SetMinimum(0); fres_ynew->SetMaximum(0.25);



  /// check Core p_T Resolution
  TCanvas *plotres = new TCanvas("plotres", "plotres",1200,800);
  plotres->cd();
  fres_ynew->Draw();

  
  TLegend *leg=new TLegend(0.53,0.50,0.85,0.85);
  leg->AddEntry(fres_ynew,"PYTHIA8 CUETP8M1 @ 5.02 TeV", "");
  leg->AddEntry(fres_ynew,"ak4PFJets", "");
  leg->AddEntry(fres_ynew,"0<|y|<2.0", "l"); 
  
  leg->SetTextFont(42);
  leg->SetFillColor(kWhite);
  leg->Draw();
  
  plotres->Update();
  
  
  
  
  //////////// Apply NPs to theory ////////////////////////////////////////////////////////
  /// Cross section ynew
  cout<<"making theory cross section hist"<<endl;
  TH1D *Xsection_theory_ynew=(TH1D*)theory_ynew->Clone("xsec_theory_ynew");
  if(printBaseDebug)Xsection_theory_ynew->Print("base");
  for(int i=0; i<theory_ynew->GetNbinsX(); ++i){
    double bin_value  = theory_ynew->GetBinContent(1+i);
    double bin_center = theory_ynew->GetBinCenter(1+i);
    double new_bin_value = bin_value * NP_ynew(bin_center);
    Xsection_theory_ynew->SetBinContent(1+i,new_bin_value);
  }
  cout<<"done making theory cross section hist"<<endl;
  if(printBaseDebug)Xsection_theory_ynew->Print("base");
  
  /////////////// Fitting Spectra with spline3
  /////////////// Not worning well for y6 bin because Klaus's spectra ends < 500 GeV
  /////////////// We must fit spectra by hand for y6 bin
  
  ///Create Cubic Splines using Cross sections
  cout<<"creating TSpline of cross section hist"<<endl;
  TSpline3 *spline3_ynew = new TSpline3(Xsection_theory_ynew);
  
  spline3_ynew->SetLineColor(6);
  
  Xsection_theory_ynew->SetAxisRange(thyBins_incl[0],thyBins_incl[n_thybins_incl],"X"); 
  
  /// check spectra with splines
  TCanvas *TheorySpectra = new TCanvas("TheorySpectra", "Theory Spectra",1200,800);
  TheorySpectra->cd(0);
  gPad->SetLogy();

  cout<<"drawing Tspline and cross section on canvas"<<endl;
  Xsection_theory_ynew->Draw();
  spline3_ynew->Draw("same");
  
  
  
  //////////////////////  Start production of Smeared spectra
  
  TRandom3 *rnd = new TRandom3();
  int nEvents=10000000;  /// Increase your atistics here typical 100Mevents are enought
  
  cout<<"creating new TH1D's"<<endl<<endl;
  TH1D *True_incl_jet_ynew    = new TH1D("True_incl_jet_ynew","True_incl_jet_ynew", n_thybins_incl, thyBins_incl);   
  True_incl_jet_ynew->Sumw2();
  if(printBaseDebug)True_incl_jet_ynew->Print("base");
  
  //TH1D *Smeared_incl_jet_ynew = new TH1D("Smeared_incl_jet_ynew","Smeared_incl_jet_ynew", n_thybins_incl, thyBins_incl);   
  TH1D *Smeared_incl_jet_ynew = new TH1D("Smeared_incl_jet_ynew","Smeared_incl_jet_ynew", n_smearedbins_incl, smearedBins_incl);   
  Smeared_incl_jet_ynew->Sumw2();
  if(printBaseDebug)Smeared_incl_jet_ynew->Print("base");
  
  TH1D *NormBinWidth_True_incl_jet_ynew    = new TH1D("NormBinWidth_True_incl_jet_ynew","NormBinWidth_True_incl_jet_ynew", n_thybins_incl, thyBins_incl); 
  NormBinWidth_True_incl_jet_ynew->Sumw2();
  if(printBaseDebug)NormBinWidth_True_incl_jet_ynew->Print("base");
  
  TH1D *NormBinWidth_Smeared_incl_jet_ynew = new TH1D("NormBinWidth_Smeared_incl_jet_ynew","NormBinWidth_Smeared_incl_jet_ynew", n_smearedbins_incl, smearedBins_incl);   
  NormBinWidth_Smeared_incl_jet_ynew->Sumw2();
  if(printBaseDebug)NormBinWidth_Smeared_incl_jet_ynew->Print("base");
  
  cout<<"constructing RooUnfoldResponse "<<endl;
  RooUnfoldResponse response_ynew(Smeared_incl_jet_ynew,True_incl_jet_ynew);
  
  
  cout<<"Evaluating ynew"<<endl;  
  for(int i=0;i<nEvents;++i){  
    double ptmin_thy = thyBins_incl[0];
    double ptmax_thy = thyBins_incl[n_thybins_incl];

    double ptmin_smeared = smearedBins_incl[0];
    double ptmax_smeared = smearedBins_incl[n_smearedbins_incl];
    
    double ptTrue  = rnd->Uniform(ptmin_thy,ptmax_thy);
    double sigma   = gres_ynew(ptTrue);
    
    //sigma=sigma*1.079; //JER Scaling factor 8 TeV  Acoording to Mikko no scaling for the moment
    
    double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
    double pt_w      =  spline3_ynew->Eval(ptTrue);
    
    True_incl_jet_ynew    -> Fill(ptTrue,pt_w);
    Smeared_incl_jet_ynew -> Fill(ptSmeared,pt_w);
    
    //printf("ptTrue=%f  ptSmeared=%f sigma=%f  ptTrue/ptSmeared=%f  pt_w=%f\n", ptTrue, ptSmeared, sigma, ptTrue/ptSmeared ,pt_w);
    
    if (((ptSmeared>=ptmin_smeared)&&(ptSmeared<=ptmax_smeared))&&
	((ptTrue>=ptmin_thy)&&(ptTrue<=ptmax_thy))){
      response_ynew.Fill(ptSmeared,ptTrue,pt_w);    }
    //else if ((ptTrue>=ptmin)&&(ptTrue<=ptmax)){
    //  response_ynew.Miss(ptTrue,pt_w);    }
    //else if ((ptSmeared>=ptmin)&&(ptSmeared<=ptmax)){
    //  response_ynew.Fake(ptSmeared,pt_w);    }
  }
  cout<<"done evaluating ynew"<<endl;
  
  NormBinWidth_True_incl_jet_ynew    = cross_section(True_incl_jet_ynew);
  NormBinWidth_True_incl_jet_ynew->SetLineColor(2);
  if(printBaseDebug)NormBinWidth_True_incl_jet_ynew->Print("base");
  
  NormBinWidth_Smeared_incl_jet_ynew = cross_section(Smeared_incl_jet_ynew);
  NormBinWidth_Smeared_incl_jet_ynew->SetLineColor(2);
  if(printBaseDebug)NormBinWidth_Smeared_incl_jet_ynew->Print("base");
  

  
  cout<<"creating ratio hist of true/smeared"<<endl;
  TH1D *rat_ynew_true_smeared=(TH1D*)True_incl_jet_ynew->Clone("True_incl_jet_ynew_clone_ratio"); 
  if(printBaseDebug)rat_ynew_true_smeared->Print("base");
  rat_ynew_true_smeared->SetTitle("ynew True/Smeared;p_T (GeV); ynew True/Smeared          ");
  
  //rebin 4 ratio only
  cout<<"rebinning..."<<endl;
  rat_ynew_true_smeared=(TH1D*)rat_ynew_true_smeared->TH1::Rebin(n_smearedbins_incl, "True_incl_jet_ynew_clone_ratio_rebin", smearedBins_incl)    ;
  if(printBaseDebug)rat_ynew_true_smeared->Print("base");

  cout<<"dividing..."<<endl;
  //rat_ynew_true_smeared->Divide(True_incl_jet_ynew,Smeared_incl_jet_ynew,1.,1.,"B");  
  //rat_ynew_true_smeared->Divide(Smeared_incl_jet_ynew);
  rat_ynew_true_smeared->Divide(rat_ynew_true_smeared,Smeared_incl_jet_ynew,1.,1.,"B");  
  


  //rat_ynew_true_smeared->SetMaximum(1.1); 
  rat_ynew_true_smeared->SetAxisRange(smearedBins_incl[0],smearedBins_incl[n_smearedbins_incl],"X"); 
  
  ///// check plots True/Smeared 
  TCanvas *trueOverSmeared = new TCanvas("trueOverSmeared", "trueOverSmeared",1600,800);
  
  trueOverSmeared->cd();
  gPad->SetGrid(); 
  rat_ynew_true_smeared->Draw();     
  
  trueOverSmeared->Update();
  
  
  
  ////// Write out histos
  outf->cd();
  
  //hists opened
  theory_y0->Write();
  theory_y1->Write();
  theory_y2->Write();
  theory_y3->Write();

  theory_ynew->Write();//made from the 4 above so it goes in the "opened" section
  
  //fits used
  fNP_ynew->Write();
  fres_ynew->Write();
  
  
  //hists made

  Xsection_theory_ynew->Write();  

  True_incl_jet_ynew->Write();
  NormBinWidth_True_incl_jet_ynew->Write("NormBinWidth_True_incl_jet_ynew");      

  Smeared_incl_jet_ynew->Write(); 
  NormBinWidth_Smeared_incl_jet_ynew->Write("NormBinWidth_Smeared_incl_jet_ynew"); 

  response_ynew.Write("response_ynew");

  rat_ynew_true_smeared->Write();

  //Cavasses created
  plotNPs->Write();
  plotres->Write();
  TheorySpectra->Write();
  trueOverSmeared->Write();
  //splines created
  spline3_ynew->Write();
  
  outf->Write();
  return;
  
  
  
}


//  steering ---------------------------------------------------------------------------------
int main(int argc, char* argv[]){  
  int rStatus = -1;   
  
  //if(argc==1) smearTheorySpectra_gaussCoreJER();
  if(argc==2) {
    smearTheorySpectra_gaussCoreJER( (string) argv[1] );  
    rStatus=0;
  }
  else {
    rStatus=1;
  }
  
  cout<<"rStatus="<<rStatus<<endl;
  return rStatus;
}
