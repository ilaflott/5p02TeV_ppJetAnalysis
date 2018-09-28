#include <cstdlib>
#include <iostream>
#include <cassert>

#include "TH2.h"
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"

#include "unf_hstrings.h"

const int CANVXW=1800;
const int CANVYW=(int)
  (
   ((float)CANVXW) * 0.70
   );
//hist-name array pattern is always
// // std::string namearray={ MCprior, MCmeas, datameas, dataunfolded}

void checkHistVector(std::vector<TH1D*> th1vect){
  //if(!th1vect){
  //  std::cout<<"*** WARNING *** vector pointer DNE or is NULL *** WARNING ***"<<std::endl;
  //  std::cout<<"*** WARNING *** i may crash!                  *** WARNING ***"<<std::endl; 
  //  std::cout<<"*** WARNING *** vector pointer DNE or is NULL *** WARNING ***"<<std::endl;
  //  std::cout<<"*** WARNING *** i may crash!                  *** WARNING ***"<<std::endl;
  //  std::cout<<"*** WARNING *** vector pointer DNE or is NULL *** WARNING ***"<<std::endl;
  //  std::cout<<"*** WARNING *** i may crash!                  *** WARNING ***"<<std::endl;
  //}
  for (int i=0;i<(int)th1vect.size();i++){
    TH1D* theHist=(TH1D*) th1vect[i];    
    if(!theHist){
      std::cout<<"*** WARNING *** theHist pointer DNE or is NULL *** WARNING ***"<<std::endl;
      std::cout<<"*** WARNING *** i="<<i<<" may crash!           *** WARNING ***"<<std::endl;      
    }
  }
  return;
}

// general stuff --------------
void writeNsaveCanvas(TFile* fout, TCanvas* canv){
  fout->cd();
  canv->Write();
  canv->SaveAs( ( ((std::string)canv->GetName())+".png").c_str());
  return;
}
// general draw stuff --------------
void sethArrayTitles(std::vector<TH1D*>harray, std::string basetitle, std::vector<std::string> htitles){
  for(int i=0; i<((int)harray.size()); i++){
    std::string htitlefinal=basetitle+" "+htitles[i];
    harray[i]->SetTitle( htitlefinal.c_str());
  }
  return;
}


// ratio draw stuff --------------
void setupRatioCanv(TCanvas* canv){
  canv->SetLogx(1);
  canv->SetLogy(0);  
  return;
}

TLegend* buildRatioLegend(std::vector<TH1D*> harray){
  TLegend* theleg= new TLegend(0.1,0.7,0.5,0.9);
  for(int i=0; i<((int)harray.size()); i++)
    theleg->AddEntry(harray[i], NULL, "lpe");
  return theleg;
}

void drawhArrayOnRatioCanv(TCanvas* canv, std::vector<TH1D*> harray, std::string plot_title, std::string yaxtitle=""){
  canv->cd();
  
  float xlo=harray[0]->GetBinLowEdge(1);
  float xhi=harray[0]->GetBinLowEdge(harray[0]->GetNbinsX()) + harray[0]->GetBinWidth(harray[0]->GetNbinsX());
  
  TLine* thelineatone=new TLine(xlo,1.,xhi,1.);
  thelineatone->SetLineStyle(5);
  TLine* thelineat1p1=(TLine*)thelineatone->Clone();
  thelineat1p1->SetY1(1.1);  thelineat1p1->SetY2(1.1);
  TLine* thelineat0p9=(TLine*)thelineatone->Clone();
  thelineat0p9->SetY1(0.9);  thelineat0p9->SetY2(0.9);
  
  TH1D* emptyh=(TH1D*)harray[0]->Clone("4drawonly");
  emptyh->Reset("M ICES");
  emptyh->SetAxisRange(0.6,1.4,"Y");
  emptyh->GetXaxis()->SetMoreLogLabels(true);
  emptyh->GetXaxis()->SetNoExponent(true);
  emptyh->GetXaxis()->SetTitle("Jet p_{T} (GeV)");  
  if(yaxtitle.length()>0)
    emptyh->GetYaxis()->SetTitle(yaxtitle.c_str());
  else
    emptyh->GetYaxis()->SetTitle("Spectra Ratios");      
  emptyh->SetTitle(plot_title.c_str());
  emptyh->DrawClone("HIST E");    
  
  thelineatone->Draw();
  thelineat1p1->Draw();
  thelineat0p9->Draw();
  
  for(int i =0; i<((int)harray.size()); i++)
    harray[i]->DrawClone("HIST E SAME");
  
  return;
}

std::vector<TH1D*> makeCommonhRatios(std::vector<TH1D*> commonh, std::string prior_name){
  TH1D* MCmeas=(TH1D*)commonh[0]  ;//->Clone("");
  TH1D* MCprior=(TH1D*)commonh[1] ;//->Clone("");
  TH1D* datameas=(TH1D*)commonh[2];//->Clone("");
  
  TH1D* hrat_MCprior_MCmeas=(TH1D*)MCprior->Clone("");
  hrat_MCprior_MCmeas->Divide(MCmeas);
  hrat_MCprior_MCmeas->SetTitle((prior_name+".MC Prior/"+prior_name+".MC Meas.").c_str());
  
  TH1D* hrat_MCprior_datameas=(TH1D*)MCprior->Clone("");
  hrat_MCprior_datameas->SetLineColor(datameas->GetLineColor());
  hrat_MCprior_datameas->SetMarkerStyle(datameas->GetMarkerStyle());
  hrat_MCprior_datameas->SetMarkerColor(datameas->GetMarkerColor());
  hrat_MCprior_datameas->Divide(datameas);
  hrat_MCprior_datameas->SetTitle((prior_name+".MC Prior/Data Meas.").c_str());
  
  TH1D* hrat_MCmeas_datameas=(TH1D*)MCmeas->Clone("");
  hrat_MCmeas_datameas->Divide(datameas);
  hrat_MCmeas_datameas->SetTitle((prior_name+".MC Meas/Data Meas.").c_str());
  
  std::vector<TH1D*>commonhrat_vector;
  commonhrat_vector.push_back(hrat_MCprior_MCmeas);
  commonhrat_vector.push_back(hrat_MCprior_datameas);
  commonhrat_vector.push_back(hrat_MCmeas_datameas);
  return commonhrat_vector;
}



// spectra draw stuff --------------
void setupSpectraCanv(TCanvas* canv){
  canv->SetLogx(1);
  canv->SetLogy(1);  
  return;
}

TLegend* buildSpectraLegend(std::vector<TH1D*> harray){
  TLegend* theleg= new TLegend(0.7,0.7,0.9,0.9);
  for(int i=0; i<((int)harray.size()); i++)
    theleg->AddEntry(harray[i], NULL, "lpe");
  return theleg;
}

void drawhArrayOnSpectraCanv(TCanvas* canv, std::vector<TH1D*> harray, std::string plot_title){
  canv->cd();
  TH1D* emptyh=(TH1D*)harray[0]->Clone("4drawonly");
  emptyh->Reset("ICES");
  emptyh->SetAxisRange((harray[0]->GetMinimum())/2.,(harray[0]->GetMaximum()*2.),"Y");
  emptyh->GetXaxis()->SetMoreLogLabels(true);
  emptyh->GetXaxis()->SetNoExponent(true);
  emptyh->GetXaxis()->SetTitle("Jet p_{T} (GeV)");
  emptyh->GetYaxis()->SetTitle(yaxtitle_xsec.c_str());      
  emptyh->SetTitle(plot_title.c_str());
  emptyh->DrawClone("HIST E");    

  for(int i =0; i<((int)harray.size()); i++)    
    harray[i]->DrawClone("HIST E SAME");
  
  return;
}











void print_compareMethods( TFile* fin_meth1, std::vector<std::string> hnames_meth1,  std::string meth1_name,
			   TFile* fin_meth2, std::vector<std::string> hnames_meth2,  std::string meth2_name,
			   std::string prior_nameA, //for file I/O
			   std::string prior_nameB, //for hist titling
			   //std::string methcomp_title   , 
			   std::vector<std::string> names){//, TFile* fout){
  
  if(!fin_meth1 || !fin_meth2){
    std::cout<<"file not found. skipping comparison between "<<meth1_name<<" and "<<meth2_name<<std::endl;
    return;
  }
  
  std::string fout_name=prior_nameA+"_"+meth1_name+".v."+meth2_name+"_"+"compareMethods.root";  
  TFile* fout=new TFile(fout_name.c_str(),"RECREATE");
  
  
  std::vector<TH1D*> meth1_hvector;
  for(int i =0; i<((int)hnames_meth1.size()); i++)
    meth1_hvector.push_back((TH1D*)fin_meth1->Get(  ((std::string)hnames_meth1[i]).c_str()));
  for(int i =0; i<((int)hnames_meth1.size()); i++)
    std::cout<<meth1_hvector[i]->GetName()<<std::endl;
  
  std::cout<<"drawing method 1's spectra"<<std::endl;
  TCanvas* meth1_spectra_canv=new TCanvas((prior_nameA+"."+meth1_name+"_spectra_canv").c_str(),
					  (prior_nameA+"."+meth1_name+"_spectra_canv").c_str(),CANVXW,CANVYW);
  setupSpectraCanv(meth1_spectra_canv);
  drawhArrayOnSpectraCanv(meth1_spectra_canv, meth1_hvector, (prior_nameB+", "+meth1_name+", Spectra"));
  sethArrayTitles(meth1_hvector,prior_nameB+"."+meth1_name, names);
  TLegend* meth1_spectra_leg=buildSpectraLegend(meth1_hvector);
  meth1_spectra_canv->cd();
  meth1_spectra_leg->Draw();
  writeNsaveCanvas(fout,meth1_spectra_canv);
  //fout->cd();  //meth1_spectra_canv->Write();
  
  
  std::vector<TH1D*> meth2_hvector;
  for(int i =0; i<((int)hnames_meth2.size()); i++)
    meth2_hvector.push_back((TH1D*)fin_meth2->Get(  ((std::string)hnames_meth2[i]).c_str()));
  
  std::cout<<"drawing method 2's spectra"<<std::endl;
  TCanvas* meth2_spectra_canv=new TCanvas((prior_nameA+"."+meth2_name+"_spectra_canv").c_str(),
				     (prior_nameA+"."+meth2_name+"_spectra_canv").c_str(),CANVXW,CANVYW);
  setupSpectraCanv(meth2_spectra_canv);
  drawhArrayOnSpectraCanv(meth2_spectra_canv, meth2_hvector, (prior_nameB+", "+meth2_name+", Spectra"));
  sethArrayTitles(meth2_hvector,prior_nameB+"."+meth2_name, names);

  TLegend* meth2_spectra_leg=buildSpectraLegend(meth2_hvector);
  meth2_spectra_canv->cd();
  meth2_spectra_leg->Draw();
  writeNsaveCanvas(fout,meth2_spectra_canv);
  //fout->cd();  //meth2_spectra_canv->Write();
  

  std::vector<TH1D*> common_hvector;
  common_hvector.push_back((TH1D*)fin_meth1->Get(  ((std::string)hnames_meth1[0]).c_str()));  
  common_hvector.push_back((TH1D*)fin_meth1->Get(  ((std::string)hnames_meth1[1]).c_str()));
  common_hvector.push_back((TH1D*)fin_meth1->Get(  ((std::string)hnames_meth1[2]).c_str()));
  sethArrayTitles(common_hvector,prior_nameB,names);

  std::cout<<"drawing spectra common to each method (prior hist + meas hists)"<<std::endl;
  TCanvas* common_ratiocanv=new TCanvas((prior_nameA+"."+meth1_name+"_"+meth2_name+"_commonhrats").c_str(),
					(prior_nameA+"."+meth1_name+"_"+meth2_name+"_commonhrats").c_str(),CANVXW,CANVYW);
  common_ratiocanv->cd();
  setupRatioCanv(common_ratiocanv);
  
  std::vector<TH1D*> common_hratvector=makeCommonhRatios(common_hvector,prior_nameB);
  TLegend* commonhrat_leg=buildRatioLegend(common_hratvector);
  drawhArrayOnRatioCanv(common_ratiocanv, common_hratvector, (prior_nameB+", Ratios, Hists common between "+meth1_name+" "+meth2_name));
  common_ratiocanv->cd();
  commonhrat_leg->Draw();
  writeNsaveCanvas(fout,common_ratiocanv);
  //fout->cd();  //common_ratiocanv->Write();
  
  
  
  








  
  std::vector<TH1D*> hvector_unf;//(common_hvector);
  hvector_unf.push_back(meth1_hvector[(meth1_hvector.size()-1)]); //bayes 
  hvector_unf.push_back(meth2_hvector[(meth2_hvector.size()-1)]); //SVD //meth2 is in the back of the vector
  hvector_unf[hvector_unf.size()-1]->SetLineColor(kOrange+8);
  
  TCanvas* unf_spectracanv=new TCanvas((prior_nameA+"."+meth1_name+"_"+meth2_name+"_unfh").c_str(),
				       (prior_nameA+"."+meth1_name+"_"+meth2_name+"_unfh").c_str(),CANVXW,CANVYW);
  unf_spectracanv->cd();
  setupSpectraCanv(unf_spectracanv);
  drawhArrayOnSpectraCanv(unf_spectracanv, hvector_unf,(prior_nameB+"."+meth1_name+"."+meth2_name+", Unf. Spectra"));
  TLegend* meth12_spectra_leg=buildSpectraLegend(hvector_unf);
  unf_spectracanv->cd();
  meth12_spectra_leg->Draw();
  writeNsaveCanvas(fout,unf_spectracanv);
  //fout->cd();  unf_spectracanv->Write();
    

  TCanvas* unf_ratiocanv=new TCanvas((prior_nameA+"."+meth1_name+"_"+meth2_name+"_unfhrats").c_str(),
				     (prior_nameA+"."+meth1_name+"_"+meth2_name+"_unfhrats").c_str(),CANVXW,CANVYW);
  unf_ratiocanv->cd();
  setupRatioCanv(unf_ratiocanv);

  float xlo=hvector_unf[0]->GetBinLowEdge(1);
  float xhi=hvector_unf[0]->GetBinLowEdge(hvector_unf[0]->GetNbinsX()) + hvector_unf[0]->GetBinWidth(hvector_unf[0]->GetNbinsX());
  TLine* thelineatone=new TLine(xlo,1.,xhi,1.);
  thelineatone->SetLineStyle(5);
  TLine* thelineat1p1=(TLine*)thelineatone->Clone();
  thelineat1p1->SetY1(1.1);  thelineat1p1->SetY2(1.1);
  TLine* thelineat0p9=(TLine*)thelineatone->Clone();
  thelineat0p9->SetY1(0.9);  thelineat0p9->SetY2(0.9);
  
  TH1D* hunf_rat=(TH1D*)hvector_unf[0]->Clone();  //hunf_rat->Reset("MIES");
  hunf_rat->Divide(hvector_unf[1]);
  hunf_rat->SetTitle(((prior_nameB+", "+meth1_name+"/"+meth2_name+",Unf. Data Ratio")).c_str());
  hunf_rat->GetXaxis()->SetMoreLogLabels(true);
  hunf_rat->GetXaxis()->SetNoExponent(true);
  hunf_rat->SetAxisRange(0.6,1.4,"Y");
  hunf_rat->GetYaxis()->SetTitle((meth1_name+"/"+meth2_name).c_str());
  
  unf_ratiocanv->cd();
  hunf_rat->Draw("HIST E");
  thelineatone->Draw();	
  thelineat1p1->Draw();	
  thelineat0p9->Draw();	  
  hunf_rat->Draw("HIST E SAME");
  
  writeNsaveCanvas(fout,unf_ratiocanv);

  ////////////
  // TO DO ?//
  ////////////
  // // compare meth1 unf data w/ meas data
  // // compare meth2 unf data w/ meas data

  // // compare meth1 unf data w/ meas MC
  // // compare meth2 unf data w/ meas MC   

  // // compare meth1 unf data w/ prior  
  // // compare meth2 unf data w/ prior
  ////////////
  // TO DO ?//
  ////////////
  
  //fout->cd();
  //hunf_rat->Write((prior_nameA+"."+meth1_name+"_"+meth2_name+"_unfth1rats").c_str());
  //unf_ratiocanv->Write();
  //unf_ratiocanv->SaveAs((prior_nameA+"."+meth1_name+"_"+meth2_name+"_unfhrats.png").c_str());
  //std::vector<TH1D*> commonwunf_hratvector=makeCommonhRatios_wunf(commonwunf_hvector,prior_name);
  //TLegend* commonhratwunf_leg=buildRatioLegendwunf(commonwunf_hratvector);
  //drawhArrayOnRatioCanv(commonwunf_ratiocanv, commonwunf_hratvector, (prior_name+"."+meth2_name+"."+meth1_name+", Unf. Ratios"));
  //commonwunf_ratiocanv->cd();
  //commonhratwunf_leg->Draw();
  //fout->cd();
  //commonwunf_ratiocanv->Write();
  
  return;
}

//i.e. match the priors, but make sure diff methods
void quick_compareMethods(TFile* fin_p8byes, TFile* fin_NLObyes, TFile* fin_NLOwNPbyes,	 
			  TFile* fin_p8svd,  TFile* fin_NLOsvd,   TFile* fin_NLOwNPsvd ){//, TFile* fout   ) {
  std::cout<<"comparing data unfolded w/ same priors yet diff methods"<<std::endl;  
  
  if( !((bool) fin_p8byes) || !((bool) fin_NLObyes) || !((bool) fin_NLOwNPbyes) ||	 
      !((bool) fin_p8svd ) || !((bool) fin_NLOsvd ) || !((bool) fin_NLOwNPsvd )  ){
    std::cout<<"file not found. exit."<<std::endl;
    assert(false);  }      
  
  //-----------------------------------------//
  //compare py8byes file w/ p8svd file (prior in common is p8 prior)
  std::cout<<"comparing PY8 bayes w/ PY8 SVD"<<std::endl;
  print_compareMethods( (TFile*)fin_p8byes, PY8_bayes_NAMES,  (std::string)"Bayes",
			(TFile*)fin_p8svd,  PY8_SVD_NAMES,    (std::string)"SVD",
			(std::string) "PY8_LO", (std::string) "PY8 LO",
			(std::vector<std::string>)HNAMES );//,
  //(TFile*) fout) ;
  //-----------------------------------------//
  
  //////-----------------------------------------//
  std::cout<<"comparing NLO bayes w/ NLO SVD"<<std::endl;
  print_compareMethods( (TFile*)fin_NLObyes, NLO_bayes_NAMES,  (std::string)"Bayes",
  			(TFile*)fin_NLOsvd,  NLO_SVD_NAMES,    (std::string)"SVD",
  			(std::string) "NNPDF_NNLO", (std::string) "NNPDF NNLO",
  			(std::vector<std::string>)HNAMES); //);,
  //(TFile*) fout) ;
  
  ////-----------------------------------------//
  
  ////-----------------------------------------//
  std::cout<<"comparing NLO+NP bayes w/ NLO+NP SVD"<<std::endl;
  print_compareMethods( (TFile*)fin_NLOwNPbyes, NLOwNP_bayes_NAMES,  (std::string)"Bayes",
			(TFile*)fin_NLOwNPsvd,  NLOwNP_SVD_NAMES,    (std::string)"SVD",
			(std::string) "NNPDF_NNLOwNP",	(std::string) "NNPDF NNLO+NP",
			(std::vector<std::string>)HNAMES);//,
  //(TFile*) fout) ;
  ////-----------------------------------------//
  
  return;
}



//i.e. match the priors, but make sure diff methods
void print_comparePriors(TFile* fin_prior1, std::vector<std::string> hnames_prior1, std::string prior1_name,std::string prior1_desc,
			 TFile* fin_prior2, std::vector<std::string> hnames_prior2, std::string prior2_name,std::string prior2_desc,
			 TFile* fin_prior3, std::vector<std::string> hnames_prior3, std::string prior3_name,std::string prior3_desc,
			 std::string meth_name, std::vector<std::string> names			 ){
  std::cout<<std::endl<<"running print_comparePriors"<<std::endl;
  
  
  if(!fin_prior1 || !fin_prior2 || !fin_prior3){
    std::cout<<"file not found. skipping comparison between priors using "<<meth_name<<" unfolding."<<std::endl;
    return;
  }
  
  //std::string fout_name=prior_nameA+"_"+meth1_name+".v."+meth2_name+"_"+"compareMethods.root";  
  std::string fout_name=meth_name+"_"+prior1_name+"."+prior2_name+"."+prior3_name+"_"+"comparePriors.root";  
  TFile* fout=new TFile(fout_name.c_str(),"RECREATE");
  

  std::vector<TH1D*> prior1_hvector;
  for(int i =0; i<((int)hnames_prior1.size()); i++)
    prior1_hvector.push_back((TH1D*)fin_prior1->Get(  ((std::string)hnames_prior1[i]).c_str()));
  checkHistVector(prior1_hvector);

  std::vector<TH1D*> prior2_hvector;
  for(int i =0; i<((int)hnames_prior2.size()); i++)
    prior2_hvector.push_back((TH1D*)fin_prior2->Get(  ((std::string)hnames_prior2[i]).c_str()));
  checkHistVector(prior2_hvector);

  std::vector<TH1D*> prior3_hvector;
  for(int i =0; i<((int)hnames_prior3.size()); i++)
    prior3_hvector.push_back((TH1D*)fin_prior3->Get(  ((std::string)hnames_prior3[i]).c_str()));
  checkHistVector(prior3_hvector);
  
  



  //PRIORX OVER MEASURED DATA
  TCanvas* priorX_measdata_canv=new TCanvas((meth_name+"_"+prior1_name+"."+prior2_name+"."+prior3_name+"_measdataratprior").c_str(),
					    (meth_name+"_"+prior1_name+"."+prior2_name+"."+prior3_name+"_measdataratprior").c_str(),
					    CANVXW, CANVYW);
  setupRatioCanv(priorX_measdata_canv);
  
  TH1D* h_measdata=(TH1D*)prior1_hvector[2]->Clone();//clone meas data 
  std::vector<TH1D*> prior_measdata_rat_vector;
  prior_measdata_rat_vector.push_back((TH1D*)prior1_hvector[1]->Clone());
  prior_measdata_rat_vector.push_back((TH1D*)prior2_hvector[1]->Clone());
  prior_measdata_rat_vector.push_back((TH1D*)prior3_hvector[1]->Clone());
  
  prior_measdata_rat_vector[0]->Divide(h_measdata);
  prior_measdata_rat_vector[1]->Divide(h_measdata);
  prior_measdata_rat_vector[2]->Divide(h_measdata);
  
  prior_measdata_rat_vector[0]->SetLineColor(kRed);      prior_measdata_rat_vector[0]->SetMarkerColor(kRed);    
  prior_measdata_rat_vector[1]->SetLineColor(kBlue);     prior_measdata_rat_vector[1]->SetMarkerColor(kBlue);   
  prior_measdata_rat_vector[2]->SetLineColor(kMagenta);  prior_measdata_rat_vector[2]->SetMarkerColor(kMagenta);
  
  drawhArrayOnRatioCanv(priorX_measdata_canv,prior_measdata_rat_vector, meth_name+" Prior MC, Ratio w/ Meas. Data", "Prior MC / Meas. Data");  
  
  prior_measdata_rat_vector[0]->SetTitle( (prior1_desc+" MC Prior").c_str() );
  prior_measdata_rat_vector[1]->SetTitle( (prior2_desc+" MC Prior").c_str() );
  prior_measdata_rat_vector[2]->SetTitle( (prior3_desc+" MC Prior").c_str() );
  
  TLegend* priorX_measdata_legend=buildRatioLegend(prior_measdata_rat_vector);
  priorX_measdata_canv->cd();
  priorX_measdata_legend->Draw();
  
  writeNsaveCanvas(fout, priorX_measdata_canv);
  
  
  //PRIORX UNFOLDED OVER MEASURED DATA
  TCanvas* priorXunf_measdata_canv=new TCanvas((meth_name+"_"+prior1_name+"."+prior2_name+"."+prior3_name+"_measdatarat").c_str(),
					       (meth_name+"_"+prior1_name+"."+prior2_name+"."+prior3_name+"_measdatarat").c_str(),
					       CANVXW, CANVYW);
  setupRatioCanv(priorXunf_measdata_canv);
  
  //TH1D* h_measdata=(TH1D*)prior1_hvector[2]->Clone();//clone meas data 
  std::vector<TH1D*> priorunf_measdata_rat_vector;



  priorunf_measdata_rat_vector.push_back((TH1D*)prior1_hvector[3]->Clone()); 
  priorunf_measdata_rat_vector.push_back((TH1D*)prior2_hvector[3]->Clone()); 
  priorunf_measdata_rat_vector.push_back((TH1D*)prior3_hvector[3]->Clone());
  

  priorunf_measdata_rat_vector[0]->Divide(h_measdata);
  priorunf_measdata_rat_vector[1]->Divide(h_measdata);
  priorunf_measdata_rat_vector[2]->Divide(h_measdata);
  
  priorunf_measdata_rat_vector[0]->SetLineColor(kRed);      priorunf_measdata_rat_vector[0]->SetMarkerColor(kRed);    
  priorunf_measdata_rat_vector[1]->SetLineColor(kBlue);     priorunf_measdata_rat_vector[1]->SetMarkerColor(kBlue);   
  priorunf_measdata_rat_vector[2]->SetLineColor(kMagenta);  priorunf_measdata_rat_vector[2]->SetMarkerColor(kMagenta);
  
  drawhArrayOnRatioCanv(priorXunf_measdata_canv,priorunf_measdata_rat_vector, meth_name+" Unf. Data, Ratio w/ Meas. Data", "Unf. Data / Meas. Data");  
  
  priorunf_measdata_rat_vector[0]->SetTitle( (prior1_desc+" Unf. Data").c_str() );
  priorunf_measdata_rat_vector[1]->SetTitle( (prior2_desc+" Unf. Data").c_str() );
  priorunf_measdata_rat_vector[2]->SetTitle( (prior3_desc+" Unf. Data").c_str() );
  
  TLegend* priorunf_measdata_legend=buildRatioLegend(priorunf_measdata_rat_vector);
  priorXunf_measdata_canv->cd();
  priorunf_measdata_legend->Draw();
  
  writeNsaveCanvas(fout, priorXunf_measdata_canv);
  

  
  
  
  //PRIOR X OVER PRIORX UNF
  TCanvas* priorXunf_priorX_canv=new TCanvas((meth_name+"_"+prior1_name+"."+prior2_name+"."+prior3_name+"_priorrat").c_str(),
					     (meth_name+"_"+prior1_name+"."+prior2_name+"."+prior3_name+"_priorrat").c_str(),
					     CANVXW, CANVYW);
  setupRatioCanv(priorXunf_priorX_canv);
  
  std::vector<TH1D*> prior_rat_vector;
  prior_rat_vector.push_back( (TH1D*)prior1_hvector[1]->Clone());
  prior_rat_vector.push_back( (TH1D*)prior2_hvector[1]->Clone());
  prior_rat_vector.push_back( (TH1D*)prior3_hvector[1]->Clone());
  
  prior_rat_vector[0]->Divide((TH1D*)prior1_hvector[3]->Clone());
  prior_rat_vector[1]->Divide((TH1D*)prior2_hvector[3]->Clone());
  prior_rat_vector[2]->Divide((TH1D*)prior3_hvector[3]->Clone());
  
  prior_rat_vector[0]->SetLineColor(kRed);      prior_rat_vector[0]->SetMarkerColor(kRed);    
  prior_rat_vector[1]->SetLineColor(kBlue);     prior_rat_vector[1]->SetMarkerColor(kBlue);   
  prior_rat_vector[2]->SetLineColor(kMagenta);  prior_rat_vector[2]->SetMarkerColor(kMagenta);
  
  drawhArrayOnRatioCanv(priorXunf_priorX_canv, prior_rat_vector, meth_name+" MC Priors, Ratio w/ Unf. Data", "MC Prior / Unf. Data");  
  
  prior_rat_vector[0]->SetTitle( (prior1_desc+" MC, Prior / Unf. Data").c_str() );
  prior_rat_vector[1]->SetTitle( (prior2_desc+" MC, Prior / Unf. Data").c_str() );
  prior_rat_vector[2]->SetTitle( (prior3_desc+" MC, Prior / Unf. Data").c_str() );
  
  TLegend* priorXunf_priorX_legend=buildRatioLegend(prior_rat_vector);
  priorXunf_priorX_canv->cd();
  priorXunf_priorX_legend->Draw();
  
  writeNsaveCanvas(fout, priorXunf_priorX_canv);
  
  
  
  //assert(false);  
  
  
  
  //compare prior1 unf data w/ prior2 unf data
  //compare prior2 unf data w/ prior3 unf data
  //compare prior1 unf data w/ prior3 unf data
  
  
  
  
 
  //std::cout<<"helloworld"<<std::endl;
  return;
}



//i.e. match the priors, but make sure diff methods
void quick_comparePriors(TFile* fin_p8byes, TFile* fin_NLObyes, TFile* fin_NLOwNPbyes,
			 TFile* fin_p8svd, TFile* fin_NLOsvd, TFile* fin_NLOwNPsvd){
  std::cout<<std::endl<<"running quick_comparePriors"<<std::endl;  
  
  //std::cout<<"comparing data unfolded w/ same methods yet diff priors"<<std::endl;    
  
  ////-----------------------------------------//
  std::cout<<std::endl<<"comparing SVD unfolded data across PY8, NLO, and NLO+NP priors"<<std::endl;
  print_comparePriors(     (TFile*)fin_p8svd,     PY8_SVD_NAMES,    "PY8_LO",        "PY8 LO",       
			   (TFile*)fin_NLOsvd,     NLO_SVD_NAMES,    "NNPDF_NNLO",    "NNPDF NNLO",   
			   (TFile*)fin_NLOwNPsvd,  NLOwNP_SVD_NAMES,    "NNPDF_NNLOwNP", "NNPDF NNLO+NP",
			   "SVD", 
			   (std::vector<std::string>)HNAMES );
  ////-----------------------------------------//


  //assert(false);
  //-----------------------------------------//
  std::cout<<std::endl<<"comparing bayes unfolded data across PY8, NLO, and NLO+NP priors"<<std::endl;
  print_comparePriors( (TFile*)fin_p8byes,   PY8_bayes_NAMES,          "PY8_LO",        "PY8 LO",       
		       (TFile*)fin_NLObyes,  NLO_bayes_NAMES,          "NNPDF_NNLO",    "NNPDF NNLO",   
		       (TFile*)fin_NLOwNPbyes,  NLOwNP_bayes_NAMES,    "NNPDF_NNLOwNP", "NNPDF NNLO+NP",
		       "Bayes", 
		       (std::vector<std::string>)HNAMES );
  //-----------------------------------------//


  
  return;
}






//void quickCompare_unfoldedData(bool compareMethods=false, bool comparePriors=true){
int main (int argc, char* argv[]){
  if(argc!=3)
    return 0;
  
  bool compareMethods=(bool)std::atoi(argv[1]);
  bool comparePriors=(bool)std::atoi(argv[2]);
  
  
  int rStatus=1;
  
  std::cout<<std::endl<<"running quickCompare_unfoldedData"<<std::endl;

  TFile* fin_PY8_bayes   =TFile::Open(FILE_PY8_bayes.c_str(),"READ");
  TFile* fin_NLO_bayes   =TFile::Open(FILE_NLO_bayes.c_str(),"READ");
  TFile* fin_NLOwNP_bayes=TFile::Open(FILE_NLOwNP_bayes.c_str(),"READ");
  
  TFile* fin_PY8_SVD   =TFile::Open(FILE_PY8_SVD.c_str(),"READ");
  TFile* fin_NLO_SVD   =TFile::Open(FILE_NLO_SVD.c_str(),"READ");
  TFile* fin_NLOwNP_SVD=TFile::Open(FILE_NLOwNP_SVD.c_str(),"READ");
  
  //TFile* fout_compmeth=new TFile("compare_methods.root","RECREATE");
  if(compareMethods)
    quick_compareMethods( fin_PY8_bayes, fin_NLO_bayes, fin_NLOwNP_bayes,
			  fin_PY8_SVD,   fin_NLO_SVD,   fin_NLOwNP_SVD);
  
  
  if(comparePriors)
    quick_comparePriors( fin_PY8_bayes, fin_NLO_bayes, fin_NLOwNP_bayes,
			 fin_PY8_SVD,   fin_NLO_SVD,   fin_NLOwNP_SVD  );    
  std::cout<<"done w/ all comparisons. exit."<<std::endl;
  
  return rStatus;
}

