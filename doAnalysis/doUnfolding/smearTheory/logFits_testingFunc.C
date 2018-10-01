#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TSystem.h"
#include "TF1.h"
#include "TH1.h"
#include <TVirtualFitter.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>


#include <cstdlib>
#include <cstring>
const double anabins_pt[] = {//analysis pt bins
  43.,   49.,   56.,   64.,   74.,   84.,   97., 
  114.,   133.,   153.,   174.,   196.,   220., 
  245.,   272.,   300.,   330.,     362.,   395., 
  430.,   468.,   507.,   548.,   592.,   638., 
  686. 
};
const int n_anabins_pt=sizeof(anabins_pt)/sizeof(double)-1;

//hacky fix, avoid in future.
void specialRebin(TH1D* hist ){
  std::string histname=(std::string)hist->GetName();
  hist=(TH1D*)hist->TH1::Rebin( n_anabins_pt, histname.c_str(), anabins_pt);
  return;
}

TH1D* safeClone(TH1D* hist, int nameint){
  std::cout<<"running safeClone"<<std::endl;
  
  std::string h_safeclone_name=((std::string)hist->GetName())+"_clone"+std::to_string(nameint);
  TH1D* h_safeclone=(TH1D*)hist->Clone(h_safeclone_name.c_str());
  h_safeclone->Reset("M ICES");
  
  for(int i=1; i<= hist->GetNbinsX();i++){
    h_safeclone->SetBinContent(i,hist->GetBinContent(i));
    h_safeclone->SetBinError(i,hist->GetBinError(i));
  }
  
  return (TH1D*)h_safeclone;
}


int main(int argc, char* argv[]){
  int rstatus=1;
  
  if(argc!=2) {
    std::cout<<"not enough arguments. do ./temp.exe <param bin>"<<std::endl;
    return -1;
  }
  
  TFile* _file0=TFile::Open("fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_gaussSmear_00eta20.root");  
  if(!_file0){
    std::cout<<"file not found. exit."<<std::endl;
    return -1;
  }
  
  //warning, hist in this/these files misnamed slightly (clean up in smearTheory code when you can)
  //TH1D* hist2Fit=(TH1D*)_file0->Get("h00eta20_wNP_specialrebin;2"); //w NP
  
  TH1D* hist2Fit=(TH1D*)_file0->Get("h00eta20_wNP;1"); //no NP, despite the hist object name...
  std::cout<<"hist2Fit nbinsx="<<hist2Fit->GetNbinsX()<<std::endl;
  hist2Fit=(TH1D*)hist2Fit->TH1::Rebin( n_anabins_pt, "h00eta20_rebin", anabins_pt);
  std::cout<<"hist2Fit nbinsx="<<hist2Fit->GetNbinsX()<<std::endl;
  std::cout<<"compare to n_anabins_pt="<<n_anabins_pt<<std::endl;
  //  assert(false);
  if(!hist2Fit){
    std::cout<<"hist not found. exit."<<std::endl;
    return -1;
  }  
  
  gStyle->SetOptFit(101);
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();    
  
  TVirtualFitter::SetDefaultFitter("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04);    
  
  int nbins=hist2Fit->GetNbinsX();
  int parambin=std::atoi(argv[1]);
  if(parambin <1 || parambin > nbins){
    std::cout<<"compare parambin w/ nbins="<<nbins<<std::endl;
    std::cout<<"input parambin="<<parambin<<" too large or small. exit."<<std::endl;
    return -1;    
  }
  
  

  //  TH1D* hist2Fit_clone1=(TH1D*)hist2Fit->Clone( ((std::string)hist2Fit->GetName()+"_clone1").c_str() );
  //  TH1D* hist2Fit_clone2=(TH1D*)hist2Fit->Clone( ((std::string)hist2Fit->GetName()+"_clone2").c_str() );
  TH1D* hist2Fit_clone1=(TH1D*)safeClone(hist2Fit,1);
  TH1D* hist2Fit_clone2=(TH1D*)safeClone(hist2Fit,2);    
  TH1D* hist2Fit_clone3=(TH1D*)safeClone(hist2Fit,3);    
  
  int startbin=1, endbin=hist2Fit->GetNbinsX();
  float xlo=hist2Fit->GetBinLowEdge(startbin);
  float xhi=hist2Fit->GetBinLowEdge(endbin) + hist2Fit->GetBinWidth(endbin);
  std::cout<<"fit range is from xlo="<<xlo<<" to xhi="<<xhi<<std::endl;        
  
  
  //  std::string fout_new_str="logFits_NNPDF_NNLO_POWPY8_NPs_param"+std::to_string(parambin)+".root";
  
  std::string fout_new_str="logFits_NNPDF_NNLO_param"+std::to_string(parambin)+".root";
  
  TFile* fout_new=new TFile(fout_new_str.c_str(),"RECREATE");  
  std::cout<<"grabbing some initial fit parameters from bin #"<<parambin<<std::endl;  
  float logFitp0=hist2Fit->GetBinContent(parambin);  std::cout<<"fit param 0 being set to="<<logFitp0<<std::endl;  //parameter 0 val of xsec
  float logFitp1=hist2Fit->GetBinCenter(parambin);  std::cout<<"fit param 1 being set to="<<logFitp1<<std::endl;  //parameter 1 val bin center of param 0's bin
 //basically, using my hist, in logx/logy coordinates, this is a linear fit
  TF1* logFit=new TF1("logfit",    "[0]*( TMath::Power( ( [1]/x ) , [2] ) )",xlo,xhi);
  logFit->SetParameter(0,logFitp0);
  logFit->SetParameter(1,logFitp1);
  hist2Fit_clone1->Fit("logfit","VMRE");  


  
  float modLogFitp0=logFit->GetParameter(0);
  float modLogFitp1=logFit->GetParameter(1);
  float modLogFitp2=logFit->GetParameter(2);
  //thus, this is a parabolic fit in logx/logy
  TF1* modLogFit=new TF1("modlogfit","[0]*( TMath::Power( ( [1]/x ) , [2]+[3]*TMath::Log10( x/[1] ) ) )",xlo,xhi);
  modLogFit->SetParameter(0,modLogFitp0);
  modLogFit->SetParameter(1,modLogFitp1);
  modLogFit->SetParameter(2,modLogFitp2);
  hist2Fit_clone2->Fit("modlogfit","VMRE");  
  
  
  
  float modLog2Fitp0=modLogFit->GetParameter(0);
  float modLog2Fitp1=modLogFit->GetParameter(1);
  float modLog2Fitp2=modLogFit->GetParameter(2);
  float modLog2Fitp3=modLogFit->GetParameter(3);
  //thus, this is a cubic fit fit in logx/logy
  TF1* modLog2Fit=new TF1("modlog2fit", "[0] * ( TMath::Power( ( [1]/x ) , ( [2] + [3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )  ) ) )",xlo,xhi);
  modLog2Fit->SetParameter(0,modLog2Fitp0);
  modLog2Fit->SetParameter(1,modLog2Fitp1);
  modLog2Fit->SetParameter(2,modLog2Fitp2);
  modLog2Fit->SetParameter(3,modLog2Fitp3);
  hist2Fit_clone3->Fit("modlog2fit","VMRE");  
  
  
  
  //TCanvas* theFIRSTcanv=new TCanvas("theFIRSTcanv","theFIRSTcanv",800,800);
  //theFIRSTcanv->SetLogx(1);
  //theFIRSTcanv->SetLogy(1);
  //theFIRSTcanv->cd();  
  
  //TH1D* hmodLogFit=(TH1D*)modLogFit->GetHistogram();
  //hmodLogFit->SetLineColor(kRed);
  //hmodLogFit->Draw("HIST");   
  //hist2Fit_clone1->Draw("HIST E");


  //TCanvas* theSECONDcanv=new TCanvas("theSECONDcanv","theSECONDcanv",800,800);
  //theSECONDcanv->SetLogx(1);
  //theSECONDcanv->SetLogy(1);
  //theSECONDcanv->cd();  
  
  //TH1D* hmodLogFit_4pms=(TH1D*)modLogFit_4pms->GetHistogram();
  //hmodLogFit_4pms->SetLineColor(kBlue);
  //hmodLogFit_4pms->Draw("HIST");    
  //hist2Fit_clone2->Draw("HIST E");
  
  fout_new->cd();
  
  
//  std::string logfit_name=        "h00eta20_logfit_NNPDF_NNLO_POWPY8_NPs_param"+std::to_string(parambin);
//  std::string modlogfit_name=  "h00eta20_modlogfit_NNPDF_NNLO_POWPY8_NPs_param"+std::to_string(parambin);  
//  std::string modlog2fit_name="h00eta20_modlog2fit_NNPDF_NNLO_POWPY8_NPs_param"+std::to_string(parambin);  

  std::string logfit_name=        "h00eta20_logfit_NNPDF_NNLO_param"+std::to_string(parambin);
  std::string modlogfit_name=  "h00eta20_modlogfit_NNPDF_NNLO_param"+std::to_string(parambin);  
  std::string modlog2fit_name="h00eta20_modlog2fit_NNPDF_NNLO_param"+std::to_string(parambin);  
  
  TCanvas* logfit_canv=new TCanvas("logfit_canv","logfit_canv",800,800);
  logfit_canv->cd();
  logfit_canv->SetLogx(1);  logfit_canv->SetLogy(1);
  hist2Fit_clone1->Draw();
  logfit_canv->SaveAs( (logfit_name + ".png").c_str() ) ;
  logfit_canv->Write( (logfit_name + "_canv").c_str());
  hist2Fit_clone1->Write(logfit_name.c_str());
  logFit->Write();
  

  

  
  TCanvas* modlogfit_canv=new TCanvas("modlogfit_canv","modlogfit_canv",800,800);
  modlogfit_canv->cd();
  modlogfit_canv->SetLogx(1);  modlogfit_canv->SetLogy(1);
  hist2Fit_clone2->Draw();
  modlogfit_canv->SaveAs( (modlogfit_name + ".png").c_str() ) ;
  modlogfit_canv->Write( (modlogfit_name + "_canv").c_str());
  hist2Fit_clone2->Write(modlogfit_name.c_str());
  modLogFit->Write();








  
  TCanvas* modlog2fit_canv=new TCanvas("modlog2fit_canv","modlog2fit_canv",800,800);
  modlog2fit_canv->cd();
  modlog2fit_canv->SetLogx(1);  modlog2fit_canv->SetLogy(1);
  hist2Fit_clone3->Draw();
  modlog2fit_canv->SaveAs( (modlog2fit_name + ".png").c_str() ) ;
  modlog2fit_canv->Write( (modlog2fit_name + "_canv").c_str());
  hist2Fit_clone3->Write(modlog2fit_name.c_str());
  modLog2Fit->Write();
  




  //hist2Fit_clone2->SaveAs( (modlogfit_name + ".png").c_str());
  //hist2Fit_clone2->Write(   modlogfit_name.c_str());
  //modLogFit->Write();
  

  fout_new->Close();
  return rstatus;
}
