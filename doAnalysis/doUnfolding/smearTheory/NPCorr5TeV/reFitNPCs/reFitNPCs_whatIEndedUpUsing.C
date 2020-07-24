#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <iostream>
#include <ostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstring>
// ROOTSYS
#include <TSystem.h>
#include <TProfile.h>
// I/O
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TEventList.h>
#include <TCut.h>
// utilities
#include <TMath.h>
#include <TRandom3.h>
#include <TStopwatch.h>
// plotting
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TROOT.h>
// histos
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
//
#include <TVirtualFitter.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>
#include <Math/MinimizerOptions.h>
#include <Math/IOptions.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMatrix.h>
#include <TMatrixT.h>
#include <TMatrixD.h>

const std::string absy_str[]={
  "0.0 < #||{y} < 0.5",
  "0.5 < #||{y} < 1.0",
  "1.0 < #||{y} < 1.5",
  "1.5 < #||{y} < 2.0"
};
const int Nybins=4;


std::string safeFtostr(double val){
  std::stringstream val_sstream;
  val_sstream<<val;
  std::string val_str=val_sstream.str();
  return val_str;
}

// I put this here to keep track of what the bins should be given the NLO calculations from Joao
std::vector<std::vector<double>>SMP_ptbins={
  {43,49,56,64,74,
   84,97,114,133,153,
   174,196,220,245,272,
   300,330,362,395,430,
   468,507,548,592,638,
   686,737,790,846,905,
   967,1032,1101},
  {43,49,56,64,74,
   84,97,114,133,153,
   174,196,220,245,272,
   300,330,362,395,430,
   468,507,548,592,638,
   686,737,790,846,905,
   967,1032},
  {43,49,56,64,74,
   84,97,114,133,153,
   174,196,220,245,272,
   300,330,362,395,430,
   468,507,548,592,638,
   686,737,790,846},
  {43,49,56,64,74,
   84,97,114,133,153,
   174,196,220,245,272,
   300,330,362,395,430,
   468,507,548,592,638,
   686,737,790}
};




void reFitNPCs(){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gROOT->ForceStyle();
  TVirtualFitter::SetDefaultFitter("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(100);
  ROOT::Math::MinimizerOptions opt;
  opt.Print();
  
  TFile* fin=TFile::Open("NLOpNP_InclusiveJets5TeV_fromJohn.root","READ");
  TH1* data[Nybins]={};
  
  TFile* fout=TFile::Open("NLOpNP_InclusiveJets5TeV_fromJohn_refit.root","RECREATE");
  TF1* nufits[Nybins]={};
  TF1* nufits_up[Nybins]={};
  TF1* nufits_down[Nybins]={};
  std::string fitfunc="[0]+[1]*pow(x,[2])";

  for(int i=0; i<Nybins; i++){
    
    data[i]=(TH1*)fin->Get(("h"+std::to_string(i)).c_str());

    //int lowestbin=7;//even higher confidence, even worse error bars
    //int lowestbin=6;//about the same as 5
    //int lowestbin=5;//higher confidence, worse error band
    int lowestbin=4;//best results so far
    //int lowestbin=3;//worse than 4
    //int lowestbin=2;//meh, no improvement

    int binstocut=0;
    if(i==0)binstocut=0;
    else if(i==1)binstocut=0;
    else if(i==2)binstocut=1;
    else if(i==3)binstocut=6;
    
    // i do this rebinning in order to cut off limited statistics at high pT and the region at low pT we don't need.
    // no bins are actually being combined here. i'm just concatenating the range for the fit.
    float xmin=data[i]->GetBinLowEdge(lowestbin);
    float xmax=-1.;
    std::vector<double>ptbins;
    for(int j=1;j<=(data[i]->GetNbinsX()-binstocut);j++){
      if(j<lowestbin){
	std::cout<<"bin j="<<j<<" below lowest pt bin to use, continue."<<std::endl;
	continue;
      }
      else if(data[i]->GetBinContent(j)>0.){
	xmax=data[i]->GetBinLowEdge(j)+data[i]->GetBinWidth(j);
	ptbins.push_back(data[i]->GetBinLowEdge(j));
	std::cout<<"ptbins.at("<<j-lowestbin<<")="<<ptbins.at(j-lowestbin)<<std::endl;
      }
      else {
	std::cout<<"bin "<<data[i]->GetBinLowEdge(j)<<" < pt < " << data[i]->GetBinLowEdge(j) +data[i]->GetBinWidth(j)<< " has no bin content, continue."<<std::endl;
	continue;
      }
    }
    ptbins.push_back(xmax);
    std::cout<<"ptbins.at("<<ptbins.size()-1<<")="<<ptbins.at(ptbins.size()-1)<<std::endl;  
    data[i]=(TH1*)data[i]->TH1::Rebin( (ptbins.size()-1) ,((std::string)data[i]->GetName()+"_rebin").c_str(),(Double_t*)ptbins.data());
    


    
    
    //this worked pretty good!
    nufits[i]=new TF1(("f"+std::to_string(i)).c_str(),fitfunc.c_str(), xmin,xmax);
    nufits[i]->SetParameter(0, 1.0);
    nufits[i]->SetParameter(1,-0.1);
    nufits[i]->SetParameter(2,-0.1);
    std::cout<<"par0 init="<<nufits[i]->GetParameter(0)<<std::endl;
    std::cout<<"par1 init="<<nufits[i]->GetParameter(1)<<std::endl;
    std::cout<<"par2 init="<<nufits[i]->GetParameter(2)<<std::endl;
    nufits[i]->SetParLimits(0,   0.5,2.);
    nufits[i]->SetParLimits(1,-100000.,0.);
    nufits[i]->SetParLimits(2,-100000.,0.);


    TFitResultPtr resultptr=(TFitResultPtr)(data[i]->Fit(nufits[i],"RES"));
    resultptr=(TFitResultPtr)(data[i]->Fit(nufits[i],"MRES"));
    resultptr->Print("V");

    // NOTE TO JOHN: All of this stuff involving the parameters and putting together a new formula string was to try to do 
    // a genuine error calculation that accounted for correlations between the variables.
    // this ended up convoluted and complicated and not worth it, since the error band it created was even *more* optimistic 
    // than the original error band you came up with. Because of this, i'm guessing it's not the way to go
    // BUT if you're curious what i did here, all I did was use the full error propagation formula;
    // let f be the fit function, i.e. f = A + B x^C
    // delta(NP) = sqrt(  (df/dA*delta(A))^2 +
    //                    (df/dB*delta(B))^2 +
    //                    (df/dC*delta(C))^2 +
    //                    2(df/dA)(df/dB)cov(A,B) +
    //                    2(df/dB)(df/dC)cov(B,C) +
    //                    2(df/dC)(df/dA)cov(C,A) +
    //                 )
    // the other mess here was that the usual conversion of a # to a string using std::to_string doesn't work well for floats and/or doubles.
    // to get around that, one needs to use the std::stringstream class, the function definition is above at the top of the script to neatness.
    TMatrixD covmat=(TMatrixD)resultptr->GetCovarianceMatrix();
    double errfact=3.;//errfact=Nsigma
    double parA   =nufits[i]->GetParameter(0);
    double parAerr=errfact*nufits[i]->GetParError( 0);
    double parB   =nufits[i]->GetParameter(1);
    double parBerr=errfact*nufits[i]->GetParError( 1);
    double parC   =nufits[i]->GetParameter(2);
    double parCerr=errfact*nufits[i]->GetParError( 2);

    std::string parA_str   =safeFtostr(parA   );
    std::string parAerr_str=safeFtostr(parAerr);
    std::string parB_str   =safeFtostr(parB   );
    std::string parBerr_str=safeFtostr(parBerr);
    std::string parC_str   =safeFtostr(parC   );
    std::string parCerr_str=safeFtostr(parCerr);
    //////DEBUG
    //std::cout<<"               par0     ="<<               par0     <<std::endl;
    //std::cout<<"std::to_string(par0)    ="<<std::to_string(par0)    <<std::endl;
    //std::cout<<"               par0_str ="<<               par0_str <<std::endl;
    //assert(false);
    //////DEBUG 
    

    double parABcov=errfact*errfact*covmat(0,1);
    double parBCcov=errfact*errfact*covmat(1,2);
    double parCAcov=errfact*errfact*covmat(2,0);
    std::string parABcov_str=safeFtostr(parABcov);
    std::string parBCcov_str=safeFtostr(parBCcov);
    std::string parCAcov_str=safeFtostr(parCAcov);
    //////DEBUG 
    //std::cout<<"               parABcov     ="<<               parABcov     <<std::endl;
    //std::cout<<"std::to_string(parABcov)    ="<<std::to_string(parABcov)    <<std::endl;
    //std::cout<<"               parABcov_str ="<<               parABcov_str <<std::endl;
    //std::cout<<"               parACcov     ="<<               parACcov     <<std::endl;
    //std::cout<<"std::to_string(parACcov)    ="<<std::to_string(parACcov)    <<std::endl;
    //std::cout<<"               parACcov_str ="<<               parACcov_str <<std::endl;
    //////DEBUG 



    std::string dfdA="(1.)";
    std::string dfdB="(pow(x,"+parC_str+"))";
    std::string dfdC="("+parB_str+"*log(x)*pow(x,"+parC_str+"))";
    std::cout<<"dfdA = "<< dfdA<<std::endl<<std::endl;		 
    std::cout<<"dfdB = "<< dfdB<<std::endl<<std::endl;		 
    std::cout<<"dfdC = "<< dfdC<<std::endl<<std::endl<<std::endl;
    std::string dfdA_dfdB="("+dfdA+"*"+dfdB+")";
    std::string dfdB_dfdC="("+dfdB+"*"+dfdC+")";
    std::string dfdC_dfdA="("+dfdC+"*"+dfdA+")";
    std::cout<<"dfdA_dfdB = "<< dfdA_dfdB<<std::endl<<std::endl;		 
    std::cout<<"dfdB_dfdC = "<< dfdB_dfdC<<std::endl<<std::endl;		 
    std::cout<<"dfdC_dfdA = "<< dfdC_dfdA<<std::endl<<std::endl<<std::endl;
    
    std::string delA_dfdA_sq="pow(("+dfdA+"*"+parAerr_str+"),2)";
    std::string delB_dfdB_sq="pow(("+dfdB+"*"+parBerr_str+"),2)";
    std::string delC_dfdC_sq="pow(("+dfdC+"*"+parCerr_str+"),2)";
    std::cout<< "delA_dfdA_sq=" << delA_dfdA_sq<< std::endl<<std::endl;		 
    std::cout<< "delB_dfdB_sq=" << delB_dfdB_sq<< std::endl<<std::endl;		 
    std::cout<< "delC_dfdC_sq=" << delC_dfdC_sq<< std::endl<<std::endl<<std::endl;
    
    std::string dfdA_dfdB_ABcov="(2.*"+dfdA+"*"+dfdB+"*"+parABcov_str+")";
    std::string dfdB_dfdC_BCcov="(2.*"+dfdB+"*"+dfdC+"*"+parBCcov_str+")";
    std::string dfdC_dfdA_CAcov="(2.*"+dfdC+"*"+dfdA+"*"+parCAcov_str+")";
    std::cout<< "dfdA_dfdB_ABcov=" << dfdA_dfdB_ABcov<< std::endl<<std::endl;		 
    std::cout<< "dfdB_dfdC_BCcov=" << dfdB_dfdC_BCcov<< std::endl<<std::endl;		 
    std::cout<< "dfdC_dfdA_CAcov=" << dfdC_dfdA_CAcov<< std::endl<<std::endl<<std::endl;
    //assert(false);
    
    //this approach seems to underestimate the error pretty hard, unless i've made an error somewhere (possible...)
    std::string totalerr_str="sqrt("+
      delA_dfdA_sq+"+"+
      delB_dfdB_sq+"+"+
      delC_dfdC_sq+"+"+
      dfdA_dfdB_ABcov+"+"+
      dfdB_dfdC_BCcov+"+"+
      dfdC_dfdA_CAcov+")";
    std::cout<<"totalerr_str="<<totalerr_str<<std::endl;

    
    //nufits_up[i]=new TF1(("f"+std::to_string(i)+"_up").c_str(), (fitfunc+"+"+totalerr_str).c_str(),xmin,xmax);
    //nufits_up[i]->SetParameter(0,nufits[i]->GetParameter(0));
    //nufits_up[i]->SetParameter(1,nufits[i]->GetParameter(1));
    //nufits_up[i]->SetParameter(2,nufits[i]->GetParameter(2));  

    nufits_up[i]=(TF1*)nufits[i]->Clone((("f"+std::to_string(i)+"_up").c_str()));
    nufits_up[i]->SetParameter(0,nufits[i]->GetParameter(0)+nufits[i]->GetParError(0));
    nufits_up[i]->SetParameter(1,nufits[i]->GetParameter(1)+nufits[i]->GetParError(1));
    nufits_up[i]->SetParameter(2,nufits[i]->GetParameter(2)-nufits[i]->GetParError(2));
    nufits_up[i]->SetLineColorAlpha(kBlue,0.6);
    nufits_up[i]->SetLineWidth(1);
    nufits_up[i]->SetLineStyle(8);
    
    //nufits_down[i]=new TF1(("f"+std::to_string(i)+"_down").c_str(), (fitfunc+"-"+totalerr_str).c_str(),xmin,xmax);
    //nufits_down[i]->SetParameter(0,nufits[i]->GetParameter(0));
    //nufits_down[i]->SetParameter(1,nufits[i]->GetParameter(1));
    //nufits_down[i]->SetParameter(2,nufits[i]->GetParameter(2));  

    nufits_down[i]=(TF1*)nufits[i]->Clone((("f"+std::to_string(i)+"_down").c_str()));
    nufits_down[i]->SetParameter(0,nufits[i]->GetParameter(0)-nufits[i]->GetParError(0));
    nufits_down[i]->SetParameter(1,nufits[i]->GetParameter(1)-nufits[i]->GetParError(1));
    nufits_down[i]->SetParameter(2,nufits[i]->GetParameter(2)+nufits[i]->GetParError(2));
    nufits_down[i]->SetLineColorAlpha(kBlue,0.6);
    nufits_down[i]->SetLineWidth(1);
    nufits_down[i]->SetLineStyle(8);
    
    // sanity check
    std::cout<<"nufits[i]     ->Eval(100)="<< nufits[i]     ->Eval(100.)<<std::endl;
    std::cout<<"nufits_up[i]  ->Eval(100)="<< nufits_up[i]  ->Eval(100.)<<std::endl;
    std::cout<<"nufits_down[i]->Eval(100)="<< nufits_down[i]->Eval(100.)<<std::endl;
    
    //if(fitfailed)std::cout<<"WARNING FIT FAILED!!!!"<<std::endl;
    //else         std::cout<<"FIT SUCCESSFUL!!!!!"<<std::endl;
    //if(fitfailed)std::cout<<""<<std::endl;
    //else         std::cout<<""<<std::endl;
    ptbins.clear();
    
    data[i]->SetTitle(absy_str[i].c_str());
    data[i]->GetXaxis()->SetMoreLogLabels(true);
    data[i]->GetXaxis()->SetNoExponent(true);
    data[i]->GetXaxis()->SetTitle("Jet p_{T} (GeV)");      
    data[i]->SetMinimum(0.90);
    data[i]->SetMaximum(1.05);

    nufits[i]->SetLineColorAlpha(kRed,0.8);
    nufits[i]->SetLineWidth(1);

  }
    
  TCanvas *canv=new TCanvas("canv", "canv", 1200, 900);
  canv->Divide(2,2);
  for(int i=0; i<Nybins;i++){
    canv->cd(i+1)->SetLogx(1);
    canv->cd(i+1)->SetLogy(0);
    canv->cd(i+1);        
    data[i]->Draw("HIST E ][");
    nufits[i]->Draw("SAME");
    nufits_up[i]->Draw("SAME");
    nufits_down[i]->Draw("SAME");
  }  
  canv->SaveAs( ("NLOpNP_InclusiveJets5TeV_fromJohn_refit.pdf"));
  
  fout->cd();
  canv->Write();
  for(int i=0; i<Nybins; i++){
    nufits[i]->Write();
    nufits_up[i]->Write();
    nufits_down[i]->Write();
    data[i]->Write(("h"+std::to_string(i)).c_str());
  }
  fout->Close();


  return;
}

// NOTE TO JOHN: i do an int main function because i'm generally using these scripts by compiling with g++ and linking to the ROOT libraries.
// this means this script won't work by simply doing 'root -l -q reFitNPCs.C++' since ROOT doesn't know what to do with int main.
// That being said, if you comment this out, and do 'root -l -q reFitNPCs.C++', it should work as expected.
int main(int argc, char*argv[]){
  reFitNPCs(); 
  return 0;
}









