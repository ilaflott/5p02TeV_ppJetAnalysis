#include <string>
#include <iostream>
#include <cstdio>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
// QUICK AND DIRTY ROOT USAGE NOTES

////for a C++ function
// quickCheck_gausSmThy(int input=defValue)
////written in the file with the same name
// quickCheck_gausSmThy.C
////the following applies...

////recompile + run w/ input, usage:
//root -l 'quickCheck_gausSmThy.C++(1)'
////recompile + run w/ default input, usage:
//root -l 'quickCheck_gausSmThy.C++'
////compile only if not done so + run w/ default input, usage:
//root -l 'quickCheck_gausSmThy.C+'



void safe_getNdraw(TFile* f, std::string plotname ){//, std::string dname="default_gausSmThy_hists"){
  TCanvas* canv=(TCanvas*)f->Get(plotname.c_str());
  if(((bool)canv)){
    canv->Draw();
  }
  else{
    cout<<"PLOT NOT FOUND; SKIPPING PLOT "<<plotname<<endl;
    return;
  }
  std::string fname=(std::string)f->GetName();
  std::string dname="default_gausSmThy_hists";
  if(fname.find(     "CT10" )!=std::string::npos )  {canv->SaveAs( ("CT10nlo_"      +plotname  +".pdf" ).c_str());     }
  else if(fname.find("CT14" )!=std::string::npos )  {canv->SaveAs( ("CT14nlo_"      +plotname  +".pdf" ).c_str());     }
  else if(fname.find("HERAPDF")!=std::string::npos ){canv->SaveAs( ("HERAPDF15NLO_"   +plotname+".pdf").c_str());    }
  else if(fname.find("MMHT"   )!=std::string::npos ){canv->SaveAs( ("MMHT2014nlo68cl_"+plotname+".pdf"   ).c_str()); }
  else if(fname.find("NNPDF"  )!=std::string::npos ){
    
    if(fname.find("spl3wgts")!=std::string::npos )dname="NNPDF_NNLO_00eta20_spl3wgts_gausSmThy_plots/";
    else if(fname.find("fitwgts")!=std::string::npos )dname="NNPDF_NNLO_00eta20_fitwgts_gausSmThy_plots/";
    else    dname="NNPDF_NNLO_00eta20_gausSmThy_plots/";
    
    //canv->SaveAs( (dname+"NNPDF30_nnlo_"   +plotname+".png"  ).c_str());  
    canv->SaveAs( (dname+"NNPDF30_nnlo_"   +plotname+".pdf"  ).c_str());  

  }
  
  return;
}

//-----------------------------------------------
std::string plotnames[]={
  //NPcorr
  "plot_y0_to_y3_NPs",
  "plot_totNPs",
  "plot_allNPs",
  //JER
  "plotJER",
  //NLOxsec
  "plot_y0_to_y3_NLOxsec",
  "plot_totNLOxsec",
  "plot_allNLOxsec",
  //NPNLOxsec
  "plot_y0_to_y3_NPNLOxsec",
  "plot_totNPNLOxsec",
  "plot_allNPNLOxsec",
  //TSplines  
  "plot_splines",
  //TF1s
  "plot_logfits_ynew",
  "plot_logfits_NPynew",
  //NLOsmearing
  "plot_response_ynew_th2",
  "plot_ynew_true_smeared_rat",
  //NP+NLOsmearing
  "plot_response_NPynew_th2",
  "plot_NPynew_true_smeared_rat"  ,
  //compare rnd w/ true hist
  "plot_true_rnd_rat"
};
const int N_plotnames=sizeof(plotnames)/sizeof(std::string);
//-----------------------------------------------
void print_gausSmThy(std::string rootFile){
  
  //std::string rootFile="fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_00eta20.root";
  
  TFile* file=TFile::Open(rootFile.c_str());
  if( !((bool)file) ){
    cout<<"TFILE "<<rootFile<<" NOT OPEN/FOUND. EXIT"<<endl;
    return;  }
  
  file->cd();
  file->ls();
  std::cout<<"opening file :" <<rootFile<<std::endl;
  std::cout<<"drawing "<<N_plotnames<<" plots from file "<<rootFile<<std::endl;
  
  for(int i=0; i<N_plotnames;i++){
    std::cout<<"drawing plot #"<<i<<std::endl;
    safe_getNdraw(file,plotnames[i].data());
  }
  
  
  return;
}
//-----------------------------------------------




//-----------------------------------------------
std::string plotnames_singleBin[]={
  //NPcorr
  "plot_NP",
  //JER
  "plotJER",
  //NLOxsec
  "plot_NLOxsec",
  //NPNLOxsec
  "plot_NPNLOxsec",
  //TSplines  
  "plot_splines",
  //NLOsmearing
  "plot_response_th2",
  "plot_true_smeared_rat",
  //NP+NLOsmearing
  "plot_response_NP_th2",
  "plot_NP_true_smeared_rat"  
};
const int N_plotnames_singleBin=sizeof(plotnames_singleBin)/sizeof(std::string);
//-----------------------------------------------
void print_gausSmThy_singleBin(int etaBin){
  std::string etaBin_str=std::to_string(etaBin);
  std::string rootFile="fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_etaBin"+etaBin_str+".root";
  
  TFile* file=TFile::Open(rootFile.c_str());
  if( !((bool)file) ){
    cout<<"TFILE "<<rootFile<<" NOT OPEN/FOUND. EXIT"<<endl;
    return;  }
  
  file->cd();
  file->ls();
  std::cout<<"opening file :" <<rootFile<<std::endl;
  std::cout<<"drawing "<<N_plotnames_singleBin<<" plots from file "<<rootFile<<std::endl;
  
  for(int i=0; i<N_plotnames_singleBin;i++){
    std::cout<<"drawing plot #"<<i<<std::endl;
    std::string plotname=(std::string)plotnames_singleBin[i].data();
    if(plotname!="plotJER")
      plotname+=etaBin_str;
    
    safe_getNdraw(file,plotname);
  }
  return;
}
//-----------------------------------------------




void quickCheck_gausSmThy( bool print00eta20=true, bool printXXetaYY=false, int etabin=0){

  std::cout<<"quickly printing gausSmThy hists" <<std::endl;
  
  if(print00eta20){
    std::cout<<"printing 00eta20 gausSmThy hists"<<std::endl;
    
    //print_gausSmThy("fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_00eta20.root"             );
    //print_gausSmThy("fnl5020_LO2_R04Jets_modify_CT14nlo_HS_gaussSmear_00eta20.root"		);
    //print_gausSmThy("fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS_gaussSmear_00eta20.root"	);
    //print_gausSmThy("fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS_gaussSmear_00eta20.root"	);
    print_gausSmThy("fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_spl3wgts_gaussSmear_00eta20.root");
    print_gausSmThy("fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_gaussSmear_00eta20.root");
  }
  if(printXXetaYY){
    std::cout<<"printing gausSmThy hists for single eta bin..."<<std::endl;
    if(etabin==0){
      std::cout<<"ERROR give me an etabin to use."<<std::endl;
      std::cout<<"not printing hists for single eta bin!!"<<std::endl;
    }
    else{
      std::cout<<"etabin="<<etabin<<std::endl;
      print_gausSmThy_singleBin(etabin);
    }
  }
  
  return;
}
//  
//
//
//
//
//
//
//
//  
//  std::string rootFile="fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_00eta20.root";
//  
//  TFile* file=TFile::Open(rootFile.c_str());
//  if( !((bool)file) ){
//    cout<<"TFILE "<<rootFile<<" NOT OPEN/FOUND. EXIT"<<endl;
//    return;  }
//  
//  file->cd();
//  file->ls();
//  std::cout<<"opening file :" <<rootFile<<std::endl;
//  std::cout<<"drawing "<<N_plotnames<<" plots from file "<<rootFile<<std::endl;
//  
//  for(int i=0; i<N_plotnames;i++){
//    std::cout<<"drawing plot #"<<i<<std::endl;
//    safe_getNdraw(file,plotnames[i].data());
//  }
//  
//  
//  return;


