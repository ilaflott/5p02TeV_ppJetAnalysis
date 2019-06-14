#include <string>
#include <iostream>
#include <cstdio>

#include "TROOT.h"
#include "TKey.h"
#include "TClass.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"

const bool DO_PNG=false, DO_PDF=true, DO_EPS=false, DO_C=false;
const bool debugMode=true;

void safe_getNdraw(TFile* f=NULL, std::string dirname="", int etabin=0){  
  
  TIter next(f->GetListOfKeys());
  TKey *key;
  
  std::string outputTag="";
  if(etabin==0)       outputTag="_00eta05";
  else if(etabin==1 ) outputTag="_05eta10";
  else if(etabin==2 ) outputTag="_10eta15";
  else if(etabin==3 ) outputTag="_15eta20";
  else if(etabin==4 ) outputTag="_20eta25";
  else if(etabin==5 ) outputTag="_25eta30";
  else outputTag="_NULLTAG";
  
  
  //bool dir_created=false;
  while ( (key = (TKey*)next()) ) {//loop over all objects in the root file
    
    if(debugMode)std::cout<<"found key! TKey->GetName()="<<key->GetName()<<std::endl;//debug
    
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TCanvas")) continue;//if its not a canvas forget it
    
    TCanvas* canv = (TCanvas*)key->ReadObj();
    if((bool)canv){if(debugMode)std::cout<<"found TCanvas! name is "<<key->GetName()<<std::endl;    }
    else {std::cout<<"canvas not opening for some reason. skip."<<std::endl;continue;}
    
    canv->Draw();
    
    if(DO_PNG){
      std::cout<<"saving .png of "<< key->GetName()<<std::endl;
      std::string png_name=dirname+((std::string)key->GetName())+outputTag+".png";
      std::cout<<"png_name="<<png_name<<std::endl;
      canv->SaveAs( png_name.c_str() );}
    if(DO_EPS){
      std::cout<<"saving .eps of "<< key->GetName()<<std::endl;
      std::string eps_name=dirname+((std::string)key->GetName())+outputTag+".eps";
      canv->SaveAs( eps_name.c_str() );}
    if(DO_PDF){
      std::cout<<"saving .pdf of "<< key->GetName()<<std::endl;
      std::string pdf_name=dirname+((std::string)key->GetName())+outputTag+".pdf";
      canv->SaveAs( pdf_name.c_str() );}
    if(DO_C){
      std::cout<<"saving .C of   "<< key->GetName()<<std::endl;
      std::string C_name  =dirname+((std::string)key->GetName())+outputTag+".C";
      canv->SaveAs( C_name.c_str() );}
    
    //return; //use for DEBUG if running by root -l <name>.C
    //assert(false);//use for DEBUG if running by ./<name>.exe
    
  }
  if(debugMode)std::cout<<"done w/ key loop"<<std::endl;
  
  return;
}

//-----------------------------------------------
void print_gausSmThy_etabin(std::string nlodir="", std::string nlofile="", int etaBin=0){
  
  std::string rootFile=nlodir+nlofile;

  std::cout<<"opening file :" <<rootFile<<std::endl;
  TFile* file=TFile::Open(rootFile.c_str());
  if( !((bool)file) ){
    cout<<"TFILE "<<rootFile<<" NOT OPEN/FOUND. EXIT"<<endl;
    return;  }
  
  file->cd();
  file->ls();      
  
  safe_getNdraw(file,nlodir, etaBin);
  
  
  
  return;
}
//-----------------------------------------------


void quickCheck_gausSmThy_etabin( std::string dir="", std::string file="", int etaBin=0){
  
  std::cout<<"quickly printing gausSmThy hists" <<std::endl;
  print_gausSmThy_etabin(dir, file, etaBin);
  
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


