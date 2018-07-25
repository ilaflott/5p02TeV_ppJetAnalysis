


#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>

#include <TROOT.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>

#include <TH1.h>
#include <TH2.h>
#include <TH1D.h>
#include <TH2D.h>






int main(int argc, char* argv[]){  int rStatus = -1;

  if(argc!=3){
    std::cout<<"do ./quickPrintNLOSpectra.exe <targetFile> <etaBin>"<<std::endl;        
    
    return rStatus;  }
  
  const std::string NLOfilename=(const std::string)argv[1];
  std::cout<<"NLOfilename="<<NLOfilename<<std::endl;
  const int etaBin=(int)std::atoi(argv[2]);
  std::cout<<"etaBin="<<etaBin<<std::endl;
  
  
  if(etaBin > 7 || etaBin < 1)    {
    std::cout<<"<etaBin> argument must be int between 1 and 7. exit."<<std::endl;        
    return rStatus;
  }
  
  TFile* NLOfile=TFile::Open(NLOfilename.c_str());
  
  if(!NLOfile){
    std::cout<<"file doesn't exist or is not open. exit."<<std::endl;
    return rStatus;
  }
  
  
  std::string outputPdf=(std::string)NLOfilename;
  outputPdf.replace( (outputPdf.length())-5,(outputPdf.length()), "");
  outputPdf+="_bin"+std::to_string(etaBin);
  outputPdf+=".pdf";
  std::cout<<"outputPdf="<<outputPdf<<std::endl;
  
  std::string outputPdf_open=outputPdf+"[";
  std::string outputPdf_close=outputPdf+"]";
  TCanvas* canvForPrint= new TCanvas("theCanv","",1400,1400);  
  canvForPrint->Print(outputPdf_open.c_str());  
  
  for (int i = 0; i <=9; i++){
    if(i==3 || i==4 || i== 5)continue;
    
    std::string plotName="h0100"+std::to_string(etaBin)+"0"+std::to_string(i);
    std::cout<<"plotName="<<plotName<<std::endl; 
    canvForPrint->SetLogx(1);
    if(i==0)
      canvForPrint->SetLogy(1);
    else 
      canvForPrint->SetLogy(0);
    canvForPrint->cd();
    TH1D* thePlotToPrint=(TH1D*)NLOfile->Get(plotName.c_str());
    thePlotToPrint->SetTitle(plotName.c_str());
    thePlotToPrint->Draw("HIST E");
    canvForPrint->Print(outputPdf.c_str());
    
  }
  
  for (int i = 0; i <=9; i++){
    if(i==3 || i==4 || i== 5)continue;
    
    std::string plotName="h1100"+std::to_string(etaBin)+"0"+std::to_string(i);
    std::cout<<"plotName="<<plotName<<std::endl; 
    canvForPrint->SetLogx(1);
    if(i==0)
      canvForPrint->SetLogy(1);
    else 
      canvForPrint->SetLogy(0);
    canvForPrint->cd();
    TH1D* thePlotToPrint=(TH1D*)NLOfile->Get(plotName.c_str());
    thePlotToPrint->SetTitle(plotName.c_str());
    thePlotToPrint->Draw("HIST E");
    canvForPrint->Print(outputPdf.c_str());
    
  }
  
  
  canvForPrint->Print(outputPdf_close.c_str());
  std::cout<<"done"<<std::endl;
  
  
  

  return rStatus;
}
