//#include "plotNfit.h"
#include "printPlots.h"
#include "printPlots_fitsNfuncs.h"
#include "printPlots_binsNconsts.h"

const bool debugMode=true;
//const std::string defOutFilename="printPlots_MCJEC_ppMC_defOut";

//int printPlots_MCJEC(const std::string inFile_MC_dir,const std::string baseName){
int printPlots_MCJEC(){

  // root style settings.
  std::cout<<std::endl<<"loading style..." <<std::endl;
  LoadStyle();  
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);  
  
  
  //figure out what radius/jetcollection we are looking at using the ppData filename
  std::size_t radPos=inFile_MC_dir.find("_ak")+strlen("_ak");
  const std::string radiusInt= inFile_MC_dir.substr( radPos,1 );
  const std::string radius="R"+radiusInt+"_";
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;
  std::size_t jetTypePos=radPos+1;
  std::size_t jetsPos=inFile_MC_dir.find("Jets");;
  
  //make final jetType strings
  const std::string jetType=inFile_MC_dir.substr( jetTypePos,(jetsPos-jetTypePos) );
  if(debugMode)std::cout<<"jetType string is = "<<jetType<<std::endl;
  const std::string fullJetType="ak"+radiusInt+jetType;//"ak3PFJets";//find in dir name
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  std::string didJetID=inFile_MC_dir.substr( inFile_MC_dir.find("MCJEC_jtID")+strlen("MCJEC_jtID"), 1 );
  //if(debugMode)std::cout<<"doJetID="<<doJetID<<std::endl;
  
  std::string doJetID;
  if(didJetID=="0") doJetID=didJetID;
  else if (didJetID=="L"||didJetID=="T") doJetID="1";
  
  std::string jobType="_MCJEC_jtID"+didJetID;
  if(debugMode)std::cout<<"jobType="<<jobType<<std::endl;

  // input file
  const std::string inFile_MC_name="/Py8_CUETP8M1_QCDjetAllPtBins_"+fullJetType+"-allFiles.root";
  std::cout<<std::endl<<"opening input file in dir "<< inFile_MC_dir <<std::endl;
  std::cout<<"input file name in dir "<< inFile_MC_name <<std::endl;
  TFile *finPP=new TFile( (inputDir+inFile_MC_dir+inFile_MC_name).c_str() );
  
  std::cout<<"closing output root file"<<std::endl<<std::endl;
  rootfout->Close();
  return 0;
}


int main(int argc, char*argv[]){
  //int rStatus=-1;
  //if(argc!=3){
  //  std::cout<<"no defaults. Do...."<<std::endl;
  //  std::cout<<"./printPlots_MCJEC.exe "<<
  //    "<target_ppMC_dir> <outputNameBase>"<<std::endl<<std::endl;
  //  return rStatus;  }
  //
  //rStatus=1;
  rStatus=printPlots_MCJEC((const std::string)argv[1], (const std::string)argv[2] );
  rStatus=printPlots_MCJEC();
  return rStatus;
}
