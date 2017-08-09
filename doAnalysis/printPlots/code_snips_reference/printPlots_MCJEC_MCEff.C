//#include "plotNfit.h"
#include "printPlots.h"
#include "printPlots_fitsNfuncs.h"
#include "printPlots_binsNconsts.h"

const bool drawProfiles=true;
const bool drawTH2s=true;
const bool debugMode=true;
const bool doJetIDPlots=true;
//const std::string defOutFilename="printPlots_jtCalibHists_ppMC_defOut";

int printPlots_jtCalibHists(const std::string inFile_MC_dir,const std::string baseName){

  // root style settings.
  std::cout<<std::endl<<"loading style..." <<std::endl;
  LoadStyle();  
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);  
  
  
  // input stuff. //

  //figure out what radius/jetcollection we are looking at using the ppData filename
  std::size_t radPos=inFile_MC_dir.find("_ak")+3;
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
  //std::string doJetID=inFile_MC_dir.substr( inFile_MC_dir.find("MCJEC_JetID")+11 );
//if(debugMode)std::cout<<"doJetID="<<doJetID<<std::endl;
  std::string doJetID="1";

  // input file
  const std::string inFile_MC_name="/Py8_CUETP8M1_QCDjetAllPtBins_"+fullJetType+"-allFiles.root";
  std::cout<<std::endl<<"opening input file in dir "<< inFile_MC_dir <<std::endl;
  std::cout<<"input file name in dir "<< inFile_MC_name <<std::endl;
  TFile *finPP=new TFile( (inputDir+inFile_MC_dir+inFile_MC_name).c_str() );
 
  std::string jobType="_MCJEC_JetID"+doJetID;
  if(debugMode)std::cout<<"jobType="<<jobType<<std::endl;

  { std::cout<<" drawing MC Eff. QA Plots..."<<std::endl;
    
    std::string thePDFFileName=outputDir+fullJetType+jobType+"_MCEff_"+baseName+".pdf";
    std::string open_thePDFFileName=thePDFFileName+"[";
    std::string close_thePDFFileName=thePDFFileName+"]";
    std::cout<<std::endl<<"creating temporary canvas for printing MCEff plots..."<<std::endl;

    TCanvas *temp_canvMCEff = new TCanvas("tempMCEff", "tempMCEff", 1200, 600);
    temp_canvMCEff->Print(open_thePDFFileName.c_str());
    temp_canvMCEff->cd();
  

    if(drawProfiles){     
      for(int j=0; j<(N_genVars); j++){     
  	if(debugMode)std::cout<<std::endl<<"j= "<<j<<std::endl;
  
  	  
  	std::string inHistName="hpp_mceff_"+genVars[j];
  	if(doJetIDPlots)inHistName+="_wJetID";
  	inHistName+="_"+radius+etaWidth;
  	  
  	std::cout<<"opening TH2F "<<inHistName << std::endl;
  	TH2F* the2DMCEffQAHist= (TH2F*)finPP->Get( inHistName.c_str() );
  	if(!the2DMCEffQAHist) {std::cout<<"no MCEff plot, continuing..."<<std::endl; continue;}
  	the2DMCEffQAHist->Print("base");
  	  
  	std::string h_Title   ="MC Eff. QA, TH2 profile";
  	if(doJetIDPlots)h_Title+=", w/ JetIDCut";      
  	  
  	//std::string h_ZAx_Title="Entries"; 
  	std::string stringRatio="(rec"+genVars[j]+"/gen"+genVars[j]+")";
  	std::string h_YAx_Title="avg "+stringRatio+"/bin";
  	  
  	std::string h_XAx_Title="";	    
  	if(j==0)h_XAx_Title="gen Jet Pt (GeV)"  ;
  	if(j==1)h_XAx_Title="gen Jet Eta"       ;
  	if(j==2)h_XAx_Title="gen Jet Phi (rad)" ;
  	  
  	std::cout<<"taking the profile of "<<inHistName << std::endl;
  	TH1F* theMCEffQAHist= (TH1F*)the2DMCEffQAHist->TH2::ProfileX("MCEffProfile",1,-1, "o");
  	  
  	theMCEffQAHist->SetTitle (    h_Title.c_str() );
  	theMCEffQAHist->SetXTitle( h_XAx_Title.c_str() );
  	theMCEffQAHist->SetYTitle( h_YAx_Title.c_str() );
  	  
  	theMCEffQAHist->SetAxisRange(0.5,2.0,"Y");
  	if(j==0)theMCEffQAHist->SetAxisRange(  0., 200., "X");//genpt
  	if(j==1)theMCEffQAHist->SetAxisRange(-3.0, 3.0 , "X");//eta
  	if(j==2)theMCEffQAHist->SetAxisRange(-4.0, 4.0 , "X");//phi
  	  
  	theMCEffQAHist->Draw(); 
  	 
  	temp_canvMCEff->Print(thePDFFileName.c_str());   
      }// gen genvar loop
      
      for(int j=0; j<(N_genVars_ptrat); j++){     
  	if(debugMode)std::cout<<std::endl<<"j= "<<j<<std::endl;      
  	
  	std::string inHistName="hpp_mceff_ptrat_"+genVars_ptrat[j];
  	if(doJetIDPlots)inHistName+="_wJetID";
  	inHistName+="_"+radius+etaWidth;
  	  
  	TH2F* the2DMCEffQAHist= (TH2F*)finPP->Get( inHistName.c_str() );
  	
  	std::string h_Title   ="MC Eff. QA, TH2 profile";
  	if(doJetIDPlots)h_Title+=", w/ JetIDCut";      
  	    
  	std::string h_YAx_Title="avg (recpt)/(genpt)/bin";	    
  	std::string h_XAx_Title="";	    
  	if(j==0)h_XAx_Title="gen Jet Eta (GeV)"  ;
  	if(j==1)h_XAx_Title="gen Jet Phi (rad)"  ;
  	if(j==2)h_XAx_Title="gen Jet dRJet"      ; 
  	
  	TH1F* theMCEffQAHist= (TH1F*)the2DMCEffQAHist->TH2::ProfileX("MCEffProfile",1,-1, "o");
  	
  	theMCEffQAHist->SetAxisRange(0.5,2.0,"Y");
  	if(j==0)theMCEffQAHist->SetAxisRange( -3.0, 3.0, "X");
  	if(j==1)theMCEffQAHist->SetAxisRange( -4.0, 4.0, "X");
  	if(j==2)theMCEffQAHist->SetAxisRange(   0., 0.5, "X");
  	
  	theMCEffQAHist->SetTitle (    h_Title.c_str() );
  	theMCEffQAHist->SetXTitle( h_XAx_Title.c_str() );
  	theMCEffQAHist->SetYTitle( h_YAx_Title.c_str() );
  	
  	theMCEffQAHist->Draw(); 
  	
  	temp_canvMCEff->Print(thePDFFileName.c_str());   	
      }//end genvar ptrat loops      
    }//end drawProfles



    if(drawTH2s){     
      for(int j=0; j<(N_genVars); j++){     
  	if(debugMode)std::cout<<std::endl<<"j= "<<j<<std::endl;
  
  	  
  	std::string inHistName="hpp_mceff_"+genVars[j];
  	if(doJetIDPlots)inHistName+="_wJetID";
  	inHistName+="_"+radius+etaWidth;
  	  
  	std::cout<<"opening TH2F "<<inHistName << std::endl;
  	TH2F* the2DMCEffQAHist= (TH2F*)finPP->Get( inHistName.c_str() );
  	if(!the2DMCEffQAHist) {std::cout<<"no MCEff plot, continuing..."<<std::endl; continue;}
  	  
  	std::string h_Title   ="MC Eff. QA, TH2";
  	if(doJetIDPlots)h_Title+=", w/ JetIDCut";      
  	  
  	//std::string h_ZAx_Title="Entries"; 
  	std::string stringRatio="(rec"+genVars[j]+"/gen"+genVars[j]+")";
  	std::string h_YAx_Title=stringRatio;
  	  
  	std::string h_XAx_Title="";	    
  	if(j==0)h_XAx_Title="gen Jet Pt (GeV)"  ;
  	if(j==1)h_XAx_Title="gen Jet Eta"       ;
  	if(j==2)h_XAx_Title="gen Jet Phi (rad)" ;
  	  
  	  
  	the2DMCEffQAHist->SetTitle (    h_Title.c_str() );
  	the2DMCEffQAHist->SetXTitle( h_XAx_Title.c_str() );
  	the2DMCEffQAHist->SetYTitle( h_YAx_Title.c_str() );
  	  
  	the2DMCEffQAHist->SetAxisRange(0.5,2.0,"Y");
  	if(j==0)the2DMCEffQAHist->SetAxisRange(  0., 200., "X");//genpt
  	if(j==1)the2DMCEffQAHist->SetAxisRange(-3.0, 3.0 , "X");//eta
  	if(j==2)the2DMCEffQAHist->SetAxisRange(-4.0, 4.0 , "X");//phi
  	  
  	the2DMCEffQAHist->Draw("COLZ"); 
  	 
  	temp_canvMCEff->Print(thePDFFileName.c_str());   
      }// gen genvar loop
      
      for(int j=0; j<(N_genVars_ptrat); j++){     
  	if(debugMode)std::cout<<std::endl<<"j= "<<j<<std::endl;      
  	
  	std::string inHistName="hpp_mceff_ptrat_"+genVars_ptrat[j];
  	if(doJetIDPlots)inHistName+="_wJetID";
  	inHistName+="_"+radius+etaWidth;
  	  
  	TH2F* the2DMCEffQAHist= (TH2F*)finPP->Get( inHistName.c_str() );
  	
  	std::string h_Title   ="MC Eff. QA, TH2 profile";
  	if(doJetIDPlots)h_Title+=", w/ JetIDCut";      
  	    
  	std::string h_YAx_Title="(recpt)/(genpt)";	    
  	std::string h_XAx_Title="";	    
  	if(j==0)h_XAx_Title="gen Jet Eta (GeV)"  ;
  	if(j==1)h_XAx_Title="gen Jet Phi (rad)"  ;
  	if(j==2)h_XAx_Title="gen Jet dRJet"      ; 
  	
  	the2DMCEffQAHist->SetAxisRange(0.5,2.0,"Y");
  	if(j==0)the2DMCEffQAHist->SetAxisRange( -3.0, 3.0, "X");
  	if(j==1)the2DMCEffQAHist->SetAxisRange( -4.0, 4.0, "X");
  	if(j==2)the2DMCEffQAHist->SetAxisRange(   0., 0.5, "X");
  	
  	the2DMCEffQAHist->SetTitle (    h_Title.c_str() );
  	the2DMCEffQAHist->SetXTitle( h_XAx_Title.c_str() );
  	the2DMCEffQAHist->SetYTitle( h_YAx_Title.c_str() );
  	
  	the2DMCEffQAHist->Draw("COLZ"); 
  	
  	temp_canvMCEff->Print(thePDFFileName.c_str());   	
      }//end genvar ptrat loops      
    }//end draw TH2s
    temp_canvMCEff->Print(close_thePDFFileName.c_str());   	
    temp_canvMCEff->Close();        
  }//end  MCEff Plots else




  

  //// save output root file. //
  //
  //{
  //
  //  std::string outFileName=outputDir+fullJetType+jobType+"_"+baseName+".root";
  //  std::cout<<std::endl<<"opening output root file "<< outFileName<<std::endl;
  //  TFile *rootfout = new TFile(outFileName.c_str(),"RECREATE"); 
  //  rootfout->cd();  
  //  std::cout<<"writing output to root file..."<<std::endl;
  //  for(int i=0;i<Nrad;++i){
  //    hMean[i]->Write();
  //    hSigma[i]->Write();    
  //    for(int j=0;j<nbins_eta;++j) {
  //	if( !(hrsp[i][j]) ) continue;
  //	else hrsp[i][j]->Write();	
  //    }
  //  }
  //  std::cout<<"closing output root file"<<std::endl<<std::endl;
  //  rootfout->Close();
  //}
  
  return 0;
}


int main(int argc, char*argv[]){
  int rStatus=-1;
  if(argc!=3){
    std::cout<<"no defaults. Do...."<<std::endl;
    std::cout<<"./printPlots_jtCalibHists.exe "<<
      "<target_ppMC_dir> <outputNameBase>"<<std::endl<<std::endl;
    return rStatus;  }
  
  rStatus=1;
  rStatus=printPlots_jtCalibHists((const std::string)argv[1], (const std::string)argv[2] );
  return rStatus;
}


//// ALT METRICS
// mean  = (fitstatus!=0) ? hrsp[nj][ip]->GetMean()     :  0.5*(fgaus->GetParameter(1) + hrsp[nj][ip]->GetMean());
// mean  = (fitstatus!=0) ? hrsp[nj][ip]->GetMean()     : fgaus->GetParameter(1);
// emean = (fitstatus!=0) ? hrsp[nj][ip]->GetMeanError(): hrsp[nj][ip]->GetMeanError(1);
// emean = (fitstatus!=0) ? hrsp[nj][ip]->GetMeanError(): sqrt(pow(fgaus->GetParError(1),2) + pow(hrsp[nj][ip]->GetMeanError(),2));      
// sig   = (fitstatus!=0) ? hrsp[nj][ip]->GetRMS()/mean : fgaus->GetParameter(2)/fgaus->GetParameter(1);
// esig  = (fitstatus!=0) 
// ? sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))) 
// : sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))); 
// esig  = (fitstatus!=0) 
// ? sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))) 
// : sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))); 
//       
      
//hrsp[nj][ip]->Write();
