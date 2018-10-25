#include "printPlots.h"

const bool debugMode=true;
const bool draw_hJER=true;
const bool draw_MCEff=false;
const bool drawProfiles = false;
const bool rebinJER=false;
const int rebinFactor=5;

//other options
const bool draw_hJERRapBins=false, doGenBinsToo=false;//RapBins -> dual-diff xsec bins, GenBins -> variable, depends on readForests
const bool draw_JERgen150to200=false, draw_JERgen30to50=false;
//const std::string defOutFilename="printPlots_MCJEC_ppMC_defOut";
const std::string didJetID="1";
const std::string doJetID="1";

int printPlots_ppMC_JERS(std::string inFile_MC_dir,const std::string outputTag){
  
  // root style settings.
  std::cout<<"forcing style"<<std::endl;
  LoadStyle();
  
  //figure out what radius/jetcollection we are looking at using the ppData filename
  std::size_t radPos=inFile_MC_dir.find("_ak")+strlen("_ak");
  const std::string radiusInt= inFile_MC_dir.substr( radPos,1 );
  const std::string radius="R"+radiusInt+"_";
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;
  
  //make final jetType strings
  std::size_t jetTypePos=radPos+1;
  std::size_t jetsPos=inFile_MC_dir.find("Jets");;  
  const std::string jetType=inFile_MC_dir.substr( jetTypePos,(jetsPos-jetTypePos) );
  if(debugMode)std::cout<<"jetType string is = "<<jetType<<std::endl;
  
  const std::string fullJetType="ak"+radiusInt+jetType;//"ak3PFJets";//find in dir name
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  
  // input file
  const std::string inFile_MC_name="/Py8_CUETP8M1_QCDjetAllPtBins_"+fullJetType+"-allFiles.root";  
  std::cout<<std::endl<<"input file dir "<< inputDir+inFile_MC_dir <<std::endl;
  std::cout<<"input file name in dir "<< inFile_MC_name <<std::endl;
  TFile *finPP=new TFile( (inputDir+inFile_MC_dir+inFile_MC_name).c_str() );
  
  // output file
  std::string outFileName=outputDir+fullJetType+"_MCJEC_"+outputTag+".root";
  std::cout<<std::endl<<"opening output root file "<< outFileName<<std::endl;
  TFile *rootfout = new TFile(outFileName.c_str(),"RECREATE");   
  rootfout->cd();    
  
  // fit + fit settings
  std::cout<<"setting up fitters..."<< std::endl;
  TVirtualFitter::SetDefaultFitter("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04); 
  
  TH1F *hrsp[nbins_pt_debug]={};  //input
  
  TH1F *hMean=NULL, *hSigma=NULL; //output    
  TH1F *hChi2NDF=NULL;
  Double_t fit_ptLo=-1.;//, fit_ptHi;
  TF1 *fgaus=NULL;     //for later

  if(!draw_hJER)
    std::cout<<std::endl<<std::endl<<"skipping hJER hists"<<std::endl<<std::endl;
  else {
    
    const std::string algname=algo+radiusInt+jetType; 
    std::cout<<"making hists for algorithm "<< algname<<"..."<<std::endl;
    
    // initialize mean hist    
    hMean = new TH1F( "hMean",    
		      Form( "Mean %s", algname.c_str()),    
		      nbins_pt_debug, ptbins_debug);
    
    hSigma = new TH1F( "hSigma",    
		       Form( "Sigma %s", algname.c_str()),    
		       nbins_pt_debug, ptbins_debug);
    
    hChi2NDF = new TH1F("hChi2NDF",    
			Form( "Chi2NDF %s", algname.c_str()),    
			nbins_pt_debug, ptbins_debug); 
    

    
    // pt bin hJER-fit loop    
    if(debugMode)std::cout<<"nbins_pt="<<nbins_pt_debug<<std::endl;    
    bool ptLoSet=false;
    for(int ip=0; ip<nbins_pt_debug; ip++){    
      
      //// reco/gen pt in gen pt bins for fits
      int ptbin_ip=(int)ptbins_debug[ip], ptbin_ip1=(int)ptbins_debug[ip+1];
      if(debugMode)std::cout<<"pt range for bin: "<<ptbin_ip<<" - "<<ptbin_ip1<< " GeV "<<std::endl;    


      // open the input hist    
      std::string inputHistName="hJER_"+doJetID+"wJetID_ptbin"+std::to_string(ip);
      hrsp[ip] = (TH1F*)finPP->Get( inputHistName.c_str() );        
      
      if(rebinJER)
	hrsp[ip]= (TH1F*)hrsp[ip]->TH1::Rebin( rebinFactor, (inputHistName+"_rebin5").c_str());	
      
      if(!hrsp[ip]){ 
	std::cout<<"no input hist named " <<  inputHistName<< ", exiting..."<<std::endl;
	assert(false);}          
      else if(hrsp[ip]->GetEntries()<3)
	continue;	
      if(debugMode)std::cout<<std::endl;  
      if(debugMode)hrsp[ip]->Print("base");    
      hrsp[ip]->Scale( 1./ hrsp[ip]->Integral() );    //scale to area of 1
      
      //fgaus = new TF1("fgaus","gaus", 0.10,1.90);        //fgaus = new TF1("fgaus","gaus", 0.70,1.30);          //fgaus = new TF1("fgaus","gaus", 0.90,1.10);    
      fgaus = new TF1("fgaus","gaus", 0.80,1.20);          
      //fgaus->SetParLimits(0,0.98,1.02);    //normalization
      //fgaus->SetParLimits(1,0.95,1.05);   // mean
      //fgaus->SetParLimits(2,0.0,0.5);     // width
      
      std::cout<<"fitting inputHist "<<inputHistName<<std::endl;    
      if(debugMode)std::cout<< "Mean= "<< hrsp[ip]->GetMean()<<" +/- " <<hrsp[ip]->GetMeanError()<< std::endl;    	


      //int fitstatus = 1;    	//fit status=0 if fit successful... annoying
      //fitstatus = hrsp[ip]->Fit(fgaus,"RQ");          //fitstatus = hrsp[ip]->Fit(fgaus,"RS");          //fitstatus = hrsp[ip]->Fit(fgaus);    
      bool fitfailed=true;
      fitfailed = ((bool)(hrsp[ip]->Fit(fgaus,"R")));    
      //hrsp[ip]->Write();

      //dont use the first few low pt bins for the fit; just write them to file for record keeping
      if(ip<3)continue;      
      if(!ptLoSet){
	fit_ptLo=ptbins_debug[ip];
	ptLoSet=true;
      }
      

      if(fitfailed)
	std::cout<<std::endl<<"!!!WARNING FIT FAILED!!!"<<std::endl<<std::endl;
      else 
	std::cout<<std::endl<<"fit success"<<std::endl<<std::endl;
      
      // get some quick stats    
      //float norm  = hrsp[ip]->GetMaximumStored();    
      float mean  = hrsp[ip]->GetMean();    
      float emean = hrsp[ip]->GetMeanError();    
      float sig   = hrsp[ip]->GetStdDev()/mean;    
      float esig  = sig* 
	sqrt (  pow( emean/mean,2) + pow( hrsp[ip]->GetStdDevError()/hrsp[ip]->GetStdDev(),2) ) ; 
      float chi2NDF = -1.;
      
      if(fitfailed)

      
      mean  = (fitfailed) ? hrsp[ip]->GetMean()        : fgaus->GetParameter(1); // z = if(condition) then(?) <do this> else(:) <do this>   
      emean = (fitfailed) ? hrsp[ip]->GetMeanError()   : fgaus->GetParError(1) ; 
      sig   = (fitfailed) ? hrsp[ip]->GetStdDev()/mean : fgaus->GetParameter(2)/mean;    
      esig  = (fitfailed) ? 
	sig * sqrt( pow( (emean/mean),2) + 
		    pow( ( hrsp[ip]->GetStdDevError()/hrsp[ip]->GetStdDev() ),2) 
		    ) : 
	sig * sqrt( pow( (emean/mean), 2) + 
		    pow( ( fgaus->GetParError(2)/fgaus->GetParameter(2) ), 2) 
		    )   ;
      chi2NDF = (fitfailed) ? -1. : fgaus->GetChisquare()/((float)fgaus->GetNDF());
      
      // set contents+errors in specific bin    
      std::cout<<"setting bin content for hMean[ip]=["<< ip<<"] "<<std::endl;    
      if(debugMode)std::cout<<"mean ="<< mean<<std::endl;    
      if(debugMode)std::cout<<"emean ="<< emean<<std::endl;    
      if(debugMode)std::cout<<"sig ="<< sig<<std::endl;    
      if(debugMode)std::cout<<"esig ="<< esig<<std::endl;    
      if(debugMode)std::cout<<"chi2NDF ="<< chi2NDF<<std::endl;    
      
      if((ip+1)==nbins_pt_debug)continue;
      else if((ip+1)<=2)continue;
      hMean->SetBinContent (ip+1, mean);    
      hMean->SetBinError   (ip+1, emean);    
      
      hSigma->SetBinContent (ip+1, sig);    
      hSigma->SetBinError   (ip+1, esig);    
      
      hChi2NDF->SetBinContent (ip+1, chi2NDF );
      hChi2NDF->SetBinError (ip+1, 0. );
      
    }// end fit-loop over ptbins    
    
    //if(debugMode)std::cout<<"closing gPad..."<<std::endl<<std::endl;    
    //gPad->Close();        
    
    // for mu v. gen pt and sigma/mu v. genpt hists_fit
    std::cout<<std::endl<<"fitting JER v genpt hists ..."<<std::endl;        
    TF1 *hSigmaFit=new TF1("hSigmaFit","[0]+[1]/(pow(x,[2])+[3]*x)", 
			   fit_ptLo,
			   ptbins_debug[nbins_pt_debug-1] );       // I want this fit to end at 967 GeV (for entire range)    
    
    float fitParam_0    =hSigmaFit->GetParameter(0);
    float fitParam_0_err=hSigmaFit->GetParError(0);
    float fitParam_1    =hSigmaFit->GetParameter(1);
    float fitParam_1_err=hSigmaFit->GetParError(1);
    float fitParam_2    =hSigmaFit->GetParameter(2);
    float fitParam_2_err=hSigmaFit->GetParError(2); 
    
    std::cout << "fitParam_0     = " << fitParam_0    << std::endl;
    std::cout << "fitParam_0_err = " << fitParam_0_err<< std::endl;
    std::cout << "%err fitParam0 = " << (fitParam_0_err/fitParam_0)*100.    << std::endl;
    std::cout << "fitParam_1     = " << fitParam_1    << std::endl;
    std::cout << "fitParam_1_err = " << fitParam_1_err<< std::endl;
    std::cout << "%err fitParam1 = " << (fitParam_1_err/fitParam_1)*100.    << std::endl;
    std::cout << "fitParam_2     = " << fitParam_2    << std::endl;
    std::cout << "fitParam_2_err =   " << fitParam_2_err<< std::endl;
    std::cout << "%err fitParam2 = " << (fitParam_2_err/fitParam_2)*100.    << std::endl;
    
    
    bool sigmaFitStatus= true;
    std::cout<<std::endl<<"fitting JER width v. pt... init. sigmaFitStatus="<<sigmaFitStatus<<std::endl<<std::endl; 
    //sigmaFitStatus=hSigma->Fit(hSigmaFit,"MEVR");    //sigmaFitStatus=hSigma->Fit(hSigmaFit,"MR");
    sigmaFitStatus=hSigma->Fit(hSigmaFit,"R");	
    
    fitParam_0    =hSigmaFit->GetParameter(0);
    fitParam_0_err=hSigmaFit->GetParError(0);
    fitParam_1    =hSigmaFit->GetParameter(1);
    fitParam_1_err=hSigmaFit->GetParError(1);
    fitParam_2    =hSigmaFit->GetParameter(2);
    fitParam_2_err=hSigmaFit->GetParError(2); 
    
    if(!sigmaFitStatus){	
      std::cout<<std::endl<<"fit success!, sigmaFitStatus="<<sigmaFitStatus<<std::endl<<std::endl; 
      std::cout << "fitParam_0     = " << fitParam_0    << std::endl;
      std::cout << "fitParam_0_err = " << fitParam_0_err<< std::endl;
      std::cout << "%err fitParam0 = " << (fitParam_0_err/fitParam_0)*100.    << std::endl;
      std::cout << "fitParam_1     = " << fitParam_1    << std::endl;
      std::cout << "fitParam_1_err = " << fitParam_1_err<< std::endl;
      std::cout << "%err fitParam1 = " << (fitParam_1_err/fitParam_1)*100.    << std::endl;
      std::cout << "fitParam_2     = " << fitParam_2    << std::endl;
      std::cout << "fitParam_2_err =   " << fitParam_2_err<< std::endl;
      std::cout << "%err fitParam2 = " << (fitParam_2_err/fitParam_2)*100.    << std::endl;
      //	assert(false);
    }
    else {
      std::cout<<std::endl<<"fit failed, sigmaFitStatus="<<sigmaFitStatus<<std::endl<<std::endl; 
      //	assert(false);
    }
    







    TCanvas* pdfoutCanv_muSigma=new TCanvas("outputPdfwLog_muSigma","outputPdfwLog", 900, 700);    
    pdfoutCanv_muSigma->Divide(1,3);
    pdfoutCanv_muSigma->cd();
    
    TPad* p1=(TPad*)pdfoutCanv_muSigma->cd(1);        
    p1->SetLogx(1);    
    p1->SetGridx(0);    
    p1->cd();
    
    TPad* p2=(TPad*)pdfoutCanv_muSigma->cd(2);        
    p2->SetLogx(1);    	
    p2->SetGridx(0);    
    p2->cd();   
    
    TPad* p3=(TPad*)pdfoutCanv_muSigma->cd(3);        
    p3->SetLogx(1);    	
    p3->SetGridx(0);
    p3->cd();        
    
    
    //draw pad 1
    p1->cd();    
    MakeHistMean( (TH1F*)hMean,ptbins_debug[0],ptbins_debug[nbins_pt_debug]);     
    hMean->Draw("HIST E1"); 
    
    TLine* meanLine=new TLine(ptbins_debug[0],1.,ptbins_debug[nbins_pt_debug],1.);    
    meanLine->SetLineStyle(2);  
    meanLine->SetLineColor(kBlue);    
    meanLine->Draw();            
    
    p2->cd();    
    MakeHistRMS( (TH1F*) hSigma, ptbins_debug[0],ptbins_debug[nbins_pt_debug]);   
    hSigma->Draw("HIST E1");
    hSigmaFit->Draw("SAME");
    
    //draw pad 3
    p3->cd();    
    MakeHistChi2NDF( (TH1F*)hChi2NDF,ptbins_debug[0],ptbins_debug[nbins_pt_debug]);           
    hChi2NDF->Draw("HIST E1");
    









    // DRAW THOSE PDFS //        
    std::string thePDFFileName=outputDir+fullJetType+"_MCJEC_"+outputTag+".pdf";    
    std::string open_thePDFFileName=thePDFFileName+"[";    
    std::string close_thePDFFileName=thePDFFileName+"]";    
    
    std::cout<<"printing JER reco/gen hists"<<std::endl;
    TCanvas* pdfoutCanv_wLogy=new TCanvas("outputPdfwLogx","output Pdf wLogx", 800, 800);    
    pdfoutCanv_wLogy->SetLogy(1);
    pdfoutCanv_wLogy->cd();            
    pdfoutCanv_wLogy->Print( open_thePDFFileName.c_str() );        
    
    for(int j=0;j<nbins_pt_debug;++j){       //if(j<3)continue;	
      int ptbin_j=(int)ptbins_debug[j];
      int ptbin_j1=(int)ptbins_debug[j+1];
      std::string hrspTitle=std::to_string(ptbin_j)+" GeV < gen jet p_{T} < "+std::to_string(ptbin_j1)+" GeV";    
      std::string hrsp_XAxTitle="reco p_{T}/gen p_{T}";      
      makeJERSHists( ( TCanvas*    ) pdfoutCanv_wLogy ,  
		     ( std::string ) thePDFFileName   ,  
		     ( TH1F*       ) hrsp[j]          , 
		     ( std::string ) hrspTitle        ,
		     ( std::string ) hrsp_XAxTitle    ,
		     ( TFile*      ) rootfout           );

    }
    //assert(false);         



         
    

    
    pdfoutCanv_muSigma->Print(thePDFFileName.c_str());        

    std::cout<<"closing output pdf file "<<thePDFFileName<<std::endl<<std::endl;                
    pdfoutCanv_muSigma->Print( close_thePDFFileName.c_str() );        
    
    rootfout->cd();
    
    hMean->Write();               
    hSigma->Write();               
    hChi2NDF->Write();                     
    hSigmaFit->Write();
    pdfoutCanv_muSigma->SetTitle("Mean, Significance, #chi^{2}/n.d.o.f. Canvas");        pdfoutCanv_muSigma->Write("mean_sigma_chi2ndf_canv");    

    
    
    //pdfoutCanv_muSigma->Close();
    //pdfoutCanv_wLogy->Close();
    
    
  }
  
  
  
  
  
  
  
  /////////////////////////////////////////////  /////////////////////////////////////////////  /////////////////////////////////////////////  /////////////////////////////////////////////
  /////////////////////////////////////////////  /////////////////////////////////////////////  /////////////////////////////////////////////  /////////////////////////////////////////////
  /////////////////////////////////////////////  /////////////////////////////////////////////  /////////////////////////////////////////////  /////////////////////////////////////////////
  if(!draw_MCEff)  std::cout<<std::endl<<std::endl<<"skipping MCEff hists"<<std::endl<<std::endl;
  else{
    
    
    
    std::cout<<" drawing MC Eff. QA Plots..."<<std::endl;
    
    //std::string thePDFFileName=outputDir+fullJetType+jobType+"_MCEff_"+outputTag+".pdf";
    std::string thePDFFileName=outputDir+fullJetType+"_MCJEC_"+outputTag+"_MCEff.pdf";
    std::string open_thePDFFileName=thePDFFileName+"[";
    std::string close_thePDFFileName=thePDFFileName+"]";
    std::cout<<std::endl<<"creating temporary canvas for printing MCEff plots..."<<std::endl;
    
    TCanvas *temp_canvMCEff = new TCanvas("tempMCEff", "tempMCEff", 1200, 600);
    temp_canvMCEff->Print(open_thePDFFileName.c_str());
    temp_canvMCEff->SetLogz(1);
    temp_canvMCEff->SetLogx(1);
    temp_canvMCEff->cd();
    
    
    for(int j=0; j<(N_genVars); j++){     
      if(debugMode)std::cout<<std::endl<<"j= "<<j<<std::endl;	
      
      
      std::string inHistName="hpp_mceff_"+genVars[j];
      if(doJetID=="1")inHistName+="_wJetID";
      inHistName+="_"+radius+etaWidth;
      
      std::cout<<"opening TH2F "<<inHistName << std::endl;
      TH2F* the2DMCEffQAHist= (TH2F*)finPP->Get( inHistName.c_str() );
      if(!the2DMCEffQAHist) {std::cout<<"no MCEff plot, continuing..."<<std::endl; continue;}
      
      if(j>2){ 
	temp_canvMCEff->SetLogx(0);
	temp_canvMCEff->cd();	}
      else {
	temp_canvMCEff->SetLogx(1);
	temp_canvMCEff->cd();
	the2DMCEffQAHist->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");	
	if(genVars_xAx[j]=="rawpt")	  the2DMCEffQAHist->SetAxisRange(0.5 , 1.5,"Y");	
      }
      the2DMCEffQAHist->SetAxisRange(1.*pow(10.,-13.), 1.*pow(10.,-4.),"Z");	
      
      std::string h_Title   ="MC QA, TH2,"+genVars[j];
      if(doJetID=="1")h_Title+=", w/ JetIDCut";      
      
      //std::string h_ZAx_Title="Entries"; 
      //std::string stringRatio="(rec"+genVars[j]+"/gen"+genVars[j]+")";
      //std::string h_YAx_Title=stringRatio;
      
      //std::string h_XAx_Title="";    
      //if(j==0)h_XAx_Title="gen Jet Pt (GeV)"  ;
      //if(j==1)h_XAx_Title="gen Jet Eta"       ;
      //if(j==2)h_XAx_Title="gen Jet Phi (rad)" ;
      
      
      the2DMCEffQAHist->SetTitle (    h_Title.c_str() );
      the2DMCEffQAHist->SetXTitle( genVars_xAx[j].c_str() );
      the2DMCEffQAHist->SetYTitle( genVars_yAx[j].c_str() );
      
      //the2DMCEffQAHist->SetAxisRange(0.5,2.0,"Y");
      //if(j==0)the2DMCEffQAHist->SetAxisRange(  0., 200., "X");//genpt
      //if(j==1)the2DMCEffQAHist->SetAxisRange(-3.0, 3.0 , "X");//eta
      //if(j==2)the2DMCEffQAHist->SetAxisRange(-4.0, 4.0 , "X");//phi
      
      the2DMCEffQAHist->Draw("COLZ"); 
      
      temp_canvMCEff->Print(thePDFFileName.c_str());   
      
      if(drawProfiles){     
	
	//std::string h_ZAx_Title="Entries"; 
	//std::string stringRatio="(rec"+genVars[j]+"/gen"+genVars[j]+")";
	//std::string h_YAx_Title="avg "+stringRatio+"/bin";
	
	//std::string h_XAx_Title="";    
	//if(j==0)h_XAx_Title="gen Jet Pt (GeV)"  ;
	//if(j==1)h_XAx_Title="gen Jet Eta"       ;
	//if(j==2)h_XAx_Title="gen Jet Phi (rad)" ;
	
	std::cout<<"taking the profile of "<<inHistName << std::endl;
	inHistName+="_TH2Prof";
	TH1F* theMCEffQAHist= (TH1F*)the2DMCEffQAHist->TH2::ProfileX( inHistName.c_str(), 1, -1, "o");
	if(  genVars[j]=="pt" ||
	       genVars[j]=="pt2" ||
	     genVars[j]=="pt3" )  {
	  theMCEffQAHist->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");	
	  theMCEffQAHist->SetAxisRange( 0.90, 1.10,"Y");	 
	}
	else if (genVars[j]=="ptrat_phi" || 
		 genVars[j]=="ptrat_eta" 	
		 ) {
	  theMCEffQAHist->SetAxisRange( 0.90, 1.10,"Y");	 }
	
	std::cout<<std::endl<<"rebinning TH2 Profile"<<std::endl;
	
	//inHistName+="_rebin2";
	//theMCEffQAHist=(TH1F*)theMCEffQAHist->Rebin(2,  inHistName.c_str());
	//theMCEffQAHist->Scale(0.5);
	
	h_Title+=" TH2Profile";
	theMCEffQAHist->SetTitle (    h_Title.c_str() );
	theMCEffQAHist->SetXTitle( genVars_xAx[j].c_str() );
	theMCEffQAHist->SetYTitle( ("bin-avg of "+genVars_yAx[j]).c_str() );
	
	//theMCEffQAHist->SetXTitle( h_XAx_Title.c_str() );
	//theMCEffQAHist->SetYTitle( h_YAx_Title.c_str() );
	
	//theMCEffQAHist->SetAxisRange(0.5,2.0,"Y");
	//if(j==0)theMCEffQAHist->SetAxisRange(  0., 200., "X");//genpt
	//if(j==1)theMCEffQAHist->SetAxisRange(-3.0, 3.0 , "X");//eta
	//if(j==2)theMCEffQAHist->SetAxisRange(-4.0, 4.0 , "X");//phi
	
	
	theMCEffQAHist->Draw("HIST E"); 
	
	temp_canvMCEff->Print(thePDFFileName.c_str());   
      }
      
      
    }// gen genvar loop
    
    temp_canvMCEff->Print(close_thePDFFileName.c_str());   
    temp_canvMCEff->Close();
  } // end  MCEff Plots
  
  

  if(debugMode)std::cout<<"closing output root file"<<std::endl<<std::endl;
  rootfout->Close();
  return 0;
}





int main(int argc, char*argv[]){
  int rStatus=-1;
  if(argc!=3){
    std::cout<<"no defaults. Do...."<<std::endl;
    std::cout<<"./printPlots_ppMC_JERS.exe <target_ppMC_dir> <output tag>"<<std::endl<<std::endl;
    return rStatus;  }
  
  rStatus=1;
  rStatus=printPlots_ppMC_JERS((const std::string)argv[1], (const std::string)argv[2] );
  return rStatus;
}

    //      TF1 *hSigmaFit=new TF1("hSigmaFit","[0]+[1]*pow(x,[2])", 
    //			     ptbins_debug[2],                     // I want this fit to start at 56 GeV or so
    //			     ptbins_debug[nbins_pt_debug]);       // I want this fit to end at 967 GeV (for entire range)



    //TH1F *hrawrsp[nbins_pt_debug]={};  //input
    //TH1F *hcorr[nbins_pt_debug]={};  //input




      /*
      // open + grab raw/gen hists
      std::string inputHistNameRaw="hJER_raw_"+doJetID+"wJetID_ptbin"+std::to_string(ip);
      hrawrsp[ip] = (TH1F*)finPP->Get( inputHistNameRaw.c_str() );        
      if(!hrawrsp[ip])
      std::cout<<"WARNING, no input hist named " <<  inputHistNameRaw<< ", moving on..."<<std::endl;	
      if(debugMode)std::cout<<std::endl;  
      if(debugMode)hrawrsp[ip]->Print("base");    	
      hrawrsp[ip]->Scale( 1./ hrawrsp[ip]->Integral() );    
      
      
      // open + grab reco/raw 
      std::string inputHistNameJES="hJES_"+doJetID+"wJetID_ptbin"+std::to_string(ip);
      hcorr[ip] = (TH1F*)finPP->Get( inputHistNameJES.c_str() );        
      if(!hcorr[ip])
      std::cout<<"WARNING, no input hist named " <<  inputHistNameJES<< ", moving on..."<<std::endl;
      if(debugMode)std::cout<<std::endl;  
      if(debugMode)hcorr[ip]->Print("base");    	
      hcorr[ip]->Scale( 1./ hcorr[ip]->Integral() );    
      */



      /*
      std::cout<<"printing JER raw/gen hists"<<std::endl;
      for(int j=0;j<nbins_pt_debug;++j){    
	
	int ptbin_j=(int)ptbins_debug[j];
	int ptbin_j1=(int)ptbins_debug[j+1];
	std::string hrspTitle=std::to_string(ptbin_j)+" GeV < gen jet p_{T} < "+std::to_string(ptbin_j1)+" GeV";    
	std::string hrsp_XAxTitle="raw p_{T}/gen p_{T}";	
	makeJERSHists(  ( TCanvas*    ) pdfoutCanv_wLogy ,  
		        ( std::string ) thePDFFileName	 ,  
		        ( TH1F*       ) hrawrsp[j]    , 
 		        ( std::string ) hrspTitle        ,
	 	        ( std::string ) hrsp_XAxTitle      );      }          
      std::cout<<"printing JES reco/gen hists"<<std::endl;
      for(int j=0;j<nbins_pt_debug;++j){    
	
	int ptbin_j=(int)ptbins_debug[j];
	int ptbin_j1=(int)ptbins_debug[j+1];
	std::string hrspTitle=std::to_string(ptbin_j)+" GeV < gen jet p_{T} < "+std::to_string(ptbin_j1)+" GeV";    
	std::string hrsp_XAxTitle="reco p_{T}/raw p_{T}";
	makeJERSHists(  ( TCanvas*    ) pdfoutCanv_wLogy ,  
		        ( std::string ) thePDFFileName	 ,  
		        ( TH1F*       ) hcorr[j]      , 
 		        ( std::string ) hrspTitle        ,
	 	        ( std::string ) hrsp_XAxTitle      );      }    
      */      
