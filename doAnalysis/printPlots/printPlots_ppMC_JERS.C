#include "printPlots.h"

const bool debugMode=true;
const bool draw_hJER=true;
const bool draw_MCEff=true;
const bool drawProfiles = true;
const bool rebinJER=false;

//other options
const bool draw_hJERRapBins=false, doGenBinsToo=false;//RapBins -> dual-diff xsec bins, GenBins -> variable, depends on readForests
const bool draw_JERgen150to200=false, draw_JERgen30to50=false;
//const std::string defOutFilename="printPlots_MCJEC_ppMC_defOut";

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
  
  std::string didJetID="1";
  std::string doJetID="1";
  
  
  
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
  TF1 *fgaus=NULL;   
  
  // input hists
  if(draw_hJER){
    
    // array of pointers for input/output hists
    TH1F *hrsp[Nrad][nbins_pt_debug]={};  //input
    //TH1F *hrawrsp[Nrad][nbins_pt_debug]={};  //input
    //TH1F *hcorr[Nrad][nbins_pt_debug]={};  //input
    
    TH1F *hMean[Nrad]={}, *hSigma[Nrad]={}; //output    
    TH1F *hChi2NDF[Nrad]={};
    
    // loop over available radii (one radii for now)
    for(int nj=0; nj<Nrad; nj++){
      
      const std::string algname=algo+srad[nj]+jetType; 
      std::cout<<"making hists for algorithm "<< algname<<"..."<<std::endl;

      // initialize mean hist    
      hMean[nj] = new TH1F(Form( "hMean%d", nj),    
			   Form( "Mean %s %s", algname.c_str(), ccent[nj]),    
			   nbins_pt_debug, ptbins_debug); //nbins_pt, ptbins_bound);    
      MakeHistMean( (TH1F*)hMean[nj],ptbins_debug[0],ptbins_debug[nbins_pt_debug]);     
      
      // initialize RMS/Sigma hist    
      hSigma[nj] = new TH1F(Form( "hSigma%d", nj),    
			    Form( "Sigma %s %s", algname.c_str(), ccent[nj]),    
			    nbins_pt_debug, ptbins_debug); //nbins_pt, ptbins_bound);    
      MakeHistRMS( (TH1F*) hSigma[nj], ptbins_debug[0],ptbins_debug[nbins_pt_debug]);   

      hChi2NDF[nj] = new TH1F(Form( "hChi2NDF%d", nj),    
			      Form( "Chi2NDF %s %s", algname.c_str(), ccent[nj]),    
			   nbins_pt_debug, ptbins_debug); //nbins_pt, ptbins_bound);    
      MakeHistChi2NDF( (TH1F*)hChi2NDF[nj],ptbins_debug[0],ptbins_debug[nbins_pt_debug]);       
      
      // pt bin hJER-fit loop    
      if(debugMode)std::cout<<"nbins_pt="<<nbins_pt_debug<<std::endl;    
      for(int ip=0; ip<nbins_pt_debug; ip++){    
	
	int ptbin_ip=(int)ptbins_debug[ip], ptbin_ip1=(int)ptbins_debug[ip+1];
	if(debugMode)std::cout<<"pt range for bin: "<<ptbin_ip<<" - "<<ptbin_ip1<< " GeV "<<std::endl;    


	/*
	// open + grab raw/gen hists
	std::string inputHistNameRaw="hJER_raw_"+doJetID+"wJetID_ptbin"+std::to_string(ip);
	hrawrsp[nj][ip] = (TH1F*)finPP->Get( inputHistNameRaw.c_str() );        
	if(!hrawrsp[nj][ip])
	  std::cout<<"WARNING, no input hist named " <<  inputHistNameRaw<< ", moving on..."<<std::endl;	
	if(debugMode)std::cout<<std::endl;  
	if(debugMode)hrawrsp[nj][ip]->Print("base");    	
	hrawrsp[nj][ip]->Scale( 1./ hrawrsp[nj][ip]->Integral() );    
	
	
	// open + grab reco/raw 
	std::string inputHistNameJES="hJES_"+doJetID+"wJetID_ptbin"+std::to_string(ip);
	hcorr[nj][ip] = (TH1F*)finPP->Get( inputHistNameJES.c_str() );        
	if(!hcorr[nj][ip])
	  std::cout<<"WARNING, no input hist named " <<  inputHistNameJES<< ", moving on..."<<std::endl;
	if(debugMode)std::cout<<std::endl;  
	if(debugMode)hcorr[nj][ip]->Print("base");    	
	hcorr[nj][ip]->Scale( 1./ hcorr[nj][ip]->Integral() );    
	*/

	//// reco/gen pt in gen pt bins for fits
	// open the input hist    
	std::string inputHistName="hJER_"+doJetID+"wJetID_ptbin"+std::to_string(ip);
	hrsp[nj][ip] = (TH1F*)finPP->Get( inputHistName.c_str() );        
	if(rebinJER){	  hrsp[nj][ip]= (TH1F*)hrsp[nj][ip]->TH1::Rebin(5, (inputHistName+"_rebin5").c_str());	}
	
	if(!hrsp[nj][ip]){ 
	  std::cout<<"no input hist named " <<  inputHistName<< ", exiting..."<<std::endl;
	  assert(false);}          
	else if(hrsp[nj][ip]->GetEntries()<3)
	  continue;	
	if(debugMode)std::cout<<std::endl;  
	if(debugMode)hrsp[nj][ip]->Print("base");    
	hrsp[nj][ip]->Scale( 1./ hrsp[nj][ip]->Integral() );    
	
	// get some quick stats    
	//float norm  = hrsp[nj][ip]->GetMaximumStored();    
	float mean  = hrsp[nj][ip]->GetMean();    
	float emean = hrsp[nj][ip]->GetMeanError();    
	float sig   = hrsp[nj][ip]->GetStdDev()/mean;    
	float esig  = sig* 
	  sqrt (  pow( emean/mean,2) + pow( hrsp[nj][ip]->GetStdDevError()/hrsp[nj][ip]->GetStdDev(),2) ) ; 
	float chi2NDF = -1.;
	

	//fgaus = new TF1("fgaus","gaus", 0.10,1.90);    
	//fgaus = new TF1("fgaus","gaus", 0.70,1.30);    
	fgaus = new TF1("fgaus","gaus", 0.80,1.20);    
	//fgaus = new TF1("fgaus","gaus", 0.90,1.10);    
	//fgaus->SetParLimits(0,0.98,1.02);    //normalization
	//fgaus->SetParLimits(1,0.95,1.05);   // mean
	//fgaus->SetParLimits(2,0.0,0.5);     // width
	
	std::cout<<"fitting..."<<std::endl;    
	std::cout<<"inputHist "<<inputHistName<<std::endl;    
	if(debugMode)std::cout<< "Mean= "<< hrsp[nj][ip]->GetMean()<<" +/- " <<hrsp[nj][ip]->GetMeanError()<< std::endl;    	
	
	int fitstatus = 1;    	//fit status=0 if fit successful... annoying
	//fitstatus = hrsp[nj][ip]->Fit(fgaus,"RQ");    
	//fitstatus = hrsp[nj][ip]->Fit(fgaus,"RS");    
	fitstatus = hrsp[nj][ip]->Fit(fgaus,"R");    
	//fitstatus = hrsp[nj][ip]->Fit(fgaus);    
	std::cout<< "Fit Status: "<< fitstatus<< std::endl;    
	
	mean  = (fitstatus!=0) ? hrsp[nj][ip]->GetMean()        : fgaus->GetParameter(1); // z = if(condition) then(?) <do this> else(:) <do this>   
	emean = (fitstatus!=0) ? hrsp[nj][ip]->GetMeanError()   : fgaus->GetParError(1) ; 
	sig   = (fitstatus!=0) ? hrsp[nj][ip]->GetStdDev()/mean : fgaus->GetParameter(2)/mean;    
	esig  = (fitstatus!=0) 
	  ? sig * sqrt( pow( (emean/mean),2) + 
			pow( ( hrsp[nj][ip]->GetStdDevError()/hrsp[nj][ip]->GetStdDev() ),2) 
			) 
	  : sig * sqrt( pow( (emean/mean), 2) + 
			pow( ( fgaus->GetParError(2)/fgaus->GetParameter(2) ), 2) 
			)   ;
	chi2NDF = (fitstatus!=0) ? -1. : fgaus->GetChisquare()/((float)fgaus->GetNDF());

	// set contents+errors in specific bin    
	std::cout<<"setting bin content for hMean[nj][ip]=["<< nj<<"]["<< ip<<"] "<<std::endl;    
	if(debugMode)std::cout<<"mean ="<< mean<<std::endl;    
	if(debugMode)std::cout<<"emean ="<< emean<<std::endl;    
	if(debugMode)std::cout<<"sig ="<< sig<<std::endl;    
	if(debugMode)std::cout<<"esig ="<< esig<<std::endl;    
	if(debugMode)std::cout<<"chi2NDF ="<< chi2NDF<<std::endl;    
	
	hMean[nj]->SetBinContent (ip+1, mean);    
	hMean[nj]->SetBinError   (ip+1, emean);    

	hSigma[nj]->SetBinContent (ip+1, sig);    
	hSigma[nj]->SetBinError   (ip+1, esig);    

	hChi2NDF[nj]->SetBinContent (ip+1, chi2NDF );
	hChi2NDF[nj]->SetBinError (ip+1, 0. );

	hrsp[nj][ip]->Write();
	
      }// end fit-loop over ptbins    
     
    }// end loop over available jet radii    
    
    
    //if(debugMode)std::cout<<"closing gPad..."<<std::endl<<std::endl;    
    //gPad->Close();    
    
    
    
    // DRAW THOSE PDFS //    
    std::string thePDFFileName=outputDir+fullJetType+"_MCJEC_"+outputTag+".pdf";    
    std::string open_thePDFFileName=thePDFFileName+"[";    
    std::string close_thePDFFileName=thePDFFileName+"]";    
    
    // for mu v. gen pt and sigma/mu v. genpt hists_fit
    TCanvas* pdfoutCanv=new TCanvas("outputPdf","outputPdf", 800, 800);    
    pdfoutCanv->cd();        
    std::cout<<"opening output pdf file "<< thePDFFileName<<std::endl<<std::endl;    
    pdfoutCanv->Print( open_thePDFFileName.c_str() );    
    

    TCanvas* pdfoutCanv_muSigma=new TCanvas("outputPdfwLog_muSigma","outputPdfwLog", 900, 700);    
    pdfoutCanv_muSigma->Divide(1,3);
    pdfoutCanv_muSigma->cd();
    
    TPad* p1=(TPad*)pdfoutCanv_muSigma->cd(1);        
    p1->SetLogx(1);    
    p1->SetGridx(1);    
    p1->cd();
    
    TPad* p2=(TPad*)pdfoutCanv_muSigma->cd(2);        
    p2->SetLogx(1);    	
    p2->SetGridx(1);    
    p2->cd();    

    TPad* p3=(TPad*)pdfoutCanv_muSigma->cd(3);        
    p3->SetLogx(1);    	
    p3->SetGridx(1);    
    p3->cd();    
    
    
    // draw and print pdfs. //    
    std::cout<<std::endl<<"printing and fitting mu/sigma hists ..."<<std::endl;    
    for(int i=0;i<Nrad;++i){    
      
      p1->cd();
      
      hMean[i]->Draw("HIST E1"); 
      
      TLine* meanLine=new TLine(ptbins_debug[0],1.,ptbins_debug[nbins_pt_debug],1.);    
      meanLine->SetLineStyle(2);  
      meanLine->SetLineColor(kBlue);    
      meanLine->Draw();    
      
      p2->cd();
      
      
      
//      TF1 *hSigmaFit=new TF1("hSigmaFit","[0]+[1]*pow(x,[2])", 
//			     ptbins_debug[2],                     // I want this fit to start at 56 GeV or so
//			     ptbins_debug[nbins_pt_debug]);       // I want this fit to end at 967 GeV (for entire range)

      TF1 *hSigmaFit=new TF1("hSigmaFit","[0]+[1]/(pow(x,[2])+[3]*x)", 
			     ptbins_debug[2],                     // I want this fit to start at 56 GeV or so
			     ptbins_debug[nbins_pt_debug-1]);       // I want this fit to end at 967 GeV (for entire range)

      
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
      //sigmaFitStatus=hSigma[i]->Fit(hSigmaFit,"MEVR");
      //sigmaFitStatus=hSigma[i]->Fit(hSigmaFit,"MR");
      sigmaFitStatus=hSigma[i]->Fit(hSigmaFit,"R");	

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
      
      hSigma[i]->Draw("HIST E1");
      hSigmaFit->Draw("SAME");
      

      
      p3->cd();

      hChi2NDF[i]->Draw("HIST E1");

      pdfoutCanv_muSigma->Print(thePDFFileName.c_str());        

      hMean[i]->Write();               
      hSigma[i]->Write();               
      hSigmaFit->Write();
      hChi2NDF[i]->Write();                     
      pdfoutCanv_muSigma->Write();



      




      TCanvas* pdfoutCanv_wLogy=new TCanvas("outputPdfwLogx","output Pdf wLogx", 800, 800);    
      pdfoutCanv_wLogy->SetLogy(1);
      pdfoutCanv_wLogy->cd();            
      




  
      
      std::cout<<"printing JER reco/gen hists"<<std::endl;

      
      for(int j=0;j<nbins_pt_debug;++j){    
	
	int ptbin_j=(int)ptbins_debug[j];
	int ptbin_j1=(int)ptbins_debug[j+1];
	std::string hrspTitle=std::to_string(ptbin_j)+" GeV < gen jet p_{T} < "+std::to_string(ptbin_j1)+" GeV";    
	std::string hrsp_XAxTitle="reco p_{T}/gen p_{T}";

	makeJERSHists(  ( TCanvas*    ) pdfoutCanv_wLogy ,  
		        ( std::string ) thePDFFileName	 ,  
		        ( TH1F*       ) hrsp[i][j]      , 
 		        ( std::string ) hrspTitle        ,
	 	        ( std::string ) hrsp_XAxTitle      );
	
      }    
      

      /*
      std::cout<<"printing JER raw/gen hists"<<std::endl;
      for(int j=0;j<nbins_pt_debug;++j){    
	
	int ptbin_j=(int)ptbins_debug[j];
	int ptbin_j1=(int)ptbins_debug[j+1];
	std::string hrspTitle=std::to_string(ptbin_j)+" GeV < gen jet p_{T} < "+std::to_string(ptbin_j1)+" GeV";    
	std::string hrsp_XAxTitle="raw p_{T}/gen p_{T}";	
	makeJERSHists(  ( TCanvas*    ) pdfoutCanv_wLogy ,  
		        ( std::string ) thePDFFileName	 ,  
		        ( TH1F*       ) hrawrsp[i][j]    , 
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
		        ( TH1F*       ) hcorr[i][j]      , 
 		        ( std::string ) hrspTitle        ,
	 	        ( std::string ) hrsp_XAxTitle      );      }    
      */      
      
      std::cout<<"closing output pdf file "<<thePDFFileName<<std::endl<<std::endl;                
      
      pdfoutCanv->Print( close_thePDFFileName.c_str() );    
      
      pdfoutCanv->Close();
      pdfoutCanv_muSigma->Close();
      pdfoutCanv_wLogy->Close();
      

    }              
    
  }
  else std::cout<<std::endl<<std::endl<<"skipping hJER hists"<<std::endl<<std::endl;
                       
  


  /////////////////////////////////////////////  /////////////////////////////////////////////  /////////////////////////////////////////////  /////////////////////////////////////////////
  /////////////////////////////////////////////  /////////////////////////////////////////////  /////////////////////////////////////////////  /////////////////////////////////////////////
  /////////////////////////////////////////////  /////////////////////////////////////////////  /////////////////////////////////////////////  /////////////////////////////////////////////

  
  
  if(!draw_MCEff)  std::cout<<std::endl<<std::endl<<"skipping MCEff hists"<<std::endl<<std::endl;
  else{
    
    { 
      
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
  }
  
  if(debugMode)std::cout<<"closing output root file"<<std::endl<<std::endl;
  rootfout->Close();
  return 0;
}





int main(int argc, char*argv[]){
  int rStatus=-1;
  if(argc!=3){
    std::cout<<"no defaults. Do...."<<std::endl;
    std::cout<<"./printPlots_ppMC_JERS.exe "<<
      "<target_ppMC_dir> <output tag>"<<std::endl<<std::endl;
    return rStatus;  }
  
  rStatus=1;
  rStatus=printPlots_ppMC_JERS((const std::string)argv[1], (const std::string)argv[2] );
  return rStatus;
}
