#include "printPlots.h"


const bool debugMode=true;
const bool drawProfiles = false;
const bool rebinJER=true;
const int rebinFactor=1;//5

//const int fit_ptlo_bin=3, fit_pthi_bin=nbins_pt_debug-1;
const int fit_ptlo_bin=1, fit_pthi_bin=nbins_pt_debug-5;
//const int N_JERfits=nbins_pt_debug-fit_ptlo_bin;
const int N_JERfits=fit_pthi_bin-fit_ptlo_bin;
const double fgaus_xlo=0.80, fgaus_xhi=1.11;
const double defbinsize=0.01;

//other options
const bool draw_hJERRapBins=false, doGenBinsToo=false;//RapBins -> dual-diff xsec bins, GenBins -> variable, depends on readForests
const bool draw_JERgen150to200=false, draw_JERgen30to50=false;
const std::string didJetID="1";
const std::string doJetID="1";
const int CANVX=1200, CANVY=1200;
const bool doabsetabinJER=true;

int printPlots_ppMC_JERS(std::string inFile_MC_dir,const std::string outputTag, const bool draw_hJER=true, const bool draw_MCEff=false, const int absetabin=-1){
  
  //if(absetabin!=-1) doabsetabinJER=true;
  
  // root style settings.
  std::cout<<"forcing style"<<std::endl;
  LoadStyle();
  
  float difffromone=1.-fgaus_xlo;
  int rebin_hrsp_arr_N=( (int) (100.*( ( (fgaus_xhi+difffromone) - (fgaus_xlo-difffromone) ) ) )+ 1 );
  if(debugMode)std::cout<<"rebin_hrsp_arr_N="<<rebin_hrsp_arr_N<<std::endl;
  double rebin_hrsp_arr[rebin_hrsp_arr_N]={0.};
  for(int i=0; i<rebin_hrsp_arr_N;i++) rebin_hrsp_arr[i]=(fgaus_xlo-difffromone)+i*defbinsize;
  if(debugMode) for(int i=0; i<rebin_hrsp_arr_N;i++) std::cout<<"rebin_hrsp_arr["<<i<<"]="<<rebin_hrsp_arr[i]<<std::endl;
  
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
  TDirectory* fout_ptbin_JER=rootfout->mkdir("all_ptbins_JER");
  
  // fit + fit settings
  std::cout<<"setting up fitters..."<< std::endl;
  TVirtualFitter::SetDefaultFitter("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04); 
  
  TH1F *hrsp[nbins_pt_debug]={};  //input
  
  TH1F *hMean=NULL, *hSigma=NULL; //output     
  TH1F *hMean_fit=NULL, *hSigma_fit=NULL; 
  TH1F *hChi2NDF=NULL, *hChi2Prob=NULL;

  Double_t fit_ptlo=ptbins_debug[fit_ptlo_bin];
  Double_t fit_pthi=ptbins_debug[fit_pthi_bin];
  TF1 *fgaus=NULL; 
  
  if(!draw_hJER)
    std::cout<<std::endl<<std::endl<<"skipping hJER hists"<<std::endl<<std::endl;
  else {
    
    const std::string algname=algo+radiusInt+jetType; 
    std::cout<<"making hists for algorithm "<< algname<<"..."<<std::endl;
    
    // initialize mean hist    
    hMean = new TH1F( "hMean", Form( "Mean %s", algname.c_str()),    
		      nbins_pt_debug, ptbins_debug);
    
    hSigma = new TH1F( "hSigma", Form( "Sigma %s", algname.c_str()),    
		       nbins_pt_debug, ptbins_debug);
    
    hMean_fit = new TH1F( "hMean_fit", Form( "Fit Mean %s", algname.c_str()),    
			  nbins_pt_debug, ptbins_debug);
    
    hSigma_fit = new TH1F( "hSigma_fit", Form( "Fit Sigma %s", algname.c_str()),    
			   nbins_pt_debug, ptbins_debug);
    
    hChi2NDF = new TH1F("hChi2NDF", Form( "Fit Chi2NDF %s", algname.c_str()),    
			nbins_pt_debug, ptbins_debug); 
    
    hChi2Prob = new TH1F("hChi2Prob", Form( "Fit Chi2Prob %s", algname.c_str()),    
			 nbins_pt_debug, ptbins_debug); 
    
    // pt bin hJER-fit loop    
    if(debugMode)std::cout<<"nbins_pt="<<nbins_pt_debug<<std::endl;    
    for(int ip=0; ip<nbins_pt_debug; ip++){    
      int ptbin_ip=(int)ptbins_debug[ip], ptbin_ip1=(int)ptbins_debug[ip+1];
      
      // open the input hist    
      std::string inputHistName="";
      if(doabsetabinJER)
	inputHistName="hJER_"+doJetID+"wJetID_etabin"+std::to_string(absetabin)+"_"+std::to_string(ptbin_ip)+"ptbin"+std::to_string(ptbin_ip1);      
      else 
	inputHistName="hJER_"+doJetID+"wJetID_ptbin"+std::to_string(ip);

      hrsp[ip] = (TH1F*)finPP->Get( inputHistName.c_str() );        
     
      if(!hrsp[ip]){ 
	std::cout<<"no input hist named " <<  inputHistName<< ", continue..."<<std::endl;
	continue;
	//assert(false);
      }          
      else if(hrsp[ip]->GetEntries()<3)
	continue;	
      

      if(debugMode)std::cout<<std::endl;  
      if(debugMode)std::cout<<"________________________________"<<std::endl;  
      if(debugMode)std::cout<<"pt range for bin: "<<ptbin_ip<<" - "<<ptbin_ip1<< " GeV "<<std::endl;    
      if(debugMode)std::cout<<"hist title = "<<hrsp[ip]->GetTitle()<<std::endl;
      if(doabsetabinJER)
	hrsp[ip]->SetName(("hJER_1wJetID_receta_"+
			   absetabins_str[absetabin]+"_to_"+absetabins_str[absetabin+1]+
			   "_genpt_"+std::to_string(ptbin_ip) + "_to_" + std::to_string(ptbin_ip1) ).c_str());
      else
	hrsp[ip]->SetName(("hJER_1wJetID_genpt_"+std::to_string(ptbin_ip) + "_to_" + std::to_string(ptbin_ip1) ).c_str());
      if(debugMode)hrsp[ip]->Print("base");    
      if(debugMode)std::cout<<std::endl;  
      
      if(rebinJER&&rebinFactor!=1)
	hrsp[ip]= (TH1F*)hrsp[ip]->TH1::Rebin( rebinFactor, (inputHistName+"_rebin5").c_str());	
      
      hrsp[ip]->Scale( 1./ hrsp[ip]->Integral() );    //scale to area of 1
      
      // get some quick stats    
      TH1F* hrsp_ip_rebin=(TH1F*)hrsp[ip]->Clone();
      hrsp_ip_rebin=(TH1F*)(hrsp_ip_rebin->TH1::Rebin( (rebin_hrsp_arr_N-1), 
						       ((std::string)hrsp[ip]->GetName()+"_rebin").c_str() , 
						       rebin_hrsp_arr));
      hrsp_ip_rebin->TH1::ResetStats();
      if(debugMode)hrsp_ip_rebin->Print("base");
      
      //float h_norm  = hrsp_ip_rebin->GetMaximumStored();          
      float h_mean  = hrsp_ip_rebin->GetMean();    
      float h_emean = hrsp_ip_rebin->GetMeanError();    
      float h_sig   = hrsp_ip_rebin->GetStdDev()/h_mean;    
      float h_esig  = h_sig * sqrt (  
				    pow( h_emean/h_mean,2) + pow( hrsp_ip_rebin->GetStdDevError()/hrsp_ip_rebin->GetStdDev(),2) ) ; 
      
      // set contents+errors in specific bin    
      std::cout<<"setting bin content for hMean[ip]=["<< ip<<"] "<<std::endl;    
      std::cout<<"h_mean ="   << h_mean<<std::endl;    
      std::cout<<"h_emean ="  << h_emean<<std::endl;    
      std::cout<<"h_sig ="    << h_sig<<std::endl;    
      std::cout<<"h_esig ="   << h_esig<<std::endl;    
      ////h_norm  = hrsp[ip]->GetMaximumStored();          
      //h_mean  = hrsp[ip]->GetMean();    
      //h_emean = hrsp[ip]->GetMeanError();    
      //h_sig   = hrsp[ip]->GetStdDev()/h_mean;    
      //h_esig  = h_sig * sqrt (  
      //			      pow( h_emean/h_mean,2) + pow( hrsp[ip]->GetStdDevError()/hrsp[ip]->GetStdDev(),2) ) ; 
      //std::cout<<"setting bin content for hMean[ip]=["<< ip<<"] "<<std::endl;    
      //std::cout<<"h_mean ="   << h_mean<<std::endl;    
      //std::cout<<"h_emean ="  << h_emean<<std::endl;    
      //std::cout<<"h_sig ="    << h_sig<<std::endl;    
      //std::cout<<"h_esig ="   << h_esig<<std::endl;    
      
      hrsp[ip]->Write();
      if(ip<fit_ptlo_bin )continue;
      if(ip>=fit_pthi_bin )continue;
      
      hMean->SetBinContent (ip+1, h_mean);    
      hMean->SetBinError   (ip+1, h_emean);    
      
      hSigma->SetBinContent (ip+1, h_sig);    
      hSigma->SetBinError   (ip+1, h_esig);    
      
      fout_ptbin_JER->cd();
      
      
      std::cout<<std::endl;
      std::cout<<"fitting gauss core for "<<hrsp[ip]->GetName()<<std::endl;   
      fgaus = new TF1("fgaus","gaus", fgaus_xlo,fgaus_xhi);  //what's been typicall used.      
      //fgaus->SetParLimits(0,0.98,1.02);    //normalization
      //fgaus->SetParLimits(1,0.95,1.05);   // mean
      //fgaus->SetParLimits(2,0.0,0.5);     // width
      
      bool fitfailed=true;
      fitfailed = ((bool)(hrsp[ip]->Fit(fgaus,"R")));    
      
      
      
      if(fitfailed)
	std::cout<<std::endl<<"!!!WARNING FIT FAILED!!!"<<std::endl<<std::endl;
      else 
	std::cout<<std::endl<<"fit success"<<std::endl<<std::endl;
      
      float f_mean  = -1.;
      float f_emean = -1.;
      float f_sig   = -1.;
      float f_esig  = -1.;      
      float chi2NDF = -1.;
      float chi2prob= -1.;      
      
      f_mean  = (fitfailed)  ? -1. : fgaus->GetParameter(1); // z = if(condition) then(?) <do this> else(:) <do this>   
      f_emean = (fitfailed)  ? -1. : fgaus->GetParError(1) ; 
      f_sig   = (fitfailed)  ? -1. : fgaus->GetParameter(2)/f_mean;    
      f_esig  = (fitfailed)  ? -1. : f_sig * sqrt( pow( (f_emean/f_mean), 2) + 
						   pow( ( fgaus->GetParError(2)/fgaus->GetParameter(2) ), 2) 
						   )   ;
      chi2NDF = (fitfailed)  ? -1. : fgaus->GetChisquare()/((float)fgaus->GetNDF());
      chi2prob = (fitfailed) ? -1. : fgaus->GetProb();
      
      std::cout<<"setting bin content for fit hists ip="<< ip<<" "<<std::endl;    
      std::cout<<"f_mean ="   << f_mean<<std::endl;    
      std::cout<<"f_emean ="  << f_emean<<std::endl;    
      std::cout<<"f_sig ="    << f_sig<<std::endl;    
      std::cout<<"f_esig ="   << f_esig<<std::endl;    
      std::cout<<"Chi2NDF ="   << chi2NDF<<std::endl;    
      std::cout<<"Chi2Prob ="   << chi2prob<<std::endl;          
      
      hMean_fit->SetBinContent (ip+1, f_mean);    
      hMean_fit->SetBinError   (ip+1, f_emean);    
      
      hSigma_fit->SetBinContent (ip+1, f_sig);    
      hSigma_fit->SetBinError   (ip+1, f_esig);    
      
      hChi2NDF->SetBinContent (ip+1, chi2NDF );
      hChi2NDF->SetBinError (ip+1, 1. );
      
      hChi2Prob->SetBinContent (ip+1, chi2prob );      
    }// end fit-loop over ptbins    
    
    rootfout->cd();//don't wanna be in directory anymore
    
    
    // for mu v. gen pt and sigma/mu v. genpt hists_fit    
    std::cout<<std::endl<<"fitting JER v genpt hists ..."<<std::endl;        
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-02); 
     
    TF1 *SigmaFit_h=new TF1("SigmaFit_h","[0]+[1]/(pow(x,[2])+[3]*x)", 
			    fit_ptlo,
			    fit_pthi );       // I want this fit to end at 967 GeV (for entire range)    
    
    bool hsigmaFitStatus= true;
    std::cout<<std::endl<<"fitting JER width v. pt... init. sigmaFitStatus="<<hsigmaFitStatus<<std::endl<<std::endl; 
    //sigmaFitStatus=hSigma->Fit(SigmaFit_h,"MEVR");    //sigmaFitStatus=hSigma->Fit(SigmaFit_h,"MR");
    hsigmaFitStatus=hSigma->Fit(SigmaFit_h,"MEVR");	
    hsigmaFitStatus=hSigma->Fit(SigmaFit_h,"EVR");	


    TF1 *SigmaFit_f=new TF1("SigmaFit_f","[0]+[1]/(pow(x,[2])+[3]*x)", 
			   fit_ptlo,
			   fit_pthi );       // I want this fit to end at 967 GeV (for entire range)    

    bool fsigmaFitStatus= true;
    std::cout<<std::endl<<"fitting JER width v. pt... init. sigmaFitStatus="<<fsigmaFitStatus<<std::endl<<std::endl; 
    //sigmaFitStatus=hSigma->Fit(SigmaFit_h,"MEVR");    //sigmaFitStatus=hSigma->Fit(SigmaFit_h,"MR");
    fsigmaFitStatus=hSigma->Fit(SigmaFit_f,"MEVR");	
    fsigmaFitStatus=hSigma_fit->Fit(SigmaFit_f,"EVR");	


    
    float fitParam_0    =SigmaFit_h->GetParameter(0);
    float fitParam_0_err=SigmaFit_h->GetParError(0);
    float fitParam_1    =SigmaFit_h->GetParameter(1);
    float fitParam_1_err=SigmaFit_h->GetParError(1);
    float fitParam_2    =SigmaFit_h->GetParameter(2);
    float fitParam_2_err=SigmaFit_h->GetParError(2); 
    float fitParam_3    =SigmaFit_h->GetParameter(3);
    float fitParam_3_err=SigmaFit_h->GetParError(3); 

    std::cout << "fitParam_0     = " << fitParam_0    << std::endl;
    std::cout << "fitParam_0_err = " << fitParam_0_err<< std::endl;
    std::cout << "%err fitParam0 = " << (fitParam_0_err/fitParam_0)*100.    << std::endl;
    std::cout << "fitParam_1     = " << fitParam_1    << std::endl;
    std::cout << "fitParam_1_err = " << fitParam_1_err<< std::endl;
    std::cout << "%err fitParam1 = " << (fitParam_1_err/fitParam_1)*100.    << std::endl;
    std::cout << "fitParam_2     = " << fitParam_2    << std::endl;
    std::cout << "fitParam_2_err =   " << fitParam_2_err<< std::endl;
    std::cout << "%err fitParam2 = " << (fitParam_2_err/fitParam_2)*100.    << std::endl;
    std::cout << "fitParam_3     = " << fitParam_3    << std::endl;
    std::cout << "fitParam_3_err =   " << fitParam_3_err<< std::endl;
    std::cout << "%err fitParam3 = " << (fitParam_3_err/fitParam_3)*100.    << std::endl;        
    
    
    TCanvas* pdfoutCanv_muSigma=new TCanvas("outputPdfwLog_muSigma","outputPdfwLog", CANVX, CANVY);    
    pdfoutCanv_muSigma->Divide(1,3);    //    pdfoutCanv_muSigma->Divide(1,4);
    pdfoutCanv_muSigma->cd();
    
    float xmin=ptbins_debug[0], xmax=ptbins_debug[fit_pthi_bin];
    if(debugMode){
      std::cout<<"xmin="<<xmin<<std::endl;
      std::cout<<"xmax="<<xmax<<std::endl;    }

    TLine* lineAtOne=new TLine(xmin,1.,ptbins_debug[fit_pthi_bin+1],1.);    
    lineAtOne->SetLineStyle(2);  
    lineAtOne->SetLineColor(kBlack);
    TLine* lineAtTen=new TLine(xmin,10.,ptbins_debug[fit_pthi_bin+1],10.);    
    lineAtTen->SetLineStyle(2);  
    lineAtTen->SetLineColor(kBlack);
    

    
    TPad* p1=(TPad*)pdfoutCanv_muSigma->cd(1);        
    p1->SetLogx(1);    
    p1->SetGridx(0);    
    p1->SetTopMargin(0.05);
    p1->cd();
    
    TPad* p2=(TPad*)pdfoutCanv_muSigma->cd(2);        
    p2->SetLogx(1);    	
    p2->SetGridx(0);    
    p2->SetTopMargin(0.);
    p2->cd();   
    
    TPad* p3=(TPad*)pdfoutCanv_muSigma->cd(3);        
    p3->SetLogx(1);    	
    p3->SetLogy(1);    	
    p3->SetGridx(0);
    p3->SetTopMargin(0.);
    p3->cd();        

//    TPad* p4=(TPad*)pdfoutCanv_muSigma->cd(4);        
//    p4->SetLogx(1);    	
//    p4->SetLogy(1);    	
//    p4->SetGridx(0);
//    p4->cd();        
    
    
    //draw pad 1
    p1->cd();    
    MakeHistMean( (TH1F*)hMean,xmin,xmax);
    hMean->DrawClone("HIST ][E"); 
    
    MakeHistMean( (TH1F*)hMean_fit, xmin,xmax);   
    hMean_fit->SetMarkerStyle(kOpenStar);
    hMean_fit->SetMarkerColor(kBlack);
    hMean_fit->SetLineColor(kMagenta);
    hMean_fit->DrawClone("HIST ][E SAME");

    TLegend* hMean_Leg=new TLegend( 0.73, 0.75 ,0.9 , 0.9, NULL,"brNDC");
    hMean_Leg->SetBorderSize(0);
    hMean_Leg->SetFillStyle(0);

    hMean_Leg->AddEntry( hMean, "#mu from Hists", "lp");
    hMean_Leg->AddEntry( hMean_fit, "#mu from Gauss Fits", "lp");
    hMean_Leg->Draw();

    lineAtOne->Draw();   
    
    //draw pad 2
    p2->cd();    
    
    MakeHistRMS( (TH1F*) hSigma ,xmin,xmax); 
    hSigma->DrawClone("HIST ][E");
    
    //this is a hist created from the gauss fits in each bin
    hSigma_fit->SetMarkerStyle(kOpenStar);
    hSigma_fit->SetMarkerColor(kBlack);
    hSigma_fit->SetLineColor(kMagenta);
    hSigma_fit->DrawClone("HIST ][E SAME");
    
    //this is a TF1; this is the fit to hSigma_fit
    //SigmaFit_h->SetLineWidth(1);
    SigmaFit_h->SetLineColor(kBlue-7);
    SigmaFit_h->Draw("SAME");
    
    //SigmaFit_f->SetLineWidth(1);
    SigmaFit_f->SetLineColor(kMagenta-2);
    SigmaFit_f->Draw("SAME");

    TLegend* hSigma_Leg=new TLegend( 0.65, 0.50 , 0.95 , 0.97, NULL,"brNDC");

    std::string SigmaFit_h_prob_str;
    makeFitProbStr(SigmaFit_h, &SigmaFit_h_prob_str);
    if(debugMode)std::cout<<"SigmaFit_h_prob_str="<<SigmaFit_h_prob_str<<std::endl;    
    std::string SigmaFit_h_chi2_str;
    makeFitChi2Str(SigmaFit_h, &SigmaFit_h_chi2_str);
    if(debugMode)std::cout<<"SigmaFit_h_chi2_str="<<SigmaFit_h_chi2_str<<std::endl;    
    std::string SigmaFit_h_ndf_str=std::to_string( SigmaFit_h->GetNDF() );
    if(debugMode)std::cout<<"SigmaFit_h_ndf_str="<<SigmaFit_h_ndf_str<<std::endl;    
    
    hSigma_Leg->AddEntry( hSigma, "#sigma/#mu from Hists", "lp");
    hSigma_Leg->AddEntry( SigmaFit_h, "Fit of #sigma/#mu from Hists", "l");   
    hSigma_Leg->AddEntry( (TObject*)0, ("#chi^{2}/NDF="+SigmaFit_h_chi2_str+"/"+SigmaFit_h_ndf_str).c_str() , "");   
    hSigma_Leg->AddEntry( (TObject*)0, ("Prob="+SigmaFit_h_prob_str).c_str() , "");   
    
    
    std::string SigmaFit_f_prob_str;
    makeFitProbStr(SigmaFit_f, &SigmaFit_f_prob_str);
    if(debugMode)std::cout<<"SigmaFit_f_prob_str="<<SigmaFit_f_prob_str<<std::endl;    
    std::string SigmaFit_f_chi2_str;
    makeFitChi2Str(SigmaFit_f, &SigmaFit_f_chi2_str);
    if(debugMode)std::cout<<"SigmaFit_f_chi2_str="<<SigmaFit_f_chi2_str<<std::endl;    
    std::string SigmaFit_f_ndf_str=std::to_string( SigmaFit_f->GetNDF() );
    if(debugMode)std::cout<<"SigmaFit_f_ndf_str="<<SigmaFit_f_ndf_str<<std::endl;    

    hSigma_Leg->AddEntry( (TObject*)0, "" , "");   
    hSigma_Leg->AddEntry( hSigma_fit, "#sigma/#mu from Gauss Fits", "lp");
    hSigma_Leg->AddEntry( SigmaFit_f, "Fit of #sigma/#mu from Gauss Fits", "l");
    hSigma_Leg->AddEntry( (TObject*)0, ("#chi^{2}/NDF="+SigmaFit_f_chi2_str+"/"+SigmaFit_f_ndf_str).c_str() , "");   
    hSigma_Leg->AddEntry( (TObject*)0, ("Prob="+SigmaFit_f_prob_str).c_str() , "");   


    hSigma_Leg->SetBorderSize(0);    
    hSigma_Leg->SetFillStyle(0);        
    hSigma_Leg->Draw();
    
    //draw pad 3
    p3->cd();    
    MakeHistChi2NDF( (TH1F*)hChi2NDF,xmin,xmax);           
    hChi2NDF->DrawClone("HIST ][E");    
    lineAtOne->Draw();                
    lineAtTen->Draw();            
    //hChi2NDF->DrawClone("HIST ][E SAME");    
    
    TLegend* hChi2NDF_Leg=new TLegend(0.65,0.85 ,0.89 , 1.0, NULL,"brNDC");
    hChi2NDF_Leg->AddEntry(hChi2NDF,"#chi^{2}/NDF from Gauss Core Fits");
    hChi2NDF_Leg->SetBorderSize(0);
    hChi2NDF_Leg->SetFillStyle(0);
    hChi2NDF_Leg->Draw();
    
    ////draw pad 4
    //p4->cd();    
    //MakeHistChi2Prob( (TH1F*)hChi2Prob,xmin,xmax);           
    //hChi2Prob->DrawClone("HIST ][E");    
    
    
    // DRAW THOSE PDFS //        
    std::string thePDFFileName=outputDir+fullJetType+"_MCJEC_"+outputTag+".pdf";    
    std::string open_thePDFFileName=thePDFFileName+"[";    
    std::string close_thePDFFileName=thePDFFileName+"]";    
    pdfoutCanv_muSigma->Print(open_thePDFFileName.c_str());            
    pdfoutCanv_muSigma->Print(thePDFFileName.c_str());            

    



    






    
    std::cout<<"printing JER reco/gen hists"<<std::endl;

    TCanvas* pdfoutCanv_wLogy=new TCanvas("JERfitCanv","Canv for JER fit", CANVX, CANVY);    
    pdfoutCanv_wLogy->SetLogy(1);
    pdfoutCanv_wLogy->cd();            
    //pdfoutCanv_wLogy->Print( open_thePDFFileName.c_str() );        
    
    
    
    //for(int j=fit_ptlo_bin;j<fit_pthi_bin;++j){       //if(j<3)continue;	
    for(int j=0;j<nbins_pt_debug;++j){       //if(j<3)continue;	
      
      int ptbin_j=(int)ptbins_debug[j];
      int ptbin_j1=(int)ptbins_debug[j+1];
      std::string hrspTitle=std::to_string(ptbin_j)+" GeV < GEN Jet p_{T} < "+std::to_string(ptbin_j1)+" GeV";    
      std::string hrsp_XAxTitle="RECO p_{T}/GEN p_{T}";  
      //      		     ( std::string ) thePDFFileName   ,      
      makeJERSHists( ( TCanvas*    ) pdfoutCanv_wLogy ,  
		     ( TH1F*       ) hrsp[j]          , 
		     ( std::string ) hrspTitle        ,
		     ( std::string ) hrsp_XAxTitle    ,
		     ( TFile*      ) rootfout         ,
		     ( TDirectory*  ) fout_ptbin_JER );
    }

    
    std::cout<<"N_JERfits="<<N_JERfits<<std::endl;
    int canv_count=(N_JERfits/9);
    std::cout<<"canv_count="<<canv_count<<std::endl;
    if(N_JERfits%9>0)canv_count++;
    std::cout<<"canv_count="<<canv_count<<std::endl;
    const int canv_count_const=canv_count;
    
    TCanvas* pdfoutCanv_3x3_arr[canv_count_const]={};
    for(int i=0;i<canv_count_const;i++){
      std::cout<<"i="<<i<<std::endl;
      pdfoutCanv_3x3_arr[i]=new TCanvas(("JERfitCanv_3x3_"+std::to_string(i)).c_str(), ("Canv #"+std::to_string(i)+" for JERfits").c_str(), CANVX, CANVY);
      pdfoutCanv_3x3_arr[i]->Divide(3,3);
      int fit_ptbin_start=i*9+fit_ptlo_bin;
      int fit_ptbin_end=fit_ptbin_start+9;
      if(fit_ptbin_end>fit_pthi_bin)fit_ptbin_end=fit_pthi_bin;
      //if(fit_ptbin_end>nbins_pt_debug)fit_ptbin_end=nbins_pt_debug;
      int drawcount=0;
      
      for(int j=fit_ptbin_start;j<fit_ptbin_end;++j){    
	std::cout<<"j="<<j<<std::endl;
	std::cout<<"drawing pad for "<<ptbins_debug[j]<<" < jet p_T < "<<ptbins_debug[j+1]<<std::endl;
	if(debugMode)hrsp[j]->Print("base");
	pdfoutCanv_3x3_arr[i]->cd(drawcount+1);
	pdfoutCanv_3x3_arr[i]->cd(drawcount+1)->SetLogy(1);
	hrsp[j]->Draw("E1")    ;	
	drawcount++;	
      }//end loop over pads of one of the 3x3 canvases
      
      pdfoutCanv_3x3_arr[i]->Print(thePDFFileName.c_str());
      rootfout->cd();
      std::string fit_pt_start_str= std::to_string(((int) ptbins_debug[fit_ptbin_start]));
      std::string fit_pt_end_str= std::to_string(((int) ptbins_debug[fit_ptbin_end]));
      if(debugMode){
	std::cout<<"fit_ptbin_start ="<< fit_ptbin_start     <<std::endl;
	std::cout<<"fit_pt_start_str ="<< fit_pt_start_str <<std::endl;
	std::cout<<"fit_ptbin_end   ="<< fit_ptbin_end       <<std::endl;
	std::cout<<"fit_pt_end_str   ="<< fit_pt_end_str   <<std::endl;
      }
      pdfoutCanv_3x3_arr[i]->SetTitle(("3x3 JER fit canv for bin"+ fit_pt_start_str + "to bin"+fit_pt_end_str ).c_str());
      pdfoutCanv_3x3_arr[i]->Write(("JERfitCanv_3x3_bins_"       + fit_pt_start_str + "_to_"  +fit_pt_end_str ).c_str());
    }//end loop over 3x3canvases

    


    std::cout<<"closing output pdf file "<<thePDFFileName<<std::endl<<std::endl;                
    pdfoutCanv_muSigma->Print( close_thePDFFileName.c_str() );        
    
    rootfout->cd();
    
    hMean->Write();               
    hSigma->Write();               
    hChi2NDF->Write();                     
    SigmaFit_h->Write();
    SigmaFit_f->Write();
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
      //the2DMCEffQAHist->SetAxisRange(1.*pow(10.,-13.), 1.*pow(10.,-4.),"Z");	
      
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
	the2DMCEffQAHist->Rebin2D(2,2, (inHistName+"_rebin2").c_str());
	TH1F* theMCEffQAHist= (TH1F*)the2DMCEffQAHist->TH2::ProfileX( inHistName.c_str(), 1, -1, "o");
	if(  genVars[j]=="pt" ||
	     genVars[j]=="pt2" ||
	     genVars[j]=="pt3" )  {
	  theMCEffQAHist->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");	
	  theMCEffQAHist->SetAxisRange( 0.90, 1.10,"Y");	 
	}
	else if (genVars[j]=="ptrat_phi" || 
		 genVars[j]=="ptrat_eta" ||
		 genVars[j]=="ptrat_drjt" 
		 ) {
	  theMCEffQAHist->SetAxisRange( 0.90, 1.15,"Y");	 
	  std::cout<<std::endl<<"rebinning TH2 Profile"<<std::endl;
	  
	  inHistName+="_rebin2";
	  //theMCEffQAHist=(TH1F*)theMCEffQAHist->Rebin(2,  inHistName.c_str());
	  //theMCEffQAHist->Rebin(2,  inHistName.c_str());
	  //theMCEffQAHist->Scale(0.5);
	  if(genVars[j]=="ptrat_drjt")
	    theMCEffQAHist->SetAxisRange(0.0,0.11,"X");
	}

	
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
  if(argc<5){
    std::cout<<"no defaults. Do...."<<std::endl;
    std::cout<<"./printPlots_ppMC_JERS.exe <target_ppMC_dir> <output tag> <draw hJER> <draw MCEFf> <absetabin>"<<std::endl<<std::endl;
    return rStatus;  }
  
  rStatus=1;
  if(argc==6)
    rStatus=printPlots_ppMC_JERS((const std::string)argv[1], (const std::string)argv[2] , (const bool)std::atoi(argv[3]), (const bool)std::atoi(argv[4]), (const int) std::atoi(argv[5]) );
  if(argc==5)
    rStatus=printPlots_ppMC_JERS((const std::string)argv[1], (const std::string)argv[2] , (const bool)std::atoi(argv[3]), (const bool)std::atoi(argv[4]), -1);
  return rStatus;
}
