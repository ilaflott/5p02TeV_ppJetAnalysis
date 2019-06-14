#include "printPlots.h"


//entire macro settings
const bool debugMode=true;
const int CANVX=1200, CANVY=1200;
const std::string doJetID="1";

//// JER settings
//hrsp
const bool doGenYBins=true;//RapBins -> dual-diff xsec bins, GenBins -> variable, depends on readForests
const bool rebinJER=true;
const int rebinFactor=1;// test
const double hgaus_xlo=0.70, hgaus_xhi=1.+(1.-hgaus_xlo);//test
const bool useconstxrange=false;
const double fgaus_xlo=0.80, fgaus_xhi=1.+(1.-fgaus_xlo);//test
const double defbinsize=0.01;//dont touch
//resolution v genpt fit options
const double errfact=1.0;// test



//MCEff settings
const bool drawProfiles = false;//rarely used




int printPlots_ppMC_JERS_etabin(std::string inFile_MC_dir,const std::string outputTag, const bool draw_hJER=true, const bool draw_MCEff=false, const int absetabin=0, const std::string fitquantity="sigma"){
  
  //if(absetabin!=-1) doabsetabinJER=true;
  
  // root style settings.
  std::cout<<"forcing style"<<std::endl;
  LoadStyle();

  bool fitsigma=false, fitsigmu=false;
  if(fitquantity=="sigma"){
    fitsigma=true; fitsigmu=false;}
  else if(fitquantity=="sigmu"){
    fitsigma=false; fitsigmu=true;}
  else assert(false);
  //const bool fitsigma=false, fitsigmu=!fitsigma;//fit either sigma, or sigma/mu


  int fit_ptlo_bin=-1, fit_pthi_bin=nbins_pt_debug+1;//test
  if(absetabin==0)      fit_ptlo_bin=fit_ptlo_bin_00y05;
  else if(absetabin==1) fit_ptlo_bin=fit_ptlo_bin_05y10;
  else if(absetabin==2) fit_ptlo_bin=fit_ptlo_bin_10y15;
  else if(absetabin==3) fit_ptlo_bin=fit_ptlo_bin_15y20;
  
  if(absetabin==0)      fit_pthi_bin=fit_pthi_bin_00y05;
  else if(absetabin==1) fit_pthi_bin=fit_pthi_bin_05y10;
  else if(absetabin==2) fit_pthi_bin=fit_pthi_bin_10y15;
  else if(absetabin==3) fit_pthi_bin=fit_pthi_bin_15y20;
    
  const int N_JERfits=fit_pthi_bin-fit_ptlo_bin;//dont touch  

  double* fgaus_width_arr=NULL;
  int nwidths=-1;
  if(!useconstxrange){
    fgaus_width_arr=setWidthArr(absetabin, &nwidths);
    assert(N_JERfits==nwidths);
  }



  //int rebin_hrsp_arr_N=( (int) (100.*( ( (hgaus_xhi+difffromone) - (hgaus_xlo-difffromone) ) ) )+ 1 );//2x size of fit range.
  int rebin_hrsp_arr_N=( (int) (100.*( ( (hgaus_xhi) - (hgaus_xlo) ) ) )+ 1 );// size of fit range; resolution is way too confident
  if(debugMode)std::cout<<"rebin_hrsp_arr_N="<<rebin_hrsp_arr_N<<std::endl;
  double rebin_hrsp_arr[rebin_hrsp_arr_N]={0.};
  //for(int i=0; i<rebin_hrsp_arr_N;i++) rebin_hrsp_arr[i]=(hgaus_xlo-difffromone)+i*defbinsize;//2x size of fit range
  for(int i=0; i<rebin_hrsp_arr_N;i++) rebin_hrsp_arr[i]=(hgaus_xlo)+i*defbinsize;//size of fit range; resulution is way too confident
  if(debugMode) for(int i=0; i<rebin_hrsp_arr_N;i++) std::cout<<"rebin_hrsp_arr["<<i<<"]="<<rebin_hrsp_arr[i]<<std::endl;
  
  //figure out what radius/jetcollection we are looking at using the ppData filename
  std::size_t radPos=inFile_MC_dir.find("_ak")+strlen("_ak");
  const std::string radiusInt= inFile_MC_dir.substr( radPos,1 );
  const std::string radius="R"+radiusInt+"_";
  std::cout<<"radius string is = "<<radius<<std::endl;
  
  //make final jetType strings
  std::size_t jetTypePos=radPos+1;
  std::size_t jetsPos=inFile_MC_dir.find("Jets");;  
  const std::string jetType=inFile_MC_dir.substr( jetTypePos,(jetsPos-jetTypePos) );
  std::cout<<"jetType string is = "<<jetType<<std::endl;
  
  const std::string fullJetType="ak"+radiusInt+jetType;//"ak3PFJets";//find in dir name
  std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  //assert(false);

  // input file
  const std::string inFile_MC_name="/Py8_CUETP8M1_QCDjetAllPtBins_"+fullJetType+"-allFiles.root";  
  std::cout<<std::endl<<"input file dir "<< inputDir+inFile_MC_dir <<std::endl;
  std::cout<<"input file name in dir "<< inFile_MC_name <<std::endl;
  TFile *finPP=new TFile( (inputDir+inFile_MC_dir+inFile_MC_name).c_str() );
  
  // output file
  std::string out_name=outputDir+fullJetType+"_PY8JER_"+outputTag;
  std::string outFileName=out_name+".root";
  std::cout<<std::endl<<"opening output root file "<< outFileName<<std::endl;
  TFile *rootfout = new TFile(outFileName.c_str(),"RECREATE");   
  rootfout->cd();    
  TDirectory* fout_ptbin_JER=rootfout->mkdir("all_ptbins_JER");
  
  // fit + fit settings
  std::cout<<"setting up fitters..."<< std::endl;
  TVirtualFitter::SetDefaultFitter("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04); 
  
  TH1F *hrsp[nbins_pt_debug]={};  //input
  
  TH1F *hMean=NULL, *hSigma=NULL, *hSigmaMu=NULL; 
  TH1F *hMean_fit=NULL, *hSigma_fit=NULL, *hSigmaMu_fit=NULL;
  TH1F *hChi2NDF=NULL, *hChi2Prob=NULL;

  Double_t fit_ptlo=ptbins_debug[fit_ptlo_bin];
  Double_t fit_pthi=ptbins_debug[fit_pthi_bin];
  TF1 *fgaus=NULL; 
  
  if(!draw_hJER)
    std::cout<<std::endl<<std::endl<<"skipping hJER hists"<<std::endl<<std::endl;
  else {
    
    const std::string algname=algo+radiusInt+jetType; 
    std::cout<<"making hists for algorithm "<< algname<<"..."<<std::endl;
    
    // init mu, sigma/mu hists (params fromhistograms themselves)
    hMean = new TH1F( "hMean", Form( "#mu %s", algname.c_str()),    
		      nbins_pt_debug, ptbins_debug);
    
    hSigma = new TH1F( "hSigma", Form( "#sigma %s", algname.c_str()),    
		       nbins_pt_debug, ptbins_debug);
    
    hSigmaMu = new TH1F( "hSigmaMu", Form( "#sigma/#mu %s", algname.c_str()),    
			 nbins_pt_debug, ptbins_debug);
    
    // init mu, sigma/mu hists (params from gauss fits)    
    hMean_fit = new TH1F( "hMean_fit", Form( "Gauss Core Fit #mu %s", algname.c_str()),    
			  nbins_pt_debug, ptbins_debug);
    
    hSigma_fit = new TH1F( "hSigma_fit", Form( "Gauss Core Fit #sigma %s", algname.c_str()),    
			   nbins_pt_debug, ptbins_debug);

    hSigmaMu_fit = new TH1F( "hSigmaMu_fit", Form( "Gauss Core Fit #sigma/#mu %s", algname.c_str()),    
			   nbins_pt_debug, ptbins_debug);

//    hSigmaUp_fit = new TH1F( "hSigmaUp_fit", Form( "Fit Sigma %s Up Unc.", algname.c_str()),    
//			   nbins_pt_debug, ptbins_debug);
//
//    hSigmaDown_fit = new TH1F( "hSigmaDown_fit", Form( "Fit Sigma %s Down Unc.", algname.c_str()),    
//			   nbins_pt_debug, ptbins_debug);

    // QA for gauss core fits
    hChi2NDF = new TH1F("hChi2NDF", Form( "Fit Chi2NDF %s", algname.c_str()),    
			nbins_pt_debug, ptbins_debug); 
    
    hChi2Prob = new TH1F("hChi2Prob", Form( "Fit Chi2Prob %s", algname.c_str()),    
			 nbins_pt_debug, ptbins_debug); 
    
    // pt bin hJER-fit loop    
    if(debugMode)std::cout<<"nbins_pt="<<nbins_pt_debug<<std::endl;    
    int width_arr_ind=0;
    for(int ip=0; ip<nbins_pt_debug; ip++){    
      
      int ptbin_ip=(int)ptbins_debug[ip], ptbin_ip1=(int)ptbins_debug[ip+1];      
      std::string inputHistName="";

      if(doGenYBins)
	inputHistName="hJER_"+doJetID+"wJetID_genetabin"+std::to_string(absetabin)+"_"+std::to_string(ptbin_ip)+"ptbin"+std::to_string(ptbin_ip1);      
      else
	inputHistName="hJER_"+doJetID+"wJetID_etabin"+std::to_string(absetabin)+"_"+std::to_string(ptbin_ip)+"ptbin"+std::to_string(ptbin_ip1);            
      
      // open the input hist    
      hrsp[ip] = (TH1F*)finPP->Get( inputHistName.c_str() );              
      if(!hrsp[ip]){ 
	std::cout<<"no input hist named " <<  inputHistName<< ", continue..."<<std::endl;
	continue;}          
      else if(hrsp[ip]->GetEntries()<50){
      	std::cout<<"#entries= " <<  hrsp[ip]->GetEntries()<< ", too few, continue..."<<std::endl;
	continue;	      }
            
      //set name for later
      if(doGenYBins){
	hrsp[ip]->SetName(("hJER_1wJetID_geneta_"+
			   absetabins_str[absetabin]+"_to_"+absetabins_str[absetabin+1]+
			   "_genpt_"+std::to_string(ptbin_ip) + "_to_" + std::to_string(ptbin_ip1) ).c_str());}
      else {
	hrsp[ip]->SetName(("hJER_1wJetID_receta_"+
			   absetabins_str[absetabin]+"_to_"+absetabins_str[absetabin+1]+
			   "_genpt_"+std::to_string(ptbin_ip) + "_to_" + std::to_string(ptbin_ip1) ).c_str());}      
      
      //FYI debug
      if(debugMode)std::cout<<std::endl;  
      if(debugMode)std::cout<<"________________________________"<<std::endl;  
      if(debugMode)std::cout<<"pt range for bin: "<<ptbin_ip<<" - "<<ptbin_ip1<< " GeV "<<std::endl;    
      if(debugMode)std::cout<<"hist title = "<<hrsp[ip]->GetTitle()<<std::endl;
      if(debugMode)hrsp[ip]->Print("base");          
      if(debugMode)std::cout<<std::endl;  
      
      //rebin if desired
      if(rebinJER&&rebinFactor>1.)
	hrsp[ip]= (TH1F*)hrsp[ip]->TH1::Rebin( rebinFactor, (inputHistName+"_rebin5").c_str());	
      
      //scaling
      std::cout<<"h_integral ="   << hrsp[ip]->Integral()<<std::endl;    
      hrsp[ip]->Scale( 1./ hrsp[ip]->Integral() );    //scale to area of 1
      
      
      //rebin again to get mu, sigma/mu from hists directly
      TH1F* hrsp_ip_rebin=(TH1F*)hrsp[ip]->Clone();
      hrsp_ip_rebin=(TH1F*)(hrsp_ip_rebin->TH1::Rebin( (rebin_hrsp_arr_N-1), 
						       ((std::string)hrsp[ip]->GetName()+"_rebin").c_str() , 
						       rebin_hrsp_arr));
      hrsp_ip_rebin->TH1::ResetStats();//reset stats after rebinning; root won't do it automatically...
      if(debugMode)hrsp_ip_rebin->Print("base");
      
      //get some numbers
      float h_norm  = hrsp_ip_rebin->GetBinContent(hrsp_ip_rebin->GetMaximumBin());
      float h_enorm = hrsp_ip_rebin->GetBinError(hrsp_ip_rebin->GetMaximumBin());
      float h_mean  = hrsp_ip_rebin->GetMean();    //mu
      float h_emean = hrsp_ip_rebin->GetMeanError()*errfact;    //mu err 
      float h_sig   = hrsp_ip_rebin->GetStdDev();     // sigma
      float h_esig  = hrsp_ip_rebin->GetStdDevError()*errfact;     // sigma
      float h_sigmu   = h_sig/h_mean;     // sigma/mu
      float h_esigmu  = h_sigmu * sqrt (  
					pow( h_emean/h_mean,2) + pow( h_esig/h_sig,2) ) ; //sigma/mu err
      
      std::cout<<"h_norm ="   << h_norm<<std::endl;    
      std::cout<<"h_enorm ="  << h_enorm<<std::endl;    
      std::cout<<"h_mean ="   << h_mean<<std::endl;    
      std::cout<<"h_emean ="  << h_emean<<std::endl;    
      std::cout<<"h_sig ="    << h_sig<<std::endl;    
      std::cout<<"h_esig ="   << h_esig<<std::endl;    
      std::cout<<"h_sigmu ="    << h_sigmu<<std::endl;    
      std::cout<<"h_esigmu ="   << h_esigmu<<std::endl;    
      //      std::cout<<"h_integral ="   << hrsp[ip]->Integral()<<std::endl;    
      std::cout<<"h_integral ="   << hrsp_ip_rebin->Integral()<<std::endl;    

      
      //why write here...? where does this write to??!
      hrsp[ip]->Write();
      if(ip<fit_ptlo_bin )continue;
      if(ip>=fit_pthi_bin )continue;
      
      //set mu, sig/mu as calc form histogram bin content
      std::cout<<"setting bin content for hMean[ip]=["<< ip<<"] "<<std::endl;    
      hMean->SetBinContent (ip+1, h_mean);    
      hMean->SetBinError   (ip+1, h_emean);          
      hSigma->SetBinContent (ip+1, h_sig);    
      hSigma->SetBinError   (ip+1, h_esig);    
      hSigmaMu->SetBinContent (ip+1, h_sigmu);    
      hSigmaMu->SetBinError   (ip+1, h_esigmu);    
      
      //
      fout_ptbin_JER->cd();            //tdirectory... why here?
      //
      
      //now get mu, sigma/mu from the gauss core fits
      std::cout<<std::endl;
      std::cout<<"fitting gauss core for "<<hrsp[ip]->GetName()<<std::endl;  

      if(useconstxrange){
	fgaus = new TF1("fgaus","gaus", fgaus_xlo,fgaus_xhi);  //what's been typicall used.      
      }
      else {
	fgaus = new TF1("fgaus","gaus", 1.-fgaus_width_arr[width_arr_ind],1.+fgaus_width_arr[width_arr_ind]);  //what's been typicall used.      
	width_arr_ind++;
      }
      //fgaus->SetParLimits(0,0.98,1.02);    //normalization
      //fgaus->SetParLimits(1,0.95,1.05);   // mean
      //fgaus->SetParLimits(2,0.0,0.5);     // width
      
      bool fitfailed=true;
      fitfailed = ((bool)(hrsp[ip]->Fit(fgaus,"R")));                      
      if(fitfailed)
	std::cout<<std::endl<<"!!!WARNING FIT FAILED!!!"<<std::endl<<std::endl;
      else 
	std::cout<<std::endl<<"fit success"<<std::endl<<std::endl;
      
      float f_norm  = -1.;
      float f_enorm = -1.;
      float f_mean  = -1.;
      float f_emean = -1.;
      float f_sig   = -1.;
      float f_esig  = -1.;      
      float f_sigmu   = -1.;
      float f_esigmu  = -1.;      
      float chi2NDF = -1.;
      float chi2prob= -1.;      
      
      f_norm  = (fitfailed)  ? -1. : fgaus->GetParameter(0); // z = if(condition) then(?) <do this> else(:) <do this>   
      f_enorm = (fitfailed)  ? -1. : fgaus->GetParError(0)*errfact ; 
      f_mean  = (fitfailed)  ? -1. : fgaus->GetParameter(1); // z = if(condition) then(?) <do this> else(:) <do this>   
      f_emean = (fitfailed)  ? -1. : fgaus->GetParError(1)*errfact ; 
      f_sig   = (fitfailed)  ? -1. : fgaus->GetParameter(2);
      f_esig  = (fitfailed)  ? -1. : fgaus->GetParError(2)*errfact ;      
      f_sigmu   = (fitfailed)  ? -1. : f_sig/f_mean;    
      f_esigmu  = (fitfailed)  ? -1. : f_sigmu * sqrt( pow( (f_emean/f_mean), 2) + 
						       pow( (f_esig/f_sig  ), 2) 
						       )   ;

      std::cout<<"f_norm ="   << f_norm<<std::endl;    
      std::cout<<"f_enorm ="  << f_enorm<<std::endl;    
      std::cout<<"f_mean ="   << f_mean<<std::endl;    
      std::cout<<"f_emean ="  << f_emean<<std::endl;    
      std::cout<<"f_sig ="    << f_sig<<std::endl;    
      std::cout<<"f_esig ="   << f_esig<<std::endl;    
      std::cout<<"f_sigmu ="  << f_sigmu<<std::endl;    
      std::cout<<"f_esigmu =" << f_esigmu<<std::endl;    
      if(!useconstxrange)std::cout<<"TF1->Integral="<<fgaus->Integral(1.-fgaus_width_arr[width_arr_ind],1.+fgaus_width_arr[width_arr_ind],1.e-9)<<std::endl;
      
      chi2NDF = (fitfailed)  ? -1. : fgaus->GetChisquare()/((float)fgaus->GetNDF());
      chi2prob = (fitfailed) ? -1. : fgaus->GetProb();
      
      if(   (f_mean<0.)||
	    (f_emean<0.)||
	    (f_sig<0.)||
	    (f_esig<0.)||
	    (chi2NDF<0.)||
	    (chi2prob<0.) ) std::cout<<"WARNING mu and/or sigma/mu are negative because FIT FAILED!!!!"<<std::endl;
      


      std::cout<<"Chi2NDF ="   << chi2NDF<<std::endl;    
      std::cout<<"Chi2Prob ="   << chi2prob<<std::endl;          
      
      std::cout<<"setting bin content for fit hists ip="<< ip<<" "<<std::endl;    
      hMean_fit->SetBinContent (ip+1, f_mean);    
      hMean_fit->SetBinError   (ip+1, f_emean);    
      
      hSigma_fit->SetBinContent (ip+1, f_sig);    
      hSigma_fit->SetBinError   (ip+1, f_esig);    
      
      hSigmaMu_fit->SetBinContent (ip+1, f_sigmu);    
      hSigmaMu_fit->SetBinError   (ip+1, f_esigmu);    
      
      hChi2NDF->SetBinContent (ip+1, chi2NDF );
      hChi2NDF->SetBinError (ip+1, 1.e-05 );
      
      hChi2Prob->SetBinContent (ip+1, chi2prob );      
      
    }// end loop over ptbins to yield mu, sigma/mu values
    
    rootfout->cd();//don't wanna be in directory anymore
    
    
    std::cout<<std::endl<<"fitting JER v genpt hists ..."<<std::endl;        
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-02); //diff tolerance setting for sigma/mu v gen pT fit
    
    // fit, mu v. gen pt, sigma/mu v. genpt hists, where vals are from histgram bin contents    
    bool hsigmaFitStatus= true;    
    TF1 *SigmaFit_h=new TF1("SigmaFit_h","[0]+[1]/(pow(x,[2])+[3]*x)", 
			    fit_ptlo,
			    fit_pthi );       // I want this fit to end at 967 GeV (for entire range)        
    std::cout<<std::endl<<"fitting JER width v. pt... init. sigmaFitStatus="<<hsigmaFitStatus<<std::endl<<std::endl; 
    
    if(fitsigma){
      if(debugMode)hsigmaFitStatus=hSigma->Fit(SigmaFit_h,"MEVR");	    
      else hsigmaFitStatus=hSigma->Fit(SigmaFit_h,"MER");	
      if(debugMode)hsigmaFitStatus=hSigma->Fit(SigmaFit_h,"EVR");	    
      else hsigmaFitStatus=hSigma->Fit(SigmaFit_h,"ER");	
    }
    else if(fitsigmu){
      if(debugMode)hsigmaFitStatus=hSigmaMu->Fit(SigmaFit_h,"MEVR");	    
      else         hsigmaFitStatus=hSigmaMu->Fit(SigmaFit_h,"MER");	
      if(debugMode)hsigmaFitStatus=hSigmaMu->Fit(SigmaFit_h,"EVR");	    
      else         hsigmaFitStatus=hSigmaMu->Fit(SigmaFit_h,"ER");	
    }

    float h_fitParam_0    =SigmaFit_h->GetParameter(0);
    float h_fitParam_0_err=SigmaFit_h->GetParError(0);
    float h_fitParam_1    =SigmaFit_h->GetParameter(1);
    float h_fitParam_1_err=SigmaFit_h->GetParError(1);
    float h_fitParam_2    =SigmaFit_h->GetParameter(2);
    float h_fitParam_2_err=SigmaFit_h->GetParError(2); 
    float h_fitParam_3    =SigmaFit_h->GetParameter(3);
    float h_fitParam_3_err=SigmaFit_h->GetParError(3); 

    std::cout << "h_fitParam_0     = " << h_fitParam_0    << std::endl;
    std::cout << "h_fitParam_0_err = " << h_fitParam_0_err<< std::endl;
    std::cout << "%err h_fitParam0 = " << (h_fitParam_0_err/h_fitParam_0)*100.    << std::endl;
    std::cout << "h_fitParam_1     = " << h_fitParam_1    << std::endl;
    std::cout << "h_fitParam_1_err = " << h_fitParam_1_err<< std::endl;
    std::cout << "%err h_fitParam1 = " << (h_fitParam_1_err/h_fitParam_1)*100.    << std::endl;
    std::cout << "h_fitParam_2     = " << h_fitParam_2    << std::endl;
    std::cout << "h_fitParam_2_err =   " << h_fitParam_2_err<< std::endl;
    std::cout << "%err h_fitParam2 = " << (h_fitParam_2_err/h_fitParam_2)*100.    << std::endl;
    std::cout << "h_fitParam_3     = " << h_fitParam_3    << std::endl;
    std::cout << "h_fitParam_3_err =   " << h_fitParam_3_err<< std::endl;
    std::cout << "%err h_fitParam3 = " << (h_fitParam_3_err/h_fitParam_3)*100.    << std::endl;        
    
    //make sys err envelope TF1s
    TF1 *SigmaFit_sysup_h=new TF1("SigmaFit_sysup_h","[0]+[1]/(pow(x,[2])+[3]*x)", 
				  fit_ptlo,
				  fit_pthi );       
    SigmaFit_sysup_h->SetParameter(0, SigmaFit_h->GetParameter(0)+SigmaFit_h->GetParError(0));
    SigmaFit_sysup_h->SetParameter(1, SigmaFit_h->GetParameter(1)+SigmaFit_h->GetParError(1));
    SigmaFit_sysup_h->SetParameter(2, SigmaFit_h->GetParameter(2)-SigmaFit_h->GetParError(2));//par 2 negative, carefi;
    SigmaFit_sysup_h->SetParameter(3, SigmaFit_h->GetParameter(3)-SigmaFit_h->GetParError(3));
    
    
    TF1 *SigmaFit_sysdown_h=new TF1("SigmaFit_sysdown_h","[0]+[1]/(pow(x,[2])+[3]*x)", 
				    fit_ptlo,
				    fit_pthi );       
    SigmaFit_sysdown_h->SetParameter(0, SigmaFit_h->GetParameter(0)-SigmaFit_h->GetParError(0));
    SigmaFit_sysdown_h->SetParameter(1, SigmaFit_h->GetParameter(1)-SigmaFit_h->GetParError(1));
    SigmaFit_sysdown_h->SetParameter(2, SigmaFit_h->GetParameter(2)+SigmaFit_h->GetParError(2));//par 2 negative, carefi;
    SigmaFit_sysdown_h->SetParameter(3, SigmaFit_h->GetParameter(3)+SigmaFit_h->GetParError(3));
    
    //get histogram versions
    TH1F *SigmaFit_h_hist= ((TH1F*) ((TH1F*)SigmaFit_h->GetHistogram())->Clone("SigmaFit_h_hist"));
    TH1F *SigmaFit_sysup_h_hist= ((TH1F*) ((TH1F*)SigmaFit_sysup_h->GetHistogram())->Clone("SigmaFit_sysup_h_hist"));
    TH1F *SigmaFit_sysdown_h_hist= ((TH1F*) ((TH1F*)SigmaFit_sysdown_h->GetHistogram())->Clone("SigmaFit_sysdown_h_hist"));
    


    // fit, mu v. gen pt, sigma/mu v. genpt hists, where vals are from gauss core fits
    bool fsigmaFitStatus= true;
    TF1 *SigmaFit_f=new TF1("SigmaFit_f","[0]+[1]/(pow(x,[2])+[3]*x)",     // this fit func from raghav/HIN/2.76 TeV/5 TeV
			    fit_ptlo,
			    fit_pthi );       // I want this fit to end at 967 GeV (for entire range)   
    
    std::cout<<std::endl<<"fitting JER width v. pt using HIN/2.76 TeV fits... init. sigmaFitStatus="<<fsigmaFitStatus<<std::endl<<std::endl; 
    if(fitsigma){
      if(debugMode)fsigmaFitStatus = hSigma_fit->Fit(SigmaFit_f,"MEVR");	
      else fsigmaFitStatus         = hSigma_fit->Fit(SigmaFit_f,"MER");	
      if(debugMode)fsigmaFitStatus = hSigma_fit->Fit(SigmaFit_f,"EVR");	
      else fsigmaFitStatus         = hSigma_fit->Fit(SigmaFit_f,"ER");	
    }
    else if(fitsigmu){
      if(debugMode)fsigmaFitStatus = hSigmaMu_fit->Fit(SigmaFit_f,"MEVR");	
      else fsigmaFitStatus         = hSigmaMu_fit->Fit(SigmaFit_f,"MER");	
      if(debugMode)fsigmaFitStatus = hSigmaMu_fit->Fit(SigmaFit_f,"EVR");	
      else fsigmaFitStatus         = hSigmaMu_fit->Fit(SigmaFit_f,"ER");	
    }
    else assert(false);

    bool fsigmaFit13tevStatus=true;
    TF1 *SigmaFit_f_13tev=new TF1("SigmaFit_f_13tev","[0]+([1]/(pow(x,[2])+[3]*x))+[4]*x",     // this fit func 13 TeV
				  fit_ptlo,
				  fit_pthi );       // I want this fit to end at 967 GeV (for entire range)    
    

    std::cout<<std::endl<<"fitting JER width v. pt using SMP/13 TeV fits... init. sigmaFitStatus="<<fsigmaFit13tevStatus<<std::endl<<std::endl; 
    if(fitsigma){
      if(debugMode)fsigmaFit13tevStatus = hSigma_fit->Fit(SigmaFit_f_13tev,"MEVR");	
      else fsigmaFit13tevStatus         = hSigma_fit->Fit(SigmaFit_f_13tev,"MER");	
      if(debugMode)fsigmaFit13tevStatus = hSigma_fit->Fit(SigmaFit_f_13tev,"EVR");	
      else fsigmaFit13tevStatus         = hSigma_fit->Fit(SigmaFit_f_13tev,"ER");	
    }
    else if (fitsigmu){
      if(debugMode)fsigmaFit13tevStatus = hSigmaMu_fit->Fit(SigmaFit_f_13tev,"MEVR");	
      else fsigmaFit13tevStatus         = hSigmaMu_fit->Fit(SigmaFit_f_13tev,"MER");	
      if(debugMode)fsigmaFit13tevStatus = hSigmaMu_fit->Fit(SigmaFit_f_13tev,"EVR");	
      else fsigmaFit13tevStatus         = hSigmaMu_fit->Fit(SigmaFit_f_13tev,"ER");	
    }
    else assert(false);
  

    bool fsigmaFit8tevStatus=true;
    TF1 *SigmaFit_f_8tev=new TF1("SigmaFit_f_8tev","sqrt( pow( [0]/x, 2) + pow([1],2)/x + pow([0],2)) ",     // this fit func 8 TeV
				 fit_ptlo,
				 fit_pthi );       // I want this fit to end at 967 GeV (for entire range)    
    
    std::cout<<std::endl<<"fitting JER width v. pt using SMP/13 TeV fits... init. sigmaFitStatus="<<fsigmaFit8tevStatus<<std::endl<<std::endl; 
    if(fitsigma){
      if(debugMode)fsigmaFit8tevStatus = hSigma_fit->Fit(SigmaFit_f_8tev,"MEVR");	
      else fsigmaFit8tevStatus         = hSigma_fit->Fit(SigmaFit_f_8tev,"MER");	
      if(debugMode)fsigmaFit8tevStatus = hSigma_fit->Fit(SigmaFit_f_8tev,"EVR");	
      else fsigmaFit8tevStatus         = hSigma_fit->Fit(SigmaFit_f_8tev,"ER");	
    }
    else if (fitsigmu) {
      if(debugMode)fsigmaFit8tevStatus = hSigmaMu_fit->Fit(SigmaFit_f_8tev,"MEVR");	
      else fsigmaFit8tevStatus         = hSigmaMu_fit->Fit(SigmaFit_f_8tev,"MER");	
      if(debugMode)fsigmaFit8tevStatus = hSigmaMu_fit->Fit(SigmaFit_f_8tev,"EVR");	
      else fsigmaFit8tevStatus         = hSigmaMu_fit->Fit(SigmaFit_f_8tev,"ER");	
    }
    else assert(false);
    
    //make sigmaup/sigmadown hists here
    
    //    TF1 *SigmaFit_sysupv2_f=new TF1("SigmaFit_sysupv2_f","[0]+[1]/(pow(x,[2])+[3]*x)", 
    //				    fit_ptlo,
    //				    fit_pthi );       // I want this fit to end at 967 GeV (for entire range)    
    //    hSigmaUp_fit->Fit(SigmaFit_sysupv2_f, "EVR");
    //    hSigmaUp_fit->Fit(SigmaFit_sysupv2_f, "ER");
    //    
    //    TF1 *SigmaFit_sysdownv2_f=new TF1("SigmaFit_sysdownv2_f","[0]+[1]/(pow(x,[2])+[3]*x)", 
    //				      fit_ptlo,
    //				      fit_pthi );       // I want this fit to end at 967 GeV (for entire range)    
    //    hSigmaDown_fit->Fit(SigmaFit_sysdownv2_f, "EVR");
    //    hSigmaDown_fit->Fit(SigmaFit_sysdownv2_f, "ER");
    
    TF1 *SigmaFit_sysup_f=new TF1("SigmaFit_sysup_f","[0]+[1]/(pow(x,[2])+[3]*x)", 
				  fit_ptlo,
				  fit_pthi );       
    SigmaFit_sysup_f->SetParameter(0, SigmaFit_f->GetParameter(0)+SigmaFit_f->GetParError(0));
    SigmaFit_sysup_f->SetParameter(1, SigmaFit_f->GetParameter(1)+SigmaFit_f->GetParError(1));
    SigmaFit_sysup_f->SetParameter(2, SigmaFit_f->GetParameter(2)-SigmaFit_f->GetParError(2));//par 2 negative, carefi;
    SigmaFit_sysup_f->SetParameter(3, SigmaFit_f->GetParameter(3)-SigmaFit_f->GetParError(3));

    
    TF1 *SigmaFit_sysdown_f=new TF1("SigmaFit_sysdown_f","[0]+[1]/(pow(x,[2])+[3]*x)", 
				    fit_ptlo,
				    fit_pthi );       
    SigmaFit_sysdown_f->SetParameter(0, SigmaFit_f->GetParameter(0)-SigmaFit_f->GetParError(0));
    SigmaFit_sysdown_f->SetParameter(1, SigmaFit_f->GetParameter(1)-SigmaFit_f->GetParError(1));
    SigmaFit_sysdown_f->SetParameter(2, SigmaFit_f->GetParameter(2)+SigmaFit_f->GetParError(2));//par 2 negative, carefi;
    SigmaFit_sysdown_f->SetParameter(3, SigmaFit_f->GetParameter(3)+SigmaFit_f->GetParError(3));

    TH1F *SigmaFit_f_hist= ((TH1F*) ((TH1F*)SigmaFit_f->GetHistogram())->Clone("SigmaFit_f_hist"));
    TH1F *SigmaFit_sysup_f_hist= ((TH1F*) ((TH1F*)SigmaFit_sysup_f->GetHistogram())->Clone("SigmaFit_sysup_f_hist"));
    TH1F *SigmaFit_sysdown_f_hist= ((TH1F*) ((TH1F*)SigmaFit_sysdown_f->GetHistogram())->Clone("SigmaFit_sysdown_f_hist"));
    //    TH1F *SigmaFit_sysupv2_f_hist= ((TH1F*) ((TH1F*)SigmaFit_sysupv2_f->GetHistogram())->Clone("SigmaFit_sysupv2_f_hist"));
    //    TH1F *SigmaFit_sysdownv2_f_hist= ((TH1F*) ((TH1F*)SigmaFit_sysdownv2_f->GetHistogram())->Clone("SigmaFit_sysdownv2_f_hist"));
    
    
    float f_fitParam_0    =SigmaFit_f->GetParameter(0);
    float f_fitParam_0_err=SigmaFit_f->GetParError(0);
    float f_fitParam_1    =SigmaFit_f->GetParameter(1);
    float f_fitParam_1_err=SigmaFit_f->GetParError(1);
    float f_fitParam_2    =SigmaFit_f->GetParameter(2);
    float f_fitParam_2_err=SigmaFit_f->GetParError(2); 
    float f_fitParam_3    =SigmaFit_f->GetParameter(3);
    float f_fitParam_3_err=SigmaFit_f->GetParError(3); 

    std::cout << "f_fitParam_0     = " << f_fitParam_0    << std::endl;
    std::cout << "f_fitParam_0_err = " << f_fitParam_0_err<< std::endl;
    std::cout << "%err f_fitParam0 = " << (f_fitParam_0_err/f_fitParam_0)*100.    << std::endl;
    std::cout << "f_fitParam_1     = " << f_fitParam_1    << std::endl;
    std::cout << "f_fitParam_1_err = " << f_fitParam_1_err<< std::endl;
    std::cout << "%err f_fitParam1 = " << (f_fitParam_1_err/f_fitParam_1)*100.    << std::endl;
    std::cout << "f_fitParam_2     = " << f_fitParam_2    << std::endl;
    std::cout << "f_fitParam_2_err =   " << f_fitParam_2_err<< std::endl;
    std::cout << "%err f_fitParam2 = " << (f_fitParam_2_err/f_fitParam_2)*100.    << std::endl;
    std::cout << "f_fitParam_3     = " << f_fitParam_3    << std::endl;
    std::cout << "f_fitParam_3_err =   " << f_fitParam_3_err<< std::endl;
    std::cout << "%err f_fitParam3 = " << (f_fitParam_3_err/f_fitParam_3)*100.    << std::endl;        

    
    //__________________________
    // DRAW THOSE PDFS //        
    std::string thePDFFileName=out_name+".pdf";    
    std::string open_thePDFFileName=thePDFFileName+"[";    
    std::string close_thePDFFileName=thePDFFileName+"]";    
    TCanvas* pdfoutCanv=new TCanvas("pdfoutCanv","pdfoutCanv", CANVX, CANVY);    
    pdfoutCanv->Print(open_thePDFFileName.c_str());            
    
    
    
    //___________________________________________________________________________
    //drawing the mu, sigma, sigma/mu, canvas
    TCanvas* muSigmaCanv=new TCanvas("muSigmaCanv","muSigmaCanv", CANVX, CANVY);    
    //muSigmaCanv->Divide(1,2);
    muSigmaCanv->Divide(1,3);
    muSigmaCanv->cd();    

    
    //float xmin=ptbins_debug[0], xmax=ptbins_debug[fit_pthi_bin];
    float xmin=ptbins_debug[fit_ptlo_bin-1], xmax=ptbins_debug[fit_pthi_bin];
    if(debugMode){
      std::cout<<"xmin="<<xmin<<std::endl;
      std::cout<<"xmax="<<xmax<<std::endl;    }

    TLine* lineAtOne=new TLine(xmin,1.,ptbins_debug[fit_pthi_bin+1],1.);    
    lineAtOne->SetLineStyle(2);  
    lineAtOne->SetLineColor(kBlack);
    TLine* lineAtTen=new TLine(xmin,10.,ptbins_debug[fit_pthi_bin+1],10.);    
    lineAtTen->SetLineStyle(2);  
    lineAtTen->SetLineColor(kBlack);
    

    //get pointers to panels on canvs as TPads, draw on TPads
    TPad* p1=(TPad*)muSigmaCanv->cd(1);        
    p1->SetLogx(1);    
    p1->SetGridx(0);    
    p1->SetTopMargin(0.05);
    //p1->SetBottomMargin();
    p1->cd();
    
    TPad* p2=(TPad*)muSigmaCanv->cd(2);        
    p2->SetLogx(1);    	
    p2->SetGridx(0);    
    p2->SetTopMargin(0.);
    //p2->SetBottomMargin();
    p2->cd();   
    
    TPad* p3=(TPad*)muSigmaCanv->cd(3);        
    p3->SetLogx(1);    	
    p3->SetGridx(0);    
    p3->SetTopMargin(0.);
    //p3->SetBottomMargin();
    p3->cd();   
    
    
    //draw pad 1
    p1->cd();    
    MakeHistMu( (TH1F*)hMean,xmin,xmax);
    hMean->DrawClone("HIST ][E"); 
    
    MakeHistMu( (TH1F*)hMean_fit, xmin,xmax);   
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
    
    float posx=0.40, posy=0.90, size=0.045;
    TLatex* jetdesc=new TLatex(posx, posy, fullJetType.c_str()); 
    jetdesc->SetTextSize(size);
    jetdesc->SetNDC();
    jetdesc->Draw();
    TLatex* absetarange=new TLatex(posx, posy-0.05, (absetabins_str[absetabin] + " < #||{y} < " + absetabins_str[absetabin+1]).c_str());
    absetarange->SetTextSize(size);
    absetarange->SetNDC();
    absetarange->Draw();

   
    /////////////////draw pad 1


    //draw pad 2
    p2->cd();    
    
    MakeHistSigma( (TH1F*) hSigma ,xmin,xmax); 
    //MakeHistSigma( (TH1F*) hSigma_fit ,xmin,xmax); 
    hSigma->DrawClone("HIST ][E");
    
    //this is a hist created from the gauss fits in each bin
    hSigma_fit->SetMarkerStyle(kOpenStar);
    hSigma_fit->SetMarkerColor(kBlack);
    hSigma_fit->SetLineColor(kMagenta);
    hSigma_fit->DrawClone("HIST ][E SAME");
    



    if(fitsigma){
      //this is a TF1 of sigma from hists
      SigmaFit_h->SetLineColor(kBlue-7);
      //SigmaFit_h->Draw("SAME");
      
      //this is a TF1 of sigma from fits of hists
      SigmaFit_f->SetLineColor(kMagenta-2);
      SigmaFit_f->Draw("SAME");
      
      SigmaFit_f_13tev->SetLineColor(kGreen-2);
      //SigmaFit_f_13tev->Draw("SAME");
      
      SigmaFit_f_8tev->SetLineColor(kTeal);
      //SigmaFit_f_8tev->Draw("SAME");
    }
    
    TLegend* hSigma_Leg=new TLegend( 0.60, 0.50 , 0.90 , 0.97, NULL,"brNDC");
    hSigma_Leg->AddEntry( hSigma, "#sigma from Hists", "lp");
    
//    if(fitsigma){
//      
//      std::string SigmaFit_h_prob_str;
//      makeFitProbStr(SigmaFit_h, &SigmaFit_h_prob_str);
//      if(debugMode)std::cout<<"SigmaFit_h_prob_str="<<SigmaFit_h_prob_str<<std::endl;    
//      std::string SigmaFit_h_chi2_str;
//      makeFitChi2Str(SigmaFit_h, &SigmaFit_h_chi2_str);
//      if(debugMode)std::cout<<"SigmaFit_h_chi2_str="<<SigmaFit_h_chi2_str<<std::endl;    
//      std::string SigmaFit_h_ndf_str=std::to_string( SigmaFit_h->GetNDF() );
//      if(debugMode)std::cout<<"SigmaFit_h_ndf_str="<<SigmaFit_h_ndf_str<<std::endl;    
//      
//      hSigma_Leg->AddEntry( SigmaFit_h, "Fit of #sigma from Hists", "l");   
//      hSigma_Leg->AddEntry( (TObject*)0, ("#chi^{2}/NDF="+SigmaFit_h_chi2_str+"/"+SigmaFit_h_ndf_str).c_str() , "");   
//      hSigma_Leg->AddEntry( (TObject*)0, ("Prob="+SigmaFit_h_prob_str).c_str() , "");   
//    }
    
//    hSigma_Leg->AddEntry( (TObject*)0, "" , "");   
    hSigma_Leg->AddEntry( hSigma_fit, "#sigma from Gauss Fits", "lp");
    
    
    if(fitsigma){
      
      std::string SigmaFit_f_prob_str;
      makeFitProbStr(SigmaFit_f, &SigmaFit_f_prob_str);
      if(debugMode)std::cout<<"SigmaFit_f_prob_str="<<SigmaFit_f_prob_str<<std::endl;    
      std::string SigmaFit_f_chi2_str;
      makeFitChi2Str(SigmaFit_f, &SigmaFit_f_chi2_str);
      if(debugMode)std::cout<<"SigmaFit_f_chi2_str="<<SigmaFit_f_chi2_str<<std::endl;    
      std::string SigmaFit_f_ndf_str=std::to_string( SigmaFit_f->GetNDF() );
      if(debugMode)std::cout<<"SigmaFit_f_ndf_str="<<SigmaFit_f_ndf_str<<std::endl;    
      
      
      hSigma_Leg->AddEntry( SigmaFit_f, "Fit of #sigma from Gauss Fits", "l");
      hSigma_Leg->AddEntry( (TObject*)0, ("#chi^{2}/NDF="+SigmaFit_f_chi2_str+"/"+SigmaFit_f_ndf_str).c_str() , "");   
      hSigma_Leg->AddEntry( (TObject*)0, ("Prob="+SigmaFit_f_prob_str).c_str() , "");   
    }


    hSigma_Leg->SetBorderSize(0);    
    hSigma_Leg->SetFillStyle(0);        
    hSigma_Leg->Draw();
    /////////////////draw pad2    

    //draw pad 3
    p3->cd();    
    
    MakeHistSigmaMu( (TH1F*) hSigmaMu ,xmin,xmax); 
    //MakeHistSigma( (TH1F*) hSigma_fit ,xmin,xmax); 
    hSigmaMu->DrawClone("HIST ][E");
    
    //this is a hist created from the gauss fits in each bin
    hSigmaMu_fit->SetMarkerStyle(kOpenStar);
    hSigmaMu_fit->SetMarkerColor(kBlack);
    hSigmaMu_fit->SetLineColor(kMagenta);
    hSigmaMu_fit->DrawClone("HIST ][E SAME");
    
    if(fitsigmu){
      
      //this is a TF1 of sigma/mu from hist
      SigmaFit_h->SetLineColor(kBlue-7);
      //SigmaFit_h->Draw("SAME");
      
      //this is a TF1 of sigma/mu from fits       
      SigmaFit_f->SetLineColor(kMagenta-2);
      SigmaFit_f->Draw("SAME");
      
      SigmaFit_f_13tev->SetLineColor(kGreen-2);
      //SigmaFit_f_13tev->Draw("SAME");
      
      SigmaFit_f_8tev->SetLineColor(kTeal);
      //SigmaFit_f_8tev->Draw("SAME");
    }
    
    TLegend* hSigmaMu_Leg=new TLegend( 0.60, 0.50 , 0.90 , 0.97, NULL,"brNDC");
    hSigmaMu_Leg->AddEntry( hSigma, "#sigma/#mu from Hists", "lp");
    
//    if(fitsigmu){
//      std::string SigmaFit_h_prob_str;
//      makeFitProbStr(SigmaFit_h, &SigmaFit_h_prob_str);
//      if(debugMode)std::cout<<"SigmaFit_h_prob_str="<<SigmaFit_h_prob_str<<std::endl;    
//      std::string SigmaFit_h_chi2_str;
//      makeFitChi2Str(SigmaFit_h, &SigmaFit_h_chi2_str);
//      if(debugMode)std::cout<<"SigmaFit_h_chi2_str="<<SigmaFit_h_chi2_str<<std::endl;    
//      std::string SigmaFit_h_ndf_str=std::to_string( SigmaFit_h->GetNDF() );
//      if(debugMode)std::cout<<"SigmaFit_h_ndf_str="<<SigmaFit_h_ndf_str<<std::endl;          
//      hSigmaMu_Leg->AddEntry( SigmaFit_h, "Fit of #sigma/#mu from Hists", "l");   
//      hSigmaMu_Leg->AddEntry( (TObject*)0, ("#chi^{2}/NDF="+SigmaFit_h_chi2_str+"/"+SigmaFit_h_ndf_str).c_str() , "");   
//      hSigmaMu_Leg->AddEntry( (TObject*)0, ("Prob="+SigmaFit_h_prob_str).c_str() , "");   
//    }    

//    hSigmaMu_Leg->AddEntry( (TObject*)0, "" , "");       
    hSigmaMu_Leg->AddEntry( hSigma_fit, "#sigma/#mu from Gauss Fits", "lp");
    
    if(fitsigmu){
      std::string SigmaFit_f_prob_str;
      makeFitProbStr(SigmaFit_f, &SigmaFit_f_prob_str);
      if(debugMode)std::cout<<"SigmaFit_f_prob_str="<<SigmaFit_f_prob_str<<std::endl;    
      std::string SigmaFit_f_chi2_str;
      makeFitChi2Str(SigmaFit_f, &SigmaFit_f_chi2_str);
      if(debugMode)std::cout<<"SigmaFit_f_chi2_str="<<SigmaFit_f_chi2_str<<std::endl;    
      std::string SigmaFit_f_ndf_str=std::to_string( SigmaFit_f->GetNDF() );
      if(debugMode)std::cout<<"SigmaFit_f_ndf_str="<<SigmaFit_f_ndf_str<<std::endl;          

      hSigmaMu_Leg->AddEntry( SigmaFit_f, "Fit of #sigma/#mu from Gauss Fits", "l");
      hSigmaMu_Leg->AddEntry( (TObject*)0, ("#chi^{2}/NDF="+SigmaFit_f_chi2_str+"/"+SigmaFit_f_ndf_str).c_str() , "");   
      hSigmaMu_Leg->AddEntry( (TObject*)0, ("Prob="+SigmaFit_f_prob_str).c_str() , "");   
    }
    
    
    hSigmaMu_Leg->SetBorderSize(0);    
    hSigmaMu_Leg->SetFillStyle(0);        
    hSigmaMu_Leg->Draw();
    ///////////////////draw pad3
    
    muSigmaCanv->Print(thePDFFileName.c_str());            
    
    rootfout->cd();
    muSigmaCanv->SetTitle("#mu, #sigma, #sigma/#mu Canvas");            muSigmaCanv->Write("mean_sigma_sigmu_canv");            
    


    //____________________________________________________________
    //FITS WITH FIT UP/DOWN ERR BAND
    TCanvas* fits_wSystCanv=new TCanvas("fits_wSystCanv", "sigma/mu fits w/ syst errs", CANVX, CANVY);
    fits_wSystCanv->Divide(1,2);
    
    //pad 1
    fits_wSystCanv->cd(1)->SetLogx(1);
    fits_wSystCanv->cd(1)->SetLogy(0);    
    fits_wSystCanv->cd(1);

    if      (fitsigma)    MakeHistSigma(  SigmaFit_sysup_h_hist, xmin, xmax);
    else if (fitsigmu)    MakeHistSigmaMu(SigmaFit_sysup_h_hist, xmin, xmax);

    SigmaFit_sysup_h_hist->SetAxisRange(0.04, 0.16,"Y");

    if(fitsigma) SigmaFit_sysup_h_hist->SetTitle("#sigma [Bin Content] v. GEN Jet p_{T} Fit");
    else if(fitsigmu) SigmaFit_sysup_h_hist->SetTitle("(#sigma/#mu) [Bin Content] v. GEN Jet p_{T} Fit");
    SigmaFit_sysup_h_hist->GetYaxis()->SetLabelSize(0.04);
    SigmaFit_sysup_h_hist->GetXaxis()->SetLabelSize(0.04);
    SigmaFit_h_hist->SetLineColor(kBlue-7);
    
    SigmaFit_sysup_h_hist->Draw("HIST ][");
    SigmaFit_h_hist->Draw("HIST ][ SAME");
    SigmaFit_sysdown_h_hist->Draw("HIST ][ SAME");
    
    TLegend* legend_h_fitsys=new TLegend(0.6,0.6,0.9,0.8);
    legend_h_fitsys->SetFillStyle(0);
    legend_h_fitsys->SetBorderSize(0);
    legend_h_fitsys->AddEntry(SigmaFit_h_hist,"Fit from Hist Statistics", "l");
    legend_h_fitsys->AddEntry(SigmaFit_sysup_h_hist,"Upper/Lower Unc.", "l");
    legend_h_fitsys->Draw();
    
    //pad 2
    fits_wSystCanv->cd(2)->SetLogx(1);
    fits_wSystCanv->cd(2)->SetLogy(0);
    fits_wSystCanv->cd(2);    
    
    if     (fitsigma)MakeHistSigma(  SigmaFit_sysup_f_hist, xmin, xmax);
    else if(fitsigmu)MakeHistSigmaMu(SigmaFit_sysup_f_hist, xmin, xmax);
    
    SigmaFit_sysup_f_hist->SetAxisRange(0.04, 0.20,"Y");

    if(fitsigma)SigmaFit_sysup_f_hist->SetTitle("#sigma [Gauss Core] v. GEN Jet p_{T} Fit");
    else if(fitsigmu)SigmaFit_sysup_f_hist->SetTitle("(#sigma/#mu) [Gauss Core] v. GEN Jet p_{T} Fit");
    SigmaFit_sysup_f_hist->GetYaxis()->SetLabelSize(0.04);
    SigmaFit_sysup_f_hist->GetXaxis()->SetLabelSize(0.04);
    SigmaFit_f_hist->SetLineColor(kMagenta-2);
    
    SigmaFit_sysup_f_hist->Draw("HIST ][");    //    SigmaFit_sysupv2_f_hist->Draw("HIST ][ SAME");
    SigmaFit_f_hist->Draw("HIST ][ SAME");
    SigmaFit_sysdown_f_hist->Draw("HIST ][ SAME");    //    SigmaFit_sysdownv2_f_hist->Draw("HIST ][ SAME");

    TLegend* legend_f_fitsys=new TLegend(0.6,0.6,0.9,0.8);
    legend_f_fitsys->SetFillStyle(0);
    legend_f_fitsys->SetBorderSize(0);
    legend_f_fitsys->AddEntry(SigmaFit_f_hist,"Fit from Gaussian Cores", "l");
    legend_f_fitsys->AddEntry(SigmaFit_sysup_f_hist,"Upper/Lower Unc.", "l");
    legend_f_fitsys->Draw();
    
    //print+write
    fits_wSystCanv->Print(thePDFFileName.c_str());
    rootfout->cd();
    fits_wSystCanv->SetTitle("Jet Resolution Fits w/ Sys Err Band"); 
    fits_wSystCanv->Write("resolutionTF1_hists_wsyst_canv");
    
      
    
    //___________________________________________
    // GAUSS CORE CHI2 QA
    TCanvas* gaussCoreChi2Canv=new TCanvas("gaussCoreChi2Canv", "gaussCoreChi2Canv",CANVX,CANVY);
    gaussCoreChi2Canv->Divide(1,2);
    //draw pad 1
    TPad* chi2p1=(TPad*)gaussCoreChi2Canv->cd(1);        
    chi2p1->SetLogx(1);    	
    chi2p1->SetLogy(1);    	
    chi2p1->SetGridx(0);
    chi2p1->SetTopMargin(0.);
    chi2p1->cd();        
    
    TPad* chi2p2=(TPad*)gaussCoreChi2Canv->cd(2);        
    chi2p2->SetLogx(1);    	
    chi2p2->SetLogy(1);    	
    chi2p2->SetGridx(0);
    chi2p2->cd();        

    ////draw pad 1
    chi2p1->cd();
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
    
    ////draw pad 2
    chi2p2->cd();    
    MakeHistChi2Prob( (TH1F*)hChi2Prob,xmin,xmax);           
    hChi2Prob->DrawClone("HIST ][E");    
    
    gaussCoreChi2Canv->Print(thePDFFileName.c_str());
    rootfout->cd();
    gaussCoreChi2Canv->SetTitle("Gauss Core Fit QA, #Chi^{2} and Probability");
    gaussCoreChi2Canv->Write("gausscore_chi2_canv");





    //___________________________________________
    //CANVASES OF INDIVIDUAL JER FITS
    std::cout<<"making+writing canvases of JER reco/gen hists"<<std::endl;
    
    TCanvas* JERfitsCanv=new TCanvas("JERfitsCanv","Canv for JER fit", CANVX, CANVY);    
    JERfitsCanv->SetLogy(1);
    JERfitsCanv->cd();                    
    for(int j=0;j<nbins_pt_debug;++j){ 
      int ptbin_j=(int)ptbins_debug[j];
      int ptbin_j1=(int)ptbins_debug[j+1];
      std::string hrspTitle=std::to_string(ptbin_j)+" GeV < GEN Jet p_{T} < "+std::to_string(ptbin_j1)+" GeV";    
      std::string hrsp_XAxTitle="RECO p_{T}/GEN p_{T}";        
      makeJERSHists( ( TCanvas*    ) JERfitsCanv ,  
		     ( TH1F*       ) hrsp[j]          , 
		     ( std::string ) hrspTitle        ,
		     ( std::string ) hrsp_XAxTitle    ,
		     ( TFile*      ) rootfout         ,
		     ( TDirectory*  ) fout_ptbin_JER );    }

    
    //___________________________________________
    //MULTIPANEL CANVASES OF JER HISTS
    std::cout<<"making+writing multi-panel canvases of JER reco/gen hists"<<std::endl;
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
    pdfoutCanv->Print( close_thePDFFileName.c_str() );        
    
    rootfout->cd();
    
    
    
    
    //FIT + HIST WRITING GOES HERE; NO CANVASES
    
    //_______________________________________________    
    //mu + sigma/mu plots, vals from histogram contents
    hMean->Write();                   hSigma->Write();               
    //fits
    SigmaFit_h->Write();    SigmaFit_sysup_h->Write();    SigmaFit_sysdown_h->Write();
    //histogram versions of fits for viewing on my local machine, which CRASHES b.c. of TF1 for some reason
    SigmaFit_h_hist->Write();    SigmaFit_sysup_h_hist->Write();    SigmaFit_sysdown_h_hist->Write();
    

    //_______________________________________________
    //mu + sigma/mu plots + chi2/ndf, vals from gauss core fits 
    hMean_fit->Write();    hSigma_fit->Write();    hChi2NDF->Write();                         
    //fits w/ sys err envelope
    SigmaFit_f->Write();    SigmaFit_sysup_f->Write();    SigmaFit_sysdown_f->Write();    
    //histogram versions of fits for viewing on my local machine, which CRASHES b.c. of TF1 for some reason
    SigmaFit_f_hist->Write();    SigmaFit_sysup_f_hist->Write();    SigmaFit_sysdown_f_hist->Write();
    
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
  
  
  
  if(debugMode)std::cout<<"closing output root file "<<rootfout->GetName()<<std::endl;
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
  if(argc==7)
    rStatus=printPlots_ppMC_JERS_etabin((const std::string)argv[1], (const std::string)argv[2] , (const bool)std::atoi(argv[3]), (const bool)std::atoi(argv[4]), (const int) std::atoi(argv[5]) , (const std::string)argv[6] );
  if(argc==6)
    rStatus=printPlots_ppMC_JERS_etabin((const std::string)argv[1], (const std::string)argv[2] , (const bool)std::atoi(argv[3]), (const bool)std::atoi(argv[4]), (const int) std::atoi(argv[5]) );
  if(argc==5)
    rStatus=printPlots_ppMC_JERS_etabin((const std::string)argv[1], (const std::string)argv[2] , (const bool)std::atoi(argv[3]), (const bool)std::atoi(argv[4]), -1);
  return rStatus;
}
