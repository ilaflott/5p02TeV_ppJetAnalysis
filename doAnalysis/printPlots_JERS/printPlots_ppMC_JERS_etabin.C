#include "printPlots.h"

const bool fitJohnsDSCB=true;//usually true
const bool drawJohnsDSCB=fitJohnsDSCB&&false;
//const bool drawAltDSCBFits=fitJohnsDSCB&&false;
//const std::string dscb_dir="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_10_2_0_pre6/src/doAnalysis/printPlots_JERS/dbc5tev_myver/output_07.22.21";
//const std::string dscb_finstr_start="/ak4PF_PY8JER_";//??eta??
//const std::string dscb_finstr_end  ="_07.22.21_minloptcuts_wPFMET_DCSB_geny.root";//_07.22.21_sigmu_wPFMET_minlopTcuts_genbinsonly_semifinalv4.root";
const std::string dscb_dir="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_10_2_0_pre6/src/doAnalysis/printPlots_JERS/dbc5tev_myver/output_09.27.21";
const std::string dscb_finstr_start="/ak4PF_PY8JER_";//??eta??
//const std::string dscb_finstr_end  ="_09.27.21_RECOptLo10_GENptLo3_TESTFIT2_GENYbinEQYbin_wPFMET_DCSB_geny.root";//_07.22.21_sigmu_wPFMET_minlopTcuts_genbinsonly_semifinalv4.root";
const std::string dscb_finstr_end  ="_09.27.21_RECOptLo10_GENptLo3_TESTFIT_GENYbinEQYbin_wPFMET_DCSB_geny.root";//_07.22.21_sigmu_wPFMET_minlopTcuts_genbinsonly_semifinalv4.root";
const std::string hjohn_dscb_mu   ="hJohnDSCBMean";
const std::string hjohn_dscb_sigma="hJohnDSCBSigma";

//entire macro settings
const bool debugMode=false;
const int CANVX=1200, CANVY=1200;
const std::string doJetID="1";

//// JER settings
//hrsp
//const bool doGenYBins=true;//RapBins -> dual-diff xsec bins, GenBins -> variable, depends on readForests
//const bool doGenEtaBins=false;//RapBins -> dual-diff xsec bins, GenBins -> variable, depends on readForests
const bool rebinJER=true;
const int rebinFactor=2;// test
const double hgaus_xlo=0.70, hgaus_xhi=1.+(1.-hgaus_xlo);//test
const bool useconstxrange=false;
const double fgaus_xlo=0.80, fgaus_xhi=1.+(1.-fgaus_xlo);//test
const double defbinsize=0.01;//dont touch
//resolution v genpt fit options
const double errfact=1.0;// test



//MCEff settings
const bool drawProfiles = false;//rarely used




int printPlots_ppMC_JERS_etabin(std::string inFile_MC_dir,const std::string outputTag, const bool draw_hJER=true, const bool draw_MCEff=false, 
				const int absetabin=0, const std::string fitquantity="sigma", const std::string bintype="y"){
  
  if      (bintype=="geny"  )std::cout<<"fitting JER as a function of gen jet rapidity."      <<std::endl;
  else if (bintype=="geneta")std::cout<<"fitting JER as a function of gen jet pseudorapidity."<<std::endl;
  else if (bintype=="y"  )std::cout<<"fitting JER as a function of reco jet rapidity."      <<std::endl;
  else if (bintype=="eta")std::cout<<"fitting JER as a function of reco jet pseudorapidity."<<std::endl;
  else assert(false);
  
  std::cout<<"absetabin=="<<absetabin<<std::endl;

  //if(absetabin!=-1) doabsetabinJER=true;
  if(draw_MCEff)std::cout<<"draw_MCEff is true but still not drawing them in this macro because that's stupid. "<<std::endl;
  else if(!draw_MCEff)std::cout<<"draw_MCEff is false and you should take the MCEff plots out of this macro and put it into a separate macro some time ian. "<<std::endl;
  if(draw_hJER)std::cout<<"draw_hJER is true. duh."<<std::endl;
  else if(!draw_hJER){ std::cout<<"draw_hJER is false. dumb. exit"<<std::endl; assert(false);}


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

  const std::string algname=algo+radiusInt+jetType; 
  std::cout<<"making hists for algorithm "<< algname<<"..."<<std::endl;

  // root style settings.
  std::cout<<"forcing style"<<std::endl;
  LoadStyle();

  // general fit setting stuff for gauss core fits
  std::cout<<"setting up fitters..."<< std::endl;
  TVirtualFitter::SetDefaultFitter("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04); 
 
  bool fitsigma=false, fitsigmu=false;
  if(fitquantity=="sigma"){
    fitsigma=true; fitsigmu=false;}
  else if(fitquantity=="sigmu"){
    fitsigma=false; fitsigmu=true;}
  else assert(false);
  //const bool fitsigma=false, fitsigmu=!fitsigma;//fit either sigma, or sigma/mu

  //which bins to loop over for gauss cores and with what width arrays
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
  std::cout<<"fit_pthi_bin="<<fit_pthi_bin<<std::endl;
  std::cout<<"fit_ptlo_bin="<<fit_ptlo_bin<<std::endl;
  std::cout<<"N_JERfits="<<N_JERfits<<std::endl;
  
  double* fgaus_width_arr=NULL;
  //double* fgaus_width_y_arr=NULL;
  int nwidths=-1;
  if(!useconstxrange){
    if      (bintype=="eta" || bintype=="geneta") fgaus_width_arr = setWidthArr_eta(absetabin, &nwidths);//, N_JERfits);
    else if (bintype=="y"   || bintype=="geny"  ) fgaus_width_arr = setWidthArr_y  (absetabin, &nwidths);//, N_JERfits);
    if(N_JERfits!=nwidths){
      std::cout<<"error! N_JERfits="<<N_JERfits<<std::endl;
      std::cout<<"error! nwidths="<<nwidths<<std::endl;
      assert(false);
    }     
    
  }

  // fit + fit settings
  Double_t fit_ptlo=ptbins_debug[fit_ptlo_bin];
  Double_t fit_pthi=ptbins_debug[fit_pthi_bin];

  //int rebin_hrsp_arr_N=( (int) (100.*( ( (hgaus_xhi+difffromone) - (hgaus_xlo-difffromone) ) ) )+ 1 );//2x size of fit range.
  int rebin_hrsp_arr_N=( (int) (100.*( ( (hgaus_xhi) - (hgaus_xlo) ) ) )+ 1 );// size of fit range; resolution is way too confident
  if(debugMode)std::cout<<"rebin_hrsp_arr_N="<<rebin_hrsp_arr_N<<std::endl;
  double rebin_hrsp_arr[rebin_hrsp_arr_N]={0.};
  //for(int i=0; i<rebin_hrsp_arr_N;i++) rebin_hrsp_arr[i]=(hgaus_xlo-difffromone)+i*defbinsize;//2x size of fit range
  for(int i=0; i<rebin_hrsp_arr_N;i++) rebin_hrsp_arr[i]=(hgaus_xlo)+i*defbinsize;//size of fit range; resulution is way too confident
  if(debugMode)   for(int i=0; i<rebin_hrsp_arr_N;i++) std::cout<<"rebin_hrsp_arr["<<i<<"]="<<rebin_hrsp_arr[i]<<std::endl;
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

  //
  ////y bins
  //TF1 *fgaus_y=NULL; 
  //TH1F *hrsp_y[nbins_pt_debug]={};  //input
  //TH1F *hMean_y=NULL, *hSigma_y=NULL, *hSigmaMu_y=NULL; 
  //TH1F *hMean_fit_y=NULL, *hSigma_fit_y=NULL, *hSigmaMu_fit_y=NULL;
  //TH1F *hChi2NDF_y=NULL, *hChi2Prob_y=NULL;
  //
  //// init mu, sigma/mu hists (params fromhistograms themselves)
  //hMean_y = new TH1F( "hMean_y", Form( "#mu %s", algname.c_str()),    
  //		      nbins_pt_debug, ptbins_debug);
  //
  //hSigma_y = new TH1F( "hSigma_y", Form( "#sigma %s", algname.c_str()),    
  //		       nbins_pt_debug, ptbins_debug);
  //
  //hSigmaMu_y = new TH1F( "hSigmaMu_y", Form( "#sigma/#mu %s", algname.c_str()),    
  //			 nbins_pt_debug, ptbins_debug);
  //
  //// init mu, sigma/mu hists (params from gauss fits)    
  //hMean_fit_y = new TH1F( "hMean_fit_y", Form( "Gauss Core Fit #mu %s", algname.c_str()),    
  //			  nbins_pt_debug, ptbins_debug);
  //
  //hSigma_fit_y = new TH1F( "hSigma_fit_y", Form( "Gauss Core Fit #sigma %s", algname.c_str()),    
  //			   nbins_pt_debug, ptbins_debug);
  //
  //hSigmaMu_fit_y = new TH1F( "hSigmaMu_fit_y", Form( "Gauss Core Fit #sigma/#mu %s", algname.c_str()),    
  //			     nbins_pt_debug, ptbins_debug);
  //
  //// QA for gauss core fits
  //hChi2NDF_y = new TH1F("hChi2NDF_y", Form( "Fit Chi2NDF %s", algname.c_str()),    
  //			nbins_pt_debug, ptbins_debug); 
  //
  //hChi2Prob_y = new TH1F("hChi2Prob_y", Form( "Fit Chi2Prob %s", algname.c_str()),    
  //			 nbins_pt_debug, ptbins_debug); 
  //
  //// pt bin hJER-fit loop    
  //if(debugMode)std::cout<<"nbins_pt="<<nbins_pt_debug<<std::endl;    
  //int width_arr_y_ind=0;
  //
  ////|y| bin version
  //for(int ip=0; ip<nbins_pt_debug; ip++){    
  //  int ptbin_ip=(int)ptbins_debug[ip], ptbin_ip1=(int)ptbins_debug[ip+1];          
  //  
  //  std::string inputHistName_y="";
  //  inputHistName_y="hJER_"+doJetID+"wJetID_ybin"+std::to_string(absetabin)+"_"+std::to_string(ptbin_ip)+"ptbin"+std::to_string(ptbin_ip1);            
  //
  //  hrsp_y[ip] = (TH1F*)finPP->Get( inputHistName_y.c_str() );              
  //  if(!hrsp_y[ip]){ 
  //    std::cout<<"no input hist named " <<  inputHistName_y<< ", continue..."<<std::endl;
  //    continue;}          
  //  else if(hrsp_y[ip]->GetEntries()<50){
  //    std::cout<<"#entries= " <<  hrsp_y[ip]->GetEntries()<< ", too few, continue..."<<std::endl;
  //    continue;	      }
  //
  //  hrsp_y[ip]->SetName(("hJER_1wJetID_recy_"+
  //			 absetabins_str[absetabin]+"_to_"+absetabins_str[absetabin+1]+
  //			 "_genpt_"+std::to_string(ptbin_ip) + "_to_" + std::to_string(ptbin_ip1) ).c_str());      
  //  
  //  if(rebinJER&&rebinFactor>1.)
  //    hrsp_y[ip]= (TH1F*)hrsp_y[ip]->TH1::Rebin( rebinFactor, (inputHistName_y+"_rebin"+std::to_string(rebinFactor)).c_str());	
  //  
  //  
  //  std::cout<<"h_integral ="   << hrsp_y[ip]->Integral()<<std::endl;    
  //  hrsp_y[ip]->Scale( 1./ hrsp_y[ip]->Integral() );    //scale to area of 1
  //  
  //  
  //  TH1F* hrsp_ip_y_rebin=(TH1F*)hrsp_y[ip]->Clone();
  //  hrsp_ip_y_rebin=(TH1F*)(hrsp_ip_y_rebin->TH1::Rebin( (rebin_hrsp_arr_N-1), 
  //							 ((std::string)hrsp_y[ip]->GetName()+"_rebin").c_str() , 
  //							 rebin_hrsp_arr));
  //  hrsp_ip_y_rebin->TH1::ResetStats();//reset stats after rebinning; root won't do it automatically...
  //  if(debugMode)hrsp_ip_y_rebin->Print("base");
  //  hrsp_y[ip]->Write();
  //
  //  if(ip<fit_ptlo_bin )continue;
  //  if(ip>=fit_pthi_bin )continue;
  //
  //  //get some numbers
  //  float h_y_norm  = hrsp_ip_y_rebin->GetBinContent(hrsp_ip_y_rebin->GetMaximumBin());
  //  float h_y_enorm = hrsp_ip_y_rebin->GetBinError(hrsp_ip_y_rebin->GetMaximumBin());
  //  float h_y_mean  = hrsp_ip_y_rebin->GetMean();    //mu
  //  float h_y_emean = hrsp_ip_y_rebin->GetMeanError()*errfact;    //mu err 
  //  float h_y_sig   = hrsp_ip_y_rebin->GetStdDev();     // sigma
  //  float h_y_esig  = hrsp_ip_y_rebin->GetStdDevError()*errfact;     // sigma
  //  float h_y_sigmu   = h_y_sig/h_y_mean;     // sigma/mu
  //  float h_y_esigmu  = h_y_sigmu * sqrt (  
  //					  pow( h_y_emean/h_y_mean,2) + pow( h_y_esig/h_y_sig,2) ) ; //sigma/mu err
  //
  //  std::cout<<"h_y_norm ="     << h_y_norm<<std::endl;    
  //  std::cout<<"h_y_enorm ="    << h_y_enorm<<std::endl;    
  //  std::cout<<"h_y_mean ="     << h_y_mean<<std::endl;    
  //  std::cout<<"h_y_emean ="    << h_y_emean<<std::endl;    
  //  std::cout<<"h_y_sig ="      << h_y_sig<<std::endl;    
  //  std::cout<<"h_y_esig ="     << h_y_esig<<std::endl;    
  //  std::cout<<"h_y_sigmu ="    << h_y_sigmu<<std::endl;    
  //  std::cout<<"h_y_esigmu ="   << h_y_esigmu<<std::endl;    
  //  std::cout<<"h_y_integral =" << hrsp_ip_y_rebin->Integral()<<std::endl;    
  //
  //  //set mu, sig/mu as calc form histogram bin content
  //  std::cout<<"setting bin content for hMean[ip]=["<< ip<<"] "<<std::endl;    
  //  hMean_y   ->SetBinContent (ip+1, h_y_mean);    
  //  hMean_y   ->SetBinError   (ip+1, h_y_emean);          
  //  hSigma_y  ->SetBinContent (ip+1, h_y_sig);    
  //  hSigma_y  ->SetBinError   (ip+1, h_y_esig);    
  //  hSigmaMu_y->SetBinContent (ip+1, h_y_sigmu);    
  //  hSigmaMu_y->SetBinError   (ip+1, h_y_esigmu);    
  //  
  //  //
  //  fout_ptbin_JER->cd();            //tdirectory... why here?
  //  //
  //
  //  //now get mu, sigma/mu from the gauss core fits
  //  std::cout<<std::endl;
  //  std::cout<<"fitting gauss core for "<<hrsp_y[ip]->GetName()<<std::endl;  
  //  
  //  if(useconstxrange){
  //    fgaus_y = new TF1("fgaus_y","gaus_y", fgaus_xlo,fgaus_xhi);  //what's been typicall used.      
  //  }
  //  else {
  //    fgaus_y = new TF1("fgaus_y","gaus_y", 1.-fgaus_width_y_arr[width_arr_y_ind],1.+fgaus_width_y_arr[width_arr_y_ind]);  //what's been typicall used.      
  //    width_arr_y_ind++;
  //  }
  //  
  //  bool fitfailed_y=true;
  //  fitfailed_y = ((bool)(hrsp_y[ip]->Fit(fgaus_y,"R")));                      
  //  if(fitfailed_y)
  //    std::cout<<std::endl<<"!!!WARNING FIT FAILED FOR |y| BIN!!!"<<std::endl<<std::endl;
  //  else 
  //    std::cout<<std::endl<<"fit success for |y| bin"<<std::endl<<std::endl;
  //  
  //  float f_y_norm   = -1.;
  //  float f_y_enorm  = -1.;
  //  float f_y_mean   = -1.;
  //  float f_y_emean  = -1.;
  //  float f_y_sig    = -1.;
  //  float f_y_esig   = -1.;      
  //  float f_y_sigmu  = -1.;
  //  float f_y_esigmu = -1.;      
  //  float y_chi2NDF  = -1.;
  //  float y_chi2prob = -1.;      
  //  
  //  f_y_norm    = (fitfailed_y)  ? -1. : fgaus_y->GetParameter(0); // z = if(condition) then(?) <do this> else(:) <do this>   
  //  f_y_enorm   = (fitfailed_y)  ? -1. : fgaus_y->GetParError(0)*errfact ; 
  //  f_y_mean    = (fitfailed_y)  ? -1. : fgaus_y->GetParameter(1); // z = if(condition) then(?) <do this> else(:) <do this>   
  //  f_y_emean   = (fitfailed_y)  ? -1. : fgaus_y->GetParError(1)*errfact ; 
  //  f_y_sig     = (fitfailed_y)  ? -1. : fgaus_y->GetParameter(2);
  //  f_y_esig    = (fitfailed_y)  ? -1. : fgaus_y->GetParError(2)*errfact ;      
  //  f_y_sigmu   = (fitfailed_y)  ? -1. : f_y_sig/f_y_mean;    
  //  f_y_esigmu  = (fitfailed_y)  ? -1. : f_y_sigmu * sqrt( pow( (f_y_emean/f_y_mean), 2) + 
  //							   pow( (f_y_esig/f_y_sig  ), 2) 
  //							   )   ;
  //  
  //  std::cout<<"f_y_norm ="   << f_y_norm<<std::endl;    
  //  std::cout<<"f_y_enorm ="  << f_y_enorm<<std::endl;    
  //  std::cout<<"f_y_mean ="   << f_y_mean<<std::endl;    
  //  std::cout<<"f_y_emean ="  << f_y_emean<<std::endl;    
  //  std::cout<<"f_y_sig ="    << f_y_sig<<std::endl;    
  //  std::cout<<"f_y_esig ="   << f_y_esig<<std::endl;    
  //  std::cout<<"f_y_sigmu ="  << f_y_sigmu<<std::endl;    
  //  std::cout<<"f_y_esigmu =" << f_y_esigmu<<std::endl;    
  //  if(!useconstxrange)std::cout<<"TF1->Integral="<<fgaus_y->Integral(1.-fgaus_width_y_arr[width_arr_y_ind],1.+fgaus_width_y_arr[width_arr_y_ind],1.e-9)<<std::endl;
  //  
  //  y_chi2NDF  = (fitfailed_y) ? -1. : fgaus_y->GetChisquare()/((float)fgaus_y->GetNDF());
  //  y_chi2prob = (fitfailed_y) ? -1. : fgaus_y->GetProb();
  //  
  //  if(   (f_y_mean<0.)||
  //	  (f_y_emean<0.)||
  //	  (f_y_sig<0.)||
  //	  (f_y_esig<0.)||
  //	  (y_chi2NDF<0.)||
  //	  (y_chi2prob<0.) ) std::cout<<"WARNING mu and/or sigma/mu are negative because FIT FAILED FOR |y| BIN!!!!"<<std::endl;
  //  
  //  
  //  
  //  std::cout<<"y_Chi2NDF ="   << y_chi2NDF<<std::endl;    
  //  std::cout<<"y_Chi2Prob ="   << y_chi2prob<<std::endl;          
  //  
  //  std::cout<<"setting bin content for fit hists ip="<< ip<<" "<<std::endl;    
  //  hMean_fit_y->SetBinContent (ip+1, f_y_mean);    
  //  hMean_fit_y->SetBinError   (ip+1, f_y_emean);    
  //  
  //  hSigma_fit_y->SetBinContent (ip+1, f_y_sig);    
  //  hSigma_fit_y->SetBinError   (ip+1, f_y_esig);    
  //  
  //  hSigmaMu_fit_y->SetBinContent (ip+1, f_y_sigmu);    
  //  hSigmaMu_fit_y->SetBinError   (ip+1, f_y_esigmu);    
  //  
  //  hChi2NDF_y->SetBinContent (ip+1, y_chi2NDF );
  //  hChi2NDF_y->SetBinError (ip+1, 1.e-05 );
  //  
  //  hChi2Prob_y->SetBinContent (ip+1, y_chi2prob );      
  //  
  //}// end loop over ptbins to yield mu, sigma/mu values
  
  //eta bins
  //TF1 *fgaus=NULL; 
  TF1 *fgaus[nbins_pt_debug]={};
  //double fgausMAX[nbins_pt_debug]={0.};
  TH1F *hrsp[nbins_pt_debug]={};  //input
  TH1F *hMean=NULL, *hSigma=NULL, *hSigmaMu=NULL; 
  TH1F *hMean_fit=NULL, *hSigma_fit=NULL, *hSigmaMu_fit=NULL;
  TH1F *hChi2NDF=NULL, *hChi2Prob=NULL;

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
  
  // QA for gauss core fits
  hChi2NDF = new TH1F("hChi2NDF_y", Form( "Fit Chi2NDF %s", algname.c_str()),    
		      nbins_pt_debug, ptbins_debug); 
  
  hChi2Prob = new TH1F("hChi2Prob_y", Form( "Fit Chi2Prob %s", algname.c_str()),    
		       nbins_pt_debug, ptbins_debug); 

  // pt bin hJER-fit loop    
  if(debugMode)std::cout<<"nbins_pt="<<nbins_pt_debug<<std::endl;    
  int width_arr_ind=0;

  for(int ip=0; ip<nbins_pt_debug; ip++){    
    

    int ptbin_ip=(int)ptbins_debug[ip], ptbin_ip1=(int)ptbins_debug[ip+1];      
    
    
    
    std::string inputHistName="";
    inputHistName="hJER_"+doJetID+"wJetID_"+bintype+"bin"+std::to_string(absetabin)+"_"+std::to_string(ptbin_ip)+"ptbin"+std::to_string(ptbin_ip1);      
    
    // open the input hist    
    hrsp[ip] = (TH1F*)finPP->Get( inputHistName.c_str() );              
    if(!hrsp[ip]){ 
      std::cout<<"no input hist named " <<  inputHistName<< ", continue..."<<std::endl;
      continue;}          
    else if(hrsp[ip]->GetEntries()<50){
      std::cout<<"inputHistName="<<inputHistName<<std::endl;
      std::cout<<"#entries= " <<  hrsp[ip]->GetEntries()<< ", too few, continue..."<<std::endl;
      continue;	      }
    
    //set name for later
    hrsp[ip]->SetName((
		       "hJER_1wJetID_"+bintype+"_"+
		       absetabins_str[absetabin]+"_to_"+absetabins_str[absetabin+1]+
		       "_genpt_"+std::to_string(ptbin_ip) + "_to_" + std::to_string(ptbin_ip1) 
		       ).c_str());
    
    //FYI debug
    if(debugMode||true)std::cout<<std::endl;  
    if(debugMode||true)std::cout<<"________________________________"<<std::endl;  
    if(debugMode||true)std::cout<<"pt range for bin: "<<ptbin_ip<<" - "<<ptbin_ip1<< " GeV "<<std::endl;    
    if(debugMode||true)std::cout<<"hist title = "<<hrsp[ip]->GetTitle()<<std::endl;
    if(debugMode||true)hrsp[ip]->Print("base");          
    if(debugMode||true)std::cout<<std::endl;  
    
    //rebin if desired
    if(rebinJER&&rebinFactor>1.){
      std::cout<<"WARNING: rebinning \n";
      hrsp[ip]= (TH1F*)hrsp[ip]->TH1::Rebin( rebinFactor, (inputHistName+"_rebin"+std::to_string(rebinFactor)).c_str());	
      //assert(false);
    }
    
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
    
    if(ip<fit_ptlo_bin )continue;
    if(ip>=fit_pthi_bin )continue;
    
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
    std::cout<<"h_integral ="   << hrsp_ip_rebin->Integral()<<std::endl;    
    
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
    std::string fgaus_name=   "fJER_1wJetID_"+bintype+"_"+
      absetabins_str[absetabin]+"_to_"+absetabins_str[absetabin+1]+
      "_genpt_"+std::to_string(ptbin_ip) + "_to_" + std::to_string(ptbin_ip1) ;

    if(useconstxrange){
      fgaus[ip] = new TF1(fgaus_name.c_str(),"gaus", fgaus_xlo,fgaus_xhi);  //what's been typicall used.      
    }
    else {
      fgaus[ip] = new TF1(fgaus_name.c_str(),"gaus", 1.-fgaus_width_arr[width_arr_ind],1.+fgaus_width_arr[width_arr_ind]);  //what's been typicall used.      
      width_arr_ind++;
    }
    //fgaus->SetParLimits(0,0.98,1.02);    //normalization
    //fgaus->SetParLimits(1,0.95,1.05);   // mean
    //fgaus->SetParLimits(2,0.0,0.5);     // width
    
    bool fitfailed=true;
    fitfailed = ((bool)(hrsp[ip]->Fit(fgaus[ip],"R")));                      
    if(fitfailed)
      std::cout<<std::endl<<"!!!WARNING FIT FAILED!!!"<<std::endl<<std::endl;
    else 
      std::cout<<std::endl<<"fit success"<<std::endl<<std::endl;

    std::cout<<"now writing histname "<<hrsp[ip]->GetName()<<std::endl;    
    fout_ptbin_JER->cd();
    fgaus[ip]->Write();
    hrsp[ip]->Write();
    rootfout->cd();

    
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
    
    f_norm  = (fitfailed)  ? -1. : fgaus[ip]->GetParameter(0); // z = if(condition) then(?) <do this> else(:) <do this>   
    f_enorm = (fitfailed)  ? -1. : fgaus[ip]->GetParError(0)*errfact ; 
    f_mean  = (fitfailed)  ? -1. : fgaus[ip]->GetParameter(1); // z = if(condition) then(?) <do this> else(:) <do this>   
    f_emean = (fitfailed)  ? -1. : fgaus[ip]->GetParError(1)*errfact ; 
    f_sig   = (fitfailed)  ? -1. : fgaus[ip]->GetParameter(2);
    f_esig  = (fitfailed)  ? -1. : fgaus[ip]->GetParError(2)*errfact ;      
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
    if(!useconstxrange)std::cout<<"TF1->Integral="<<fgaus[ip]->Integral(1.-fgaus_width_arr[width_arr_ind],1.+fgaus_width_arr[width_arr_ind],1.e-9)<<std::endl;
    
    chi2NDF = (fitfailed)  ? -1. : fgaus[ip]->GetChisquare()/((float)fgaus[ip]->GetNDF());
    chi2prob = (fitfailed) ? -1. : fgaus[ip]->GetProb();
    
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
    
    //fgausMAX[ip]=fgaus[ip]->GetMaximum(0.,2.);
    
    
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
    if(debugMode)hsigmaFitStatus=hSigma->Fit(SigmaFit_h,"QMER");	    
    else hsigmaFitStatus        =hSigma->Fit(SigmaFit_h,"QMER");	
    if(debugMode)hsigmaFitStatus=hSigma->Fit(SigmaFit_h,"QER");	    
    else hsigmaFitStatus        =hSigma->Fit(SigmaFit_h,"QER");	
  }
  else if(fitsigmu){
    if(debugMode)hsigmaFitStatus=hSigmaMu->Fit(SigmaFit_h,"QMER");	    
    else         hsigmaFitStatus=hSigmaMu->Fit(SigmaFit_h,"QMER");	
    if(debugMode)hsigmaFitStatus=hSigmaMu->Fit(SigmaFit_h,"QER");	    
    else         hsigmaFitStatus=hSigmaMu->Fit(SigmaFit_h,"QER");	
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
  
  


  std::string polNstr2;
  if(     absetabin==0)polNstr2="[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)+[5]*pow(log(x),5)";
  else if(absetabin==1)polNstr2="[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)+[5]*pow(log(x),5)";
  else if(absetabin==2)polNstr2="[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)+[5]*pow(log(x),5)";
  else if(absetabin==3)polNstr2="[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)+[5]*pow(log(x),5)+[6]*pow(log(x),6)";
  TF1 *MuFit_f= new TF1("MuFit_f",polNstr2.c_str(), fit_ptlo, fit_pthi);
  hMean_fit->Fit(MuFit_f,"MER");
  TFitResultPtr hMean_fit_rp=(TFitResultPtr)hMean_fit->Fit(MuFit_f,"ERS");
  
  std::string mufit_uncstring_dterms=TF1_mufituncstr_CovMat(hMean_fit_rp, MuFit_f->GetNpar(),true);//diagonal covariance terms only
  TF1 *MuFit_f_unc_dterms= new TF1("MuFit_f_unc_dterms",mufit_uncstring_dterms.c_str(), fit_ptlo, fit_pthi);
  std::cout<<MuFit_f_unc_dterms->GetNpar()<<"\n";
  
  std::string mufit_uncstring=TF1_mufituncstr_CovMat(hMean_fit_rp, MuFit_f->GetNpar(),false);//all terms
  TF1 *MuFit_f_unc= new TF1("MuFit_f_unc",mufit_uncstring.c_str(), fit_ptlo, fit_pthi);
  std::cout<<MuFit_f_unc->GetNpar()<<"\n";
  //assert(false);
  
  TF1* MuFit_f_sysup=new TF1("MuFit_f_sysup", (polNstr2+"+"+mufit_uncstring).c_str(),fit_ptlo,fit_pthi);//full unc calculation/propagation
  for(int i=0; i<MuFit_f->GetNpar();i++)
    MuFit_f_sysup->SetParameter(i,MuFit_f->GetParameter(i));
  MuFit_f_sysup->SetLineColor(kMagenta-8);
  MuFit_f_sysup->SetLineStyle(1);
  MuFit_f_sysup->SetLineWidth(1);

  TF1* MuFit_f_sysdown=new TF1("MuFit_f_sysdown", (polNstr2+"-"+mufit_uncstring).c_str(),fit_ptlo,fit_pthi);
  for(int i=0; i<MuFit_f->GetNpar();i++)
    MuFit_f_sysdown->SetParameter(i,MuFit_f->GetParameter(i));
  MuFit_f_sysdown->SetLineColor(kMagenta-8);
  MuFit_f_sysdown->SetLineStyle(1);
  MuFit_f_sysdown->SetLineWidth(1);
  
  TF1* MuFit_f_relsys=new TF1("MuFit_f_relsys", (mufit_uncstring+"/("+polNstr2+")").c_str(),fit_ptlo,fit_pthi);//full unc calculation/propagation
  for(int i=0; i<MuFit_f->GetNpar();i++)
    MuFit_f_relsys->SetParameter(i,MuFit_f->GetParameter(i));
  MuFit_f_relsys->SetLineColor(kMagenta-8);
  MuFit_f_relsys->SetLineStyle(1);
  MuFit_f_relsys->SetLineWidth(1);

  TF1* MuFit_f_2xsysup=new TF1("MuFit_f_2xsysup", (polNstr2+"+2.*"+mufit_uncstring).c_str(),fit_ptlo,fit_pthi);//full unc calculation/propagation
  for(int i=0; i<MuFit_f->GetNpar();i++)
    MuFit_f_2xsysup->SetParameter(i,MuFit_f->GetParameter(i));
  MuFit_f_2xsysup->SetLineColor(kMagenta-6);
  MuFit_f_2xsysup->SetLineStyle(3);
  MuFit_f_2xsysup->SetLineWidth(1);

  TF1* MuFit_f_2xsysdown=new TF1("MuFit_f_2xsysdown", (polNstr2+"-2.*"+mufit_uncstring).c_str(),fit_ptlo,fit_pthi);
  for(int i=0; i<MuFit_f->GetNpar();i++)
    MuFit_f_2xsysdown->SetParameter(i,MuFit_f->GetParameter(i));
  MuFit_f_2xsysdown->SetLineColor(kMagenta-6);
  MuFit_f_2xsysdown->SetLineStyle(2);
  MuFit_f_2xsysdown->SetLineWidth(1);

  TF1* MuFit_f_3xsysup=new TF1("MuFit_f_3xsysup", (polNstr2+"+3.*"+mufit_uncstring).c_str(),fit_ptlo,fit_pthi);//full unc calculation/propagation
  for(int i=0; i<MuFit_f->GetNpar();i++)
    MuFit_f_3xsysup->SetParameter(i,MuFit_f->GetParameter(i));
  MuFit_f_3xsysup->SetLineColor(kMagenta-4);
  MuFit_f_3xsysup->SetLineStyle(3);
  MuFit_f_3xsysup->SetLineWidth(1);

  TF1* MuFit_f_3xsysdown=new TF1("MuFit_f_3xsysdown", (polNstr2+"-3.*"+mufit_uncstring).c_str(),fit_ptlo,fit_pthi);
  for(int i=0; i<MuFit_f->GetNpar();i++)
    MuFit_f_3xsysdown->SetParameter(i,MuFit_f->GetParameter(i));
  MuFit_f_3xsysdown->SetLineColor(kMagenta-4);
  MuFit_f_3xsysdown->SetLineStyle(3);
  MuFit_f_3xsysdown->SetLineWidth(1);

  TF1* MuFit_f_10xsysup=new TF1("MuFit_f_10xsysup", (polNstr2+"+10.*"+mufit_uncstring).c_str(),fit_ptlo,fit_pthi);//full unc calculation/propagation
  for(int i=0; i<MuFit_f->GetNpar();i++)
    MuFit_f_10xsysup->SetParameter(i,MuFit_f->GetParameter(i));
  MuFit_f_10xsysup->SetLineColor(kMagenta-4);
  MuFit_f_10xsysup->SetLineStyle(3);
  MuFit_f_10xsysup->SetLineWidth(1);

  TF1* MuFit_f_10xsysdown=new TF1("MuFit_f_10xsysdown", (polNstr2+"-10.*"+mufit_uncstring).c_str(),fit_ptlo,fit_pthi);
  for(int i=0; i<MuFit_f->GetNpar();i++)
    MuFit_f_10xsysdown->SetParameter(i,MuFit_f->GetParameter(i));
  MuFit_f_10xsysdown->SetLineColor(kMagenta-4);
  MuFit_f_10xsysdown->SetLineStyle(3);
  MuFit_f_10xsysdown->SetLineWidth(1);

  TF1* MuFit_f_dsysup=new TF1("MuFit_f_dsysup", (polNstr2+"+"+mufit_uncstring_dterms).c_str(),fit_ptlo,fit_pthi);//typical unc. propagation w/o covariances accounted for (d stands for diagonal entries of covmatrix only)
  for(int i=0; i<MuFit_f->GetNpar();i++)
    MuFit_f_dsysup->SetParameter(i,MuFit_f->GetParameter(i));
  MuFit_f_dsysup->SetLineColor(kMagenta-3);
  MuFit_f_dsysup->SetLineStyle(9);
  MuFit_f_dsysup->SetLineWidth(1);


  TF1* MuFit_f_dsysdown=new TF1("MuFit_f_dsysdown", (polNstr2+"-"+mufit_uncstring_dterms).c_str(),fit_ptlo,fit_pthi);
  for(int i=0; i<MuFit_f->GetNpar();i++)
    MuFit_f_dsysdown->SetParameter(i,MuFit_f->GetParameter(i));
  MuFit_f_dsysdown->SetLineColor(kMagenta-3);
  MuFit_f_dsysdown->SetLineStyle(9);
  MuFit_f_dsysdown->SetLineWidth(1);

  TF1* MuFit_f_altsysup=new TF1("MuFit_f_altsysup", (polNstr2).c_str(),fit_ptlo,fit_pthi); //just shift fit params up/down, overestimate/simple band.
  for(int i=0; i<MuFit_f->GetNpar();i++)    
    MuFit_f_altsysup->SetParameter(i,MuFit_f->GetParameter(i) + MuFit_f->GetParError(i) );
  MuFit_f_altsysup->SetLineColor(kMagenta-9);
  MuFit_f_altsysup->SetLineStyle(10);
  MuFit_f_altsysup->SetLineWidth(1);

  
  TF1* MuFit_f_altsysdown=new TF1("MuFit_f_altsysdown", (polNstr2).c_str(),fit_ptlo,fit_pthi);
  for(int i=0; i<MuFit_f->GetNpar();i++)
    MuFit_f_altsysdown->SetParameter(i,MuFit_f->GetParameter(i) - MuFit_f->GetParError(i));
  MuFit_f_altsysdown->SetLineColor(kMagenta-9);
  MuFit_f_altsysdown->SetLineStyle(10);
  MuFit_f_altsysdown->SetLineWidth(1);


  
  // fit, mu v. gen pt, sigma/mu v. genpt hists, where vals are from gauss core fits
  //bool fsigmaFitStatus= true;
  std::string sigmafit_string="[0]+[1]/(pow(x,[2])+[3]*x)";
  TF1 *SigmaFit_f=new TF1("SigmaFit_f",sigmafit_string.c_str(),     // this fit func from raghav/HIN/2.76 TeV/5 TeV
			  fit_ptlo, fit_pthi );       // I want this fit to end at 967 GeV (for entire range)   
  TFitResultPtr hSigma_fit_rp;

  std::cout<<std::endl<<"fitting JER width v. pt using HIN/2.76 TeV fits... \n \n";
  if(fitsigma){
    if(debugMode)hSigma_fit_rp = (TFitResultPtr)hSigma_fit->Fit(SigmaFit_f,"MEVR");	
    else hSigma_fit_rp         = (TFitResultPtr)hSigma_fit->Fit(SigmaFit_f,"MER");	
    if(debugMode)hSigma_fit_rp = (TFitResultPtr)hSigma_fit->Fit(SigmaFit_f,"SMEVR");	
    else hSigma_fit_rp         = (TFitResultPtr)hSigma_fit->Fit(SigmaFit_f,"SMER");	
  }
  else if(fitsigmu){
    if(debugMode)hSigma_fit_rp = (TFitResultPtr)hSigmaMu_fit->Fit(SigmaFit_f,"MEVR");	
    else hSigma_fit_rp         = (TFitResultPtr)hSigmaMu_fit->Fit(SigmaFit_f,"MER");	
    if(debugMode)hSigma_fit_rp = (TFitResultPtr)hSigmaMu_fit->Fit(SigmaFit_f,"SMEVR");	
    else hSigma_fit_rp         = (TFitResultPtr)hSigmaMu_fit->Fit(SigmaFit_f,"SMER");	
  }
  else assert(false);
  std::string sigmafit_uncstring=TF1_sigmafituncstr_CovMat(hSigma_fit_rp,SigmaFit_f->GetNpar());

  //assert(false);
  
  bool fsigmaFit13tevStatus=true;
  TF1 *SigmaFit_f_13tev=new TF1("SigmaFit_f_13tev","[0]+([1]/(pow(x,[2])+[3]*x))+[4]*x",     // this fit func 13 TeV
				fit_ptlo,
				fit_pthi );       // I want this fit to end at 967 GeV (for entire range)    
  
  
  std::cout<<std::endl<<"fitting JER width v. pt using SMP/13 TeV fits... init. sigmaFitStatus="<<fsigmaFit13tevStatus<<std::endl<<std::endl; 
  if(fitsigma){
    if(debugMode)fsigmaFit13tevStatus = hSigma_fit->Fit(SigmaFit_f_13tev,"QMER");	
    else fsigmaFit13tevStatus         = hSigma_fit->Fit(SigmaFit_f_13tev,"QMER");	
    if(debugMode)fsigmaFit13tevStatus = hSigma_fit->Fit(SigmaFit_f_13tev,"QER");	
    else fsigmaFit13tevStatus         = hSigma_fit->Fit(SigmaFit_f_13tev,"QER");	
  }
  else if (fitsigmu){
    if(debugMode)fsigmaFit13tevStatus = hSigmaMu_fit->Fit(SigmaFit_f_13tev,"QMER");	
    else fsigmaFit13tevStatus         = hSigmaMu_fit->Fit(SigmaFit_f_13tev,"QMER");	
    if(debugMode)fsigmaFit13tevStatus = hSigmaMu_fit->Fit(SigmaFit_f_13tev,"QEVR");	
    else fsigmaFit13tevStatus         = hSigmaMu_fit->Fit(SigmaFit_f_13tev,"QER");	
  }
  else assert(false);


  bool fsigmaFit8tevStatus=true;
  TF1 *SigmaFit_f_8tev=new TF1("SigmaFit_f_8tev","sqrt( pow( [0]/x, 2) + pow([1],2)/x + pow([0],2)) ",     // this fit func 8 TeV
			       fit_ptlo,
			       fit_pthi );       // I want this fit to end at 967 GeV (for entire range)    
  
  std::cout<<std::endl<<"fitting JER width v. pt using SMP/13 TeV fits... init. sigmaFitStatus="<<fsigmaFit8tevStatus<<std::endl<<std::endl; 
  if(fitsigma){
    if(debugMode)fsigmaFit8tevStatus = hSigma_fit->Fit(SigmaFit_f_8tev,"QMER");	
    else fsigmaFit8tevStatus         = hSigma_fit->Fit(SigmaFit_f_8tev,"QMER");	
    if(debugMode)fsigmaFit8tevStatus = hSigma_fit->Fit(SigmaFit_f_8tev,"QER");	
    else fsigmaFit8tevStatus         = hSigma_fit->Fit(SigmaFit_f_8tev,"QER");	
  }
  else if (fitsigmu) {
    if(debugMode)fsigmaFit8tevStatus = hSigmaMu_fit->Fit(SigmaFit_f_8tev,"QMER");	
    else fsigmaFit8tevStatus         = hSigmaMu_fit->Fit(SigmaFit_f_8tev,"QMER");	
    if(debugMode)fsigmaFit8tevStatus = hSigmaMu_fit->Fit(SigmaFit_f_8tev,"QER");	
    else fsigmaFit8tevStatus         = hSigmaMu_fit->Fit(SigmaFit_f_8tev,"QER");	
  }
  else assert(false);
  
  TF1 *SigmaFit_sysup_f=new TF1("SigmaFit_sysup_f", (sigmafit_string+"+"+sigmafit_uncstring).c_str(),
				fit_ptlo,
				fit_pthi );       
  for(int i=0; i<SigmaFit_sysup_f->GetNpar();i++)
    SigmaFit_sysup_f->SetParameter(i, SigmaFit_f->GetParameter(i));
  SigmaFit_sysup_f->SetLineWidth(1);
  SigmaFit_sysup_f->SetLineColor(kMagenta-8);
  SigmaFit_sysup_f->SetLineStyle(2);
  
  TF1 *SigmaFit_sysdown_f=new TF1("SigmaFit_sysdown_f",(sigmafit_string+"-"+sigmafit_uncstring).c_str(),
				  fit_ptlo,
				  fit_pthi );       
  for(int i=0; i<SigmaFit_sysdown_f->GetNpar();i++)
    SigmaFit_sysdown_f->SetParameter(i, SigmaFit_f->GetParameter(i));
  SigmaFit_sysdown_f->SetLineWidth(1);
  SigmaFit_sysdown_f->SetLineColor(kMagenta-8);
  SigmaFit_sysdown_f->SetLineStyle(2);


  TF1 *SigmaFit_10xsysup_f=new TF1("SigmaFit_10xsysup_f", (sigmafit_string+"+10.*"+sigmafit_uncstring).c_str(),
				fit_ptlo,
				fit_pthi );       
  for(int i=0; i<SigmaFit_10xsysup_f->GetNpar();i++)
    SigmaFit_10xsysup_f->SetParameter(i, SigmaFit_f->GetParameter(i));
  SigmaFit_10xsysup_f->SetLineWidth(1);
  SigmaFit_10xsysup_f->SetLineColor(kMagenta-4);
  SigmaFit_10xsysup_f->SetLineStyle(3);
  
  TF1 *SigmaFit_10xsysdown_f=new TF1("SigmaFit_10xsysdown_f",(sigmafit_string+"-10.*"+sigmafit_uncstring).c_str(),
				  fit_ptlo,
				  fit_pthi );       
  for(int i=0; i<SigmaFit_10xsysdown_f->GetNpar();i++)
    SigmaFit_10xsysdown_f->SetParameter(i, SigmaFit_f->GetParameter(i));
  SigmaFit_10xsysdown_f->SetLineWidth(1);
  SigmaFit_10xsysdown_f->SetLineColor(kMagenta-4);
  SigmaFit_10xsysdown_f->SetLineStyle(3);

  TF1 *SigmaFit_relsys_f=new TF1("SigmaFit_relsys_f",(sigmafit_uncstring+"/("+sigmafit_string+")").c_str(),
				  fit_ptlo,
				  fit_pthi );       
  for(int i=0; i<SigmaFit_relsys_f->GetNpar();i++)
    SigmaFit_relsys_f->SetParameter(i, SigmaFit_f->GetParameter(i));
  SigmaFit_relsys_f->SetLineWidth(1);
  SigmaFit_relsys_f->SetLineStyle(3);

  
  TH1F *SigmaFit_f_hist= ((TH1F*) ((TH1F*)SigmaFit_f->GetHistogram())->Clone("SigmaFit_f_hist"));
  TH1F *SigmaFit_sysup_f_hist= ((TH1F*) ((TH1F*)SigmaFit_sysup_f->GetHistogram())->Clone("SigmaFit_sysup_f_hist"));
  TH1F *SigmaFit_sysdown_f_hist= ((TH1F*) ((TH1F*)SigmaFit_sysdown_f->GetHistogram())->Clone("SigmaFit_sysdown_f_hist"));
  


  TF1 *SigmaFit_altsysup_f=new TF1("SigmaFit_altsysup_f","[0]+[1]/(pow(x,[2])+[3]*x)", 
				fit_ptlo,
				fit_pthi );       
  SigmaFit_altsysup_f->SetParameter(0, SigmaFit_f->GetParameter(0)+SigmaFit_f->GetParError(0));
  SigmaFit_altsysup_f->SetParameter(1, SigmaFit_f->GetParameter(1)+SigmaFit_f->GetParError(1));
  SigmaFit_altsysup_f->SetParameter(2, SigmaFit_f->GetParameter(2)-SigmaFit_f->GetParError(2));//par 2 negative, carefi;
  SigmaFit_altsysup_f->SetParameter(3, SigmaFit_f->GetParameter(3)-SigmaFit_f->GetParError(3));
  
  TF1 *SigmaFit_altsysdown_f=new TF1("SigmaFit_altsysdown_f","[0]+[1]/(pow(x,[2])+[3]*x)", 
				  fit_ptlo,
				  fit_pthi );       
  SigmaFit_altsysdown_f->SetParameter(0, SigmaFit_f->GetParameter(0)-SigmaFit_f->GetParError(0));
  SigmaFit_altsysdown_f->SetParameter(1, SigmaFit_f->GetParameter(1)-SigmaFit_f->GetParError(1));
  SigmaFit_altsysdown_f->SetParameter(2, SigmaFit_f->GetParameter(2)+SigmaFit_f->GetParError(2));//par 2 negative, carefi;
  SigmaFit_altsysdown_f->SetParameter(3, SigmaFit_f->GetParameter(3)+SigmaFit_f->GetParError(3));
  
  
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


  //assert(false);

  TFile *findscb=NULL;

  TH1F *hMu_dscb=NULL;
  TH1F *hMu_dscb_clone=NULL;
  TH1F *hMu_dscb_exclude_clone=NULL;
  TF1 *MuFit_dscb=NULL;    
  TF1 *MuFit_dscb_unc_dterms=NULL;  TF1 *MuFit_dscb_unc	    =NULL;
  TF1* MuFit_dscb_sysup	    =NULL;  TF1* MuFit_dscb_sysdown   =NULL;
  TF1* MuFit_dscb_dsysup    =NULL;  TF1* MuFit_dscb_dsysdown  =NULL;
  TF1* MuFit_dscb_altsysup  =NULL;  TF1* MuFit_dscb_altsysdown=NULL;



  TH1F *hSigma_dscb=NULL;
  TH1F *hSigma_dscb_clone=NULL;
  TH1F *hSigma_dscb_exclude_clone=NULL;
  TF1 *SigmaFit_dscb=NULL;


  TH1F *hSigmaMu_dscb=NULL;
  TH1F *hSigmaMu_dscb_clone=NULL;
  TH1F *hSigmaMu_dscb_exclude_clone=NULL;
  TF1 *SigmaMuFit_dscb=NULL;

  
  if(fitJohnsDSCB){

    std::string findscb_etabin_str;
    if(     absetabin==0)findscb_etabin_str="00eta05";
    else if(absetabin==1)findscb_etabin_str="05eta10";
    else if(absetabin==2)findscb_etabin_str="10eta15";
    else if(absetabin==3)findscb_etabin_str="15eta20";
    else assert(false);
    std::string findscb_name=dscb_dir+dscb_finstr_start+ findscb_etabin_str + dscb_finstr_end;
    findscb=TFile::Open(findscb_name.c_str(),"READ");


    float dscbfit_ptlo;
    if(absetabin==0){
      dscbfit_ptlo=28.;
      //dscbfit_ptlo=32.;
    }
    else
      dscbfit_ptlo=43.;

    //assert(false);
    hMu_dscb=(TH1F*)findscb->Get(hjohn_dscb_mu.c_str());    
    hMu_dscb_clone=(TH1F*)hMu_dscb->Clone(((std::string)hMu_dscb->GetName()+"_clone").c_str());
    hMu_dscb_clone->Reset("MICES");
    //assert(false);

    //assert(false);


    hMu_dscb_exclude_clone=(TH1F*)hMu_dscb->Clone(((std::string)hMu_dscb->GetName()+"_exclude_clone").c_str());
    hMu_dscb_exclude_clone->Reset("MICES");

    bool dontFitThis_Mu=false;
    for(int i=1; i<=hMu_dscb->GetNbinsX(); i++){
      if(     absetabin==0){
	if(i==999){
	//if( i==12   ) { 	  
	//if( i==9 || i==32 ) { 	  
	  std::cout<<"hMu_dscb->GetBinLowEdge("<<i<<")="<<hMu_dscb->GetBinLowEdge(i)<<std::endl;
	  dontFitThis_Mu=true;
	} }
      else if(absetabin==1){
	if(i==999){
	  //if(i==8 || i==9){
	  std::cout<<"hMu_dscb->GetBinLowEdge("<<i<<")="<<hMu_dscb->GetBinLowEdge(i)<<std::endl;
	  dontFitThis_Mu=true;
	} }
      else if(absetabin==2){
	if(i==999){
	  //if(i==11 || i==12){
	  std::cout<<"hMu_dscb->GetBinLowEdge("<<i<<")="<<hMu_dscb->GetBinLowEdge(i)<<std::endl;
	  dontFitThis_Mu=true;
	} }
      else if(absetabin==3){
	if(i==999){
	  //if(i==10){
	  std::cout<<"hMu_dscb->GetBinLowEdge("<<i<<")="<<hMu_dscb->GetBinLowEdge(i)<<std::endl;
	  dontFitThis_Mu=true;
	} }
      
      if(dontFitThis_Mu){	
	hMu_dscb_exclude_clone->SetBinContent(i, hMu_dscb->GetBinContent(i));
	//hMu_dscb_exclude_clone->SetBinError(i, 1.e-10*hMu_dscb->GetBinContent(i));
	hMu_dscb_exclude_clone->SetBinError(i, hMu_dscb->GetBinError(i));
      }
      else{
	hMu_dscb_clone->SetBinContent(i, hMu_dscb->GetBinContent(i));
	//hMu_dscb_clone->SetBinError(i, 1.e-10*hMu_dscb->GetBinContent(i));
	hMu_dscb_clone->SetBinError(i, hMu_dscb->GetBinError(i));
      }
      dontFitThis_Mu=false;
    }

    //assert(false);
    hSigma_dscb=(TH1F*)findscb->Get((hjohn_dscb_sigma).c_str());
    hSigma_dscb_clone=(TH1F*)hSigma_dscb->Clone(((std::string)hSigma_dscb->GetName()+"_clone").c_str());
    hSigma_dscb_clone->Reset("MICES");

    hSigma_dscb_exclude_clone=(TH1F*)hSigma_dscb->Clone(((std::string)hSigma_dscb->GetName()+"_exclude_clone").c_str());
    hSigma_dscb_exclude_clone->Reset("MICES");
    
    bool dontFitThis_Sigma=false;
    for(int i=1; i<=hSigma_dscb->GetNbinsX(); i++){
      if(     absetabin==0){
	if(i==999){
	  //if(i == 12  ){
	  //if(i == 9 || i==32 || i==13){
	  std::cout<<"hSigma_dscb->GetBinLowEdge("<<i<<")="<<hSigma_dscb->GetBinLowEdge(i)<<std::endl;
	  dontFitThis_Sigma=true;
	} }
      else if(absetabin==1){
	if(i==999){
	  //if(i==8 || i==9){
	  std::cout<<"hSigma_dscb->GetBinLowEdge("<<i<<")="<<hSigma_dscb->GetBinLowEdge(i)<<std::endl;
	  dontFitThis_Sigma=true;
	} }
      else if(absetabin==2){
	if(i==999){
	//if( i==13){
	  std::cout<<"hSigma_dscb->GetBinLowEdge("<<i<<")="<<hSigma_dscb->GetBinLowEdge(i)<<std::endl;
	  dontFitThis_Sigma=true;
	} }
      else if(absetabin==3){
	if(i==999){
	  //if(i==10){
	  std::cout<<"hSigma_dscb->GetBinLowEdge("<<i<<")="<<hSigma_dscb->GetBinLowEdge(i)<<std::endl;
	  dontFitThis_Sigma=true;
	} }
      
      if(dontFitThis_Sigma){	
	hSigma_dscb_exclude_clone->SetBinContent(i, hSigma_dscb->GetBinContent(i));
	//hSigma_dscb_exclude_clone->SetBinError(i, 1.e-10*hSigma_dscb->GetBinContent(i));
	hSigma_dscb_exclude_clone->SetBinError(i, hSigma_dscb->GetBinError(i));
      }
      else{
	hSigma_dscb_clone->SetBinContent(i, hSigma_dscb->GetBinContent(i));
	//hSigma_dscb_clone->SetBinError(i, 1.e-10*hSigma_dscb->GetBinContent(i));
	hSigma_dscb_clone->SetBinError(i, hSigma_dscb->GetBinError(i));
      }
      dontFitThis_Sigma=false;
    }
    //assert(false);




    std::string hSigmaMu_dscb_name=(std::string)hSigma_dscb->GetName();
    hSigmaMu_dscb=(TH1F*)hSigma_dscb->Clone( hSigmaMu_dscb_name.c_str());
    for(int i=1; i<=hSigmaMu_dscb->GetNbinsX();i++)
      if(hSigmaMu_dscb->GetBinContent(i)>0. && hMu_dscb->GetBinContent(i)>0.)
	hSigmaMu_dscb->SetBinContent( i , hSigmaMu_dscb->GetBinContent(i)/hMu_dscb->GetBinContent(i));
    hSigmaMu_dscb_clone=(TH1F*)hSigmaMu_dscb->Clone(((std::string)hSigmaMu_dscb->GetName()+"_clone").c_str());
    hSigmaMu_dscb_clone->Reset("MICES");
    
    hSigmaMu_dscb_exclude_clone=(TH1F*)hSigmaMu_dscb->Clone(((std::string)hSigmaMu_dscb->GetName()+"_exclude_clone").c_str());
    hSigmaMu_dscb_exclude_clone->Reset("MICES");
    
    bool dontFitThis_SigmaMu=false;
    for(int i=1; i<=hSigmaMu_dscb->GetNbinsX(); i++){
      if(     absetabin==0){
	if(i==999){
	  //if( i==7){
	  //if(i ==9 || i==32 || i==13){
	  std::cout<<"hSigmaMu_dscb->GetBinLowEdge("<<i<<")="<<hSigmaMu_dscb->GetBinLowEdge(i)<<std::endl;
	  dontFitThis_SigmaMu=true;}
      }
      else if(absetabin==1){
	if(i==999){
	  //if(i==8 || i==9){
	  std::cout<<"hSigmaMu_dscb->GetBinLowEdge("<<i<<")="<<hSigmaMu_dscb->GetBinLowEdge(i)<<std::endl;
	  dontFitThis_SigmaMu=true;}
      } 
      else if(absetabin==2){
	if(i==999){
	  //if(i==30 || i==11){
	  std::cout<<"hSigmaMu_dscb->GetBinLowEdge("<<i<<")="<<hSigmaMu_dscb->GetBinLowEdge(i)<<std::endl;
	  dontFitThis_SigmaMu=true;}
      } 
      else if(absetabin==3){
	if(i==999){
	  //if(i==10){
	  std::cout<<"hSigmaMu_dscb->GetBinLowEdge("<<i<<")="<<hSigmaMu_dscb->GetBinLowEdge(i)<<std::endl;
	  dontFitThis_SigmaMu=true;}
      } 
      
      if(dontFitThis_SigmaMu){	
	hSigmaMu_dscb_exclude_clone->SetBinContent(i, hSigmaMu_dscb->GetBinContent(i));
	hSigmaMu_dscb_exclude_clone->SetBinError(i, 1.e-10*hSigmaMu_dscb->GetBinContent(i));
      }
      else{
	hSigmaMu_dscb_clone->SetBinContent(i, hSigmaMu_dscb->GetBinContent(i));
	hSigmaMu_dscb_clone->SetBinError(i, 1.e-10*hSigmaMu_dscb->GetBinContent(i));
      }
      dontFitThis_SigmaMu=false;
    }
    //assert(false);


    //assert(false);




    std::string polNstr;
    //if(     absetabin==0)polNstr="pol4";
    //else if(absetabin==1)polNstr="pol4";
    //else if(absetabin==2)polNstr="pol4";
    //else if(absetabin==3)polNstr="pol3";
    if(     absetabin==0)polNstr="[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)+[5]*pow(log(x),5)";//+[8]*pow(log(x),8)";
    else if(absetabin==1)polNstr="[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)+[5]*pow(log(x),5)";//+[8]*pow(log(x),8)";
    else if(absetabin==2)polNstr="[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)+[5]*pow(log(x),5)";//+[8]*pow(log(x),8)";
    else if(absetabin==3)polNstr="[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)+[5]*pow(log(x),5)+[6]*pow(log(x),6)";//6'th order seems good here..

    bool hmudscbFitStatus= true;    
    std::string dscb_fitOpts;
    //if(absetabin==0)
    //dscb_fitOpts="MER";//this followed by "SER" will give the fit of john's params with uncertainties, after i widened the pT limits on JER measurements (currently only works for 00y05, see "TESTFIT2" in OUTDESC in run.sh script
    //else
    dscb_fitOpts="WMER";//this followed by "SWER" will give the fit of johns' params before i widened pT limits on JER measurements
    
    MuFit_dscb=new TF1("MuFit_dscb", polNstr.c_str() ,
		       dscbfit_ptlo, //fit_ptlo,			  
		       fit_pthi );
    
    std::cout<<std::endl<<"fitting JER DSCB mu v. pt... init. sigmaFitStatus="<<hmudscbFitStatus<<std::endl<<std::endl; 
    hMu_dscb_clone->Fit(MuFit_dscb, dscb_fitOpts.c_str());
    TFitResultPtr hMu_dscb_rp;
    
    //if(absetabin==0)   {
    //hMu_dscb_rp=(TFitResultPtr)hMu_dscb_clone->Fit(MuFit_dscb, "SER"); 
    //}
    //else            {
    hMu_dscb_rp=(TFitResultPtr)hMu_dscb_clone->Fit(MuFit_dscb, "SWER");
    //}
    
    std::string dscb_uncstring_dterms=TF1_mufituncstr_CovMat(hMu_dscb_rp, MuFit_dscb->GetNpar(),true);//diagonal covariance terms only
    MuFit_dscb_unc_dterms= new TF1("MuFit_dscb_f_unc_dterms",dscb_uncstring_dterms.c_str(), dscbfit_ptlo, fit_pthi);
    std::cout<<MuFit_dscb_unc_dterms->GetNpar()<<"\n";
    
    std::string dscb_uncstring=TF1_mufituncstr_CovMat(hMu_dscb_rp, MuFit_dscb->GetNpar(),false);//all terms
    MuFit_dscb_unc= new TF1("MuFit_dscb_unc",dscb_uncstring.c_str(), dscbfit_ptlo, fit_pthi);
    std::cout<<MuFit_dscb_unc->GetNpar()<<"\n";
    //assert(false);
    
    MuFit_dscb_sysup=new TF1("MuFit_dscb_sysup", (polNstr+"+"+dscb_uncstring).c_str(),fit_ptlo,fit_pthi);//full unc calculation/propagation
    for(int i=0; i<MuFit_dscb->GetNpar();i++)
      MuFit_dscb_sysup->SetParameter(i,MuFit_dscb->GetParameter(i));
    MuFit_dscb_sysup->SetLineColor(kBlue-8);
    MuFit_dscb_sysup->SetLineStyle(1);
    MuFit_dscb_sysup->SetLineWidth(1);
    
    MuFit_dscb_sysdown=new TF1("MuFit_dscb_sysdown", (polNstr+"-"+dscb_uncstring).c_str(),fit_ptlo,fit_pthi);
    for(int i=0; i<MuFit_dscb->GetNpar();i++)
      MuFit_dscb_sysdown->SetParameter(i,MuFit_dscb->GetParameter(i));
    MuFit_dscb_sysdown->SetLineColor(kBlue-8);
    MuFit_dscb_sysdown->SetLineStyle(1);
    MuFit_dscb_sysdown->SetLineWidth(1);
    
    MuFit_dscb_dsysup=new TF1("MuFit_dscb_dsysup", (polNstr+"+"+dscb_uncstring_dterms).c_str(),fit_ptlo,fit_pthi);//typical unc. propagation w/o covariances accounted for (d stands for diagonal entries of covmatrix only)
    for(int i=0; i<MuFit_dscb->GetNpar();i++)
      MuFit_dscb_dsysup->SetParameter(i,MuFit_dscb->GetParameter(i));
    MuFit_dscb_dsysup->SetLineColor(kBlue-3);
    MuFit_dscb_dsysup->SetLineStyle(9);
    MuFit_dscb_dsysup->SetLineWidth(1);
    
    
    MuFit_dscb_dsysdown=new TF1("MuFit_dscb_dsysdown", (polNstr+"-"+dscb_uncstring_dterms).c_str(),fit_ptlo,fit_pthi);
    for(int i=0; i<MuFit_dscb->GetNpar();i++)
      MuFit_dscb_dsysdown->SetParameter(i,MuFit_dscb->GetParameter(i));
    MuFit_dscb_dsysdown->SetLineColor(kBlue-3);
    MuFit_dscb_dsysdown->SetLineStyle(9);
    MuFit_dscb_dsysdown->SetLineWidth(1);
    
    MuFit_dscb_altsysup=new TF1("MuFit_dscb_altsysup", (polNstr).c_str(),fit_ptlo,fit_pthi); //just shift fit params up/down, overestimate/simple band.
    for(int i=0; i<MuFit_dscb->GetNpar();i++)    
      MuFit_dscb_altsysup->SetParameter(i,MuFit_dscb->GetParameter(i) + MuFit_dscb->GetParError(i) );
    MuFit_dscb_altsysup->SetLineColor(kBlue-9);
    MuFit_dscb_altsysup->SetLineStyle(10);
    MuFit_dscb_altsysup->SetLineWidth(1);
    
    
    MuFit_dscb_altsysdown=new TF1("MuFit_dscb_altsysdown", (polNstr).c_str(),fit_ptlo,fit_pthi);
    for(int i=0; i<MuFit_dscb->GetNpar();i++)
      MuFit_dscb_altsysdown->SetParameter(i,MuFit_dscb->GetParameter(i) - MuFit_dscb->GetParError(i));
    MuFit_dscb_altsysdown->SetLineColor(kBlue-9);
    MuFit_dscb_altsysdown->SetLineStyle(10);
    MuFit_dscb_altsysdown->SetLineWidth(1);




    dscb_fitOpts="WME";//for fits other than the mean
    bool hsigmadscbFitStatus= true;    
    SigmaFit_dscb=new TF1("SigmaFit_dscb","[0]+[1]/(pow(x,[2])+[3]*x)", 
			  dscbfit_ptlo, //fit_ptlo,			  
			  fit_pthi ); 
    
    //hSigma_dscb->Fit(SigmaFit_dscb,"MER");
    //std::cout<<std::endl<<"fitting JER DSCB width v. pt...[WME opts] init. sigmaFitStatus="<<hsigmadscbFitStatus<<std::endl<<std::endl; 
    //hSigma_dscb_clone->Fit(SigmaFit_dscb,"WME");
    std::cout<<std::endl<<"fitting JER DSCB width v. pt...[WE opts] init. sigmaFitStatus="<<hsigmadscbFitStatus<<std::endl<<std::endl; 
    hSigma_dscb_clone->Fit(SigmaFit_dscb,dscb_fitOpts.c_str());//"WE");

    float dscb_fitParam_0    =SigmaFit_dscb->GetParameter(0);
    float dscb_fitParam_0_err=SigmaFit_dscb->GetParError(0);
    float dscb_fitParam_1    =SigmaFit_dscb->GetParameter(1);
    float dscb_fitParam_1_err=SigmaFit_dscb->GetParError(1);
    float dscb_fitParam_2    =SigmaFit_dscb->GetParameter(2);
    float dscb_fitParam_2_err=SigmaFit_dscb->GetParError(2); 
    float dscb_fitParam_3    =SigmaFit_dscb->GetParameter(3);
    float dscb_fitParam_3_err=SigmaFit_dscb->GetParError(3); 

    std::cout << "dscb_fitParam_0     = " << dscb_fitParam_0    << std::endl;
    std::cout << "dscb_fitParam_0_err = " << dscb_fitParam_0_err<< std::endl;
    std::cout << "%err dscb_fitParam0 = " << (dscb_fitParam_0_err/dscb_fitParam_0)*100.    << std::endl;
    std::cout << "dscb_fitParam_1     = " << dscb_fitParam_1    << std::endl;
    std::cout << "dscb_fitParam_1_err = " << dscb_fitParam_1_err<< std::endl;
    std::cout << "%err dscb_fitParam1 = " << (dscb_fitParam_1_err/dscb_fitParam_1)*100.    << std::endl;
    std::cout << "dscb_fitParam_2     = " << dscb_fitParam_2    << std::endl;
    std::cout << "dscb_fitParam_2_err =   " << dscb_fitParam_2_err<< std::endl;
    std::cout << "%err dscb_fitParam2 = " << (dscb_fitParam_2_err/dscb_fitParam_2)*100.    << std::endl;
    std::cout << "dscb_fitParam_3     = " << dscb_fitParam_3    << std::endl;
    std::cout << "dscb_fitParam_3_err =   " << dscb_fitParam_3_err<< std::endl;
    std::cout << "%err dscb_fitParam3 = " << (dscb_fitParam_3_err/dscb_fitParam_3)*100.    << std::endl;        




    bool hsigmudscbFitStatus= true;    
    SigmaMuFit_dscb=new TF1("SigmaMuFit_dscb","[0]+[1]/(pow(x,[2])+[3]*x)", 
			  dscbfit_ptlo, //fit_ptlo,			  
			  fit_pthi ); 
    
    //hSigma_dscb->Fit(SigmaMuFit_dscb,"MER");
    //std::cout<<std::endl<<"fitting JER DSCB width v. pt...[WME opts] init. sigmaFitStatus="<<hsigmudscbfitstatus<<std::endl<<std::endl; 
    //hSigma_dscb_clone->Fit(SigmaMuFit_dscb,"WME");
    std::cout<<std::endl<<"fitting JER DSCB width v. pt...[WE opts] init. sigmaFitStatus="<<hsigmudscbFitStatus<<std::endl<<std::endl; 
    hSigmaMu_dscb_clone->Fit(SigmaMuFit_dscb,dscb_fitOpts.c_str());//"WE");

    float dscb_fitParam2_0    =SigmaMuFit_dscb->GetParameter(0);
    float dscb_fitParam2_0_err=SigmaMuFit_dscb->GetParError(0);
    float dscb_fitParam2_1    =SigmaMuFit_dscb->GetParameter(1);
    float dscb_fitParam2_1_err=SigmaMuFit_dscb->GetParError(1);
    float dscb_fitParam2_2    =SigmaMuFit_dscb->GetParameter(2);
    float dscb_fitParam2_2_err=SigmaMuFit_dscb->GetParError(2); 
    float dscb_fitParam2_3    =SigmaMuFit_dscb->GetParameter(3);
    float dscb_fitParam2_3_err=SigmaMuFit_dscb->GetParError(3); 

    std::cout << "dscb_fitParam2_0     = " << dscb_fitParam2_0    << std::endl;
    std::cout << "dscb_fitParam2_0_err = " << dscb_fitParam2_0_err<< std::endl;
    std::cout << "%err dscb_fitParam20 = " << (dscb_fitParam2_0_err/dscb_fitParam2_0)*100.    << std::endl;
    std::cout << "dscb_fitParam2_1     = " << dscb_fitParam2_1    << std::endl;
    std::cout << "dscb_fitParam2_1_err = " << dscb_fitParam2_1_err<< std::endl;
    std::cout << "%err dscb_fitParam21 = " << (dscb_fitParam2_1_err/dscb_fitParam2_1)*100.    << std::endl;
    std::cout << "dscb_fitParam2_2     = " << dscb_fitParam2_2    << std::endl;
    std::cout << "dscb_fitParam2_2_err =   " << dscb_fitParam2_2_err<< std::endl;
    std::cout << "%err dscb_fitParam22 = " << (dscb_fitParam2_2_err/dscb_fitParam2_2)*100.    << std::endl;
    std::cout << "dscb_fitParam2_3     = " << dscb_fitParam2_3    << std::endl;
    std::cout << "dscb_fitParam2_3_err =   " << dscb_fitParam2_3_err<< std::endl;
    std::cout << "%err dscb_fitParam23 = " << (dscb_fitParam2_3_err/dscb_fitParam2_3)*100.    << std::endl;        

    //assert(false);

    
  }
  
  
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
  //float xmin=ptbins_debug[fit_ptlo_bin-1], xmax=ptbins_debug[fit_pthi_bin];
  float xmin=ptbins_debug[fit_ptlo_bin], xmax=ptbins_debug[fit_pthi_bin];
  if(debugMode){
    std::cout<<"xmin="<<xmin<<std::endl;
    std::cout<<"xmax="<<xmax<<std::endl;    }
  
  TLine* lineAtOne=new TLine(xmin,1.,ptbins_debug[fit_pthi_bin],1.);    
  lineAtOne->SetLineStyle(2);  
  lineAtOne->SetLineColor(kBlack);
  TLine* lineAtTen=new TLine(xmin,10.,ptbins_debug[fit_pthi_bin],10.);    
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

  MakeHistMu( (TH1F*)hMean_fit, xmin,xmax);   
  hMean_fit->SetMarkerStyle(kOpenStar);
  //hMean_fit->SetMarkerSize(0);
  hMean_fit->SetMarkerColor(kBlack);
  hMean_fit->SetLineColor(kMagenta);
  hMean_fit->SetLineWidth(1);
  //hMean_fit->DrawClone("HIST ][E SAME");
  hMean_fit->DrawClone("HIST ][E");

  MakeHistMu( (TH1F*)hMean,xmin,xmax);
  //hMean->DrawClone("HIST ][E");  

  MuFit_f->SetLineColor(kMagenta+1);
  MuFit_f->SetLineWidth(1);
  MuFit_f->Draw("SAME");

  if(drawJohnsDSCB){
    hMu_dscb_clone->SetMarkerStyle(kOpenSquare);
    hMu_dscb_clone->SetMarkerColor(kBlue-1);
    hMu_dscb_clone->DrawClone("HIST ][E SAME");
    hMu_dscb_exclude_clone->SetMarkerStyle(kFullCross);
    hMu_dscb_exclude_clone->SetMarkerColor(kRed);
    hMu_dscb_exclude_clone->SetLineColor(kRed);
    hMu_dscb_exclude_clone->DrawClone("HIST ][E SAME");

    MuFit_dscb->SetLineColor(kBlue);
    MuFit_dscb->SetLineWidth(1);
    MuFit_dscb->Draw("SAME");

    MuFit_dscb_sysup->Draw("SAME");
    MuFit_dscb_dsysup->Draw("SAME");
    MuFit_dscb_altsysup->Draw("SAME");
    MuFit_dscb_sysdown->Draw("SAME");
    MuFit_dscb_dsysdown->Draw("SAME");
    MuFit_dscb_altsysdown->Draw("SAME");

    hMu_dscb_clone->DrawClone("HIST ][E SAME");

  }
  else{
    MuFit_f_sysup->Draw("SAME");
    //MuFit_f_2xsysup->Draw("SAME");
    //MuFit_f_3xsysup->Draw("SAME");
    MuFit_f_10xsysup->Draw("SAME");
    //MuFit_f_dsysup->Draw("SAME");
    //MuFit_f_altsysup->Draw("SAME");
    MuFit_f_sysdown->Draw("SAME");
    //MuFit_f_2xsysdown->Draw("SAME");
    ///MuFit_f_3xsysdown->Draw("SAME");
    MuFit_f_10xsysdown->Draw("SAME");
    //MuFit_f_dsysdown->Draw("SAME");
    //MuFit_f_altsysdown->Draw("SAME");
    
    hMean_fit->DrawClone("HIST ][E SAME");  

  }
  
  TLegend* hMean_Leg=new TLegend( 0.60, 0.60 ,0.9 , 0.97, NULL,"brNDC");
  hMean_Leg->SetBorderSize(0);
  hMean_Leg->SetFillStyle(0);
  
  //hMean_Leg->AddEntry( hMean, "#mu from Hists", "lp");
  hMean_Leg->AddEntry( hMean_fit, "#mu from Gauss Fits", "lpe");
  hMean_Leg->AddEntry( MuFit_f, "Fit", "l");
  if(drawJohnsDSCB){
    hMean_Leg->AddEntry( hMu_dscb_clone, "#mu from DSCB Fits", "lpe");
    hMean_Leg->AddEntry( hMu_dscb_exclude_clone, "excluded #mu from DSCB Fits", "lpe");
    hMean_Leg->AddEntry( MuFit_dscb, "Fit of #mu from DSCB Fits", "l");
    
    //std::string MuFit_dscb_prob_str;
    //makeFitProbStr(MuFit_dscb, &MuFit_dscb_prob_str);
    //if(debugMode)std::cout<<"MuFit_dscb_prob_str="<<MuFit_dscb_prob_str<<std::endl;    
    //std::string MuFit_dscb_chi2_str;
    //makeFitChi2Str(MuFit_dscb, &MuFit_dscb_chi2_str);
    //if(debugMode)std::cout<<"MuFit_dscb_chi2_str="<<MuFit_dscb_chi2_str<<std::endl;    
    //std::string MuFit_dscb_ndf_str=std::to_string( MuFit_dscb->GetNDF() );
    //if(debugMode)std::cout<<"MuFit_dscb_ndf_str="<<MuFit_dscb_ndf_str<<std::endl;    
    
    //hMean_Leg->AddEntry( MuFit_dscb, "Fit of #mu from DSCB Fits", "l"); 
    //hMean_Leg->AddEntry( (TObject*)0, ("#chi^{2}/NDF="+MuFit_dscb_chi2_str+"/"+MuFit_dscb_ndf_str).c_str() , "");   
    //hMean_Leg->AddEntry( (TObject*)0, ("Prob="+MuFit_dscb_prob_str).c_str() , "");   
    
  }
  else{

    hMean_Leg->AddEntry( MuFit_f_sysup, "Fit Unc.", "l");
    //hMean_Leg->AddEntry( MuFit_f_2xsysup, "2#times(Fit Unc.)", "l");
    //hMean_Leg->AddEntry( MuFit_f_3xsysup, "3#times(Fit Unc.)", "l");
    hMean_Leg->AddEntry( MuFit_f_10xsysup, "10#times(Fit Unc.)", "l");
    //hMean_Leg->AddEntry( MuFit_f_dsysup, "Fit Unc. (no covariance)", "l");
    //hMean_Leg->AddEntry( MuFit_f_altsysup, "Param. Shift up/down by unc.", "l");

  }
    

  hMean_Leg->Draw();
  
  lineAtOne->Draw();
  
  float posx=0.40, posy=0.90, size=0.045;
  TLatex* jetdesc=new TLatex(posx, posy, fullJetType.c_str()); 
  jetdesc->SetTextSize(size);
  jetdesc->SetNDC();
  jetdesc->Draw();
  
  std::string absetarange_str=absetabins_str[absetabin];
  if     (bintype=="y")      absetarange_str+=" < #||{y^{RECO}} < ";
  else if(bintype=="eta")    absetarange_str+=" < #||{#eta^{RECO}} < ";
  else if(bintype=="geny")   absetarange_str+=" < #||{y^{GEN}} < ";
  else if(bintype=="geneta") absetarange_str+=" < #||{#eta^{GEN}} < ";
  absetarange_str+= absetabins_str[absetabin+1];
  TLatex* absetarange=new TLatex(posx, posy-0.05, 
				 absetarange_str.c_str());
  absetarange->SetTextSize(size);
  absetarange->SetNDC();
  absetarange->Draw();
  
  
  /////////////////draw pad 1
  
  
  //draw pad 2
  p2->cd();    
  
  MakeHistSigma( (TH1F*) hSigma_fit ,xmin,xmax); 
  //this is a hist created from the gauss fits in each bin
  hSigma_fit->SetMarkerStyle(kOpenStar);
  hSigma_fit->SetMarkerColor(kBlack);
  hSigma_fit->SetMarkerSize(0);
  hSigma_fit->SetLineColor(kMagenta);
  hSigma_fit->SetLineWidth(1);  
  hSigma_fit->DrawClone("HIST ][ E");
  
  
  if(drawJohnsDSCB){
    //MakeHistSigma( (TH1F*) hSigma_dscb_clone ,xmin,xmax);   
    hSigma_dscb_clone->SetMarkerStyle(kFullSquare);
    hSigma_dscb_clone->DrawClone("HIST ][ E SAME");
    hSigma_dscb_exclude_clone->SetMarkerStyle(kFullCross);
    hSigma_dscb_exclude_clone->SetMarkerColor(kRed);
    hSigma_dscb_exclude_clone->SetLineColor(kRed);
    hSigma_dscb_exclude_clone->DrawClone("HIST ][ E SAME");
  }
  else{
    MakeHistSigma( (TH1F*) hSigma ,xmin,xmax);   
    //hSigma->DrawClone("HIST ][ E SAME");
  }
  
      
  
  if(fitsigma){
    //this is a TF1 of sigma from hists
    SigmaFit_h->SetLineColor(kBlue-7);
    //SigmaFit_h->Draw("SAME");
    
    //this is a TF1 of sigma from fits of hists
    SigmaFit_f->SetLineColor(kMagenta-2);
    SigmaFit_f->Draw("SAME");
    SigmaFit_sysup_f->Draw("SAME");
    SigmaFit_sysdown_f->Draw("SAME");
    SigmaFit_10xsysup_f->Draw("SAME");
    SigmaFit_10xsysdown_f->Draw("SAME");
    
    SigmaFit_f_13tev->SetLineColor(kGreen-2);
    //SigmaFit_f_13tev->Draw("SAME");
    
    SigmaFit_f_8tev->SetLineColor(kTeal);
    //SigmaFit_f_8tev->Draw("SAME");
    
    if(drawJohnsDSCB){
      SigmaFit_dscb->SetLineColor(kBlue);
      SigmaFit_dscb->SetLineWidth(1);
      SigmaFit_dscb->Draw("SAME");
    }
  }
  
  TLegend* hSigma_Leg=new TLegend( 0.60, 0.50 , 0.90 , 0.97, NULL,"brNDC");
  hSigma_Leg->AddEntry( hSigma_fit, "#sigma from Gauss Fits", "lp");
  if(drawJohnsDSCB)  {
    hSigma_Leg->AddEntry( hSigma_dscb_clone, "#sigma from DSCB Fits", "lp");
    hSigma_Leg->AddEntry( hSigma_dscb_exclude_clone, "excluded #sigma from DSCB Fits", "lp");
  }
  //else hSigma_Leg->AddEntry( hSigma, "#sigma from Hists", "lp");
 
  
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
    
  
  if(fitsigma){
    
    //std::string SigmaFit_f_prob_str;
    //makeFitProbStr(SigmaFit_f, &SigmaFit_f_prob_str);
    //if(debugMode)std::cout<<"SigmaFit_f_prob_str="<<SigmaFit_f_prob_str<<std::endl;    
    //std::string SigmaFit_f_chi2_str;
    //makeFitChi2Str(SigmaFit_f, &SigmaFit_f_chi2_str);
    //if(debugMode)std::cout<<"SigmaFit_f_chi2_str="<<SigmaFit_f_chi2_str<<std::endl;    
    //std::string SigmaFit_f_ndf_str=std::to_string( SigmaFit_f->GetNDF() );
    //if(debugMode)std::cout<<"SigmaFit_f_ndf_str="<<SigmaFit_f_ndf_str<<std::endl;    
    
    
    hSigma_Leg->AddEntry( SigmaFit_f, "Fit", "l");
    hSigma_Leg->AddEntry( SigmaFit_sysup_f, "Fit Unc.", "l");
    hSigma_Leg->AddEntry( SigmaFit_10xsysup_f, "10#times(Fit Unc.)", "l");
    //hSigma_Leg->AddEntry( (TObject*)0, ("#chi^{2}/NDF="+SigmaFit_f_chi2_str+"/"+SigmaFit_f_ndf_str).c_str() , "");   
    //hSigma_Leg->AddEntry( (TObject*)0, ("Prob="+SigmaFit_f_prob_str).c_str() , "");   
    
    //std::string SigmaFit_dscb_prob_str;
    //makeFitProbStr(SigmaFit_dscb, &SigmaFit_dscb_prob_str);
    //if(debugMode)std::cout<<"SigmaFit_dscb_prob_str="<<SigmaFit_dscb_prob_str<<std::endl;    
    //std::string SigmaFit_dscb_chi2_str;
    //makeFitChi2Str(SigmaFit_dscb, &SigmaFit_dscb_chi2_str);
    //if(debugMode)std::cout<<"SigmaFit_dscb_chi2_str="<<SigmaFit_dscb_chi2_str<<std::endl;    
    //std::string SigmaFit_dscb_ndf_str=std::to_string( SigmaFit_dscb->GetNDF() );
    //if(debugMode)std::cout<<"SigmaFit_dscb_ndf_str="<<SigmaFit_dscb_ndf_str<<std::endl;    
    //    
    if(drawJohnsDSCB)hSigma_Leg->AddEntry( SigmaFit_dscb, "Fit of #sigma from DSCB Fits", "l"); 
    //hSigma_Leg->AddEntry( (TObject*)0, ("#chi^{2}/NDF="+SigmaFit_dscb_chi2_str+"/"+SigmaFit_dscb_ndf_str).c_str() , "");   
    //hSigma_Leg->AddEntry( (TObject*)0, ("Prob="+SigmaFit_dscb_prob_str).c_str() , "");   
 }
  
  
  hSigma_Leg->SetBorderSize(0);    
  hSigma_Leg->SetFillStyle(0);        
  hSigma_Leg->Draw();
  /////////////////draw pad2    
  
  //draw pad 3
  p3->cd();    
  
  MakeHistSigmaMu( (TH1F*) hSigmaMu ,xmin,xmax); 
  //MakeHistSigma( (TH1F*) hSigma_fit ,xmin,xmax); 
  //hSigmaMu->DrawClone("HIST ][E");
  
  //this is a hist created from the gauss fits in each bin
  MakeHistSigmaMu( (TH1F*) hSigmaMu_fit, xmin, xmax);
  hSigmaMu_fit->SetMarkerStyle(kOpenStar);
  hSigmaMu_fit->SetMarkerColor(kBlack);
  hSigmaMu_fit->SetLineColor(kMagenta);
  //hSigmaMu_fit->DrawClone("HIST ][E SAME");
  hSigmaMu_fit->DrawClone("HIST ][E");


  if(drawJohnsDSCB){
    hSigmaMu_dscb_clone->SetMarkerStyle(kFullSquare);
    hSigmaMu_dscb_clone->DrawClone("HIST ][E SAME");
    hSigmaMu_dscb_exclude_clone->SetMarkerStyle(kFullCross);
    hSigmaMu_dscb_exclude_clone->SetMarkerColor(kRed);
    hSigmaMu_dscb_exclude_clone->SetLineColor(kRed);
    hSigmaMu_dscb_exclude_clone->DrawClone("HIST ][E SAME");
  }
  
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

    if(drawJohnsDSCB)
      SigmaMuFit_dscb->Draw("SAME");
    
  }
  
  TLegend* hSigmaMu_Leg=new TLegend( 0.60, 0.50 , 0.90 , 0.97, NULL,"brNDC");
  //hSigmaMu_Leg->AddEntry( hSigma, "#sigma/#mu from Hists", "lp");
  
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
  if(drawJohnsDSCB)
    hSigmaMu_Leg->AddEntry( hSigma_dscb, "#sigma/#mu from DSCB Fits", "lp");
  
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
    std::cout<<"j="<<j<<std::endl;
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
  //assert(false);
  
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
  //TF1* patclone4Draw[nbins_pt_debug]={};
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


      //
      ////if(drawAltDSCBFits){
      ////if(drawAltDSCBFits && ptbins_debug[j]>50. && ptbins_debug[j]<58.){
      //if(drawAltDSCBFits){
      //	std::string patDSCB_name="dJER_1wJetID_geny_0.0_to_0.5_genpt_" + std::to_string((int)ptbins_debug[j])+ "_to_"+std::to_string((int)ptbins_debug[j+1])+"/f";
      //	//std::string johnDSCB_name="dJER_1wJetID_geny_0.0_to_0.5_genpt_"+ std::to_string((int)ptbins_debug[j])+ "_to_"+std::to_string((int)ptbins_debug[j+1])+"/johnf";
      //	
      //	std::cout<<"patDSCB_name="<<patDSCB_name<<"\n";
      //	TF1* patDSCB=(TF1*)findscb->Get(patDSCB_name.c_str());
      //	//if(!patDSCB)std::cout<<"this fit doesnt exist (i don't think) \n";
      //	//else        std::cout<<"this fit   DOES exist (i think) \n";
      //
      //	//std::cout<<"johnDSCB_name="<<johnDSCB_name<<"\n";
      //	//TF1* johnDSCB=(TF1*)findscb->Get(johnDSCB_name.c_str());
      //	//if(!johnDSCB)std::cout<<"this fit doesnt exist (i don't think) \n";
      //	//else         std::cout<<"this fit   DOES exist (i think) \n";
      //	//assert(false);
      //	
      //	if(patDSCB) {
      //
      //
      //	  double patnorm=fgausMAX[j]/
      //	    patDSCB->GetMaximum(0., 2.);
      //	  std::cout<<"patnorm="<<patnorm<<"\n";
      //	  
      //	  //neat trick; TF1 as functor, using lambda functions https://root.cern.ch/doc/master/classTF1.html#F5
      //	  patclone4Draw[j] = new TF1(("patclone_drawonly_"+std::to_string((int)ptbins_debug[j])+"_to_"+std::to_string((int)ptbins_debug[j+1])).c_str(),
      //				     [&](double*x, double *p){ return p[0]*patDSCB->Eval(x[0]); },
      //				     0., 2., 1);
      //	  patclone4Draw[j]->SetParameter(0., patnorm);
      //	  
      //	  std::cout<<"test Eval patclone @ 1 = "<<patclone4Draw[j]->Eval(1.);
      //	  
      //	  //patDSCB->SetLineColor(kGreen);
      //	  //patDSCB->Draw("SAME");
      //	  patclone4Draw[j]->SetLineColor(kGreen);
      //	  patclone4Draw[j]->Draw("SAME");
      //	  //assert(false);
      //	}
      //	//if(johnDSCB){
      //	//  johnDSCB->SetLineColor(kBlue);
      //	//  johnDSCB->Draw("SAME");
      //	//}
      //	
      //}//end draw altDSCB fits
      
      drawcount++;	
    }//end loop over pads of one of the 3x3 canvases
    
    
    pdfoutCanv_3x3_arr[i]->Print(thePDFFileName.c_str());
    //assert(false);
    rootfout->cd();
    //assert(false);
    std::string fit_pt_start_str= std::to_string(((int) ptbins_debug[fit_ptbin_start]));
    std::string fit_pt_end_str= std::to_string(((int) ptbins_debug[fit_ptbin_end]));
    if(debugMode){
      std::cout<<"fit_ptbin_start ="<< fit_ptbin_start     <<std::endl;
      std::cout<<"fit_pt_start_str ="<< fit_pt_start_str <<std::endl;
      std::cout<<"fit_ptbin_end   ="<< fit_ptbin_end       <<std::endl;
      std::cout<<"fit_pt_end_str   ="<< fit_pt_end_str   <<std::endl;
    }
    pdfoutCanv_3x3_arr[i]->SetTitle(("3x3 JER fit canv for bin"+ fit_pt_start_str + "to bin"+fit_pt_end_str ).c_str());
    //assert(false);
    //pdfoutCanv_3x3_arr[i]->Write(("JERfitCanv_3x3_bins_"       + fit_pt_start_str + "_to_"  +fit_pt_end_str ).c_str());
    //assert(false);
  }//end loop over 3x3canvases
  
  std::cout<<"closing output pdf file "<<thePDFFileName<<std::endl<<std::endl;                
  pdfoutCanv->Print( close_thePDFFileName.c_str() );        
  
  rootfout->cd();
  
  
  
  
  //FIT + HIST WRITING GOES HERE; NO CANVASES
  
  //_______________________________________________    
  //mu + sigma/mu plots, vals from histogram contents
  hMean->Write();                   hSigma->Write();               
  //fits of vals from hist contents
  SigmaFit_h->Write();    SigmaFit_sysup_h->Write();    SigmaFit_sysdown_h->Write();
  //histogram versions of fits for viewing on my local machine, which CRASHES b.c. of TF1 for some reason
  SigmaFit_h_hist->Write();    SigmaFit_sysup_h_hist->Write();    SigmaFit_sysdown_h_hist->Write();
  
  
  //_______________________________________________
  //mu + sigma/mu plots + chi2/ndf, vals from gauss core fits 
  hMean_fit->Write();    hSigma_fit->Write();    hChi2NDF->Write();                         
  //fits w/ sys err envelope of the gauss core fit extracted values
  MuFit_f->Write();  MuFit_f_relsys->Write();
  MuFit_f_sysup->Write(); MuFit_f_sysdown->Write();
  MuFit_f_2xsysup->Write(); MuFit_f_2xsysdown->Write();
  MuFit_f_3xsysup->Write(); MuFit_f_3xsysdown->Write();
  MuFit_f_10xsysup->Write(); MuFit_f_10xsysdown->Write();
  MuFit_f_dsysup->Write(); MuFit_f_dsysdown->Write();
  MuFit_f_altsysup->Write(); MuFit_f_altsysdown->Write();
  MuFit_f_unc->Write(); MuFit_f_unc_dterms->Write();
  SigmaFit_f->Write();  SigmaFit_relsys_f->Write();    
  SigmaFit_sysup_f->Write();    SigmaFit_sysdown_f->Write();    
  SigmaFit_10xsysup_f->Write();    SigmaFit_10xsysdown_f->Write();    
  //histogram versions of fits for viewing on my local machine, which CRASHES b.c. of TF1 for some reason
  SigmaFit_f_hist->Write();    SigmaFit_sysup_f_hist->Write();    SigmaFit_sysdown_f_hist->Write();
  
  if(fitJohnsDSCB){
    
    
    //mu + sigma vals from histogram contents
    hSigma_dscb_clone->Write("hSigma_dscbfit");    
    hMu_dscb_clone->Write("hMu_dscbfit");

    //fits of extracted dscb parameters
    MuFit_dscb->Write(); 

    MuFit_dscb_sysup->Write(); MuFit_dscb_sysdown->Write();
    MuFit_dscb_dsysup->Write(); MuFit_dscb_dsysdown->Write();
    MuFit_dscb_altsysup->Write(); MuFit_dscb_altsysdown->Write();
    MuFit_dscb_unc->Write(); MuFit_dscb_unc_dterms->Write();

    SigmaFit_dscb->Write();
    SigmaMuFit_dscb->Write();
  }
  
  
  
  
  
  
  
  
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
  if(argc==8)
    rStatus=printPlots_ppMC_JERS_etabin((const std::string)argv[1], (const std::string)argv[2] , (const bool)std::atoi(argv[3]), (const bool)std::atoi(argv[4]), (const int) std::atoi(argv[5]) , (const std::string)argv[6] , (const std::string)argv[7]);
  if(argc==7)
    rStatus=printPlots_ppMC_JERS_etabin((const std::string)argv[1], (const std::string)argv[2] , (const bool)std::atoi(argv[3]), (const bool)std::atoi(argv[4]), (const int) std::atoi(argv[5]) , (const std::string)argv[6] );
  if(argc==6)
    rStatus=printPlots_ppMC_JERS_etabin((const std::string)argv[1], (const std::string)argv[2] , (const bool)std::atoi(argv[3]), (const bool)std::atoi(argv[4]), (const int) std::atoi(argv[5]) );
  if(argc==5)
    rStatus=printPlots_ppMC_JERS_etabin((const std::string)argv[1], (const std::string)argv[2] , (const bool)std::atoi(argv[3]), (const bool)std::atoi(argv[4]), -1);
  return rStatus;
}


























/*
  
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
  } // end  MCEff 
*/
