#include "unfoldSpectra.h"

//other settings
const bool drawPDFs=true; 
const bool debugMode=false, debugWrite=false;
const bool drawRespMatrix=true;
const bool dokIterQA=true;
const bool doJetID=true;

// CODE --------------------------------------------------
int bayesUnfoldNLOMCSpectra(  std::string inFile_MC_dir , std::string inFile_MC_name, //input 
			      std::string baseName ,  //output 
			      bool useNPCorrSpectra=true   , bool useSimpBins=false, const int kIterInput=4 ){
  
  // BINNING -----------  
  if(!useSimpBins)std::cout<<"using analysis pt bins"<<std::endl;
  else std::cout<<"using simple pt bins"<<std::endl<<std::endl;
  
  double* boundaries_pt_gen=setBinning(  useSimpBins , "gen"  );
  int     nbins_pt_gen=setNBins(  useSimpBins , "gen"  );  
  double* boundaries_pt_reco      = setBinning( useSimpBins , "reco" );
  int     nbins_pt_reco =   setNBins(  useSimpBins , "reco"  );
  double* boundaries_pt_gen_mat   = setBinning( useSimpBins , "gen"  );
  int     nbins_pt_gen_mat =   setNBins(  useSimpBins , "gen"  );
  double* boundaries_pt_reco_mat  = setBinning( useSimpBins , "reco" );
  int     nbins_pt_reco_mat =   setNBins(  useSimpBins , "reco"  );
  
  // ERROR/WEIGHTS/STATS HANDLING ------------------
  RooUnfold::ErrorTreatment errorTreatment;
  if(!doToyErrs) errorTreatment = RooUnfold::kCovariance;
  else errorTreatment = RooUnfold::kCovToy; 
  if(debugMode)std::cout<<"doToyErrs="<<doToyErrs<<std::endl; 
  
  if(debugMode)std::cout<<"TH2 GetDefaultSumw2="<<TH2::GetDefaultSumw2()<<std::endl;
  if(debugMode)std::cout<<"setting TH2DefSumw2 to true..."<<std::endl; 
  TH2::SetDefaultSumw2(true);
  
  if(debugMode)std::cout<<"TH1 GetDefaultSumw2="<<TH1::GetDefaultSumw2()<<std::endl;
  if(debugMode)std::cout<<"setting TH1DefSumw2 to true..."<<std::endl; 
  TH1::SetDefaultSumw2(true);
  
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  
  
  

  // STRINGS -----------
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl; 
  // inFile_MC_dir  =SCRATCH_BASE+inFile_MC_dir;
  
  //std::size_t radPos=inFile_MC_dir.find("_ak")+3;  
  //const std::string radius="R"+inFile_MC_dir.substr( radPos,1 );
  const std::string radiusint="4";
  const std::string radius="R"+radiusint;
  const std::string RandEtaRange="_"+radius+"_20_eta_20";
  const std::string RandEtaRange_plotTitle=" "+radius+" 20eta20";// R4 20eta20
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;
  if(debugMode)std::cout<<"RandEtaRange string is = "<<RandEtaRange<<std::endl;  
  
  //std::size_t jetTypePos=radPos+1, jetsPos=inFile_MC_dir.find("Jets");
  const std::string jetType="PFJets";//inFile_MC_dir.substr( jetTypePos,(jetsPos-jetTypePos) );//PFJets
  const std::string fullJetType="ak"+radiusint+jetType;//inFile_MC_dir.substr( radPos,1 )+jetType;//"ak4PFJets"
  if(debugMode)std::cout<<"jetType string is = "<<jetType<<std::endl;
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  
  std::string NLOMCtitle_str;
  if(     inFile_MC_name.find("CT10")!=std::string::npos)   NLOMCtitle_str="CT10 NLO";
  else if(inFile_MC_name.find("CT14")!=std::string::npos)   NLOMCtitle_str="CT14 NLO";
  else if(inFile_MC_name.find("HERAPDF")!=std::string::npos)NLOMCtitle_str="HERAPDF15 NLO";
  else if(inFile_MC_name.find("MMHT")!=std::string::npos)   NLOMCtitle_str="MMHT2014 NLO";
  else if(inFile_MC_name.find("NNPDF")!=std::string::npos)  NLOMCtitle_str="NNPDF30 NNLO";
  if(useNPCorrSpectra)NLOMCtitle_str+="+NP";

  // OUTPUT FILE, NAME(S) -----------
  std::string outFileName=unfoldMCSpectra_outdir+fullJetType;
  if(doJetID)outFileName+="_wjtID";
  if(!useSimpBins)outFileName+="_anabins";
  else outFileName+="_simpbins";
  if(useNPCorrSpectra) outFileName+="_NPNLO";
  else outFileName+="_NLO";
  outFileName+=baseName;
  
  std::string outRespMatPdfFile =  outFileName+"_respMat.pdf";
  std::string outBayesPdfFile =  outFileName+".pdf";
  std::string outRootFile     =  outFileName+".root";  
  
  checkNRenameFiles ( (const std::string) outFileName, 
		      &outRespMatPdfFile, 
		      &outBayesPdfFile, 
		      &outRootFile);
  
  
  
  if(debugMode)std::cout<<"opening output file: "<<outRootFile<<std::endl;  
  TFile* fout = new TFile(outRootFile.c_str(),"RECREATE");   
  if(debugMode)fout->cd();  
  
  
  // ppMC input histos -------------------------
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from MC file"<<std::endl; 
  std::cout<<"input MC dir : "<< (inFile_MC_dir)  <<std::endl; 
  std::cout<<"MC file name : "<< (inFile_MC_name)<<std::endl;   std::cout<<std::endl<<std::endl;  
  TFile *fin_ppMC = TFile::Open( (inFile_MC_dir+inFile_MC_name).c_str());
  
  if(debugWrite)
    fout->cd();
  
  // ---------- open MC "response" matrix
  std::string TH2_title="response";//"_ynew_th2";//"hpp_mcclosure_matrix";//+RandEtaRange;
  if(useNPCorrSpectra)TH2_title+="_NPynew_th2";
  else TH2_title+="_ynew_th2";

  //get the response matrix made by readforests
  TH2D* hmat = (TH2D*)fin_ppMC->Get( TH2_title.c_str() );
  hmat->Scale(1./NLOMCscaling);//05/09/18
  if(debugWrite)hmat->Write(TH2_title.c_str());
  if(debugMode)hmat->Print("base");
  
  // rebinned matrix ---------------
  TH2_title+="_clone";
  TH2D* hmat_rebin = (TH2D*)hmat->Clone( (TH2_title).c_str() );
  if(debugWrite)hmat_rebin->Write(TH2_title.c_str());
  if(debugMode)hmat_rebin->Print("base"); 
  
  TH2_title+="_rebins";
  hmat_rebin=(TH2D*) reBinTH2(hmat_rebin, (TH2_title).c_str(), 
  			      (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
			      (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );  
  if(debugWrite)hmat_rebin->Write(TH2_title.c_str());
  if(debugMode)hmat_rebin->Print("base"); 
  
  if(clearOverUnderflows){
    TH2_title+="_noOverUnderFlows";
    TH2clearOverUnderflows((TH2D*)hmat_rebin);
    if(debugWrite)hmat_rebin->Write(TH2_title.c_str());
    if(debugMode)hmat_rebin->Print("base");  }
  
  
  // error and %error for response matrix ---------------
  std::string errTH2_title="hmat_errors_rebins";
  TH2D* hmat_errors=makeRespMatrixErrors( (TH2D*) hmat,
  					  (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
  					  (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  ); 
  
  //TH2D* hmat_errors=makeNLORespMatrixErrors( (TH2D*) hmat,
  //					     (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
  //					     (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
  
  //TH2D* hmat_errors= (TH2D*)hmat_rebin->Clone(errTH2_title.c_str());
  //hmat_errors->Reset("M ICES");
  //for(int i=1; i<=(hmat_errors->GetNbinsX());i++)
  //for(int j=1; j<=(hmat_errors->GetNbinsY());j++)
  //hmat_errors->SetBinContent(i,j,hmat_rebin->GetBinError(i,j));
  
  if(debugWrite)hmat_errors->Write(errTH2_title.c_str());
  if(debugMode)hmat_errors->Print("base");
  
  if(clearOverUnderflows){
    errTH2_title+="_noOverUnderFlows";
    TH2clearOverUnderflows((TH2D*)hmat_errors);
    if(debugWrite)hmat_errors->Write(errTH2_title.c_str());
    if(debugMode)hmat_errors->Print("base");  }
  
  //// give response matrix the correct errors
  setRespMatrixErrs( (TH2D*)hmat_rebin, (TH2D*) hmat_errors , (bool)zeroBins);  
  TH2_title+="_wseterrs";
  if(debugWrite)hmat_rebin->Write(TH2_title.c_str());
  if(debugMode)hmat_rebin->Print("base");
  
  // ---------- reco, measured same-side spectra used to create response matrix, and for "sameside" unfolding test
  std::string histTitle2="smeared_rnd";//"_NPynew";//"hpp_mcclosure_reco";//    hpp_reco_wJetID_R4_20_eta_20
  if(useNPCorrSpectra) histTitle2+="_NPynew";
  else histTitle2+="_ynew";
  //  if(doJetID)histTitle2+="_wJetID";
  //  histTitle2+=RandEtaRange;
  
  TH1D*  hrec_sameside = (TH1D*)fin_ppMC->Get( histTitle2.c_str() );   
  hrec_sameside->Scale(1./NLOMCscaling);
  hrec_sameside->Scale(4.);
  multiplyBinWidth(hrec_sameside);
  if(debugWrite)hrec_sameside->Write(histTitle2.c_str());
  if(debugMode)hrec_sameside->Print("base");
  
  histTitle2+="_clone";
  TH1D *hrec_sameside_rebin = (TH1D*)hrec_sameside->Clone( (histTitle2).c_str() );
  std::cout<<"rebinning hrec_sameside..."<<std::endl;
  histTitle2+="_rebins";
  hrec_sameside_rebin = (TH1D*)hrec_sameside_rebin->Rebin( nbins_pt_reco, (histTitle2).c_str() , boundaries_pt_reco);
  if(debugWrite)hrec_sameside_rebin->Write( histTitle2.c_str() );   
  if(debugMode)hrec_sameside_rebin->Print("base");  

  
  if(clearOverUnderflows){
    histTitle2+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hrec_sameside_rebin);
    if(debugWrite)hrec_sameside_rebin->Write( histTitle2.c_str() );
    if(debugMode)hrec_sameside_rebin->Print("base");    }
  
  //cosmetics
  hrec_sameside_rebin->SetMarkerStyle(kOpenSquare);
  hrec_sameside_rebin->SetMarkerColor(kBlue-3);     
  hrec_sameside_rebin->SetLineColor(kBlue-3);     
  hrec_sameside_rebin->SetMarkerSize(1.02);     
  
  
  
  
  // ---------- gen, MC truth spectra
  std::string genHistTitle="theory_rnd";//"hpp_mcclosure_gen"
  if(useNPCorrSpectra) genHistTitle+="_NPynew";
  else genHistTitle+="_ynew";
 
  TH1D* hgen = (TH1D*)fin_ppMC->Get( genHistTitle.c_str() );
  hgen->Scale(1./NLOMCscaling);//05/09/18
  hgen->Scale(4.);
  multiplyBinWidth(hgen);
  if(debugWrite)hgen->Write(genHistTitle.c_str());
  if(debugMode)hgen->Print("base");    
  
  genHistTitle+="_clone";
  TH1D* hgen_rebin = (TH1D*)hgen->Clone( (genHistTitle).c_str() );
  genHistTitle+="_rebins";
  hgen_rebin = (TH1D*)hgen_rebin->Rebin(nbins_pt_gen, (genHistTitle).c_str() , boundaries_pt_gen);
  if(debugWrite)hgen_rebin->Write(genHistTitle.c_str());
  if(debugMode)hgen_rebin->Print("base"); 
  
  if(clearOverUnderflows){
    genHistTitle+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hgen_rebin);
    if(debugWrite)hgen_rebin->Write(genHistTitle.c_str());
    if(debugMode)hgen_rebin->Print("base");    }

  //cosmetics
  hgen_rebin->SetMarkerStyle(kOpenStar);
  hgen_rebin->SetMarkerColor(kMagenta);
  hgen_rebin->SetLineColor(kMagenta);
  hgen_rebin->SetMarkerSize(1.02);     	



  
  // reco/gen response hist

  std::string histTitle3="hpp_rec_sameside_response_rebin"+RandEtaRange;
  TH1D* hrec_sameside_resp_rebin=NULL;  
  std::string histTitle3gen="hpp_gen_response_rebin"+RandEtaRange;
  TH1D* hgen_resp_rebin=NULL;

  if(fillRespHists) {
    //if(useTH2ProjRespHist){
    //  
    //  //reco response histo
    //  if(debugMode)std::cout<<"boundaries_pt_gen[0]  = "<< ((int)boundaries_pt_gen[0]) << std::endl;
    //  if(debugMode)std::cout<<"boudnaries_pt_gen["<<nbins_pt_gen<<"] = "<< (int)boundaries_pt_gen[nbins_pt_gen] <<std::endl;
    //  if(debugMode)std::cout<<"bin 1 pt low = "<<
    //		     (int)hmat->GetYaxis()->GetBinLowEdge((int)(boundaries_pt_gen[0]+1.))<<std::endl;
    //  
    //  if(debugMode)std::cout<<"bin "<<nbins_pt_gen<<" pt hi = "<<
    //		     hmat->GetYaxis()->GetBinLowEdge(((int)boundaries_pt_gen[nbins_pt_gen]))
    //		     + hmat->GetYaxis()->GetBinWidth(((int)boundaries_pt_gen[nbins_pt_gen]))<<std::endl;    
    //
    //  //profile error types; default(sum/sqrt(N)), s, i, *g*
    //  hrec_sameside_resp_rebin = (TH1D*)hmat->TH2::ProjectionX( 
    //							    (histTitle3+"_profile").c_str() , 
    //							    (int)(boundaries_pt_gen[0]+1.), 
    //							    ((int)boundaries_pt_gen[nbins_pt_gen]), 	"g");
    //  hrec_sameside_resp_rebin=(TH1D*)hrec_sameside_resp_rebin->Rebin( nbins_pt_reco, (histTitle3+"_profile_rebin").c_str() , boundaries_pt_reco);
    //  hrec_sameside_resp_rebin->Write( (histTitle3+"_profile_rebin").c_str());      
    //  
    //  //divideBinWidth(hrec_sameside_resp_rebin);    
    //  hrec_sameside_resp_rebin->Write( (histTitle3+"_profile_rebin_divBinWidth").c_str());
    //  if(debugMode)std::cout<<"profile bin #1 pt low = "<<
    //		     (int)hrec_sameside_resp_rebin->GetBinLowEdge(1)<<std::endl;
    //  if(debugMode)std::cout<<"profile bin #"<<nbins_pt_gen<<" pt hi = "<<
    //		     hrec_sameside_resp_rebin->GetBinLowEdge(nbins_pt_reco)
    //		     + hrec_sameside_resp_rebin->GetBinWidth(nbins_pt_reco)<<std::endl;    
    //
    //  if(clearOverUnderflows){
    //	TH1clearOverUnderflows((TH1*)hrec_sameside_resp_rebin);
    //	hrec_sameside_resp_rebin->Write((histTitle3+"_profile_rebin_divBinWidth_noOverUnderflow").c_str());
    //	if(debugMode)hrec_sameside_resp_rebin->Print("base");    }
    //  
    //  //gen version here
    //  if(debugMode)std::cout<<"boundaries_pt_reco[0]  = "<< ((int)boundaries_pt_reco[0]) << std::endl;
    //  if(debugMode)std::cout<<"boudnaries_pt_reco["<<nbins_pt_reco<<"] = "<< (int)boundaries_pt_reco[nbins_pt_reco] <<std::endl;
    //  if(debugMode)std::cout<<"bin 1 pt low = "<<
    //		     (int)hmat->GetXaxis()->GetBinLowEdge((int)(boundaries_pt_reco[0]+1.))<<std::endl;
    //  
    //  if(debugMode)std::cout<<"bin "<<nbins_pt_reco<<" pt hi = "<<
    //		     hmat->GetXaxis()->GetBinLowEdge(((int)boundaries_pt_reco[nbins_pt_reco]))
    //		     + hmat->GetXaxis()->GetBinWidth(((int)boundaries_pt_reco[nbins_pt_reco]))<<std::endl;    
    //  
    //  //profile error types; default(sum/sqrt(N)), s, i, *g*
    //  hgen_resp_rebin = (TH1D*)hmat->TH2::ProjectionY( 
    //						   (histTitle3gen+"_profile").c_str() , 
    //						   (int)(boundaries_pt_reco[0]+1.), 
    //						   ((int)boundaries_pt_reco[nbins_pt_reco]), 	"g");
    //  hgen_resp_rebin=(TH1D*)hgen_resp_rebin->Rebin( nbins_pt_gen, (histTitle3gen+"_profile_rebin").c_str() , boundaries_pt_gen);
    //  hgen_resp_rebin->Write( (histTitle3gen+"_profile_rebin").c_str());      
    //  
    //  //divideBinWidth(hgen_resp_rebin);    
    //  hgen_resp_rebin->Write( (histTitle3gen+"_profile_rebin_divBinWidth").c_str());
    //  if(debugMode)std::cout<<"profile bin #1 pt low = "<<
    //		     (int)hgen_resp_rebin->GetBinLowEdge(1)<<std::endl;
    //  if(debugMode)std::cout<<"profile bin #"<<nbins_pt_gen<<" pt hi = "<<
    //		     hgen_resp_rebin->GetBinLowEdge(nbins_pt_gen)
    //		     + hgen_resp_rebin->GetBinWidth(nbins_pt_gen)<<std::endl;    
    //
    //  if(clearOverUnderflows){
    //	TH1clearOverUnderflows((TH1*)hgen_resp_rebin);
    //	hgen_resp_rebin->Write((histTitle3gen+"_profile_rebin_divBinWidth_noOverUnderflow").c_str());
    //	if(debugMode)hgen_resp_rebin->Print("base");    }
    //}
    

    //else{
    hrec_sameside_resp_rebin = (TH1D*)hrec_sameside_rebin->Clone( (histTitle3+"_clone").c_str() );      
    if(debugWrite)hrec_sameside_resp_rebin->Write( (histTitle3+"_clone").c_str() );    
    
    hgen_resp_rebin = (TH1D*)hgen_rebin->Clone( (histTitle3gen+"_clone").c_str() );      
    if(debugWrite)hgen_resp_rebin->Write( (histTitle3gen+"_clone").c_str() );       
    //}
    
  }
  else{
    hrec_sameside_resp_rebin = new TH1D( (histTitle3+"_empty").c_str() ,"", 
					 nbins_pt_reco, boundaries_pt_reco);  
    
    hgen_resp_rebin = new TH1D( (histTitle3gen+"_empty").c_str() ,"", 
				nbins_pt_gen, boundaries_pt_gen);     }

  //cosmetics
  hrec_sameside_resp_rebin->SetMarkerStyle(kOpenSquare);
  hrec_sameside_resp_rebin->SetMarkerColor(kBlue-3);     
  hrec_sameside_resp_rebin->SetLineColor(kBlue-3);     
  hrec_sameside_resp_rebin->SetMarkerSize(1.02);       

  if(debugWrite)hrec_sameside_resp_rebin->Write((histTitle3+"_final").c_str());
  if(debugMode)hrec_sameside_resp_rebin->Print(" base");  
  
  //cosmetics
  hgen_resp_rebin->SetMarkerStyle(kOpenStar);
  hgen_resp_rebin->SetMarkerColor(kMagenta);
  hgen_resp_rebin->SetLineColor(kMagenta);
  hgen_resp_rebin->SetMarkerSize(1.02);   
  
  if(debugWrite)hgen_resp_rebin->Write((histTitle3gen+"_final").c_str());
  if(debugMode)hgen_resp_rebin->Print("base");    
  
  
  
  
  
  // input histos -------------------------
  
  //for output
  if(debugWrite)fout->cd();    
  
  //HERE
  //ISSUE: NEED INDEPENDENT RECO SPECTRA w/ SEPARATE RANDOM # GEN TO PUT IN AS TEST SPECTRA
  //std::string histTitle="hpp_mcclosure_reco_test";
  std::string histTitle="smeared_rnd";//"hpp_mcclosure_reco_test";
  if(useNPCorrSpectra) histTitle+="_NPynew_test";
  else histTitle+="_ynew_test";

  
  TH1D*  hrec = (TH1D*)fin_ppMC->Get( histTitle.c_str() ); 
  hrec->Scale(1./NLOMCscaling);
  hrec->Scale(4.);
  multiplyBinWidth(hrec);
  if(debugWrite)hrec->Write();
  if(debugMode)hrec->Print("base");
  
  histTitle+="_clone";
  TH1D *hrec_rebin = (TH1D*)hrec->Clone( (histTitle).c_str() );
  std::cout<<"rebinning hrec..."<<std::endl;
  histTitle+="_rebins";
  hrec_rebin = (TH1D*)hrec_rebin->Rebin( nbins_pt_reco, (histTitle).c_str() , boundaries_pt_reco);
  if(debugWrite)hrec_rebin->Write(histTitle.c_str());   
  if(debugMode)hrec_rebin->Print("base");  
  
  
  if(clearOverUnderflows){
    histTitle+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hrec_rebin);
    if(debugWrite)hrec_rebin->Write(histTitle.c_str());
    if(debugMode)hrec_rebin->Print("base");  
  }
  
  //cosmetics
  hrec_rebin->SetMarkerStyle(kOpenCircle);
  hrec_rebin->SetMarkerColor(kBlue);     
  hrec_rebin->SetLineColor(kBlue);     
  hrec_rebin->SetMarkerSize(1.02);     

  
  
  
  // thy spectra  
  TH1D* CT10nlo  =(TH1D*)makeThyHist_00eta20_v2((fNLOFile_R04_CT10nlo  ).c_str(),useNPCorrSpectra);
  CT10nlo->SetMarkerSize(0);
  CT10nlo->SetLineColor(kBlack);  
  TH1D* CT14nlo  =(TH1D*)makeThyHist_00eta20_v2((fNLOFile_R04_CT14nlo  ).c_str(),useNPCorrSpectra);
  CT14nlo->SetMarkerSize(0);
  CT14nlo->SetLineColor(kGreen);  
  TH1D* HERAPDF  =(TH1D*)makeThyHist_00eta20_v2((fNLOFile_R04_HERAPDF  ).c_str(),useNPCorrSpectra);
  HERAPDF->SetMarkerSize(0);
  HERAPDF->SetLineColor(kViolet-5);  
  TH1D* MMHTnlo  =(TH1D*)makeThyHist_00eta20_v2((fNLOFile_R04_MMHTnlo  ).c_str(),useNPCorrSpectra);
  MMHTnlo->SetMarkerSize(0);
  MMHTnlo->SetLineColor(kOrange+7);  
  TH1D* NNPDFnnlo=(TH1D*)makeThyHist_00eta20_v2((fNLOFile_R04_NNPDFnnlo).c_str(),useNPCorrSpectra);
  NNPDFnnlo->SetMarkerSize(0);
  NNPDFnnlo->SetLineColor(kCyan-6);  

  
  
  
  
  
  // Bayesian unfolding -------------------------   
  // RooUnfoldResponse, Bayes --------------------
  std::cout<<"calling RooUnfoldResponse "<<std::endl;
  
  RooUnfoldResponse roo_resp( hrec_sameside_resp_rebin, hgen_resp_rebin, hmat_rebin, ("Response_matrix"+RandEtaRange).c_str()) ;
  //RooUnfoldResponse roo_resp( 0, 0, hmat_rebin, ("Response_matrix"+RandEtaRange).c_str()) ;
  roo_resp.UseOverflow(doOverUnderflows);    
  if(debugWrite)roo_resp.Write();

  TH1D* hfak=  (TH1D*) roo_resp.Hfakes() ;
  hfak->Print("base");
  //cosmetics
  hfak->SetMarkerColor(kGreen);
  hfak->SetLineColor(kGreen);
  
  
  std::cout<<"calling RooUnfoldBayes..."<<std::endl;
  RooUnfoldBayes unf_bayes( &roo_resp, hrec_rebin, kIterInput );
  unf_bayes.SetVerbose(verbosity);
  if(doToyErrs){
    std::cout<<"using toy errors, suppressing text output"<<std::endl;
    unf_bayes.SetNToys(10000);
    unf_bayes.SetVerbose(1);
  }
  std::cout<<"RooUnfoldBayes Overflow Status: " << unf_bayes.Overflow()<<std::endl;
  
  TH1D *hunf = (TH1D*)unf_bayes.Hreco(errorTreatment);     std::cout<<std::endl; 
  hunf->SetName("ppMC_OS_BayesUnf_Spectra");
  hunf->SetTitle(("OS "+NLOMCtitle_str+" Toy MC Unf.").c_str());
  if(debugMode)hunf->Print("base");
  
  //cosmetics
  hunf->SetMarkerStyle(kOpenCircle);
  hunf->SetMarkerColor(kRed);
  hunf->SetLineColor(kRed);
  hunf->SetMarkerSize(1.02);     
  
  std::cout<<"folding unfolded MC histogram!!"<<std::endl;
  TH1D* hfold=(TH1D*)roo_resp.ApplyToTruth(hunf);          
  hfold->SetName("ppMC_OS_BayesFold_Spectra");
  hfold->SetTitle(("OS "+NLOMCtitle_str+" Toy MC Fold(Unf.)").c_str());
  if(debugMode)hfold->Print("base");
  
  
  //cosmetics
  hfold->SetMarkerStyle(kOpenCircle);
  hfold->SetMarkerColor(kGreen-5);
  hfold->SetLineColor(  kGreen-5);
  hfold->SetMarkerSize(1.02);     




  
  std::cout<<"calling RooUnfoldBayes FOR SAME SIDE..."<<std::endl;
  RooUnfoldBayes unf_ss_bayes( &roo_resp, hrec_sameside_rebin, kIterInput );
  unf_ss_bayes.SetVerbose(verbosity);
  if(doToyErrs){
    std::cout<<"using toy errors, suppressing text output"<<std::endl;
    unf_ss_bayes.SetNToys(10000);
    unf_ss_bayes.SetVerbose(1);
  }
  std::cout<<"RooUnfoldBayes SAMESIDE Overflow Status: " << unf_ss_bayes.Overflow()<<std::endl;
  
  TH1D *hunf_ss = (TH1D*)unf_ss_bayes.Hreco(errorTreatment);     std::cout<<std::endl; 
  hunf_ss->SetName("ppMC_SS_BayesUnf_Spectra");
  //hunf_ss->SetTitle("SS Toy NLO MC Unf.");
  hunf_ss->SetTitle(("SS "+NLOMCtitle_str+" Toy MC Unf.").c_str());
  if(debugMode)hunf_ss->Print("base");
  
  //cosmetics
  hunf_ss->SetMarkerStyle(kOpenSquare);
  hunf_ss->SetMarkerColor(kRed);
  hunf_ss->SetLineColor(kRed);
  hunf_ss->SetMarkerSize(1.02);
     
  

  std::cout<<"folding unfolded same side Toy NLO MC histogram!!"<<std::endl;
  TH1D* hfold_ss=(TH1D*)roo_resp.ApplyToTruth(hunf_ss);          
  hfold_ss->SetName("ppMC_SS_BayesFold_Spectra");
  hfold_ss->SetTitle(("SS "+NLOMCtitle_str+" Toy MC Fold(Unf.)").c_str());
  //hfold_ss->SetTitle("SS Toy NLO MC Fold(Unf.)");
  if(debugMode)hfold_ss->Print("base");
  
  //cosmetics
  hfold_ss->SetMarkerStyle(kOpenSquare);
  hfold_ss->SetMarkerColor(kGreen-5);
  hfold_ss->SetLineColor(  kGreen-5);
  hfold_ss->SetMarkerSize(1.02);     
  




  
  
  std::cout<<"folding MC Truth histogram!!"<<std::endl;
  TH1D* hfold_truth=(TH1D*)roo_resp.ApplyToTruth(hgen_rebin);          
  hfold_truth->SetName("ppMC_SS_BayesFoldMCTruth_Spectra");
  hfold_truth->SetTitle(("Truth "+NLOMCtitle_str+" Toy MC Folded").c_str());
  hfold_truth->SetTitle("MC Truth Folded");
  if(debugMode)hfold_truth->Print("base");
  
  
  //cosmetics
  hfold_truth->SetMarkerStyle(kOpenStar);
  hfold_truth->SetMarkerColor(kGreen+3);
  hfold_truth->SetLineColor(  kGreen+3);
  hfold_truth->SetMarkerSize(1.02);     
  
  

  std::cout<<"calling Ereco... getting Cov Mat"<<std::endl;
  TMatrixD covmat = unf_bayes.Ereco(errorTreatment);
  
  std::cout<<"converting TMatrixD covmat to TH2D* covmat_TH2"<<std::endl;
  TH2D* covmat_TH2= new TH2D(covmat);  //this gives the TH2 a name of "TMatrixDBase" by default.
  covmat_TH2->SetName("Bayes_covarianceMatrix_TH2");
  if(debugMode)covmat_TH2->Print("base");
  if(debugWrite)covmat_TH2->Write();
  
  TH2D* covmatabsval_TH2=(TH2D*)absVal_TH2Content((TH2D*)covmat_TH2);
  
  std::cout<<"calculating pearson coefficients"<<std::endl;
  TH2D* PearsonBayes = (TH2D*)CalculatePearsonCoefficients( &covmat, false,"Bayes_pearson");
  if(debugMode)PearsonBayes->Print("base");
  if(debugWrite)PearsonBayes->Write();

  std::cout<<"getting unfolding matrix"<<std::endl;
  TH2D* unfmat_TH2=new TH2D(unf_bayes.UnfoldingMatrix());
  unfmat_TH2->SetName("Bayes_unfoldingMatrix_TH2");  
  if(debugMode)unfmat_TH2->Print("base");
  if(debugWrite)unfmat_TH2->Write();
  
  



  // --------- RATIOS WITH MC TRUTH ----------------
  TH1D *hgen_rebin_ratiodenom=(TH1D*)hgen_rebin->Clone("ppMC_Gen_Ratio_denom");
  if(debugMode)hgen_rebin_ratiodenom->Print("base");
  
  TH1D *h_genratio_oppunf = (TH1D*)hunf->Clone( "ppMC_Gen_Ratio_OppUnf" );
  h_genratio_oppunf->SetTitle( "OS Toy NLO MC Unf./SS Toy NLO MC Truth" );
  h_genratio_oppunf->Divide(hgen_rebin_ratiodenom);
  if(debugMode)h_genratio_oppunf->Print("base");

  TH1D *h_genratio_ssunf = (TH1D*)hunf_ss->Clone( "ppMC_Gen_Ratio_SSUnf" );
  h_genratio_ssunf->SetTitle( "SS Toy NLO MC Unf./SS Toy NLO MC Truth" );
  h_genratio_ssunf->Divide(hgen_rebin_ratiodenom);
  if(debugMode)h_genratio_ssunf->Print("base");
  
  TH1D *h_genratio_oppfold = (TH1D*)hfold->Clone( "ppMC_Gen_Ratio_OppFold" );
//  if(nbins_pt_reco>nbins_pt_gen)
//    h_genratio_oppfold=(TH1D*)h_genratio_oppfold->Rebin(nbins_pt_gen,"ppMC_Gen_Ratio_OppFold_rebin",boundaries_pt_gen);
  h_genratio_oppfold->SetTitle( "OS Toy NLO MC Fold(Unf.)/SS Toy NLO MC Truth" );
  h_genratio_oppfold->Divide(hgen_rebin_ratiodenom);
  if(debugMode)h_genratio_oppfold->Print("base");

  TH1D *h_genratio_ssfold = (TH1D*)hfold_ss->Clone( "ppMC_Gen_Ratio_SSFold" );
//  if(nbins_pt_reco>nbins_pt_gen)
//    h_genratio_oppfold=(TH1D*)h_genratio_oppfold->Rebin(nbins_pt_gen,"ppMC_Gen_Ratio_OppFold_rebin",boundaries_pt_gen);
  h_genratio_ssfold->SetTitle( "SS Toy NLO MC Fold(Unf.)/SS Toy NLO MC Truth" );
  h_genratio_ssfold->Divide(hgen_rebin_ratiodenom);
  if(debugMode)h_genratio_ssfold->Print("base");
  
  TH1D *h_genratio_oppmeas = (TH1D*)hrec_rebin->Clone( "ppMC_Gen_Ratio_OppMeas" );
//  if(nbins_pt_reco>nbins_pt_gen)
//    h_genratio_oppmeas=(TH1D*)h_genratio_oppmeas->Rebin(nbins_pt_gen,"ppMC_Gen_Ratio_Meas_rebin",boundaries_pt_gen);
  h_genratio_oppmeas->SetTitle( "OS Toy NLO MC Meas./SS Toy NLO MC Truth" );
  h_genratio_oppmeas->Divide(hgen_rebin_ratiodenom);
  if(debugMode)h_genratio_oppmeas->Print("base");

  
  TH1D *h_genratio_ssmeas = (TH1D*)hrec_sameside_rebin->Clone( "ppMC_Gen_Ratio_SSMeas" );
  h_genratio_ssmeas->SetTitle( "SS Toy NLO MC Meas./SS Toy NLO MC Truth" );
  h_genratio_ssmeas->Divide(hgen_rebin_ratiodenom);
  if(debugMode)h_genratio_ssmeas->Print("base");
  
  




  // ---------------- RATIOS WITH MC OS MEAS -----------------
  Double_t fac_min1=1.;
  Double_t fac=(hrec_rebin->Integral()/hrec_sameside_rebin->Integral());

  TH1D *hrec_rebin_fakecorr=(TH1D*)hrec_rebin->Clone("ppMC_OS_Meas_Ratio_denom");
  //hrec_sameside_rebin_fakecorr=(TH1D*)hrec_sameside_rebin_fakecorr->Rebin(nbins_pt_gen,"ppMC_OS_Meas_Ratio_denom_rebin",boundaries_pt_gen);
  if(debugMode)hrec_rebin_fakecorr->Print("base");
  
  TH1D *hrec_sameside_rebin_fakecorr=(TH1D*)hrec_sameside_rebin->Clone("ppMC_SS_Meas_Ratio_denom");
  //hrec_rebin_fakecorr=(TH1D*)hrec_rebin_fakecorr->Rebin(nbins_pt_gen,"ppMC_OS_Meas_Ratio_denom_rebin",boundaries_pt_gen);
  if(debugMode)hrec_sameside_rebin_fakecorr->Print("base");
  
  //fake correction
  int nbins_fakcorr=hrec_rebin_fakecorr->GetNbinsX();
  for(int i=1;i<=nbins_fakcorr; i++){
    Double_t hrec_i     = hrec_rebin_fakecorr->GetBinContent(i); Double_t hrec_i_err = hrec_rebin_fakecorr->GetBinError(i);
    Double_t hrecss_i     = hrec_sameside_rebin_fakecorr->GetBinContent(i); Double_t hrecss_i_err = hrec_sameside_rebin_fakecorr->GetBinError(i);
    Double_t hfak_i     = hfak->GetBinContent(i); Double_t hfak_i_err = hfak->GetBinError(i);
    hrec_rebin_fakecorr->SetBinContent( i, hrec_i-fac*hfak_i );     
    hrec_rebin_fakecorr->SetBinError(   i, sqrt(hrec_i_err*hrec_i_err + fac*hfak_i_err*fac*hfak_i_err));  
    hrec_sameside_rebin_fakecorr->SetBinContent( i, hrecss_i-fac_min1*hfak_i );     
    hrec_sameside_rebin_fakecorr->SetBinError(   i, sqrt(hrecss_i_err*hrecss_i_err + fac_min1*hfak_i_err*fac_min1*hfak_i_err));  
  }
  
  
  TH1D *h_recratio_oppunf = (TH1D*)hunf->Clone( "ppMC_OS_Meas_Ratio_OppUnf" );
  h_recratio_oppunf->SetTitle( "OS Toy NLO MC Unf./OS Toy NLO MC Meas." );
  h_recratio_oppunf->Divide(hrec_rebin);
  if(debugMode)h_recratio_oppunf->Print("base");
  
  TH1D *h_recratio_oppfold = (TH1D*)hfold->Clone( "ppMC_OS_Meas_Ratio_OppFold" );
  h_recratio_oppfold->SetTitle( "OS Toy NLO MC Fold(Unf.)/OS Toy NLO MC Meas." );
  h_recratio_oppfold->Divide(hrec_rebin_fakecorr);
  if(debugMode)h_recratio_oppunf->Print("base");
  
  TH1D *h_recratio_ssunf = (TH1D*)hunf_ss->Clone( "ppMC_OS_Meas_Ratio_SSUnf" );
  h_recratio_ssunf->SetTitle( "SS Toy NLO MC Unf./SS Toy NLO MC Meas." );
  h_recratio_ssunf->Divide(hrec_rebin);
  if(debugMode)h_recratio_oppunf->Print("base");
  
  TH1D *h_recratio_ssfold = (TH1D*)hfold->Clone( "ppMC_OS_Meas_Ratio_SSFold" );
  h_recratio_ssfold->SetTitle( "SS Toy NLO MC Fold(Unf.)/OS Toy NLO MC Meas." );
  h_recratio_ssfold->Divide(hrec_rebin_fakecorr);
  if(debugMode)h_recratio_ssunf->Print("base");

  TH1D *h_recratio_ssmeas = (TH1D*)hrec_sameside_rebin->Clone( "ppMC_OS_Meas_Ratio_SSMeas" );
  h_recratio_ssmeas->SetTitle( "SS Toy NLO MC Meas./OS Toy NLO MC Meas." );
  h_recratio_ssmeas->Divide(hrec_rebin);
  if(debugMode)h_recratio_ssmeas->Print("base");
 
  TH1D *h_recratio_ssgen = (TH1D*)hgen_rebin->Clone( "ppMC_OS_Meas_Ratio_Truth" );    
  h_recratio_ssgen->SetTitle( "SS Toy NLO MC Truth/OS Toy NLO MC Meas." );
  h_recratio_ssgen->Divide(hrec_rebin);
  if(debugMode)h_recratio_ssgen->Print("base");
  
  // ---------------- FOLDED RATIOS ----------------- //
  
  TH1D *h_foldratio_datafold=(TH1D*)hfold->Clone("ppMC_OS_Meas_Ratio_OppFold");
  h_foldratio_datafold->SetTitle( "OS Toy NLO MC Fold(Unf.)/OS Toy NLO MC Meas." );
  h_foldratio_datafold->Divide(hrec_rebin_fakecorr);
  if(debugMode)h_foldratio_datafold->Print("base");
    
  TH1D *h_foldratio_mcfold=(TH1D*)hfold_truth->Clone("ppMC_OS_Meas_Ratio_TruthFold");
  h_foldratio_mcfold->SetTitle( "SS Toy NLO MC Truth Fold/SS Toy NLO MC Meas." );
  h_foldratio_mcfold->Divide(hrec_sameside_rebin_fakecorr);
  if(debugMode)h_foldratio_mcfold->Print("base");

  TH1D *h_foldratio_ssfold=(TH1D*)hfold_ss->Clone("ppMC_OS_Meas_Ratio_SSFold");
  h_foldratio_ssfold->SetTitle( "SS Toy NLO MC Fold(Unf.)/SS Toy NLO MC Meas." );
  h_foldratio_ssfold->Divide(hrec_sameside_rebin_fakecorr);
  if(debugMode)h_foldratio_mcfold->Print("base");
  
  


  TH1D* hunf_bayes[nKiterMax]={};
  TH1D* hfold_bayes[nKiterMax]={};
  TH1D* hgen_unf_ratio[nKiterMax]={};
  TH1D* hrec_unf_ratio[nKiterMax]={};
  
  TH1D* hgen_fold_ratio[nKiterMax]={};
  TH1D* hrec_fold_ratio[nKiterMax]={};
  
  int kIter_start=kIterInput-kIterRange;
  int kIter_end=kIterInput+kIterRange;
  
  if(dokIterQA){
    if(debugMode)std::cout<<"unfolding across diff kIter values"<<std::endl;
    if(debugMode)std::cout<<kIter_start<<" <= kIter <= "<<kIter_end<<std::endl;
    for(int ki=0; ki<nKiterMax;ki++){
      int current_kIter=kIter_start+ki;
      if(debugMode)std::cout<<"ki="<<ki<<", current_kIter="<<current_kIter<<std::endl;
      
      RooUnfoldBayes unf_bayes_kIterQA(&roo_resp, hrec_rebin, current_kIter);        
      
      hunf_bayes[ki]=(TH1D*)unf_bayes_kIterQA.Hreco(errorTreatment);
      hunf_bayes[ki]->SetName(("ppMC_OS_unf_kIter"+std::to_string(current_kIter)).c_str());
      hunf_bayes[ki]->SetTitle(("OS MC Unf., kIter="+std::to_string(current_kIter)+";Jet p_{T}; #sigma" ).c_str());
      //cosmetics, keep here so they propagate to the ratio hists
      hunf_bayes[ki]->SetMarkerStyle(kOpenCircle);
      hunf_bayes[ki]->SetMarkerColor(kRed);
      hunf_bayes[ki]->SetLineColor(kRed);
      hunf_bayes[ki]->SetMarkerSize(1.02);     
      if(debugMode)hunf_bayes[ki]->Print("base");
      
      hfold_bayes[ki]=(TH1D*)roo_resp.ApplyToTruth(hunf_bayes[ki]);
      hfold_bayes[ki]->SetName(("ppMC_OS_fold_kIter"+std::to_string(current_kIter)).c_str());
      hfold_bayes[ki]->SetTitle(("OS MC Fold(Unf.), kIter="+std::to_string(current_kIter)+"; Jet p_{T}; #sigma").c_str());
      if(debugMode)hfold_bayes[ki]->Print("base");
      //cosmetics
      hfold_bayes[ki]->SetMarkerStyle(kOpenCircle);
      hfold_bayes[ki]->SetMarkerColor(kGreen-5);
      hfold_bayes[ki]->SetLineColor(  kGreen-5);
      hfold_bayes[ki]->SetMarkerSize(1.02);     
      
      if(debugMode)hfold_bayes[ki]->Print("base");
      
      hgen_unf_ratio[ki]=(TH1D*)hunf_bayes[ki]->Clone(("ppMC_OS_unf_genratio_kIter"+std::to_string(current_kIter)).c_str());
      hgen_unf_ratio[ki]->SetTitle(("OS MC Unf./SS MC Truth, kIter="+std::to_string(current_kIter)+"; Jet p_{T}; ratio").c_str());
      hgen_unf_ratio[ki]->Divide(hgen_rebin);
      
      hrec_unf_ratio[ki]=(TH1D*)hunf_bayes[ki]->Clone(("ppMC_OS_unf_recratio_kIter"+std::to_string(current_kIter)).c_str());
      //hrec_unf_ratio[ki]->SetTitle(("OS MC Unf./SS MC Fake Corr. Meas., kIter="+std::to_string(current_kIter)+"; Jet p_{T}; ratio").c_str());
      //hrec_unf_ratio[ki]->Divide(hrec_sameside_rebin_fakecorr);
      hrec_unf_ratio[ki]->SetTitle(("OS MC Unf./SS MC Meas., kIter="+std::to_string(current_kIter)+"; Jet p_{T}; ratio").c_str());
      hrec_unf_ratio[ki]->Divide(hrec_sameside_rebin);
      
      hgen_fold_ratio[ki]=(TH1D*)hfold_bayes[ki]->Clone(("ppMC_OS_fold_genratio_kIter"+std::to_string(current_kIter)).c_str());
      hgen_fold_ratio[ki]->SetTitle(("OS MC Fold./SS MC Truth, kIter="+std::to_string(current_kIter)+"; Jet p_{T}; ratio").c_str());
      hgen_fold_ratio[ki]->Divide(hgen_rebin);
      
      hrec_fold_ratio[ki]=(TH1D*)hfold_bayes[ki]->Clone(("ppMC_OS_fold_recratio_kIter"+std::to_string(current_kIter)).c_str());
      hrec_fold_ratio[ki]->SetTitle(("OS MC Fold./SS MC Fake Corr. Meas., kIter="+std::to_string(current_kIter)+"; Jet p_{T}; ratio").c_str());
      hrec_fold_ratio[ki]->Divide(hrec_sameside_rebin_fakecorr);
      
      //do bin width normalizations here, so the ratios in this loop arent affected by hgen/hrec hists not having been bin-width normalized yet
      hunf_bayes[ki]->Scale(1./etaBinWidth);
      divideBinWidth(hunf_bayes[ki]);
      
      hfold_bayes[ki]->Scale(1./etaBinWidth);
      divideBinWidth(hfold_bayes[ki]);
      
    }//end kIterLoop
  }//end if(dokIterQA)
  //assert(false);







  // ----- PUT BIN WIDTH(s) DIVISIONS + NORMALIZING HERE ----- //
  
  // -- MC SS RECO -- //  
  histTitle2+="_normbinwidth";
  hrec_sameside_rebin->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hrec_sameside_rebin); 
  if(debugWrite)hrec_sameside_rebin->Write( histTitle2.c_str() );
  if(debugMode)hrec_sameside_rebin->Print("base");  
  
  // -- MC SS GEN -- //  
  genHistTitle+="_normbinwidth";
  hgen_rebin->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hgen_rebin);
  if(debugWrite)hgen_rebin->Write(genHistTitle.c_str());
  if(debugMode)hgen_rebin->Print("base");  

  // -- MC SS GEN RESP -- //  
  hgen_resp_rebin->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hgen_resp_rebin);
  if(debugWrite)hgen_resp_rebin->Write(genHistTitle.c_str());
  if(debugMode)hgen_resp_rebin->Print("base");  

  // -- MC SS RESP FAKES -- //  
  hfak->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hfak);
  if(debugWrite)hfak->Write();
  if(debugMode)hfak->Print("base");  
  
  // -- MC SS REC RESP -- //  
  hrec_sameside_resp_rebin->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hrec_sameside_resp_rebin);
  if(debugWrite)hrec_sameside_resp_rebin->Write();
  if(debugMode)hrec_sameside_resp_rebin->Print("base");  
  
  // -- MC FOLD(MC SS TRUTH) -- //
  hfold_truth->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hfold_truth); 
  if(debugWrite)hfold_truth->Write();
  if(debugMode)hfold_truth->Print("base");  
  
  // -- MC OS RECO -- //  
  histTitle+="_normbinwidth";
  hrec_rebin->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hrec_rebin); 
  if(debugWrite)hrec_rebin->Write(histTitle.c_str());
  if(debugMode)hrec_rebin->Print("base");  
  
  // -- MC OS UNF -- //
  hunf->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hunf); 
  if(debugWrite)hunf->Write();
  if(debugMode)hunf->Print("base");  
  
  // -- MC OS FOLD(UNF) -- //
  hfold->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hfold); 
  if(debugWrite)hfold->Write();
  if(debugMode)hfold->Print("base");  

  // -- MC SS UNF -- //
  hunf_ss->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hunf_ss); 
  if(debugWrite)hunf_ss->Write();
  if(debugMode)hunf_ss->Print("base");  
  
  // -- MC SS FOLD(UNF) -- //
  hfold_ss->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hfold_ss); 
  if(debugWrite)hfold_ss->Write();
  if(debugMode)hfold_ss->Print("base");  



  // -- MC OS UNF x2 (for thy comparison) -- //
  TH1D* hunf_x2=(TH1D*)hunf->Clone( ( ((std::string)hunf->GetName()) + "_2xScaled").c_str() );
  //hunf_x2->Scale(2.);
  if(debugWrite)hunf_x2->Write();
  if(debugMode)hunf_x2->Print("base");
  
  // -- MC OS UNF x2 (for thy comparison) -- //
  TH1D* hgen_rebin_x2=(TH1D*)hgen_rebin->Clone( ( ((std::string)hgen_rebin->GetName()) + "_2xScaled").c_str() );
  //hgen_rebin_x2->Scale(2.);
  if(debugWrite)hgen_rebin_x2->Write();
  if(debugMode)hgen_rebin_x2->Print("base");


  // ----- END BIN WIDTH DIVISIONS + NORMALIZING ----- //
  
  // ---------------- THY RATIOS (THY / DATA) ----------------- //
  // must be after the binwidth divisions + normalization (the thy hists are made this way by default)
  
  TH1D* h_thyratio_CT10nlo  =(TH1D*)CT10nlo  ->Clone("");
  h_thyratio_CT10nlo = (TH1D*)h_thyratio_CT10nlo->Rebin(nbins_pt_gen,"pp_CT10Thy_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_CT10nlo  ->SetTitle("CT10 NLO/OS Toy NLO MC Unf.");  
  h_thyratio_CT10nlo  ->Divide(hunf_x2);
  if(debugMode)h_thyratio_CT10nlo  ->Print("base");
  
  TH1D* h_thyratio_CT14nlo  =(TH1D*)CT14nlo  ->Clone("");  
  h_thyratio_CT14nlo=(TH1D*)h_thyratio_CT14nlo->Rebin(nbins_pt_gen,"pp_CT14Thy_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_CT14nlo  ->SetTitle("CT14 NLO/OS Toy NLO MC Unf.");
  h_thyratio_CT14nlo  ->Divide(hunf_x2);
  if(debugMode)h_thyratio_CT14nlo  ->Print("base");
  
  TH1D* h_thyratio_HERAPDF  =(TH1D*)HERAPDF  ->Clone("");  
  h_thyratio_HERAPDF=(TH1D*)h_thyratio_HERAPDF->Rebin(nbins_pt_gen,"pp_HERAPDF_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_HERAPDF  ->SetTitle("HERAPDF15 NLO/OS Toy NLO MC Unf.");
  h_thyratio_HERAPDF  ->Divide(hunf_x2);
  if(debugMode)h_thyratio_HERAPDF  ->Print("base");
  
  TH1D* h_thyratio_MMHTnlo  =(TH1D*)MMHTnlo  ->Clone("");
  h_thyratio_MMHTnlo=(TH1D*)h_thyratio_MMHTnlo->Rebin(nbins_pt_gen,"pp_MMHT_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_MMHTnlo  ->SetTitle("MMHT14 NLO/OS MC Unf.");
  h_thyratio_MMHTnlo  ->Divide(hunf_x2);
  if(debugMode)h_thyratio_MMHTnlo  ->Print("base");
  
  TH1D* h_thyratio_NNPDFnnlo=(TH1D*)NNPDFnnlo->Clone("");
  h_thyratio_NNPDFnnlo=(TH1D*)h_thyratio_NNPDFnnlo->Rebin(nbins_pt_gen,"pp_NNPDFnnlo_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_NNPDFnnlo->SetTitle("NNPDF NNLO/OS MC Unf.");
  h_thyratio_NNPDFnnlo->Divide(hunf_x2);
  if(debugMode)h_thyratio_NNPDFnnlo->Print("base");

  //  TH1D* h_thyratio_mctruth=(TH1D*)NNPDFnnlo->Clone("");
  TH1D* h_thyratio_mctruth=(TH1D*)hgen_rebin->Clone("");
  h_thyratio_mctruth=(TH1D*)h_thyratio_mctruth->Rebin(nbins_pt_gen,"pp_NNPDFnnlo_ToyMCTruth_Ratio_rebin",boundaries_pt_gen);
  //h_thyratio_mctruth->SetTitle("PY8 GEN/Data Unf.");  
  h_thyratio_mctruth->SetTitle("NNPDF NNLO Toy MC Truth/OS MC Unf.");
  h_thyratio_mctruth->Divide(hunf_x2);
  if(debugMode)h_thyratio_mctruth->Print("base");

  

  
  
  
    

  
  
  //  drawPDFS -------------------------------------------------- 

  //TCanvas pointers for writing canvs to file
  TCanvas *canv_spectra=NULL, *canv_mc_fakes_spectra=NULL, *canv_thy_spectra_1=NULL, *canv_thy_spectra_2=NULL;
  TCanvas *canv_gen_ratio=NULL, *canv_rec_ratio=NULL, *canv_fold_ratio=NULL, *canv_thy_ratio=NULL; 
  TCanvas *canv_covmat=NULL, *canv_absval_covmat=NULL, *canv_pearson=NULL, *canv_unfmat=NULL, *canv_mat_rebin=NULL, *canv_mat_percerrs=NULL;
  
  TH2D* hmat_percenterrs=NULL;
  TCanvas *canv_3x3spectra=NULL, *canv_3x3genratio=NULL, *canv_3x3recratio=NULL;

  //misc printing prep
  TLine* theLineAtOne= new TLine( boundaries_pt_gen_mat[0], 1., 
				  (boundaries_pt_gen_mat[nbins_pt_gen_mat]), 1.);
  theLineAtOne->SetLineWidth(1);
  theLineAtOne->SetLineStyle(2);
  theLineAtOne->SetLineColor(36);    
  
  TLine* theLineAtp9= (TLine*)theLineAtOne->Clone();
  theLineAtp9->SetY1(0.9);    theLineAtp9->SetY2(0.9);       
  
  TLine* theLineAt1p1= (TLine*)theLineAtOne->Clone();
  theLineAt1p1->SetY1(1.1);    theLineAt1p1->SetY2(1.1);
  
  
  if(!drawPDFs)std::cout<<std::endl<<"NOT drawing PDFs for Bayesian Unfolding!"<<std::endl<<std::endl;
  else{ std::cout<<std::endl<<"drawing output PDFs for Bayesian Unfolding..."<<std::endl;
    std::string outPdfFile=outBayesPdfFile;
    std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";
    
    // canvas for printing
    TCanvas* canvForPrint=new TCanvas("printCanvas","for printing",           1400,1000);    
    canvForPrint->cd();
    //canvForPrint->SetGridx(1);    //canvForPrint->SetGridy(1);
    canvForPrint->Print(open_outPdfFile.c_str()); 
    
    

    // gen/meas/unf spectra hists ---------------    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(1);    
    
    setupSpectraHist(hgen_rebin          , useSimpBins);    
    setupSpectraHist(hrec_rebin	       , useSimpBins);
    setupSpectraHist(hrec_sameside_rebin , useSimpBins);
    setupSpectraHist(hunf		       , useSimpBins);    
    setupSpectraHist(hunf_ss		       , useSimpBins);    
    
    //hrec_rebin->SetTitle("Bayesian, OS and SS Toy NLO MC Spectra");
    hrec_rebin->SetTitle(("Bayesian, OS and SS "+NLOMCtitle_str+"Toy MC Spectra").c_str());
    
    hrec_rebin->DrawClone("P E");           
    hrec_sameside_rebin->DrawClone("P E SAME");           
    hgen_rebin->DrawClone("P E SAME");                 
    hunf->DrawClone("P E SAME");               
    hunf_ss->DrawClone("P E SAME");
    
    TLegend* legend_in1 = new TLegend( 0.7,0.7,0.9,0.9 );
    legend_in1->AddEntry(hrec_sameside_rebin, "SS Toy NLO MC Meas."  ,  "lp");
    legend_in1->AddEntry(hrec_rebin,          "OS Toy NLO MC Meas." , "lp");	
    legend_in1->AddEntry(hunf,                "OS Toy NLO MC Unf." ,  "lp");
    legend_in1->AddEntry(hunf_ss,             "SS Toy NLO MC Unf." ,  "lp");
    legend_in1->AddEntry(hgen_rebin,          "SS Toy NLO MC Truth"   , "lp");
    
    legend_in1->Draw();
    
    canvForPrint->Print(outPdfFile.c_str());      
        canv_spectra=(TCanvas*)canvForPrint->DrawClone();    
    
    // ratios w MC Truth------------    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(0);
    
    setupRatioHist(h_genratio_oppunf , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_genratio_ssunf , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_genratio_oppmeas, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_genratio_ssmeas , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    //setupRatioHist(h_genratio_oppmeas, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    //setupRatioHist(h_genratio_oppfold, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    //setupRatioHist(h_genratio_ssmeas, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);    
    
    //h_genratio_oppunf->SetTitle( "Bayesian, Ratios w/ SS Toy NLO MC Truth Spectra" );
    h_genratio_oppunf->SetTitle(("Bayes, "+NLOMCtitle_str+" Toy MC, Ratio w/ SS Toy NLO MC Truth").c_str());
    h_genratio_oppunf->GetYaxis()->SetTitle("Ratio w/ SS Toy NLO MC Truth");
    

    h_genratio_oppunf->DrawClone("P E");            // os mc unf/mc truth
    h_genratio_oppmeas->DrawClone("P E SAME");      // os mc meas/mc truth
    h_genratio_ssmeas->DrawClone("P E SAME");       // ss mc meas/mc truth
    h_genratio_ssunf->DrawClone("P E SAME");       // ss mc meas/mc truth
    //h_genratio_oppfold->DrawClone("P E SAME");
    
    TLegend* legend2 = new TLegend( 0.1,0.8,0.2,0.9 );
    legend2->AddEntry(h_genratio_oppunf,  "OS Toy NLO MC Unf." ,  "lp");
    legend2->AddEntry(h_genratio_ssunf,  "SS Toy NLO MC Unf." ,  "lp");
    legend2->AddEntry(h_genratio_oppmeas, "OS Toy NLO MC Meas." , "lp"); 
    legend2->AddEntry(h_genratio_ssmeas, "SS Toy NLO MC Meas.", "lp");
    //legend2->AddEntry(h_genratio_oppfold, "OS Toy NLO MC Fold(Unf.)", "lp");
    
    legend2->Draw();
    

    theLineAtp9 ->Draw();
    theLineAtOne->Draw();
    theLineAt1p1->Draw();
    
    canvForPrint->Print(outPdfFile.c_str());
    canv_gen_ratio=(TCanvas*)canvForPrint->DrawClone();

    // ratio w/ data meas -------------------
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(0);      
    
    //h_recratio_oppunf->SetTitle("Bayesian, Ratios w/ OS Toy NLO MC Meas. Spectra");
    //h_recratio_oppunf->GetYaxis()->SetTitle("Ratio w/ OS Toy NLO MC Meas.");

    //h_recratio_ssmeas->SetTitle("Bayesian, Ratios w/ OS Toy NLO MC Meas. Spectra");
    h_recratio_ssmeas->SetTitle(("Bayesian, "+ NLOMCtitle_str+" Toy MC, Ratios w/ OS Toy NLO MC Meas. Spectra").c_str());
    

    h_recratio_ssmeas->GetYaxis()->SetTitle("Ratio w/ OS Toy NLO MC Meas.");
    
    setupRatioHist(h_recratio_ssmeas , useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);    
    setupRatioHist(h_recratio_oppunf , useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);    
    setupRatioHist(h_recratio_ssunf , useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);    
    setupRatioHist(h_recratio_ssgen  , useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);    
    
    h_recratio_ssmeas->DrawClone("P E");
    h_recratio_oppunf->DrawClone("P E SAME");
    h_recratio_ssunf->DrawClone("P E SAME");
    h_recratio_ssgen->DrawClone("P E SAME");
    //h_recratio_oppfold->DrawClone("P E SAME");  
    //h_recratio_ssfold->DrawClone("P E SAME");
    
    TLegend* legend4 = new TLegend( 0.1,0.8,0.2,0.9 );
    legend4->AddEntry(h_recratio_ssmeas, "SS Toy NLO MC Meas.", "lp");
    legend4->AddEntry(h_recratio_oppunf,  "OS Toy NLO MC Unf.", "lp");
    legend4->AddEntry(h_recratio_ssgen,   "SS Toy NLO MC Truth", "lp");
    legend4->AddEntry(h_recratio_ssunf, "SS Toy NLO MC Unf.", "p");
    //legend4->AddEntry(h_recratio_oppfold, "OS Toy NLO MC Fold(Unf.)", "lp");
    //legend4->AddEntry(h_recratio_ssfold, "MC Fold(Unf.)", "p");
    legend4->Draw();


    theLineAtp9 ->Draw();
    theLineAtOne->Draw();
    theLineAt1p1->Draw();
    

    
    canvForPrint->Print(outPdfFile.c_str());
    canv_rec_ratio=(TCanvas*)canvForPrint->DrawClone();

    // Folding Ratios -------------------
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(0);
    
    setupRatioHist(h_foldratio_datafold, useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);
    setupRatioHist(h_foldratio_mcfold, useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);
    

    h_foldratio_datafold->SetTitle( (NLOMCtitle_str+" Toy MC, Folded Ratios w/ (Corrected) Toy NLO MC Meas. Spectra" ).c_str());
    h_foldratio_datafold->GetYaxis()->SetTitle("Ratio w/ Toy NLO MC Meas.");
    
    h_foldratio_datafold->DrawClone("P E");
    h_foldratio_mcfold->DrawClone("P E SAME");
    h_foldratio_ssfold->DrawClone("P E SAME");
    
    TLegend* legendfold = new TLegend( 0.1,0.8,0.4,0.9, NULL, "NBNDC" );
    
    legendfold->AddEntry(h_foldratio_datafold ,  "OS Toy NLO MC Fold(Unf.)/OS Toy NLO MC Meas." , "lp");
    legendfold->AddEntry(h_foldratio_mcfold  ,   "SS Toy NLO MC Truth Fold/SS Toy NLO MC Meas.",  "lp");
    legendfold->AddEntry(h_foldratio_ssfold  ,   "SS Toy NLO MC Fold(Unf.)/SS Toy NLO MC Meas.",  "lp");
    
    legendfold->Draw();
    
    theLineAtp9 ->Draw();
    theLineAtOne->Draw();
    theLineAt1p1->Draw();


    canvForPrint->Print(outPdfFile.c_str());    
    canv_fold_ratio=(TCanvas*)canvForPrint->DrawClone();

    
    // gen/reco response spectra hists ---------------    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(1);    
    
    setupSpectraHist(hrec_sameside_resp_rebin , useSimpBins);
    setupSpectraHist(hrec_rebin , useSimpBins);
    //setupSpectraHist(hgen_resp_rebin          , useSimpBins);
    setupSpectraHist(hfak                     , useSimpBins);
    
    //histograms given to roo_resp + the fakes histo made from these spectra + response matrix projection
    //hrec_sameside_resp_rebin->SetTitle("MC Only, Response Spectra and Fakes");
    hfak->SetTitle((NLOMCtitle_str+" Toy MC Only, Spectra and Kinematic Fakes").c_str());

    //hrec_sameside_resp_rebin->GetYaxis()->SetTitle("N_{Jets}/L_{int}");    
    
    hfak->DrawClone("P E");           
    hrec_sameside_resp_rebin->DrawClone("P E SAME");           
    hrec_rebin->DrawClone("P E SAME");
    hgen_resp_rebin->DrawClone("P E SAME");                 
    
    TLegend* legend_resp = new TLegend( 0.6,0.7,0.9,0.9 , NULL, "NBNDC");
    legend_resp->AddEntry(hfak,          "SS Toy NLO MC Kinematic Fakes" , "lp");    
    legend_resp->AddEntry(hgen_resp_rebin,          "SS Toy NLO MC Truth" , "lp");
    legend_resp->AddEntry(hrec_sameside_resp_rebin, "SS Toy NLO MC Meas." , "lp");
    legend_resp->AddEntry(hrec_rebin, "OS Toy NLO MC Meas." , "lp");
    
    legend_resp->Draw();
    
    canv_mc_fakes_spectra=(TCanvas*)canvForPrint->DrawClone();    
    canvForPrint->Print(outPdfFile.c_str());       


    //// gen/reco response spectra hists ---------------    
    //canvForPrint->cd();
    //if(!useSimpBins)canvForPrint->SetLogx(1);
    //canvForPrint->SetLogy(1);    
    //
    //setupSpectraHist(hfold       , useSimpBins);
    //setupSpectraHist(hfold_truth , useSimpBins);
    //
    ////histograms given to roo_resp + the fakes histo made from these spectra + response matrix projection
    //hfold->SetTitle("Folded Spectra, Data and MC");
    //
    //hfold->DrawClone("P E");           
    //hfold_truth->DrawClone("P E SAME");           
    //
    //TLegend* legend_fold = new TLegend( 0.6,0.7,0.9,0.9 );
    //legend_fold->AddEntry(hfold,       "Data Fold(Unf.)" , "lp");    
    //legend_fold->AddEntry(hfold_truth, "MC Fold(Truth)"        , "lp");    
    //legend_fold->Draw();
    //
    //canvForPrint->Print(outPdfFile.c_str());      
   
    
    //    // MC gen/meas/unf spectra hists ---------------
    //    
    //    canvForPrint->cd();
    //    if(!useSimpBins)canvForPrint->SetLogx(1);
    //    canvForPrint->SetLogy(1);    
    //    
    //    setupSpectraHist(hgen_rebin, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);    
    //    setupSpectraHist(hrec_sameside_rebin , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    //    setupSpectraHist(hfold_sameside_clone	  , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    //    setupSpectraHist(hunf_sameside	  , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    //    
    //    hrec_sameside_rebin->SetTitle("MC Spectra Only");
    //    
    //    hrec_sameside_rebin->DrawClone("P E");           
    //    hgen_rebin->DrawClone("P E SAME");                 
    //    hunf_sameside->DrawClone("P E SAME");               
    //    //hfold_sameside_clone->DrawClone("P E SAME");               
    //    
    //    TLegend* legend_in2 = new TLegend( 0.6,0.7,0.9,0.9 );
    //    
    //    legend_in2->AddEntry(hunf_sameside,        "MC Unf." , "lp");
    //    legend_in2->AddEntry(hgen_rebin,          "MC Truth" , "lp");
    //    legend_in2->AddEntry(hrec_sameside_rebin, "MC Meas." , "lp");
    //    //legend_in2->AddEntry(hfold_sameside_clone, "Py8 Fold.", "lp");    
    //    
    //    legend_in2->Draw();
    //    
    //    canvForPrint->Print(outPdfFile.c_str());      
    
    
    // thy spectra CT10/14 NNPDF NLO---------------------------
    TLegend* legendThy1 =new TLegend( 0.7,0.7,0.9,0.9 );    
    
    //CT10nlo->SetLineColor(kRed);
    setupSpectraHist(CT10nlo  ,useSimpBins);
    legendThy1->AddEntry(CT10nlo  ,"CT10 NLO","l");
    if(debugWrite){fout->cd(); CT10nlo->Write("CT10_NLO_R04_jtpt");}
    
    //CT14nlo->SetLineColor(kGreen);
    setupSpectraHist(CT14nlo  ,useSimpBins);
    legendThy1->AddEntry(CT14nlo  ,"CT14 NLO","l");
    if(debugWrite){fout->cd(); CT14nlo->Write("CT14_NLO_R04_jtpt");}
    
    if(debugWrite)fout->cd();//makeThyHist will make me cd to last thy file opened...
    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(1);                  
    
    legendThy1->AddEntry(hunf_x2,"(OS Toy NLO MC Unf.)","lp");
    legendThy1->AddEntry(hgen_rebin_x2,"(SS Toy NLO MC Truth)", "lp");
        
    CT10nlo->SetAxisRange(1e-07,1e+03,"Y");//for y axis in nanbarns
    CT10nlo->SetTitle("NLO Thy w/ MC, OS Unf., SS Truth");

    CT10nlo  ->DrawClone("][HIST E");
    CT14nlo  ->DrawClone("][HIST E SAME");
    
    hgen_rebin_x2->DrawClone("P E SAME");
    hunf_x2->DrawClone("P E SAME");   //just for axis range
    

    legendThy1->Draw();

    canv_thy_spectra_1=(TCanvas*)canvForPrint->DrawClone();
    canvForPrint->Print(outPdfFile.c_str());

    
    // thy spectra HERA/MMHT/NNPDF ---------------------------

    TLegend* legendThy =new TLegend( 0.7,0.7,0.9,0.9 );    

    //HERAPDF->SetLineColor(kBlue);
    setupSpectraHist(HERAPDF  ,useSimpBins);
    legendThy->AddEntry(HERAPDF  ,"HERAPDF 2015 NLO","l");
    if(debugWrite){fout->cd(); HERAPDF->Write("HERAPDF105_NLO_R04_jtpt");}
    
    //MMHTnlo->SetLineColor(kMagenta);
    setupSpectraHist(MMHTnlo  ,useSimpBins);
    legendThy->AddEntry(MMHTnlo  ,"MMHT 2014 NLO   ","l");
    if(debugWrite){fout->cd(); MMHTnlo->Write("MMHT2014_NLO_R04_jtpt"); }
    
    //NNPDFnnlo->SetLineColor(kBlack);
    setupSpectraHist(NNPDFnnlo,useSimpBins);
    legendThy->AddEntry(NNPDFnnlo,"NNPDF NNLO","l");
    if(debugWrite){fout->cd(); NNPDFnnlo->Write("NNPDF_NLO_R04_jtpt");    }
    
    if(debugWrite)fout->cd();//makeThyHist will make me cd to last thy file opened...
    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(1);                  
    
    legendThy->AddEntry(hunf_x2,"(OS Toy NLO MC Unf.)","lp");
    legendThy->AddEntry(hgen_rebin_x2,"(SS Toy NLO MC Truth)", "lp");
    
    HERAPDF->SetAxisRange(1e-07,1e+03,"Y");//for y axis in nanbarns
    HERAPDF->SetTitle("NLO Thy w/ MC, OS Unf., SS Truth");
    
    HERAPDF  ->DrawClone("][HIST E");
    MMHTnlo  ->DrawClone("][HIST E SAME");
    NNPDFnnlo->DrawClone("][HIST E SAME");
    
    hgen_rebin_x2->DrawClone("P E SAME");
    hunf_x2->DrawClone("P E SAME");   //just for axis range
    

    legendThy->Draw();

    canvForPrint->Print(outPdfFile.c_str());
    canv_thy_spectra_2=(TCanvas*)canvForPrint->DrawClone();

    // thy ratios w hunf------------    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(0);
    
    setupRatioHist(h_thyratio_CT10nlo , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_thyratio_CT14nlo , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_thyratio_HERAPDF , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_thyratio_MMHTnlo , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_thyratio_NNPDFnnlo, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_thyratio_mctruth, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);

    h_thyratio_CT10nlo->SetTitle( "Thy Ratios w/ (Bayes Unf. OS MC)" );
    h_thyratio_CT10nlo->GetYaxis()->SetTitle("Thy / (Bayes Unf. OS MC)");
    
    h_thyratio_CT10nlo ->DrawClone( "][HIST ");      
    h_thyratio_CT14nlo ->DrawClone( "][HIST SAME"); 
    //h_thyratio_HERAPDF ->DrawClone( "][HIST SAME"); 
    //h_thyratio_MMHTnlo ->DrawClone( "][HIST SAME"); 
    h_thyratio_NNPDFnnlo->DrawClone("][HIST E SAME"); 
    h_thyratio_mctruth->DrawClone("P E SAME"); 
    
    
    TLegend* legendthyrat = new TLegend( 0.1,0.7,0.3,0.9 );
    legendthyrat->AddEntry(h_thyratio_CT10nlo ,  "CT10 PDF NLO" ,    "l");
    legendthyrat->AddEntry(h_thyratio_CT14nlo ,  "CT14 PDF NLO" ,    "l"); 
    //legendthyrat->AddEntry(h_thyratio_HERAPDF ,  "HERAPDF 2015 NLO", "l");
    //legendthyrat->AddEntry(h_thyratio_MMHTnlo ,  "MMHT 2014 NLO",    "l");
    legendthyrat->AddEntry(h_thyratio_NNPDFnnlo, "NNPDF NNLO",       "l");
    legendthyrat->AddEntry(h_thyratio_mctruth, "NNPDF NNLO Toy MC Truth",       "l");
    
    legendthyrat->Draw();
    
    theLineAtp9 ->Draw();
    theLineAtOne->Draw();
    theLineAt1p1->Draw();


    canvForPrint->Print(outPdfFile.c_str());
    canv_thy_ratio=(TCanvas*)canvForPrint->DrawClone();

    
    

    


    
    // NJets v. Jet p_T hist (all weight=1) ---------------------
    
//    if( (bool)hJetQA_jtptEntries )      {
//      TLine* theLineAt100_reco= new TLine( boundaries_pt_reco_mat[0]  
//					   ,100.   
//					   ,(boundaries_pt_reco_mat[nbins_pt_reco_mat])
//					   ,100.);
//      if(debugMode)std::cout<<"theLineAt1p1_reco starts at "<<boundaries_pt_reco_mat[0]  <<std::endl;
//      if(debugMode)std::cout<<"theLineAt1p1_reco ends at   "<<(boundaries_pt_reco_mat[nbins_pt_reco_mat])  <<std::endl;
//      theLineAt100_reco->SetLineWidth(1);
//      theLineAt100_reco->SetLineStyle(2);
//      theLineAt100_reco->SetLineColor(36);
//      
//      canvForPrint->cd();
//      if(!useSimpBins)canvForPrint->SetLogx(1);
//      canvForPrint->SetLogy(1);
//      
//      setupSpectraHist(hJetQA_jtptEntries, useSimpBins);
//      hJetQA_jtptEntries->GetYaxis()->SetTitle("NJets/Bin [weights=1]");
//      
//      hJetQA_jtptEntries->DrawClone("HIST E");
//
//      theLineAt100_reco->Draw();
//
//      canvForPrint->Print(outPdfFile.c_str());
//      
//    }
    

    // covariance matrix ---------------      
    canvForPrint->cd();
    canvForPrint->SetLogx(0);
    canvForPrint->SetLogy(0);
    canvForPrint->SetLogz(1);
      
    //always use simp bins for covmat to avoid log scaling the x/y axes
    matStylePrint( (TH2D*)covmat_TH2, "Bayesian Unfolding Covariance Matrix", canvForPrint, outPdfFile, true);
    canv_covmat=(TCanvas*)canvForPrint->DrawClone(); 
    
    //always use simp bins for covmat to avoid log scaling the x/y axes
    matStylePrint( (TH2D*)covmatabsval_TH2, "Bayesian Unfolding |Covariance Matrix|", canvForPrint, outPdfFile, true);
    canv_absval_covmat=(TCanvas*)canvForPrint->DrawClone(); 
    
    // pearson matrix ---------------      
    canvForPrint->cd();
    canvForPrint->SetLogx(0);
    canvForPrint->SetLogy(0);
    canvForPrint->SetLogz(0);
    
    //always use simp bins for covmat to avoid log scaling the x/y axes
    matStylePrint( (TH2D*)PearsonBayes, "Bayesian Unfolding Pearson Matrix", canvForPrint, outPdfFile, true);
    canv_pearson=(TCanvas*)canvForPrint->DrawClone(); 
    
    // unfolding matrix ---------------      
    canvForPrint->cd();
    canvForPrint->SetLogx(0);
    canvForPrint->SetLogy(0);
    canvForPrint->SetLogz(1);
    
    //always use simp bins for covmat to avoid log scaling the x/y axes
    matStylePrint( (TH2D*)unfmat_TH2, "Bayesian Unfolding Matrix", canvForPrint, outPdfFile, true);
    canv_unfmat=(TCanvas*)canvForPrint->DrawClone(); 


    // matrix rebinned ---------------    
    canvForPrint->cd();
    canvForPrint->SetLogz(1);
    
    if(!useSimpBins){
      canvForPrint->SetLogx(1);    
      canvForPrint->SetLogy(1);}
    else{
      canvForPrint->SetLogx(0);    
      canvForPrint->SetLogy(0);          }
    
    matStylePrint(hmat_rebin, "MC Response Matrix, Rebinned", canvForPrint, outPdfFile, useSimpBins);
    canv_mat_rebin=(TCanvas*)canvForPrint->DrawClone(); 
    


    // percent error matrix in binning of interest ---------------    
    hmat_percenterrs= makeRespMatrixPercentErrs( (TH2D*) hmat_errors, (TH2D*) hmat_rebin,
						       (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
						       (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );		     
    hmat_percenterrs->SetName("hmat_percenterrs_foroutputpdf");
    if(debugWrite)hmat_percenterrs->Write();
    if(debugMode) hmat_percenterrs->Print("base");
    
    canvForPrint->cd();
    canvForPrint->SetLogz(1);
    
    if(!useSimpBins){
      canvForPrint->SetLogx(1);    
      canvForPrint->SetLogy(1);}
    else{
      canvForPrint->SetLogx(0);    
      canvForPrint->SetLogy(0);          }
    
    matStylePrint(hmat_percenterrs, "MC Response Matrix, Bin % Errors", canvForPrint, outPdfFile, useSimpBins);
    canv_mat_percerrs=(TCanvas*)canvForPrint->DrawClone();

    //--------------- 3x3 canvases
    if(dokIterQA){
      if(debugMode)std::cout<<"drawing plots across diff kIter values"<<std::endl;
      if(debugMode)std::cout<<kIter_start<<" <= kIter <= "<<kIter_end<<std::endl;
      canv_3x3spectra  =new TCanvas("canv_3x3spectra","canv 3x3 spectra", 1500, 1500);   canv_3x3spectra ->Divide(3,3);
      canv_3x3genratio =new TCanvas("canv_3x3genratio","canv 3x3 genratio", 1500, 1500); canv_3x3genratio->Divide(3,3);
      canv_3x3recratio =new TCanvas("canv_3x3recratio","canv 3x3 recratio", 1500, 1500); canv_3x3recratio->Divide(3,3);            
      TLegend* leg_3x3=new TLegend(0.62,0.75,0.9,0.9,NULL,"NBNDC");
      
      for(int ki=0; ki<nKiterMax;ki++){
	int current_kIter=kIter_start+ki;
	int ki_canv=ki+1;
	if(debugMode)std::cout<<"ki="<<ki<<", current_kIter="<<current_kIter<<std::endl;
 
	canv_3x3spectra->cd(ki_canv)->SetLogx(1);
	canv_3x3spectra->cd(ki_canv)->SetLogy(1);
	canv_3x3spectra->cd(ki_canv);
 
	hunf_bayes[ki]->DrawClone("HIST E");
	hfold_bayes[ki]->DrawClone("HIST E SAME");
 
	leg_3x3->Clear();//removes any entries in legend
	leg_3x3->AddEntry(hunf_bayes[ki] ,"Unf."      ,"lp");
	leg_3x3->AddEntry(hfold_bayes[ki],"Fold(Unf.)","lp");
	leg_3x3->DrawClone();
 
	canv_3x3genratio->cd(ki_canv)->SetLogx(1);
	canv_3x3genratio->cd(ki_canv)->SetLogy(0);
	canv_3x3genratio->cd(ki_canv);
	hgen_unf_ratio[ki]->SetAxisRange(0.5,1.5,"Y");
 
	hgen_unf_ratio[ki]->DrawClone("HIST E");
	hgen_fold_ratio[ki]->DrawClone("HIST E SAME");
 
	leg_3x3->Clear();//removes any entries in legend
	leg_3x3->AddEntry(hgen_unf_ratio[ki] ,"Unf."      ,"lp");
	leg_3x3->AddEntry(hgen_fold_ratio[ki],"Fold(Unf.)","lp");
	leg_3x3->DrawClone();
 
	canv_3x3recratio->cd(ki_canv)->SetLogx(1);
	canv_3x3recratio->cd(ki_canv)->SetLogy(0);
	canv_3x3recratio->cd(ki_canv);
 
	hrec_unf_ratio[ki]->SetAxisRange(0.5,1.5,"Y");
 
	hrec_unf_ratio[ki]->DrawClone("HIST E");
	hrec_fold_ratio[ki]->DrawClone("HIST E SAME"); 
 
	leg_3x3->Clear();//removes any entries in legend
	leg_3x3->AddEntry(hrec_unf_ratio[ki] ,"Unf."      ,"lp");
	leg_3x3->AddEntry(hrec_fold_ratio[ki],"Fold(Unf.)","lp");
	leg_3x3->DrawClone();
 
      }
      canv_3x3spectra ->Print(outPdfFile.c_str());
      canv_3x3genratio->Print(outPdfFile.c_str());
      canv_3x3recratio->Print(outPdfFile.c_str());
      
    }
    canvForPrint->Print(close_outPdfFile.c_str());      
    //assert(false);


    //---------------          
    //canvForPrint->Print(close_outPdfFile.c_str());      
     
    //---------------          
    if(drawRespMatrix){    
      //assert(false);
      drawRespMatrixFile((TH2D*)hmat, (TH2D*)hmat_rebin, (TH2D*)hmat_errors,
			 //(TH2D*)covmat_TH2, (TH2D*)PearsonBayes, (TH2D*)unfmat_TH2,
			 (double*)boundaries_pt_reco_mat, (int)nbins_pt_reco_mat,
			 (double*)boundaries_pt_gen_mat , (int)nbins_pt_gen_mat,
			 (std::string)outRespMatPdfFile , (bool)useSimpBins, (TFile*)fin_ppMC);    }//end draw response matrix
    


  std::cout<<std::endl<<"draw pdfs for Bayesian Unfolding done!"<<std::endl;


  }// end draw pdfs

  std::cout<<"writing output hists to file... "<<std::endl;

  fout->cd();
  
    
  
  // input data ------------------
  hrec_rebin->SetTitle("Data Meas.");hrec_rebin->Write("OS_MC_meas");
  //  if((bool)hJetQA_jtptEntries){ 
  //    hJetQA_jtptEntries->SetTitle("Data N_{Jets}");hJetQA_jtptEntries->Write("OS_MC_njets");}
  
  // input SS MC ---------------------
  hgen_rebin->SetTitle("SS MC Truth");  hgen_rebin->Write("SS_MC_truth");
  hrec_sameside_rebin->SetTitle("SS MC Meas.");  hrec_sameside_rebin->Write("SS_MC_meas");
  
  hmat->SetTitle("SS MC Response Matrix");  hmat->Write("SS_MC_mat");
  hmat_rebin->SetTitle("SS MC Response Matrix Rebinned");  hmat_rebin->Write("SS_MC_mat_rebin"); 
  hmat_errors->SetTitle("SS MC Response Matrix Errors");  hmat_errors->Write("SS_MC_mat_rebin_errors");
  if((bool)hmat_percenterrs){
    hmat_percenterrs->SetTitle("SS MC Response Matrix Percent Errors");hmat_percenterrs->Write("SS_MC_mat_rebin_percerrors");  }
  
  // input thy ----------------  
  CT10nlo  ->SetTitle("CT10 NLO Spectra");         CT10nlo  ->Write("NLO_CT10_NLO_R04_jtpt");	      
  CT14nlo  ->SetTitle("CT14 NLO Spectra");	   CT14nlo  ->Write("NLO_CT14_NLO_R04_jtpt");	      
  HERAPDF  ->SetTitle("HERAPDF NLO Spectra");	   HERAPDF  ->Write("NLO_HERAPDF105_NLO_R04_jtpt");  
  MMHTnlo  ->SetTitle("MMHT NLO Spectra");	   MMHTnlo  ->Write("NLO_MMHT2014_NLO_R04_jtpt");    
  NNPDFnnlo->SetTitle("NNPDF NNLO Spectra");	   NNPDFnnlo->Write("NLO_NNPDF_NLO_R04_jtpt");       
  
  // output hists -------------
  hunf->SetTitle("OS MC Unf."); hunf->Write("OS_MC_unf");    
  hunf_ss->SetTitle("SS MC Unf."); hunf_ss->Write("SS_MC_unf");    
  hfold->SetTitle("OS MC Fold(Unf.)"); hfold->Write("OS_MC_fold");        
  hfold_ss->SetTitle("SS MC Fold(Unf.)"); hfold_ss->Write("SS_MC_fold");        
  hrec_rebin_fakecorr->SetTitle("OS MC Fake Corr. Meas.");hrec_rebin_fakecorr->Write("OS_MC_measfakcorr");
  
  hfold_truth->SetTitle("SS MC Fold(Truth)");hfold_truth->Write("SS_MC_truth_fold");
  hfak->SetTitle("SS MC Meas. Fakes");hfak->Write("SS_MC_meas_fakes");
  hrec_sameside_rebin_fakecorr->SetTitle("SS MC Fake Corr. Meas."); hrec_sameside_rebin_fakecorr->Write("SS_MC_measfakcorr");
  
  covmat_TH2->Write("covmat");
  covmatabsval_TH2->Write("covmatabsval");
  PearsonBayes->Write("pearson");
  unfmat_TH2->Write("unfmat");

  // output ratio comparisons -------------
  // gen ratios (denom=mc truth)
  h_genratio_oppunf ->SetTitle("OS MC Unf./SS MC Truth");h_genratio_oppunf ->Write("ratio_OS_MC_unf_SS_MC_truth");  //data unf/mc truth	   
  h_genratio_oppfold->SetTitle("OS MC Fold(Unf.)/SS MC Truth");h_genratio_oppfold->Write("ratio_OS_MC_fold_SS_MC_truth"); //data fold(unf)/mc truth 
  h_genratio_oppmeas->SetTitle("OS MC Meas./SS MC Truth");h_genratio_oppmeas->Write("ratio_OS_MC_meas_SS_MC_truth"); //data meas/mc truth	   
  h_genratio_ssmeas ->SetTitle("SS MC Meas./SS MC Truth");h_genratio_ssmeas ->Write("ratio_SS_MC_meas_SS_MC_truth");  //mc meas/mc truth        
  
  // rec ratios (denom=data meas)
  h_recratio_oppunf ->SetTitle("OS MC Unf./OS MC Meas.");h_recratio_oppunf ->Write("ratio_OS_MC_unf_OS_MC_meas");  //data unf       / data meas	   
  h_recratio_oppfold->SetTitle("OS MC Fold(Unf.)/OS MC Meas.");h_recratio_oppfold->Write("ratio_OS_MC_fold_OS_MC_meas"); //data fold(unf) / data meas     
  h_recratio_ssmeas ->SetTitle("SS MC Meas./OS MC Meas.");h_recratio_ssmeas ->Write("ratio_SS_MC_meas_OS_MC_meas");  //mc meas        / data meas	   
  h_recratio_ssgen  ->SetTitle("SS MC Truth/OS MC Meas.");h_recratio_ssgen  ->Write("ratio_SS_MC_truth_OS_MC_meas");   //mc truth       / data meas        
  
  // fold ratio test
  h_foldratio_datafold->SetTitle("OS MC Fold(Unf.)/OS MC Fake Corr. Meas.");h_foldratio_datafold->Write("ratio_OS_MC_fold_OS_MC_measfakcorr");
  h_foldratio_mcfold  ->SetTitle("SS MC Fold(Truth)/SS MC Fake Corr. Meas.");   h_foldratio_mcfold  ->Write("ratio_SS_MC_fold_OS_MC_measfakcorr");
  
  // thy ratios w/ unfolded data
  h_thyratio_CT10nlo  ->Write("ratio_CT10_NLO_OS_MC_unf");
  h_thyratio_CT14nlo  ->Write("ratio_CT14_NLO_OS_MC_unf");
  h_thyratio_HERAPDF  ->Write("ratio_HERAPDF_NLO_OS_MC_unf");
  h_thyratio_MMHTnlo  ->Write("ratio_MMHTnlo_NLO_OS_MC_unf");    
  h_thyratio_NNPDFnnlo->Write("ratio_NNPDFnnlo_NNLO_OS_MC_unf");
  h_thyratio_mctruth->Write("ratio_NNPDFnnlo_ToyMCTruth_OS_MC_unf");


  if(dokIterQA){
    
    for(int ki=0; ki<nKiterMax;ki++) hunf_bayes[ki]->Write();
    for(int ki=0; ki<nKiterMax;ki++) hfold_bayes[ki]->Write();
    for(int ki=0; ki<nKiterMax;ki++) hgen_unf_ratio[ki]->Write();
    for(int ki=0; ki<nKiterMax;ki++) hgen_fold_ratio[ki]->Write();
    for(int ki=0; ki<nKiterMax;ki++) hrec_unf_ratio[ki]->Write();
    for(int ki=0; ki<nKiterMax;ki++) hrec_fold_ratio[ki]->Write();
      
    
  }

  // canvases ----------------
  if(drawPDFs){
    canv_spectra          ->SetTitle("I/O Spectra Canvas");        canv_spectra           ->Write("canv_spectra");
    canv_mc_fakes_spectra ->SetTitle("SS MC Fakes Spectra Canvas");   canv_mc_fakes_spectra  ->Write("canv_mc_fakes_spectra");
    canv_thy_spectra_1    ->SetTitle("NLO Thy Spectra 1 Canvas");  canv_thy_spectra_1     ->Write("canv_thy_spectra_1");
    canv_thy_spectra_2    ->SetTitle("NLO Thy Spectra 2 Canvas");  canv_thy_spectra_2     ->Write("canv_thy_spectra_2");                                
    
    canv_gen_ratio        ->SetTitle("SS MC Truth Ratios Canvas");   canv_gen_ratio          ->Write("canv_gen_ratio");
    canv_rec_ratio        ->SetTitle("OS MC Meas Ratios Canvas");  canv_rec_ratio          ->Write("canv_meas_ratio");
    canv_fold_ratio       ->SetTitle("Fold Test Ratios Canvas");  canv_fold_ratio         ->Write("canv_fakcorr_meas_ratio");
    canv_thy_ratio        ->SetTitle("NLO Thy Ratios Canvas");    canv_thy_ratio          ->Write("canv_thy_ratio");                                

    canv_covmat           ->SetTitle("Covariance Matrix Canvas");           canv_covmat        ->Write("canv_covmat");
    canv_absval_covmat    ->SetTitle("Abs Val. Covariance Matrix Canvas");  canv_absval_covmat ->Write("canv_covmatabsval");
    canv_pearson          ->SetTitle("Pearson Matrix Canvas");      canv_pearson               ->Write("canv_pearson");
    //  assert(false);      
    canv_unfmat           ->SetTitle("Unfolding Matrix Canvas");    canv_unfmat                ->Write("canv_unfmat");
    canv_mat_rebin        ->SetTitle("SS MC Response Matrix Canvas");  canv_mat_rebin             ->Write("canv_mat_rebin");
    canv_mat_percerrs     ->SetTitle("SS MC Response Matrix % Errors Canvas");  canv_mat_percerrs ->Write("canv_mat_percerrors");

    if(dokIterQA){ canv_3x3spectra->SetTitle("3x3 Unf. Spectra kIter QA Canvas"); canv_3x3spectra ->Write("canv_3x3spectra");}
    if(dokIterQA){canv_3x3genratio->SetTitle("3x3 Gen Ratio kIter QA Canvas");    canv_3x3genratio->Write("canv_3x3genratio");}
    if(dokIterQA){canv_3x3recratio->SetTitle("3x3 Rec Ratio kIter QA Canvas");    canv_3x3recratio->Write("canv_3x3recratio");}

  }

  
  

  

  
  
  
  


  
  std::cout<<"closing input/output root files..."<<std::endl;
  fout->Close();
  fin_ppMC->Close();
  return 0;
  
} // end unfoldMCSpectra



//  steering ---------------------------------------------------------------------------------
int main(int argc, char* argv[]){  int rStatus = -1;
  
  if( argc!=7 ){
    std::cout<<"do ./bayesUnfoldMCSpectra.exe <targMCDir> <baseOutputName> <doJetID> <useSimpleBins> <kIterInput>"<<std::endl;
    std::cout<<"actually... just open the damn code and look"<<std::endl;
    
    return rStatus;  }
  
  rStatus=1; // runtime error
  
  rStatus=bayesUnfoldNLOMCSpectra(  (const std::string)argv[1], (const std::string)argv[2], (const std::string)argv[3], 
				    (bool)std::atoi(argv[4]) ,   (bool)std::atoi(argv[5])  , (const int)std::atoi(argv[6]) ); 
  
  std::cout<<std::endl<<"done!"<<std::endl<<" return status: "<<rStatus<<std::endl<<std::endl;
  return rStatus;
}
