#include "unfoldSpectra.h"

//other settings
const bool drawPDFs=true; 
const bool debugMode=false, debugWrite=false;
const bool drawRespMatrix=true;
const bool dokIterQA=true;
//const bool applyNPCorrs=true; 


// CODE --------------------------------------------------
int bayesUnfoldDataSpectra_wNLO( std::string inFile_Data_dir , std::string inFile_MC_dir   , std::string inFile_MC_name, std::string baseName , 
			    const bool applyNPCorrs=false,
			    const bool doJetID=true     , const bool useSimpBins=false, const int kIterInput=4 ){
  
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
  inFile_Data_dir=SCRATCH_BASE+inFile_Data_dir;
  //inFile_MC_dir  =SCRATCH_BASE+inFile_MC_dir;
  
  //std::size_t radPos=inFile_MC_dir.find("_ak")+3;  
  //  const std::string radius="R"+inFile_MC_dir.substr( radPos,1 );
  const std::string radiusint="4";
  const std::string radius="R"+radiusint;
  const std::string RandEtaRange="_"+radius+"_20_eta_20";
  const std::string RandEtaRange_plotTitle=" "+radius+" 20eta20";// R4 20eta20
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;
  if(debugMode)std::cout<<"RandEtaRange string is = "<<RandEtaRange<<std::endl;  
  
  //  std::size_t jetTypePos=radPos+1, jetsPos=inFile_MC_dir.find("Jets");
  //  const std::string jetType=inFile_MC_dir.substr( jetTypePos,(jetsPos-jetTypePos) );//PFJets
  const std::string jetType="PFJets";
  //  const std::string fullJetType="ak"+inFile_MC_dir.substr( radPos,1 )+jetType;//"ak4PFJets"
  const std::string fullJetType="ak"+radiusint+jetType;
  if(debugMode)std::cout<<"jetType string is = "<<jetType<<std::endl;
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  
  
  // INFILE NAME(S) -----------
  //const std::string inFile_MC_name="/Py8_CUETP8M1_QCDjetAllPtBins_"+fullJetType+"-allFiles.root";
  //const std::string inFile_Data_name="/HighPtJetTrig_"+fullJetType+"-allFiles.root";
  const std::string inFile_Data_name="/HighPtJetTrig_ak4PF-allFiles.root";
  
  // OUTPUT FILE, NAME(S) -----------
  std::string outFileName=unfoldDataSpectra_outdir+fullJetType;
  if(doJetID)outFileName+="_wjtID";
  if(!useSimpBins)outFileName+="_anabins";
  else outFileName+="_simpbins";
  outFileName+="_"+baseName;
  
  if(applyNPCorrs)
    outFileName+="_NLOMC_wNP";
  else
    outFileName+="_NLOMC";
  if(debugMode)std::cout<<"baseName="<<baseName<<std::endl;
  if(debugMode)std::cout<<"outFileName="<<outFileName<<std::endl;

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
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from NLO file"<<std::endl; 
  std::cout<<"input NLO dir : "<< (inFile_MC_dir)  <<std::endl; 
  std::cout<<"NLO file name : "<< (inFile_MC_name)<<std::endl;   std::cout<<std::endl<<std::endl;  
  TFile *fpp_MC = TFile::Open( (inFile_MC_dir+inFile_MC_name).c_str());
  
  // ---------- open MC "response" matrix 
  std::string TH2_title="response";//"_ynew_th2"
  if(applyNPCorrs)TH2_title+="_NPynew_th2";
  else TH2_title+="_ynew_th2";
  
  //get the response matrix made by readforests
  TH2D* hmat = (TH2D*)fpp_MC->Get( TH2_title.c_str() );
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
  if(debugWrite)hmat_errors->Write(errTH2_title.c_str());
  if(debugMode)hmat_errors->Print("base");
  
  if(clearOverUnderflows){
    errTH2_title+="_noOverUnderFlows";
    TH2clearOverUnderflows((TH2D*)hmat_errors);
    if(debugWrite)hmat_errors->Write(errTH2_title.c_str());
    if(debugMode)hmat_errors->Print("base");  }
  
  // give response matrix the correct errors
  setRespMatrixErrs( (TH2D*)hmat_rebin, (TH2D*) hmat_errors , (bool)zeroBins);  
  TH2_title+="_wseterrs";
  if(debugWrite)hmat_rebin->Write(TH2_title.c_str());
  if(debugMode)hmat_rebin->Print("base");
  
    
  // ---------- reco, measured same-side spectra used to create response matrix, and for "sameside" unfolding test
  std::string histTitle2="smeared_rnd";
  if(applyNPCorrs) histTitle2+="_NPynew";
  else histTitle2+="_ynew";
  
  TH1D*  hrec_sameside = (TH1D*)fpp_MC->Get( histTitle2.c_str() ); 
  hrec_sameside->Scale(1./NLOMCscaling);
  hrec_sameside->Scale(etaBinWidth);//cause of the way the thy hists are made
  multiplyBinWidth(hrec_sameside);
  if(debugWrite)hrec_sameside->Write(histTitle2.c_str());
  if(debugMode)hrec_sameside->Print("base");
  
  
  
  histTitle2+="_clone";
  TH1D *hrec_sameside_rebin = (TH1D*)hrec_sameside->Clone( (histTitle2).c_str() );
  if(debugWrite)hrec_sameside_rebin->Write( histTitle2.c_str() );
  if(debugMode)hrec_sameside_rebin->Print("base");  
  
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
  

  // ---------- gen, NLO truth spectra
  std::string genHistTitle="theory_rnd";//"hpp_mcclosure_gen"
  if(applyNPCorrs) genHistTitle+="_NPynew";
  else genHistTitle+="_ynew";
  
  TH1D* hgen = (TH1D*)fpp_MC->Get( genHistTitle.c_str() );
  hgen->Scale(1./NLOMCscaling);//05/09/18
  hgen->Scale(etaBinWidth);
  multiplyBinWidth(hgen);  
  if(debugWrite)hgen->Write(genHistTitle.c_str());
  if(debugMode)hgen->Print("base");    
  
  genHistTitle+="_clone";
  TH1D* hgen_rebin = (TH1D*)hgen->Clone( (genHistTitle).c_str() );
  if(debugWrite)hgen_rebin->Write(genHistTitle.c_str());
  if(debugMode)hgen_rebin->Print("base");
  
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
  

  // ppData input histos -------------------------
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from Data file"<<std::endl; 
  std::cout<<"input Data dir : "<< (inFile_Data_dir)  <<std::endl; 
  std::cout<<"Data file name : "<< (inFile_Data_name)<<std::endl;   std::cout<<std::endl;//<<std::endl;  
  TFile *fpp_Data = TFile::Open( (inFile_Data_dir+inFile_Data_name).c_str());
  
  //for output
  if(debugWrite)fout->cd();    
  
  TH1D* hJetQA_jtptEntries=(TH1D*)fpp_Data->Get("hJetQA_1wJetID_jtptEntries");
  if(hJetQA_jtptEntries){
    hJetQA_jtptEntries = (TH1D*)hJetQA_jtptEntries->Rebin(nbins_pt_reco, hJetQA_jtptEntries->GetTitle(), boundaries_pt_reco);  
  }
  
  // ---------- reco, measured spectra to unfold
  std::string histTitle="hJetQA";
  if(doJetID)histTitle+="_1wJetID";
  else histTitle+="_0wJetID";
  histTitle+="_jtpt";
  
  TH1D*  hrec = (TH1D*)fpp_Data->Get( histTitle.c_str() ); 
  if(debugWrite)hrec->Write();
  if(debugMode)hrec->Print("base");
  
  histTitle+="_divBylumietabin";
  double effIntgrtdLumi=computeEffLumi(fpp_Data);
  hrec->Scale(1./effIntgrtdLumi); // lumi
  if(debugWrite)hrec->Write( (histTitle).c_str() );
  if(debugMode)hrec->Print("base");
  
  histTitle+="_clone";
  TH1D *hrec_rebin = (TH1D*)hrec->Clone( (histTitle).c_str() );
  if(debugWrite)hrec_rebin->Write(histTitle.c_str());
  if(debugMode)hrec_rebin->Print("base");
  
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
  std::string NLOMCtitle_str="";//this is for the Toy NLO used specifically. the TH1s below are pure theory, no toy procedure on top.
  std::string CT10NPs ="" ; 
  std::string CT14NPs ="" ; 
  std::string HERANPs ="" ; 
  std::string MMHTNPs ="" ; 
  std::string NNPDFNPs="" ; 

  TH1D* CT10nlo  =(TH1D*)makeThyHist_00eta20_v2((fNLOFile_R04_CT10nlo  ).c_str(),applyNPCorrs, &CT10NPs);
  CT10nlo->SetMarkerSize(0);
  CT10nlo->SetLineColor(kBlack);  
  CT10nlo = (TH1D*)CT10nlo->Rebin(nbins_pt_gen,"pp_CT10Thy_rebin",boundaries_pt_gen);  

  TH1D* CT14nlo  =(TH1D*)makeThyHist_00eta20_v2((fNLOFile_R04_CT14nlo  ).c_str(),applyNPCorrs, &CT14NPs);
  CT14nlo->SetMarkerSize(0);
  CT14nlo->SetLineColor(kGreen);  
  CT14nlo=(TH1D*)CT14nlo->Rebin(nbins_pt_gen,"pp_CT14Thy_rebin",boundaries_pt_gen);

  TH1D* HERAPDF  =(TH1D*)makeThyHist_00eta20_v2((fNLOFile_R04_HERAPDF  ).c_str(),applyNPCorrs, &HERANPs);
  HERAPDF->SetMarkerSize(0);
  HERAPDF->SetLineColor(kViolet-5);  
  HERAPDF=(TH1D*)HERAPDF->Rebin(nbins_pt_gen,"pp_HERAPDF_rebin",boundaries_pt_gen);
  
  TH1D* MMHTnlo  =(TH1D*)makeThyHist_00eta20_v2((fNLOFile_R04_MMHTnlo  ).c_str(),applyNPCorrs, &MMHTNPs);
  MMHTnlo->SetMarkerSize(0);
  MMHTnlo->SetLineColor(kOrange+7);  
  MMHTnlo=(TH1D*)MMHTnlo->Rebin(nbins_pt_gen,"pp_MMHT_rebin",boundaries_pt_gen);
  
  TH1D* NNPDFnnlo=(TH1D*)makeThyHist_00eta20_v2((fNLOFile_R04_NNPDFnnlo).c_str(),applyNPCorrs, &NNPDFNPs);
  NNPDFnnlo->SetMarkerSize(0);
  NNPDFnnlo->SetLineColor(kCyan-6); 
  NNPDFnnlo=(TH1D*)NNPDFnnlo->Rebin(nbins_pt_gen,"pp_NNPDFnlo_rebin",boundaries_pt_gen);
  
  if(     inFile_MC_name.find("CT10")!=std::string::npos)   NLOMCtitle_str="CT10 NLO"      + CT10NPs ;
  else if(inFile_MC_name.find("CT14")!=std::string::npos)   NLOMCtitle_str="CT14 NLO"      + CT14NPs ;
  else if(inFile_MC_name.find("HERAPDF")!=std::string::npos)NLOMCtitle_str="HERAPDF15 NLO" + HERANPs ;
  else if(inFile_MC_name.find("MMHT")!=std::string::npos)   NLOMCtitle_str="MMHT2014 NLO"  + MMHTNPs ;
  else if(inFile_MC_name.find("NNPDF")!=std::string::npos)  NLOMCtitle_str="NNPDF30 NNLO"  + NNPDFNPs;

  // ----- unfolding setup+use below here

  // RooUnfold Response Matrix 
  std::cout<<"calling RooUnfoldResponse "<<std::endl;
  RooUnfoldResponse roo_resp( hrec_sameside_resp_rebin, hgen_resp_rebin, hmat_rebin, ("Response_matrix"+RandEtaRange).c_str()) ;
  roo_resp.UseOverflow(doOverUnderflows);    
  if(debugWrite)roo_resp.Write();
  
  TH1D* hfak=  (TH1D*) roo_resp.Hfakes() ;
  hfak->Print("base");

  //cosmetics
  hfak->SetMarkerColor(kGreen);
  hfak->SetLineColor(kGreen);
  
  
  // Bayesian unfolding -------------------------   
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
  hunf->SetName("ppData_BayesUnf_Spectra");
  hunf->SetTitle( ("Unf. Data, kIter="+std::to_string(kIterInput)).c_str());
  if(debugMode)hunf->Print("base");

  std::string descString="";
  getUnfDetails(hunf,&descString);
  std::string methodString=", Bayesian ";      

  //cosmetics
  hunf->SetMarkerStyle(kOpenCircle);
  hunf->SetMarkerColor(kRed);
  hunf->SetLineColor(kRed);
  hunf->SetMarkerSize(1.02);     
  
  std::cout<<"folding unfolded data histogram!!"<<std::endl;
  TH1D* hfold=(TH1D*)roo_resp.ApplyToTruth(hunf);          
  hfold->SetName("ppData_BayesFold_Spectra");
  hfold->SetTitle(("Fold. Data, kIter="+std::to_string(kIterInput)).c_str());
  if(debugMode)hfold->Print("base");
  
  //cosmetics
  hfold->SetMarkerStyle(kOpenCircle);
  hfold->SetMarkerColor(kGreen-5);
  hfold->SetLineColor(  kGreen-5);
  hfold->SetMarkerSize(1.02);     
  
  std::cout<<"folding NLO Truth histogram!!"<<std::endl;
  TH1D* hfold_truth=(TH1D*)roo_resp.ApplyToTruth(hgen_rebin);          
  hfold_truth->SetName("ppMC_BayesFoldMCTruth_Spectra");
  hfold_truth->SetTitle(("Fold. Toy NLO Truth, kIter="+std::to_string(kIterInput)).c_str());
  if(debugMode)hfold_truth->Print("base");
  
  //cosmetics
  hfold_truth->SetMarkerStyle(kOpenStar);
  hfold_truth->SetMarkerColor(kGreen+3);
  hfold_truth->SetLineColor(  kGreen+3);
  hfold_truth->SetMarkerSize(1.02);     
  
  //i should probably do this differently, but for now ignore
//  std::cout<<"folding thy histogram!!"<<std::endl;
//  TH1D* CT10nlo_rebin=(TH1D*)CT10nlo->Rebin(nbins_pt_gen,"MC_CT10nlo_fold_rebin",boundaries_pt_gen);
//  CT10nlo_rebin->Scale(etaBinWidth);
//  multiplyBinWidth(CT10nlo_rebin);  
//  TH1clearOverUnderflows((TH1*)CT10nlo_rebin);
//  if(debugMode) CT10nlo_rebin->Print("base");
//  CT10nlo_rebin->Print("base");
//  //assert(false);
//
//  TH1D* hfold_CT10nlo=(TH1D*)roo_resp.ApplyToTruth(CT10nlo_rebin);          
//  hfold_CT10nlo->SetName("ppMC_BayesFoldCT10nlo_Spectra");
//  hfold_CT10nlo->SetTitle("Fold. CT10nlo");
//  if(debugMode)hfold_CT10nlo->Print("base");

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
  




  // --------- RATIOS WITH TOY NLO TRUTH ----------------
  TH1D *hgen_rebin_ratiobin=(TH1D*)hgen_rebin->Clone("ppData_Gen_Ratio_denom");
  if(debugMode)hgen_rebin_ratiobin->Print("base");
  
  
  TH1D *h_genratio_oppunf = (TH1D*)hunf->Clone( "ppData_Gen_Ratio_OppUnf" );
  h_genratio_oppunf->SetTitle( "Unf. Data/Toy NLO Truth" );
  h_genratio_oppunf->Divide(hgen_rebin_ratiobin);
  if(debugMode)h_genratio_oppunf->Print("base");
  
  TH1D *h_genratio_oppfold = (TH1D*)hfold->Clone( "ppData_Gen_Ratio_OppFold" );
  h_genratio_oppfold->SetTitle( "Fold. Data/Toy NLO Truth" );
  h_genratio_oppfold->Divide(hgen_rebin_ratiobin);
  if(debugMode)h_genratio_oppfold->Print("base");
  
  TH1D *h_genratio_oppmeas = (TH1D*)hrec_rebin->Clone( "ppData_Gen_Ratio_Meas" );
  h_genratio_oppmeas->SetTitle( "Meas. Data/Toy NLO Truth" );
  h_genratio_oppmeas->Divide(hgen_rebin_ratiobin);
  if(debugMode)h_genratio_oppmeas->Print("base");
  
  TH1D *h_genratio_ssmeas = (TH1D*)hrec_sameside_rebin->Clone( "ppMC_Gen_Ratio4_SSMeas" );
  h_genratio_ssmeas->SetTitle( "JER Smeared NLO/Toy NLO Truth" );
  h_genratio_ssmeas->Divide(hgen_rebin_ratiobin);
  if(debugMode)h_genratio_ssmeas->Print("base");
  
  

  // ---------------- RATIOS WITH MEAS -----------------
  TH1D *h_recratio_oppunf = (TH1D*)hunf->Clone( "ppData_Meas_Ratio_OppUnf" );
  h_recratio_oppunf->SetTitle( "Unf. Data/RECO Data" );
  h_recratio_oppunf->Divide(hrec_rebin);
  if(debugMode)h_recratio_oppunf->Print("base");
  
  TH1D *h_recratio_ssmeas = (TH1D*)hrec_sameside_rebin->Clone( "ppMC_Meas_Ratio_SSMeas" );
  h_recratio_ssmeas->SetTitle( "JER Smeared NLO/RECO Data" );
  h_recratio_ssmeas->Divide(hrec_rebin);
  if(debugMode)h_recratio_ssmeas->Print("base");  
  
  TH1D *h_recratio_ssgen = (TH1D*)hgen_rebin->Clone( "ppMC_Meas_Ratio_SSTruth" );    
  h_recratio_ssgen->SetTitle( "Toy NLO Truth/RECO Data" );
  h_recratio_ssgen->Divide(hrec_rebin);
  if(debugMode)h_recratio_ssgen->Print("base");

  // ---------------- FAKE CORRECTION -----------------
  //Double_t fac=(hrec_rebin->Integral()/hrec_sameside_rebin->Integral());
  
  TH1D *hrec_rebin_fakecorr=(TH1D*)hrec_rebin->Clone("ppData_Meas_Ratio_denom");
  if(debugMode)hrec_rebin_fakecorr->Print("base");  
  
  TH1D *hrec_sameside_rebin_fakecorr=(TH1D*)hrec_sameside_rebin->Clone("ppMC_Meas_Ratio_denom");
  if(debugMode)hrec_sameside_rebin_fakecorr->Print("base");  

  TH1D *hfold_fakecorr=(TH1D*)hfold->Clone("ppData_Folded_wFakes");
  if(debugMode)hrec_rebin_fakecorr->Print("base");  
  
  TH1D *hfold_truth_fakecorr=(TH1D*)hfold_truth->Clone("ppMC_TruthFold_wFakes");
  if(debugMode)hrec_sameside_rebin_fakecorr->Print("base");
  
  
  applyFakeCorrection(hrec_rebin_fakecorr ,  hrec_sameside_rebin_fakecorr,
		      hfold_fakecorr      ,  hfold_truth_fakecorr        ,
		      hfak      );
  

  // ---------------- FOLDED RATIOS (TOY NLO AND DATA) ----------------- //
  
  TH1D *h_recratio_oppfold = (TH1D*)hfold->Clone( "ppData_MeasmFakes_Ratio_DataFold" );
  h_recratio_oppfold->SetTitle( "Data Fold(Unf.)/(RECO Data - Fakes)" );
  h_recratio_oppfold->Divide(hrec_rebin_fakecorr);
  if(debugMode)h_recratio_oppfold->Print("base");  
  
  TH1D *h_recratio_truthfold = (TH1D*)hfold_truth->Clone( "ppData_MeasmFakes_Ratio_TruthFold" );
  h_recratio_truthfold->SetTitle( "Fold Toy NLO Truth/(JER Smeared NLO - Fakes)" );
  h_recratio_truthfold->Divide(hrec_sameside_rebin_fakecorr);
  if(debugMode)h_recratio_truthfold->Print("base");  
  
  TH1D *h_foldratio_datafold=(TH1D*)hfold_fakecorr->Clone("ppData_Meas_Ratio_DataFoldpFakes");
  h_foldratio_datafold->SetTitle( "(Data Fold(Unf.) + Fakes)/RECO Data" );
  h_foldratio_datafold->Divide(hrec_rebin);
  if(debugMode)h_foldratio_datafold->Print("base");  
  
  TH1D *h_foldratio_mcfold=(TH1D*)hfold_truth_fakecorr->Clone("ppData_Meas_Ratio_TruthFoldpFakes");
  h_foldratio_mcfold->SetTitle( "(Fold Toy NLO Truth + Fakes)/JER Smeared NLO" );
  h_foldratio_mcfold->Divide(hrec_sameside_rebin);
  if(debugMode)h_foldratio_mcfold->Print("base");



  
  TH1D* hunf_bayes[nKiterMax]={};
  TH1D* hfold_bayes[nKiterMax]={};
  TH1D* hgen_unf_ratio[nKiterMax]={};
  TH1D* hrec_unf_ratio[nKiterMax]={};
  
  TH1D* hgen_fold_ratio[nKiterMax]={};
  TH1D* hrec_fold_ratio[nKiterMax]={};
  TH2D* hcovmat_bayes[nKiterMax]={};
  TH2D* hcovmatabsval_bayes[nKiterMax]={};
  TH2D* hpearson_bayes[nKiterMax]={};
  
  int kIter_start=kIterInput-kIterRange;
  int kIter_end=kIterInput+kIterRange;
  
  if(dokIterQA){
    if(debugMode)std::cout<<"unfolding across diff kIter values"<<std::endl;
    if(debugMode)std::cout<<kIter_start<<" <= kIter <= "<<kIter_end<<std::endl;
    for(int ki=0; ki<nKiterMax;ki++){
      int current_kIter=kIter_start+ki;
      std::cout<<"ki="<<ki<<", current_kIter="<<current_kIter<<std::endl;      
      RooUnfoldBayes unf_bayes_kIterQA(&roo_resp, hrec_rebin, current_kIter);        
      std::cout<<"unf Chi2="<<unf_bayes_kIterQA.RooUnfold::Chi2((const TH1*)hgen_rebin)<<std::endl;
      
      hunf_bayes[ki]=(TH1D*)unf_bayes_kIterQA.Hreco(errorTreatment);
      hunf_bayes[ki]->SetName(("Data_unf_kIter"+std::to_string(current_kIter)).c_str());
      hunf_bayes[ki]->SetTitle(("Data Unf., kIter="+std::to_string(current_kIter)+";Jet p_{T}; #sigma" ).c_str());
      //cosmetics, keep here so they propagate to the ratio hists
      hunf_bayes[ki]->SetMarkerStyle(kOpenCircle);
      hunf_bayes[ki]->SetMarkerColor(kRed);
      hunf_bayes[ki]->SetLineColor(kRed);
      hunf_bayes[ki]->SetMarkerSize(1.02);     

      if(debugMode)hunf_bayes[ki]->Print("base");
      
      hfold_bayes[ki]=(TH1D*)roo_resp.ApplyToTruth(hunf_bayes[ki]);
      hfold_bayes[ki]->SetName(("Data_fold_kIter"+std::to_string(current_kIter)).c_str());
      hfold_bayes[ki]->SetTitle(("Data Fold(Unf.), kIter="+std::to_string(current_kIter)+"; Jet p_{T}; #sigma").c_str());
      if(debugMode)hfold_bayes[ki]->Print("base");
      //cosmetics
      hfold_bayes[ki]->SetMarkerStyle(kOpenCircle);
      hfold_bayes[ki]->SetMarkerColor(kGreen-5);
      hfold_bayes[ki]->SetLineColor(  kGreen-5);
      hfold_bayes[ki]->SetMarkerSize(1.02);     
      
      if(debugMode)hfold_bayes[ki]->Print("base");
      
      hgen_unf_ratio[ki]=(TH1D*)hunf_bayes[ki]->Clone(("Data_unf_genratio_kIter"+std::to_string(current_kIter)).c_str());
      hgen_unf_ratio[ki]->SetTitle(("Data Unf./Toy NLO Truth Truth, kIter="+std::to_string(current_kIter)+"; Jet p_{T}; ratio").c_str());
      hgen_unf_ratio[ki]->Divide(hgen_rebin);
      
      hrec_unf_ratio[ki]=(TH1D*)hunf_bayes[ki]->Clone(("Data_unf_recratio_kIter"+std::to_string(current_kIter)).c_str());
      hrec_unf_ratio[ki]->SetTitle(("Data Unf./Data Meas., kIter="+std::to_string(current_kIter)+"; Jet p_{T}; ratio").c_str());
      hrec_unf_ratio[ki]->Divide(hrec_rebin);
      //hrec_unf_ratio[ki]->SetTitle(("Data Unf./Data Fake Corr. Meas., kIter="+std::to_string(current_kIter)+"; Jet p_{T}; ratio").c_str());
      //hrec_unf_ratio[ki]->Divide(hrec_rebin_fakecorr);
      
      hgen_fold_ratio[ki]=(TH1D*)hfold_bayes[ki]->Clone(("Data_fold_genratio_kIter"+std::to_string(current_kIter)).c_str());
      hgen_fold_ratio[ki]->SetTitle(("Data Fold./Toy NLO Truth Truth, kIter="+std::to_string(current_kIter)+"; Jet p_{T}; ratio").c_str());
      hgen_fold_ratio[ki]->Divide(hgen_rebin);
      
      hrec_fold_ratio[ki]=(TH1D*)hfold_bayes[ki]->Clone(("Data_fold_recratio_kIter"+std::to_string(current_kIter)).c_str());
      hrec_fold_ratio[ki]->SetTitle(("Data Fold./Data Fake Corr. Meas., kIter="+std::to_string(current_kIter)+"; Jet p_{T}; ratio").c_str());
      hrec_fold_ratio[ki]->Divide(hrec_rebin_fakecorr);
      

      //do bin width normalizations here, so the ratios in this loop arent affected by hgen/hrec hists not having been bin-width normalized yet
      hunf_bayes[ki]->Scale(1./etaBinWidth);
      divideBinWidth(hunf_bayes[ki]);
      
      hfold_bayes[ki]->Scale(1./etaBinWidth);
      divideBinWidth(hfold_bayes[ki]);


      if(debugMode)std::cout<<"getting covariance matrix, kIter="<<current_kIter<<std::endl;
      TMatrixD covmat = unf_bayes_kIterQA.Ereco(errorTreatment);
      hcovmat_bayes[ki]=new TH2D(covmat);
      hcovmat_bayes[ki]->SetName( ("Bayes_covmat_"+std::to_string(current_kIter)).c_str() );      
      hcovmat_bayes[ki]->SetTitle( ("Bayes Covariance Matrix, kIter="+std::to_string(current_kIter)).c_str()  );
      hcovmat_bayes[ki]->GetXaxis()->SetTitle("RECO Jet p_{T} Bin #");
      hcovmat_bayes[ki]->GetYaxis()->SetTitle("GEN Jet p_{T} Bin #");
      hcovmat_bayes[ki]->GetZaxis()->SetLabelSize(0.035);
       
      if(debugMode)std::cout<<"calling absval_th2content"<<std::endl;
      hcovmatabsval_bayes[ki]=(TH2D*)absVal_TH2Content((TH2D*)hcovmat_bayes[ki]);
      hcovmatabsval_bayes[ki]->SetName( ("Bayes_covmatabsval_"+std::to_string(current_kIter)).c_str() );      
      hcovmatabsval_bayes[ki]->SetTitle( ("Bayes |Covariance| Matrix, kIter="+std::to_string(current_kIter)).c_str()  );
      hcovmatabsval_bayes[ki]->GetXaxis()->SetTitle("RECO Jet p_{T} Bin #");
      hcovmatabsval_bayes[ki]->GetYaxis()->SetTitle("GEN Jet p_{T} Bin #");
      hcovmatabsval_bayes[ki]->GetZaxis()->SetLabelSize(0.035);
      
      if(debugMode)std::cout<<"calling CalculatePearsonCoefficients... getting pearson matrix"<<std::endl;
      hpearson_bayes[ki] = (TH2D*) CalculatePearsonCoefficients(&covmat, false, ("Bayes_pearson_"+std::to_string(current_kIter)) );   
      hpearson_bayes[ki]->SetName( ("Bayes_pearson_"+std::to_string(current_kIter)).c_str() );      
      hpearson_bayes[ki]->SetTitle( ("Bayes Pearson Matrix, kIter="+std::to_string(current_kIter)).c_str()  );
      hpearson_bayes[ki]->GetXaxis()->SetTitle("RECO Jet p_{T} Bin #");
      hpearson_bayes[ki]->GetYaxis()->SetTitle("GEN Jet p_{T} Bin #");
      hpearson_bayes[ki]->GetZaxis()->SetLabelSize(0.035);
      hpearson_bayes[ki]->SetAxisRange(-1., 1., "Z");  
      
    }//end kIterLoop
  }//end if(dokIterQA)
  //assert(false);

  // ----- PUT BIN WIDTH(s) DIVISIONS + NORMALIZING HERE ----- //  
  // -- JER Smeared NLO -- //  
  histTitle2+="_normbinwidth";
  hrec_sameside_rebin->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hrec_sameside_rebin); 
  if(debugWrite)hrec_sameside_rebin->Write( histTitle2.c_str() );
  if(debugMode)hrec_sameside_rebin->Print("base");  

  // -- JER Smeared NLO -- //  
  histTitle2+="_fakecorr";
  hrec_sameside_rebin_fakecorr->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hrec_sameside_rebin_fakecorr); 
  if(debugWrite)hrec_sameside_rebin_fakecorr->Write( histTitle2.c_str() );
  if(debugMode)hrec_sameside_rebin_fakecorr->Print("base");  
  
  // -- Toy NLO Truth -- //  
  genHistTitle+="_normbinwidth";
  hgen_rebin->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hgen_rebin);
  if(debugWrite)hgen_rebin->Write(genHistTitle.c_str());
  if(debugMode)hgen_rebin->Print("base");  

  // -- Toy NLO Truth RESP -- //  
  hgen_resp_rebin->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hgen_resp_rebin);
  if(debugWrite)hgen_resp_rebin->Write(genHistTitle.c_str());
  if(debugMode)hgen_resp_rebin->Print("base");  

  // -- Toy NLO RESP FAKES -- //  
  hfak->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hfak);
  if(debugWrite)hfak->Write();
  if(debugMode)hfak->Print("base");  
  
  // -- JER Smeared NLO RESP -- //  
  hrec_sameside_resp_rebin->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hrec_sameside_resp_rebin);
  if(debugWrite)hrec_sameside_resp_rebin->Write();
  if(debugMode)hrec_sameside_resp_rebin->Print("base");  
  
  // -- FOLD(Toy NLO TRUTH) -- //
  hfold_truth->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hfold_truth); 
  if(debugWrite)hfold_truth->Write();
  if(debugMode)hfold_truth->Print("base");  

  // -- FOLD(Toy NLO TRUTH) + FAKES -- //
  hfold_truth_fakecorr->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hfold_truth_fakecorr); 
  if(debugWrite)hfold_truth_fakecorr->Write();
  if(debugMode)hfold_truth_fakecorr->Print("base");    
  
  // -- DATA RECO -- //  
  histTitle+="_normbinwidth";
  hrec_rebin->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hrec_rebin); 
  if(debugWrite)hrec_rebin->Write(histTitle.c_str());
  if(debugMode)hrec_rebin->Print("base");  

  // -- DATA RECO + FAKES -- //  
  histTitle+="_fakecorr";
  hrec_rebin_fakecorr->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hrec_rebin_fakecorr); 
  if(debugWrite)hrec_rebin_fakecorr->Write(histTitle.c_str());
  if(debugMode)hrec_rebin_fakecorr->Print("base");  
  
  // -- DATA UNF -- //
  hunf->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hunf); 
  if(debugWrite)hunf->Write();
  if(debugMode)hunf->Print("base");  
  
  // -- DATA FOLD(UNF) -- //
  hfold->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hfold); 
  if(debugWrite)hfold->Write();
  if(debugMode)hfold->Print("base");  
  // ----- END BIN WIDTH DIVISIONS + NORMALIZING ----- //
  
  // -- DATA FOLD(UNF) + FAKES-- //
  hfold_fakecorr->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hfold_fakecorr); 
  if(debugWrite)hfold_fakecorr->Write();
  if(debugMode) hfold_fakecorr->Print("base");  

  // ---------------- THY RATIOS (THY / DATA) ----------------- //
  // must be after the binwidth divisions + normalization (the thy hists are made this way by default)
  
  TH1D* h_thyratio_CT10nlo  =(TH1D*)CT10nlo  ->Clone("");
  //  h_thyratio_CT10nlo = (TH1D*)h_thyratio_CT10nlo->Rebin(nbins_pt_gen,"pp_CT10Thy_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_CT10nlo  ->SetTitle("CT10 NLO/Data Unf.");  
  h_thyratio_CT10nlo  ->Divide(hunf);
  if(debugMode)h_thyratio_CT10nlo  ->Print("base");
  
  TH1D* h_thyratio_CT14nlo  =(TH1D*)CT14nlo  ->Clone("");  
  //  h_thyratio_CT14nlo=(TH1D*)h_thyratio_CT14nlo->Rebin(nbins_pt_gen,"pp_CT14Thy_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_CT14nlo  ->SetTitle("CT14 NLO/Data Unf.");
  h_thyratio_CT14nlo  ->Divide(hunf);
  if(debugMode)h_thyratio_CT14nlo  ->Print("base");
  
  TH1D* h_thyratio_HERAPDF  =(TH1D*)HERAPDF  ->Clone("");  
  //  h_thyratio_HERAPDF=(TH1D*)h_thyratio_HERAPDF->Rebin(nbins_pt_gen,"pp_HERAPDF_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_HERAPDF  ->SetTitle("HERAPDF15 NLO/Data Unf.");
  h_thyratio_HERAPDF  ->Divide(hunf);
  if(debugMode)h_thyratio_HERAPDF  ->Print("base");
  
  TH1D* h_thyratio_MMHTnlo  =(TH1D*)MMHTnlo  ->Clone("");
  //  h_thyratio_MMHTnlo=(TH1D*)h_thyratio_MMHTnlo->Rebin(nbins_pt_gen,"pp_MMHT_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_MMHTnlo  ->SetTitle("MMHT14 NLO/Data Unf.");
  h_thyratio_MMHTnlo  ->Divide(hunf);
  if(debugMode)h_thyratio_MMHTnlo  ->Print("base");
  
  TH1D* h_thyratio_NNPDFnnlo=(TH1D*)NNPDFnnlo->Clone("");
  //  h_thyratio_NNPDFnnlo=(TH1D*)h_thyratio_NNPDFnnlo->Rebin(nbins_pt_gen,"pp_NNPDFnlo_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_NNPDFnnlo->SetTitle("NNPDF NLO/Data Unf.");
  h_thyratio_NNPDFnnlo->Divide(hunf);
  if(debugMode)h_thyratio_NNPDFnnlo->Print("base");
  
  TH1D* h_thyratio_mctruth=(TH1D*)hgen_rebin->Clone("");
  //  h_thyratio_mctruth=(TH1D*)h_thyratio_mctruth->Rebin(nbins_pt_gen,"pp_MCTruth_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_mctruth->SetTitle("Toy NLO Truth/Data Unf.");
  h_thyratio_mctruth->Divide(hunf);
  if(debugMode)h_thyratio_NNPDFnnlo->Print("base");
  
  

  
  
  
    

  
  
  //  drawPDFS -------------------------------------------------- 
  //TCanvas pointer for writing the pdf file
  TCanvas* canvForPrint=NULL;
  std::string outPdfFile=outBayesPdfFile;
  std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";
  
  //TCanvas pointers for writing canvs to file
  TCanvas *canv_spectra=NULL, *canv_mc_fakes_spectra=NULL, *canv_thy_spectra_1=NULL, *canv_thy_spectra_2=NULL;
  TCanvas *canv_gen_ratio=NULL, *canv_rec_ratio=NULL, *canv_fold_ratio=NULL, *canv_fold_ratio2=NULL, *canv_thy_ratio=NULL,*canv_thy_ratio2=NULL; 
  TCanvas *canv_covmat=NULL, *canv_absval_covmat=NULL, *canv_pearson=NULL, *canv_unfmat=NULL, *canv_mat_rebin=NULL, *canv_mat_percerrs=NULL;
  TCanvas *canv_3x3spectra=NULL, *canv_3x3genratio=NULL, *canv_3x3recratio=NULL, *canv_kIterRatio;
  TCanvas *canv_3x3covmat=NULL, *canv_3x3covmatabsval=NULL, *canv_3x3pearson=NULL;
  
  TH2D* hmat_percenterrs=NULL;
  TH1D* hkiter_ratio=NULL;
  if(!drawPDFs)std::cout<<std::endl<<"NOT drawing PDFs for Bayesian Unfolding!"<<std::endl<<std::endl;
  else{ std::cout<<std::endl<<"drawing output PDFs for Bayesian Unfolding..."<<std::endl;
    
    // canvas for printing
    //TCanvas* canvForPrint=new TCanvas("printCanvas","for printing",           1400,1000);    
    canvForPrint=new TCanvas("printCanvas","for printing",           1400,1000);    
    canvForPrint->cd();
    canvForPrint->Print(open_outPdfFile.c_str()); 
    
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

    std::vector<std::string> desclines={ "#sqrt{s} = 5.02 TeV",
					 "ak4PF Jets",
					 "#||{y} < 2.0"};
    
    // gen/meas/unf spectra hists ---------------    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(1);    
    
    setupSpectraHist(hgen_rebin          , useSimpBins);    
    setupSpectraHist(hrec_rebin	       , useSimpBins);
    setupSpectraHist(hrec_sameside_rebin , useSimpBins);
    setupSpectraHist(hunf		       , useSimpBins);    
    
    hrec_rebin->SetTitle(("Jet Spectra"+methodString+descString).c_str());
    
    
    hrec_rebin->DrawClone("P E");           
    hgen_rebin->DrawClone("P E SAME");                 
    hunf->DrawClone("P E SAME");               
    hrec_sameside_rebin->DrawClone("P E SAME");           
    
    TLegend* legend_in1 = new TLegend( 0.55,0.7,0.9,0.9 );
    legend_in1->AddEntry(hunf,                "Data Unf." ,  "lp");
    legend_in1->AddEntry(hrec_rebin,          "Data Meas." , "lp");	
    legend_in1->AddEntry(hrec_sameside_rebin, ("Smeared "+NLOMCtitle_str).c_str()  ,  "lp");	 
    legend_in1->AddEntry(hgen_rebin,          ("Toy "+NLOMCtitle_str+" Truth").c_str()   , "lp");
    legend_in1->SetBorderSize(0);
    legend_in1->SetFillStyle(0);
    legend_in1->Draw();    
    drawTLatex( 0.42, 0.83, 0.035, canvForPrint, desclines);
    
    canv_spectra=(TCanvas*)canvForPrint->DrawClone();
    canvForPrint->Print(outPdfFile.c_str());      
    
    
    
    // ratios w PY8 Truth------------    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(0);
    
    setupRatioHist(h_genratio_oppunf , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_genratio_oppmeas, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_genratio_ssmeas , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_genratio_oppfold, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    
    h_genratio_oppunf->SetTitle(("Toy NLO Truth Ratio"+methodString+descString).c_str());
    h_genratio_oppunf->GetYaxis()->SetTitle("Ratio w/ Toy NLO Truth");    
    
    h_genratio_oppunf->DrawClone("P E");            //data unf/mc truth
    h_genratio_oppmeas->DrawClone("P E SAME");      // data meas/mc truth
    h_genratio_ssmeas->DrawClone("P E SAME");       // mc meas/mc truth
    //h_genratio_oppfold->DrawClone("P E SAME");
    
    TLegend* legend2 = new TLegend( 0.55,0.75,0.95,0.9 );
    legend2->AddEntry(h_genratio_oppunf,  "Data Unf." ,  "lp");
    legend2->AddEntry(h_genratio_oppmeas, "Data Meas." , "lp"); 
    legend2->AddEntry(h_genratio_ssmeas, ("Smeared "+NLOMCtitle_str).c_str(), "lp");

    //legend2->AddEntry(h_genratio_oppfold, "Data Fold(Unf.)", "lp");
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->Draw();    
    drawTLatex( 0.42, 0.83, 0.035, canvForPrint, desclines);
    
    theLineAtp9 ->Draw();
    theLineAtOne->Draw();
    theLineAt1p1->Draw();
    
    canv_gen_ratio=(TCanvas*)canvForPrint->DrawClone();
    canvForPrint->Print(outPdfFile.c_str());





    // ratio w/ data meas -------------------
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(0);      
    
    //h_recratio_oppunf->SetTitle("Bayesian, Ratios w/ Data Meas. Spectra");
    //h_recratio_oppunf->GetYaxis()->SetTitle("Ratio w/ Data Meas.");
    h_recratio_ssmeas->SetTitle(("Data Meas. Ratio"+methodString+descString).c_str());
    h_recratio_ssmeas->GetYaxis()->SetTitle("Ratio w/ Data Meas.");
    
    setupRatioHist(h_recratio_ssmeas , useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);    
    setupRatioHist(h_recratio_oppunf , useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);    
    setupRatioHist(h_recratio_ssgen  , useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);    
    
    h_recratio_ssmeas->DrawClone("P E");
    h_recratio_oppunf->DrawClone("P E SAME");
    h_recratio_ssgen->DrawClone("P E SAME");
    ////h_recratio_ssunf->DrawClone("P E SAME");
    //h_recratio_oppfold->DrawClone("P E SAME");  
    //h_recratio_ssfold->DrawClone("P E SAME");
    
    TLegend* legend4 = new TLegend( 0.55,0.75,0.9,0.9);
    legend4->AddEntry(h_recratio_oppunf,  "Data Unf.", "lp");
    legend4->AddEntry(h_recratio_ssmeas, ("Smeared "+NLOMCtitle_str).c_str()  ,  "lp");	 
    legend4->AddEntry(h_recratio_ssgen,  ("Toy "+NLOMCtitle_str+" Truth").c_str()   , "lp");
    ////legend4->AddEntry(h_recratio_ssunf, NULL, "p");
    //legend4->AddEntry(h_recratio_oppfold, "Data Fold(Unf.)", "lp");
    //legend4->AddEntry(h_recratio_ssfold, "PY8 Fold(Unf.)", "p");    
    legend4->SetBorderSize(0);
    legend4->SetFillStyle(0);
    legend4->Draw();    
    drawTLatex( 0.42, 0.83, 0.035, canvForPrint, desclines);

    theLineAtp9->Draw();
    theLineAtOne->Draw();
    theLineAt1p1->Draw();
    
    canv_rec_ratio=(TCanvas*)canvForPrint->DrawClone();
    canvForPrint->Print(outPdfFile.c_str());



    // Folding Ratios (fold+fakes) -------------------
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(0);
    
    setupRatioHist(h_foldratio_datafold, useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);
    setupRatioHist(h_foldratio_mcfold, useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);
    
    h_foldratio_datafold->SetTitle(("Folding Ratios"+methodString+descString).c_str());// "Folded Ratios w/ (Corrected) Meas. Spectra" );
    h_foldratio_datafold->GetYaxis()->SetTitle("Ratio w/ Meas.");
    
    h_foldratio_datafold->DrawClone("P E");
    h_foldratio_mcfold->DrawClone("P E SAME");
    
    TLegend* legendfold = new TLegend( 0.5,0.8,0.9,0.9);
    legendfold->AddEntry(h_foldratio_datafold ,  "Fold(Data Unf.) + Fakes/Data Meas." , "lp");
    legendfold->AddEntry(h_foldratio_mcfold  ,   "Fold(Toy NLO Truth) + Fakes/Smeared Toy NLO",  "lp");
    legendfold->SetBorderSize(0);
    legendfold->SetFillStyle(0);
    legendfold->Draw();    
    drawTLatex( 0.37, 0.85, 0.035, canvForPrint, desclines);
    
    theLineAtp9->Draw();
    theLineAtOne->Draw();
    theLineAt1p1->Draw();

    canv_fold_ratio=(TCanvas*)canvForPrint->DrawClone();
    canvForPrint->Print(outPdfFile.c_str());


    // Folding Ratios (meas-fakes) -------------------
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(0);
    
    setupRatioHist(h_recratio_oppfold, useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);
    setupRatioHist(h_recratio_truthfold, useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);
    
    h_recratio_oppfold->SetTitle(("Folding Ratios"+methodString+descString).c_str());// "Folded Ratios w/ (Corrected) Meas. Spectra" );
    h_recratio_oppfold->GetYaxis()->SetTitle("Ratio w/ (Meas. - Fakes)");
    
    h_recratio_oppfold->DrawClone("P E");
    h_recratio_truthfold->DrawClone("P E SAME");
    
    TLegend* legendfold2 = new TLegend( .5,.8,.9,.9);
    legendfold2->AddEntry(h_recratio_oppfold ,  "Fold(Data Unf.)/(Data Meas. - Fakes)" , "lp");
    legendfold2->AddEntry(h_recratio_truthfold  ,   "Fold(Toy NLO Truth)/(Smeared Toy NLO - Fakes)",  "lp");
    legendfold2->SetBorderSize(0);
    legendfold2->SetFillStyle(0);
    legendfold2->Draw();    
    drawTLatex( 0.37, 0.85, 0.035, canvForPrint, desclines);
    
    theLineAtp9->Draw();
    theLineAtOne->Draw();
    theLineAt1p1->Draw();

    canv_fold_ratio2=(TCanvas*)canvForPrint->DrawClone();
    canvForPrint->Print(outPdfFile.c_str());
    
    
    // gen/reco response spectra hists ---------------    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(1);    
    
    setupSpectraHist(hrec_sameside_resp_rebin , useSimpBins);
    setupSpectraHist(hgen_resp_rebin          , useSimpBins);
    setupSpectraHist(hfak                     , useSimpBins);
    
    //histograms given to roo_resp + the fakes histo made from these spectra + response matrix projection
    hfak->SetTitle( (NLOMCtitle_str+" Only, Kinematic Fakes Spectra").c_str());
    hfak->SetMaximum(2.*(hrec_sameside_resp_rebin->GetMaximum()));

    hfak->DrawClone("P E");           
    hrec_sameside_resp_rebin->DrawClone("P E SAME");           
    hgen_resp_rebin->DrawClone("P E SAME");                 

    TLegend* legend_resp = new TLegend( 0.65,0.75,0.9,0.9 );
    legend_resp->AddEntry(hfak,          "Toy NLO Kinematic Fakes" , "lp");    
    legend_resp->AddEntry(hgen_resp_rebin,          "Toy NLO Truth Response" , "lp");
    legend_resp->AddEntry(hrec_sameside_resp_rebin, "Smeared Toy NLO Response" , "lp");    
    legend_resp->SetBorderSize(0);
    legend_resp->SetFillStyle(0);
    legend_resp->Draw();    
    drawTLatex( 0.54, 0.84, 0.035, canvForPrint, desclines);
    
    canv_mc_fakes_spectra=(TCanvas*)canvForPrint->DrawClone();
    canvForPrint->Print(outPdfFile.c_str());           
    
    
    // thy spectra CT10/14 NNPDF NLO---------------------------
    TLegend* legendThy1 =new TLegend( 0.50,0.65,0.90,0.9 );    
    
    setupSpectraHist(CT10nlo  ,useSimpBins);
    if(debugWrite){fout->cd(); CT10nlo->Write("CT10_NLO_R04_jtpt");}    
    setupSpectraHist(CT14nlo  ,useSimpBins);
    if(debugWrite){fout->cd(); CT14nlo->Write("CT14_NLO_R04_jtpt");}
    setupSpectraHist(NNPDFnnlo,useSimpBins);        
    if(debugWrite){fout->cd(); NNPDFnnlo->Write("NNPDF_NLO_R04_jtpt"); }    
    if(debugWrite)fout->cd();//makeThyHist will make me cd to last thy file opened...
    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(1);                  
    

    //CT10nlo->SetAxisRange(1e-07,1e+03,"Y");//for y axis in nanbarns
    if(applyNPCorrs)    CT10nlo->SetTitle(("NLO+NPs Jet Spectra"+methodString+descString).c_str());
    else    CT10nlo->SetTitle(("NLO Jet Spectra"+methodString+descString).c_str());
    
    CT10nlo  ->DrawClone("][HIST E");
    CT14nlo  ->DrawClone("][HIST E SAME");    
    NNPDFnnlo->DrawClone("][HIST E SAME");    
    hgen_rebin->DrawClone("P E SAME");
    hunf->DrawClone("P E SAME");   //just for axis range
    
    legendThy1->AddEntry(CT10nlo  ,(  "CT10 NLO"+CT10NPs).c_str(),"l");
    legendThy1->AddEntry(CT14nlo  ,(  "CT14 NLO"+CT14NPs).c_str(),"l");
    legendThy1->AddEntry(NNPDFnnlo,("NNPDF NNLO"+NNPDFNPs).c_str(),"l");
    legendThy1->AddEntry(hunf,"Data Unf.","lp");
    legendThy1->AddEntry(hgen_rebin,("Toy "+NLOMCtitle_str+" Truth").c_str(), "lp");
    legendThy1->SetBorderSize(0);
    legendThy1->SetFillStyle(0);
    legendThy1->Draw();    
    drawTLatex( 0.40, 0.84, 0.035, canvForPrint, desclines);
    
    canv_thy_spectra_1=(TCanvas*)canvForPrint->DrawClone();
    canvForPrint->Print(outPdfFile.c_str());

    // thy spectra HERA/MMHT/NNPDF ---------------------------

    setupSpectraHist(HERAPDF  ,useSimpBins);
    if(debugWrite){fout->cd(); HERAPDF->Write("HERAPDF105_NLO_R04_jtpt");}    
    setupSpectraHist(MMHTnlo  ,useSimpBins);
    if(debugWrite){fout->cd(); MMHTnlo->Write("MMHT2014_NLO_R04_jtpt"); }
    if(debugWrite)fout->cd();//makeThyHist will make me cd to last thy file opened...
    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(1);                  
    
    if(applyNPCorrs)    MMHTnlo->SetTitle(("NLO+NPs Jet Spectra"+methodString+descString).c_str());
    else    MMHTnlo->SetTitle(("NLO Jet Spectra"+methodString+descString).c_str());    
    MMHTnlo  ->DrawClone("][HIST E");
    HERAPDF  ->DrawClone("][HIST E SAME");

    
    hgen_rebin->DrawClone("P E SAME");
    hunf->DrawClone("P E SAME");   //just for axis range
    
    TLegend* legendThy =new TLegend( 0.50,0.70,0.90,0.9 );        
    legendThy->AddEntry(MMHTnlo  ,(   "MMHT 2014 NLO"+MMHTNPs).c_str(),"l");    
    legendThy->AddEntry(HERAPDF  ,("HERAPDF 2015 NLO"+HERANPs).c_str(),"l");
    legendThy->AddEntry(hunf,"Data Unf.","lp");
    legendThy->AddEntry(hgen_rebin,("Toy "+NLOMCtitle_str+" Truth").c_str(), "lp");    
    legendThy->SetBorderSize(0);
    legendThy->SetFillStyle(0);
    legendThy->Draw();    
    drawTLatex( 0.40, 0.84, 0.035, canvForPrint, desclines);
    

    canv_thy_spectra_2=(TCanvas*)canvForPrint->DrawClone();
    canvForPrint->Print(outPdfFile.c_str());

    // thy ratios w hunf------------    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(0);
    
    setupRatioHist(h_thyratio_CT10nlo , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_thyratio_CT14nlo , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_thyratio_NNPDFnnlo, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_thyratio_mctruth, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    
    if(applyNPCorrs) h_thyratio_CT10nlo->SetTitle(("NLO+NPs, Ratios w/ Unf. Data"+methodString+descString).c_str()); 
    else h_thyratio_CT10nlo->SetTitle(("NLO, Ratios w/ Unf. Data"+methodString+descString).c_str()); 
    
    h_thyratio_CT10nlo->GetYaxis()->SetTitle("Thy / Unf. Data");

    h_thyratio_CT10nlo ->DrawClone( "][HIST ");      
    h_thyratio_CT14nlo ->DrawClone( "][HIST SAME"); 
    h_thyratio_NNPDFnnlo->DrawClone("][HIST E SAME"); 
    h_thyratio_mctruth->DrawClone("P E SAME");     
    
    TLegend* legendthyrat = new TLegend( 0.1,0.7,0.5,0.9 );
    legendthyrat->AddEntry(h_thyratio_CT10nlo ,    ("CT10 PDF NLO"    +CT10NPs).c_str(),    "l");
    legendthyrat->AddEntry(h_thyratio_CT14nlo ,    ("CT14 PDF NLO"    +CT14NPs).c_str(),    "l"); 
    legendthyrat->AddEntry(h_thyratio_NNPDFnnlo,   ("NNPDF NNLO"      +NNPDFNPs).c_str(),       "l");
    legendthyrat->AddEntry(h_thyratio_mctruth, ("Toy "+NLOMCtitle_str+" Truth").c_str(),       "lp");    
    legendthyrat->SetBorderSize(0);
    legendthyrat->SetFillStyle(0);
    legendthyrat->Draw();
    drawTLatex( 0.50, 0.84, 0.035, canvForPrint, desclines);    
    
    theLineAtp9 ->Draw();
    theLineAtOne->Draw();
    theLineAt1p1->Draw();
    
    canv_thy_ratio=(TCanvas*)canvForPrint->DrawClone();
    canvForPrint->Print(outPdfFile.c_str());
    
    // thy ratios w hunf------------    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(0);
    
    setupRatioHist(h_thyratio_HERAPDF , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_thyratio_MMHTnlo , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    
    if(applyNPCorrs) h_thyratio_HERAPDF->SetTitle(("NLO+NPs, Ratios w/ Unf. Data"+methodString+descString).c_str()); 
    else h_thyratio_HERAPDF->SetTitle(("NLO, Ratios w/ Unf. Data"+methodString+descString).c_str()); 
    
    h_thyratio_HERAPDF->GetYaxis()->SetTitle("Thy / Unf. Data");
    
    h_thyratio_HERAPDF ->DrawClone( "][HIST E"); 
    h_thyratio_MMHTnlo ->DrawClone( "][HIST E SAME"); 
    h_thyratio_mctruth->DrawClone("P E SAME");     
    
    TLegend* legendthyrat2 = new TLegend( 0.1,0.7,0.5,0.9 );
    legendthyrat2->AddEntry(h_thyratio_HERAPDF ,  ("HERAPDF 2015 NLO"+HERANPs).c_str(), "l");
    legendthyrat2->AddEntry(h_thyratio_MMHTnlo ,  ("MMHT 2014 NLO"   +MMHTNPs).c_str(),    "l");
    legendthyrat2->AddEntry(h_thyratio_mctruth, ("Toy "+NLOMCtitle_str+" Truth").c_str(),       "lp");    
    legendthyrat2->SetBorderSize(0);
    legendthyrat2->SetFillStyle(0);
    legendthyrat2->Draw();
    drawTLatex( 0.50, 0.84, 0.035, canvForPrint, desclines);    
    
    theLineAtp9 ->Draw();
    theLineAtOne->Draw();
    theLineAt1p1->Draw();
    
    canv_thy_ratio2=(TCanvas*)canvForPrint->DrawClone();
    canvForPrint->Print(outPdfFile.c_str());
    
    
    // NJets v. Jet p_T hist (all weight=1) ---------------------    
    if( (bool)hJetQA_jtptEntries )      {
      TLine* theLineAt100_reco= new TLine( boundaries_pt_reco_mat[0]  
					   ,100.   
					   ,(boundaries_pt_reco_mat[nbins_pt_reco_mat])
					   ,100.);
      if(debugMode)std::cout<<"theLineAt1p1_reco starts at "<<boundaries_pt_reco_mat[0]  <<std::endl;
      if(debugMode)std::cout<<"theLineAt1p1_reco ends at   "<<(boundaries_pt_reco_mat[nbins_pt_reco_mat])  <<std::endl;
      theLineAt100_reco->SetLineWidth(1);
      theLineAt100_reco->SetLineStyle(2);
      theLineAt100_reco->SetLineColor(36);
      
      canvForPrint->cd();
      if(!useSimpBins)canvForPrint->SetLogx(1);
      canvForPrint->SetLogy(1);
      
      setupSpectraHist(hJetQA_jtptEntries, useSimpBins);

      hJetQA_jtptEntries->SetTitle("Jet Statistics, Weights=1");
      hJetQA_jtptEntries->GetYaxis()->SetTitle("N_Jets [weight=1]");      
      hJetQA_jtptEntries->DrawClone("HIST E");
      
      theLineAt100_reco->Draw();
      
      canvForPrint->Print(outPdfFile.c_str());    }
    

    // covariance matrix ---------------      //always use simp bins for covmat to avoid log scaling the x/y axes
    canvForPrint->cd();
    canvForPrint->SetLogx(0);
    canvForPrint->SetLogy(0);
    canvForPrint->SetLogz(1);
        
    matStylePrint( (TH2D*)covmat_TH2, ("Covariance Matrix"+methodString+descString), canvForPrint, outPdfFile, true);        
    canv_covmat=(TCanvas*)canvForPrint->DrawClone(); 
    

    
    matStylePrint( (TH2D*)covmatabsval_TH2, ("|Covariance Matrix|"+methodString+descString), canvForPrint, outPdfFile, true);
    canv_absval_covmat=(TCanvas*)canvForPrint->DrawClone(); 
    
    // pearson matrix ---------------      //always use simp bins for covmat to avoid log scaling the x/y axes
    canvForPrint->cd();
    canvForPrint->SetLogx(0);
    canvForPrint->SetLogy(0);
    canvForPrint->SetLogz(0);
    
    matStylePrint( (TH2D*)PearsonBayes, ("Pearson Matrix"+methodString+descString), canvForPrint, outPdfFile, true);
    canv_pearson=(TCanvas*)canvForPrint->DrawClone(); 
    
    // unfolding matrix ---------------      //always use simp bins for covmat to avoid log scaling the x/y axes
    canvForPrint->cd();
    canvForPrint->SetLogx(0);
    canvForPrint->SetLogy(0);
    canvForPrint->SetLogz(1);

    matStylePrint( (TH2D*)unfmat_TH2, ("Unfolding Matrix"+methodString+descString), canvForPrint, outPdfFile, true);
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
    matStylePrint(hmat_rebin, (NLOMCtitle_str+" Response Matrix, Rebinned"), canvForPrint, outPdfFile, useSimpBins);
    canv_mat_rebin=(TCanvas*)canvForPrint->DrawClone(); 
    
    
    // percent error matrix in binning of interest ---------------    
    hmat_percenterrs= makeRespMatrixPercentErrs( (TH2D*) hmat_errors, (TH2D*) hmat_rebin,
						 (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
						 (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );		     
    hmat_percenterrs->SetName("hmat_percenterrs_foroutputpdf");
    if(debugWrite) hmat_percenterrs->Write();
    if(debugMode)  hmat_percenterrs->Print("base");
    
    canvForPrint->cd();
    canvForPrint->SetLogz(1);
    
    if(!useSimpBins){
      canvForPrint->SetLogx(1);    
      canvForPrint->SetLogy(1);}
    else{
      canvForPrint->SetLogx(0);    
      canvForPrint->SetLogy(0);          }
    
    matStylePrint(hmat_percenterrs, (NLOMCtitle_str+" Response Matrix, Bin % Errors"), canvForPrint, outPdfFile, useSimpBins);
    canv_mat_percerrs=(TCanvas*)canvForPrint->DrawClone(); 
    



    if(dokIterQA){
      if(debugMode)std::cout<<"drawing plots across diff kIter values"<<std::endl;
      if(debugMode)std::cout<<kIter_start<<" <= kIter <= "<<kIter_end<<std::endl;
      canv_3x3spectra  =new TCanvas("canv_3x3spectra","canv 3x3 spectra", 1500, 1500);   canv_3x3spectra ->Divide(3,3);
      canv_3x3genratio =new TCanvas("canv_3x3genratio","canv 3x3 genratio", 1500, 1500); canv_3x3genratio->Divide(3,3);
      canv_3x3recratio =new TCanvas("canv_3x3recratio","canv 3x3 recratio", 1500, 1500); canv_3x3recratio->Divide(3,3);                  
      canv_3x3covmat       =new TCanvas("canv_3x3covmat","canv 3x3 covmat", 1500, 1500); canv_3x3covmat      ->Divide(3,3);                  
      canv_3x3covmatabsval =new TCanvas("canv_3x3covmatabsval","canv 3x3 covmatabsval", 1500, 1500); canv_3x3covmatabsval->Divide(3,3);                  
      canv_3x3pearson      =new TCanvas("canv_3x3pearson","canv 3x3 pearson", 1500, 1500); canv_3x3pearson     ->Divide(3,3);                  
      TLegend* leg_3x3=new TLegend(0.62,0.75,0.9,0.9,NULL,"NBNDC");
      
      for(int ki=0; ki<nKiterMax;ki++){
	int current_kIter=kIter_start+ki;
	int ki_canv=ki+1;
	if(debugMode)std::cout<<"ki="<<ki<<", current_kIter="<<current_kIter<<std::endl;
 
	canv_3x3spectra->cd(ki_canv)->SetLogx(1);
	canv_3x3spectra->cd(ki_canv)->SetLogy(1);
	canv_3x3spectra->cd(ki_canv);
	
	hunf_bayes[ki]->GetXaxis()->SetNoExponent(true);
	hunf_bayes[ki]->GetXaxis()->SetMoreLogLabels(true);
	
	hunf_bayes[ki]->DrawClone("HIST E");
	//hfold_bayes[ki]->DrawClone("HIST E SAME");
 
	leg_3x3->Clear();//removes any entries in legend
	leg_3x3->AddEntry(hunf_bayes[ki] ,"Unf."      ,"lp");
	//leg_3x3->AddEntry(hfold_bayes[ki],"Fold(Unf.)","lp");
	//leg_3x3->DrawClone();
 
	canv_3x3genratio->cd(ki_canv)->SetLogx(1);
	canv_3x3genratio->cd(ki_canv)->SetLogy(0);
	canv_3x3genratio->cd(ki_canv);
	hgen_unf_ratio[ki]->SetAxisRange(0.5,1.5,"Y");
 

	hgen_unf_ratio[ki]->GetXaxis()->SetNoExponent(true);
	hgen_unf_ratio[ki]->GetXaxis()->SetMoreLogLabels(true);

	hgen_unf_ratio[ki]->DrawClone("HIST");
	theLineAtp9 ->Draw();
	theLineAtOne->Draw();
	theLineAt1p1->Draw();	
	hgen_unf_ratio[ki]->DrawClone("HIST E SAME");
	//hgen_fold_ratio[ki]->DrawClone("HIST E SAME");
	
	
	leg_3x3->Clear();//removes any entries in legend
	leg_3x3->AddEntry(hgen_unf_ratio[ki] ,"Unf."      ,"lp");
	//leg_3x3->AddEntry(hgen_fold_ratio[ki],"Fold(Unf.)","lp");
	//leg_3x3->DrawClone();
 
	canv_3x3recratio->cd(ki_canv)->SetLogx(1);
	canv_3x3recratio->cd(ki_canv)->SetLogy(0);
	canv_3x3recratio->cd(ki_canv);
	
	hrec_unf_ratio[ki]->GetXaxis()->SetNoExponent(true);
	hrec_unf_ratio[ki]->GetXaxis()->SetMoreLogLabels(true);

	hrec_unf_ratio[ki]->SetAxisRange(0.5,1.5,"Y");
	
	hrec_unf_ratio[ki]->DrawClone("HIST");
	theLineAtp9 ->Draw();
	theLineAtOne->Draw();
	theLineAt1p1->Draw();
	hrec_unf_ratio[ki]->DrawClone("HIST E SAME");
	//hrec_fold_ratio[ki]->DrawClone("HIST E SAME"); 
 
	leg_3x3->Clear();//removes any entries in legend
	leg_3x3->AddEntry(hrec_unf_ratio[ki] ,"Unf."      ,"lp");
	//leg_3x3->AddEntry(hrec_fold_ratio[ki],"Fold(Unf.)","lp");
	//leg_3x3->DrawClone();

	//canv_3x3covmat      
	canv_3x3covmat->cd(ki_canv)->SetLogx(0);
	canv_3x3covmat->cd(ki_canv)->SetLogy(0);
	canv_3x3covmat->cd(ki_canv)->SetLogz(1);
	canv_3x3covmat->cd(ki_canv);	
	hcovmat_bayes[ki]->Draw("COLZ");
	
	//canv_3x3covmatabsval
	canv_3x3covmatabsval->cd(ki_canv)->SetLogx(0);
	canv_3x3covmatabsval->cd(ki_canv)->SetLogy(0);
	canv_3x3covmatabsval->cd(ki_canv)->SetLogz(1);
	canv_3x3covmatabsval->cd(ki_canv);	
	hcovmatabsval_bayes[ki]->Draw("COLZ");
	//canv_3x3pearson     	
	canv_3x3pearson->cd(ki_canv)->SetLogx(0);
	canv_3x3pearson->cd(ki_canv)->SetLogy(0);
	canv_3x3pearson->cd(ki_canv)->SetLogz(0);
	canv_3x3pearson->cd(ki_canv);	
	hpearson_bayes[ki]->Draw("COLZ");	
      }
      
      
      canv_3x3spectra ->Print(outPdfFile.c_str());
      canv_3x3genratio->Print(outPdfFile.c_str());
      canv_3x3recratio->Print(outPdfFile.c_str());
      canv_3x3covmat->Print(outPdfFile.c_str());
      canv_3x3covmatabsval->Print(outPdfFile.c_str());
      canv_3x3pearson->Print(outPdfFile.c_str());
      
      
      
      canv_kIterRatio=new TCanvas("canv_kiter_ratio","canv kiter ratio" , 1500, 1500);

      canv_kIterRatio->SetName(("canv_hunf_ratio_kIter"+std::to_string(kIter_start)+"_v_kIter"+std::to_string(kIter_end)).c_str());
      canv_kIterRatio->SetTitle(("Unf. Ratio, kIter="+std::to_string(kIter_start)+"/kIter="+std::to_string(kIter_end)+" Canvas").c_str());
      canv_kIterRatio->cd()->SetLogx(1);
      canv_kIterRatio->cd()->SetLogy(0);
      canv_kIterRatio->cd();
      

      hkiter_ratio=(TH1D*)hunf_bayes[0]->Clone();
      hkiter_ratio->Divide(hunf_bayes[nKiterMax-1]);
      hkiter_ratio->SetName(("hunf_ratio_kIter"+std::to_string(kIter_start)+"_v_"+std::to_string(kIter_end)).c_str());
      hkiter_ratio->SetTitle(("Unf. Ratio, kIter="+std::to_string(kIter_start)+"/kIter="+std::to_string(kIter_end)+";Jet p_{T};ratio").c_str());
      hkiter_ratio->GetXaxis()->SetNoExponent(true);
      hkiter_ratio->GetXaxis()->SetMoreLogLabels(true);
      hkiter_ratio->SetAxisRange(0.5,1.5,"Y");
      
      hkiter_ratio->DrawClone("HIST E");
      theLineAtp9 ->Draw();
      theLineAtOne->Draw();
      theLineAt1p1->Draw();
      hkiter_ratio->DrawClone("HIST E SAME");

      
      canv_kIterRatio->Print(outPdfFile.c_str());
      
      
    }
    
    
    // for the resp matrix file drawing only
    if(drawRespMatrix){    
      drawRespMatrixFile((TH2D*)hmat, (TH2D*)hmat_rebin, (TH2D*)hmat_errors,
			 //(TH2D*)covmat_TH2, (TH2D*)PearsonBayes, (TH2D*)unfmat_TH2,
			 (double*)boundaries_pt_reco_mat, (int)nbins_pt_reco_mat,
			 (double*)boundaries_pt_gen_mat , (int)nbins_pt_gen_mat,
			 (std::string)outRespMatPdfFile , (bool)useSimpBins, (TFile*)fpp_MC);    }//end draw response matrix
    
    std::cout<<std::endl<<"draw pdfs for Bayesian Unfolding done!"<<std::endl;
  


  }// end draw pdfs
  
  std::cout<<"writing output hists to file... "<<std::endl;

  fout->cd();
  

  
// input data ------------------
  hrec_rebin->SetTitle( "Data Meas.");  hrec_rebin->Write("Data_meas");
  if((bool)hJetQA_jtptEntries){ 
    hJetQA_jtptEntries->SetTitle("Data N_{Jets}");hJetQA_jtptEntries->Write("Data_njets");}
  
  // input PY8 ---------------------
  hgen_rebin->SetTitle("PY8 Truth");  hgen_rebin->Write("MC_truth");
  hrec_sameside_rebin->SetTitle("PY8 Meas.");  hrec_sameside_rebin->Write("MC_meas");
  
  hmat->SetTitle("PY8 Response Matrix");  hmat->Write("MC_mat");
  hmat_rebin->SetTitle("PY8 Response Matrix Rebinned");  hmat_rebin->Write("MC_mat_rebin"); 
  hmat_errors->SetTitle("PY8 Response Matrix Errors");  hmat_errors->Write("MC_mat_rebin_errors");
  if((bool)hmat_percenterrs){
    hmat_percenterrs->SetTitle("PY8 Response Matrix Percent Errors");hmat_percenterrs->Write("MC_mat_rebin_percerrors");  }
  
  // input thy ---------------- 
  CT10nlo  ->SetTitle("CT10 NLO Spectra");         CT10nlo  ->Write("NLO_CT10_NLO_R04_jtpt");	      
  CT14nlo  ->SetTitle("CT14 NLO Spectra");	   CT14nlo  ->Write("NLO_CT14_NLO_R04_jtpt");	      
  HERAPDF  ->SetTitle("HERAPDF NLO Spectra");	   HERAPDF  ->Write("NLO_HERAPDF105_NLO_R04_jtpt");  
  MMHTnlo  ->SetTitle("MMHT NLO Spectra");	   MMHTnlo  ->Write("NLO_MMHT2014_NLO_R04_jtpt");    
  NNPDFnnlo->SetTitle("NNPDF NNLO Spectra");	   NNPDFnnlo->Write("NLO_NNPDF_NLO_R04_jtpt");       
  
  // output hists -------------
  hfak->SetTitle("PY8 Meas. Fakes");hfak->Write("MC_meas_fakes");
  
  hunf->SetTitle(("Data Unf. kIter="+std::to_string(kIterInput)).c_str() );hunf->Write("Data_unf");      
  hfold->SetTitle(("Data Fold(Unf.) kIter="+std::to_string(kIterInput)).c_str() ); hfold->Write("Data_fold");        
  hfold_fakecorr->SetTitle(("Data Fold(Unf.) + Fakes, kIter="+std::to_string(kIterInput)).c_str() ); hfold_fakecorr->Write("Data_foldfakcorr");        
  hrec_rebin_fakecorr->SetTitle("Data Fake Corr. Meas.");hrec_rebin_fakecorr->Write("Data_measfakcorr");  
  
  hfold_truth->SetTitle(("PY8 Fold(Truth), kIter="+std::to_string(kIterInput)).c_str() );hfold_truth->Write("MC_truth_fold");
  hfold_truth_fakecorr->SetTitle(("PY8 Fold(Truth) + Fakes, kIter="+std::to_string(kIterInput)).c_str() );hfold_truth_fakecorr->Write("MC_truth_foldfakcorr");
  hrec_sameside_rebin_fakecorr->SetTitle("PY8 Fake Corr. Meas."); hrec_sameside_rebin_fakecorr->Write("MC_measfakcorr");
  
  covmat_TH2->Write("covmat");
  covmatabsval_TH2->Write("covmatabsval");
  PearsonBayes->Write("pearson");
  unfmat_TH2->Write("unfmat");
  
  // output ratio comparisons -------------
  // gen ratios (denom=mc truth)
  h_genratio_oppunf ->SetTitle("Data Unf./PY8 Truth");h_genratio_oppunf ->Write("ratio_Data_unf_MC_truth");  //data unf/mc truth	   
  h_genratio_oppfold->SetTitle("Data Fold(Unf.)/PY8 Truth");h_genratio_oppfold->Write("ratio_Data_fold_MC_truth"); //data fold(unf)/mc truth 
  h_genratio_oppmeas->SetTitle("Data Meas./PY8 Truth");h_genratio_oppmeas->Write("ratio_Data_meas_MC_truth"); //data meas/mc truth	   
  h_genratio_ssmeas ->SetTitle("PY8 Meas./PY8 Truth");h_genratio_ssmeas ->Write("ratio_MC_meas_MC_truth");  //mc meas/mc truth        
  
  // rec ratios (denom=data meas)
  h_recratio_oppunf ->SetTitle("Data Unf./Data Meas.");h_recratio_oppunf ->Write("ratio_Data_unf_Data_meas");  //data unf       / data meas	   
  h_recratio_ssmeas ->SetTitle("PY8 Meas./Data Meas.");h_recratio_ssmeas ->Write("ratio_MC_meas_Data_meas");  //mc meas        / data meas	   
  h_recratio_ssgen  ->SetTitle("PY8 Truth/Data Meas.");h_recratio_ssgen  ->Write("ratio_MC_truth_Data_meas");   //mc truth       / data meas        
  
  // fold ratio test
  h_recratio_oppfold   ->SetTitle("Data Fold(Unf.)/(Data Meas. - Fakes)");  h_recratio_oppfold   ->Write(); //data fold(unf) / data meas - fakes    
  h_recratio_truthfold ->SetTitle("PY8 Fold(Truth)/(PY8 Meas. - Fakes)");      h_recratio_truthfold  ->Write(); //PY8 fold(truth) / data meas - fakes    
  h_foldratio_datafold ->SetTitle("(Data Fold(Unf.) + Fakes)/Data Meas."); h_foldratio_datafold ->Write();
  h_foldratio_mcfold   ->SetTitle("(PY8 Fold(Truth) + Fakes)/PY8 Meas.");    h_foldratio_mcfold   ->Write();
  
  // thy ratios w/ unfolded data
  h_thyratio_CT10nlo  ->Write("ratio_CT10_NLO_Data_unf");
  h_thyratio_CT14nlo  ->Write("ratio_CT14_NLO_Data_unf");
  h_thyratio_HERAPDF  ->Write("ratio_HERAPDF_NLO_Data_unf");
  h_thyratio_MMHTnlo  ->Write("ratio_MMHTnlo_NLO_Data_unf");    
  h_thyratio_NNPDFnnlo->Write("ratio_NNPDFnnlo_NNLO_Data_unf");
  h_thyratio_mctruth->Write("ratio_PY8_MCTruth_Data_unf");
  
  if(dokIterQA){
    hkiter_ratio->Write();
    
    TDirectory* fout_kiter_dir=fout->mkdir("all_kiter_plots");
    fout_kiter_dir->cd();
    for(int ki=0; ki<nKiterMax;ki++) hunf_bayes[ki]->Write();
    for(int ki=0; ki<nKiterMax;ki++) hfold_bayes[ki]->Write();
    for(int ki=0; ki<nKiterMax;ki++) hgen_unf_ratio[ki]->Write();
    for(int ki=0; ki<nKiterMax;ki++) hgen_fold_ratio[ki]->Write();
    for(int ki=0; ki<nKiterMax;ki++) hrec_unf_ratio[ki]->Write();
    for(int ki=0; ki<nKiterMax;ki++) hrec_fold_ratio[ki]->Write();
    for(int ki=0; ki<nKiterMax;ki++) hcovmat_bayes[ki]->Write();
    for(int ki=0; ki<nKiterMax;ki++) hcovmatabsval_bayes[ki]->Write();
    for(int ki=0; ki<nKiterMax;ki++) hpearson_bayes[ki]->Write();
    fout->cd();//return to base dir of file
  }



  
  // canvases ----------------
  if(drawPDFs){
    canv_spectra          ->SetTitle("I/O Spectra Canvas");        canv_spectra           ->Write("canv_spectra");
    canv_mc_fakes_spectra ->SetTitle("PY8 Fakes Spectra Canvas");   canv_mc_fakes_spectra  ->Write("canv_mc_fakes_spectra");
    canv_thy_spectra_1    ->SetTitle("NLO Thy Spectra 1 Canvas");  canv_thy_spectra_1     ->Write("canv_thy_spectra_1");
    canv_thy_spectra_2    ->SetTitle("NLO Thy Spectra 2 Canvas");  canv_thy_spectra_2     ->Write("canv_thy_spectra_2");                                
    
    canv_gen_ratio        ->SetTitle("PY8 Truth Ratios Canvas");   canv_gen_ratio          ->Write("canv_gen_ratio");
    canv_rec_ratio        ->SetTitle("Data Meas Ratios Canvas");  canv_rec_ratio          ->Write("canv_meas_ratio");
    canv_fold_ratio       ->SetTitle("Fold Test Ratios Canvas");  canv_fold_ratio         ->Write("canv_fold_ratio");
    canv_fold_ratio2       ->SetTitle("Fold Test Ratios v2 Canvas");  canv_fold_ratio2         ->Write("canv_fold_ratio2");
    canv_thy_ratio        ->SetTitle("NLO Thy Ratios Canvas");    canv_thy_ratio          ->Write("canv_thy_ratio");                                
    canv_thy_ratio2        ->SetTitle("NLO Thy Ratios Canvas");    canv_thy_ratio2          ->Write("canv_thy_ratio2");                                
    
    canv_covmat           ->SetTitle("Covariance Matrix Canvas");           canv_covmat        ->Write("canv_covmat");
    canv_absval_covmat    ->SetTitle("Abs Val. Covariance Matrix Canvas");  canv_absval_covmat ->Write("canv_covmatabsval");
    canv_pearson          ->SetTitle("Pearson Matrix Canvas");      canv_pearson               ->Write("canv_pearson");
    canv_unfmat           ->SetTitle("Unfolding Matrix Canvas");    canv_unfmat                ->Write("canv_unfmat");
    canv_mat_rebin        ->SetTitle("PY8 Response Matrix Canvas");  canv_mat_rebin             ->Write("canv_mat_rebin");
    canv_mat_percerrs     ->SetTitle("PY8 Response Matrix % Errors Canvas");  canv_mat_percerrs ->Write("canv_mat_percerrors");
    if(dokIterQA){ canv_3x3spectra->SetTitle("3x3 Unf. Spectra kIter QA Canvas");  canv_3x3spectra ->Write("canv_3x3spectra");}
    if(dokIterQA){ canv_3x3genratio->SetTitle("3x3 Gen Ratio kIter QA Canvas");    canv_3x3genratio->Write("canv_3x3genratio");}
    if(dokIterQA){ canv_3x3recratio->SetTitle("3x3 Rec Ratio kIter QA Canvas");    canv_3x3recratio->Write("canv_3x3recratio");}
    if(dokIterQA){ canv_3x3covmat->SetTitle("3x3 Covariance Matrix kIter QA Canvas");    canv_3x3covmat->Write("canv_3x3covmat");}
    if(dokIterQA){ canv_3x3covmatabsval->SetTitle("3x3 |Covariance| Matrix kIter QA Canvas");    canv_3x3covmatabsval->Write("canv_3x3covmatabsval");}
    if(dokIterQA){ canv_3x3pearson->SetTitle("3x3 PearsonMatrix kIter QA Canvas");    canv_3x3pearson->Write("canv_3x3pearson");}
    if(dokIterQA) canv_kIterRatio->Write();

  }
  
  std::cout<<"closing input/output root files..."<<std::endl;
  
  fout->Close();
  if(drawPDFs){
    makeCombinedPlots(outRootFile, canvForPrint, outPdfFile,applyNPCorrs);  
    canvForPrint->Print(close_outPdfFile.c_str());          }
  fpp_MC->Close();
  fpp_Data->Close();
  return 0;
  
} // end unfoldDataSpectra



//  steering ---------------------------------------------------------------------------------
int main(int argc, char* argv[]){  
  int rStatus = -1;  
  if( argc!=9 ){
    std::cout<<"argc!=9 error, argc=="<<argc<<std::endl;
    std::cout<<"do ./bayesUnfoldDataSpectra_wNLO.exe <targDataDir> <targMCDir> <targMCName> <baseOutputName> <applyNPCorrs> <doJetID> <useSimpleBins> <kIterInput>"<<std::endl;
    return rStatus;  }  
  rStatus=bayesUnfoldDataSpectra_wNLO(  (const std::string)argv[1], (const std::string)argv[2], (const std::string)argv[3], (const std::string)argv[4], 
					(int)std::atoi(argv[5]),
					(int)std::atoi(argv[6]) , (int)std::atoi(argv[7]), (const int)std::atoi(argv[8])  );   
  std::cout<<std::endl<<"done!"<<std::endl<<" return status: "<<rStatus<<std::endl<<std::endl;
  return rStatus;
}




    
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
    //legend_fold->AddEntry(hfold_truth, "PY8 Fold(Truth)"        , "lp");    
    //legend_fold->Draw();
    //
    //canvForPrint->Print(outPdfFile.c_str());      
    //    // PY8 gen/meas/unf spectra hists ---------------
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
    //    hrec_sameside_rebin->SetTitle("PY8 Spectra Only");
    //    
    //    hrec_sameside_rebin->DrawClone("P E");           
    //    hgen_rebin->DrawClone("P E SAME");                 
    //    hunf_sameside->DrawClone("P E SAME");               
    //    //hfold_sameside_clone->DrawClone("P E SAME");               
    //    
    //    TLegend* legend_in2 = new TLegend( 0.6,0.7,0.9,0.9 );
    //    
    //    legend_in2->AddEntry(hunf_sameside,        "PY8 Unf." , "lp");
    //    legend_in2->AddEntry(hgen_rebin,          "PY8 Truth" , "lp");
    //    legend_in2->AddEntry(hrec_sameside_rebin, "PY8 Meas." , "lp");
    //    //legend_in2->AddEntry(hfold_sameside_clone, "PY8 Fold.", "lp");    
    //    
    //    legend_in2->Draw();
    //    
    //    canvForPrint->Print(outPdfFile.c_str());      

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
