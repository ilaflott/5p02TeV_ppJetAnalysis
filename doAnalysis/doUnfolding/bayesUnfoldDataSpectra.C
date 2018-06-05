#include "unfoldSpectra.h"

//other settings
const bool drawPDFs=true; 
const bool debugMode=false;
const bool drawRespMatrix=true;

const float etaBinWidth=4.;//for |y|<2.0

// CODE --------------------------------------------------
int bayesUnfoldDataSpectra( std::string inFile_Data_dir , std::string inFile_MC_dir , std::string baseName , 
			    const bool doJetID=true , const bool useSimpBins=false){
  
  
  
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
  inFile_MC_dir  =SCRATCH_BASE+inFile_MC_dir;
  
  std::size_t radPos=inFile_MC_dir.find("_ak")+3;  
  const std::string radius="R"+inFile_MC_dir.substr( radPos,1 );
  const std::string RandEtaRange="_"+radius+"_20_eta_20";
  const std::string RandEtaRange_plotTitle=" "+radius+" 20eta20";// R4 20eta20
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;
  if(debugMode)std::cout<<"RandEtaRange string is = "<<RandEtaRange<<std::endl;  
  
  std::size_t jetTypePos=radPos+1, jetsPos=inFile_MC_dir.find("Jets");
  const std::string jetType=inFile_MC_dir.substr( jetTypePos,(jetsPos-jetTypePos) );//PFJets
  const std::string fullJetType="ak"+inFile_MC_dir.substr( radPos,1 )+jetType;//"ak4PFJets"
  if(debugMode)std::cout<<"jetType string is = "<<jetType<<std::endl;
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  
  
  
  // INFILE NAME(S) -----------
  const std::string inFile_MC_name="/Py8_CUETP8M1_QCDjetAllPtBins_"+fullJetType+"-allFiles.root";
  const std::string inFile_Data_name="/HighPtJetTrig_"+fullJetType+"-allFiles.root";
  
  
  // OUTPUT FILE, NAME(S) -----------
  std::string outFileName=unfoldDataSpectra_outdir+fullJetType;
  if(doJetID)outFileName+="_wjtID";
  if(!useSimpBins)outFileName+="_anabins";
  else outFileName+="_simpbins";
  outFileName+="_"+baseName;
  
  
  std::string outRespMatPdfFile =  outFileName+"_respMat.pdf";
  std::string outBayesPdfFile =  outFileName+".pdf";
  std::string outRootFile     =  outFileName+".root";  
  
  checkNRenameFiles ( (const std::string) outFileName, 
		      &outRespMatPdfFile, 
		      &outBayesPdfFile, 
		      &outRootFile);
  
  
  
  if(debugMode)std::cout<<"opening output file: "<<outRootFile<<std::endl;  
  TFile* fout = new TFile(outRootFile.c_str(),"RECREATE");   

  //for output
  fout->cd();  
  
  
  // ppMC input histos -------------------------
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from MC file"<<std::endl; 
  std::cout<<"input MC dir : "<< (inFile_MC_dir)  <<std::endl; 
  std::cout<<"MC file name : "<< (inFile_MC_name)<<std::endl;   std::cout<<std::endl<<std::endl;  
  TFile *fpp_MC = TFile::Open( (inFile_MC_dir+inFile_MC_name).c_str());
  
  //for output
  fout->cd();  

  // ---------- open MC "response" matrix
  std::string TH2_title="hpp_matrix";//+RandEtaRange;
  if(doJetID)TH2_title+="_wJetID";
  TH2_title+=RandEtaRange;
  
  //get the response matrix made by readforests
  TH2F* hmat = (TH2F*)fpp_MC->Get( TH2_title.c_str() );
  hmat->Write(TH2_title.c_str());
  if(debugMode)hmat->Print("base");
  
  //TH2_title+="_divByetabin";
  //hmat->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  //hmat->Write( TH2_title.c_str());
  //if(debugMode)hmat->Print("base");
  
  // rebinned matrix ---------------
  TH2_title+="_clone";
  TH2F* hmat_rebin = (TH2F*)hmat->Clone( (TH2_title).c_str() );
  hmat_rebin->Write(TH2_title.c_str());
  if(debugMode)hmat_rebin->Print("base"); 
  
  TH2_title+="_rebins";
  hmat_rebin=(TH2F*) reBinTH2(hmat_rebin, (TH2_title).c_str(), 
			      (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
			      (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );  
  hmat_rebin->Write(TH2_title.c_str());
  if(debugMode)hmat_rebin->Print("base"); 
  
  //TH2_title+="_normbinwidth";
  //divideBinWidth_TH2(hmat_rebin);
  //hmat_rebin->Write(TH2_title.c_str());
  //if(debugMode)hmat_rebin->Print("base"); 
  
  if(clearOverUnderflows){
    TH2_title+="_noOverUnderFlows";
    TH2clearOverUnderflows((TH2F*)hmat_rebin);
    hmat_rebin->Write(TH2_title.c_str());
    if(debugMode)hmat_rebin->Print("base");  }
  
  
  // error and %error for response matrix ---------------
  std::string errTH2_title="hmat_errors_rebins";
  TH2F* hmat_errors=makeRespMatrixErrors( (TH2F*) hmat,
					  (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
					  (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
  hmat_errors->Write(errTH2_title.c_str());
  if(debugMode)hmat_errors->Print("base");
  
  //errTH2_title+="_normbinwidth";
  //divideBinWidth_TH2(hmat_errors);
  //hmat_errors->Write(errTH2_title.c_str());
  //if(debugMode)hmat_errors->Print("base");
  
  if(clearOverUnderflows){
    errTH2_title+="_noOverUnderFlows";
    TH2clearOverUnderflows((TH2F*)hmat_errors);
    hmat_errors->Write(errTH2_title.c_str());
    if(debugMode)hmat_errors->Print("base");  }
  
  
  // give response matrix the correct errors
  setRespMatrixErrs( (TH2F*)hmat_rebin, (TH2F*) hmat_errors , (bool)zeroBins);

  TH2_title+="_wseterrs";
  hmat_rebin->Write(TH2_title.c_str());
  if(debugMode)hmat_rebin->Print("base");
  
  




  
  // ---------- reco, measured same-side spectra used to create response matrix, and for "sameside" unfolding test
  std::string histTitle2="hpp_reco";//    hpp_reco_wJetID_R4_20_eta_20
  if(doJetID)histTitle2+="_wJetID";
  histTitle2+=RandEtaRange;
  
  TH1F*  hrec_sameside = (TH1F*)fpp_MC->Get( histTitle2.c_str() ); 
  hrec_sameside->Write(histTitle2.c_str());
  if(debugMode)hrec_sameside->Print("base");
  
  histTitle2+="_clone";
  TH1F *hrec_sameside_rebin = (TH1F*)hrec_sameside->Clone( (histTitle2).c_str() );
  hrec_sameside_rebin->Write( histTitle2.c_str() );
  if(debugMode)hrec_sameside_rebin->Print("base");


  std::cout<<"rebinning hrec_sameside..."<<std::endl;
  histTitle2+="_rebins";
  hrec_sameside_rebin = (TH1F*)hrec_sameside_rebin->Rebin( nbins_pt_reco, (histTitle2).c_str() , boundaries_pt_reco);
  hrec_sameside_rebin->Write( histTitle2.c_str() );   
  if(debugMode)hrec_sameside_rebin->Print("base");  
  
  if(clearOverUnderflows){
    histTitle2+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hrec_sameside_rebin);
    hrec_sameside_rebin->Write( histTitle2.c_str() );
    if(debugMode)hrec_sameside_rebin->Print("base");    }
  
  //cosmetics
  hrec_sameside_rebin->SetMarkerStyle(kOpenSquare);
  hrec_sameside_rebin->SetMarkerColor(kBlue-3);     
  hrec_sameside_rebin->SetMarkerSize(1.02);     
  

  
  
  // ---------- gen, MC truth spectra
  std::string genHistTitle="hpp_gen";
  if(doJetID)genHistTitle+="_wJetID";
  genHistTitle+=RandEtaRange;
  
  TH1F* hgen = (TH1F*)fpp_MC->Get( genHistTitle.c_str() );
  hgen->Write(genHistTitle.c_str());
  if(debugMode)hgen->Print("base");    
  
  genHistTitle+="_clone";
  TH1F* hgen_rebin = (TH1F*)hgen->Clone( (genHistTitle).c_str() );
  hgen_rebin->Write(genHistTitle.c_str());
  if(debugMode)hgen_rebin->Print("base");
  
  genHistTitle+="_rebins";
  hgen_rebin = (TH1F*)hgen_rebin->Rebin(nbins_pt_gen, (genHistTitle).c_str() , boundaries_pt_gen);
  hgen_rebin->Write(genHistTitle.c_str());
  if(debugMode)hgen_rebin->Print("base"); 
  
  if(clearOverUnderflows){
    genHistTitle+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hgen_rebin);
    hgen_rebin->Write(genHistTitle.c_str());
    if(debugMode)hgen_rebin->Print("base");    }

  //cosmetics
  hgen_rebin->SetMarkerStyle(kOpenStar);
  hgen_rebin->SetMarkerColor(kMagenta);
  hgen_rebin->SetMarkerSize(1.02);     	



  
  // reco/gen response hist

  std::string histTitle3="hpp_rec_sameside_response_rebin"+RandEtaRange;
  TH1F* hrec_sameside_resp_rebin=NULL;  
  std::string histTitle3gen="hpp_gen_response_rebin"+RandEtaRange;
  TH1F* hgen_resp_rebin=NULL;

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
    //  hrec_sameside_resp_rebin = (TH1F*)hmat->TH2::ProjectionX( 
    //							    (histTitle3+"_profile").c_str() , 
    //							    (int)(boundaries_pt_gen[0]+1.), 
    //							    ((int)boundaries_pt_gen[nbins_pt_gen]), 	"g");
    //  hrec_sameside_resp_rebin=(TH1F*)hrec_sameside_resp_rebin->Rebin( nbins_pt_reco, (histTitle3+"_profile_rebin").c_str() , boundaries_pt_reco);
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
    //  hgen_resp_rebin = (TH1F*)hmat->TH2::ProjectionY( 
    //						   (histTitle3gen+"_profile").c_str() , 
    //						   (int)(boundaries_pt_reco[0]+1.), 
    //						   ((int)boundaries_pt_reco[nbins_pt_reco]), 	"g");
    //  hgen_resp_rebin=(TH1F*)hgen_resp_rebin->Rebin( nbins_pt_gen, (histTitle3gen+"_profile_rebin").c_str() , boundaries_pt_gen);
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
    hrec_sameside_resp_rebin = (TH1F*)hrec_sameside_rebin->Clone( (histTitle3+"_clone").c_str() );      
    hrec_sameside_resp_rebin->Write( (histTitle3+"_clone").c_str() );    
    
    hgen_resp_rebin = (TH1F*)hgen_rebin->Clone( (histTitle3gen+"_clone").c_str() );      
    hgen_resp_rebin->Write( (histTitle3gen+"_clone").c_str() );       
    //}
    
  }
  else{
    hrec_sameside_resp_rebin = new TH1F( (histTitle3+"_empty").c_str() ,"", 
					 nbins_pt_reco, boundaries_pt_reco);  
    
    hgen_resp_rebin = new TH1F( (histTitle3gen+"_empty").c_str() ,"", 
				nbins_pt_gen, boundaries_pt_gen);     }


  hrec_sameside_resp_rebin->SetMarkerStyle(kOpenSquare);
  hrec_sameside_resp_rebin->SetMarkerColor(kBlue-3);     
  hrec_sameside_resp_rebin->SetMarkerSize(1.02);       
  hrec_sameside_resp_rebin->Write((histTitle3+"_final").c_str());
  if(debugMode)hrec_sameside_resp_rebin->Print(" base");  
  
  hgen_resp_rebin->SetMarkerStyle(kFullStar);
  hgen_resp_rebin->SetMarkerColor(kMagenta);
  hgen_resp_rebin->SetMarkerSize(1.02);     	
  hgen_resp_rebin->Write((histTitle3gen+"_final").c_str());
  if(debugMode)hgen_resp_rebin->Print("base");    
  

  
  
  
  // ppData input histos -------------------------
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from Data file"<<std::endl; 
  std::cout<<"input Data dir : "<< (inFile_Data_dir)  <<std::endl; 
  std::cout<<"Data file name : "<< (inFile_Data_name)<<std::endl;   std::cout<<std::endl;//<<std::endl;  
  TFile *fpp_Data = TFile::Open( (inFile_Data_dir+inFile_Data_name).c_str());

  //for output
  fout->cd();    
  
  
  
  // ---------- reco, measured spectra to unfold
  std::string histTitle="hJetQA";//hJetQA_1wJetID_jtpt
  if(doJetID)histTitle+="_1wJetID";
  else histTitle+="_0wJetID";
  histTitle+="_jtpt";
  //histTitle+=RandEtaRange;
  
  TH1F*  hrec = (TH1F*)fpp_Data->Get( histTitle.c_str() ); 
  hrec->Write();
  if(debugMode)hrec->Print("base");
  
  histTitle+="_divBylumietabin";
  float effIntgrtdLumi=computeEffLumi(fpp_Data);
  hrec->Scale(1./effIntgrtdLumi); // lumi
  hrec->Write( (histTitle).c_str() );
  if(debugMode)hrec->Print("base");
  
  //assert(false);
  
  histTitle+="_clone";
  TH1F *hrec_rebin = (TH1F*)hrec->Clone( (histTitle).c_str() );
  hrec_rebin->Write(histTitle.c_str());
  if(debugMode)hrec_rebin->Print("base");
  
  std::cout<<"rebinning hrec..."<<std::endl;
  histTitle+="_rebins";
  hrec_rebin = (TH1F*)hrec_rebin->Rebin( nbins_pt_reco, (histTitle).c_str() , boundaries_pt_reco);
  hrec_rebin->Write(histTitle.c_str());   
  if(debugMode)hrec_rebin->Print("base");  
  
  
  if(clearOverUnderflows){
    histTitle+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hrec_rebin);
    hrec_rebin->Write(histTitle.c_str());
    if(debugMode)hrec_rebin->Print("base");  
  }
  
  //cosmetics
  hrec_rebin->SetMarkerStyle(kOpenCircle);
  hrec_rebin->SetMarkerColor(kBlue);     
  hrec_rebin->SetMarkerSize(1.02);     

  
  
  
  
  

  
  
  
  
  
  // Bayesian unfolding -------------------------   
  // RooUnfoldResponse, Bayes --------------------
  std::cout<<"calling RooUnfoldResponse "<<std::endl;
  
  RooUnfoldResponse roo_resp( hrec_sameside_resp_rebin, hgen_resp_rebin, hmat_rebin, ("Response_matrix"+RandEtaRange).c_str()) ;
  //RooUnfoldResponse roo_resp( 0, 0, hmat_rebin, ("Response_matrix"+RandEtaRange).c_str()) ;
  roo_resp.UseOverflow(doOverUnderflows);    
  roo_resp.Write();

  TH1F* hfak=  (TH1F*) roo_resp.Hfakes() ;
  hfak->Print("base");
  hfak->SetMarkerColor(kGreen);
  //hrec_sameside_resp_rebin->Print("base");

  //assert(false);

  
  
  std::cout<<"calling RooUnfoldBayes..."<<std::endl;
  RooUnfoldBayes unf_bayes( &roo_resp, hrec_rebin, kIter );
  unf_bayes.SetVerbose(2);
  if(doToyErrs){
    std::cout<<"using toy errors, suppressing text output"<<std::endl;
    unf_bayes.SetNToys(10000);
    unf_bayes.SetVerbose(1);
  }
  std::cout<<"RooUnfoldBayes Overflow Status: " << unf_bayes.Overflow()<<std::endl;
  
  TH1F *hunf = (TH1F*)unf_bayes.Hreco(errorTreatment);     std::cout<<std::endl; 
  hunf->SetName("ppData_BayesUnf_Spectra");
  hunf->SetTitle("Unf. Data");
  if(debugMode)hunf->Print("base");

  //cosmetics
  hunf->SetMarkerStyle(kOpenCircle);
  hunf->SetMarkerColor(kRed);
  hunf->SetMarkerSize(1.02);     
  
  std::cout<<"reFolding unfolded histogram!!"<<std::endl;
  TH1F* hfold=(TH1F*)roo_resp.ApplyToTruth(hunf);        

  hfold->SetName("ppData_BayesFold_Spectra");
  hfold->SetTitle("Fold. Data");
  //hfold->Scale(hrec_rebin->Integral()/hfold->Integral());
  if(debugMode)hfold->Print("base");


  //cosmetics
  hfold->SetMarkerStyle(kOpenCircle);
  hfold->SetMarkerColor(kGreen);
  hfold->SetMarkerSize(1.02);     



  std::cout<<"calling Ereco... getting Cov Mat"<<std::endl;
  TMatrixD covmat = unf_bayes.Ereco(errorTreatment);
  
  std::cout<<"converting TMatrixD covmat to TH2D* covmat_TH2"<<std::endl;
  TH2D* covmat_TH2= new TH2D(covmat);  //this gives the TH2 a name of "TMatrixDBase" by default.
  covmat_TH2->SetName("Bayes_covarianceMatrix_TH2");
  if(debugMode)covmat_TH2->Print("base");
  covmat_TH2->Write();
  
  std::cout<<"calculating pearson coefficients"<<std::endl;
  TH2D* PearsonBayes = (TH2D*)CalculatePearsonCoefficients( &covmat, false,"Bayes_pearson");
  if(debugMode)PearsonBayes->Print("base");
  //PearsonBayes->SetName("covarianceMatrix_TH2");
  PearsonBayes->Write();

  std::cout<<"getting unfolding matrix"<<std::endl;
  TH2D* unfmat_TH2=new TH2D(unf_bayes.UnfoldingMatrix());
  unfmat_TH2->SetName("Bayes_unfoldingMatrix_TH2");  
  if(debugMode)unfmat_TH2->Print("base");
  unfmat_TH2->Write();
  
  //std::cout<<"reg param="<<unf_bayes.GetRegParm()<<std::endl;
  
  // RATIOS WITH OPP SIDE GEN
  TH1F *hgen_rebin_ratiobin=(TH1F*)hgen_rebin->Clone("ppData_Gen_Ratio_denom");
  hgen_rebin_ratiobin=(TH1F*)hgen_rebin_ratiobin->Rebin(nbins_pt_reco,"ppData_Gen_Ratio_denom_rebin",boundaries_pt_reco);
  if(debugMode)hgen_rebin_ratiobin->Print("base");
  
  TH1F *h_genratio_oppunf = (TH1F*)hunf->Clone( "ppData_Gen_Ratio_OppUnf" );
  h_genratio_oppunf->SetTitle( "Unf. Data/GEN Py8" );
  h_genratio_oppunf->Divide(hgen_rebin);
  if(debugMode)h_genratio_oppunf->Print("base");
  
  TH1F *h_genratio_oppfold = (TH1F*)hfold->Clone( "ppData_Gen_Ratio_OppFold" );
  h_genratio_oppfold->SetTitle( "Fold. Data/GEN Py8" );
  h_genratio_oppfold->Divide(hgen_rebin_ratiobin);
  if(debugMode)h_genratio_oppfold->Print("base");
  
  TH1F *h_genratio_oppmeas = (TH1F*)hrec_rebin->Clone( "ppData_Gen_Ratio_Meas" );
  h_genratio_oppmeas->SetTitle( "Meas. Data/GEN Py8" );
  h_genratio_oppmeas->Divide(hgen_rebin_ratiobin);
  if(debugMode)h_genratio_oppmeas->Print("base");
  
  
  // RATIOS WITH OPP SIDE MEAS
  TH1F *h_recratio_oppunf = (TH1F*)hunf->Clone( "ppData_Meas_Ratio_OppUnf" );
  h_recratio_oppunf=(TH1F*)h_recratio_oppunf->Rebin(nbins_pt_reco, "ppData_Meas_Ratio_OppUnf_rebin" , boundaries_pt_reco);
  h_recratio_oppunf->SetTitle( "Unf. Data/RECO Data" );
  h_recratio_oppunf->Divide(hrec_rebin);
  if(debugMode)h_recratio_oppunf->Print("base");
  
  TH1F *h_recratio_oppfold = (TH1F*)hfold->Clone( "ppData_Meas_Ratio_OppFold" );
  h_recratio_oppfold=(TH1F*)h_recratio_oppfold->Rebin(nbins_pt_reco, "ppData_Meas_Ratio_OppFold_rebin" , boundaries_pt_reco);
  h_recratio_oppfold->SetTitle( "Fold. Data/RECO Data" );
  h_recratio_oppfold->Divide(hrec_rebin);
  if(debugMode)h_recratio_oppfold->Print("base");
  
  
  
  
  
  //std::cout<<"calling sameside RooUnfoldBayes..."<<std::endl;
  RooUnfoldBayes unf_sameside_bayes( &roo_resp, hrec_sameside_rebin, kIter );
  //RooUnfoldBayes unf_sameside_bayes( &roo_resp, hrec_sameside_resp_rebin, kIter );
  
  TH1F *hunf_sameside = (TH1F*)unf_sameside_bayes.Hreco(errorTreatment);    std::cout<<std::endl; 
  hunf_sameside->SetName("ppMC_BayesUnf_sameSideSpectra");
  hunf_sameside->SetTitle("Unf. Py8");
  if(debugMode)hunf_sameside->Print("base");
  //cosmetics
  hunf_sameside->SetMarkerStyle(kOpenSquare);
  hunf_sameside->SetMarkerColor(kRed-3);
  hunf_sameside->SetMarkerSize(1.02);     
  
  
  std::cout<<"reFolding unfolded sameside histogram!!"<<std::endl;
  TH1F* hfold_sameside=(TH1F*)roo_resp.ApplyToTruth(hunf_sameside);    
  hfold_sameside->SetName("ppMC_BayesFold_sameSide_Spectra");
  hfold_sameside->SetTitle("Fold. Py8");
  if(debugMode)hfold_sameside->Print("base");
  
  
  //cosmetics
  hfold_sameside->SetMarkerStyle(kOpenSquare);
  hfold_sameside->SetMarkerColor(kGreen-3);
  hfold_sameside->SetMarkerSize(1.02);     
  
  //hfold_sameside->Scale(hrec_sameside_rebin->Integral()/hfold_sameside->Integral());
  if(debugMode)hfold_sameside->Print("base");
  
  
  
  // SAMESIDE RATIO W/ OPPSIDE MEAS
  TH1F *h_recratio_ssunf = (TH1F*)hunf_sameside->Clone( "ppMC_Meas_Ratio_SSUnf" );
  h_recratio_ssunf=(TH1F*)h_recratio_ssunf->Rebin(nbins_pt_reco, "ppMC_Meas_Ratio_SSUnf_rebin" , boundaries_pt_reco);
  h_recratio_ssunf->SetTitle( "Unf. Py8/RECO Data" );
  h_recratio_ssunf->Divide(hrec_rebin);
  if(debugMode)h_recratio_ssunf->Print("base");

  TH1F *h_recratio_ssmeas = (TH1F*)hrec_sameside_rebin->Clone( "ppMC_Meas_Ratio_SSMeas" );
  h_recratio_ssmeas=(TH1F*)h_recratio_ssmeas->Rebin(nbins_pt_reco, "ppMC_Meas_Ratio_SSMeas_rebin" , boundaries_pt_reco);
  h_recratio_ssmeas->SetTitle( "RECO Py8/RECO Data" );
  h_recratio_ssmeas->Divide(hrec_rebin);
  if(debugMode)h_recratio_ssmeas->Print("base");
  
  
  TH1F *h_recratio_ssfold = (TH1F*)hfold_sameside->Clone( "ppMC_Meas_Ratio_SSFold" );
  h_recratio_ssfold=(TH1F*)h_recratio_ssfold->Rebin(nbins_pt_reco, "ppMC_Meas_Ratio_SSFold_rebin" , boundaries_pt_reco);
  h_recratio_ssfold->SetTitle( "Fold. Py8/RECO Data" );
  h_recratio_ssfold->Divide(hrec_rebin);
  if(debugMode)h_recratio_ssfold->Print("base");
  
  
  TH1F *h_recratio_ssgen = (TH1F*)hgen_rebin->Clone( "ppMC_Meas_Ratio_SSTruth" );    
  h_recratio_ssgen=(TH1F*)h_recratio_ssgen->Rebin(nbins_pt_reco, "ppMC_Meas_Ratio_SSTruth_rebin" , boundaries_pt_reco);
  h_recratio_ssgen->SetTitle( "GEN Py8/RECO Data" );
  h_recratio_ssgen->Divide(hrec_rebin);
  if(debugMode)h_recratio_ssgen->Print("base");
  
  
  // SAMESIDE RATIO W/ OPPSIDE GEN
  TH1F *h_genratio_ssunf = (TH1F*)hunf_sameside->Clone( "ppMC_Gen_Ratio_SSUnf" );
  h_genratio_ssunf->SetTitle( "Unf. Py8/GEN Py8" );
  h_genratio_ssunf->Divide(hgen_rebin);
  if(debugMode)h_genratio_ssunf->Print("base");
  
  
  TH1F *h_genratio_ssfold = (TH1F*)hfold_sameside->Clone( "ppMC_Gen_Ratio_SSFold" );
  h_genratio_ssfold->SetTitle( "Fold. Py8/GEN Py8" );
  h_genratio_ssfold->Divide(hgen_rebin_ratiobin);
  if(debugMode)h_genratio_ssfold->Print("base");
  
  TH1F *h_genratio_ssmeas = (TH1F*)hrec_sameside_rebin->Clone( "ppMC_Gen_Ratio4_SSMeas" );
  h_genratio_ssmeas->SetTitle( "RECO Py8/GEN Py8" );
  h_genratio_ssmeas->Divide(hgen_rebin_ratiobin);
  if(debugMode)h_genratio_ssmeas->Print("base");
  
  




  // ----- PUT BIN WIDTH(s) DIVISIONS + NORMALIZING HERE ----- //

  // -- MC RECO -- //  
  histTitle2+="_normbinwidth";
  hrec_sameside_rebin->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hrec_sameside_rebin); 
  hrec_sameside_rebin->Write( histTitle2.c_str() );
  if(debugMode)hrec_sameside_rebin->Print("base");  
  
  // -- MC GEN -- //  
  genHistTitle+="_normbinwidth";
  hgen_rebin->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hgen_rebin);
  hgen_rebin->Write(genHistTitle.c_str());
  if(debugMode)hgen_rebin->Print("base");  

  // -- MC GEN RESP -- //  
  hgen_resp_rebin->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hgen_resp_rebin);
  hgen_resp_rebin->Write(genHistTitle.c_str());
  if(debugMode)hgen_resp_rebin->Print("base");  

  // -- MC RESP FAKES -- //  
  hfak->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hfak);
  hfak->Write();
  if(debugMode)hfak->Print("base");  
  
  // -- MC REC RESP -- //  
  hrec_sameside_resp_rebin->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hrec_sameside_resp_rebin);
  hrec_sameside_resp_rebin->Write();
  if(debugMode)hrec_sameside_resp_rebin->Print("base");  
  
  // -- DATA RECO -- //  
  histTitle+="_normbinwidth";
  hrec_rebin->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hrec_rebin); 
  hrec_rebin->Write(histTitle.c_str());
  if(debugMode)hrec_rebin->Print("base");  
  
  // -- DATA UNF -- //
  hunf->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hunf); 
  hunf->Write();
  if(debugMode)hunf->Print("base");  
  
  


  // ----- END BIN WIDTH DIVISIONS + NORMALIZING ----- //
  

  
  
  //  drawPDFS -------------------------------------------------- 
  if(!drawPDFs)std::cout<<std::endl<<"NOT drawing PDFs for Bayesian Unfolding!"<<std::endl<<std::endl;
  else{ std::cout<<std::endl<<"drawing output PDFs for Bayesian Unfolding..."<<std::endl;
    std::string outPdfFile=outBayesPdfFile;
    std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";
    
    // canvas for printing
    TCanvas* canvForPrint=new TCanvas("printCanvas","for printing",           1400,1000);    
    canvForPrint->cd();
    //canvForPrint->SetGridx(1);    //canvForPrint->SetGridy(1);
    canvForPrint->Print(open_outPdfFile.c_str()); 
    

    //misc printing prep
    TLine* theLineAtOne_gen= new TLine( boundaries_pt_gen_mat[0]  
					,1.   
					,(boundaries_pt_gen_mat[nbins_pt_reco_mat])
					,1.);
    if(debugMode)std::cout<<"theLineAtOne_gen starts at "<<boundaries_pt_gen_mat[0]  <<std::endl;
    if(debugMode)std::cout<<"theLineAtOne_gen ends at   "<<(boundaries_pt_gen_mat[nbins_pt_reco_mat])  <<std::endl;
    theLineAtOne_gen->SetLineWidth(1);
    theLineAtOne_gen->SetLineStyle(2);
    theLineAtOne_gen->SetLineColor(36);

    
    TLine* theLineAtp9_gen= new TLine( boundaries_pt_gen_mat[0]  
				       ,0.9   
				       ,(boundaries_pt_gen_mat[nbins_pt_reco_mat])
				       ,0.9);
    if(debugMode)std::cout<<"theLineAtp9_gen starts at "<<boundaries_pt_gen_mat[0]  <<std::endl;
    if(debugMode)std::cout<<"theLineAtp9_gen ends at   "<<(boundaries_pt_gen_mat[nbins_pt_reco_mat])  <<std::endl;
    theLineAtp9_gen->SetLineWidth(1);
    theLineAtp9_gen->SetLineStyle(2);
    theLineAtp9_gen->SetLineColor(36);
    
    TLine* theLineAt1p1_gen= new TLine( boundaries_pt_gen_mat[0]  
					,1.1   
					,(boundaries_pt_gen_mat[nbins_pt_reco_mat])
					,1.1);
    if(debugMode)std::cout<<"theLineAt1p1_gen starts at "<<boundaries_pt_gen_mat[0]  <<std::endl;
    if(debugMode)std::cout<<"theLineAt1p1_gen ends at   "<<(boundaries_pt_gen_mat[nbins_pt_reco_mat])  <<std::endl;
    theLineAt1p1_gen->SetLineWidth(1);
    theLineAt1p1_gen->SetLineStyle(2);
    theLineAt1p1_gen->SetLineColor(36);
    
    
    TLine* theLineAtOne_reco= new TLine( boundaries_pt_reco_mat[0]  
					 ,1.   
					 ,(boundaries_pt_reco_mat[nbins_pt_reco_mat])
					 ,1.);    
    if(debugMode)std::cout<<"theLineAtOne_reco starts at "<<boundaries_pt_reco_mat[0]  <<std::endl;
    if(debugMode)std::cout<<"theLineAtOne_reco ends at   "<<boundaries_pt_reco_mat[nbins_pt_reco_mat]  <<std::endl;
    theLineAtOne_reco->SetLineWidth(1);
    theLineAtOne_reco->SetLineStyle(2);
    theLineAtOne_reco->SetLineColor(36);
            

    //TH1F* hfold_clone=(TH1F*)hfold->Clone("hfold_clone");
    
    // gen/meas/unf spectra hists ---------------    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(1);    
    
    setupSpectraHist(hgen_rebin          , useSimpBins);    
    setupSpectraHist(hrec_rebin	       , useSimpBins);
    setupSpectraHist(hrec_sameside_rebin , useSimpBins);
    //setupSpectraHist(hfold_clone	       , useSimpBins);
    setupSpectraHist(hunf		       , useSimpBins);    
    
    hgen_rebin->SetTitle("MC, Data, and Bayes Unf. Spectra");

    
    hgen_rebin->Draw("P E");                 
    hrec_rebin->Draw("P E SAME");           
    //hfold_clone->Draw("P E SAME");           
    hunf->Draw("P E SAME");               
    //hrec_sameside_rebin->Draw("P E SAME");           
    
    TLegend* legend_in1 = new TLegend( 0.7,0.7,0.9,0.9 );
    
    //legend_in1->AddEntry(hrec_sameside_rebin, "MC Meas."  ,  "lp");
    legend_in1->AddEntry(hunf,                "Data Unf." ,  "lp");
    legend_in1->AddEntry(hrec_rebin,          "Data Meas." , "lp");	
    legend_in1->AddEntry(hgen_rebin,          "MC Truth"   , "lp");
    //legend_in1->AddEntry(hfold_clone,          "Data Fold(Unf.)", "lp");    
    
    legend_in1->Draw();
    
    canvForPrint->Print(outPdfFile.c_str());      
    
    
    
    // gen/reco response spectra hists ---------------    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(1);    
    
    setupSpectraHist(hrec_sameside_resp_rebin , useSimpBins);
    setupSpectraHist(hgen_resp_rebin          , useSimpBins);
    setupSpectraHist(hfak                     , useSimpBins);
    
    //histograms given to roo_resp + the fakes histo made from these spectra + response matrix projection
    hrec_sameside_resp_rebin->SetTitle("MC Only, Response Spectra and Fakes");
    //hrec_sameside_resp_rebin->GetYaxis()->SetTitle("N_{Jets}/L_{int}");    
    
    hrec_sameside_resp_rebin->Draw("P E");           
    hgen_resp_rebin->Draw("P E SAME");                 
    hfak->Draw("P E SAME");           

    TLegend* legend_resp = new TLegend( 0.6,0.7,0.9,0.9 );
    legend_resp->AddEntry(hfak,          "MC Response Matrix Fakes" , "lp");    
    legend_resp->AddEntry(hgen_resp_rebin,          "MC Truth Response" , "lp");
    legend_resp->AddEntry(hrec_sameside_resp_rebin, "MC Meas. Response" , "lp");
    
    legend_resp->Draw();
    
    canvForPrint->Print(outPdfFile.c_str());      
    
    
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
//    hrec_sameside_rebin->Draw("P E");           
//    hgen_rebin->Draw("P E SAME");                 
//    hunf_sameside->Draw("P E SAME");               
//    //hfold_sameside_clone->Draw("P E SAME");               
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
    


    // ratios w MC Truth------------    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(0);
    
    setupRatioHist(h_genratio_oppunf, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_genratio_oppmeas, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_genratio_oppfold, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_genratio_ssmeas, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);    
    
    h_genratio_oppunf->SetTitle( "Bayesian, Ratios w/ MC Truth Spectra" );
    h_genratio_oppunf->GetYaxis()->SetTitle("Ratio w/ MC Truth");

    h_genratio_oppunf->Draw("P E");            //data unf/mc truth
    h_genratio_oppmeas->Draw("P E SAME");      // data meas/mc truth
    h_genratio_ssmeas->Draw("P E SAME");       // mc meas/mc truth
    //h_genratio_oppfold->Draw("P E SAME");
    
    TLegend* legend2 = new TLegend( 0.1,0.8,0.2,0.9 );
    legend2->AddEntry(h_genratio_oppunf,  "Data Unf." ,  "lp");
    legend2->AddEntry(h_genratio_oppmeas, "Data Meas." , "lp"); 
    legend2->AddEntry(h_genratio_ssmeas, "MC Meas.", "lp");
    //legend2->AddEntry(h_genratio_oppfold, "Data Fold(Unf.)", "lp");
    
    legend2->Draw();
    
    theLineAtOne_gen->Draw();
    theLineAtp9_gen->Draw();
    theLineAt1p1_gen->Draw();
    
    canvForPrint->Print(outPdfFile.c_str());




    // SS unfolding, ratio w/ MC truth -------------------
    //    canvForPrint->cd();
    //    if(!useSimpBins)canvForPrint->SetLogx(1);
    //    canvForPrint->SetLogy(0);
    //    
    //    setupRatioHist(h_genratio_ssunf, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    //    setupRatioHist(h_genratio_ssmeas, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    //    setupRatioHist(h_genratio_ssfold, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    //    
    //    h_genratio_ssunf->SetTitle( "MC Ratios w/ MC Truth Spectra" );
    //    
    //    h_genratio_ssunf->Draw("P E");
    //    h_genratio_ssmeas->Draw("P E SAME");
    //    //h_genratio_ssfold->Draw("P E SAME");
    //    
    //    TLegend* legend3 = new TLegend( 0.1,0.8,0.4,0.9 );
    //    
    //    legend3->AddEntry(h_genratio_ssmeas ,  "MC Meas./MC Truth" , "lp");
    //    legend3->AddEntry(h_genratio_ssunf  ,  "MC Unf./MC Truth",  "lp");
    //    //legend3->AddEntry(h_genratio_ssfold ,  "Fold. Py8/GEN Py8", "p");
    //    
    //    legend3->Draw();
    //    
    //    theLineAtOne_gen->Draw();
    //    
    //    canvForPrint->Print(outPdfFile.c_str());
    


    // ratio w/ data meas -------------------
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(0);      
    
    h_recratio_oppunf->SetTitle("Ratios w/ Data Meas. Spectra");
    h_recratio_oppunf->GetYaxis()->SetTitle("Ratio w/ Data Meas.");

    setupRatioHist(h_recratio_oppunf, useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);
    
    
    h_recratio_oppunf->Draw("P E");
    h_recratio_ssmeas->Draw("P E SAME");
    h_recratio_ssgen->Draw("P E SAME");
    ////h_recratio_ssunf->Draw("P E SAME");
    //h_recratio_oppfold->Draw("P E SAME");  
    //h_recratio_ssfold->Draw("P E SAME");
    
    TLegend* legend4 = new TLegend( 0.1,0.8,0.2,0.9 );
    legend4->AddEntry(h_recratio_oppunf,  "Data Unf.", "lp");
    legend4->AddEntry(h_recratio_ssmeas, "MC Meas.", "lp");
    legend4->AddEntry(h_recratio_ssgen,   "MC Truth", "lp");
    ////legend4->AddEntry(h_recratio_ssunf, NULL, "p");
    //legend4->AddEntry(h_recratio_oppfold, "Data Fold(Unf.)", "lp");
    //legend4->AddEntry(h_recratio_ssfold, "MC Fold(Unf.)", "p");
    legend4->Draw();
    
    theLineAtOne_reco->Draw();
    theLineAtp9_gen->Draw();
    theLineAt1p1_gen->Draw();
    
    canvForPrint->Print(outPdfFile.c_str());
    



    
    
    // covariance matrix ---------------      
    canvForPrint->cd();
    canvForPrint->SetLogx(0);
    canvForPrint->SetLogy(0);
    canvForPrint->SetLogz(1);
    
    //always use simp bins for covmat to avoid log scaling the x/y axes
    matStylePrint( (TH2*)covmat_TH2, "Bayesian Unfolding Covariance Matrix", canvForPrint, outPdfFile, true);
    
    // pearson matrix ---------------      
    canvForPrint->cd();
    canvForPrint->SetLogx(0);
    canvForPrint->SetLogy(0);
    canvForPrint->SetLogz(0);
    
    //always use simp bins for covmat to avoid log scaling the x/y axes
    matStylePrint( (TH2*)PearsonBayes, "Bayesian Unfolding Pearson Matrix", canvForPrint, outPdfFile, true);
    
    
    // unfolding matrix ---------------      
    canvForPrint->cd();
    canvForPrint->SetLogx(0);
    canvForPrint->SetLogy(0);
    canvForPrint->SetLogz(1);
      
    //always use simp bins for covmat to avoid log scaling the x/y axes
    matStylePrint( (TH2*)unfmat_TH2, "Bayesian Unfolding Matrix", canvForPrint, outPdfFile, true);



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
    
    


    // percent error matrix in binning of interest ---------------    
    TH2F* hmat_percenterrs= makeRespMatrixPercentErrs( (TH2F*) hmat_errors, (TH2F*) hmat_rebin,
						       (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
						       (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );		     
    hmat_percenterrs->SetName("hmat_percenterrs_foroutputpdf");
    hmat_percenterrs->Write();
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
    
    //---------------          
    canvForPrint->Print(close_outPdfFile.c_str());      
     
    //---------------          
    if(drawRespMatrix){    
      drawRespMatrixFile((TH2*)hmat, (TH2*)hmat_rebin, (TH2*)hmat_errors,
			 //(TH2*)covmat_TH2, (TH2*)PearsonBayes, (TH2*)unfmat_TH2,
			 (double*)boundaries_pt_reco_mat, (int)nbins_pt_reco_mat,
			 (double*)boundaries_pt_gen_mat , (int)nbins_pt_gen_mat,
			 (std::string)outRespMatPdfFile , (bool)useSimpBins);    }//end draw response matrix
    





  }// end draw pdfs

  
  std::cout<<"writing output hists to file... "<<std::endl;
  
  hunf->Write();    
  hunf_sameside->Write();
  
  hfold->Write();    
  hfold_sameside->Write();
  
  h_genratio_oppunf->Write();  
  h_genratio_oppfold->Write();  
  h_genratio_oppmeas->Write();  
  
  h_genratio_ssunf->Write();  
  h_genratio_ssfold->Write();  
  h_genratio_ssmeas->Write();  
  
  h_recratio_oppunf->Write();  
  h_recratio_oppfold->Write();  
  
  h_recratio_ssunf->Write();  
  h_recratio_ssmeas->Write();  
  h_recratio_ssfold->Write();  
  h_recratio_ssgen->Write();  
  
  
  
  
  

  std::cout<<std::endl<<"draw pdfs for Bayesian Unfolding done!"<<std::endl;
  
  if(debugMode) std::cout<<"writing output root file..."<<std::endl;
  fout->Write();
  
  std::cout<<"closing input/output root files..."<<std::endl;
  fout->Close();
  fpp_MC->Close();
  fpp_Data->Close();
  return 0;
  
} // end unfoldDataSpectra



//  steering ---------------------------------------------------------------------------------
int main(int argc, char* argv[]){  int rStatus = -1;
  
  if( argc!=6 ){
    std::cout<<"do ./bayesUnfoldDataSpectra.exe <targDataDir> <targMCDir> <baseOutputName> <doJetID> <useSimpleBins>"<<std::endl;
    std::cout<<"actually... just open the damn code and look"<<std::endl;
    
    return rStatus;  }
  
  rStatus=1; // runtime error
  
  rStatus=bayesUnfoldDataSpectra(  (const std::string)argv[1], (const std::string)argv[2], (const std::string)argv[3], 
				   (int)std::atoi(argv[4]) , (int)std::atoi(argv[5]) 				   
				   ); 
  
  std::cout<<std::endl<<"done!"<<std::endl<<" return status: "<<rStatus<<std::endl<<std::endl;
  return rStatus;
}












