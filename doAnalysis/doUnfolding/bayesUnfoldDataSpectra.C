#include "unfoldSpectra.h"

//other settings
const bool drawPDFs=true; 
const bool debugMode=false;
const bool drawRespMatrix=true;

const double etaBinWidth=4.;//for |y|<2.0

//const double MCscaling=1e-06;//05/09/18
//const double MCscaling=1e+01;//05/29/18
const double MCscaling=1e+00;//06/18/18

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
  if(debugMode)fout->cd();  
  
  
  // ppMC input histos -------------------------
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from MC file"<<std::endl; 
  std::cout<<"input MC dir : "<< (inFile_MC_dir)  <<std::endl; 
  std::cout<<"MC file name : "<< (inFile_MC_name)<<std::endl;   std::cout<<std::endl<<std::endl;  
  TFile *fpp_MC = TFile::Open( (inFile_MC_dir+inFile_MC_name).c_str());
  
  // ---------- open MC "response" matrix
  std::string TH2_title="hpp_matrix";//+RandEtaRange;
  if(doJetID)TH2_title+="_wJetID";
  TH2_title+=RandEtaRange;
  
  //get the response matrix made by readforests
  TH2D* hmat = (TH2D*)fpp_MC->Get( TH2_title.c_str() );
  hmat->Scale(1./MCscaling);//05/09/18
  if(debugMode)hmat->Write(TH2_title.c_str());
  if(debugMode)hmat->Print("base");
  
  // rebinned matrix ---------------
  TH2_title+="_clone";
  TH2D* hmat_rebin = (TH2D*)hmat->Clone( (TH2_title).c_str() );
  if(debugMode)hmat_rebin->Write(TH2_title.c_str());
  if(debugMode)hmat_rebin->Print("base"); 
  
  TH2_title+="_rebins";
  hmat_rebin=(TH2D*) reBinTH2(hmat_rebin, (TH2_title).c_str(), 
			      (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
			      (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );  
  if(debugMode)hmat_rebin->Write(TH2_title.c_str());
  if(debugMode)hmat_rebin->Print("base"); 
  
  if(clearOverUnderflows){
    TH2_title+="_noOverUnderFlows";
    TH2clearOverUnderflows((TH2D*)hmat_rebin);
    if(debugMode)hmat_rebin->Write(TH2_title.c_str());
    if(debugMode)hmat_rebin->Print("base");  }
  
  
  // error and %error for response matrix ---------------
  std::string errTH2_title="hmat_errors_rebins";
  TH2D* hmat_errors=makeRespMatrixErrors( (TH2D*) hmat,
					  (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
					  (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
  if(debugMode)hmat_errors->Write(errTH2_title.c_str());
  if(debugMode)hmat_errors->Print("base");
  
  if(clearOverUnderflows){
    errTH2_title+="_noOverUnderFlows";
    TH2clearOverUnderflows((TH2D*)hmat_errors);
    if(debugMode)hmat_errors->Write(errTH2_title.c_str());
    if(debugMode)hmat_errors->Print("base");  }
  
  // give response matrix the correct errors
  setRespMatrixErrs( (TH2D*)hmat_rebin, (TH2D*) hmat_errors , (bool)zeroBins);
  
  TH2_title+="_wseterrs";
  if(debugMode)hmat_rebin->Write(TH2_title.c_str());
  if(debugMode)hmat_rebin->Print("base");
  
  




  
  // ---------- reco, measured same-side spectra used to create response matrix, and for "sameside" unfolding test
  std::string histTitle2="hpp_reco";//    hpp_reco_wJetID_R4_20_eta_20
  if(doJetID)histTitle2+="_wJetID";
  histTitle2+=RandEtaRange;
  
  TH1D*  hrec_sameside = (TH1D*)fpp_MC->Get( histTitle2.c_str() ); 
  hrec_sameside->Scale(1./MCscaling);//05/09/18
  if(debugMode)hrec_sameside->Write(histTitle2.c_str());
  if(debugMode)hrec_sameside->Print("base");
  
  //assert(false);
  
  histTitle2+="_clone";
  TH1D *hrec_sameside_rebin = (TH1D*)hrec_sameside->Clone( (histTitle2).c_str() );
  if(debugMode)hrec_sameside_rebin->Write( histTitle2.c_str() );
  if(debugMode)hrec_sameside_rebin->Print("base");
  
  
  std::cout<<"rebinning hrec_sameside..."<<std::endl;
  histTitle2+="_rebins";
  hrec_sameside_rebin = (TH1D*)hrec_sameside_rebin->Rebin( nbins_pt_reco, (histTitle2).c_str() , boundaries_pt_reco);
  if(debugMode)hrec_sameside_rebin->Write( histTitle2.c_str() );   
  if(debugMode)hrec_sameside_rebin->Print("base");  
  
  if(clearOverUnderflows){
    histTitle2+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hrec_sameside_rebin);
    if(debugMode)hrec_sameside_rebin->Write( histTitle2.c_str() );
    if(debugMode)hrec_sameside_rebin->Print("base");    }
  
  //cosmetics
  hrec_sameside_rebin->SetMarkerStyle(kOpenSquare);
  hrec_sameside_rebin->SetMarkerColor(kBlue-3);     
  hrec_sameside_rebin->SetLineColor(kBlue-3);     
  hrec_sameside_rebin->SetMarkerSize(1.02);     
  
  
  
  
  // ---------- gen, MC truth spectra
  std::string genHistTitle="hpp_gen";
  if(doJetID)genHistTitle+="_wJetID";
  genHistTitle+=RandEtaRange;
  
  TH1D* hgen = (TH1D*)fpp_MC->Get( genHistTitle.c_str() );
  hgen->Scale(1./MCscaling);//05/09/18
  if(debugMode)hgen->Write(genHistTitle.c_str());
  if(debugMode)hgen->Print("base");    
  
  genHistTitle+="_clone";
  TH1D* hgen_rebin = (TH1D*)hgen->Clone( (genHistTitle).c_str() );
  if(debugMode)hgen_rebin->Write(genHistTitle.c_str());
  if(debugMode)hgen_rebin->Print("base");
  
  genHistTitle+="_rebins";
  hgen_rebin = (TH1D*)hgen_rebin->Rebin(nbins_pt_gen, (genHistTitle).c_str() , boundaries_pt_gen);
  if(debugMode)hgen_rebin->Write(genHistTitle.c_str());
  if(debugMode)hgen_rebin->Print("base"); 
  
  if(clearOverUnderflows){
    genHistTitle+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hgen_rebin);
    if(debugMode)hgen_rebin->Write(genHistTitle.c_str());
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
    if(debugMode)hrec_sameside_resp_rebin->Write( (histTitle3+"_clone").c_str() );    
    
    hgen_resp_rebin = (TH1D*)hgen_rebin->Clone( (histTitle3gen+"_clone").c_str() );      
    if(debugMode)hgen_resp_rebin->Write( (histTitle3gen+"_clone").c_str() );       
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

  if(debugMode)hrec_sameside_resp_rebin->Write((histTitle3+"_final").c_str());
  if(debugMode)hrec_sameside_resp_rebin->Print(" base");  
  
  //cosmetics
  hgen_resp_rebin->SetMarkerStyle(kFullStar);
  hgen_resp_rebin->SetMarkerColor(kMagenta);
  hgen_resp_rebin->SetLineColor(kMagenta);
  hgen_resp_rebin->SetMarkerSize(1.02);   
  
  if(debugMode)hgen_resp_rebin->Write((histTitle3gen+"_final").c_str());
  if(debugMode)hgen_resp_rebin->Print("base");    
  

  
  
  
  // ppData input histos -------------------------
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from Data file"<<std::endl; 
  std::cout<<"input Data dir : "<< (inFile_Data_dir)  <<std::endl; 
  std::cout<<"Data file name : "<< (inFile_Data_name)<<std::endl;   std::cout<<std::endl;//<<std::endl;  
  TFile *fpp_Data = TFile::Open( (inFile_Data_dir+inFile_Data_name).c_str());
  
  //for output
  if(debugMode)fout->cd();    
  
  TH1D* hJetQA_jtptEntries=(TH1D*)fpp_Data->Get("hJetQA_1wJetID_jtptEntries");
  if(hJetQA_jtptEntries){
    hJetQA_jtptEntries = (TH1D*)hJetQA_jtptEntries->Rebin(nbins_pt_reco, hJetQA_jtptEntries->GetTitle(), boundaries_pt_reco);  
  }
  
  // ---------- reco, measured spectra to unfold
  std::string histTitle="hJetQA";//hJetQA_1wJetID_jtpt
  if(doJetID)histTitle+="_1wJetID";
  else histTitle+="_0wJetID";
  histTitle+="_jtpt";
  //histTitle+=RandEtaRange;
  
  TH1D*  hrec = (TH1D*)fpp_Data->Get( histTitle.c_str() ); 
  if(debugMode)hrec->Write();
  if(debugMode)hrec->Print("base");
  
  histTitle+="_divBylumietabin";
  double effIntgrtdLumi=computeEffLumi(fpp_Data);
  hrec->Scale(1./effIntgrtdLumi); // lumi
  if(debugMode)hrec->Write( (histTitle).c_str() );
  if(debugMode)hrec->Print("base");
  
  //assert(false);
  
  histTitle+="_clone";
  TH1D *hrec_rebin = (TH1D*)hrec->Clone( (histTitle).c_str() );
  if(debugMode)hrec_rebin->Write(histTitle.c_str());
  if(debugMode)hrec_rebin->Print("base");
  
  std::cout<<"rebinning hrec..."<<std::endl;
  histTitle+="_rebins";
  hrec_rebin = (TH1D*)hrec_rebin->Rebin( nbins_pt_reco, (histTitle).c_str() , boundaries_pt_reco);
  if(debugMode)hrec_rebin->Write(histTitle.c_str());   
  if(debugMode)hrec_rebin->Print("base");  
  
  
  if(clearOverUnderflows){
    histTitle+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hrec_rebin);
    if(debugMode)hrec_rebin->Write(histTitle.c_str());
    if(debugMode)hrec_rebin->Print("base");  
  }
  
  //cosmetics
  hrec_rebin->SetMarkerStyle(kOpenCircle);
  hrec_rebin->SetMarkerColor(kBlue);     
  hrec_rebin->SetLineColor(kBlue);     
  hrec_rebin->SetMarkerSize(1.02);     

  
  
  
  // thy spectra
  
  TH1D* CT10nlo  =(TH1D*)makeThyHist_00eta20((fNLOFile_R04_CT10nlo  ).c_str());
  CT10nlo->SetMarkerSize(0);
  CT10nlo->SetLineColor(kBlack);  
  TH1D* CT14nlo  =(TH1D*)makeThyHist_00eta20((fNLOFile_R04_CT14nlo  ).c_str());
  CT14nlo->SetMarkerSize(0);
  CT14nlo->SetLineColor(kGreen);  
  TH1D* HERAPDF  =(TH1D*)makeThyHist_00eta20((fNLOFile_R04_HERAPDF  ).c_str());
  HERAPDF->SetMarkerSize(0);
  HERAPDF->SetLineColor(kViolet-5);  
  TH1D* MMHTnlo  =(TH1D*)makeThyHist_00eta20((fNLOFile_R04_MMHTnlo  ).c_str());
  MMHTnlo->SetMarkerSize(0);
  MMHTnlo->SetLineColor(kOrange+7);  
  TH1D* NNPDFnnlo=(TH1D*)makeThyHist_00eta20((fNLOFile_R04_NNPDFnnlo).c_str());
  NNPDFnnlo->SetMarkerSize(0);
  NNPDFnnlo->SetLineColor(kCyan-6);  

  
  
  
  
  
  // Bayesian unfolding -------------------------   
  // RooUnfoldResponse, Bayes --------------------
  std::cout<<"calling RooUnfoldResponse "<<std::endl;
  
  RooUnfoldResponse roo_resp( hrec_sameside_resp_rebin, hgen_resp_rebin, hmat_rebin, ("Response_matrix"+RandEtaRange).c_str()) ;
  //RooUnfoldResponse roo_resp( 0, 0, hmat_rebin, ("Response_matrix"+RandEtaRange).c_str()) ;
  roo_resp.UseOverflow(doOverUnderflows);    
  if(debugMode)roo_resp.Write();

  TH1D* hfak=  (TH1D*) roo_resp.Hfakes() ;
  hfak->Print("base");
  //cosmetics
  hfak->SetMarkerColor(kGreen);
  hfak->SetLineColor(kGreen);
  
  
  std::cout<<"calling RooUnfoldBayes..."<<std::endl;
  RooUnfoldBayes unf_bayes( &roo_resp, hrec_rebin, kIter );
  unf_bayes.SetVerbose(2);
  if(doToyErrs){
    std::cout<<"using toy errors, suppressing text output"<<std::endl;
    unf_bayes.SetNToys(10000);
    unf_bayes.SetVerbose(1);
  }
  std::cout<<"RooUnfoldBayes Overflow Status: " << unf_bayes.Overflow()<<std::endl;
  
  TH1D *hunf = (TH1D*)unf_bayes.Hreco(errorTreatment);     std::cout<<std::endl; 
  hunf->SetName("ppData_BayesUnf_Spectra");
  hunf->SetTitle("Unf. Data");
  if(debugMode)hunf->Print("base");
  
  //cosmetics
  hunf->SetMarkerStyle(kOpenCircle);
  hunf->SetMarkerColor(kRed);
  hunf->SetLineColor(kRed);
  hunf->SetMarkerSize(1.02);     
  
  std::cout<<"folding unfolded data histogram!!"<<std::endl;
  TH1D* hfold=(TH1D*)roo_resp.ApplyToTruth(hunf);          
  hfold->SetName("ppData_BayesFold_Spectra");
  hfold->SetTitle("Fold. Data");
  if(debugMode)hfold->Print("base");
  
  
  //cosmetics
  hfold->SetMarkerStyle(kOpenCircle);
  hfold->SetMarkerColor(kGreen-5);
  hfold->SetLineColor(  kGreen-5);
  hfold->SetMarkerSize(1.02);     
  
  
  std::cout<<"folding MC Truth histogram!!"<<std::endl;
  TH1D* hfold_truth=(TH1D*)roo_resp.ApplyToTruth(hgen_rebin);          
  hfold_truth->SetName("ppMC_BayesFoldMCTruth_Spectra");
  hfold_truth->SetTitle("Fold. MC Truth");
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
  if(debugMode)covmat_TH2->Write();
  
  std::cout<<"calculating pearson coefficients"<<std::endl;
  TH2D* PearsonBayes = (TH2D*)CalculatePearsonCoefficients( &covmat, false,"Bayes_pearson");
  if(debugMode)PearsonBayes->Print("base");
  if(debugMode)PearsonBayes->Write();

  std::cout<<"getting unfolding matrix"<<std::endl;
  TH2D* unfmat_TH2=new TH2D(unf_bayes.UnfoldingMatrix());
  unfmat_TH2->SetName("Bayes_unfoldingMatrix_TH2");  
  if(debugMode)unfmat_TH2->Print("base");
  if(debugMode)unfmat_TH2->Write();
  
  



  // --------- RATIOS WITH MC TRUTH ----------------
  TH1D *hgen_rebin_ratiobin=(TH1D*)hgen_rebin->Clone("ppData_Gen_Ratio_denom");
  if(debugMode)hgen_rebin_ratiobin->Print("base");
  
  TH1D *h_genratio_oppunf = (TH1D*)hunf->Clone( "ppData_Gen_Ratio_OppUnf" );
  h_genratio_oppunf->SetTitle( "Unf. Data/GEN Py8" );
  h_genratio_oppunf->Divide(hgen_rebin_ratiobin);
  if(debugMode)h_genratio_oppunf->Print("base");
  
  TH1D *h_genratio_oppfold = (TH1D*)hfold->Clone( "ppData_Gen_Ratio_OppFold" );
  if(nbins_pt_reco>nbins_pt_gen)
    h_genratio_oppfold=(TH1D*)h_genratio_oppfold->Rebin(nbins_pt_gen,"ppData_Gen_Ratio_OppFold_rebin",boundaries_pt_gen);
  h_genratio_oppfold->SetTitle( "Fold. Data/GEN Py8" );
  h_genratio_oppfold->Divide(hgen_rebin_ratiobin);
  if(debugMode)h_genratio_oppfold->Print("base");
  
  TH1D *h_genratio_oppmeas = (TH1D*)hrec_rebin->Clone( "ppData_Gen_Ratio_Meas" );
  if(nbins_pt_reco>nbins_pt_gen)
    h_genratio_oppmeas=(TH1D*)h_genratio_oppmeas->Rebin(nbins_pt_gen,"ppData_Gen_Ratio_Meas_rebin",boundaries_pt_gen);
  h_genratio_oppmeas->SetTitle( "Meas. Data/GEN Py8" );
  h_genratio_oppmeas->Divide(hgen_rebin_ratiobin);
  if(debugMode)h_genratio_oppmeas->Print("base");

  
  TH1D *h_genratio_ssmeas = (TH1D*)hrec_sameside_rebin->Clone( "ppMC_Gen_Ratio4_SSMeas" );
  h_genratio_ssmeas->SetTitle( "RECO Py8/GEN Py8" );
  h_genratio_ssmeas->Divide(hgen_rebin_ratiobin);
  if(debugMode)h_genratio_ssmeas->Print("base");
  
  




  // ---------------- RATIOS WITH DATA MEAS -----------------
  TH1D *hrec_rebin_ratiobin=(TH1D*)hrec_rebin->Clone("ppData_Meas_Ratio_denom");
  hrec_rebin_ratiobin=(TH1D*)hrec_rebin_ratiobin->Rebin(nbins_pt_gen,"ppData_Meas_Ratio_denom_rebin",boundaries_pt_gen);
  if(debugMode)hgen_rebin_ratiobin->Print("base");
  
  
  TH1D *h_recratio_oppunf = (TH1D*)hunf->Clone( "ppData_Meas_Ratio_OppUnf" );
  h_recratio_oppunf->SetTitle( "Unf. Data/RECO Data" );
  h_recratio_oppunf->Divide(hrec_rebin_ratiobin);
  if(debugMode)h_recratio_oppunf->Print("base");
  
  TH1D *h_recratio_oppfold = (TH1D*)hfold->Clone( "ppData_Meas_Ratio_OppFold" );
  h_recratio_oppfold->SetTitle( "Fold. Data/RECO Data" );
  h_recratio_oppfold->Divide(hrec_rebin);
  if(debugMode)h_recratio_oppfold->Print("base");
  
  
  
  TH1D *h_recratio_ssmeas = (TH1D*)hrec_sameside_rebin->Clone( "ppMC_Meas_Ratio_SSMeas" );
  h_recratio_ssmeas->SetTitle( "RECO Py8/RECO Data" );
  h_recratio_ssmeas->Divide(hrec_rebin);
  if(debugMode)h_recratio_ssmeas->Print("base");
  
  
  TH1D *h_recratio_ssgen = (TH1D*)hgen_rebin->Clone( "ppMC_Meas_Ratio_SSTruth" );    
  h_recratio_ssgen->SetTitle( "GEN Py8/RECO Data" );
  h_recratio_ssgen->Divide(hrec_rebin_ratiobin);
  if(debugMode)h_recratio_ssgen->Print("base");
  
  // ---------------- FOLDED RATIOS (MC AND DATA) ----------------- //
  
  TH1D *h_foldratio_datafold=(TH1D*)hfold->Clone("ppData_Fold_Ratio");
  h_foldratio_datafold->SetTitle( "Data Fold(Unf.)/RECO Data" );
  h_foldratio_datafold->Divide(hrec_rebin_ratiobin);
  if(debugMode)h_foldratio_datafold->Print("base");
  
  
  TH1D *h_foldratio_mcfold=(TH1D*)hfold_truth->Clone("ppMC_Fold_Ratio");
  h_foldratio_mcfold->SetTitle( "Fold GEN Py8/RECO Py8" );
  h_foldratio_mcfold->Divide(hrec_sameside_rebin);
  if(debugMode)h_foldratio_mcfold->Print("base");
  
  


  // ----- PUT BIN WIDTH(s) DIVISIONS + NORMALIZING HERE ----- //
  
  // -- MC RECO -- //  
  histTitle2+="_normbinwidth";
  hrec_sameside_rebin->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hrec_sameside_rebin); 
  if(debugMode)hrec_sameside_rebin->Write( histTitle2.c_str() );
  if(debugMode)hrec_sameside_rebin->Print("base");  
  
  // -- MC GEN -- //  
  genHistTitle+="_normbinwidth";
  hgen_rebin->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hgen_rebin);
  if(debugMode)hgen_rebin->Write(genHistTitle.c_str());
  if(debugMode)hgen_rebin->Print("base");  

  // -- MC GEN RESP -- //  
  hgen_resp_rebin->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hgen_resp_rebin);
  if(debugMode)hgen_resp_rebin->Write(genHistTitle.c_str());
  if(debugMode)hgen_resp_rebin->Print("base");  

  // -- MC RESP FAKES -- //  
  hfak->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hfak);
  if(debugMode)hfak->Write();
  if(debugMode)hfak->Print("base");  
  
  // -- MC REC RESP -- //  
  hrec_sameside_resp_rebin->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  divideBinWidth(hrec_sameside_resp_rebin);
  if(debugMode)hrec_sameside_resp_rebin->Write();
  if(debugMode)hrec_sameside_resp_rebin->Print("base");  
  
  // -- MC FOLD(TRUTH) -- //
  hfold_truth->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hfold_truth); 
  if(debugMode)hfold_truth->Write();
  if(debugMode)hfold_truth->Print("base");  
  
  // -- DATA RECO -- //  
  histTitle+="_normbinwidth";
  hrec_rebin->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hrec_rebin); 
  if(debugMode)hrec_rebin->Write(histTitle.c_str());
  if(debugMode)hrec_rebin->Print("base");  
  
  // -- DATA UNF -- //
  hunf->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hunf); 
  if(debugMode)hunf->Write();
  if(debugMode)hunf->Print("base");  
  
  // -- DATA FOLD(UNF) -- //
  hfold->Scale(1./etaBinWidth); // |y| bin width
  divideBinWidth(hfold); 
  if(debugMode)hfold->Write();
  if(debugMode)hfold->Print("base");  


  // ----- END BIN WIDTH DIVISIONS + NORMALIZING ----- //
  
  // ---------------- THY RATIOS (THY / DATA) ----------------- //
  // must be after the binwidth divisions + normalization (the thy hists are made this way by default)
  
  TH1D* h_thyratio_CT10nlo  =(TH1D*)CT10nlo  ->Clone("");
  h_thyratio_CT10nlo = (TH1D*)h_thyratio_CT10nlo->Rebin(nbins_pt_gen,"pp_CT10Thy_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_CT10nlo  ->SetTitle("CT10 NLO/Data Unf.");  
  h_thyratio_CT10nlo  ->Divide(hunf);
  if(debugMode)h_thyratio_CT10nlo  ->Print("base");
  
  TH1D* h_thyratio_CT14nlo  =(TH1D*)CT14nlo  ->Clone("");  
  h_thyratio_CT14nlo=(TH1D*)h_thyratio_CT14nlo->Rebin(nbins_pt_gen,"pp_CT14Thy_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_CT14nlo  ->SetTitle("CT14 NLO/Data Unf.");
  h_thyratio_CT14nlo  ->Divide(hunf);
  if(debugMode)h_thyratio_CT14nlo  ->Print("base");
  
  TH1D* h_thyratio_HERAPDF  =(TH1D*)HERAPDF  ->Clone("");  
  h_thyratio_HERAPDF=(TH1D*)h_thyratio_HERAPDF->Rebin(nbins_pt_gen,"pp_HERAPDF_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_HERAPDF  ->SetTitle("HERAPDF15 NLO/Data Unf.");
  h_thyratio_HERAPDF  ->Divide(hunf);
  if(debugMode)h_thyratio_HERAPDF  ->Print("base");
  
  TH1D* h_thyratio_MMHTnlo  =(TH1D*)MMHTnlo  ->Clone("");
  h_thyratio_MMHTnlo=(TH1D*)h_thyratio_MMHTnlo->Rebin(nbins_pt_gen,"pp_MMHT_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_MMHTnlo  ->SetTitle("MMHT14 NLO/Data Unf.");
  h_thyratio_MMHTnlo  ->Divide(hunf);
  if(debugMode)h_thyratio_MMHTnlo  ->Print("base");
  
  TH1D* h_thyratio_NNPDFnnlo=(TH1D*)NNPDFnnlo->Clone("");
  h_thyratio_NNPDFnnlo=(TH1D*)h_thyratio_NNPDFnnlo->Rebin(nbins_pt_gen,"pp_NNPDFnlo_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_NNPDFnnlo->SetTitle("NNPDF NLO/Data Unf.");
  h_thyratio_NNPDFnnlo->Divide(hunf);
  if(debugMode)h_thyratio_NNPDFnnlo->Print("base");
  

  
  
  
    

  
  
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
					,(boundaries_pt_gen_mat[nbins_pt_gen_mat])
					,1.);
    if(debugMode)std::cout<<"theLineAtOne_gen starts at "<<boundaries_pt_gen_mat[0]  <<std::endl;
    if(debugMode)std::cout<<"theLineAtOne_gen ends at   "<<(boundaries_pt_gen_mat[nbins_pt_gen_mat])  <<std::endl;
    theLineAtOne_gen->SetLineWidth(1);
    theLineAtOne_gen->SetLineStyle(2);
    theLineAtOne_gen->SetLineColor(36);    
    
    TLine* theLineAtp9_gen= new TLine( boundaries_pt_gen_mat[0]  
				       ,0.9   
				       ,(boundaries_pt_gen_mat[nbins_pt_gen_mat])
				       ,0.9);
    if(debugMode)std::cout<<"theLineAtp9_gen starts at "<<boundaries_pt_gen_mat[0]  <<std::endl;
    if(debugMode)std::cout<<"theLineAtp9_gen ends at   "<<(boundaries_pt_gen_mat[nbins_pt_gen_mat])  <<std::endl;
    theLineAtp9_gen->SetLineWidth(1);
    theLineAtp9_gen->SetLineStyle(2);
    theLineAtp9_gen->SetLineColor(36);
    
    TLine* theLineAt1p1_gen= new TLine( boundaries_pt_gen_mat[0]  
					,1.1   
					,(boundaries_pt_gen_mat[nbins_pt_gen_mat])
					,1.1);
    if(debugMode)std::cout<<"theLineAt1p1_gen starts at "<<boundaries_pt_gen_mat[0]  <<std::endl;
    if(debugMode)std::cout<<"theLineAt1p1_gen ends at   "<<(boundaries_pt_gen_mat[nbins_pt_gen_mat])  <<std::endl;
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
    
    
    TLine* theLineAtp9_reco= new TLine( boundaries_pt_reco_mat[0]  
					,0.9   
					,(boundaries_pt_reco_mat[nbins_pt_reco_mat])
					,0.9);
    if(debugMode)std::cout<<"theLineAtp9_reco starts at "<<boundaries_pt_reco_mat[0]  <<std::endl;
    if(debugMode)std::cout<<"theLineAtp9_reco ends at   "<<(boundaries_pt_reco_mat[nbins_pt_reco_mat])  <<std::endl;
    theLineAtp9_reco->SetLineWidth(1);
    theLineAtp9_reco->SetLineStyle(2);
    theLineAtp9_reco->SetLineColor(36);
    
    TLine* theLineAt1p1_reco= new TLine( boundaries_pt_reco_mat[0]  
					 ,1.1   
					 ,(boundaries_pt_reco_mat[nbins_pt_reco_mat])
					 ,1.1);
    if(debugMode)std::cout<<"theLineAt1p1_reco starts at "<<boundaries_pt_reco_mat[0]  <<std::endl;
    if(debugMode)std::cout<<"theLineAt1p1_reco ends at   "<<(boundaries_pt_reco_mat[nbins_pt_reco_mat])  <<std::endl;
    theLineAt1p1_reco->SetLineWidth(1);
    theLineAt1p1_reco->SetLineStyle(2);
    theLineAt1p1_reco->SetLineColor(36);
    
    

    

    // gen/meas/unf spectra hists ---------------    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(1);    
    
    setupSpectraHist(hgen_rebin          , useSimpBins);    
    setupSpectraHist(hrec_rebin	       , useSimpBins);
    setupSpectraHist(hrec_sameside_rebin , useSimpBins);
    setupSpectraHist(hunf		       , useSimpBins);    
    
    //hgen_rebin->SetTitle("MC, Data, and Bayes Unf. Spectra");
    hrec_rebin->SetTitle("MC, Data, and Bayes Unf. Spectra");
    
    
    hrec_rebin->Draw("P E");           
    hgen_rebin->Draw("P E SAME");                 
    hunf->Draw("P E SAME");               
    //hrec_sameside_rebin->Draw("P E SAME");           
    
    TLegend* legend_in1 = new TLegend( 0.7,0.7,0.9,0.9 );
    //legend_in1->AddEntry(hrec_sameside_rebin, "MC Meas."  ,  "lp");
    legend_in1->AddEntry(hrec_rebin,          "Data Meas." , "lp");	
    legend_in1->AddEntry(hunf,                "Data Unf." ,  "lp");
    legend_in1->AddEntry(hgen_rebin,          "MC Truth"   , "lp");
    
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
    //hrec_sameside_resp_rebin->SetTitle("MC Only, Response Spectra and Fakes");
    hfak->SetTitle("MC Only, Response Spectra and Fakes");
    //hrec_sameside_resp_rebin->GetYaxis()->SetTitle("N_{Jets}/L_{int}");    
    
    hfak->Draw("P E");           
    hrec_sameside_resp_rebin->Draw("P E SAME");           
    hgen_resp_rebin->Draw("P E SAME");                 

    TLegend* legend_resp = new TLegend( 0.6,0.7,0.9,0.9 );
    legend_resp->AddEntry(hfak,          "MC Response Matrix Fakes" , "lp");    
    legend_resp->AddEntry(hgen_resp_rebin,          "MC Truth Response" , "lp");
    legend_resp->AddEntry(hrec_sameside_resp_rebin, "MC Meas. Response" , "lp");
    
    legend_resp->Draw();
    
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
    //hfold->Draw("P E");           
    //hfold_truth->Draw("P E SAME");           
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
    
    
    // thy spectra CT10/14 NNPDF NLO---------------------------
    TLegend* legendThy1 =new TLegend( 0.7,0.7,0.9,0.9 );    

    //CT10nlo->SetLineColor(kRed);
    setupSpectraHist(CT10nlo  ,useSimpBins);
    legendThy1->AddEntry(CT10nlo  ,"CT10 NLO","l");
    if(debugMode){fout->cd(); CT10nlo->Write("CT10_NLO_R04_jtpt");}

    //CT14nlo->SetLineColor(kGreen);
    setupSpectraHist(CT14nlo  ,useSimpBins);
    legendThy1->AddEntry(CT14nlo  ,"CT14 NLO","l");
    if(debugMode){fout->cd(); CT14nlo->Write("CT14_NLO_R04_jtpt");}
    
    if(debugMode)fout->cd();//makeThyHist will make me cd to last thy file opened...
    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(1);                  
    
    legendThy1->AddEntry(hunf,"Data Unf.","lp");
    legendThy1->AddEntry(hgen_rebin,"MC Truth", "lp");
        
    CT10nlo->SetAxisRange(1e-07,1e+03,"Y");//for y axis in nanbarns
    CT10nlo->SetTitle("NLO Thy w/ Data Unf., MC Truth");

    CT10nlo  ->Draw("][HIST E");
    CT14nlo  ->Draw("][HIST E SAME");
    
    hgen_rebin->Draw("P E SAME");
    hunf->Draw("P E SAME");   //just for axis range
    

    legendThy1->Draw();

    canvForPrint->Print(outPdfFile.c_str());

    
    // thy spectra HERA/MMHT/NNPDF ---------------------------

    TLegend* legendThy =new TLegend( 0.7,0.7,0.9,0.9 );    

    //HERAPDF->SetLineColor(kBlue);
    setupSpectraHist(HERAPDF  ,useSimpBins);
    legendThy->AddEntry(HERAPDF  ,"HERAPDF 2015 NLO","l");
    if(debugMode){fout->cd(); HERAPDF->Write("HERAPDF105_NLO_R04_jtpt");}
    
    //MMHTnlo->SetLineColor(kMagenta);
    setupSpectraHist(MMHTnlo  ,useSimpBins);
    legendThy->AddEntry(MMHTnlo  ,"MMHT 2014 NLO   ","l");
    if(debugMode){fout->cd(); MMHTnlo->Write("MMHT2014_NLO_R04_jtpt"); }
    
    //NNPDFnnlo->SetLineColor(kBlack);
    setupSpectraHist(NNPDFnnlo,useSimpBins);
    legendThy->AddEntry(NNPDFnnlo,"NNPDF NNLO","l");
    if(debugMode){fout->cd(); NNPDFnnlo->Write("NNPDF_NLO_R04_jtpt");    }
    
    if(debugMode)fout->cd();//makeThyHist will make me cd to last thy file opened...
    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(1);                  
    
    legendThy->AddEntry(hunf,"Data Unf.","lp");
    legendThy->AddEntry(hgen_rebin,"MC Truth", "lp");
    
    HERAPDF->SetAxisRange(1e-07,1e+03,"Y");//for y axis in nanbarns
    HERAPDF->SetTitle("NLO Thy w/ Data Unf., MC Truth");
    
    HERAPDF  ->Draw("][HIST E");
    MMHTnlo  ->Draw("][HIST E SAME");
    NNPDFnnlo->Draw("][HIST E SAME");
    
    hgen_rebin->Draw("P E SAME");
    hunf->Draw("P E SAME");   //just for axis range
    

    legendThy->Draw();

    canvForPrint->Print(outPdfFile.c_str());


    // ratios w MC Truth------------    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(0);
    
    setupRatioHist(h_genratio_oppunf , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_genratio_oppmeas, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_genratio_ssmeas , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    //setupRatioHist(h_genratio_oppmeas, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    //setupRatioHist(h_genratio_oppfold, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    //setupRatioHist(h_genratio_ssmeas, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);    
    
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




    // thy ratios w hunf------------    
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(0);
    
    setupRatioHist(h_thyratio_CT10nlo , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_thyratio_CT14nlo , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_thyratio_HERAPDF , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_thyratio_MMHTnlo , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_thyratio_NNPDFnnlo, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    
    h_thyratio_CT10nlo->SetTitle( "Thy Ratios w/ Bayes Unf. Data" );
    h_thyratio_CT10nlo->GetYaxis()->SetTitle("Thy / Bayes Unf. Data");
    
    h_thyratio_CT10nlo ->Draw( "][HIST ");      
    h_thyratio_CT14nlo ->Draw( "][HIST SAME"); 
    h_thyratio_HERAPDF ->Draw( "][HIST SAME"); 
    h_thyratio_MMHTnlo ->Draw( "][HIST SAME"); 
    h_thyratio_NNPDFnnlo->Draw("][HIST SAME"); 

    
    TLegend* legendthyrat = new TLegend( 0.1,0.7,0.3,0.9 );
    legendthyrat->AddEntry(h_thyratio_CT10nlo ,  "CT10 PDF NLO" ,    "l");
    legendthyrat->AddEntry(h_thyratio_CT14nlo ,  "CT14 PDF NLO" ,    "l"); 
    legendthyrat->AddEntry(h_thyratio_HERAPDF ,  "HERAPDF 2015 NLO", "l");
    legendthyrat->AddEntry(h_thyratio_MMHTnlo ,  "MMHT 2014 NLO",    "l");
    legendthyrat->AddEntry(h_thyratio_NNPDFnnlo, "NNPDF NNLO",       "l");
    
    legendthyrat->Draw();
    
    theLineAtOne_gen->Draw();
    theLineAtp9_gen->Draw();
    theLineAt1p1_gen->Draw();
 


    canvForPrint->Print(outPdfFile.c_str());


    
    

    // Folding Ratios -------------------
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(0);
    
    setupRatioHist(h_foldratio_datafold, useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);
    setupRatioHist(h_foldratio_mcfold, useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);
    
    h_foldratio_datafold->SetTitle( "Folded Spectra Ratios w/ RECO Spectra" );
    h_foldratio_datafold->GetYaxis()->SetTitle("Ratio w/ Meas.");

    h_foldratio_datafold->Draw("P E");
    h_foldratio_mcfold->Draw("P E SAME");
    
    TLegend* legendfold = new TLegend( 0.1,0.8,0.4,0.9 );
    
    legendfold->AddEntry(h_foldratio_datafold ,  "Fold(Data Unf.)/Data Meas." , "lp");
    legendfold->AddEntry(h_foldratio_mcfold  ,   "Fold(MC Truth)/MC Meas.",  "lp");
    
    legendfold->Draw();
    
    
    theLineAtOne_gen->Draw();
    theLineAt1p1_gen->Draw();
    theLineAtp9_gen->Draw();

    canvForPrint->Print(outPdfFile.c_str());
    
    
    
    // ratio w/ data meas -------------------
    canvForPrint->cd();
    if(!useSimpBins)canvForPrint->SetLogx(1);
    canvForPrint->SetLogy(0);      
    
    //h_recratio_oppunf->SetTitle("Bayesian, Ratios w/ Data Meas. Spectra");
    //h_recratio_oppunf->GetYaxis()->SetTitle("Ratio w/ Data Meas.");

    h_recratio_ssmeas->SetTitle("Bayesian, Ratios w/ Data Meas. Spectra");
    h_recratio_ssmeas->GetYaxis()->SetTitle("Ratio w/ Data Meas.");
    
    setupRatioHist(h_recratio_ssmeas , useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);    
    setupRatioHist(h_recratio_oppunf , useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);    
    setupRatioHist(h_recratio_ssgen  , useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);    
    
    h_recratio_ssmeas->Draw("P E");
    h_recratio_oppunf->Draw("P E SAME");
    h_recratio_ssgen->Draw("P E SAME");
    ////h_recratio_ssunf->Draw("P E SAME");
    //h_recratio_oppfold->Draw("P E SAME");  
    //h_recratio_ssfold->Draw("P E SAME");
    
    TLegend* legend4 = new TLegend( 0.1,0.8,0.2,0.9 );
    legend4->AddEntry(h_recratio_ssmeas, "MC Meas.", "lp");
    legend4->AddEntry(h_recratio_oppunf,  "Data Unf.", "lp");
    legend4->AddEntry(h_recratio_ssgen,   "MC Truth", "lp");
    ////legend4->AddEntry(h_recratio_ssunf, NULL, "p");
    //legend4->AddEntry(h_recratio_oppfold, "Data Fold(Unf.)", "lp");
    //legend4->AddEntry(h_recratio_ssfold, "MC Fold(Unf.)", "p");
    legend4->Draw();


    theLineAtOne_reco->Draw();
    theLineAtp9_reco->Draw();
    theLineAt1p1_reco->Draw();
    

    
    canvForPrint->Print(outPdfFile.c_str());
    


    
    // NJets v. Jet p_T hist (all weight=1) ---------------------
    
    if( (bool)hJetQA_jtptEntries )
      {
	canvForPrint->cd();
	if(!useSimpBins)canvForPrint->SetLogx(1);
	canvForPrint->SetLogy(1);
	
	setupSpectraHist(hJetQA_jtptEntries, useSimpBins);
	hJetQA_jtptEntries->GetYaxis()->SetTitle("NJets/Bin [weights=1]");
	
	hJetQA_jtptEntries->Draw("HIST E");
	canvForPrint->Print(outPdfFile.c_str());
      }

    // covariance matrix ---------------      
    canvForPrint->cd();
    canvForPrint->SetLogx(0);
    canvForPrint->SetLogy(0);
    canvForPrint->SetLogz(1);
      
      //always use simp bins for covmat to avoid log scaling the x/y axes
    matStylePrint( (TH2D*)covmat_TH2, "Bayesian Unfolding Covariance Matrix", canvForPrint, outPdfFile, true);
    
    // pearson matrix ---------------      
    canvForPrint->cd();
    canvForPrint->SetLogx(0);
    canvForPrint->SetLogy(0);
    canvForPrint->SetLogz(0);
    
    //always use simp bins for covmat to avoid log scaling the x/y axes
    matStylePrint( (TH2D*)PearsonBayes, "Bayesian Unfolding Pearson Matrix", canvForPrint, outPdfFile, true);
    
    
    // unfolding matrix ---------------      
    canvForPrint->cd();
    canvForPrint->SetLogx(0);
    canvForPrint->SetLogy(0);
    canvForPrint->SetLogz(1);
      
    //always use simp bins for covmat to avoid log scaling the x/y axes
    matStylePrint( (TH2D*)unfmat_TH2, "Bayesian Unfolding Matrix", canvForPrint, outPdfFile, true);



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
    TH2D* hmat_percenterrs= makeRespMatrixPercentErrs( (TH2D*) hmat_errors, (TH2D*) hmat_rebin,
						       (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
						       (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );		     
    hmat_percenterrs->SetName("hmat_percenterrs_foroutputpdf");
    if(debugMode)hmat_percenterrs->Write();
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
      drawRespMatrixFile((TH2D*)hmat, (TH2D*)hmat_rebin, (TH2D*)hmat_errors,
			 //(TH2D*)covmat_TH2, (TH2D*)PearsonBayes, (TH2D*)unfmat_TH2,
			 (double*)boundaries_pt_reco_mat, (int)nbins_pt_reco_mat,
			 (double*)boundaries_pt_gen_mat , (int)nbins_pt_gen_mat,
			 (std::string)outRespMatPdfFile , (bool)useSimpBins, (TFile*)fpp_MC);    }//end draw response matrix
    





  }// end draw pdfs

  fout->cd();
  
  std::cout<<"writing output hists to file... "<<std::endl;
  
  // input data ------------------
  hrec_rebin->Write();
  
  // input MC ---------------------
  hgen_rebin->Write();
  hrec_sameside_rebin->Write();
  if((bool)hJetQA_jtptEntries)      hJetQA_jtptEntries->Write();
  
  //response matrix
  hmat->Write();
  hmat_rebin->Write();
  hmat_errors->Write();    
  
  // input thy ----------------
  CT10nlo  ->Write("CT10_NLO_R04_jtpt");
  CT14nlo  ->Write("CT14_NLO_R04_jtpt");
  HERAPDF  ->Write("HERAPDF105_NLO_R04_jtpt");
  MMHTnlo  ->Write("MMHT2014_NLO_R04_jtpt");
  NNPDFnnlo->Write("NNPDF_NLO_R04_jtpt");
  
  // output data -------------
  hunf->Write();    
  hfold->Write();        
  hfold_truth->Write();
  
  // output gen ratios (denom=mc truth)
  h_genratio_oppunf->Write();  //data unf/mc truth	   
  h_genratio_oppfold->Write(); //data fold(unf)/mc truth 
  h_genratio_oppmeas->Write(); //data meas/mc truth	   
  h_genratio_ssmeas->Write();  //mc meas/mc truth        
  
  // output rec ratios (denom=data meas)
  h_recratio_oppunf->Write();  //data unf       / data meas	   
  h_recratio_oppfold->Write(); //data fold(unf) / data meas     
  h_recratio_ssmeas->Write();  //mc meas        / data meas	   
  h_recratio_ssgen->Write();   //mc truth       / data meas        
  
  //output fold ratio test
  h_foldratio_datafold->Write();
  h_foldratio_mcfold->Write();
  
  // thy ratios w/ unfolded data
  h_thyratio_CT10nlo ->Write();
  h_thyratio_CT14nlo ->Write();
  h_thyratio_HERAPDF ->Write();
  h_thyratio_MMHTnlo ->Write();    
  h_thyratio_NNPDFnnlo->Write();
  
  //old but may reappear at some point?!
  //h_recratio_ssunf->Write();  
  //h_recratio_ssfold->Write();  
  //h_genratio_ssunf->Write();  
  //h_genratio_ssfold->Write();  
  
    
  

  
  
  
  

  std::cout<<std::endl<<"draw pdfs for Bayesian Unfolding done!"<<std::endl;
  
  //if(debugMode) std::cout<<"writing output root file..."<<std::endl;
  //fout->Write();
  
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

























  //std::cout<<"calling sameside RooUnfoldBayes..."<<std::endl;
  //RooUnfoldBayes unf_sameside_bayes( &roo_resp, hrec_sameside_rebin, kIter );
  //RooUnfoldBayes unf_sameside_bayes( &roo_resp, hrec_sameside_resp_rebin, kIter );
  
  //TH1D *hunf_sameside = (TH1D*)unf_sameside_bayes.Hreco(errorTreatment);    std::cout<<std::endl; 
  //hunf_sameside->SetName("ppMC_BayesUnf_sameSideSpectra");
  //hunf_sameside->SetTitle("Unf. Py8");
  //if(debugMode)hunf_sameside->Print("base");
  //
  ////cosmetics
  //hunf_sameside->SetMarkerStyle(kOpenSquare);
  //hunf_sameside->SetMarkerColor(kRed-3);
  //hunf_sameside->SetMarkerSize(1.02);     
  //
  //
  //std::cout<<"reFolding unfolded sameside histogram!!"<<std::endl;
  //TH1D* hfold_sameside=(TH1D*)roo_resp.ApplyToTruth(hunf_sameside);    
  //hfold_sameside->SetName("ppMC_BayesFold_sameSide_Spectra");
  //hfold_sameside->SetTitle("Fold. Py8");
  //if(debugMode)hfold_sameside->Print("base");
  //
  ////cosmetics
  //hfold_sameside->SetMarkerStyle(kOpenSquare);
  //hfold_sameside->SetMarkerColor(kGreen-3);
  //hfold_sameside->SetMarkerSize(1.02);     
  
  //hfold_sameside->Scale(hrec_sameside_rebin->Integral()/hfold_sameside->Integral());
  //if(debugMode)hfold_sameside->Print("base");  
  
  
  // SAMESIDE RATIO W/ OPPSIDE MEAS
  //TH1D *h_recratio_ssunf = (TH1D*)hunf_sameside->Clone( "ppMC_Meas_Ratio_SSUnf" );
  ////h_recratio_ssunf=(TH1D*)h_recratio_ssunf->Rebin(nbins_pt_reco, "ppMC_Meas_Ratio_SSUnf_rebin" , boundaries_pt_reco);
  //h_recratio_ssunf->SetTitle( "Unf. Py8/RECO Data" );
  //h_recratio_ssunf->Divide(hrec_rebin_ratiobin);
  //if(debugMode)h_recratio_ssunf->Print("base");

  
