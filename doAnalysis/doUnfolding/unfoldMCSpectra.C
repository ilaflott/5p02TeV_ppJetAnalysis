#include "unfoldSpectra.h"

const int kRegDraw  = 4 ; // array entries w/ arguments 0-8. 4 -> middle hist on 3x3 SVDplot

const bool drawPDFs=true; 
//const bool drawPDFs_BayesInputHistos= doBayes && drawPDFs;
const bool drawPDFs_SVDInputHistos  = doSVD   && drawPDFs;

//other settings
const bool doToyErrs=true;
const bool debugMode=true;
const bool debugPearson=(false && debugMode) ;


// CODE --------------------------------------------------
int unfoldMCSpectra( std::string inFile_MC_dir , const std::string baseName , 
		     const bool doJetID=true , const int kRegCenter= 8 ){
  


  // BINNING -----------
  double* boundaries_pt_reco ;   int nbins_pt_reco ;   // reco spectra bins
  double* boundaries_pt_gen  ;   int nbins_pt_gen  ;   // gen spectra bins

  double* boundaries_pt_reco_mat ;   int nbins_pt_reco_mat ;   // matrix bins
  double* boundaries_pt_gen_mat  ;   int nbins_pt_gen_mat  ;   // 

  int Nbins2Clear_gen_lowpt=0;
  //int Nbins2Clear_gen_highpt=0;
  int Nbins2Clear_reco_lowpt=0;
  int Nbins2Clear_reco_highpt=0;




  if(useSimplePtBinning){
    std::cout<<"using simple pt bins"<<std::endl<<std::endl;
    boundaries_pt_reco = (double*)simpbins_pt_reco ; //spectra
    nbins_pt_reco      = (int)n_simpbins_pt_reco ; 
    boundaries_pt_gen = (double*)simpbins_pt_gen ;
    nbins_pt_gen      = (int)n_simpbins_pt_gen ; 

    boundaries_pt_reco_mat = (double*)simpbins_pt_reco ; //matrix
    nbins_pt_reco_mat      = (int)n_simpbins_pt_reco ; 
    boundaries_pt_gen_mat = (double*)simpbins_pt_gen ;
    nbins_pt_gen_mat      = (int)n_simpbins_pt_gen ;                    
    Nbins2Clear_gen_lowpt=0;
    //Nbins2Clear_gen_highpt=0;
    Nbins2Clear_reco_lowpt=0;
    Nbins2Clear_reco_highpt=0;

} 

  else{
    std::cout<<"using analysis pt bins"<<std::endl<<std::endl;
    boundaries_pt_reco = (double*)anabins_pt_reco ; //spectra
    nbins_pt_reco      = (int)n_anabins_pt_reco ; 
    boundaries_pt_gen = (double*)anabins_pt_gen ;
    nbins_pt_gen      = (int)n_anabins_pt_gen ; 

    boundaries_pt_reco_mat = (double*)anabins_pt_reco ; //matrix
    nbins_pt_reco_mat      = (int)n_anabins_pt_reco ; 
    boundaries_pt_gen_mat = (double*)anabins_pt_gen ;
    nbins_pt_gen_mat      = (int)n_anabins_pt_gen ;                  
    Nbins2Clear_gen_lowpt=0;
    //Nbins2Clear_gen_highpt=0;
    Nbins2Clear_reco_lowpt=0;
    Nbins2Clear_reco_highpt=0;
  }



  // STRINGS -----------
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl; 
  //inFile_Data_dir=SCRATCH_BASE+inFile_Data_dir;
  inFile_MC_dir  =SCRATCH_BASE+inFile_MC_dir;
  
  std::size_t radPos=inFile_MC_dir.find("_ak")+3;
  const std::string radiusInt= inFile_MC_dir.substr( radPos,1 );
  const std::string radius="R"+radiusInt;//R4
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;
  const std::string RandEtaRange="_"+radius+"_20_eta_20";//_R4_20_eta_20
  if(debugMode)std::cout<<"RandEtaRange string is = "<<RandEtaRange<<std::endl;
  
  const std::string Rstring="_"+radius; //_R4
  const std::string Rstring_plotTitle=" "+radius; // R4
  const std::string RandEtaRange_plotTitle=Rstring_plotTitle+" 20eta20";// R4 20eta20
  
  std::size_t jetTypePos=radPos+1, jetsPos=inFile_MC_dir.find("Jets");
  const std::string jetType=inFile_MC_dir.substr( jetTypePos,(jetsPos-jetTypePos) );//PFJets
  if(debugMode)std::cout<<"jetType string is = "<<jetType<<std::endl;

  const std::string fullJetType="ak"+radiusInt+jetType;//"ak4PFJets"
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  


  // INFILE NAME(S) -----------
  const std::string inFile_MC_name="/Py8_CUETP8M1_QCDjetAllPtBins_"+fullJetType+"-allFiles.root";
  //const std::string inFile_Data_name="/HighPtJetTrig_"+fullJetType+"-allFiles.root";

  // OUTPUT FILE, NAME(S) -----------
  std::string outFileName=unfoldMCSpectra_outdir+fullJetType;//+"_"+baseName;//+".root";  
  if(doJetID)outFileName+="_wjtID";//+".root";  
  if(!useSimplePtBinning)outFileName+="_anabins";
  else outFileName+="_simpbins";
  outFileName+="_"+baseName;

  std::string outRespMatPdfFile =  outFileName+"_respMat.pdf";
  //std::string outInputSpectraPdfFile =  outFileName+"_inputSpectra.pdf";
  std::string outBayesPdfFile =  outFileName+"_Bayes.pdf";
  std::string outSVDPdfFile   =  outFileName+"_SVD.pdf"; // see drawPDFs part for rest of string
  std::string outRootFile     =  outFileName+"_File.root";  
  
  
  if(debugMode)std::cout<<"opening output file: "<<outRootFile<<std::endl;
  TFile* fout = new TFile(outRootFile.c_str(),"RECREATE"); 



  // ERROR/WEIGHTS/STATS HANDLING ------------------
  RooUnfold::ErrorTreatment errorTreatment;
  if(!doToyErrs) errorTreatment = RooUnfold::kCovariance;
  else errorTreatment = RooUnfold::kCovToy; 

  if(debugMode)std::cout<<"doToyErrs="<<doToyErrs<<std::endl; 
   
  std::cout<<"TH2 GetDefaultSumw2="<<TH2::GetDefaultSumw2()<<std::endl;
  std::cout<<"setting TH2DefSumw2 to true..."<<std::endl; 
  TH2::SetDefaultSumw2(true);
  
  std::cout<<"TH1 GetDefaultSumw2="<<TH1::GetDefaultSumw2()<<std::endl;
  std::cout<<"setting TH1DefSumw2 to true..."<<std::endl; 
  TH1::SetDefaultSumw2(true);
  
  gStyle->SetOptStat(0);



  // ppMC input histos -------------------------
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from MC file"<<std::endl; 
  if(debugMode){
    std::cout<<"input MC dir : "<< (inFile_MC_dir)  <<std::endl; 
    std::cout<<"MC file name : "<< (inFile_MC_name)<<std::endl; 
    std::cout<<std::endl<<std::endl;  }
  TFile *fpp_MC = TFile::Open( (inFile_MC_dir+inFile_MC_name).c_str());
  


  // ---------- reco, measured spectra to unfold
  std::string histTitle="hpp_mcclosure_reco_test";
  if(doJetID)histTitle+="_wJetID";
  histTitle+=RandEtaRange;
  
  TH1F *hrec;
  hrec = (TH1F*)fpp_MC->Get( histTitle.c_str() ); 

  fout->cd();
  hrec->Write("hpp_mcclosure_reco_test_forChecking");
  hrec->Print("base");
  //std::cout<<"GetSumw2N="<<hrec->TH1::GetSumw2N()<<std::endl<<std::endl;

  TH1F *hrec_anabin;
  hrec_anabin = (TH1F*)hrec->Clone( (histTitle+"_clone").c_str() );
  //hrec_anabin->Scale(1/effIntgrtdLumi_vz);
  hrec_anabin->Print("base");
  //std::cout<<"GetSumw2N="<<hrec_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;

  if(debugMode)std::cout<<"rebinning hrec..."<<std::endl;
  hrec_anabin = (TH1F*)hrec_anabin->Rebin(nbins_pt_reco, (histTitle+"_anabins").c_str() , boundaries_pt_reco); 
  hrec_anabin->Write("hpp_mcclosure_reco_test_rebinned");
  hrec_anabin->Print("base");  
  //std::cout<<"GetSumw2N="<<hrec_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;
  
  divideBinWidth(hrec_anabin); 
  hrec_anabin->Write("hpp_mcclosure_reco_test_rebinned_divBinWidth");
  hrec_anabin->Print("base");  
  //std::cout<<"GetSumw2N="<<hrec_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;
  
  if(clearOverUnderflows){
    TH1clearOverUnderflows((TH1*)hrec_anabin);
    hrec_anabin->Write("hpp_mcclosure_reco_test_rebinned_divBinWidth_noOverUnderFlows");
    hrec_anabin->Print("base");  
    //std::cout<<"GetSumw2N="<<hrec_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;    
  }

  if(zeroBins){
    TH1zeroBins((TH1*)hrec_anabin,Nbins2Clear_reco_lowpt, Nbins2Clear_reco_highpt);
    hrec_anabin->Write("hpp_mcclosure_reco_test_rebinned_divBinWidth_noOverUnderFlows_zeroBins");
    hrec_anabin->Print("base");  
    //std::cout<<"GetSumw2N="<<hrec_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;    
  }
  
  
  
  // response hist, for output? what is this for if it's empty?
  TH1F* hrec_resp_anabin;
  if(fillRespHists) hrec_resp_anabin = (TH1F*)hrec_anabin->Clone("recanabinClone4unf");
  else{
    hrec_resp_anabin = new TH1F( ("hpp_rec_response_anabin"+RandEtaRange).c_str(),"", 
				 nbins_pt_reco, boundaries_pt_reco);}
  hrec_resp_anabin->Print(" base");  
  //std::cout<<"GetSumw2N="<<hrec_resp_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;    

  

  // ---------- reco, measured spectra used to create response matrix, and for "sameside" unfolding test
  std::string histTitle2="hpp_mcclosure_reco";
  if(doJetID)histTitle2+="_wJetID";
  histTitle2+=RandEtaRange;
  
  TH1F *hrec_sameside;
  hrec_sameside = (TH1F*)fpp_MC->Get( histTitle2.c_str() ); 
  hrec_sameside->Print("base");  
  //std::cout<<"GetSumw2N="<<hrec_sameside->TH1::GetSumw2N()<<std::endl<<std::endl;    
  
  TH1F *hrec_sameside_anabin;
  hrec_sameside_anabin = (TH1F*)hrec_sameside->Clone( (histTitle2+"_clone").c_str() );
  //hrec_sameside_anabin->Scale(1/effIntgrtdLumi_vz);
  hrec_sameside_anabin->Print("base");
  //std::cout<<"GetSumw2N="<<hrec_sameside_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;    
  
  if(debugMode)std::cout<<"rebinning hrec_sameside..."<<std::endl;
  hrec_sameside_anabin = (TH1F*)hrec_sameside_anabin->Rebin(nbins_pt_reco, (histTitle2+"_anabins").c_str() , boundaries_pt_reco); 
  hrec_sameside_anabin->Print("base");
  //std::cout<<"GetSumw2N="<<hrec_sameside_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;    
  
  divideBinWidth(hrec_sameside_anabin); 
  hrec_sameside_anabin->Print("base");
  //std::cout<<"GetSumw2N="<<hrec_sameside_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;    
  
  if(clearOverUnderflows){
    TH1clearOverUnderflows((TH1*)hrec_sameside_anabin);
    hrec_sameside_anabin->Print("base");
    //std::cout<<"GetSumw2N="<<hrec_sameside_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;   
  }

  if(zeroBins){
    TH1zeroBins((TH1*)hrec_sameside_anabin,Nbins2Clear_reco_lowpt, Nbins2Clear_reco_highpt);
    hrec_sameside_anabin->Print("base");  
    //std::cout<<"GetSumw2N="<<hrec_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;    
  }

  if(doMCIntegralScaling){
    std::cout<<"scaling hrec_anabin to hrec_sameside_anabin"<<std::endl;
    float integral_rec_anabin=hrec_anabin->Integral();
    float integral_rec_ss_anabin=hrec_sameside_anabin->Integral();
    std::cout<<"scale factor="<<integral_rec_ss_anabin/integral_rec_anabin<<std::endl;
    hrec_anabin->Scale(integral_rec_ss_anabin/integral_rec_anabin);
  }
  


  // ---------- gen, MC truth spectra
  std::string genHistTitle="hpp_mcclosure_gen";
  if(doJetID)genHistTitle+="_wJetID";
  genHistTitle+=RandEtaRange;
  
  TH1F* hgen = (TH1F*)fpp_MC->Get( genHistTitle.c_str() );
  hgen->Print("base");    
  //std::cout<<"GetSumw2N="<<hgen->TH1::GetSumw2N()<<std::endl<<std::endl;    
  
  TH1F* hgen_anabin = (TH1F*)hgen->Clone( (genHistTitle+"_clone").c_str() );
  hgen_anabin = (TH1F*)hgen_anabin->Rebin(nbins_pt_gen, (genHistTitle+"_anabins").c_str() , boundaries_pt_gen);
  hgen_anabin->Print("base"); 
  //std::cout<<"GetSumw2N="<<hgen_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;    
  
  divideBinWidth(hgen_anabin);
  hgen_anabin->Print("base");  
  //std::cout<<"GetSumw2N="<<hgen_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;    
  
  if(clearOverUnderflows){
    TH1clearOverUnderflows((TH1*)hgen_anabin);
    hgen_anabin->Print("base");  
    //std::cout<<"GetSumw2N="<<hgen_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;      
  }

  //if(zeroBins){
  //TH1zeroBins((TH1*)hgen_anabin, Nbins2Clear_gen);
  //hgen_anabin->Print("base");  
  ////std::cout<<"GetSumw2N="<<hrec_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;    
  //}
  
  TH1F* hgen_resp_anabin;
  if(fillRespHists) hgen_resp_anabin = (TH1F*)hgen_anabin->Clone("genanabinClone4unf");
  else{
    hgen_resp_anabin = new TH1F( ("hpp_gen_response_anabin"+RandEtaRange).c_str() ,"", 
				 nbins_pt_gen, boundaries_pt_gen);
    //  hgen_resp_anabin->Sumw2(); 
  }
  hgen_resp_anabin->Print("base");  
  //std::cout<<"GetSumw2N="<<hgen_resp_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;      
  



  std::cout<<std::endl<<"writing input hists to file..."<<std::endl;
  
  fout->cd();
  hgen->Write(); 
  hrec->Write();
  hgen_anabin->Write(); 
  hrec_anabin->Write(); 
  
  hrec_sameside->Write();
  hrec_sameside_anabin->Write(); 
  
  if(fillRespHists) hgen_resp_anabin->Write();
  if(fillRespHists) hrec_resp_anabin->Write();
  


  
  
  // ---------- matrix, MC "response"
  std::string TH2_title="hpp_mcclosure_matrix";//+RandEtaRange;
  if(doJetID)TH2_title+="_wJetID";
  TH2_title+=RandEtaRange;
  
  TH2F* hmat = (TH2F*)fpp_MC->Get( TH2_title.c_str() );

  fout->cd();
  hmat->Write();
  hmat->Print("base");
  //std::cout<<"GetSumw2N="<<hmat->TH1::GetSumw2N()<<std::endl<<std::endl;      


  TH2F* hmat_anabin = (TH2F*)hmat->Clone( (TH2_title+"_clone").c_str() );
  hmat_anabin->Print("base"); 
  //std::cout<<"GetSumw2N="<<hmat_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;      

  hmat_anabin=(TH2F*) reBinTH2(hmat_anabin, (TH2_title+"_anabins").c_str(), 
			       (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
			       (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
			       //			       (double*) boundaries_pt, (int) nbins_pt );
  hmat_anabin->Write("hmat_rebinned");
  hmat_anabin->Print("base"); 
  //std::cout<<"GetSumw2N="<<hmat_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;      

  
  divideBinWidth_TH2(hmat_anabin);
  hmat_anabin->Write("hmat_rebinned_binWidthDiv");
  hmat_anabin->Print("base"); 
  //std::cout<<"GetSumw2N="<<hmat_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;      

  //if(!doOverUnderflows){
  if(clearOverUnderflows){
    TH2clearOverUnderflows((TH2F*)hmat_anabin);
    hmat_anabin->Write("hmat_rebinned_binWidthDiv_noOverUnderFlow");
    hmat_anabin->Print("base");
    //std::cout<<"GetSumw2N="<<hmat_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;      
  }
  
  if(normalizedMCMatrix){
    normalizeMC_TH2(hmat_anabin);
    hmat_anabin->Write("hmat_rebinned_binWidthDiv_noOverUnderFlow_normd");
    hmat_anabin->Print("base");
    //std::cout<<"GetSumw2N="<<hmat_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;      
  }
  
  TH2F* hmat_errors=makeRespMatrixErrors( (TH2F*) hmat,
					  (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
					  (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
  hmat_errors->Print("base");

  divideBinWidth_TH2(hmat_errors);
  hmat_errors->Print("base");
  
  if(clearOverUnderflows){
    TH2clearOverUnderflows((TH2F*)hmat_errors);
    hmat_errors->Print("base");
    //std::cout<<"GetSumw2N="<<hmat_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;      
  }

  setRespMatrixErrs( (TH2F*)hmat_anabin, (TH2F*) hmat_errors );
  
  
  TH2F* hmat_percenterrs=makeRespMatrixPercentErrs( (TH2F*) hmat_errors, (TH2F*) hmat_anabin,
						    (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
						    (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
  hmat_percenterrs->Print("base");
  
  
  
  if(drawPDFs){    
    
    bool drawRespMatrix=true;
    //bool drawInputSpectra=true;
    
    if(drawRespMatrix){
      std::cout<<std::endl<<"drawing input response matrices..."<<std::endl;
      
      std::string outPdfFile=outRespMatPdfFile;
      std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";
      
      TCanvas* tempCanvForPdfPrint=new TCanvas("tempCanv_respMat","",1200,1200);    
      tempCanvForPdfPrint->cd();
      
      tempCanvForPdfPrint->Print(open_outPdfFile.c_str()); 
      

      // matrix ---------------
    
      if(normalizedMCMatrix) {
	tempCanvForPdfPrint->SetLogz(0);

	hmat_anabin->SetAxisRange(0. , 1,"Z");
	hmat_anabin->SetTitle("ppMC rebinned input, normalized");

	hmat->SetAxisRange(0. , 1,"Z");
	hmat->SetTitle("ppMC input, normalized");
      }
      else	{
	tempCanvForPdfPrint->SetLogz(1);

	hmat_anabin->SetAxisRange(0.0000000000000001,.00001,"Z");
	hmat_anabin->SetTitle("ppMC rebinned input. not normalized");

	hmat->SetAxisRange(0.0000000000000001,.00001,"Z");
	hmat->SetTitle("ppMC input. not normalized");
      }
      
      
      tempCanvForPdfPrint->cd();

      tempCanvForPdfPrint->SetLogx(0);
      tempCanvForPdfPrint->SetLogy(0);

      hmat->SetTitle("ppMC Resp Matrix, original");

      hmat->SetAxisRange(0.,1000.,"X");
      hmat->SetAxisRange(0.,1000.,"Y");

      
      hmat->GetZaxis()->SetLabelSize(0.025);
      
      hmat->GetYaxis()->SetLabelSize(0.02);
      hmat->GetYaxis()->SetTitleSize(0.023);
      hmat->GetYaxis()->SetTitle("gen p_{t}");
      
      hmat->GetXaxis()->SetLabelSize(0.02);
      hmat->GetXaxis()->SetTitleSize(0.025);
      hmat->GetXaxis()->SetTitle("reco p_{t}   ");
      hmat->Draw("COLZ");           
      
      tempCanvForPdfPrint->Print(outPdfFile.c_str());


      // matrix w/ specific range ---------------
    
      tempCanvForPdfPrint->cd();

      tempCanvForPdfPrint->SetLogx(0);
      tempCanvForPdfPrint->SetLogy(0);
      tempCanvForPdfPrint->SetLogz(1);

      hmat->SetTitle("ppMC Resp Matrix, no-log axes");

      hmat->SetAxisRange(boundaries_pt_reco_mat[0],boundaries_pt_reco_mat[nbins_pt_reco_mat],"X");
      hmat->SetAxisRange(boundaries_pt_gen_mat[0],boundaries_pt_gen_mat[nbins_pt_gen_mat],"Y");
      
      hmat->GetZaxis()->SetLabelSize(0.025);
      
      hmat->GetYaxis()->SetLabelSize(0.02);
      hmat->GetYaxis()->SetTitleSize(0.023);
      hmat->GetYaxis()->SetTitle("gen p_{t}");
      
      hmat->GetXaxis()->SetLabelSize(0.02);
      hmat->GetXaxis()->SetTitleSize(0.025);
      hmat->GetXaxis()->SetTitle("reco p_{t}   ");
      hmat->Draw("COLZ");           
      
      tempCanvForPdfPrint->Print(outPdfFile.c_str());
      

      // matrix w/ log ---------------

      tempCanvForPdfPrint->cd();

      tempCanvForPdfPrint->SetLogx(1);
      tempCanvForPdfPrint->SetLogy(1);
      tempCanvForPdfPrint->SetLogz(1);

      hmat->SetTitle("ppMC Resp Matrix, w/ log axes");

      //hmat->SetAxisRange(boundaries_pt_reco_mat[0],boundaries_pt_reco_mat[nbins_pt_reco_mat],"X");
      //hmat->SetAxisRange(boundaries_pt_gen_mat[0],boundaries_pt_gen_mat[nbins_pt_gen_mat],"Y");

      //hmat->GetZaxis()->SetLabelSize(0.025);      
      hmat->GetYaxis()->SetMoreLogLabels(true);
      hmat->GetYaxis()->SetNoExponent(true);
      //hmat->GetYaxis()->SetLabelSize(0.02);
      //hmat->GetYaxis()->SetTitleSize(0.025);
      //hmat->GetYaxis()->SetTitle("gen p_{t}");
      //
      hmat->GetXaxis()->SetMoreLogLabels(true);
      hmat->GetXaxis()->SetNoExponent(true);
      //hmat->GetXaxis()->SetLabelSize(0.02);
      //hmat->GetXaxis()->SetTitleSize(0.025);
      //hmat->GetXaxis()->SetTitle("reco p_{t}   ");

      hmat->Draw("COLZ");           
      
      tempCanvForPdfPrint->Print(outPdfFile.c_str());

      // matrix rebinned ---------------
      
      tempCanvForPdfPrint->cd();
      
      tempCanvForPdfPrint->SetLogx(0);
      tempCanvForPdfPrint->SetLogy(0);
      tempCanvForPdfPrint->SetLogz(1);

      hmat_anabin->SetTitle("ppMC Resp Matrix rebin, no log axes");
      
      hmat_anabin->GetZaxis()->SetLabelSize(0.025);
      
      hmat_anabin->GetYaxis()->SetMoreLogLabels(true);
      hmat_anabin->GetYaxis()->SetNoExponent(true);
      hmat_anabin->GetYaxis()->SetLabelSize(0.02);
      hmat_anabin->GetYaxis()->SetTitleSize(0.025);
      hmat_anabin->GetYaxis()->SetTitle("gen p_{t}");
      
      hmat_anabin->GetXaxis()->SetLabelSize(0.02);
      hmat_anabin->GetXaxis()->SetTitleSize(0.025);
      hmat_anabin->GetXaxis()->SetTitle("reco p_{t}   ");
      hmat_anabin->Draw("COLZ");           
      
      tempCanvForPdfPrint->Print(outPdfFile.c_str());


      // matrix rebinned ---------------
      
      tempCanvForPdfPrint->cd();
      
      tempCanvForPdfPrint->SetLogx(1);
      tempCanvForPdfPrint->SetLogy(1);
      tempCanvForPdfPrint->SetLogz(1);

      hmat_anabin->SetTitle("ppMC Resp Matrix rebin, w/ log axes");
      
      hmat_anabin->GetZaxis()->SetLabelSize(0.025);
      
      hmat_anabin->GetYaxis()->SetMoreLogLabels(true);
      hmat_anabin->GetYaxis()->SetNoExponent(true);
      hmat_anabin->GetYaxis()->SetLabelSize(0.02);
      hmat_anabin->GetYaxis()->SetTitleSize(0.025);
      hmat_anabin->GetYaxis()->SetTitle("gen p_{t}");
      
      hmat_anabin->GetXaxis()->SetMoreLogLabels(true);
      hmat_anabin->GetXaxis()->SetNoExponent(true);
      hmat_anabin->GetXaxis()->SetLabelSize(0.02);
      hmat_anabin->GetXaxis()->SetTitleSize(0.025);
      hmat_anabin->GetXaxis()->SetTitle("reco p_{t}   ");
      hmat_anabin->Draw("COLZ");           
      
      tempCanvForPdfPrint->Print(outPdfFile.c_str());


      // error matrix in binning of interest ---------------
      
      tempCanvForPdfPrint->cd();
      
      tempCanvForPdfPrint->SetLogx(1);
      tempCanvForPdfPrint->SetLogy(1);
      tempCanvForPdfPrint->SetLogz(1);

      hmat_errors->SetTitle("Resp Matrix Errors, log axes");

      hmat_errors->SetAxisRange(0.0000000000000001,.00001,"Z");      
      hmat_errors->GetZaxis()->SetLabelSize(0.025);
      
      hmat_errors->GetYaxis()->SetMoreLogLabels(true);
      hmat_errors->GetYaxis()->SetNoExponent(true);
      hmat_errors->GetYaxis()->SetLabelSize(0.02);
      hmat_errors->GetYaxis()->SetTitleSize(0.025);
      hmat_errors->GetYaxis()->SetTitle("gen p_{t}");
      
      hmat_errors->GetXaxis()->SetMoreLogLabels(true);
      hmat_errors->GetXaxis()->SetNoExponent(true);
      hmat_errors->GetXaxis()->SetLabelSize(0.02);
      hmat_errors->GetXaxis()->SetTitleSize(0.025);
      hmat_errors->GetXaxis()->SetTitle("reco p_{t}   ");
      hmat_errors->Draw("COLZ");           
      
      tempCanvForPdfPrint->Print(outPdfFile.c_str());

      // percent error matrix in binning of interest , no log ---------------
      
      tempCanvForPdfPrint->cd();
      
      tempCanvForPdfPrint->SetLogx(0);
      tempCanvForPdfPrint->SetLogy(0);
      tempCanvForPdfPrint->SetLogz(1);

      hmat_percenterrs->SetAxisRange(0.1,100.,"Z");      

      hmat_percenterrs->SetTitle("Resp Matrix % Errors non-log x/y axes");
      
      hmat_percenterrs->GetZaxis()->SetLabelSize(0.025);
      hmat_percenterrs->GetYaxis()->SetLabelSize(0.02);
      hmat_percenterrs->GetYaxis()->SetTitleSize(0.025);

      hmat_percenterrs->GetYaxis()->SetTitle("gen p_{t}");
      
      hmat_percenterrs->GetXaxis()->SetLabelSize(0.02);
      hmat_percenterrs->GetXaxis()->SetTitleSize(0.025);

      hmat_percenterrs->GetXaxis()->SetTitle("reco p_{t}   ");
      hmat_percenterrs->Draw("COLZ");           
      
      tempCanvForPdfPrint->Print(outPdfFile.c_str());



      // percent error matrix in binning of interest , no log ---------------
      
      tempCanvForPdfPrint->cd();
      
      tempCanvForPdfPrint->SetLogx(1);
      tempCanvForPdfPrint->SetLogy(1);
      tempCanvForPdfPrint->SetLogz(1);

      hmat_percenterrs->SetAxisRange(0.1,100.,"Z");      

      hmat_percenterrs->SetTitle("Resp Matrix % Errors w/ log-axes");
      
      hmat_percenterrs->GetZaxis()->SetLabelSize(0.025);
      
      hmat_percenterrs->GetYaxis()->SetMoreLogLabels(true);
      hmat_percenterrs->GetYaxis()->SetNoExponent(true);
      hmat_percenterrs->GetYaxis()->SetLabelSize(0.02);
      hmat_percenterrs->GetYaxis()->SetTitleSize(0.025);
      hmat_percenterrs->GetYaxis()->SetTitle("gen p_{t}");
      
      hmat_percenterrs->GetXaxis()->SetMoreLogLabels(true);
      hmat_percenterrs->GetXaxis()->SetNoExponent(true);
      hmat_percenterrs->GetXaxis()->SetLabelSize(0.02);
      hmat_percenterrs->GetXaxis()->SetTitleSize(0.025);
      hmat_percenterrs->GetXaxis()->SetTitle("reco p_{t}   ");
      hmat_percenterrs->Draw("COLZ");           
      
      tempCanvForPdfPrint->Print(outPdfFile.c_str());


      // % error of matrix rebinned ---------------
      
//      tempCanvForPdfPrint->cd();
//      
//      tempCanvForPdfPrint->SetLogx(1);
//      tempCanvForPdfPrint->SetLogy(1);
//      tempCanvForPdfPrint->SetLogz(0);
      
//TH1F* hmat_anabin_errs=(TH1F*)respMatrixError(hmat_anabin,
//						    boundaries_pt_reco_mat , nbins_pt_reco_mat,
//						    boundaries_pt_gen_mat  , nbins_pt_gen_mat);
//
//hmat_anabin->GetZaxis()->SetLabelSize(0.025);
//
//hmat_anabin->GetYaxis()->SetMoreLogLabels(true);
//hmat_anabin->GetYaxis()->SetNoExponent(true);
//hmat_anabin->GetYaxis()->SetLabelSize(0.02);
//hmat_anabin->GetYaxis()->SetTitleSize(0.025);
//hmat_anabin->GetYaxis()->SetTitle("gen p_{t}");
//
//hmat_anabin->GetXaxis()->SetMoreLogLabels(true);
//hmat_anabin->GetXaxis()->SetNoExponent(true);
//hmat_anabin->GetXaxis()->SetLabelSize(0.02);
//hmat_anabin->GetXaxis()->SetTitleSize(0.025);
//hmat_anabin->GetXaxis()->SetTitle("reco p_{t}   ");
//hmat_anabin->Draw("COLZ");           
//      
//      tempCanvForPdfPrint->Print(outPdfFile.c_str());
      

      // close file ---------------

      tempCanvForPdfPrint->Print(close_outPdfFile.c_str()); }
  }
  













  // Bayesian unfolding ------------------------- 
  if(!doBayes)  std::cout<<std::endl<<"   skipping Bayesian Unfolding..."<<std::endl<<std::endl;
  else{ std::cout<<std::endl<<std::endl<<std::endl<<"   beginning Bayesian Unfolding..."<<std::endl;
    

    // RooUnfoldResponse, Bayes --------------------
    std::cout<<"calling RooUnfoldResponse "<<std::endl;
    RooUnfoldResponse roo_resp( hrec_resp_anabin, hgen_resp_anabin, hmat_anabin, ("Response_matrix"+RandEtaRange).c_str()) ;

    roo_resp.UseOverflow(doOverUnderflows);    
        
    std::cout<<"calling RooUnfoldBayes..."<<std::endl;
    RooUnfoldBayes unf_bayes( &roo_resp, hrec_anabin, kIter );

    TH1F *hunf = (TH1F*)unf_bayes.Hreco(errorTreatment);     std::cout<<std::endl; 

    hunf->SetName("ppMC_BayesUnf_Spectra");
    hunf->SetTitle("ppMC BayesUnf Spectra");
    hunf->Print("base");
    //std::cout<<"GetSumw2N="<<hunf->TH1::GetSumw2N()<<std::endl<<std::endl ;
     
    TH1F *hfold = (TH1F*)roo_resp.ApplyToTruth(hunf);    std::cout<<std::endl; 
    hfold->SetName("ppMC_BayesFold_Spectra");
    hfold->SetTitle("ppMC BayesFold Spectra");
    hfold->Print("base");
    //std::cout<<"GetSumw2N="<<hfold->TH1::GetSumw2N()<<std::endl<<std::endl; 

    RooUnfoldBayes unf_sameside_bayes( &roo_resp, hrec_sameside_anabin, kIter );

    TH1F *hunf_sameside = (TH1F*)unf_sameside_bayes.Hreco(errorTreatment);    std::cout<<std::endl; 
    hunf_sameside->SetName("ppMC_BayesUnf_sameSideSpectra");
    hunf_sameside->SetTitle("ppMC BayesUnf sameSideSpectra");
    hunf_sameside->Print("base");
    //std::cout<<"GetSumw2N="<<hunf_sameside->TH1::GetSumw2N()<<std::endl<<std::endl; 


    TH1F *hfold_sameside = (TH1F*)roo_resp.ApplyToTruth(hunf_sameside);    std::cout<<std::endl; 
    hfold_sameside->SetName("ppMC_BayesFold_sameSide_Spectra");
    hfold_sameside->SetTitle("ppMC BayesFold sameSide Spectra");
    hfold_sameside->Print("base");
    //std::cout<<"GetSumw2N="<<hfold_sameside->TH1::GetSumw2N()<<std::endl<<std::endl; 


    // Gen Ratio Plots --------------------
    TH1F *hgen_anabin_ratiobin=(TH1F*)hgen_anabin->Clone("ppMC_Gen_Ratio_denom");
    hgen_anabin_ratiobin->Print("base");
    //std::cout<<"GetSumw2N="<<hgen_anabin_ratiobin->TH1::GetSumw2N()<<std::endl<<std::endl; 

    hgen_anabin_ratiobin=(TH1F*)hgen_anabin_ratiobin->Rebin(nbins_pt_reco,"ppMC_Gen_Ratio_denom_rebin",boundaries_pt_reco);
    hgen_anabin_ratiobin->Print("base");
    //std::cout<<"GetSumw2N="<<hgen_anabin_ratiobin->TH1::GetSumw2N()<<std::endl<<std::endl; 


    TH1F *h_genratio_oppunf = (TH1F*)hunf->Clone( "ppMC_Gen_Ratio_OppUnf" );
    h_genratio_oppunf->Print("base");
    //std::cout<<"GetSumw2N="<<h_genratio_oppunf->TH1::GetSumw2N()<<std::endl<<std::endl; 

    h_genratio_oppunf->SetTitle( "Bayes Unf./Gen. Truth" );
    h_genratio_oppunf->SetMarkerStyle(24);
    h_genratio_oppunf->SetMarkerColor(2);
    h_genratio_oppunf->Divide(hgen_anabin);

    h_genratio_oppunf->Print("base");
    //std::cout<<"GetSumw2N="<<h_genratio_oppunf->TH1::GetSumw2N()<<std::endl<<std::endl; 


    TH1F *h_genratio_oppfold = (TH1F*)hfold->Clone( "ppMC_Gen_Ratio_OppFold" );
    h_genratio_oppfold->Print("base");
    //std::cout<<"GetSumw2N="<<h_genratio_oppfold->TH1::GetSumw2N()<<std::endl<<std::endl; 

    h_genratio_oppfold->SetTitle( "Bayes Fold/Gen. Truth" );
    h_genratio_oppfold->SetMarkerStyle(24);
    h_genratio_oppfold->SetMarkerColor(3);
    h_genratio_oppfold->Divide(hgen_anabin_ratiobin);

    h_genratio_oppfold->Print("base");
    //std::cout<<"GetSumw2N="<<h_genratio_oppfold->TH1::GetSumw2N()<<std::endl<<std::endl; 


    TH1F *h_genratio_oppmeas = (TH1F*)hrec_anabin->Clone( "ppMC_Gen_Ratio_Meas" );
    h_genratio_oppmeas->Print("base");
    //std::cout<<"GetSumw2N="<<h_genratio_oppmeas->TH1::GetSumw2N()<<std::endl<<std::endl; 

    h_genratio_oppmeas->SetTitle( "Meas./Gen. Truth" );
    h_genratio_oppmeas->SetMarkerStyle(24);
    h_genratio_oppmeas->SetMarkerColor(4);
    h_genratio_oppmeas->Divide(hgen_anabin_ratiobin);

    h_genratio_oppmeas->Print("base");
    //std::cout<<"GetSumw2N="<<h_genratio_oppmeas->TH1::GetSumw2N()<<std::endl<<std::endl; 


    TH1F *h_genratio_ssunf = (TH1F*)hunf_sameside->Clone( "ppMC_Gen_Ratio_SSUnf" );

    h_genratio_ssunf->Print("base");
    //std::cout<<"GetSumw2N="<<h_genratio_ssunf->TH1::GetSumw2N()<<std::endl<<std::endl; 

    h_genratio_ssunf->SetTitle( "SameSide Bayes Unf./Gen. Truth" );
    h_genratio_ssunf->SetMarkerStyle(25);
    h_genratio_ssunf->SetMarkerColor(2);
    h_genratio_ssunf->Divide(hgen_anabin);

    h_genratio_ssunf->Print("base");
    //std::cout<<"GetSumw2N="<<h_genratio_ssunf->TH1::GetSumw2N()<<std::endl<<std::endl; 


    TH1F *h_genratio_ssfold = (TH1F*)hfold_sameside->Clone( "ppMC_Gen_Ratio_SSFold" );
    h_genratio_ssfold->Print("base");
    //std::cout<<"GetSumw2N="<<h_genratio_ssfold->TH1::GetSumw2N()<<std::endl<<std::endl; 

    h_genratio_ssfold->SetTitle( "SameSide Bayes Fold/Gen. Truth" );
    h_genratio_ssfold->SetMarkerStyle(25);
    h_genratio_ssfold->SetMarkerColor(3);
    h_genratio_ssfold->Divide(hgen_anabin_ratiobin);

    h_genratio_ssfold->Print("base");
    //std::cout<<"GetSumw2N="<<h_genratio_ssfold->TH1::GetSumw2N()<<std::endl<<std::endl; 


    TH1F *h_genratio_ssmeas = (TH1F*)hrec_sameside_anabin->Clone( "ppMC_Gen_Ratio4_SSMeas" );
    h_genratio_ssmeas->Print("base");
    //std::cout<<"GetSumw2N="<<h_genratio_ssmeas->TH1::GetSumw2N()<<std::endl<<std::endl; 

    h_genratio_ssmeas->SetTitle( "Same Side Meas./Gen. Truth" );
    h_genratio_ssmeas->SetMarkerStyle(25);
    h_genratio_ssmeas->SetMarkerColor(4);
    h_genratio_ssmeas->Divide(hgen_anabin_ratiobin);

    h_genratio_ssmeas->Print("base");
    //std::cout<<"GetSumw2N="<<h_genratio_ssmeas->TH1::GetSumw2N()<<std::endl<<std::endl; 




    // Rec Ratio Plots -------------------
    TH1F *h_recratio_oppunf = (TH1F*)hunf->Clone( "ppMC_Meas_Ratio_OppUnf" );
    h_recratio_oppunf->Print("base");
    //std::cout<<"GetSumw2N="<<h_recratio_oppunf->TH1::GetSumw2N()<<std::endl<<std::endl; 

    h_recratio_oppunf=(TH1F*)h_recratio_oppunf->Rebin(nbins_pt_reco, "ppMC_Meas_Ratio_OppUnf_rebin" , boundaries_pt_reco);

    h_recratio_oppunf->Print("base");
    //std::cout<<"GetSumw2N="<<h_recratio_oppunf->TH1::GetSumw2N()<<std::endl<<std::endl; 

    h_recratio_oppunf->SetTitle( "Bayes Unf./Meas." );
    h_recratio_oppunf->SetMarkerStyle(24);
    h_recratio_oppunf->SetMarkerColor(2);
    h_recratio_oppunf->Divide(hrec_anabin);

    h_recratio_oppunf->Print("base");
    //std::cout<<"GetSumw2N="<<h_recratio_oppunf->TH1::GetSumw2N()<<std::endl<<std::endl; 


    TH1F *h_recratio_ssunf = (TH1F*)hunf_sameside->Clone( "ppMC_Meas_Ratio_SSUnf" );
    h_recratio_ssunf->Print("base");
    //std::cout<<"GetSumw2N="<<h_recratio_ssunf->TH1::GetSumw2N()<<std::endl<<std::endl; 

    h_recratio_ssunf=(TH1F*)h_recratio_ssunf->Rebin(nbins_pt_reco, "ppMC_Meas_Ratio_SSUnf_rebin" , boundaries_pt_reco);
    h_recratio_ssunf->Print("base");
    //std::cout<<"GetSumw2N="<<h_recratio_ssunf->TH1::GetSumw2N()<<std::endl<<std::endl; 

    h_recratio_ssunf->SetTitle( "SameSide Bayes Unf./Meas." );
    h_recratio_ssunf->SetMarkerStyle(25);
    h_recratio_ssunf->SetMarkerColor(2);
    h_recratio_ssunf->Divide(hrec_anabin);

    h_recratio_ssunf->Print("base");
    //std::cout<<"GetSumw2N="<<h_recratio_ssunf->TH1::GetSumw2N()<<std::endl<<std::endl; 


    TH1F *h_recratio_oppfold = (TH1F*)hfold->Clone( "ppMC_Meas_Ratio_OppFold" );
    h_recratio_oppfold->Print("base");
    //std::cout<<"GetSumw2N="<<h_recratio_oppfold->TH1::GetSumw2N()<<std::endl<<std::endl; 

    //h_recratio_oppfold=(TH1F*)h_recratio_oppfold->Rebin(nbins_pt_reco, "ppMC_Meas_Ratio_OppFold_rebin" , boundaries_pt_reco);
    //h_recratio_oppfold->Print("base");
    ////std::cout<<"GetSumw2N="<<h_recratio_oppfold->TH1::GetSumw2N()<<std::endl<<std::endl; 

    h_recratio_oppfold->SetTitle( "Bayes Fold./Meas." );
    h_recratio_oppfold->SetMarkerStyle(24);
    h_recratio_oppfold->SetMarkerColor(3);
    h_recratio_oppfold->Divide(hrec_anabin);

    h_recratio_oppfold->Print("base");
    //std::cout<<"GetSumw2N="<<h_recratio_oppfold->TH1::GetSumw2N()<<std::endl<<std::endl; 


    TH1F *h_recratio_ssfold = (TH1F*)hfold_sameside->Clone( "ppMC_Meas_Ratio_SSFold" );
    h_recratio_ssfold->Print("base");
    //std::cout<<"GetSumw2N="<<h_recratio_ssfold->TH1::GetSumw2N()<<std::endl<<std::endl; 

    //h_recratio_ssfold=(TH1F*)h_recratio_ssfold->Rebin(nbins_pt_reco, "ppMC_Meas_Ratio_SSFold_rebin" , boundaries_pt_reco);
    //h_recratio_ssfold->Print("base");
    ////std::cout<<"GetSumw2N="<<h_recratio_ssfold->TH1::GetSumw2N()<<std::endl<<std::endl; 

    h_recratio_ssfold->SetTitle( "SameSide Bayes Fold./Meas." );
    h_recratio_ssfold->SetMarkerStyle(25);
    h_recratio_ssfold->SetMarkerColor(3);
    h_recratio_ssfold->Divide(hrec_anabin);

    h_recratio_ssfold->Print("base");
    //std::cout<<"GetSumw2N="<<h_recratio_ssfold->TH1::GetSumw2N()<<std::endl<<std::endl; 


    TH1F *h_recratio_ssgen = (TH1F*)hgen_anabin->Clone( "ppMC_Meas_Ratio_SSTruth" );
    h_recratio_ssgen->Print("base");
    //std::cout<<"GetSumw2N="<<h_recratio_ssgen->TH1::GetSumw2N()<<std::endl<<std::endl; 

    h_recratio_ssgen=(TH1F*)h_recratio_ssgen->Rebin(nbins_pt_reco, "ppMC_Meas_Ratio_SSTruth_rebin" , boundaries_pt_reco);
    h_recratio_ssgen->Print("base");
    //std::cout<<"GetSumw2N="<<h_recratio_ssgen->TH1::GetSumw2N()<<std::endl<<std::endl; 

    h_recratio_ssgen->SetTitle( "Truth/Meas." );
    h_recratio_ssgen->SetMarkerStyle(24);
    h_recratio_ssgen->SetMarkerColor(6);
    h_recratio_ssgen->Divide(hrec_anabin);

    h_recratio_ssgen->Print("base");
    //std::cout<<"GetSumw2N="<<h_recratio_ssgen->TH1::GetSumw2N()<<std::endl<<std::endl; 


    std::cout<<"writing output hists to file... "<<std::endl;

    hunf->Write();    
    hunf_sameside->Write();
    
    h_genratio_oppunf->Write();  
    h_genratio_oppfold->Write();  
    h_genratio_oppmeas->Write();  
    
    h_genratio_ssunf->Write();  
    h_genratio_ssfold->Write();  
    h_genratio_ssmeas->Write();  
    
    h_recratio_oppunf->Write();  
    h_recratio_oppfold->Write();  

    h_recratio_ssunf->Write();  
    h_recratio_ssfold->Write();  
    h_recratio_ssgen->Write();  

    //  drawPDFS -------------------------------------------------- 
    if(!drawPDFs)std::cout<<std::endl<<"NOT drawing PDFs for Bayesian Unfolding!"<<std::endl<<std::endl;
    else{ std::cout<<std::endl<<"drawing output PDFs for Bayesian Unfolding..."<<std::endl;
      std::string outPdfFile=outBayesPdfFile;
      std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";

      // temp canvas for printing
      TCanvas* tempCanvForPdfPrint_wLogy=new TCanvas("tempCanv_Logy","",1000,1000);    
      tempCanvForPdfPrint_wLogy->cd();
      tempCanvForPdfPrint_wLogy->SetLogy(1);
      tempCanvForPdfPrint_wLogy->SetLogx(1);

      TCanvas* tempCanvForPdfPrint_wLogz=new TCanvas("tempCanv_Logz","",1200,1200);    
      tempCanvForPdfPrint_wLogz->cd();
      //tempCanvForPdfPrint_wLogz->SetLogz(1);
      tempCanvForPdfPrint_wLogz->SetLogx(1);
      tempCanvForPdfPrint_wLogz->SetLogy(1);
      
      TCanvas* tempCanvForPdfPrint=new TCanvas("tempCanv","",           1000,1000);    
      tempCanvForPdfPrint->cd();
      tempCanvForPdfPrint->SetLogy(0);
      tempCanvForPdfPrint->SetLogx(1);
      tempCanvForPdfPrint->Print(open_outPdfFile.c_str()); 
      
      ///////////////

      TLine* recCutLine1= new TLine( boundaries_pt_reco[Nbins2Clear_reco_lowpt],0.,
				     boundaries_pt_reco[Nbins2Clear_reco_lowpt],2.);
      recCutLine1->SetLineWidth(1);
      recCutLine1->SetLineStyle(2);
      recCutLine1->SetLineColor(4);

      TLine* recCutLine2= new TLine( boundaries_pt_reco[nbins_pt_reco-Nbins2Clear_reco_highpt],0.,boundaries_pt_reco[nbins_pt_reco-Nbins2Clear_reco_highpt],2.);
      recCutLine2->SetLineWidth(1);
      recCutLine2->SetLineStyle(2);
      recCutLine2->SetLineColor(4);

      ///////////////

      TLine* genCutLine1= new TLine( boundaries_pt_gen[Nbins2Clear_gen_lowpt],0.,boundaries_pt_gen[Nbins2Clear_gen_lowpt],2.);
      genCutLine1->SetLineWidth(1);
      genCutLine1->SetLineStyle(2);
      genCutLine1->SetLineColor(kGreen+3);

      TLine* genCutLine2= new TLine( boundaries_pt_gen[nbins_pt_gen],0.,boundaries_pt_gen[nbins_pt_gen],2.);
      genCutLine2->SetLineWidth(1);
      genCutLine2->SetLineStyle(2);
      genCutLine2->SetLineColor(kGreen+3);

      ///////////////

      TLine* genMatCutLine1= new TLine( boundaries_pt_gen_mat[0],0.,boundaries_pt_gen_mat[0],2.);
      genMatCutLine1->SetLineWidth(1);
      genMatCutLine1->SetLineStyle(2);
      genMatCutLine1->SetLineColor(3);

      TLine* genMatCutLine2= new TLine( boundaries_pt_gen_mat[nbins_pt_gen_mat],0.,boundaries_pt_gen_mat[nbins_pt_gen_mat],2.);
      genMatCutLine2->SetLineWidth(1);
      genMatCutLine2->SetLineStyle(2);
      genMatCutLine2->SetLineColor(3);

      TLine* recMatCutLine1= new TLine( boundaries_pt_reco_mat[0],0.,boundaries_pt_reco_mat[0],2.);
      recMatCutLine1->SetLineWidth(1);
      recMatCutLine1->SetLineStyle(2);
      recMatCutLine1->SetLineColor(kBlue-7);

      TLine* recMatCutLine2= new TLine( boundaries_pt_reco_mat[nbins_pt_reco_mat],0.,boundaries_pt_reco_mat[nbins_pt_reco_mat],2.);
      recMatCutLine2->SetLineWidth(1);
      recMatCutLine2->SetLineStyle(2);
      recMatCutLine2->SetLineColor(kBlue-7);

      TLine* theLineAtOne= new TLine( boundaries_pt_gen_mat[0],1.,boundaries_pt_gen_mat[nbins_pt_gen_mat],1.);
      theLineAtOne->SetLineWidth(1);
      theLineAtOne->SetLineStyle(2);
      theLineAtOne->SetLineColor(36);

      
      // gen/meas/unf spectra ---------------
      
      tempCanvForPdfPrint_wLogy->cd();
      //tempCanvForPdfPrint_wLogy->SetLogx(1);
      
      hgen_anabin->SetTitle("ppMC Spectra, input and output");
      hgen_anabin->GetXaxis()->SetTitle("jet p_{T} (GeV)");
      hgen_anabin->GetXaxis()->SetMoreLogLabels(true);
      hgen_anabin->GetXaxis()->SetNoExponent(true);
      hgen_anabin->SetAxisRange(boundaries_pt_gen_mat[0], boundaries_pt_gen_mat[nbins_pt_gen_mat],"X");           

      hrec_anabin->SetMarkerStyle(kOpenTriangleUp);
      hrec_anabin->SetMarkerColor(kBlue);     
      hrec_anabin->SetMarkerSize(1.02);     
      
      hunf->SetMarkerStyle(kOpenCircle);
      hunf->SetMarkerColor(kRed);
      hunf->SetMarkerSize(1.02);     
      
      hrec_sameside_anabin->SetMarkerStyle(kOpenTriangleDown);
      hrec_sameside_anabin->SetMarkerColor(kBlue-3);     
      hrec_sameside_anabin->SetMarkerSize(1.02);     
      
      hunf_sameside->SetMarkerStyle(kOpenSquare);
      hunf_sameside->SetMarkerColor(kRed-3);
      hunf_sameside->SetMarkerSize(1.02);     
      
      hgen_anabin->SetMarkerStyle(kOpenCross);
      hgen_anabin->SetMarkerColor(kGreen-2);
      hgen_anabin->SetMarkerSize(1.02);     	
      
      hgen_anabin->Draw();           
      hrec_anabin->Draw("SAME");           
      hunf->Draw("SAME");           
      hunf_sameside->Draw("SAME");           
      hrec_sameside_anabin->Draw("SAME");           
      
      
      TLegend* legend_in = new TLegend( 0.6,0.8,0.9,0.9 );
      legend_in->AddEntry(hrec_anabin,          "Meas.", "p");	
      legend_in->AddEntry(hunf,                    "Unf.", "p");
      legend_in->AddEntry(hrec_sameside_anabin, "Meas., Same Side", "p");
      legend_in->AddEntry(hunf_sameside,           "Unf., Same Side", "p");
      legend_in->AddEntry(hgen_anabin,          "Gen. Truth", "p");
      legend_in->Draw();
      
      recCutLine1->Draw();
      recCutLine2->Draw();
      genCutLine1->Draw();
      genCutLine2->Draw();
      
      //recMatCutLine1->Draw();
      //recMatCutLine2->Draw();
      //genMatCutLine1->Draw();
      //genMatCutLine2->Draw();
      

      tempCanvForPdfPrint_wLogy->Print(outPdfFile.c_str());
      
      

      //---------------
      
      tempCanvForPdfPrint->cd();
      //h_genratio_oppunf->SetTitle("ratio w/ Gen. Truth");
      h_genratio_oppunf->SetAxisRange(0., 2., "Y");
      h_genratio_oppunf->GetXaxis()->SetMoreLogLabels(true);
      h_genratio_oppunf->GetXaxis()->SetNoExponent(true);
      h_genratio_oppunf->SetAxisRange(boundaries_pt_gen_mat[0], boundaries_pt_gen_mat[nbins_pt_gen_mat],"X");           

      h_genratio_oppunf->Draw();
      h_genratio_oppfold->Draw("SAME");
      h_genratio_oppmeas->Draw("SAME");
      h_genratio_ssunf->Draw("SAME");
      h_genratio_ssfold->Draw("SAME");
      h_genratio_ssmeas->Draw("SAME");

      //h_genratio4->Draw("SAME");
      
      TLegend* legend2 = new TLegend( 0.6,0.8,0.9,0.9 );
      legend2->AddEntry(h_genratio_oppunf, "Unf./Gen. Truth", "p");
      legend2->AddEntry(h_genratio_oppfold, NULL, "p");
      legend2->AddEntry(h_genratio_oppmeas, NULL, "p");
      legend2->AddEntry(h_genratio_ssunf, "Same Side Unf./Gen. Truth", "p");
      legend2->AddEntry(h_genratio_ssfold, NULL, "p");
      legend2->AddEntry(h_genratio_ssmeas, NULL, "p");

      legend2->Draw();

      theLineAtOne->Draw();


      recCutLine1->Draw();
      recCutLine2->Draw();
      genCutLine1->Draw();
      genCutLine2->Draw();

      recMatCutLine1->Draw();
      recMatCutLine2->Draw();
      genMatCutLine1->Draw();
      genMatCutLine2->Draw();

      tempCanvForPdfPrint->Print(outPdfFile.c_str());

      //---------------
      	
      tempCanvForPdfPrint->cd();
      h_recratio_oppunf->SetAxisRange(0., 2., "Y");
      h_recratio_oppunf->GetXaxis()->SetMoreLogLabels(true);
      h_recratio_oppunf->GetXaxis()->SetNoExponent(true);
      h_recratio_oppunf->SetAxisRange(boundaries_pt_gen_mat[0], boundaries_pt_gen_mat[nbins_pt_gen_mat],"X");           

      h_recratio_oppunf->Draw();
      h_recratio_oppfold->Draw("SAME");
      h_recratio_ssunf->Draw("SAME");
      h_recratio_ssfold->Draw("SAME");
      h_recratio_ssgen->Draw("SAME");
      
      TLegend* legend3 = new TLegend( 0.6,0.8,0.9,0.9 );
      legend3->AddEntry(h_recratio_oppunf, NULL, "p");
      legend3->AddEntry(h_recratio_oppfold, NULL, "p");
      legend3->AddEntry(h_recratio_ssunf, NULL, "p");
      legend3->AddEntry(h_recratio_ssfold, NULL, "p");
      legend3->AddEntry(h_recratio_ssgen, NULL, "p");
      legend3->Draw();
      
      theLineAtOne->Draw();

      recCutLine1->Draw();
      recCutLine2->Draw();
      genCutLine1->Draw();
      genCutLine2->Draw();

      recMatCutLine1->Draw();
      recMatCutLine2->Draw();
      genMatCutLine1->Draw();
      genMatCutLine2->Draw();
      
      tempCanvForPdfPrint->Print(outPdfFile.c_str());

      //---------------
	
    
      tempCanvForPdfPrint->Print(close_outPdfFile.c_str());      

    }// end draw pdfs
  }  // end doBayes

  std::cout<<std::endl<<"Bayesian Unfolding done, all output written!"<<std::endl;



  // SVD unfolding ------------------------------------------------------------------------- 

  if(!doSVD) std::cout<<std::endl<<"   skipping SVD unfolding..."<<std::endl<<std::endl;
  else{ std::cout<<std::endl<<std::endl<<std::endl<<"   beginning SVD unfolding..."<<std::endl;
    std::cout<<std::endl<<"initializing kReg parameter array w/ "<<nKregMax<<" elements"<<std::endl;
    std::cout<<kRegCenter-kRegRange<<" <= kReg[i] <= "<<kRegCenter+kRegRange<<std::endl;
    
    int kReg[nKregMax];
    for(int i=(-1*kRegRange); (i+kRegRange)<nKregMax; ++i) kReg[i+kRegRange]=kRegCenter+i;

    if(debugMode)for(int i=(-1*kRegRange); (i+kRegRange)<nKregMax; ++i) 
		   std::cout<<"kReg["<<i+kRegRange<<"]="<<kReg[i+kRegRange]<<std::endl;
    
    if(debugMode)std::cout<<std::endl<<"creating histo-arrays..."<<std::endl;
    TH1 *hrec_folded_ratio[nKregMax];     TH1 *hrec_unfolded_ratio[nKregMax];

    TH2 *hPearsonSVDPriorMeas[nKregMax];      
    TH1 *hFoldedSVDPriorMeas[nKregMax];     
    TH1F *hunf_svd[nKregMax];        
    TH1F *hratio_svd[nKregMax];      
    
    
    std::cout<<"calling RooUnfoldResponse"<<std::endl;
    RooUnfoldResponse roo_resp(hrec_resp_anabin, hgen_resp_anabin, hmat_anabin, ("Response_matrix_anabin"+Rstring).c_str());
    //RooUnfoldResponse roo_resp( 0 , 0, hmat_anabin, ("SVD_hmat_anabin_resp"+RandEtaRange).c_str(), "SVD response matrix") ;
    //RooUnfoldResponse roo_resp( hrec_anabin_clone , hgen_anabin_clone, hmat_anabin, ("SVD_hmat_anabin_resp"+RandEtaRange).c_str(), "SVD response matrix") ;
    
    //roo_resp.UseOverflow(doOverUnderflows);

    if(debugMode)std::cout<<"opening TCanvases..."<<std::endl;
    TCanvas *cPearsonMatrixIter = new TCanvas("cPearsonMatrixIter","", 1500, 1500);      cPearsonMatrixIter->Divide(3,3);	
    TCanvas *cRatio             = new TCanvas("cRatio","",             1500, 1500);      cRatio->Divide(3,3);      		
    TCanvas *cSpectra           = new TCanvas("cSpectra","",           1500, 1500);      cSpectra->Divide(3,3);	   
    TCanvas *cRatioCheck        = new TCanvas("cRatioCheck","",        1500, 1500);    
    TCanvas *c11                = new TCanvas("c11"," Singular Values and divectors", 1500, 1500);      c11->Divide(2);                    

    TLegend *leg[nKregMax],*leg1[nKregMax];      

    std::cout<<"entering SVD Unfolding Loop..."<<std::endl;
    for(int kr = 0; kr < nKregMax; ++kr){
      std::string kRegRandEtaRange_plotTitle=" kReg"+std::to_string(kReg[kr])+RandEtaRange_plotTitle;
      std::string kRegRandEtaRange="_kReg"+std::to_string(kReg[kr])+RandEtaRange;
      
      std::cout<<std::endl<<std::endl<<"kr="<<kr<<" , kReg = "<<kReg[kr]<<std::endl<<std::endl;
   
      // unfold the MC dist, then create a ratio with the gen-level MC dist
      if(debugMode)std::cout<<std::endl<<"calling RooUnfoldSvd and Hreco..."<<std::endl<<std::endl;
      RooUnfoldSvd unf_svd(&roo_resp, hrec_anabin, kReg[kr]);
      hunf_svd[kr] = (TH1F*)unf_svd.Hreco(errorTreatment);
      
      std::cout<<std::endl;   
      hunf_svd[kr]->SetName( ("ppMC_SVDUnf_Spectra"+kRegRandEtaRange).c_str() );
      hunf_svd[kr]->SetLineStyle(33);
      hunf_svd[kr]->SetLineColor(kBlue);
      hunf_svd[kr]->Print("base");
  
//      TH1F* hgen_anabin_forDiv=(TH1F*)hgen_anabin->Clone("hgen_anabin_clone_forSVDdiv");
//      hgen_anabin_forDiv=(TH1F*)hgen_anabin_forDiv->Rebin(nbins_pt_reco, "hgen_anabin_clone_forSVDdiv_recobins" , boundaries_pt_reco);
//
//      hratio_svd[kr]->Divide(hgen_anabin_forDiv);

      hratio_svd[kr] = (TH1F*)hunf_svd[kr]->Clone( ("ppMC_SVDUnf_Ratio"+kRegRandEtaRange).c_str());
      hratio_svd[kr]->Divide(hgen_anabin);
      hratio_svd[kr]->SetMarkerStyle(33);
      hratio_svd[kr]->SetMarkerColor(kRed);
      hratio_svd[kr]->Print("base");

      // get covariance matrix and calculate pearson coefficients
      if(debugMode)std::cout<<std::endl<<"calling Ereco..."<<std::endl;
      TMatrixD covmat = unf_svd.Ereco(errorTreatment);


      if(debugPearson)std::cout<<"debugPearson is ON!"<<std::endl;

      if(debugPearson)std::cout<<std::endl;
      if(debugMode)std::cout<<"calling CalculatePearsonCoefficients..."<<std::endl;
      TMatrixD *pearson = CalculatePearsonCoefficients(&covmat, debugPearson );
      if(debugPearson)std::cout<<std::endl;

      if(debugMode)std::cout<<"creating TH2 for pearson matrix..."<<std::endl;
      hPearsonSVDPriorMeas[kr] = new TH2D (*pearson);
      //hPearsonSVDPriorMeas[kr] = reBinPearsonTH2( pearson );
      
      
      // drawing/visuals
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cPearsonMatrixIter canvas..."<<std::endl;
      cPearsonMatrixIter->cd(kr+1);

      //cPearsonMatrixIter->cd(kr+1)->SetLogz(1);
      //int kr_add1=kr+1;
      //i'm sorry, programming gods
      //TPad * thePad;
      //if       (kr<3)   thePad = new TPad( ("pad"+std::to_string(kr)).c_str(), "thepad" , 0.33*kr    , 0.67 , 0.33*(kr+1)     , 1.00 );
      //else if  (kr<6)   thePad = new TPad( ("pad"+std::to_string(kr)).c_str(), "thepad" , 0.33*(kr%3), 0.34 , 0.33*((kr%3)+1) , 0.67 );
      //else /*if(kr<9)*/ thePad = new TPad( ("pad"+std::to_string(kr)).c_str(), "thepad" , 0.33*(kr%3), 0.00 , 0.33*((kr%3)+1) , 0.34 );
      //cPearsonMatrixIter->SetLogz(1);
      //thePad->SetLogz(1);
      //thePad->Draw();
      //thePad->cd();
      
      //float ZminNmax=1.0;             //int XYmin=1; int XYmax=99;
      //int XYmin=1; int XYmax=1000;
      //hPearsonSVDPriorMeas[kr]->SetMinimum(-1*ZminNmax);  
      //hPearsonSVDPriorMeas[kr]->SetMaximum(1.0);

      hPearsonSVDPriorMeas[kr]->GetZaxis()->SetLabelSize(0.035);
      hPearsonSVDPriorMeas[kr]->GetXaxis()->SetTitle("reco bin num");
      hPearsonSVDPriorMeas[kr]->GetYaxis()->SetTitle("gen bin num");

      //hPearsonSVDPriorMeas[kr]->SetAxisRange(0, nbins_pt, "X");
      //hPearsonSVDPriorMeas[kr]->SetAxisRange(0, nbins_pt, "Y");
      //hPearsonSVDPriorMeas[kr]->SetAxisRange(.0000000001, 1, "Z");

      hPearsonSVDPriorMeas[kr]->SetAxisRange(-1., 1., "Z");
      hPearsonSVDPriorMeas[kr]->SetName( ("hPearsonSVDPriorMeas"+kRegRandEtaRange).c_str() );
      hPearsonSVDPriorMeas[kr]->SetTitle( ("pearsonMatrix SVDPriorMeas"+kRegRandEtaRange_plotTitle).c_str()  );
      hPearsonSVDPriorMeas[kr]->Draw("COLZ");

      if(debugMode)std::cout<<"applying roo_resp to \"truth\" histo hunf_svd[kr="<<kr<<"]..."<<std::endl;
      hFoldedSVDPriorMeas[kr] = roo_resp.ApplyToTruth(hunf_svd[kr]);
      hFoldedSVDPriorMeas[kr]->SetName( ("hFoldedSVDPriorMeas"+kRegRandEtaRange).c_str() );
      //hFoldedSVDPriorMeas[kr]->SetLineStyle(33);
      hFoldedSVDPriorMeas[kr]->SetLineColor(kRed);
      



      if(debugMode)std::cout<<std::endl<<"drawing stuff on cSpectra canvas..."<<std::endl<<std::endl;
      cSpectra->cd(kr+1);
      cSpectra->cd(kr+1)->SetLogy();
      cSpectra->cd(kr+1)->SetLogx();

      // what's happening here exactly?
      //std::cout <<"CHECK: kr="<<kr<<"  and kReg[kr]="<<kReg[kr]<<std::endl<<std::endl;
      TH1F* hrec_anabin_clone=(TH1F*)hrec_anabin->Clone("hrec_anabin_clone");
      hrec_anabin_clone->SetTitle( ("Reco Jet spectra,"+kRegRandEtaRange_plotTitle).c_str() );
      hrec_anabin_clone->SetXTitle("Jet p_{T} (GeV/c)");
      hrec_anabin_clone->SetMarkerStyle(24);
      hrec_anabin_clone->SetMarkerColor(kBlack);
      hrec_anabin_clone->SetAxisRange(boundaries_pt_gen[0], boundaries_pt_gen[nbins_pt_gen], "X");
      hrec_anabin_clone->Print("base");

      hrec_anabin_clone->Draw();
      hunf_svd[kr]->Draw("same");
      hFoldedSVDPriorMeas[kr]->Draw("SAME");


      leg[kr] = new TLegend(0.57, 0.75, 0.9, 0.9, NULL,"NBNDC");//x1,y1,x2,y2,header,option
      leg[kr]->AddEntry(hrec_anabin_clone,"Measured","pl");
      leg[kr]->AddEntry(hunf_svd[kr],"Unfolded","pl");
      leg[kr]->AddEntry(hFoldedSVDPriorMeas[kr],"Folded","pl");
      leg[kr]->SetTextSize(0.04);
      leg[kr]->Draw();

      // draw on cRatio canvas
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio canvas..."<<std::endl;
      cRatio->cd(kr+1);
      cRatio->cd(kr+1)->SetLogx();

      hrec_folded_ratio[kr] = (TH1F*)hFoldedSVDPriorMeas[kr]->Clone( ("ppMC_SVDUnf_reFold_Ratio"+kRegRandEtaRange).c_str() );
      hrec_folded_ratio[kr]->SetTitle( ("Ratios w/ Meas.,"+kRegRandEtaRange_plotTitle).c_str() );
      hrec_folded_ratio[kr]->SetMarkerStyle(27);
      hrec_folded_ratio[kr]->SetMarkerColor(kRed);
      hrec_folded_ratio[kr]->SetXTitle("Jet p_{T} (GeV/c)");
      hrec_folded_ratio[kr]->SetAxisRange(boundaries_pt_gen[0], boundaries_pt_gen[nbins_pt_gen], "X");
      hrec_folded_ratio[kr]->SetAxisRange(0.0, 2.0, "Y");
      hrec_folded_ratio[kr]->Divide(hrec_anabin);
      //hrec_folded_ratio[kr]->SetAxisRange(boundaries_pt_reco[0], boundaries_pt_reco[nbins_pt_reco]);
      //hrec_folded_ratio[kr]->Print("base");
      hrec_folded_ratio[kr]->Draw();

      hrec_unfolded_ratio[kr] = (TH1F*)hunf_svd[kr]->Clone( ("ppMC_SVDUnf_Ratio"+kRegRandEtaRange).c_str());
      hrec_unfolded_ratio[kr]->Print("base");
      hrec_unfolded_ratio[kr]->SetTitle( ("SVD Unf./Meas.,"+kRegRandEtaRange_plotTitle).c_str() );
      hrec_unfolded_ratio[kr]->SetMarkerStyle(27);
      hrec_unfolded_ratio[kr]->SetMarkerColor(kBlue);      
      hrec_unfolded_ratio[kr] = (TH1F*)hrec_unfolded_ratio[kr]->Rebin(nbins_pt_reco, ("ppMC_SVDUnf_Ratio_rebin4div"+kRegRandEtaRange).c_str() , boundaries_pt_reco);
      //if(clearOverUnderflows)TH1clearOverUnderflows((TH1*)hrec_unfolded_ratio[kr]);
      hrec_unfolded_ratio[kr]->Divide(hrec_anabin);
      hrec_unfolded_ratio[kr]->Print("base");
      //hrec_unfolded_ratio[kr]->SetAxisRange(boundaries_pt_reco[0], boundaries_pt_reco[nbins_pt_reco])
      hrec_unfolded_ratio[kr]->Draw("same");

      TLine* theLine= new TLine(boundaries_pt_gen[0],1.,boundaries_pt_gen[nbins_pt_gen],1.);
      theLine->SetLineWidth(1);
      theLine->SetLineStyle(2);
      theLine->SetLineColor(36);
      //theLine->Draw();


      leg1[kr] = new TLegend(0.57, 0.75, 0.9, 0.9, NULL,"BRNDC");
      leg1[kr]->AddEntry(hrec_unfolded_ratio[kr],"Unf./Meas.","pl");
      leg1[kr]->AddEntry(hrec_folded_ratio[kr],"reFold/Meas.","pl");
      //leg->AddEntry(hSVD_prior,"Prior, normalized to data","pl");
      leg1[kr]->SetTextSize(0.04); 
      leg1[kr]->Draw();
  	
      //  singular values and d_i vector ---------------------------
      //Note that these do not depend on the regularization.
      //The opposite: they tell you which regularization to use! (ian note: how?)
      if(kr == kRegDraw){      
	TSVDUnfold *svdUnfold = unf_svd.Impl();
	//svdUnfold->SetNormalize(true);
        //if(debugMode) std::cout<<std::endl<<"    !!! kr==0 !!!"<<std::endl<<std::endl;

        // singular values
        if(debugMode)std::cout << "  getting singular values... " << std::endl;
        TH1 *hSVal = (TH1*)svdUnfold->GetSV();
        if(debugMode)for(int bin=1; bin<=hSVal->GetNbinsX(); bin++)
		       std::cout<<"bin: "<<bin<<",  SV: "<<hSVal->GetBinContent(bin)<< std::endl;

        if(debugMode)std::cout<<std::endl<<"drawing singular values on c11 canvas.."<<std::endl<<std::endl;
        c11->cd(1);
        gPad->SetLogy();        //c11->SetLogy();
        hSVal->SetTitle(" singular values ");
        hSVal->SetXTitle(" singular values ");
        //hSVal->SetAxisRange(0,35,"X");
        hSVal->DrawCopy("E");
	//hSVal->Draw("E");
        drawText( "5.02 TeV pp, ak4PF Jets",0.358173, 0.8459761, 19);
	drawText( "Prompt-Reco, Jet80+LowerJets",0.358173, 0.8159761, 19);
	drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
	drawText( ("kReg="+std::to_string(kReg[kr])).c_str() ,0.408173, 0.7359761, 19);

        // di vector values
        if(debugMode)std::cout << "  getting di vector... " <<  std::endl;
        TH1D *hdi = (TH1D*)svdUnfold->GetD();
        if(debugMode)for(int bin=1; bin<=hdi->GetNbinsX(); bin++)
		       std::cout<<"i: "<<bin<<",  di: "<<hdi->GetBinContent(bin)<<std::endl;

        if(debugMode)std::cout<<std::endl<<"drawing di vectors on c11 canvas.."<<std::endl;
        c11->cd(2);
        gPad->SetLogy();        //cdi->SetLogy();
        hdi->SetTitle(" di vectors ");
        hdi->SetXTitle(" |d_{i}^{kreg}| ");
        //hdi->SetAxisRange(0,35,"X");
	//hdi->Draw("E");
        hdi->DrawCopy("E");
        //drawText( "5.02 TeV ak4PFJets",0.358173, 0.8459761, 19);
	//drawText( "ppData, Prompt Reco HighPtJets",0.358173, 0.8159761, 19);
	//drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
	//drawText( ("kReg="+std::to_string(kReg[kr])).c_str() ,0.408173, 0.7359761, 19);
	const std:: string MCdesc= "QCD Dijets, PY8 Tune CUETP8M1"; 
        if(debugMode)std::cout<<std::endl<<"done with kr==0 specifics"<<std::endl<<std::endl;
      }

      hrec_unfolded_ratio[kr]->Print("base");
      hrec_folded_ratio[kr]->Print("base");
    }// kReg loop ends
    
    // loop over histos pointers in arrays to write
    fout->cd();

    std::cout<<std::endl<<"writing SVD histo arrays to file..."<<std::endl;    fout->cd();
    for(int kr = 0; kr<nKregMax; ++kr){
      hunf_svd[kr]->Write();            
      hratio_svd[kr]->Write();          
      hPearsonSVDPriorMeas[kr]->Write();
      hFoldedSVDPriorMeas[kr]->Write(); 
    }
    
    if(!drawPDFs){ if(debugMode)std::cout<<std::endl<<"NOT drawing PDFs for SVD Unfolding"<<std::endl<<std::endl; }
    else{ std::cout<<std::endl<<"drawing PDFs for SVD Unfolding..."<<std::endl;
      
      // form filename string, open pdf file and draw the canvases we have so far
      std::string outPdfFile=outSVDPdfFile;
      std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";
      cPearsonMatrixIter->Print(open_outPdfFile.c_str());


      //  already-drawn canvases --------------------------------------------------
      cSpectra->cd()           ;   cSpectra->Print(outPdfFile.c_str());
      cRatio->cd()             ;   cRatio->Print(outPdfFile.c_str());
      cPearsonMatrixIter->cd() ;   cPearsonMatrixIter->Print(outPdfFile.c_str());
      c11->cd()                ;   c11->Print(outPdfFile.c_str());


      //  cRatio Check --------------------------------------------------
      cRatioCheck->cd();
      cRatioCheck->SetLogx(1);

      hrec_folded_ratio[kRegDraw]->SetAxisRange(0.5, 1.5, "Y");
      hrec_folded_ratio[kRegDraw]->SetAxisRange(boundaries_pt_reco[0], boundaries_pt_gen[nbins_pt_gen], "X");
      hrec_folded_ratio[kRegDraw]->SetTitle("SVD Unf./Meas.");
      hrec_folded_ratio[kRegDraw]->Draw();

      TLegend * leg2 = new TLegend(0.14, 0.79, 0.34, 0.87, NULL,"NBNDC");
      leg2->AddEntry(hrec_unfolded_ratio[kRegDraw],"Unf./Meas.","pl");
      leg2->AddEntry(hrec_folded_ratio[kRegDraw],"Refold/Meas.","pl");
      //leg->AddEntry(hSVD_prior,"Prior, normalized to data","pl");
      leg2->SetTextSize(0.02);
      leg2->Draw();

      hrec_unfolded_ratio[kRegDraw]->Draw("same");
      drawText( "5.02 TeV pp, ak4PF Jets",          0.14, 0.75, 22);
      drawText( "Prompt-Reco, Jet80+LowerJets",     0.14, 0.72, 22);
      drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.14, 0.69, 22);

      TLine* theLine= new TLine(boundaries_pt_gen[0],1.,boundaries_pt_gen[nbins_pt_gen],1.);
      theLine->SetLineWidth(1);
      theLine->SetLineStyle(2);
      theLine->SetLineColor(36);
      theLine->Draw();

      cRatioCheck->Print(outPdfFile.c_str());
      cRatioCheck->Print(close_outPdfFile.c_str());

      std::cout<<std::endl<<"done drawing SVD PDFs!"<<std::endl<<std::endl;
    }// end drawPDFs
  }// end SVD specific

  if(debugMode) std::cout<<"writing output root file..."<<std::endl;
  fout->Write();

  std::cout<<"closing input/output root files..."<<std::endl;
  fout->Close();
  fpp_MC->Close();
  return 0;

} // end unfoldDataSpectra



//  steering ---------------------------------------------------------------------------------
int main(int argc, char* argv[]){  int rStatus = -1;
  
  if( argc!=5 ){
    std::cout<<"do ./unfoldMCSpectra.exe <targMCDir> <baseOutputName> <doJetID> <kRegCenter>"<<std::endl;
    std::cout<<"actually... just open the damn code and look"<<std::endl;

    return rStatus;  }
  
  rStatus=1; // runtime error
  
  rStatus=unfoldMCSpectra(  (const std::string)argv[1], (const std::string)argv[2], 
			    (int)std::atoi(argv[3]) , (int)std::atoi(argv[4]) 
			    ); 
  
  std::cout<<std::endl<<"done!"<<std::endl<<" return status: "<<rStatus<<std::endl<<std::endl;
  return rStatus;
}














  //std::cout<<"underflow="<<hrec_anabin->GetBinContent(0)<<std::endl;
  //std::cout<<"overflow="<<hrec_anabin->GetBinContent(1001)<<std::endl;
  //float _56to64GeVSum = 0.;
  //for(int i=57 ; i < 65; i++) {
  //  std::cout<<"i="<<i<<std::endl;
  //  std::cout<<"bin center ="<< hrec_anabin->GetBinCenter(i) <<std::endl;
  //  std::cout<<"bin content ="<< hrec_anabin->GetBinContent(i) <<std::endl;
  //  _56to64GeVSum += hrec_anabin->GetBinContent(i);
  //  std::cout<<"Sum now ="<< _56to64GeVSum<<std::endl;
  //
  //}
  //std::cout<<"before rebinning: 56-64GeV's contents are:"<< _56to64GeVSum <<std::endl;//57 because 1 GeV bins-> 56-57 is 57'th bin because 0-1 is the first bin
  //std::cout<<"before rebinning: first bin's width is:"    << hrec_anabin->GetBinWidth(57)<<std::endl;
  //std::cout<<"after rebinning: first bin's center is    :"<< hrec_anabin->GetBinCenter(1)<<std::endl;
  //std::cout<<"after rebinning: first bin's contents are :"<< hrec_anabin->GetBinContent(1)<<std::endl;
  //std::cout<<"after rebinning: first bin's width is     :"<< hrec_anabin->GetBinWidth(1)<<std::endl;
