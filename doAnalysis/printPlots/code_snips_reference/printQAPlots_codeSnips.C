	//{//pthat weighted only MC, can use this for Vz weighting
	//  TH1F* theMCEvtQAHist2= (TH1F*)finMC->Get( "hpthatWeightedVz" );
	//  theMCEvtQAHist2->Scale( 1./theMCEvtQAHist2->GetBinWidth(0) );
	//  theMCEvtQAHist2->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist2->Integral() );
	//  theMCEvtQAHist2->SetMarkerStyle(kDot);
	//  theMCEvtQAHist2->SetMarkerSize(1.3);
	//  theMCEvtQAHist2->SetMarkerColor( kBlack );
	//  theMCEvtQAHist2->SetLineColor( altOverlayLineColor );
	//  theMCEvtQAHist2->Draw("SAME");            
	//  theEvtQALeg->AddEntry(theMCEvtQAHist2,"MC w/o vz-weights","lp");	
	//}
	

	////for ratio w/ pre-vz-weights MC, use this for Vz weighting
	//TH1F* theDataEvtQAHist2=(TH1F*)theDataEvtQAHist->Clone("DataVzClone4Ratio2");
	//
	//TH1F* theMCEvtQAHist2= (TH1F*)finMC->Get( "hpthatWeightedVz" );
	//theMCEvtQAHist2->Scale( 1./theMCEvtQAHist2->GetBinWidth(0) );
	//theMCEvtQAHist2->Scale( theDataEvtQAHist2->Integral()/theMCEvtQAHist2->Integral() );
	//
	//TH1F *theRatio2=theDataEvtQAHist2;
	//theRatio2->SetLineColor( altRatioLineColor1 );	
	//theRatio2->Divide(theMCEvtQAHist2);	
	//theEvtQALeg->AddEntry(theRatio2,"MC not vz-weighted","lp"); 		
	//
	//theRatio2->Draw("SAME");	
	
	////for ratio w/ vzweights
	//TH1F* theDataEvtQAHist3=(TH1F*)theDataEvtQAHist->Clone("DataVzClone4Ratio3");
	//
	//TH1F* theMCEvtQAHist3= (TH1F*)finMC->Get( "hWeightedAnnaVz" );
	//theMCEvtQAHist3->Scale( 1./theMCEvtQAHist3->GetBinWidth(0) );
	//theMCEvtQAHist3->Scale( theDataEvtQAHist3->Integral()/theMCEvtQAHist3->Integral() );
   	//
	//TH1F *theRatio3=theDataEvtQAHist3;
	//theRatio3->SetLineColor( altRatioLineColor2 );	
	//theRatio3->Divide(theMCEvtQAHist3);	
	//theEvtQALeg->AddEntry(theRatio3,"Annas VzWeights","lp"); 	
	//theRatio3->Draw("SAME");
	


      /*

      */
      






















/*
  // ratio plots----------------------
  {   //vz, MC v Data
  std::string Datafilename=fullDataFilename; 
    TFile *finData = new TFile(Datafilename.c_str());

    std::string theHistName_den="hVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName_den<<std::endl;
    TH1F* theJetQAHist_den= (TH1F*)finData->Get( theHistName_den.c_str() );
    theJetQAHist_den->Scale(1/integratedLuminosity);
    theJetQAHist_den->Scale(1/theJetQAHist_den->GetBinWidth(1));

    std::string theHistName_den2="hWeightedVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName_den<<std::endl;
    TH1F* theJetQAHist_den2= (TH1F*)finData->Get( theHistName_den.c_str() );
    theJetQAHist_den2->Scale(1/integratedLuminosity);
    theJetQAHist_den2->Scale(1/theJetQAHist_den2->GetBinWidth(1));
    
    {//MC, norm, no weights
      
      std::string theHistName_num="hVz";
      if(debugMode)std::cout<<"theHistName="<<theHistName_num<<std::endl;
      TH1F* theJetQAHist_num = (TH1F*)finMC->Get( theHistName_num.c_str() );      
      theJetQAHist_num->Scale(1/theJetQAHist_num->GetBinWidth(1));
      theJetQAHist_num->Scale(theJetQAHist_den->Integral()/theJetQAHist_num->Integral());
      
      TH1F* theRatio=theJetQAHist_num;
      theRatio->Divide(theJetQAHist_den);
      theRatio->SetTitle("evtQA, ratio, MC norm., data");
      theRatio->SetYTitle("MC/Data");
      theRatio->SetXTitle("vz (cm)");
      theRatio->SetAxisRange(0.,2.5,"Y");
      theRatio->SetMarkerColor(kRed+2);
      theRatio->Draw();


      TH1F* theJetQAHist_num2= (TH1F*)theJetQAHist_num->Clone("cloneHist");
      theJetQAHist_num2->Scale(theJetQAHist_den2->Integral()/theJetQAHist_num2->Integral());

      TH1F* theRatio2=theJetQAHist_num2;
      theRatio2->Divide(theJetQAHist_den2);
      theRatio2->SetMarkerColor(kAzure);
      theRatio2->Draw("same");

      TLegend* theLeg = new TLegend(0.1,0.1,0.3,0.3,NULL,"NBNDC");
      theLeg->AddEntry(theRatio,"MC/data w/o trigWghts");
      theLeg->AddEntry(theRatio2,"MC/data w/ trig weights");
      theLeg->Draw();

      temp_canv->Print( thePDFFileName.c_str() );
    }

    {//MC, norm, vz weights
      std::string theHistName_num="hvzWeightedVz";
      if(debugMode)std::cout<<"theHistName="<<theHistName_num<<std::endl;
      TH1F* theJetQAHist_num = (TH1F*)finMC->Get( theHistName_num.c_str() );
      theJetQAHist_num->Scale(1/theJetQAHist_num->GetBinWidth(1));
      theJetQAHist_num->Scale(theJetQAHist_den->Integral()/theJetQAHist_num->Integral());
      
      TH1F* theRatio=theJetQAHist_num;
      theRatio->Divide(theJetQAHist_den);
      theRatio->SetTitle("evtQA, ratio, MC norm. w/ vz weights, data");
      theRatio->SetYTitle("MC/Data");
      theRatio->SetXTitle("vz (cm)");
      theRatio->SetAxisRange(0.,2.5,"Y");
      theRatio->SetMarkerColor(kRed+2);
      theRatio->Draw();


      TH1F* theJetQAHist_num2= (TH1F*)theJetQAHist_num->Clone("cloneHist");
      theJetQAHist_num2->Scale(theJetQAHist_den2->Integral()/theJetQAHist_num2->Integral());

      TH1F* theRatio2=theJetQAHist_num2;
      theRatio2->Divide(theJetQAHist_den2);
      theRatio2->SetMarkerColor(kAzure);
      theRatio2->Draw("same");

      TLegend* theLeg = new TLegend(0.1,0.1,0.3,0.3,NULL,"NBNDC");
      theLeg->AddEntry(theRatio,"MC/data w/o trigWghts");
      theLeg->AddEntry(theRatio2,"MC/data w/ trig weights");
      theLeg->Draw();

      temp_canv->Print( thePDFFileName.c_str() );
    }

    {//MC, norm, pthat weights
      std::string theHistName_num="hpthatWeightedVz";
      if(debugMode)std::cout<<"theHistName="<<theHistName_num<<std::endl;
      TH1F* theJetQAHist_num = (TH1F*)finMC->Get( theHistName_num.c_str() );
      theJetQAHist_num->Scale(1/theJetQAHist_num->GetBinWidth(1));
      theJetQAHist_num->Scale(theJetQAHist_den->Integral()/theJetQAHist_num->Integral());
      
      TH1F* theRatio=theJetQAHist_num;
      theRatio->Divide(theJetQAHist_den);
      theRatio->SetTitle("evtQA, ratio, MC norm. w/ pthat weights, data");
      theRatio->SetYTitle("MC/Data");
      theRatio->SetXTitle("vz (cm)");
      theRatio->SetAxisRange(0.,2.5,"Y");
      theRatio->SetMarkerColor(kRed+2);
      theRatio->Draw();

      TH1F* theJetQAHist_num2= (TH1F*)theJetQAHist_num->Clone("cloneHist");
      theJetQAHist_num2->Scale(theJetQAHist_den2->Integral()/theJetQAHist_num2->Integral());

      TH1F* theRatio2=theJetQAHist_num2;
      theRatio2->Divide(theJetQAHist_den2);
      theRatio2->SetMarkerColor(kAzure);
      theRatio2->Draw("same");

      TLegend* theLeg = new TLegend(0.1,0.1,0.25,0.25,NULL,"NBNDC");
      theLeg->AddEntry(theRatio,"MC/data w/o trigWghts");
      theLeg->AddEntry(theRatio2,"MC/data w/ trig weights");
      theLeg->Draw();

      temp_canv->Print( thePDFFileName.c_str() );
    }

    {//MC, norm, all weights
      std::string theHistName_num="hWeightedVz";
      if(debugMode)std::cout<<"theHistName="<<theHistName_num<<std::endl;
      TH1F* theJetQAHist_num = (TH1F*)finMC->Get( theHistName_num.c_str() );
      theJetQAHist_num->Scale(1/theJetQAHist_num->GetBinWidth(1));
      theJetQAHist_num->Scale(theJetQAHist_den->Integral()/theJetQAHist_num->Integral());
      
      TH1F* theRatio=theJetQAHist_num;
      theRatio->Divide(theJetQAHist_den);
      theRatio->SetTitle("evtQA, ratio, MC norm. w/ all weights, data");
      theRatio->SetYTitle("MC/Data");
      theRatio->SetXTitle("vz (cm)");
      theRatio->SetAxisRange(0.,2.5,"Y");
      theRatio->SetMarkerColor(kRed+2);
      theRatio->Draw();

      TH1F* theJetQAHist_num2= (TH1F*)theJetQAHist_num->Clone("cloneHist");
      theJetQAHist_num2->Scale(theJetQAHist_den2->Integral()/theJetQAHist_num2->Integral());

      TH1F* theRatio2=theJetQAHist_num2;
      theRatio2->Divide(theJetQAHist_den2);
      theRatio2->SetMarkerColor(kAzure);
      theRatio2->Draw("same");

      TLegend* theLeg = new TLegend(0.1,0.1,0.25,0.25,NULL,"NBNDC");
      theLeg->AddEntry(theRatio,"MC/data w/o trigWghts");
      theLeg->AddEntry(theRatio2,"MC/data w/ trig weights");
      theLeg->Draw();

      temp_canv->Print( thePDFFileName.c_str() );
    }

    }*/
  






  //// unfolding plots----------------------
  //
  ////unfolding inputs
  //{  
  //  std::string theHistName="hpp_gen_"+pfRad_etaWidth;
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}
  //{  
  //  std::string theHistName="hpp_reco_"+pfRad_etaWidth;
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}
  //{  
  //  std::string theHistName="hpp_matrix_"+pfRad_etaWidth;
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH2F* theJetQAHist= (TH2F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}
  //
  ////MCClosure inputs
  //{  
  //  std::string theHistName="hpp_mcclosure_gen_"+pfRad_etaWidth;
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}
  //{  
  //  std::string theHistName="hpp_mcclosure_data_"+pfRad_etaWidth;
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}
  //{  
  //  std::string theHistName="hpp_mcclosure_data_train_"+pfRad_etaWidth;
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}
  //{  
  //  std::string theHistName="hpp_mcclosure_matrix_"+pfRad_etaWidth;
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH2F* theJetQAHist= (TH2F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}

  //// Jet Corrections----------------------
  ////JEC
  //{  
  //  std::string theHistName="";
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}
  ////JER
  //{  
  //  std::string theHistName="";
  //  if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
  //  TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
  //  theJetQAHist->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}





  ////data only ratio, vz
  //{
  //  std::string theHistName_num="hWeightedVz";
  //  if(debugMode)std::cout<<"theHistName="<<theHistName_num<<std::endl;
  //  TH1F* theJetQAHist_num= (TH1F*)fin->Get( theHistName_num.c_str() );
  //
  //  std::string theHistName_den="hVz";
  //  if(debugMode)std::cout<<"theHistName="<<theHistName_den<<std::endl;
  //  TH1F* theJetQAHist_den= (TH1F*)fin->Get( theHistName_den.c_str() );
  //
  //  TH1F* theRatio=theJetQAHist_num;
  //  theRatio->Divide(theJetQAHist_den);
  //  theRatio->SetTitle("data Event QA, vz w/ over vz w/o TrigWeight");
  //  theRatio->SetYTitle("weighted/unweighted");
  //  theRatio->SetXTitle("vz (cm)");
  //  theRatio->SetAxisRange(0., 18., "Y");
  //
  //  temp_canv->Draw();
  //  temp_canv->Print( thePDFFileName.c_str() );
  //}
  ////----------------------





    /*
    //else if(drawMCPlots)
    // MC Only evt hists
    {//unweighted Vz
    std::string theHistName="hVz";
    
    TH1F* theEvtQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
    
    theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));
    theJetQAHist->SetTitle("MC evt QA, vz w/o weights");
    theJetQAHist->SetYTitle("Num. MC Evts/bin");
    theJetQAHist->SetXTitle("vz (cm)");
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
    }
    {//vzweighted vz
    std::string theHistName="hvzWeightedVz";
    
    TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
    theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));
    theJetQAHist->SetTitle("MC event QA, vz w/ vz weights only");
    theJetQAHist->SetYTitle("millibarns");
    theJetQAHist->SetXTitle("vz (cm)");
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
    }
    {//pthat
    std::string theHistName="hpthat";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
    theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));
    
    theJetQAHist->SetTitle("MC event QA, pthat, not weighted");
    theJetQAHist->SetYTitle("Num. Evts/bin");
    theJetQAHist->SetXTitle("pthat (GeV)");
    
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
    }
    {//pthat weighted pthat
    std::string theHistName="hWeightedpthat";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)finMC->Get( theHistName.c_str() );
    theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));
    
    theJetQAHist->SetTitle("MC event QA, pthat, weighted");
    theJetQAHist->SetYTitle("Cross-section (millibarn)");
    theJetQAHist->SetXTitle("pthat (GeV)");
    
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
    }

    //Data Only Hists
    //else if(drawDataPlots)
    {//Data unweighted Vz
    std::string theHistName="hVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    
    TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
    theJetQAHist->SetTitle("data Event QA, vz, without TrigWeight");
    theJetQAHist->SetYTitle("millibarn");
    theJetQAHist->SetXTitle("vz (cm)");
    
    theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));
    theJetQAHist->Scale(1/integratedLuminosity);
    
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
    }
    //----------------------
    */  

