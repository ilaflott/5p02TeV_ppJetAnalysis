	else if(var[j]=="jteta"||
		var[j]=="jtphi"||
		var[j]=="jty"){ 
		
		
	  temp_canvJet->cd();								//draw command for eta/phi/y
	  jetpad1->SetLogy(0);
	  jetpad1->Draw();
	  jetpad1->cd();
	  theDataJetQAHist->SetAxisRange(0,0.45,"Y"); //This line breaks the plot if it's not in logy(0)
	  theDataJetQAHist->Draw(); 
	  theMCJetQAHist->Draw("SAME"); 
	  
	  float t1Loc1=0.66, t1Loc2=0.74; 	    
	  TLatex *t1=new TLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());
	  t1->SetTextFont(63);
	  t1->SetTextColor(kBlack);
	  t1->SetTextSize(10);
	  t1->SetLineWidth(1);
	  t1->SetNDC();
	  t1->Draw();
	  
	  TLatex *t2=new TLatex((t1Loc1),(t1Loc2-.05),"p_{t}>50 GeV, |#eta|<2.0");
	  t2->SetTextFont(63);
	  t2->SetTextColor(kBlack);
	  t2->SetTextSize(10);
	  t2->SetLineWidth(1);
	  t2->SetNDC();
	  t2->Draw();	

	  float legx1=0.78, legx2=legx1+0.11;
	  float legy1=0.68, legy2=legy1+0.09;
	  TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
	  theJetQALeg->AddEntry(theDataJetQAHist,"data","lp");
	  theJetQALeg->AddEntry(theMCJetQAHist,  "MC"  ,"lp");
	  theJetQALeg->Draw();

	  //Ratio Plots
	temp_canvJet->cd();	  
	jetpad2->Draw();
	jetpad2->cd();
	TH1F* theRatio=(TH1F*)theDataJetQAHist->Clone("dataJetHistClone4Ratio");
   
	theRatio->SetMarkerColor( kRed );
	theRatio->SetLineColor( theRatioLineColor );
	
	theRatio->SetAxisRange(0,2.,"Y");
	theRatio->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	

	theRatio->GetYaxis()->SetTitleSize(15);
    theRatio->GetYaxis()->SetTitleFont(43);
    theRatio->GetYaxis()->SetTitleOffset(2);
    theRatio->GetYaxis()->SetLabelFont(43); 
    theRatio->GetYaxis()->SetLabelSize(13);
    theRatio->GetXaxis()->SetTitleSize(20);
    theRatio->GetXaxis()->SetTitleFont(43);
    theRatio->GetXaxis()->SetTitleOffset(4.);
    theRatio->GetXaxis()->SetLabelFont(43); 
    theRatio->GetXaxis()->SetLabelSize(15);
	
	theRatio->SetTitle("");
	theRatio->SetXTitle( h_XAx_Title.c_str() );
	theRatio->SetYTitle("Data/MC");
	theRatio->Divide(theMCJetQAHist);
	theRatio->Draw("h");
	
	TLine* lineAtOne          = new TLine(jetQAxmin[j],1.0,jetQAxmax[j],1.0); 
	lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
	TLine* lineAtOneHalf      = new TLine(jetQAxmin[j],0.5,jetQAxmax[j],0.5); 
	lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
	TLine* lineAtOneEtOneHalf = new TLine(jetQAxmin[j],1.5,jetQAxmax[j],1.5); 
	lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
	
	lineAtOne->Draw("same");
	lineAtOneHalf->Draw("same");
	lineAtOneEtOneHalf->Draw("same");
	  