int main (int argc, char *argv[]){
  
  std::cout<<" now opening Data File "<<std::endl<<ppData_fullFilename<<std::endl<<std::endl;
  TFile* finData = new TFile(ppData_fullFilename.c_str());
  
  std::cout<<" now opening MC File "<<std::endl<<ppMC_fullFilename<<std::endl<<std::endl;
  TFile* finMC = new TFile(ppMC_fullFilename.c_str());
  
  
  TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inHistName.c_str() );
  theDataEvtQAHist->Scale( 1/theDataEvtQAHist->GetBinWidth(0) );
  theDataEvtQAHist->Scale( 1/theLumi );
  
  theDataEvtQAHist->SetTitle (    h_Title.c_str() );
  theDataEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
  theDataEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
  
  theDataEvtQAHist->SetMarkerStyle(kDot);
  theDataEvtQAHist->SetMarkerSize(1.1);
  theDataEvtQAHist->SetMarkerColor( kBlack);
  theDataEvtQAHist->SetLineColor( theRatioLineColor );
  theDataEvtQAHist->SetAxisRange(0.,1.5,"Y");
  

  TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( "hpthatWeightedVz" );
  theMCEvtQAHist->Scale( 1/theMCEvtQAHist->GetBinWidth(0) );
  theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );

  TH1F *theRatio=theDataEvtQAHist;
  theRatio->SetLineColor( altRatioLineColor1 );
  theRatio->Divide(theMCEvtQAHist);
  theEvtQALeg->AddEntry(theRatio,"MC not vz-weighted","lp");

  //theRatio->Draw("SAME");
  

  Int_t NvzWeightBins=60;
  Float_t theVzBinWidth=theRatio->TH1::GetBinWidth(0);
  Float_t xLow=-15., xHigh=15.;
  
  std::cout<<"now grabbing vzWeights for "<< NvzWeightBins<<"bins for "<<xLow<<"< vz <"<<xHigh << std::endl;
  for (int i=1;i<NvzWeightBins+1;++i){//binsX loop
    Float_t leftSideOfBin=-15.+(i-1)*theVzBinWidth;
    Float_t rightSideOfBin=-15.+(i)*theVzBinWidth;
    Float_t vzWeight = theRatio->TH1::GetBinContent(i);
    
    if(i!=NvzWeightBins)
      std::cout<<"  or i="<<i<<", "<<leftSideOfBin<<"<vz<="<<rightSideOfBin<<", vzWeight="<<vzWeight<< std::endl;
  else
    std::cout<<"  or i="<<i<<", "<<leftSideOfBin<<"<vz<="<<rightSideOfBin<<", vzWeight="<<vzWeight<< std::endl;
  }


}
