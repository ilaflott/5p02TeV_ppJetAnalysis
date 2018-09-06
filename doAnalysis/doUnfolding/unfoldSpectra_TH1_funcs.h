void divideBinWidth(TH1 *h){
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in divideBinWidth"<<std::endl<<std::endl;
  //h->Sumw2();
  if(funcDebug)h->Print("base");
  double nEntries=h->GetEntries();
  int numbins=h->GetNbinsX();

  //std::cout<<"#of bin should be = "<<nbins_pt<<std::endl;
  if(funcDebug)std::cout<<"#of bins = "<<numbins<<std::endl;
  if(funcDebug)std::cout<<"integral before = "<<h->Integral()<<std::endl;
  if(funcDebug)std::cout<<"entries before = "<<nEntries<<std::endl;

  for (int i=1;i<=numbins;i++)
    {

      if(funcDebug)std::cout<<std::endl<<"bin-i = "<< i <<std::endl;
      if(funcDebug)std::cout<<"bin-i width = "<<  h->GetBinWidth(i)  << std::endl;
      if(funcDebug)std::cout<<"bin-i pt range is "<<   h->GetBinLowEdge(i) << " to " << (h->GetBinLowEdge(i)+h->GetBinWidth(i)) << std::endl;
      if(funcDebug)std::cout<<std::endl<<"old val = "   <<      h->GetBinContent(i) << std::endl;
      if(funcDebug)std::cout<<"old valErr = "<<      h->GetBinError(i)   << std::endl;


      Double_t val    = h->GetBinContent(i);
      Double_t valErr = h->GetBinError(i);

      val    /= h->GetBinWidth(i);
      valErr /= h->GetBinWidth(i);

      h->SetBinContent(i,val);
      h->SetBinError(i,valErr);

      if(funcDebug)std::cout<<"new val = "   <<      h->GetBinContent(i) << std::endl;
      if(funcDebug)std::cout<<"new valErr = "<<      h->GetBinError(i)   << std::endl;

    }

  if(funcDebug)std::cout<<"integral after = "<<h->Integral()<<std::endl;
  if(funcDebug)std::cout<<"entries after = "<<h->GetEntries()<<std::endl;

  //h->SetEntries(nEntries);

  //std::cout<<"exiting function"<<std::endl;
  if(funcDebug)std::cout<<std::endl<<"divideBinWidth done"<<std::endl<<std::endl;
  return;
}


void multiplyBinWidth(TH1 *h){
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in divideBinWidth"<<std::endl<<std::endl;
  //h->Sumw2();
  if(funcDebug)h->Print("base");
  double nEntries=h->GetEntries();
  int numbins=h->GetNbinsX();

  //std::cout<<"#of bin should be = "<<nbins_pt<<std::endl;
  if(funcDebug)std::cout<<"#of bins = "<<numbins<<std::endl;
  if(funcDebug)std::cout<<"integral before = "<<h->Integral()<<std::endl;
  if(funcDebug)std::cout<<"entries before = "<<nEntries<<std::endl;

  for (int i=1;i<=numbins;i++)
    {

      if(funcDebug)std::cout<<std::endl<<"bin-i = "<< i <<std::endl;
      if(funcDebug)std::cout<<"bin-i width = "<<  h->GetBinWidth(i)  << std::endl;
      if(funcDebug)std::cout<<"bin-i pt range is "<<   h->GetBinLowEdge(i) << " to " << (h->GetBinLowEdge(i)+h->GetBinWidth(i)) << std::endl;
      if(funcDebug)std::cout<<std::endl<<"old val = "   <<      h->GetBinContent(i) << std::endl;
      if(funcDebug)std::cout<<"old valErr = "<<      h->GetBinError(i)   << std::endl;


      Double_t val    = h->GetBinContent(i);
      Double_t valErr = h->GetBinError(i);

      val    *= h->GetBinWidth(i);
      valErr *= h->GetBinWidth(i);

      h->SetBinContent(i,val);
      h->SetBinError(i,valErr);
      
      if(funcDebug)std::cout<<"new val = "   <<      h->GetBinContent(i) << std::endl;
      if(funcDebug)std::cout<<"new valErr = "<<      h->GetBinError(i)   << std::endl;
      
    }

  if(funcDebug)std::cout<<"integral after = "<<h->Integral()<<std::endl;
  if(funcDebug)std::cout<<"entries after = "<<h->GetEntries()<<std::endl;

  //h->SetEntries(nEntries);

  //std::cout<<"exiting function"<<std::endl;
  if(funcDebug)std::cout<<std::endl<<"divideBinWidth done"<<std::endl<<std::endl;
  return;
}





void TH1zeroBins(TH1* h, int numBins2Clear_lowpt, int numBins2Clear_highpt){
  bool funcDebug=false;
  bool clearLastBin=false;
  if(funcDebug)std::cout<<std::endl<<"in TH1zeroBins"<<std::endl<<std::endl;
  if(funcDebug&&clearLastBin)std::cout<<"WARNING clearLastBin=true"<<std::endl;
  //h->Sumw2();
  if(funcDebug)h->Print("base");
  
  //double nEntries=h->GetEntries();
  int numbins=h->GetNbinsX();
  if(funcDebug)std::cout<<"numbins="<<numbins<<std::endl;

  int numBins2Clear= numBins2Clear_lowpt + numBins2Clear_highpt;
  if(numBins2Clear>=numbins || numBins2Clear==0){
    std::cout<<"cannot zero requested # of bins... exiting"<<std::endl<<std::endl;
    return;
  }
  else{
    if(funcDebug)std::cout<<"setting bin content(s) to zero"<<std::endl;
    
    for(int i=1;i<=numBins2Clear_lowpt;i++){
      h->SetBinContent(i,0);
      h->SetBinError(i,0);    }
    
    for(int i=numbins;(numbins-i)<numBins2Clear_highpt;i--){
      h->SetBinContent(i,0);
      h->SetBinError(i,0);    }
  }
  
  
  if(funcDebug)h->Print("base");
  if(funcDebug)std::cout<<std::endl<<"TH1zeroBins done."<<std::endl<<std::endl;
  return;
}




void TH1clearOverUnderflows(TH1* h)
{
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"WARNING!!!!"<<std::endl<<std::endl;
  //if(funcDebug)std::cout<<"doOverUnderflows="<<doOverUnderflows<<std::endl;
  if(funcDebug)std::cout<<"clearing Over/Underflow Bins...."<<std::endl;
  h->TH1::ClearUnderflowAndOverflow();
  return;
}


//returns the highest bin # in histogram h that has non-zero bin content
void TH1FindHighBin(TH1* h){
  int nbins=h->GetNbinsX();
  std::cout<<"nbins="<<nbins<<std::endl;
  for(int i=0;i<=(nbins+1);i++ ){
    
    
    
    double content=h->GetBinContent(i);
    
    if(content!=0)
      std::cout<<"i="<<i<<" content="<<content<<std::endl;
    else if( i!=0 && i!=(nbins+1) )
      std::cout<<"bin i="<<i<<" binlowedge="<<h->GetBinLowEdge(i)<<" has no content"<<std::endl;
    else 
      std::cout<<"bin i="<<i<< " has no content"<<std::endl;
    
    if(i==0 && content!=0)
      std::cout<<"WARNING! Underflow bin has content="<<content<<std::endl;
    if(i==(nbins+1) && content!=0)
      std::cout<<"WARNING! Overflow bin has contents="<<content<<std::endl;
    
  }
  return;
}
