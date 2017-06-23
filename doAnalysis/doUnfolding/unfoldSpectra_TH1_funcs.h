void divideBinWidth(TH1 *h){
  std::cout<<std::endl<<"in divideBinWidth"<<std::endl<<std::endl;
  bool funcDebug=false;

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


      Float_t val    = h->GetBinContent(i);
      Float_t valErr = h->GetBinError(i);

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
  std::cout<<std::endl<<"divideBinWidth done"<<std::endl<<std::endl;
  return;
}





void TH1zeroBins(TH1* h, int numBins2Clear){
  std::cout<<std::endl<<"in TH1zeroBins"<<std::endl<<std::endl;
  bool funcDebug=true;
  bool clearLastBin=false;
  if(funcDebug&&clearLastBin)std::cout<<"WARNING clearLastBin=true"<<std::endl;
  //h->Sumw2();
  if(funcDebug)h->Print("base");

  //double nEntries=h->GetEntries();
  int numbins=h->GetNbinsX();
  if(funcDebug)std::cout<<"numbins="<<numbins<<std::endl;


  if(numBins2Clear>=numbins || numBins2Clear==0){
    std::cout<<"cannot zero requested # of bins... exiting"<<std::endl<<std::endl;
    return;
  }
  else{
    if(funcDebug)std::cout<<"setting bin content(s) to zero"<<std::endl;

    for(int i=1;i<=numBins2Clear;i++){
      h->SetBinContent(i,0);
      h->SetBinError(i,0);    }

    if(clearLastBin){
      std::cout<<"WARNING clearing last bin"<<std::endl;
      h->SetBinContent(numbins,0);
      h->SetBinError(numbins,0);
    }
  }

  if(funcDebug)h->Print("base");
  std::cout<<std::endl<<"TH1zeroBins done."<<std::endl<<std::endl;
  return;
}




void TH1clearOverUnderflows(TH1* h)
{
  //std::cout<<std::endl<<"WARNING!!!!"<<std::endl<<std::endl;
  //std::cout<<"doOverUnderflows="<<doOverUnderflows<<std::endl;
  std::cout<<"clearing Over/Underflow Bins...."<<std::endl;
  h->TH1::ClearUnderflowAndOverflow();
  return;
}
