void divideBinWidth_TH2(TH2F *h){
  std::cout<<std::endl<<"in divideBinWidth_TH2"<<std::endl<<std::endl;
  bool funcDebug=false;

  //h->Sumw2();
  if(funcDebug) h->Print("base");      std::cout<<std::endl;

  double numEntries=h->GetEntries();
  int nbinsx=h->GetNbinsX();
  int nbinsy=h->GetNbinsY();

  if(funcDebug)std::cout<<"#of xbins = "<<nbinsx<<std::endl;
  if(funcDebug)std::cout<<"#of ybins = "<<nbinsy<<std::endl;
  if(funcDebug)std::cout<<"total #of bins = "<<nbinsy*nbinsx<<std::endl;
  if(funcDebug)std::cout<<"integral before = "<<h->Integral()<<std::endl;
  if(funcDebug)std::cout<<"entries before = "<<numEntries<<std::endl;

  TAxis* the_xAxis=h->GetXaxis();
  TAxis* the_yAxis=h->GetYaxis();
  for (int i=1;i<=nbinsx;i++)    {
    for (int j=1;j<=nbinsy;j++)    {

      Float_t xWidth=h->GetXaxis()->GetBinWidth(i);
      Float_t yWidth=h->GetYaxis()->GetBinWidth(j);



      if(funcDebug)std::cout<<std::endl<<"xbin-i = "<< i <<std::endl;
      if(funcDebug)std::cout<<"xbin-i width = "<<  the_xAxis->GetBinWidth(i)  << std::endl;
      if(funcDebug)std::cout<<"xbin-i pt range is "<<   the_xAxis->GetBinLowEdge(i) << " to "
                            << (the_xAxis->GetBinLowEdge(i)+the_xAxis->GetBinWidth(i)) << std::endl;

      if(funcDebug)std::cout<<std::endl<<"ybin-j = "<< j <<std::endl;
      if(funcDebug)std::cout<<"ybin-j width = "<<  the_yAxis->GetBinWidth(j)  << std::endl;
      if(funcDebug)std::cout<<"ybin-j pt range is "<<   the_yAxis->GetBinLowEdge(j) << " to "
                            << (the_yAxis->GetBinLowEdge(j)+the_yAxis->GetBinWidth(j)) << std::endl;


      if(funcDebug)std::cout<<std::endl<<"old val = "   <<      h->GetBinContent(i,j) << std::endl;
      if(funcDebug)std::cout<<"old valErr = "<<      h->GetBinError(i,j)   << std::endl;


      Float_t val = (h->GetBinContent(i,j))/(xWidth*yWidth);
      //Float_t valErr = (h->GetBinError(i,j))/(xWidth*yWidth);

      h->SetBinContent(i,j,val);
      //h->SetBinError(i,j,valErr);
      h->SetBinError(i,j,0.0);

      if(funcDebug)std::cout<<std::endl<<"new val = "   <<      h->GetBinContent(i,j) << std::endl;
      if(funcDebug)std::cout<<"new valErr = "<<      h->GetBinError(i,j)   << std::endl;

    }
  }

  //h->SetEntries(numEntries);
  if(funcDebug)h->Print("base");      std::cout<<std::endl;

  std::cout<<std::endl<<"divideBinWidth_TH2 done"<<std::endl<<std::endl;
  return;
}




//normalizes MC response TH2 to the rowSum of the gen-axis
void normalizeMC_TH2(TH2F* inputTH2){
  std::cout<<std::endl<<"in normalizeMC_TH2"<<std::endl<<std::endl;

  inputTH2->Sumw2();
  inputTH2->Print("base");  std::cout<<std::endl;

  double numEntries=inputTH2->GetEntries();

  int colStart=1;
  TAxis *xaxis = inputTH2->GetXaxis();
  int nbins_x= xaxis->GetNbins();

  int rowStart=1;
  TAxis *yaxis = inputTH2->GetYaxis();
  int nbins_y = yaxis->GetNbins();

  //if(doOverUnderflows){
  //  nbins_x++;    nbins_y++;
  //  rowStart--;   colStart--; }

  const int rowSums_len=nbins_x;
  float rowSums[rowSums_len]={0.};
  float rowSumErrs[rowSums_len]={0.};

  //calculate sums
  std::cout<<"computing row sums..."<<std::endl;
  for(int rowNum=rowStart;rowNum<=nbins_y;rowNum++){
    float theSum=0.;
    float theSumErr=0.;
    for(int colNum=colStart;colNum<=nbins_x;colNum++){
      float theVal=inputTH2->GetBinContent(colNum,rowNum);
      float theValErr=inputTH2->GetBinError(colNum,rowNum);
      theSum+=theVal;
      theSumErr+=(theValErr*theValErr);
    }
    rowSums[rowNum]=theSum;
    rowSumErrs[rowNum]=std::sqrt(theSumErr);
  }
  std::cout<<std::endl;
  //scale each entry in each row down by that row's sum
  std::cout<<"scaling rows..."<<std::endl;
  for(int  rowNum=rowStart; rowNum<=nbins_y; rowNum++){

    float theSum=rowSums[rowNum];
    float theSumErr=rowSumErrs[rowNum];

    for( int colNum=colStart; colNum<=nbins_x; colNum++){

      float theVal=(inputTH2->GetBinContent(colNum,rowNum));
      float theValErr=(inputTH2->GetBinError(colNum,rowNum));

      float finalVal=theVal/theSum;
      float finalErr=(theVal/theSum)* std::sqrt(
                                                (theSumErr/theSum)*(theSumErr/theSum) +
                                                (theValErr/theVal)*(theValErr/theVal) )  ;


      inputTH2->SetBinContent(colNum,rowNum,finalVal);
      inputTH2->SetBinError(colNum,rowNum,finalErr);
    }
  }
  std::cout<<std::endl;

  inputTH2->SetEntries(numEntries);
  inputTH2->Print("base");

  std::cout<<std::endl<<"exiting normalizeMC_TH2"<<std::endl<<std::endl;
  return ;


}




TH2F* reBinTH2(TH2F* inputTH2, std::string rebinTH2_name,
               double* boundaries_pt_reco, int nbins_pt_reco,
               double* boundaries_pt_gen  , int nbins_pt_gen           ){
  std::cout<<std::endl<<"in reBinTH2"<<std::endl<<std::endl;
  bool funcDebug = true;

  //double numEntries=inputTH2->GetEntries();
  //inputTH2->Sumw2();
  if(funcDebug)inputTH2->Print("base");  std::cout<<std::endl;


  TAxis *xaxis = inputTH2->GetXaxis(); //reco pt axis
  int nbins_x= xaxis->GetNbins();
  int xbinstart=1;

  TAxis *yaxis = inputTH2->GetYaxis(); //gen pt axis
  int nbins_y = yaxis->GetNbins();
  int ybinstart=1;

  TH2F *reBinnedTH2 = new TH2F(rebinTH2_name.c_str(), inputTH2->GetTitle(),
                               nbins_pt_reco, boundaries_pt_reco ,
                               nbins_pt_gen, boundaries_pt_gen );

  //if(doOverUnderflows){
  //  xbinstart--;
  //  ybinstart--;
  //  nbins_x++;
  //  nbins_y++;
  //}
  //reBinnedTH2->Sumw2();

  for (  int xbin=xbinstart ; xbin <= nbins_x ; xbin++ ) {
    for (int ybin=ybinstart ; ybin <= nbins_y ; ybin++ ) {

      reBinnedTH2->Fill( xaxis->GetBinCenter(xbin) , yaxis->GetBinCenter(ybin) ,
                         inputTH2->GetBinContent(xbin,ybin) );

    }  }//end x-y loop

  //reBinnedTH2->SetEntries(numEntries);

  if(funcDebug) std::cout<<"getSumW2N="<<reBinnedTH2->GetSumw2N()<<std::endl;
  if(funcDebug) reBinnedTH2->Print("base");  std::cout<<std::endl;


  std::cout<<std::endl<<"reBinTH2 done"<<std::endl<<std::endl;
  return reBinnedTH2;
}



TH2F* reBinPearsonTH2(TMatrixD* pearson, const double* boundaries_pt, const int nbins_pt){
  std::cout<<std::endl<<"in reBinPearsonTH2"<<std::endl<<std::endl;

  bool debugPearson=false;
  //const double* the_ptBins=boundaries_pt;//ptbins
  //const int numbins=nbins_pt;//nbins

  TH2F *reBinnedTH2 = new TH2F(  "pearsonCustomBins" , "pearsonMatrix custom bins"   ,
                                 nbins_pt,boundaries_pt,
                                 nbins_pt,boundaries_pt);
  std::cout<<std::endl<<"in reBinPearsonTH2 function..."<<std::endl;
  reBinnedTH2->Print("base");

  for (int j=1;   j<=nbins_pt;j++) {
    for (int i=1; i<=nbins_pt;i++) {

      double pearson_i_j=(*pearson)(i-1,j-1);

      if(fabs(pearson_i_j) > 1.0) {
        if(debugPearson)std::cout <<std::endl<<"WARNING:: pearson["<<i-1<<"]["<<j-1<<"]="<<pearson_i_j<<std::endl<<std::endl;
        //      pearson_i_j=10.;
        //std::cout<< "pearson_i_j=10"<<std::endl
        continue; }
      else {
        if(debugPearson)std::cout <<"pearson["<<i-1<<"]["<<j-1<<"]="<<pearson_i_j<<std::endl;   }

      reBinnedTH2->Fill( i, j, fabs( (float) pearson_i_j ) );
    }
  }


  reBinnedTH2->Print("base");

  std::cout<<std::endl<<"reBinPearsonTH2 done"<<std::endl<<std::endl;
  return reBinnedTH2;
}
