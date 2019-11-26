
const int tupelvar_xAx_reBin[N_tupelvars]={
10,//  "Pt",
5,//  "Eta",
5,//  "Phi",
10,//  "E",
1,//  "Id",
10,//  "RawPt",
10,//  "RawE",
10,//  "HfHadE",
10,//  "HfEmE",
5,//  "ChHadFrac",
5,//  "NeutralHadAnHfFrac",
5,//  "ChEmFrac",
5,//  "NeutralEmFrac",
2,//  "ChMult",
2,//  "ConstCnt",
2//  "NeuMult"
};

const float tupeljetQAxmin[N_tupelvars]={
30.,//  "Pt",
-5.16,//  "Eta",
-3.15,//  "Phi",
30.,//  "E",
0.,//  "Id",
30.,//  "RawPt",
30.,//  "RawE",
0.,//  "HfHadE",
0.,//  "HfEmE",
0.,//  "ChHadFrac",
0.,//  "NeutralHadAnHfFrac",
0.,//  "ChEmFrac",
0.,//  "NeutralEmFrac",
0.,//  "ChMult",
0.,//  "ConstCnt",
0.//  "NeuMult"
};

const float tupeljetQAxmax[N_tupelvars]={
100.,//  "Pt",
5.16,//  "Eta",
3.15,//  "Phi",
1000.,//  "E",
5,//  "Id",
1000,//  "RawPt",
1000.,//  "RawE",
1000.,//  "HfHadE",
1000.,//  "HfEmE",
2.,//  "ChHadFrac",
2.,//  "NeutralHadAnHfFrac",
2.,//  "ChEmFrac",
2.,//  "NeutralEmFrac",
100.,//  "ChMult",
100.,//  "ConstCnt",
100.//  "NeuMult"
};

const int var_xAx_reBin[]={                                   
  10, 10,     //reco/rawjtpt  
  //  10, 10, 10,    //no Res, L2Res, L3Res jtpt
  
  5,  2,      //jteta, phi    

  1, 5, 5,    //trks	      
  1, 5,                       

  1, 5, 5,    //ph	      				           
  1, 5,                       

  1, 5, 5,    //ch  	      
  1, 5,                       

  1, 5, 5,    //ne	      	
  1, 5, 5,    //e	      				      
  1, 5, 5,    //mu	      				      
  
  2, 2, 2,    // N Const.     
  
  1, 1, 1,    //xj,Aj,dphi    
  20, 20        //sub/lead jttpt
};

const float jetQAxmax[]={
  (1000. ), (1000.),             //reco/rawjtpt  
  //  (1000. ), (1000.), (1000.), //jtpt forRes, L2Res, L3Res

  (4.5   ), (3.2  ),  		 //jteta, phi    

  (59) /*(60 )*/, (1.95), (1.95), 	 //trks	      //note; 60 is commented out because of the way axis range setting works
  (59) /*(60 )*/, (1.95), 		              // see: https://root.cern.ch/doc/master/TAxis_8cxx_source.html#l00405 
                                                      // the way to "catch NaN" causes the overflow bin in these hists to be drawn
  (59 ), (1.95), (1.95), 	 //ph	      
  (59 ), (1.95), 		                 

  (59 ), (1.95), (1.95),	 //ch  	      
  (59 ), (1.95), 		              
   
  (19 ), (1.95), (1.95),	 //ne	      
  (9 ), (1.95), (1.95),	 //e	      
  (9 ), (1.95), (1.95),	 //mu	      

  (98),   (98),   (98),   	 // N Const.     

  (1.   ), (1.  ),   (3.2  ), 	 //xj,Aj,dphi    
  (1000. ), (1000.) 		 //sub/lead jttpt
};


const float jetQAxmin[]={
  (50. ), (50.  ) ,           //reco/rawjtpt  
  //  (50. ), (50.  ) ,         (50.  ) ,        //jtpt forRes, L2Res, L3 Res

  (-5. ), (-3.2 ) , 	      //jteta, phi    

  (0.  ), (0.   ) , (0.),      //trks	      
  (0.  ), (0.   ) , 

  (0.  ), (0.   ) , (0.),      //ph	      
  (0.  ), (0.   ) , 

  (0.  ), (0.   ) , (0.),      //ch  	      
  (0.  ), (0.   ) , 

  (0.  ), (0.   ) , (0.),      //ne	      
  (0.  ), (0.   ) , (0.),      //e	      
  (0.  ), (0.   ) , (0.),      //mu	      

  (0.  ), (0.   ) , (0.),      // N Const.     

  (0.  ), (0.   ) , (2.0),      //xj,Aj,dphi    
  (50. ), (50.  )          //sub/lead jttpt
};		



//const int jetTrigQABinning=10;
const int jetSpectraRapBinning=20;

// hist painting ------------------------

// line colors
//const int theDataOverlayLineColor=1, theMCOverlayLineColor=1;//, altOverlayLineColor=3; 
//const int theRatioLineColor=1;//,altRatioLineColor1=8, altRatioLineColor2=7;

// marker colors
//const int theDataOverlayMarkerColor=2, theMCOverlayMarkerColor=4;//,theRatioMarkerColor=9;
//const std::string theMCColor="kBlue";
//const Color_t theMCColor="kBlue";

// hist painting ------------------------

const int   theTrigOverlayLineColor[] ={  1,  1,  1,  1,  1, 1        };//,  1 };
const int theTrigOverlayMarkerColor[] ={  1,  2,  3,  6,  7, 8        };//,  4 };
const int      theTrigOverlayMarker[] ={ 20, 20, 20, 20, 20, 20       };//, 32 };
const int  theRapOverlayMarkerColor[] ={  2,  3,  6,  7,  1,  4 , 8, 9};






void globalHistStyle(){


  //// global style settings

  gStyle->SetOptStat(0);  
  //gStyle->SetOptStat("nemruoi");
  
  //gStyle->SetLegendBorderSize(0);
  //gStyle->SetLegendFillColor(0);
  
  gStyle->SetPadTickY(1);   
  gStyle->SetTickLength(0.015,"y");
  //gStyle->SetLabelSize(0.02,"y");
  //gStyle->SetTitleOffset(0.89,"y");
  //gStyle->SetTitleSize(0.03,"y");
  
  gStyle->SetPadTickX(1);
  //gStyle->SetTickLength(0.015,"x");
  //gStyle->SetLabelSize(0.02,"x");
  //gStyle->SetTitleOffset(0.94,"x");
  //gStyle->SetTitleSize(0.03,"x");

  gROOT->ForceStyle();
  
  return;

}

void globalTrigHistStyle(){
  
  
  //// global style settings
  
  gStyle->SetOptStat(0);  
  //gStyle->SetOptStat("nemruoi");
  
  //gStyle->SetLegendBorderSize(0);
  //gStyle->SetLegendFillColor(0);
  gStyle->SetTitleFontSize(0.03);
  
  gStyle->SetPadTickY(1);   
  gStyle->SetTickLength(0.015,"y");
  
  //gStyle->SetLabelSize(0.04,"y");
  //gStyle->SetTitleOffset(0.89,"y");
  //gStyle->SetTitleSize(0.03,"y");
  
  
  gStyle->SetPadTickX(1);
  //gStyle->SetTickLength(0.015,"x");
  
  //gStyle->SetLabelSize(0.08,"x");
  //gStyle->SetTitleOffset(0.94,"x");
  //gStyle->SetTitleSize(0.04,"x");
  
  gROOT->ForceStyle();
  
  return;

}



void setupJetTrigSpectraRatioCanvas(TCanvas* canv, 
				    TPad* specpad,  TPad* ratpad){  
  canv->cd();  
  specpad->SetLogx(1);  specpad->SetLogy(1);
  //specpad->SetLogx(0);  specpad->SetLogy(1);
  //specpad->SetGridx(1);  specpad->SetGridy(1);  
  specpad->SetBottomMargin(0);
  specpad->Draw();  
  ratpad->SetLogx(1);  ratpad->SetLogy(0);
  //ratpad->SetLogx(0);  ratpad->SetLogy(0);
  //ratpad->SetGridx(1);  ratpad->SetGridy(0);  
  ratpad->SetTopMargin(0);  ratpad->SetBottomMargin(0.3);  
  ratpad->Draw();  
  return;
}


//1=kBlack
//2=kRed
//3=kGreen
//4=kBlue
const int theDataLineColor=1,theDataMarkerColor=1;
const float theDataMarkerSize=0.95;
void dataHistStyle(TH1* h ){
  
  //marker
  h->SetMarkerStyle(kFullDotLarge);
  h->SetMarkerSize(theDataMarkerSize);
  h->SetMarkerColor(theDataMarkerColor);
  
  //line
  h->SetLineColor(theDataLineColor);
  return;
}

const int theMCLineColor=2, theMCMarkerColor=1;
const float theMCMarkerSize=0.;
void MCHistStyle(TH1* h ){
  
  //marker
  h->SetMarkerColorAlpha(theMCMarkerColor,0.);
  //h->SetMarkerStyle(kFullDotLarge);
  //h->SetMarkerSize(theMCMarkerSize);
  
  //line
  h->SetLineColor(theMCLineColor);
  return;
}


void ratioHistStyle(TH1* theRatio , std::string h_XAx_Title, std::string h_YAx_Title=AUAxTitle){
  
  //theRatio->SetAxisRange(0.,2.,"Y");    
  //theRatio->SetAxisRange(0.5,1.5,"Y");    
  //theRatio->SetAxisRange(0.25,1.75,"Y");    
  //theRatio->SetAxisRange(0.4,1.1,"Y");    
  theRatio->GetYaxis()->SetLabelFont(43); 
  theRatio->GetYaxis()->SetLabelSize(13);
  
  theRatio->GetYaxis()->SetTitleSize(15);
  theRatio->GetYaxis()->SetTitleFont(43);
  theRatio->GetYaxis()->SetTitleOffset(2);
  
  theRatio->SetYTitle(h_YAx_Title.c_str() );
  
  
  theRatio->GetXaxis()->SetLabelFont(43); 
  theRatio->GetXaxis()->SetLabelSize(15);
  
  theRatio->GetXaxis()->SetTitleSize(20);
  theRatio->GetXaxis()->SetTitleFont(43);
  theRatio->GetXaxis()->SetTitleOffset(4.);
  
  theRatio->SetXTitle( h_XAx_Title.c_str() );
  
  theRatio->SetTitle("");
  theRatio->SetLineColor( theMCLineColor );    
  
  return;
}


void jetIDratioHistStyle(TH1* theRatio , std::string h_XAx_Title, std::string h_YAx_Title=AUAxTitle){
  
  //theRatio->SetAxisRange(0.,2.,"Y");    
  //theRatio->SetAxisRange(0.5,1.5,"Y");    
  //theRatio->SetAxisRange(0.25,1.75,"Y");    
  //theRatio->SetAxisRange(0.4,1.1,"Y");  
  //theRatio->GetYaxis()->SetNoExponent(true);  
  //  theRatio->GetYaxis()->SetMoreLogLabels(true);  
  theRatio->GetYaxis()->SetLabelFont(43); 
  theRatio->GetYaxis()->SetLabelSize(14);
  
  theRatio->GetYaxis()->SetTitleSize(16);
  theRatio->GetYaxis()->SetTitleFont(43);
  theRatio->GetYaxis()->SetTitleOffset(1.1);
  
  theRatio->SetYTitle(h_YAx_Title.c_str() );
  
  
  theRatio->GetXaxis()->SetLabelFont(43); 
  theRatio->GetXaxis()->SetLabelSize(16);
  
  theRatio->GetXaxis()->SetTitleSize(20);
  theRatio->GetXaxis()->SetTitleFont(43);
  theRatio->GetXaxis()->SetTitleOffset(4.);
  
  theRatio->SetXTitle( h_XAx_Title.c_str() );
  
  theRatio->SetTitle("");
  theRatio->SetLineColor( theMCLineColor );    
  
  return;
}

void trigRatioHistStyle(TH1* h, int j){
  
  h->SetMarkerStyle( theTrigOverlayMarker[j]);
  h->SetMarkerColor( theTrigOverlayMarkerColor[j] );
  h->SetMarkerSize( 0.);//0.90);

  h->SetLineColor(theTrigOverlayMarkerColor[j]);
  
  h->GetXaxis()->SetLabelSize(0.08);
  h->GetXaxis()->SetTitleSize(0.08);
  h->GetXaxis()->SetTitleOffset(1.1);
  
  h->GetYaxis()->SetLabelSize(0.08);
  h->GetYaxis()->SetTitleSize(0.08);
  h->GetYaxis()->SetTitleOffset(0.4);
  //h->GetYaxis()->SetTitle("Combo/Indiv.");
  
  
  //h->SetLineColor(   theTrigOverlayLineColor[j]   );
  
  //h->SetAxisRange(40.,1000.,"X");
  
  return;
}

void trigRatioHistStyle2(TH1* h, int j){
  
  h->SetMarkerStyle( 24 );
  h->SetMarkerColor( theTrigOverlayMarkerColor[j] );
  h->SetMarkerSize( 0.90);
  
  h->GetXaxis()->SetLabelSize(0.08);
  h->GetXaxis()->SetTitleSize(0.08);
  
  h->GetYaxis()->SetLabelSize(0.08);
  //h->GetYaxis()->SetTitleSize(0.08);
  //h->GetYaxis()->SetTitle("Combo/Indiv.");
  
  
  
  h->SetLineColor(   theTrigOverlayLineColor[j]   );
  
  //h->SetAxisRange(40.,1000.,"X");
  
  return;
}

void trigSpectraHistStyle(TH1* h, int j){
  
  h->SetMarkerStyle(theTrigOverlayMarker[j]);
  h->SetMarkerColor( theTrigOverlayMarkerColor[j] );
  //h->SetMarkerSize(0.90);
  
  h->SetLineColor(   theTrigOverlayLineColor[j]   );
  
  //h->SetAxisRange(40.,1000.,"X");
  //h->SetAxisRange(0.,50.,"Y");
  
  return;
}

void trigSpectra2HistStyle(TH1* h, int j){
  
  //h->SetMarkerStyle(theTrigOverlayMarker[j]);
  h->SetMarkerStyle(24);
  h->SetMarkerColor( theTrigOverlayMarkerColor[j] );
  h->SetMarkerSize(0.90);
  
  h->SetLineColor(   theTrigOverlayLineColor[j]   );
  
  //h->SetAxisRange(40.,1000.,"X");
  //h->SetAxisRange(0.,50.,"Y");
  
  return;
}


void trigPtHistStyle(TH1* h, int j){
  
  h->SetMarkerStyle(theTrigOverlayMarker[j]);
  h->SetMarkerColor( theTrigOverlayMarkerColor[j] );
  //h->SetMarkerSize(0.90);
  
  h->SetLineColor(   theTrigOverlayLineColor[j]   );
  
  //h->GetXaxis()->SetLabelSize(0.04);
  //h->GetXaxis()->SetTitleSize(0.04);
  //
  //h->GetYaxis()->SetLabelSize(0.04);
  //h->GetYaxis()->SetTitleSize(0.04);
  
  //h->SetAxisRange(40.,1000.,"X");
  
  return;
}




void trigEtaHistStyle(TH1* h, int j){
  
  h->SetMarkerStyle(theTrigOverlayMarker[j]);
  h->SetMarkerColor( theTrigOverlayMarkerColor[j] );
  h->SetMarkerSize(0.90);
  
  h->SetLineColor(   theTrigOverlayLineColor[j]   );
  
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitleSize(0.04);
  
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTitleSize(0.04);
  
  h->SetAxisRange(-5.1, 5.1,"X");
  
  return;
}



TLatex* makeTLatex(float x, float y, std::string inputString){
  TLatex *t = new TLatex( x,y,(inputString).c_str());
  t->SetTextFont(63);
  t->SetTextColor(kBlack);
  t->SetTextSize(18);
  t->SetLineWidth(1);
  t->SetNDC();
  return t;
}

TLatex* makeTLatex(float x, float y, std::string inputString, float size){
  TLatex *t = new TLatex( x,y,(inputString).c_str());
  t->SetTextFont(63);
  t->SetTextColor(kBlack);
  t->SetTextSize(size);
  t->SetLineWidth(1);
  t->SetNDC();
  return t;
}

