////////////////////////////
//deriveDijetResonse snips//
/////////////////////////////



  // for(int i=0; i<nbins_pt; i++){
  //   for(int j=0; j<nbins_eta; j++){
  //     avgA[i][j] = avgA[i][j]/(double)nEntries[i][j];
  //     avgB[i][j] = avgB[i][j]/(double)nEntriesB[i][j];
  //     if(nEntries[i][j]) hRelResponse[i]->SetBinContent(j+1,(1+avgA[i][j])/(1-avgA[i][j]));
  //     //if(nEntries[i][j]) hRelResponse[i]->SetBinError(j+1,hRelResponse[i]->GetBinContent(j+1)*(1./sqrt(nEntries[i][j])));
  //     if(nEntries[i][j]) hRelResponse[i]->SetBinError(j+1,avgAHisto[i][j]->GetRMS()*(1./sqrt(nEntries[i][j])));
  //     else hRelResponse[i]->SetBinContent(j+1,0);
      
  //     if(nEntriesB[i][j]) hMPFResponse[i]->SetBinContent(j+1,(1+avgB[i][j])/(1-avgB[i][j]));
  //     //if(nEntriesB[i][j]) hMPFResponse[i]->SetBinError(j+1,hMPFResponse[i]->GetBinContent(j+1)*(1./sqrt(nEntriesB[i][j])));
  //     if(nEntriesB[i][j]) hMPFResponse[i]->SetBinError(j+1,avgBHisto[i][j]->GetRMS()*(1./sqrt(nEntriesB[i][j])));
  //     else hMPFResponse[i]->SetBinContent(j+1,0);
   
  //     if(nEntriesAbs[i][j]) hMPFAbsPhoResponse[i]->SetBinContent(j+1,avgAbsPhoResponse[i][j]/nEntriesAbs[i][j]);
  //     //if(nEntriesAbs[i][j]) hMPFAbsPhoResponse[i]->SetBinError(j+1,hMPFAbsPhoResponse[i]->GetBinContent(j+1)*(1./sqrt(nEntriesAbs[i][j])));
  //     if(nEntriesAbs[i][j]) hMPFAbsPhoResponse[i]->SetBinError(j+1,hAvgAbsPhoResponse[i][j]->GetRMS()*(1./sqrt(nEntriesAbs[i][j])));
  //     else hMPFAbsPhoResponse[i]->SetBinContent(j+1,0);

  //     int totEntriesA=0, totEntriesB=0, totEntriesAbs=0;
  //     for(int j=0; j<nbins_eta; j++){ 
  // 	totEntriesA+=nEntries[i][j]; 
  // 	totEntriesB+=nEntriesB[i][j];
  // 	totEntriesAbs+=nEntriesAbs[i][j];
  //     }
  //     hRelResponse[i]->SetEntries(totEntriesA);
  //     hMPFResponse[i]->SetEntries(totEntriesB);
  //     hMPFAbsPhoResponse[i]->SetEntries(totEntriesAbs);
  //   }
  // }



    // hRelResponse[i]->SetMarkerColor(color[i]);
    // hRelResponse[i]->SetLineColor(color[i]);
    // hRelResponse[i]->SetTitle(Form("Rel, %g<p^{avg}_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    // hRelResponse[i]->Write();
    // //hRelResponse[i]->Scale(1./hRelResponse[i]->Integral());

    // hMPFResponse[i]->SetMarkerColor(color[i]);
    // hMPFResponse[i]->SetLineColor(color[i]);
    // hMPFResponse[i]->SetMarkerStyle(21);
    // hMPFResponse[i]->SetLineStyle(2);
    // hMPFResponse[i]->SetTitle(Form("MPF dijet, %g<p^{avg}_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    // hMPFResponse[i]->Write();
    // //hMPFResponse[i]->Scale(1./hMPFResponse[i]->Integral());

    // hMPFAbsPhoResponse[i]->SetMarkerColor(color[i]);
    // hMPFAbsPhoResponse[i]->SetLineColor(color[i]);
    // hMPFAbsPhoResponse[i]->SetTitle(Form("MPF #gamma+Jet, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    // hMPFAbsPhoResponse[i]->Write();



























////////////////////////////
//RAA_read_data_pp snips//
/////////////////////////////
    // still don't know why this was here in original code
    // scared to get rid of it completely
    //int maxtrg= -1;      
    //for(int ii=4; ii>=0; ii--){
    //  if(trgDec[ii]==1){
    //	maxtrg=ii;
    //	break;
    //  }
    //}

// //INITIALIZATIONS
// //switches
// const bool debugMode     = true; //int iii=0;    /*debug*/    //std::cout<<"here"<<iii<<std::endl;iii++;
// const bool fastDebugMode = false; //int iii=0;    /*debug*/    //std::cout<<"here"<<iii<<std::endl;iii++;
// 
// //for convenience during testing, coding, etc.
// const bool atMIT  = true;
// const std::string hadoopDir    = "/mnt/hadoop/cms";
// const std::string xrootdDirect = "root://cmsxrootd.fnal.gov/";
// 
// //filelist to run over
// //const std::string infile_Forest = "HighPtJet80_5p02TeV_pp_data_forests.txt";
// const std::string infile_Forest = "testFile_2015_5p02TeV_pp_data_forests.txt";//2 files long, for debugging
// 
// //default output file name
// const std::string defaultOutputName = "test_output.root"; 
// const std::string outdir="";
// 
// //pt binning
// const int ptbins[] = {15, 30, 50, 80, 120, 170, 220, 300, 500};
// const int nbins_pt = sizeof(ptbins)/sizeof(int)-1;
// 
// //eta binning
// static const char *etaWidth = (char*)"20_eta_20";
// 
// //root file directories + tree names
// const std::string defaultDirNames[]={ "hltanalysis"    ,
// 				      "skimanalysis"   ,
// 				      "GARBAGE ENTRY"  , //later replaced with Form("ak%d%sJetAnalyzer",radius, jetType.c_str()),
// 				      //"ppTrack"        ,
// 				      "hiEvtAnalyzer"  ,
// 				      "hltobject"      , //hlt40
// 				      "hltobject"      , //hlt60
// 				      "hltobject"      , //hlt80
// 				      "hltobject"      }; //hlt100
// 
// const std::string trees[]={ "HltTree"                    ,
// 			    "HltTree"                    ,
// 			    "t"                          , 
// 			    //"trackTree"                  ,
// 			    "HiTree"                     ,
// 			    "HLT_AK4CaloJet40_Eta5p1_v"  ,
// 			    "HLT_AK4CaloJet60_Eta5p1_v"  ,
// 			    "HLT_AK4CaloJet80_Eta5p1_v"  ,
// 			    "HLT_AK4CaloJet100_Eta5p1_v" };
// const int N_trees = sizeof(defaultDirNames)/sizeof(std::string);//note, for string array, don't take off one
//                                                                 //not a C++ reason, just a way the loops are organized in this routine
// 
// //Add jet variables
// const std::string var[] = {"jtpt" , "rawpt", "jteta"    , "jtphi", "trkMax", "trkSum"   , "trkHardSum", 
// 			   "chMax", "chSum", "chHardSum", "phMax", "phSum" , "phHardSum", "neMax"     , 
// 			   "neSum", "eMax" , "eSum"     , "muMax", "muSum" , "Aj"       , "xj"        };
// const int N_vars = sizeof(var)/sizeof(std::string);


// static const int nbins_pt = 39;
// static const double boundaries_pt[nbins_pt+1] = {
//   3, 4, 5, 7, 9, 12, 
//   15, 18, 21, 24, 28,
//   32, 37, 43, 49, 56,
//   64, 74, 84, 97, 114,
//   133, 153, 174, 196,
//   220, 245, 272, 300, 
//   330, 362, 395, 430,
//   468, 507, 548, 592,
//   638, 686, 1000 
// };
// these are the only radii we are interested for the RAA analysis: 2,3,4,5
// static const int no_radius = 3; 
// static const int list_radius[no_radius] = {2,3,4};
// static const int  trigValue = 4;
// static const char trigName [trigValue][256] = {"HLT40","HLT60","HLT80","Combined"};






    // vector <float> arry_TrgObj;
    // if(trgObjpt_40->size()!=0)
    //   for(unsigned i = 0; i<trgObjpt_40->size(); ++i) arry_TrgObj.push_back(trgObjpt_40->at(i));
    // if(trgObjpt_60->size()!=0) 
    //   for(unsigned i = 0; i<trgObjpt_60->size(); ++i) arry_TrgObj.push_back(trgObjpt_60->at(i));
    // if(trgObjpt_80->size()!=0) 
    //   for(unsigned i = 0; i<trgObjpt_80->size(); ++i) arry_TrgObj.push_back(trgObjpt_80->at(i));
    // if(trgObjpt_100->size()!=0)
    //   for(unsigned i = 0; i<trgObjpt_100->size(); ++i) arry_TrgObj.push_back(trgObjpt_100->at(i));

    // // find the max in the array
    // float TrgObjMax = 0.0;
    // for(unsigned i = 0; i<arry_TrgObj.size();++i){
    //   if(arry_TrgObj[i] > TrgObjMax){
    // 	TrgObjMax = arry_TrgObj[i];
    //   }
    // }

    // if(jet40_F && TrgObjMax>=40 && TrgObjMax<60) is40 = true;
    // if(jet60_F && TrgObjMax>=60 && TrgObjMax<80) is60 = true;
    // if(jet80_F && TrgObjMax>=80 && TrgObjMax<100) is80 = true;
    // if(jet100_F && TrgObjMax>=100) is100 = true;
    // int weight_jet40  = jet40_p_F * jet40_l1seed_p_F; 
    // int weight_jet60  = jet60_p_F * jet60_l1seed_p_F; 
    // int weight_jet80  = jet80_p_F * jet80_l1seed_p_F; 
    // int weight_jet100 = jet100_p_F * jet100_l1seed_p_F;


    // double weight_eS = 1.0;
    // if(is40) weight_eS = weight_jet40;
    // if(is60) weight_eS = weight_jet60;
    // if(is80) weight_eS = weight_jet80;
    // if(is100) weight_eS = weight_jet100;

    // arry_TrgObj.clear();


  // TH1F *hpp_JEC_TrgObj80;
  // TH1F *hpp_JEC_TrgObj60;
  // TH1F *hpp_JEC_TrgObj40;
  // TH1F *hpp_JEC_TrgObjComb;

  // TH1F *hpp_JEC_minus_TrgObj80;
  // TH1F *hpp_JEC_minus_TrgObj60;
  // TH1F *hpp_JEC_minus_TrgObj40;
  // TH1F *hpp_JEC_minus_TrgObjComb;

  // TH1F *hpp_JEC_gaus_TrgObj80;
  // TH1F *hpp_JEC_gaus_TrgObj60;
  // TH1F *hpp_JEC_gaus_TrgObj40;
  // TH1F *hpp_JEC_gaus_TrgObjComb;
  
  // TH1F *hpp_Smear_TrgObj80;
  // TH1F *hpp_Smear_TrgObj60;
  // TH1F *hpp_Smear_TrgObj40;
  // TH1F *hpp_Smear_TrgObjComb;
  
  // TH1F *hpp_anaBin_TrgObj80;
  // TH1F *hpp_anaBin_TrgObj60;
  // TH1F *hpp_anaBin_TrgObj40;
  // TH1F *hpp_anaBin_TrgObjComb;




  // hpp_JEC_TrgObj80 = new TH1F(Form("hpp_JEC_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_TrgObj60 = new TH1F(Form("hpp_JEC_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_TrgObj40 = new TH1F(Form("hpp_JEC_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_TrgObjComb = new TH1F(Form("hpp_JEC_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),1000, 0, 1000);

  // hpp_JEC_minus_TrgObj80 = new TH1F(Form("hpp_JEC_minus_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_minus_TrgObj60 = new TH1F(Form("hpp_JEC_minus_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_minus_TrgObj40 = new TH1F(Form("hpp_JEC_minus_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_minus_TrgObjComb = new TH1F(Form("hpp_JEC_minus_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),1000, 0, 1000);

  // hpp_JEC_gaus_TrgObj80 = new TH1F(Form("hpp_JEC_gaus_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_gaus_TrgObj60 = new TH1F(Form("hpp_JEC_gaus_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_gaus_TrgObj40 = new TH1F(Form("hpp_JEC_gaus_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_JEC_gaus_TrgObjComb = new TH1F(Form("hpp_JEC_gaus_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),1000, 0, 1000);

  
  // hpp_Smear_TrgObj80 = new TH1F(Form("hpp_Smear_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_Smear_TrgObj60 = new TH1F(Form("hpp_Smear_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_Smear_TrgObj40 = new TH1F(Form("hpp_Smear_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_Smear_TrgObjComb = new TH1F(Form("hpp_Smear_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),1000, 0, 1000);

  // hpp_anaBin_TrgObj80 = new TH1F(Form("hpp_anaBin_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  // hpp_anaBin_TrgObj60 = new TH1F(Form("hpp_anaBin_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  // hpp_anaBin_TrgObj40 = new TH1F(Form("hpp_anaBin_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  // hpp_anaBin_TrgObjComb = new TH1F(Form("hpp_anaBin_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);





      // float rawpt = rawpt_F[jet];

      // float JEC_upShift   = (float)recpt;
      // float JEC_downShift = (float)recpt;
      // float JEC_gaussmear = (float)recpt;

      // if(jet40_F == 1 && jet60_F==0 && jet80_F==0) {
      // 	hpp_TrgObj40->Fill(recpt, jet40_p_F);
      // // 	//hpp_JEC_TrgObj40->Fill(recpt * (1 + rnd.Gaus(0, 0.01)), jet40_p_F);
      // // 	hpp_JEC_TrgObj40->Fill(JEC_upShift, jet40_p_F);
      // // 	hpp_JEC_minus_TrgObj40->Fill(JEC_downShift, jet40_p_F);
      // // 	hpp_JEC_gaus_TrgObj40->Fill(JEC_gaussmear, jet40_p_F);
      // // 	hpp_Smear_TrgObj40->Fill(recpt+ rnd.Gaus(0,0.01), jet40_p_F);
      // // 	hpp_anaBin_TrgObj40->Fill(recpt, jet40_p_F);
      // }
    
      // if(jet60_F == 1 && jet80_F==0) {
      // 	hpp_TrgObj60->Fill(recpt);
      // // 	//hpp_JEC_TrgObj60->Fill(recpt * (1 + rnd.Gaus(0, 0.01)));
      // // 	hpp_JEC_TrgObj60->Fill(JEC_upShift);
      // // 	hpp_JEC_minus_TrgObj60->Fill(JEC_downShift);
      // // 	hpp_JEC_gaus_TrgObj60->Fill(JEC_gaussmear);
      // // 	hpp_Smear_TrgObj60->Fill(recpt+rnd.Gaus(0,0.01));
      // // 	hpp_anaBin_TrgObj60->Fill(recpt);
      // }








      // if(jet80_F == 1) {      
      // 	hpp_TrgObj80->Fill(recpt);
      // // 	hpp_JEC_TrgObj80->Fill(JEC_upShift);
      // // 	hpp_JEC_minus_TrgObj80->Fill(JEC_downShift);
      // // 	hpp_JEC_gaus_TrgObj80->Fill(JEC_gaussmear);
      // // 	hpp_Smear_TrgObj80->Fill(recpt+rnd.Gaus(0,0.01));
      // // 	hpp_anaBin_TrgObj80->Fill(recpt);
      // }      









  // divideBinWidth(hpp_TrgObjComb);
  // divideBinWidth(hpp_TrgObj80);
  // divideBinWidth(hpp_TrgObj60);
  // divideBinWidth(hpp_TrgObj40);
  
  // hpp_JEC_TrgObjComb->Add(hpp_JEC_TrgObj80);
  // hpp_JEC_TrgObjComb->Add(hpp_JEC_TrgObj60);
  // hpp_JEC_TrgObjComb->Add(hpp_JEC_TrgObj40);

  // divideBinWidth(hpp_JEC_TrgObjComb);
  // divideBinWidth(hpp_JEC_TrgObj80);
  // divideBinWidth(hpp_JEC_TrgObj60);
  // divideBinWidth(hpp_JEC_TrgObj40);

  // hpp_JEC_minus_TrgObjComb->Add(hpp_JEC_minus_TrgObj80);
  // hpp_JEC_minus_TrgObjComb->Add(hpp_JEC_minus_TrgObj60);
  // hpp_JEC_minus_TrgObjComb->Add(hpp_JEC_minus_TrgObj40);

  // divideBinWidth(hpp_JEC_minus_TrgObjComb);
  // divideBinWidth(hpp_JEC_minus_TrgObj80);
  // divideBinWidth(hpp_JEC_minus_TrgObj60);
  // divideBinWidth(hpp_JEC_minus_TrgObj40);

  // hpp_JEC_gaus_TrgObjComb->Add(hpp_JEC_gaus_TrgObj80);
  // hpp_JEC_gaus_TrgObjComb->Add(hpp_JEC_gaus_TrgObj60);
  // hpp_JEC_gaus_TrgObjComb->Add(hpp_JEC_gaus_TrgObj40);

  // divideBinWidth(hpp_JEC_gaus_TrgObjComb);
  // divideBinWidth(hpp_JEC_gaus_TrgObj80);
  // divideBinWidth(hpp_JEC_gaus_TrgObj60);
  // divideBinWidth(hpp_JEC_gaus_TrgObj40);
  
  // hpp_Smear_TrgObjComb->Add(hpp_Smear_TrgObj80);
  // hpp_Smear_TrgObjComb->Add(hpp_Smear_TrgObj60);
  // hpp_Smear_TrgObjComb->Add(hpp_Smear_TrgObj40);

  // divideBinWidth(hpp_Smear_TrgObjComb);
  // divideBinWidth(hpp_Smear_TrgObj80);
  // divideBinWidth(hpp_Smear_TrgObj60);
  // divideBinWidth(hpp_Smear_TrgObj40);

  // hpp_anaBin_TrgObjComb->Add(hpp_anaBin_TrgObj80);
  // hpp_anaBin_TrgObjComb->Add(hpp_anaBin_TrgObj60);
  // hpp_anaBin_TrgObjComb->Add(hpp_anaBin_TrgObj40);

  // divideBinWidth(hpp_anaBin_TrgObjComb);
  // divideBinWidth(hpp_anaBin_TrgObj80);
  // divideBinWidth(hpp_anaBin_TrgObj60);
  // divideBinWidth(hpp_anaBin_TrgObj40);
  
