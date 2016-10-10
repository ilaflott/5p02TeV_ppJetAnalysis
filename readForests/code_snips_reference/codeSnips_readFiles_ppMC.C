// FROM readFiles_ppMC.C

  // code this up if hardcoded pthatWeights become annoying, for now it should be okay. 
  //// open weight files
  //std::cout<<"MCweightsFile used is "<<MCweightsFile<<std::endl;
  //std::ifstream instr_MCweights(MCweightsFile.c_str(),std::ifstream::in);
  ////read in weights+pthat bins
  
  //// hltanalysis
  ////L1
  //int jet40_l1seed_p_F;  int jet60_l1seed_p_F;  int jet80_l1seed_p_F;  int jet100_l1seed_p_F;  
  ////HLT
  //int jet40_F;           int jet60_F;           int jet80_F;           int jet100_F;           
  //int jet40_p_F;         int jet60_p_F;         int jet80_p_F;         int jet100_p_F;         
  
  //// four trees, four specific HLT trigger pt objects
  ////currently (5/26/16): HLT_AK${radius}CaloJet${trgObjpt}_Eta${}_v
  //std::vector<double> *trgObjpt_40  ;   
  //std::vector<double> *trgObjpt_60  ;   
  //std::vector<double> *trgObjpt_80  ;   
  //std::vector<double> *trgObjpt_100 ;
  
  //// hltanalysis, I should try figuring out a better way to do this part
  //std::string HLTBranches[N_HLTBits]; for(int i=0;i<N_L1Bits;i++) HLTBranches[i]=HLTBitStrings[i]+"_v1";
  //jetpp[3]->SetBranchAddress( HLTBranches[0].c_str() , &jet40_F);
  //jetpp[3]->SetBranchAddress( HLTBranches[1].c_str() , &jet60_F);
  //jetpp[3]->SetBranchAddress( HLTBranches[2].c_str() , &jet80_F);
  //jetpp[3]->SetBranchAddress( HLTBranches[3].c_str() , &jet100_F);
  //std::string HLTPresclBranches[N_HLTBits]; for(int i=0;i<N_HLTBits;i++) HLTPresclBranches[i]=HLTBitStrings[i]+"_v1_Prescl";
  //jetpp[3]->SetBranchAddress( HLTPresclBranches[0].c_str() , &jet40_p_F);
  //jetpp[3]->SetBranchAddress( HLTPresclBranches[1].c_str() , &jet60_p_F);
  //jetpp[3]->SetBranchAddress( HLTPresclBranches[2].c_str() , &jet80_p_F);
  //jetpp[3]->SetBranchAddress( HLTPresclBranches[3].c_str() , &jet100_p_F);
  //std::string L1PresclBranches[N_L1Bits]  ; for(int i=0;i<N_HLTBits;i++) L1PresclBranches[i]=L1BitStrings[i]+"_Prescl";
  //jetpp[3]->SetBranchAddress( L1PresclBranches[0].c_str()  , &jet40_l1seed_p_F);
  //jetpp[3]->SetBranchAddress( L1PresclBranches[1].c_str()  , &jet60_l1seed_p_F);
  //jetpp[3]->SetBranchAddress( L1PresclBranches[2].c_str()  , &jet80_l1seed_p_F);
  //jetpp[3]->SetBranchAddress( L1PresclBranches[3].c_str()  , &jet100_l1seed_p_F);
  
  //// specific HLT path trees
  ////ONE path ONE tree ONE pt branch (see trees[] in header)
  ////e.g. trgObjpt_40 is filled with jet pt from the specific jet40 HLT tree/branch 
  //jetpp[4]->SetBranchAddress("pt",&trgObjpt_40);
  //jetpp[5]->SetBranchAddress("pt",&trgObjpt_60);  
  //jetpp[6]->SetBranchAddress("pt",&trgObjpt_80);  
  //jetpp[7]->SetBranchAddress("pt",&trgObjpt_100);
  


  //// for triggerpt spectra before/after jetID
  //TH1F *hpp_TrgObj40[2];
  //TH1F *hpp_TrgObj60[2];
  //TH1F *hpp_TrgObj80[2];
  //TH1F *hpp_TrgObj100[2];
  //TH1F *hpp_TrgObjComb[2];
  //TH1F *hpp_CombJetpT[2];
  
  ////40, 60, 80, 100
  //hpp_TrgObj40[0]   = new TH1F(Form("hpp_HLT40_noJetID_R%d_%s"                        , radius, etaWidth), 
  //			       Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  //hpp_TrgObj60[0]   = new TH1F(Form("hpp_HLT60_noJetID_R%d_%s"               , radius, etaWidth), 
  //			       Form("Spectra from  Jet 60 && !jet80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  //hpp_TrgObj80[0]   = new TH1F(Form("hpp_HLT80_noJetID_R%d_%s"     , radius, etaWidth), 
  //			       Form("Spectra from  Jet 80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  //hpp_TrgObj100[0]  = new TH1F(Form("hpp_HLT100_noJetID_R%d_%s"     , radius, etaWidth), 
  //			       Form("Spectra from  Jet 100 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  ////combined
  //hpp_TrgObjComb[0] = new TH1F(Form("hpp_HLTComb_noJetID_R%d_%s"    , radius, etaWidth), 
  //			       Form("Trig Combined Spectra R%d %s " , radius, etaWidth), 2000, 0, 2000);
  //hpp_CombJetpT[0]  = new TH1F(Form("hpp_TrgCombTest_noJetID_R%d_%s"           , radius, etaWidth), 
  //			       Form("Trig Combined Spectra KurtMethod R%d %s " , radius, etaWidth), 2000, 0, 2000);
  ////40, 60, 80, 100
  //hpp_TrgObj40[1]   = new TH1F(Form("hpp_HLT40_JetID_R%d_%s"                          , radius, etaWidth), 
  //			       Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  //hpp_TrgObj60[1]   = new TH1F(Form("hpp_HLT60_JetID_R%d_%s"                 , radius, etaWidth), 
  //			       Form("Spectra from  Jet 60 && !jet80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  //hpp_TrgObj80[1]   = new TH1F(Form("hpp_HLT80_JetID_R%d_%s"       , radius, etaWidth), 
  //			       Form("Spectra from  Jet 80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  //hpp_TrgObj100[1]  = new TH1F(Form("hpp_HLT100_JetID_R%d_%s"       , radius, etaWidth), 
  //			       Form("Spectra from  Jet 100 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  ////combined
  //hpp_TrgObjComb[1] = new TH1F(Form("hpp_HLTComb_JetID_R%d_%s"      , radius, etaWidth), 
  //			       Form("Trig Combined Spectra R%d %s " , radius, etaWidth), 2000, 0, 2000);
  //hpp_CombJetpT[1]  = new TH1F(Form("hpp_TrgCombTest_JetID_R%d_%s"             , radius, etaWidth), 
  //			       Form("Trig Combined Spectra KurtMethod R%d %s " , radius, etaWidth), 2000, 0, 2000);
  


    //// total prescale array
    //int treePrescl[4]={ (jet40_p_F*jet40_l1seed_p_F), (jet60_p_F*jet60_l1seed_p_F), 
    //                  (jet80_p_F*jet80_l1seed_p_F), (jet100_p_F*jet100_l1seed_p_F) };    
    
    //// HLT decision array
    //bool trgDec[4]={ (bool)jet40_F, (bool)jet60_F, (bool)jet80_F, (bool)jet100_F  };
    
    //// prefer highPt l1 trigger for jet if it matches multiple l1jet pt criteria, then compute weight
    //double triggerPt=0.;
    //unsigned int trgObj40_size =trgObjpt_40->size();
    //unsigned int trgObj60_size =trgObjpt_60->size();
    //unsigned int trgObj80_size =trgObjpt_80->size();
    //unsigned int trgObj100_size=trgObjpt_100->size();
    
    //if(jet40_F){
    //for(unsigned int itt=0; itt<trgObj40_size; ++itt){
    //	if(trgObjpt_40->at(itt) > triggerPt) triggerPt = trgObjpt_40->at(itt);
    //  }
    //}
    //if(jet60_F){
    //for(unsigned int itt=0; itt<trgObj60_size; ++itt){
    //	if(trgObjpt_60->at(itt) > triggerPt) triggerPt = trgObjpt_60->at(itt);
    //  }
    //}
    //if(jet80_F){
    //for(unsigned int itt=0; itt<trgObj80_size; ++itt){
    //	if(trgObjpt_80->at(itt) > triggerPt) triggerPt = trgObjpt_80->at(itt);
    //  }
    //}
    //if(jet100_F){
    //for(unsigned int itt=0; itt<trgObj100_size; ++itt){
    //	if(trgObjpt_100->at(itt) > triggerPt) triggerPt = trgObjpt_100->at(itt);
    //  }
    //}    
    
    //if(debugMode&&nEvt%250==0)std::cout <<"triggerPt ="<<triggerPt<<std::endl;
    //if(debugMode&&nEvt%250==0)std::cout<<"jet 40, trgDec[0]  = "<< trgDec[0]  <<std::endl;
    //if(debugMode&&nEvt%250==0)std::cout<<"jet 60, trgDec[1]  = "<< trgDec[1]  <<std::endl;
    //if(debugMode&&nEvt%250==0)std::cout<<"jet 80, trgDec[2]  = "<< trgDec[2]  <<std::endl;
    //if(debugMode&&nEvt%250==0)std::cout<<"jet 100, trgDec[3] = "<< trgDec[3]  <<std::endl<<std::endl;


    // input "offline" HLT selection using rawpt here?
    
    //if( trgDec[0] && (triggerPt>=40  &&  triggerPt<60 ) ) is40  = true;    
    //if( trgDec[1] && (triggerPt>=60  &&  triggerPt<80 ) ) is60  = true;    
    //if( trgDec[2] && (triggerPt>=80  &&  triggerPt<100) ) is80  = true;    
    //if( trgDec[3] && (triggerPt>=100                  ) ) is100 = true;    
	
    //if(debugMode&&nEvt%250==0)std::cout<<"after checking trigger pt "<<std::endl;
    //if(debugMode&&nEvt%250==0)std::cout<<"jet 40,  is40  = "<< is40  <<std::endl;
    //if(debugMode&&nEvt%250==0)std::cout<<"jet 60,  is60  = "<< is60  <<std::endl;
    //if(debugMode&&nEvt%250==0)std::cout<<"jet 80,  is80  = "<< is80  <<std::endl;
    //if(debugMode&&nEvt%250==0)std::cout<<"jet 100, is100 = "<< is100  <<std::endl<<std::endl;
	
//if(is40)  hpp_TrgObj40[0]->Fill(recpt, treePrescl[0]);//recpt here+below prescale weighted
      //if(is60)  hpp_TrgObj60[0]->Fill(recpt, treePrescl[1]);
      //if(is80)  hpp_TrgObj80[0]->Fill(recpt, treePrescl[2]);
      //if(is100) hpp_TrgObj100[0]->Fill(recpt, treePrescl[3]);
      //if(true)  hpp_CombJetpT[0]->Fill(recpt, weight_eS);//recpt here is weighted-weighted
      
      //if(is40 )  hpp_TrgObj40[1]->Fill(recpt, treePrescl[0]);
      //if(is60 )  hpp_TrgObj60[1]->Fill(recpt, treePrescl[1]);
      //if(is80 )  hpp_TrgObj80[1]->Fill(recpt, treePrescl[2]);
      //if(is100) hpp_TrgObj100[1]->Fill(recpt, treePrescl[3]);
      //if(true)  hpp_CombJetpT[1]->Fill(recpt, weight_eS);


  //// before JetID
  //hpp_TrgObjComb[0]->Add(hpp_TrgObj100[0]);
  //hpp_TrgObjComb[0]->Add(hpp_TrgObj80[0] );
  //hpp_TrgObjComb[0]->Add(hpp_TrgObj60[0] );
  //hpp_TrgObjComb[0]->Add(hpp_TrgObj40[0] );

  //// after JetID
  //hpp_TrgObjComb[1]->Add(hpp_TrgObj100[1]);
  //hpp_TrgObjComb[1]->Add(hpp_TrgObj80[1] );
  //hpp_TrgObjComb[1]->Add(hpp_TrgObj60[1] );
  //hpp_TrgObjComb[1]->Add(hpp_TrgObj40[1] );
  
  //// before JetID
  //hpp_TrgObj40[0]->Print("base");
  //hpp_TrgObj60[0]->Print("base");
  //hpp_TrgObj80[0]->Print("base");
  //hpp_TrgObj100[0]->Print("base");
  //hpp_TrgObjComb[0]->Print("base");//diff between this
  //hpp_CombJetpT[0]->Print("base");//and this?
  
  //// after JetID
  //hpp_TrgObj40[1]->Print("base");
  //hpp_TrgObj60[1]->Print("base");
  //hpp_TrgObj80[1]->Print("base");
  //hpp_TrgObj100[1]->Print("base");
  //hpp_TrgObjComb[1]->Print("base");
  //hpp_CombJetpT[1]->Print("base");
  
  //if(debugMode)std::cout<<"misc clean up.."<<std::endl;
  //trgObjpt_40->clear();
  //trgObjpt_60->clear();
  //trgObjpt_80->clear();
  //trgObjpt_100->clear();

// FROM readFiles_ppMC.h



//// string arrays
////L1
//const std::string L1BitStrings[]={//this array is a good idea
//  "L1_SingleJet28_BptxAND",
//  "L1_SingleJet40_BptxAND",
//  "L1_SingleJet48_BptxAND",
//  "L1_SingleJet52_BptxAND"
//};
//const int N_L1Bits=sizeof(L1BitStrings)/sizeof(std::string);
//
////HLT
//const std::string HLTBitStrings[N_L1Bits]={
//  "HLT_AK4CaloJet40_Eta5p1",
//  "HLT_AK4CaloJet60_Eta5p1",
//  "HLT_AK4CaloJet80_Eta5p1",
//  "HLT_AK4CaloJet100_Eta5p1"
//};
