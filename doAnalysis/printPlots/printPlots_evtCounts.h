void printDataEventCountReport(TFile* fin){
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"printDataEventCountReport Called"<<std::endl<<std::endl;
  
  if(!fin){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }

  TH1F *h_NEvents_skipped = (TH1F*)fin->Get("NEvents_skipped");  
  TH1F *h_NEvents         = (TH1F*)fin->Get("NEvents");
  TH1F *h_NEvents_read    = (TH1F*)fin->Get("NEvents_read");        
  TH1F *h_NEvents_skimCut = (TH1F*)fin->Get("NEvents_skimCut");
  TH1F *h_NEvents_vzCut   = (TH1F*)fin->Get("NEvents_vzCut");

  TH1F *h_NEvents_trigd_1=(TH1F*)fin->Get("NEvents_trigd_1");
  TH1F *h_NEvents_trigd_2=(TH1F*)fin->Get("NEvents_trigd_2");


  //  TH1F *h_NEvents_jet100=(TH1F*)fin->Get("NEvents_jet100Trigd"),    *h_NEvents_jet80=(TH1F*)fin->Get("NEvents_jet80Trigd"),    *h_NEvents_jet60=(TH1F*)fin->Get("NEvents_jet60Trigd"),    *h_NEvents_jet40=(TH1F*)fin->Get("NEvents_jet40Trigd");
  TH1F *h_NEvents_HLTMB  = (TH1F*)fin->Get("NEvents_HLTMB");
  TH1F *h_NEvents_HLT40  = (TH1F*)fin->Get("NEvents_HLT40");
  TH1F *h_NEvents_HLT60  = (TH1F*)fin->Get("NEvents_HLT60");
  TH1F *h_NEvents_HLT80  = (TH1F*)fin->Get("NEvents_HLT80");
  TH1F *h_NEvents_HLT100 = (TH1F*)fin->Get("NEvents_HLT100");

  TH1F *h_NEvents_HLT40thresh  = (TH1F*)fin->Get("NEvents_HLT40thresh");
  TH1F *h_NEvents_HLT60thresh  = (TH1F*)fin->Get("NEvents_HLT60thresh");
  TH1F *h_NEvents_HLT80thresh  = (TH1F*)fin->Get("NEvents_HLT80thresh");
  TH1F *h_NEvents_HLT100thresh = (TH1F*)fin->Get("NEvents_HLT100thresh");

  TH1F *h_NEvents_isMB  = (TH1F*)fin->Get("NEvents_isMB");
  TH1F *h_NEvents_is40  = (TH1F*)fin->Get("NEvents_is40");
  TH1F *h_NEvents_is60  = (TH1F*)fin->Get("NEvents_is60");
  TH1F *h_NEvents_is80  = (TH1F*)fin->Get("NEvents_is80");
  TH1F *h_NEvents_is100 = (TH1F*)fin->Get("NEvents_is100");

  TH1F* hTrigComb=(TH1F*)fin->Get("TrigComb");

  
  std::cout<<std::endl;  
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<"--- DATA Event Count Report ---"<<std::endl;
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<std::endl;  

  std::cout<<"--- EVENT LOOP ---"<<std::endl;
  std::cout<<"# evts (description) = count / pscl weighted count" <<std::endl<<std::endl;
  std::cout<<"in file(s) opened            = " << h_NEvents->GetBinContent(1)          << " / " << h_NEvents->GetBinContent(2)<<std::endl; 
  std::cout<<"skipped                      = " << h_NEvents_skipped->GetBinContent(1) << " / " << h_NEvents_skipped->GetBinContent(2)<<std::endl;
  std::cout<<"read = in files - skipped    = " << h_NEvents_read->GetBinContent(1)    << " / " << h_NEvents_read->GetBinContent(2)<<std::endl;  
  std::cout<<"skimCuts  = read - !skim     = " << h_NEvents_skimCut->GetBinContent(1) << " / " << h_NEvents_skimCut->GetBinContent(2) << std::endl;//same
  std::cout<<"vzCuts    = skimCuts - !vz   = " << h_NEvents_vzCut->GetBinContent(1)   << " / " << h_NEvents_vzCut->GetBinContent(2)   << std::endl; //same
  std::cout<<"excl trig w/o quality = read - !trigd        = " << h_NEvents_trigd_1->GetBinContent(1) << " / " << h_NEvents_trigd_1->GetBinContent(2)<<std::endl;//non-duplicates
  std::cout<<"excl trig w/ quality = read - !qual - !trigd = " << h_NEvents_trigd_2->GetBinContent(1) << " / " << h_NEvents_trigd_2->GetBinContent(2)<<std::endl;//close to trigd_1
  std::cout << std::endl;
  
  std::cout<<std::endl;
  std::cout<<"Trigger Counts summary, post duplicate skip, post quality (unless otherwise noted)"<<std::endl;
  std::cout<<"Trigger Combination type = " << hTrigComb->GetTitle() << std::endl;
  std::cout<<std::endl;
  std::cout<<"Pass Trigger Inclusively; Evt Counts"<<std::endl;
  std::cout<<"# evts passing HLTMB      = " << h_NEvents_HLTMB ->GetBinContent(1) << " / " << h_NEvents_HLTMB ->GetBinContent(2) <<std::endl;
  std::cout<<"# evts passing HLT40      = " << h_NEvents_HLT40 ->GetBinContent(1) << " / " << h_NEvents_HLT40 ->GetBinContent(2) <<std::endl;
  std::cout<<"# evts passing HLT60      = " << h_NEvents_HLT60 ->GetBinContent(1) << " / " << h_NEvents_HLT60 ->GetBinContent(2) <<std::endl;
  std::cout<<"# evts passing HLT80      = " << h_NEvents_HLT80 ->GetBinContent(1) << " / " << h_NEvents_HLT80 ->GetBinContent(2) <<std::endl;
  std::cout<<"# evts passing HLT100     = " << h_NEvents_HLT100->GetBinContent(1) << " / " << h_NEvents_HLT100->GetBinContent(2) <<std::endl;
  std::cout<<std::endl;    
  std::cout<<"Pass Trigger + Threshold Inclusively; Evt Counts"<<std::endl;
  //std::cout<<"Trigger Thresholds for 40/60/80/100 = " << HLTthresh[0] <<"/ "<< HLTthresh[1]<<"/ "<< HLTthresh[2]<<"/ "<< HLTthresh[3]<<std::endl;
  std::cout<<"# evts passing HLT40 +threshold     = " << h_NEvents_HLT40thresh->GetBinContent(1)  << " / " << h_NEvents_HLT40thresh->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing HLT60 +threshold     = " << h_NEvents_HLT60thresh->GetBinContent(1)  << " / " << h_NEvents_HLT60thresh->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing HLT80 +threshold     = " << h_NEvents_HLT80thresh->GetBinContent(1)  << " / " << h_NEvents_HLT80thresh->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing HLT100+threshold     = " << h_NEvents_HLT100thresh->GetBinContent(1) << " / " << h_NEvents_HLT100thresh->GetBinContent(2) <<std::endl;
  std::cout<<std::endl;    
  std::cout<<"Pass Trigger Exclusively; Evt Counts"<<std::endl;
  std::cout<<"excl trigd w/o quality   = " << h_NEvents_trigd_1->GetBinContent(1) << " / " << h_NEvents_trigd_1->GetBinContent(2)<<std::endl;//non-duplicates
  std::cout<<"excl trigd w/ quality    = " << h_NEvents_trigd_2->GetBinContent(1) << " / " << h_NEvents_trigd_2->GetBinContent(2)<<std::endl;//same
  std::cout<<"# evts passing isMB      = " << h_NEvents_isMB->GetBinContent(1)  << " / " << h_NEvents_isMB->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing is40      = " << h_NEvents_is40->GetBinContent(1)  << " / " << h_NEvents_is40->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing is60      = " << h_NEvents_is60->GetBinContent(1)  << " / " << h_NEvents_is60->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing is80      = " << h_NEvents_is80->GetBinContent(1)  << " / " << h_NEvents_is80->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing is100     = " << h_NEvents_is100->GetBinContent(1) << " / " << h_NEvents_is100->GetBinContent(2) <<std::endl;
  
  return;

}


void printMCEventCountReport(TFile* fin){
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"printMCEventCountReport Called"<<std::endl<<std::endl;
  
  if(!fin){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  
  TH1F *h_NEvents         = (TH1F*)fin->Get("NEvents");              if(!h_NEvents        ) assert(false);
  TH1F *h_NEvents_read    = (TH1F*)fin->Get("NEvents_read");         if(!h_NEvents_read   ) assert(false);
  TH1F *h_NEvents_skimCut = (TH1F*)fin->Get("NEvents_skimCut");	     if(!h_NEvents_skimCut) assert(false);
  TH1F *h_NEvents_vzCut   = (TH1F*)fin->Get("NEvents_vzCut");	     if(!h_NEvents_vzCut  ) assert(false);

  std::cout<<std::endl;  
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<"--- MC Event Count Report ---"<<std::endl;
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<std::endl;  
    
  std::cout<<"--- EVENT LOOP ---"<<std::endl;
  std::cout<<"# evts (description) = count / pt-hat weighted count" <<std::endl<<std::endl;  

  std::cout << std::endl;
  std::cout<<"in file(s) opened               = " << h_NEvents->GetBinContent(1)    << " / " << h_NEvents->GetBinContent(2)<<std::endl;
  std::cout<<"read from file(s) (eq to above) = " << h_NEvents_read->GetBinContent(1)    << " / " << h_NEvents_read->GetBinContent(2)<<std::endl;
  std::cout<<"skimCuts  = read - !skim        = " << h_NEvents_skimCut->GetBinContent(1) << " / " << h_NEvents_skimCut->GetBinContent(2) << std::endl;//same
  std::cout<<"vzCuts    = skimCuts - !vz      = " << h_NEvents_vzCut->GetBinContent(1)   << " / " << h_NEvents_vzCut->GetBinContent(2)   << std::endl; //same
  std::cout << std::endl;

  return;

}



void printDataJetCountReport(TFile* fin){
  std::cout<<std::endl;  
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<"--- DATA Jet Count Report ---"<<std::endl;
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<std::endl;  
  

  TH1F* h_NJets           =(TH1F*)fin->Get("NJets");	    
  TH1F* h_NJets_jtptCut   =(TH1F*)fin->Get("NJets_jtptCut");  
  TH1F* h_NJets_jtptCut_Hi   =(TH1F*)fin->Get("NJets_jtptCut_Hi");  
  TH1F* h_NJets_jtetaCut1 =(TH1F*)fin->Get("NJets_jtetaCut1");
  TH1F* h_NJets_jtetaCut2 =(TH1F*)fin->Get("NJets_jtetaCut2");
  TH1F* h_NJets_JetIDCut =(TH1F*)fin->Get("NJets_JetIDCut");

  std::cout<<"--- JET LOOP ---"<<std::endl;
  std::cout<<"# jets (desc) = count / pscl weighted count" <<std::endl<<std::endl;  

  std::cout<<std::endl;    
  std::cout<<"# of Jets read                   = " << h_NJets->GetBinContent(1) << " / " << h_NJets->GetBinContent(2) <<std::endl; 
  std::cout<<"# of Jets w/ jtpt > jtptCut      = " << h_NJets_jtptCut->GetBinContent(1) << " / " << h_NJets_jtptCut->GetBinContent(2) <<std::endl;   
  std::cout<<"# of Jets w/ jtpt > jtptCut_Hi   = " << h_NJets_jtptCut_Hi->GetBinContent(1) << " / " << h_NJets_jtptCut_Hi->GetBinContent(2) <<std::endl;   
  std::cout<<"# of Jets w/ |jteta|<4.7         = " << h_NJets_jtetaCut1->GetBinContent(1) << " / " << h_NJets_jtetaCut1->GetBinContent(2) <<std::endl; 
  std::cout<<"# of Jets w/ passing jtEtaCut(s) = " << h_NJets_jtetaCut2->GetBinContent(1) << " / " << h_NJets_jtetaCut2->GetBinContent(2) <<std::endl; 
  if(h_NJets_JetIDCut)
  std::cout<<"# of Jets w/ jetIDCut            = " << h_NJets_JetIDCut->GetBinContent(1) << " / " << h_NJets_JetIDCut->GetBinContent(2) <<std::endl; 
  std::cout << std::endl;
  
  return;
}

void printMCJetCountReport(TFile* fin){

  TH1F* h_NJets		      =(TH1F*)fin->Get("NJets"		    );
  TH1F* h_NJets_matched	      =(TH1F*)fin->Get("NJets_matched"	    );
  TH1F* h_NJets_deltaRCut     =(TH1F*)fin->Get("NJets_deltaRCut"    );
  TH1F* h_NJets_jtptCut	      =(TH1F*)fin->Get("NJets_jtptCut"	    );
  TH1F* h_NJets_jtptCut_Hi    =(TH1F*)fin->Get("NJets_jtptCut_Hi"   );
  TH1F* h_NJets_jtetaCut1     =(TH1F*)fin->Get("NJets_jtetaCut1"    );
  TH1F* h_NJets_jtetaCut2     =(TH1F*)fin->Get("NJets_jtetaCut2"    );
  TH1F* h_NJets_genjtptCut    =(TH1F*)fin->Get("NJets_genjtptCut"   );
  TH1F* h_NJets_genjtptCut_Hi =(TH1F*)fin->Get("NJets_genjtptCut_Hi");
  TH1F* h_NJets_JetIDCut      =(TH1F*)fin->Get("NJets_JetIDCut"     );

  std::cout<<std::endl;  
  std::cout<<"---------------------------"<<std::endl;
  std::cout<<"--- MC Jet Count Report ---"<<std::endl;
  std::cout<<"---------------------------"<<std::endl;
  std::cout<<std::endl;  
  
  std::cout<<"--- JET LOOP ---"<<std::endl;
  std::cout<<"# jets (desc) = count / pt-hat weighted count" <<std::endl<<std::endl;  

  std::cout<<std::endl;
  std::cout<<"# of Jets read                     = " << h_NJets->GetBinContent(1)               << " / " << h_NJets->GetBinContent(2) <<std::endl;
  std::cout<<"# of GEN-RECO Matched Jets         = " << h_NJets_matched->GetBinContent(1)       << " / " << h_NJets_matched->GetBinContent(2) <<std::endl;
  std::cout<<"# of Matched Jets meeting dR cut   = " << h_NJets_deltaRCut->GetBinContent(1)     << " / " << h_NJets_deltaRCut->GetBinContent(2) <<std::endl;  
  std::cout<<"# of Jets w/ jtpt > jtptCut        = " << h_NJets_jtptCut->GetBinContent(1)       << " / " << h_NJets_jtptCut->GetBinContent(2) <<std::endl;
  std::cout<<"# of Jets w/ jtpt > jtptCut_Hi     = " << h_NJets_jtptCut_Hi->GetBinContent(1)    << " / " << h_NJets_jtptCut_Hi->GetBinContent(2) <<std::endl;
  std::cout<<"# of Jets w/ |jteta|<4.7           = " << h_NJets_jtetaCut1->GetBinContent(1)     << " / " << h_NJets_jtetaCut1->GetBinContent(2) <<std::endl;
  std::cout<<"# of Jets w/ passing jtEtaCut(s)   = " << h_NJets_jtetaCut2->GetBinContent(1)     << " / " << h_NJets_jtetaCut2->GetBinContent(2) <<std::endl;
  std::cout<<"# of Jets w/ genjtpt>genjtptCut    = " << h_NJets_genjtptCut->GetBinContent(1)    << " / " << h_NJets_genjtptCut->GetBinContent(2) <<std::endl;
  std::cout<<"# of Jets w/ genjtpt<genjtptCut_Hi = " << h_NJets_genjtptCut_Hi->GetBinContent(1) << " / " << h_NJets_genjtptCut_Hi->GetBinContent(2) <<std::endl;
  if(h_NJets_JetIDCut)
  std::cout<<"# of Jets w/ jetIDCut              = " << h_NJets_JetIDCut->GetBinContent(1) << " / " << h_NJets_JetIDCut->GetBinContent(2) <<std::endl;   
  std::cout << std::endl;
  
  return;
}




long double computeEffLumi(TFile* finData){
  bool funcDebug=false;
  if(funcDebug)std::cout<<"computeEffLumi called."<<std::endl;
  long double effIntgrtdLumi_vz=1., LumiEff_vz=1.;// eff-> effective, Eff->efficiency
  
  TH1F *h_NEvents_vzCut   = (TH1F*)finData->Get("NEvents_vzCut");
  TH1F *h_NEvents_read    = (TH1F*)finData->Get("NEvents_read"); 

  std::cout<<std::endl<<"integrated luminosity     ="<<intgrtdLumi<<std::endl;               
  

  std::cout<<"total # evts post vzCut = " << h_NEvents_vzCut->GetBinContent(1) << std::endl;
  std::cout<<"total # evts read = " << h_NEvents_read->GetBinContent(1) << std::endl;

  std::cout<<"total # (weighted) evts post vzCut = " << h_NEvents_vzCut->GetBinContent(2) << std::endl;     
  std::cout<<"total # (weighted) evts read = " << h_NEvents_read->GetBinContent(2) << std::endl;     
  
  LumiEff_vz = h_NEvents_vzCut->GetBinContent(2)/h_NEvents_read->GetBinContent(2);                                    
  //LumiEff_vz = h_NEvents_vzCut->GetBinContent(1)/h_NEvents_read->GetBinContent(1);                                    
  std::cout<<"lumi efficiency, vz cuts             ="<<LumiEff_vz<<std::endl;                                   
  effIntgrtdLumi_vz=intgrtdLumi*LumiEff_vz;                                                                   
  std::cout<<"effective ingrtdLumi for dataset+cuts="<<effIntgrtdLumi_vz<<std::endl;                            
  
  return effIntgrtdLumi_vz;
}







///////////////////////////// TRK MULT VERSION(S)
void getEventCountsTrkMult(TFile* fin, bool isData){
  bool funcDebug=false;
  std::cout<<" grabbing raw event counts "<<std::endl;   
  
  
  if(!fin){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  
  if(funcDebug)std::cout<<"opening hists"<<std::endl;
  TH1F *h_NEvents         = (TH1F*)fin->Get("NEvents");
  TH1F *h_NEvents_read    = (TH1F*)fin->Get("NEvents_read");      
  TH1F *h_NEvents_skimCut = (TH1F*)fin->Get("NEvents_skimCut");
  TH1F *h_NEvents_vzCut   = (TH1F*)fin->Get("NEvents_vzCut");
  TH1F *h_NEvents_trkCuts_1   = (TH1F*)fin->Get("NEvents_trkCuts_1");
  TH1F *h_NEvents_trkCuts_2   = (TH1F*)fin->Get("NEvents_trkCuts_2");
  TH1F *h_NEvents_trkCuts_3   = (TH1F*)fin->Get("NEvents_trkCuts_3");
  
  TH1F *h_NEvents_skipped = NULL;
  TH1F *h_NEvents_trigd=NULL;
  TH1F *h_NEvents_jet100=NULL,*h_NEvents_jet80=NULL,*h_NEvents_jet60=NULL,*h_NEvents_jet40=NULL;
  
  if(isData)
    { h_NEvents_skipped=(TH1F*)fin->Get("NEvents_skipped");      
      h_NEvents_trigd=(TH1F*)fin->Get("NEvents_trigd");      
      h_NEvents_jet40=(TH1F*)fin->Get("NEvents_jet40Trigd");      
      h_NEvents_jet60=(TH1F*)fin->Get("NEvents_jet60Trigd");      
      h_NEvents_jet80=(TH1F*)fin->Get("NEvents_jet80Trigd");      
      h_NEvents_jet100=(TH1F*)fin->Get("NEvents_jet100Trigd");      
    }
  
  if(isData)std::cout<<std::endl<<"/// Data Evt-Count Summary ///"<<std::endl<<std::endl;      
  else std::cout<<std::endl<<"/// MC Evt-Count Summary ///"<<std::endl<<std::endl;      
  
  std::cout<<"Total Num of Events in file(s) opened       = " <<
    h_NEvents->GetEntries()<<std::endl;
  if(isData)
    std::cout<<"Total Num of Events skipped from those file(s) = " <<
      h_NEvents_skipped->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read from those file(s) = " <<
	h_NEvents_read->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read passing skimCuts   = " <<
    h_NEvents_skimCut->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read passing vzCuts and skimCuts    = " <<
    h_NEvents_vzCut->GetEntries()<<std::endl<<std::endl;
  
  if(isData){
    std::cout<<"Total Num of Events trigd in those file(s) = " <<
      h_NEvents_trigd->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events jet40 in those file(s) = " <<
      h_NEvents_jet40->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events jet60 in those file(s) = " <<
      h_NEvents_jet60->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events jet80 in those file(s) = " <<
      h_NEvents_jet80->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events jet100 in those file(s) = " <<
      h_NEvents_jet100->GetEntries()<<std::endl;}
    
  std::cout<<"Total Num of Events read passing trkHiPur Cut (#1) = " <<
    h_NEvents_trkCuts_1->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read passing vtx spacing Cut (#2) = " <<
    h_NEvents_trkCuts_2->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read passing oneVtxInRange Cut (#3) = " <<
    h_NEvents_trkCuts_3->GetEntries()<<std::endl;
  
  return;

}

long double computeEffLumiTrkMult(TFile* finData){

  long double effIntgrtdLumi_vz=1., LumiEff_vz=1.;// eff-> effective, Eff->efficiency
  
  std::cout<<std::endl<<"dataset integrated Luminosity (microbarns) ="<<intgrtdLumi<<std::endl;               
  TH1F *h_NEvents_vzCut   = (TH1F*)finData->Get("NEvents_trkCuts_3");
  TH1F *h_NEvents_read    = (TH1F*)finData->Get("NEvents_read");      
  
  LumiEff_vz = h_NEvents_vzCut->GetEntries()/h_NEvents_read->GetEntries();                                    
  effIntgrtdLumi_vz=intgrtdLumi*LumiEff_vz;                                                                   
  
  std::cout<<std::endl<<"lumi efficiency, vz cuts="<<LumiEff_vz<<std::endl;                                   
  std::cout<<"effective ingrtdLumi for dataset+cuts="<<effIntgrtdLumi_vz<<std::endl;                          
  //LumiEff_vz=0.999;//temp
  return effIntgrtdLumi_vz;
}














/////////////////////////#####################################/////////////////////////
/////////////////////////#####################################/////////////////////////
/////////////////////////#####################################/////////////////////////
///////////////////////////// DEPRECATED VERSION(S)

void printJetCountReport(TFile* fin, bool isData){ //deprecated
  //void printJetCountReport(TFile* fin){
  if(isData){
    std::cout<<std::endl;
    std::cout<<"------------- DATA -----------"<<std::endl;
    std::cout<<"--- Jet-Event Loop Summary ---"<<std::endl;
    std::cout<<"------------------------------"<<std::endl;
    std::cout<<std::endl;
  }
  else {
    std::cout<<std::endl;
    std::cout<<"-------------  MC  -----------"<<std::endl;
    std::cout<<"--- Jet-Event Loop Summary ---"<<std::endl;
    std::cout<<"------------------------------"<<std::endl;
    std::cout<<std::endl;
  }
  TH1F* h_NJets           =(TH1F*)fin->Get("NJets");	    
  std::cout<<"# jets read                   = / " << h_NJets->GetBinContent(1)           << " / " << h_NJets->GetBinContent(2) <<std::endl; 
  if(isData){
    TH1F* h_NJets_jtptCut   =(TH1F*)fin->Get("NJets_jtptCut");  
    TH1F* h_NJets_jtetaCut1 =(TH1F*)fin->Get("NJets_jtetaCut1");
    TH1F* h_NJets_jtetaCut2 =(TH1F*)fin->Get("NJets_jtetaCut2");
    
    std::cout<<"# jets w jtpt > jtptCut       = / " << h_NJets_jtptCut->GetBinContent(1)   << " / " << h_NJets_jtptCut->GetBinContent(2) <<std::endl;   
    std::cout<<"# jets w |jteta|<4.7          = / " << h_NJets_jtetaCut1->GetBinContent(1) << " / " << h_NJets_jtetaCut1->GetBinContent(2) <<std::endl; 
    std::cout<<"# jets w passing jtEtaCut(s)  = / " << h_NJets_jtetaCut2->GetBinContent(1) << " / " << h_NJets_jtetaCut2->GetBinContent(2) <<std::endl; 
  }
  else { //isMC
    TH1F* h_NJets_kmatCut=(TH1F*)fin->Get("NJets_kmatCut");
    std::cout<<"# jets passing kinematic cuts = / " << h_NJets_kmatCut->GetBinContent(1)   << " / " << h_NJets_kmatCut->GetBinContent(2) <<std::endl;   
  }
  
  TH1F* h_NJets_JetIDCut =(TH1F*)fin->Get("NJets_JetIDCut");
  if((bool)h_NJets_JetIDCut)
    std::cout<<"# jets w jetIDCut             = / " << h_NJets_JetIDCut->GetBinContent(1) << " / " << h_NJets_JetIDCut->GetBinContent(2) <<std::endl; 
  
  return;
}








void getEventCounts(TFile* fin, bool isData){ //deprecated
  bool funcDebug=false;
  std::cout<<" grabbing raw event counts "<<std::endl;   
  
  if(!fin){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  
  if(funcDebug)std::cout<<"opening hists"<<std::endl;
  TH1F *h_NEvents         = (TH1F*)fin->Get("NEvents");
  TH1F *h_NEvents_read    = (TH1F*)fin->Get("NEvents_read");      
  TH1F *h_NEvents_skimCut = (TH1F*)fin->Get("NEvents_skimCut");
  TH1F *h_NEvents_vzCut   = (TH1F*)fin->Get("NEvents_vzCut");
  
//  if(!h_NEvents        ||
//     !h_NEvents_read   ||
//     !h_NEvents_skimCut||
//     !h_NEvents_vzCut   ) assert(false);
  if(!h_NEvents        ) assert(false);
  if(!h_NEvents_read   ) assert(false);
  if(!h_NEvents_skimCut) assert(false);
  if(!h_NEvents_vzCut  ) assert(false);
  
  TH1F *h_NEvents_skipped = NULL;    
  if(isData){
    h_NEvents_skipped=(TH1F*)fin->Get("NEvents_skipped");      
    std::cout<<std::endl<<"/// Data Evt-Count Summary ///"<<std::endl<<std::endl;      
  }
  else std::cout<<std::endl<<"/// MC Evt-Count Summary ///"<<std::endl<<std::endl;      
  
  std::cout<<"Total Num of Events in file(s) opened       = " <<
    h_NEvents->GetEntries()<<std::endl;
  if(isData)  std::cout<<"Total Num of Events skipped from those file(s) = " <<
		h_NEvents_skipped->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read from those file(s) = " <<
    h_NEvents_read->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read passing skimCuts   = " <<
    h_NEvents_skimCut->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read passing vzCuts and skimCuts    = " <<
    h_NEvents_vzCut->GetEntries()<<std::endl<<std::endl;
  

//  TH1F *h_NEvents_trigd=NULL;
//  TH1F *h_NEvents_jet100=NULL,*h_NEvents_jet80=NULL,*h_NEvents_jet60=NULL,*h_NEvents_jet40=NULL;  
//  if(isData){
//    
//    h_NEvents_trigd=(TH1F*)fin->Get("NEvents_trigd");      
//
//    h_NEvents_jet40=(TH1F*)fin->Get("NEvents_jet40Trigd");      
//    h_NEvents_jet60=(TH1F*)fin->Get("NEvents_jet60Trigd");      
//    h_NEvents_jet80=(TH1F*)fin->Get("NEvents_jet80Trigd");      
//    h_NEvents_jet100=(TH1F*)fin->Get("NEvents_jet100Trigd");      
//    
//
//    std::cout<<"Total Num of Events trigd in those file(s) = " <<
//      h_NEvents_trigd->GetEntries()<<std::endl;
//   a ssert(false);      
//    std::cout<<"Total Num of Events jet40 in those file(s) = " <<
//      h_NEvents_jet40->GetEntries()<<std::endl;
//    std::cout<<"Total Num of Events jet60 in those file(s) = " <<
//      h_NEvents_jet60->GetEntries()<<std::endl;
//    std::cout<<"Total Num of Events jet80 in those file(s) = " <<
//      h_NEvents_jet80->GetEntries()<<std::endl;
//
//    if(h_NEvents_jet100)//not always made anymore
//      std::cout<<"Total Num of Events jet100 in those file(s) = " <<
//	h_NEvents_jet100->GetEntries()<<std::endl;
//  }
  
  return;

}
















void old_printDataEventCountReport(TFile* fin){//deprecated
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"printDataEventCountReport Called"<<std::endl<<std::endl;
  
  if(!fin){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  
  std::cout<<std::endl;  
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<"--- DATA Event Count Report ---"<<std::endl;
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<std::endl;  
  
  std::cout<<"--- EVENT LOOP ---"<<std::endl;
  std::cout<<"# evts (description) = / count / pscl weighted count" <<std::endl<<std::endl;

  TH1F *h_NEvents         = (TH1F*)fin->Get("NEvents");
  TH1F *h_NEvents_trigd_1=(TH1F*)fin->Get("NEvents_trigd_1");
  TH1F *h_NEvents_nopscl =(TH1F*)fin->Get("NEvents_nopscl");
  TH1F *h_NEvents_jet100=(TH1F*)fin->Get("NEvents_jet100Trigd"),    *h_NEvents_jet80=(TH1F*)fin->Get("NEvents_jet80Trigd"),    *h_NEvents_jet60=(TH1F*)fin->Get("NEvents_jet60Trigd"),    *h_NEvents_jet40=(TH1F*)fin->Get("NEvents_jet40Trigd");
  
  std::cout<<"# evts in file(s) opened    = / "    <<h_NEvents->GetBinContent(1)        << " / "<< h_NEvents->GetBinContent(2)<<std::endl; 
  std::cout<<"# evts trigd in those files = / "    <<h_NEvents_trigd_1->GetBinContent(1) << " / "<< h_NEvents_trigd_1->GetBinContent(2)<<std::endl;
  std::cout<<"# evts w pscl==1            = / " << h_NEvents_nopscl->GetBinContent(1) << " / " << h_NEvents_nopscl->GetBinContent(2) <<std::endl << std::endl;
  std::cout<<"# evts passing is40  = / " << h_NEvents_jet40->GetBinContent(1)  << " / " << h_NEvents_jet40->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing is60  = / " << h_NEvents_jet60->GetBinContent(1)  << " / " << h_NEvents_jet60->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing is80  = / " << h_NEvents_jet80->GetBinContent(1)  << " / " << h_NEvents_jet80->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing is100 = / " << h_NEvents_jet100->GetBinContent(1) << " / " << h_NEvents_jet100->GetBinContent(2) <<std::endl<<std::endl;  
  
  TH1F *h_NEvents_read    = (TH1F*)fin->Get("NEvents_read");        
  TH1F *h_NEvents_skipped = (TH1F*)fin->Get("NEvents_skipped");  
  TH1F *h_NEvents_trigd_2=(TH1F*)fin->Get("NEvents_trigd_2");
  
  std::cout<<"# evts in file(s) opened           = / "    <<h_NEvents->GetBinContent(1)        << " / "<< h_NEvents->GetBinContent(2)<<std::endl; //repeat on purpose
  std::cout<<"# evts read from those file(s)     = / " << h_NEvents_read->GetBinContent(1)    << " / " << h_NEvents_read->GetBinContent(2)<<std::endl;
  std::cout<<"# evts skipped from those file(s)  = / " << h_NEvents_skipped->GetBinContent(1) << " / " << h_NEvents_skipped->GetBinContent(2)<<std::endl;
  std::cout<<"# evts passing a trigger post-skip = / " << h_NEvents_trigd_2->GetBinContent(1) << " / " << h_NEvents_trigd_2->GetBinContent(2) <<std::endl << std::endl;
  
  TH1F *h_NEvents_skimCut = (TH1F*)fin->Get("NEvents_skimCut");
  TH1F *h_NEvents_vzCut   = (TH1F*)fin->Get("NEvents_vzCut");
  
  std::cout<<"# evts passing skimCuts  = / " << h_NEvents_skimCut->GetBinContent(1) << " / " << h_NEvents_skimCut->GetBinContent(2) << std::endl;
  std::cout<<"# evts passing vzCuts    = / " << h_NEvents_vzCut->GetBinContent(1)   << " / " << h_NEvents_vzCut->GetBinContent(2)   << std::endl<<std::endl;
  
  
  TH1F *hpp_HLT40trgEta=(TH1F*)fin->Get("isHLT40_trgEta"), *hpp_HLT60trgEta=(TH1F*)fin->Get("isHLT60_trgEta"), *hpp_HLT80trgEta=(TH1F*)fin->Get("isHLT80_trgEta"),   *hpp_HLT100trgEta=(TH1F*)fin->Get("isHLT100_trgEta"),*hpp_HLTCombtrgPt=(TH1F*)fin->Get("HLTComb_trgPt") ;
  TH1F *hpp_IncHLT40trgEta=(TH1F*)fin->Get("IncHLT40_trgEta"), *hpp_IncHLT60trgEta=(TH1F*)fin->Get("IncHLT60_trgEta"), *hpp_IncHLT80trgEta=(TH1F*)fin->Get("IncHLT80_trgEta"),   *hpp_IncHLT100trgEta=(TH1F*)fin->Get("IncHLT100_trgEta");
  
  std::cout<<"# evts passing trig  = / " << hpp_HLTCombtrgPt ->GetEntries() << " / " <<  hpp_HLTCombtrgPt->Integral()<<std::endl;
  std::cout<<"# evts passing is40  = / " << hpp_HLT40trgEta ->GetEntries()  << " / " <<  hpp_HLT40trgEta->Integral()<<std::endl;
  std::cout<<"# evts passing is60  = / " << hpp_HLT60trgEta ->GetEntries()  << " / " <<  hpp_HLT60trgEta->Integral()<<std::endl;
  std::cout<<"# evts passing is80  = / " << hpp_HLT80trgEta ->GetEntries()  << " / " <<  hpp_HLT80trgEta->Integral()<<std::endl;
  std::cout<<"# evts passing is100 = / " << hpp_HLT100trgEta->GetEntries()  << " / " <<  hpp_HLT100trgEta->Integral()<<std::endl<<std::endl;  
  
  std::cout<<"# evts passing each trigger with no exclusion" <<std::endl;
  std::cout<<"# evts passing HLT40  = / " << hpp_IncHLT40trgEta ->GetEntries()  << " / " <<  hpp_IncHLT40trgEta ->Integral()<<std::endl;
  std::cout<<"# evts passing HLT60  = / " << hpp_IncHLT60trgEta ->GetEntries()  << " / " <<  hpp_IncHLT60trgEta ->Integral()<<std::endl;
  std::cout<<"# evts passing HLT80  = / " << hpp_IncHLT80trgEta ->GetEntries()  << " / " <<  hpp_IncHLT80trgEta ->Integral()<<std::endl;
  std::cout<<"# evts passing HLT100 = / " << hpp_IncHLT100trgEta->GetEntries()  << " / " <<  hpp_IncHLT100trgEta->Integral()<<std::endl<<std::endl;  
  
  // don't seem to need this... but just dont delete, you never know. 
  //std::cout<<"--- EVENTS W JETS ---"<<std::endl;
  //std::cout<<"# evts (desc) = / count / pscl weighted count" <<std::endl<<std::endl;  
  //TH1F* h_NEvents_withJets           =(TH1F*)fin->Get("NEvents_withJets");	    
  //TH1F* h_NEvents_withJets_jtptCut   =(TH1F*)fin->Get("NEvents_withJets_jtptCut");  
  //TH1F* h_NEvents_withJets_jtetaCut1 =(TH1F*)fin->Get("NEvents_withJets_jtetaCut1");
  //TH1F* h_NEvents_withJets_jtetaCut2 =(TH1F*)fin->Get("NEvents_withJets_jtetaCut2");
  //TH1F* h_NEvents_withJets_JetIDCut = (TH1F*)fin->Get("NEvents_withJets_JetIDCut");
  //std::cout<<"# evts w min. 1 jet                   = / " << h_NEvents_withJets->GetBinContent(1)           << " / " << h_NEvents_withJets->GetBinContent(2) <<std::endl; 
  //std::cout<<"# evts w min. 1 jet w jtpt > jtPtCut = / " << h_NEvents_withJets_jtptCut->GetBinContent(1)   << " / " << h_NEvents_withJets_jtptCut->GetBinContent(2) << std::endl;
  //std::cout<<"# evts w min. 1 jet w |jteta| < 4.7  = / " << h_NEvents_withJets_jtetaCut1->GetBinContent(1) << " / " << h_NEvents_withJets_jtetaCut1->GetBinContent(2) << std::endl; 
  //std::cout<<"# evts w min. 1 jet w jtEtaCut(s)    = / " << h_NEvents_withJets_jtetaCut2->GetBinContent(1) << " / " << h_NEvents_withJets_jtetaCut2->GetBinContent(2) << std::endl;    
  //if((bool)h_NEvents_withJets_JetIDCut)
  //  std::cout<<"# evts w min. 1 jet w JetIDCut       = / " << h_NEvents_withJets_JetIDCut->GetBinContent(1) << " / " << h_NEvents_withJets_JetIDCut->GetBinContent(2) << std::endl;         
  //std::cout<<std::endl;  
  
  return;

}


void old_printMCEventCountReport(TFile* fin){//deprecated
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"printMCEventCountReport Called"<<std::endl<<std::endl;
  
  if(!fin){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  
  std::cout<<std::endl;  
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<"--- MC Event Count Report ---"<<std::endl;
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<std::endl;  
  
  TH1F *h_NEvents         = (TH1F*)fin->Get("NEvents");              if(!h_NEvents        ) assert(false);
  TH1F *h_NEvents_read    = (TH1F*)fin->Get("NEvents_read");         if(!h_NEvents_read   ) assert(false);
  TH1F *h_NEvents_skimCut = (TH1F*)fin->Get("NEvents_skimCut");	     if(!h_NEvents_skimCut) assert(false);
  TH1F *h_NEvents_vzCut   = (TH1F*)fin->Get("NEvents_vzCut");	     if(!h_NEvents_vzCut  ) assert(false);
  
  std::cout<<"--- EVENT LOOP ---"<<std::endl;
  //std::cout<<"# evts (description) = / count / pt-hat weighted count" <<std::endl<<std::endl;//TODO
  std::cout<<"# evts (description) =  count " <<std::endl<<std::endl;
  
  std::cout<<"Total Num of Events in file(s) opened                   = " <<    h_NEvents->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read from those file(s)             = " <<    h_NEvents_read->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read passing skimCuts               = " <<    h_NEvents_skimCut->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read passing vzCuts and skimCuts    = " <<    h_NEvents_vzCut->GetEntries()<<std::endl<<std::endl;
  
  return;

}



void old_printDataJetCountReport(TFile* fin){//deprecated
  std::cout<<std::endl;  
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<"--- DATA Jet Count Report ---"<<std::endl;
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<std::endl;  
  

  std::cout<<"--- JET LOOP ---"<<std::endl;
  std::cout<<"# jets (desc) = / count / pscl weighted count" <<std::endl<<std::endl;  

  TH1F* h_NJets           =(TH1F*)fin->Get("NJets");	    
  TH1F* h_NJets_jtptCut   =(TH1F*)fin->Get("NJets_jtptCut");  
  TH1F* h_NJets_jtetaCut1 =(TH1F*)fin->Get("NJets_jtetaCut1");
  TH1F* h_NJets_jtetaCut2 =(TH1F*)fin->Get("NJets_jtetaCut2");
  TH1F* h_NJets_JetIDCut =(TH1F*)fin->Get("NJets_JetIDCut");
  
  std::cout<<"# jets read                   = / " << h_NJets->GetBinContent(1)           << " / " << h_NJets->GetBinContent(2) <<std::endl; 
  std::cout<<"# jets w jtpt > jtptCut       = / " << h_NJets_jtptCut->GetBinContent(1)   << " / " << h_NJets_jtptCut->GetBinContent(2) <<std::endl;   
  std::cout<<"# jets w |jteta|<4.7          = / " << h_NJets_jtetaCut1->GetBinContent(1) << " / " << h_NJets_jtetaCut1->GetBinContent(2) <<std::endl; 
  std::cout<<"# jets w passing jtEtaCut(s)  = / " << h_NJets_jtetaCut2->GetBinContent(1) << " / " << h_NJets_jtetaCut2->GetBinContent(2) <<std::endl; 
  if((bool)h_NJets_JetIDCut){
  std::cout<<"# jets w jetIDCut             = / " << h_NJets_JetIDCut->GetBinContent(1) << " / " << h_NJets_JetIDCut->GetBinContent(2) <<std::endl;   }
  
  return;
}

void old_printMCJetCountReport(TFile* fin){//deprecated

  std::cout<<std::endl;  
  std::cout<<"---------------------------"<<std::endl;
  std::cout<<"--- MC Jet Count Report ---"<<std::endl;
  std::cout<<"---------------------------"<<std::endl;
  std::cout<<std::endl;  
  
  std::cout<<"--- JET LOOP ---"<<std::endl;
  std::cout<<"# jets (desc) = / count / pscl weighted count" <<std::endl<<std::endl;  
  
  TH1F* h_NJets           =(TH1F*)fin->Get("NJets");	    
  TH1F* h_NJets_kmatCut=(TH1F*)fin->Get("NJets_kmatCut");
  TH1F* h_NJets_JetIDCut =(TH1F*)fin->Get("NJets_JetIDCut");
  
  std::cout<<"# jets read                   = / " << h_NJets->GetBinContent(1)           << " / " << h_NJets->GetBinContent(2) <<std::endl; 
  std::cout<<"# jets passing kinematic cuts = / " << h_NJets_kmatCut->GetBinContent(1)   << " / " << h_NJets_kmatCut->GetBinContent(2) <<std::endl;   
  if((bool)h_NJets_JetIDCut){
  std::cout<<"# jets w jetIDCut             = / " << h_NJets_JetIDCut->GetBinContent(1) << " / " << h_NJets_JetIDCut->GetBinContent(2) <<std::endl;   }
  
  return;
}
