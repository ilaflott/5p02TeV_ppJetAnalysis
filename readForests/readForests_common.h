////////// (initializa/declara)tions //////////
// C++, C, etc.
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

// ROOTSYS
#include <TSystem.h>
#include <TProfile.h>
#include <TStyle.h>
// I/O
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
// utilities
#include <TMath.h>
#include <TStopwatch.h>
#include <TF1.h>
//#include <TClonesArray.h>
// histos
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>


void divideBinWidth(TH1 *h){
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();++i){//binsX loop
    Double_t val=h->GetBinContent(i);
    Double_t valErr=h->GetBinError(i);
    val/=h->GetBinWidth(i);
    valErr/=h->GetBinWidth(i);
    h->SetBinContent(i,val);
    h->SetBinError(i,valErr);
  }//end nbinsX loop
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  return;
}


double trigComb(bool *trgDec, int *treePrescl, double triggerPt){
  double weight_eS=0.;
  if(trgDec[0] && triggerPt>=40.  && triggerPt<60. ) weight_eS=treePrescl[0];
  if(trgDec[1] && triggerPt>=60.  && triggerPt<80. ) weight_eS=treePrescl[1];
  if(trgDec[2] && triggerPt>=80.  && triggerPt<100.) weight_eS=treePrescl[2];
  if(trgDec[3] && triggerPt>=100. )                  weight_eS=treePrescl[3];

  return weight_eS;
}

float deltaphi(float phi1, float phi2){
  float pi=TMath::Pi();
  float dphi=TMath::Abs(phi1-phi2);
  if(dphi > pi)dphi -=2.*pi;
  return TMath::Abs(dphi);
}







const float pi=TMath::Pi();
const float twoPi=2.*pi;
const float piOvTwo=pi/2.;
const float threePiOvTwo=(3.*pi/2.);
float deltaphi_jettrk(float deltaphi){
  //float pi=TMath::Pi();
  
  if(deltaphi > (threePiOvTwo) )    
    deltaphi= deltaphi - twoPi ;
  else if (deltaphi < (-1*piOvTwo) )
    deltaphi = twoPi + deltaphi;
  else 
    deltaphi = deltaphi;
  
  return deltaphi;
}



//// STRINGS
// ------------------------------------------------------------------------------------------------------------

// eta width
const std::string etaWidth="20_eta_20";

// variable names for QA Plots, based on names in ak4PFJetAnalyzer
const std::string var[]={
  //jets, 0-3=4 vars
  "jtpt"  , "rawpt" ,
  "jteta" , "jtphi" ,
  //jet constituents, 4-30=27 vars
  "trkN", "trkSum", "trkMax", //trk
  "trkHardN", "trkHardSum",
  "phN", "phSum", "phMax",    //ph
  "phHardN", "phHardSum",
  "chN", "chSum", "chMax",    //ch
  "chHardN", "chHardSum",
  "neN" ,  "neSum" , "neMax" , //ne
  "eN"  ,  "eSum"  , "eMax"  ,    //e
  "muN" ,  "muSum" , "muMax" , //mu
  "neuMult", "chMult", "numConst",
  // dijets, 31-35=5 vars
  "Aj" , "xj" , "dphi",    //JERS might not like me
  "leadJetPt", "subleadJetPt"
};
const int N_vars=sizeof(var)/sizeof(std::string);
const int jtInd=0, jtConInd=4, dijtInd=31;

const std::string tupelvar[]={
  "Pt",
  "Eta",
  "Phi",
  "E",
  "Id",
  "RawPt",
  "RawE",
  "HfHadE",
  "HfEmE",
  "ChHadFrac",
  "NeutralHadAnHfFrac",
  "ChEmFrac",
  "NeutralEmFrac",
  "ChMult",
  "ConstCnt",
  "NeuMult"
};
const int N_tupelvars=sizeof(tupelvar)/sizeof(std::string);

const int tupel_nbins[N_tupelvars]={
2000  ,//  "Pt",
517*2,//  "Eta",
315*2,//  "Phi",
2000,//  "E",
4,//  "Id",
2000,//  "RawPt",
2000,//  "RawE",
2000,//  "HfHadE",
2000,//  "HfEmE",
200,//  "ChHadFrac",
200,//  "NeutralHadAnHfFrac",
200,//  "ChEmFrac",
200,//  "NeutralEmFrac",
100,//  "ChMult",
100,//  "ConstCnt"  
100 //"NeuMult"
};
const float tupel_xlo[N_tupelvars]={
0.,//  "Pt",
-5.17,//  "Eta",
-3.15,//  "Phi",
0.,//  "E",
0.,//  "Id",
0.,//  "RawPt",
0.,//  "RawE",
0.,//  "HfHadE",
0.,//  "HfEmE",
0.,//  "ChHadFrac",
0.,//  "NeutralHadAnHfFrac",
0.,//  "ChEmFrac",
0.,//  "NeutralEmFrac",
0.,//  "ChMult",
0.,//  "ConstCnt"  
0.// "NeuMult"
};
const float tupel_xhi[N_tupelvars]={
2000.,//  "Pt",
5.17,//  "Eta",
3.15,//  "Phi",
2000.,//  "E",
4.,//  "Id",
2000.,//  "RawPt",
2000.,//  "RawE",
2000.,//  "HfHadE",
2000.,//  "HfEmE",
2.,//  "ChHadFrac",
2.,//  "NeutralHadAnHfFrac",
2.,//  "ChEmFrac",
2.,//  "NeutralEmFrac",
100.,//  "ChMult",
100.,//  "ConstCnt"    
100.//"NeuMult"
};
  

const int tupel_jtInd=0, tupel_jetID=7;


//L1
const std::string L1BitStrings[]={//this array is a good idea
  "L1_SingleJet28_BptxAND",
  "L1_SingleJet40_BptxAND",
  "L1_SingleJet48_BptxAND",
  "L1_SingleJet52_BptxAND"
};
//const int N_L1Bits=N_HLTBits;
const int N_L1Bits=sizeof(L1BitStrings)/sizeof(std::string);

// HLT, both Calo and PF available
const double HLTCalothresh[]={
  45., 65., 85., 110.   // thought HLT80 threshold was too strict for Calo, now looser
};
const std::string Calo_HLTBitStrings[]={
  "HLT_AK4CaloJet40_Eta5p1",
  "HLT_AK4CaloJet60_Eta5p1",
  "HLT_AK4CaloJet80_Eta5p1",
  "HLT_AK4CaloJet100_Eta5p1" 
};//const int N_HLTBits=sizeof(Calo_HLTBitStrings)/sizeof(std::string);
const double HLTPFthresh[]={
  //55., 75., 95., 135.  //95 GeV as min/max for HLT80/60 too low
  //55., 75., 105., 135.   //55/75/105 is john's suggestion; i feel that 75 as min/max for HLT60/40 might be slightly too low. 
  //55., 75., 105., 135.   //TODOAGAIN/JOHNS/W |Y|<3.0 FOR INCL SPECTRA COMPARISON 
  //55., 80., 105., 135.   //TODO/TRY 1 -- DONE
  //55., 85., 105., 135.   //TODO/TRY 2 -- DONE
  //55., 90., 105., 135.   //TODO/TRY 3 -- DONE
  55., 90., 110., 135.   //TODO/TRY 4 -- DONE
  //60., 90., 110., 135.   //TODO/TRY 4 -- NOT DONE
};
const std::string PF_HLTBitStrings[]={
  "HLT_AK4PFJet40_Eta5p1",
  "HLT_AK4PFJet60_Eta5p1",
  "HLT_AK4PFJet80_Eta5p1",
  "HLT_AK4PFJet100_Eta5p1" 
};
const int N_HLTBits=sizeof(PF_HLTBitStrings)/sizeof(std::string);

// data tree name array
const std::string dataTreeNames[]={
  "GARBAGE ENTRY",
  "hiEvtAnalyzer/HiTree", 
  "skimanalysis/HltTree",   
  "ppTrack/trackTree",
  "hltanalysis/HltTree",
  "tupel/EventTree",  
  "hltobject/", 
  "hltobject/", 
  "hltobject/", 
  "hltobject/"
};
const int N_dataTrees=sizeof(dataTreeNames)/sizeof(std::string);


// MC tree name array (NOT FINISHED, WILL ERROR IN CURRENT FORM)... will it though?
const std::string MCTreeNames[]={
  "GARBAGE ENTRY",
  "hiEvtAnalyzer/HiTree", 
  "skimanalysis/HltTree",
  "ppTrack/trackTree",
  "tupel/EventTree"
  // "hltanalysis/HltTree",
  // "hltobject/"+Calo_HLTBitStrings[0]+"_v" ,
  // "hltobject/"+Calo_HLTBitStrings[1]+"_v" ,
  // "hltobject/"+Calo_HLTBitStrings[2]+"_v" ,
  // "hltobject/"+Calo_HLTBitStrings[3]+"_v" 
};
const int N_MCTrees=sizeof(MCTreeNames)/sizeof(std::string);

// for 5 tev pp jets 2k17
const float absetabins[]={
  //  0.0, 0.5, 1.0, 1.5, 2.0, 2.4, 2.5, 2.7, 3.0, 3.2, 3.7, 4.2, 4.7, 5.1 //alt binning i've seen before
  0.0, 
  0.5,
  1.0, 
  1.5,
  2.0//, 
  //2.5,
  //3.0//, 
  //  3.2,
  //  4.7
};
const int nbins_abseta=sizeof(absetabins)/sizeof(float)-1;//this is the # of bins in a root histo
//const int nbins_abseta=sizeof(absetabins)/sizeof(float); //this is the # of elements in the array
std::string absetabins_str[]={  
  "0.0", 
  "0.5",
  "1.0", 
  "1.5",
  "2.0"//, 
  //"2.5",
  //"3.0"//
  //"3.2",
  //"4.7"
};


const float etabins[]={
  -5.1 ,
  -4.7 , -4.4 , -4.1 , -3.8 , -3.5 ,
  -3.2 , -3.0 , -2.7 , -2.4 , -2.1 ,
  -1.8 , -1.5 , -1.2 ,
  -0.9 , -0.6 , -0.3 ,
  0.0 ,
  0.3 ,  0.6 ,  0.9 ,
  1.2 ,  1.5 ,  1.8 ,
  2.1 ,  2.4 ,  2.7 ,  3.0 ,  3.2 ,
  3.5 ,  3.8 ,  4.1 ,  4.4 ,  4.7 ,
  5.1
};
const int nbins_eta=sizeof(etabins)/sizeof(float)-1;


////straight up SMP binning for |y| bins in steps of 0.5
const float ptbins[]={
  //  //1., 
  //  //5., 
  //  //6., 
  //  //8., 
  //  //10., 
  //  //12., 
  //  //15., 
  //  18., 
  //  21., 
  //  24., 
  28.,
  32., 
  37.,
  43., //garbage bins
  49., //typical gen pt cut
  56., //typical recopt cut
  64.,
  74.,
  84.,
  97.,
  114.,
  133.,
  153.,
  174.,
  196.,
  220.,
  245.,
  272.,
  300.,
  330.,
  362.,
  395.,
  430.,
  468.,
  507.,
  548.,
  592.,
  638.,
  686.,
  737.,
  790.,
  846.,
  905.,
  967.,//,
  //1000.,
  1032., //generally, garbage bins
  1101., 
  1172., 
  1248., 
  1327.,
  1410.//,
  //1497., 1588., 1684., 1784., 1890., 2000., 2116.
};
const int nbins_pt=sizeof(ptbins)/sizeof(float)-1;//above values define edges of bins, not centers, so subtract one





//////John's pt binning for |y| bins in steps of 0.5
//const float ptbins[]={
//  56.,  
//  64.,  
//  84., 
//  114., 
//  153., 
//  196., 
//  245.,
//  300.,
//  362.,
//  430., 
//  507., 
//  592., 
//  686., 
//  790., 
//  905.,
//  1032.,
//  1101.
//};
//const int nbins_pt=sizeof(ptbins)/sizeof(float)-1;//above values define edges of bins, not centers, so subtract one
//
//
//

////Chris's pt binning for |y|<2.0
//const float ptbins[]={
//  20.,
//  30.,
//  40.,
//  50.,
//  60.,
//  70.,
//  80.,
//  90.,
//  100.,
//  110.,
//  120.,
//  130.,
//  140.,
//  150.,
//  160.,
//  170.,
//  180.,
//  190.,
//  200.,
//  210.,
//  220.,
//  230.,
//  240.,
//  250.,
//  260.,
//  270.,
//  280.,
//  290.,
//  310.,
//  330.,
//  350.,
//  370.,
//  400.,
//  430.,
//  500.,
//  1000.,
//  1500.
//
//};
//const int nbins_pt=sizeof(ptbins)/sizeof(float)-1;//above values define edges of bins, not centers, so subtract one
































// // my APSDNP Binning
//  32.,   37., 43., 49., 56., 64., 74., 84., 97., 114., 133., 153., 174., 196., 220., 245., 272., 300., 330., 362., 395., 430., 468.,
//  507., 548., 592., 638., 686., 737., 790., 846., 905., 967., //def junk bins after this or so //  1000., 1050. 

//binning used by john for multiplicity stuff
const float ptbins2[]={
  56,80,
  100,150,200,300,400,
  500,800,1000
};
const int nbins_pt2=sizeof(ptbins2)/sizeof(float)-1;//above values define edges of bins, not centers, so subtract one
























//pthat bins and weights for semi-private MC in millibarns
const float pthatbins[]={ 15., 30., 50., 80., 120., 170., 220., 280., 370., 460., 540., 9999. };
const int nbins_pthat=sizeof(pthatbins)/sizeof(int)-1;
const double pthatWeights[]={
  /*(pthat>=15)&&(pthat<30)    , n[0]=932778   , xsDiff=0.49235 mb    , weight= */ 5.27832e-07,
  /*(pthat>=30)&&(pthat<50)    , n[1]=903567   , xsDiff=0.030482 mb   , weight= */ 3.37352e-08,
  /*(pthat>=50)&&(pthat<80)    , n[2]=983531   , xsDiff=0.0035721 mb  , weight= */ 3.63191e-09,
  /*(pthat>=80)&&(pthat<120)   , n[3]=1820782  , xsDiff=0.00042494 mb , weight= */ 2.33383e-10,
  /*(pthat>=120)&&(pthat<170)  , n[4]=1080554  , xsDiff=5.873e-05 mb  , weight= */ 5.43517e-11,
  /*(pthat>=170)&&(pthat<220)  , n[5]=836152   , xsDiff=9.199e-06 mb  , weight= */ 1.10016e-11,
  /*(pthat>=220)&&(pthat<280)  , n[6]=954396   , xsDiff=2.2564e-06 mb , weight= */ 2.36422e-12,
  /*(pthat>=280)&&(pthat<370)  , n[7]=1083994  , xsDiff=6.336e-07 mb  , weight= */ 5.84505e-13,
  /*(pthat>=370)&&(pthat<460)  , n[8]=948240   , xsDiff=1.0884e-07 mb , weight= */ 1.14781e-13,
  /*(pthat>=460)&&(pthat<540)  , n[9]=1558268  , xsDiff=2.215e-08 mb  , weight= */ 1.42145e-14,
  /*(pthat>=540)&&(pthat<9999) , n[10]=2597338 , xsDiff=1.001e-08 mb  , weight= */ 3.85395e-15
};
////pthat bins and weights for official MC in millibarns
//const float pthatbins[]={ 15., 30., 50., 80., 120., 170., 220., 280., 370., 9999. }; 
//const int nbins_pthat=sizeof(pthatbins)/sizeof(int)-1;
//const float pthatWeights[]={
//  /*(pthat>=15)&&(pthat<30)     */   2.67718E-06,
//  /*(pthat>=30)&&(pthat<50)     */   1.59143E-07,
//  /*(pthat>=50)&&(pthat<80)     */   1.69053E-08,
//  /*(pthat>=80)&&(pthat<120)    */   1.95826E-09,
//  /*(pthat>=120)&&(pthat<170)   */   2.67247E-10,
//  /*(pthat>=170)&&(pthat<220)   */   4.29455E-11,
//  /*(pthat>=220)&&(pthat<280)   */   9.59090E-12,
//  /*(pthat>=280)&&(pthat<370)   */   2.32869E-12,
//  /*(pthat>=370)&&(pthat<9999)  */   4.38436E-13 
//};


const bool jetIDDebug=false;
bool jetID_00eta24(float jetIDpt, 
		   float neSum, float phSum, float chSum, float eSum, float muSum,
		   int numConst, int chMult){
  bool funcDebug=jetIDDebug;
  if(funcDebug)
    std::cout<<"jetID_00eta24 called.";
  bool passesJetID=false;
  if( neSum/jetIDpt    < 0.99 &&
      phSum/jetIDpt    < 0.99 &&
      numConst         > 1    &&      
      muSum/jetIDpt    < 0.80 &&
      chSum/jetIDpt    > 0.   && 
      chMult           > 0    &&
      eSum/jetIDpt     < 0.99    ) passesJetID=true;	      
  if(funcDebug)
    std::cout<<" passesJetID="<<passesJetID<<std::endl;
  return passesJetID;
}

bool jetID_24eta27(float jetIDpt,
		   float neSum, float phSum, float muSum, 
		   int numConst){
  bool funcDebug=jetIDDebug;
  if(funcDebug)
    std::cout<<"jetID_24eta27 called.";
  bool passesJetID=false;
  if( neSum/jetIDpt    < 0.99 &&
      phSum/jetIDpt    < 0.99 &&
      numConst         > 1    &&
      muSum/jetIDpt    < 0.80    ) passesJetID=true;	      
  if(funcDebug)
    std::cout<<" passesJetID="<<passesJetID<<std::endl;
  return passesJetID;
}

bool jetID_27eta30(float jetIDpt,
		   float neSum, float phSum, int neuMult){
  //int numConst){ // int numConst, int neuMult){//
  bool funcDebug=jetIDDebug;
  if(funcDebug)
    std::cout<<"jetID_27eta30 called.";
  bool passesJetID=false;
  if(  phSum/jetIDpt  < .99 && 
       neSum/jetIDpt < .99 &&
       neuMult > 0 )passesJetID=true;   
  //  if(  true && 
  //       true &&
  //       numConst            > 0       ) passesJetID=true;   
  //diff versions         // CMSSW [76,80]X criterion
  //  if(  phSum/jetIDpt > 0.00 && 
  //       neSum/jetIDpt < 1.00 && 
  //       numConst            > 0       ) passesJetID=true; 
  //  if(  phSum/jetIDpt [< 0.90 ] / [ > 0.01 &&]	
  //       neSum/jetIDpt [null   ] / [ < 0.98 &&]	   
  //       neuMult            [> 2    ] / [ > 2      ] ) passesJetID=true;
  if(funcDebug)
    std::cout<<" passesJetID="<<passesJetID<<std::endl;
  return passesJetID;
}

bool jetID_32eta47(float jetIDpt, 
		   float phSum){// float phSum, float neSum, float trkSum, float neuMult){//
  bool funcDebug=jetIDDebug;
  if(funcDebug)
    std::cout<<"jetID_32eta47 called.";
  bool passesJetID=false;
  if( phSum < 0.4 &&
      true               &&
      true       ) passesJetID=true; 
  //diff version
//  if( phSum_F[jet]/jetIDpt > 0. &&                      // else if( phSum_F[jet]/jetIDpt < 0.90 &&
//      trkSum_F[jet] < 0.4       && 
//      neSum_F[jet]/jetIDpt > 0. &&                         //          neSum_F[jet]/jetIDpt < null &&
//      true       ) passesJetID=true;     //          neuMult            > 10  
  if(funcDebug)
    std::cout<<" passesJetID="<<passesJetID<<std::endl;
  return passesJetID;
}



//v2 jet ID is currently based on loose criteria of 13 TeV
bool jetID_00eta24_v2( float jetNeutralHadAndHfFrac ,
		       float jetNeutalEmFrac        ,
		       float jetConstCnt            ,
		       float jetChHadFrac           ,
		       float jetChMult              ,
		       float jetChEmFrac              ){
  bool funcDebug=jetIDDebug;
  if(funcDebug)
    std::cout<<"jetID_00eta24_v2 called.";
  bool passesJetID=false;
  if(  jetNeutralHadAndHfFrac   < 0.99 &&
       jetNeutalEmFrac          < 0.99 &&
       jetConstCnt              > 1.    &&      
       jetChHadFrac             > 0.   && 
       jetChMult                > 0.   &&
       jetChEmFrac              < 0.99    ) passesJetID=true;	      
  if(funcDebug)
    std::cout<<" passesJetID="<<passesJetID<<std::endl;
  return passesJetID;
}

bool jetID_24eta27_v2( float jetNeutralHadAndHfFrac ,	   
		       float jetNeutalEmFrac        ,		   
  		       float jetConstCnt            )          {
  bool funcDebug=jetIDDebug;
  if(funcDebug)
    std::cout<<"jetID_24eta27_v2 called.";
  bool passesJetID=false;
  if(   jetNeutralHadAndHfFrac    < 0.99 &&
        jetNeutalEmFrac           < 0.99 &&
        jetConstCnt               > 1.       ) passesJetID=true;	      
  if(funcDebug)
    std::cout<<" passesJetID="<<passesJetID<<std::endl;
  return passesJetID;
}

bool jetID_27eta30_v2 ( float jetNeutralEmFrac,
		        float jetNeuMult    ){  
  bool funcDebug=jetIDDebug;
  if(funcDebug)
    std::cout<<"jetID_27eta30 called.";
  bool passesJetID=false;
  if(   jetNeutralEmFrac   < 0.90 &&
	jetNeuMult   > 2.)passesJetID=true;
  if(funcDebug)
    std::cout<<" passesJetID="<<passesJetID<<std::endl;
  return passesJetID;
}

bool jetID_32eta47_v2(float jetNeutralEmFrac,
		      float jetNeuMult    ){
  bool funcDebug=jetIDDebug;
  if(funcDebug)
    std::cout<<"jetID_32eta47 called.";
  bool passesJetID=false;
  if( jetNeutralEmFrac <0.9 &&
      jetNeuMult >10. )passesJetID=true;
  if(funcDebug)
    std::cout<<" passesJetID="<<passesJetID<<std::endl;
  return passesJetID;
}



//void makeJetQAHists(TH1D** hJetQA, int jetIDint, bool fillDijetHists, int N_vars){
//  
//
//  for(int j = 0; j<N_vars; ++j){    
//    //jets
//    if(var[j]=="jtpt"||var[j]=="rawpt")
//      hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;", var[j].c_str()) , 2500,0,2500);
//    else if(var[j]=="jteta")
//      hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;",var[j].c_str()) , 100,-5,+5);
//    else if(var[j]=="jtphi")
//      hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;",var[j].c_str()) , 100,-4,+4);
//    //jetconst. counts
//    else if(var[j]=="trkN"|| var[j]=="phN"|| var[j]=="chN"|| var[j]=="neN"|| var[j]=="eN"|| var[j]=="muN")
//      hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 60,0,60);
//    else if(var[j]=="trkHardN"|| var[j]=="phHardN"|| var[j]=="chHardN")
//      hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 60,0,60);
//    else if(var[j]=="neuMult"|| var[j]=="chMult"|| var[j]=="numConst")
//      hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 100,0,100);
//    //dijets
//    else if (fillDijetHists){
//      if(var[j]=="dphi")
//	hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;",var[j].c_str()) , 50,0,4);
//      else if(var[j]=="leadJetPt"||var[j]=="subleadJetPt")
//	hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;", var[j].c_str()), 2500,0,2500);
//      else //xj and Aj binnings
//	hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 200,0,2);
//    }
//    //consituent binnings
//    else
//      hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 200,0,2);
//  }
//  return;
//}












const bool doVzWeights=true;


double cpuVzWeight_poly(float vz_F){
  bool funcDebug=false;
  if(funcDebug)
    std::cout<<"vz_F="<<vz_F<<std::endl;
  // new set using hVz w/ weight=1 and no trig req. / pthat weighted MC // older weights
  double p0 = 1.165;//1.159;
  double p1 = 0.009171;//0.008914;
  double p2 = -0.003257;//-0.003057;
  double p3 = -2.919e-05;//-2.287e-05;
  double p4 = -4.647e-06;//-5.287e-06;
  double p5 = 2.13e-08;//-1.582e-09;
  double p6 = 2.399e-08;//2.569e-08;
  double p7 = 4.267e-12;//2.782e-11;
  double p8 = -2.078e-11;//-2.097e-11;
  
  double vzWeight=0;
  vzWeight+=p0;
  vzWeight+=p1*pow(vz_F,1);
  vzWeight+=p2*pow(vz_F,2);
  vzWeight+=p3*pow(vz_F,3);
  vzWeight+=p4*pow(vz_F,4);
  vzWeight+=p5*pow(vz_F,5);
  vzWeight+=p6*pow(vz_F,6);
  vzWeight+=p7*pow(vz_F,7);
  vzWeight+=p8*pow(vz_F,8);  

  if(funcDebug)
    std::cout<<"vzWeight="<<vzWeight<<std::endl;

  return vzWeight;
}

//takes in a TH2, fills it with the contents of the TH1 according to the rules of an outerproduct
//binning of a axes of TH2 assumed symmetric. binning of TH1 assumed equivalent to binning of each axis of TH2
inline void fillCovMatrix(TH2D* covmat=NULL, TH1D* hist=NULL, int nbins=-1, double weight=0.){
  //bool funcDebug=true;
  for(int i=1; i<=nbins; i++){
    double val_i=hist->GetBinContent(i);
    if(!(val_i>0.))continue;
    double center_i=hist->GetBinCenter(i);
    for(int j=i; j<=nbins; j++) {
      if(i!=j) { //if(i!=j){
	double val_j=hist->GetBinContent(j);
	if( !(val_j>0) ) continue;
	double center_j=hist->GetBinCenter(j);
	covmat->Fill(center_i, center_j, weight*weight*val_i*val_j);
	covmat->Fill(center_j, center_i, weight*weight*val_i*val_j);      }      //covmat->SetBinContent(i,j,val_i*val_j);      }       
      else{// (i==j){	
	covmat->Fill(center_i, center_i, weight*weight*val_i*val_i);}      //covmat->SetBinContent(i,j, val_i*val_i);            }      
    }	
  }
  return;
}


// non-event triggered vz weights
const double vzWeights[]={           
0.0579229 ,         // vzLow=-24
0.0622996 ,
0.061942  ,
0.0902018 ,
0.0744494 ,
0.0785317 ,                // vzLow=-21.5
0.088784  ,
0.106222  ,
0.113585  ,
0.142753  ,
0.147792  ,               // vzLow=-19
0.155734  ,
0.192015  ,
0.216288  ,
0.245479  ,
0.283993  ,               // vzLow=-16.5
0.31437	  ,
0.352618  ,
0.390549  ,
0.43129	  ,
0.468473  ,               // vzLow=-14
0.501128  ,
0.539377  ,
0.582421  ,
0.647345  ,
0.643217  ,               // vzLow=-11.5
0.707504  ,
0.73216	  ,
0.768601  ,
0.819635  ,
0.836843  ,               // vzLow=-9
0.874235  ,
0.89897	  ,
0.917435  ,
0.948561  ,
0.983956  ,               // vzLow=-6.5
1.02518	  ,
1.03036	  ,
1.039	  ,
1.06393	  ,
1.0824 	  ,            // vzLow=-4
1.09788	  ,
1.1064	  ,
1.12646	  ,
1.13744	  ,
1.14794   ,              // vzLow=-1.5
1.14349	  ,
1.16373	  ,
1.16068	  ,
1.16672	  ,
1.16045   ,              // vzLow=1
1.17421	  ,
1.15623	  ,
1.1739	  ,
1.14798	  ,
1.14334   ,             // vzLow=3.5
1.13985	  ,
1.12823	  ,
1.11934	  ,
1.11583	  ,
1.09352   ,              // vzLow=6
1.05137	  ,
1.04333	  ,
1.01349	  ,
0.986353  ,
0.959904  ,              // vzLow=8.5
0.929557  ,
0.913313  ,
0.88634	  ,
0.827843  ,
0.807459  ,              // vzLow=11
0.743207  ,
0.727486  ,
0.677889  ,
0.644392  ,
0.594507  ,               // vzLow=13.5
0.56391	  ,
0.491046  ,
0.465995  ,
0.429629  ,
0.390355  ,               // vzLow=16
0.337874  ,
0.31006	  ,
0.274059  ,
0.238072  ,
0.217408  ,               // vzLow=18.5
0.189755  ,
0.169915  ,
0.164071  ,
0.137212  ,
0.136876  ,               // vzLow=21
0.0945742 ,
0.101419  ,
0.0768486 ,
0.075515  ,
0.0863076                 // vzLow=23.5
};
const int nbins_vzWeights=sizeof(vzWeights)/sizeof(float);//should be 96
const float minbinValue_vzWeights=-24., maxbinValue_vzWeights=24.;
const float binsize_vzWeights=(maxbinValue_vzWeights-minbinValue_vzWeights)/nbins_vzWeights;//should be 0.5






const UInt_t runarray[]={
  262163,
  262165,
  262167,
  262173,
  262174,
  262204,
  262205,
  262235,
  262248,
  262250,
  262252,
  262253,
  262254,
  262270,
  262271,
  262272,
  262273,
  262274,
  262275,
  262277,
  262325,
  262326,
  262327,
  262328,
  999999
};
const int Nruns=sizeof(runarray)/sizeof(UInt_t);



//use me for Jet80/LowerJets/MinBias6 Study
const UInt_t runarray_forExclRunStudy[]={
  262173,
  262204,
  262205,
  262248,
  262272,
  262273,
  262275,
  262277,
  262325,
  262326,
  262327,
  262328
};
const int NGoodRuns=sizeof(runarray_forExclRunStudy)/sizeof(UInt_t);

////use me for Jet80/LowerJets Study
//const UInt_t runarray_forExclRunStudy[]={
//  262163,
//  262165,
//  262167,
//  262173,
//  262204,
//  262205,
//  262248,
//  262250,
//  262252,
//  262253,
//  262254,
//  262270,
//  262271,
//  262272,
//  262273,
//  262275,
//  262277,
//  262325,
//  262326,
//  262327,
//  262328
//};
//const int NGoodRuns=sizeof(runarray_forExclRunStudy)/sizeof(UInt_t);
