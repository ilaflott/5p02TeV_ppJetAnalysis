//// ------ SMP SPECTRA BINS ORIGINAL ------
//// ------ DO NOT TOUCH ------
//const double anabins_pt_SMP[] = {
//  1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 
//  32, 37, 43, 49, 
//  56, 64, 74, 84,97, 114, 
//  133, 153, 174, 196, 220, 
//  245, 272, 300,   330, 
//  362, 395, 430, 
//  468,  507, 
//  548, 592, 638, 
//  686, 
//  737, 790, 846,  905,   967,  1032, 
//  1101, 1172, 1248, 1327, 
//  //1410, 1497, 1588, 1684, 1784, 1890, 2000,
//  //2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 
//  //4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000
//  };





//// ------ CHRIS MCGINNS BINNING ------ //
//// ------ GEN SPECTRA BINS, IN PROGRESS------
//const double anabins_pt_gen[] = {//analysis pt bins
//  60.,70.,80.,90.,100.,110.,
//  120.,130.,140.,150.,160.,170.,
//  180.,190.,200.,210.,220.,230.,
//  240.,250.,260.,270.,280.,290.,
//  310.,330.,350.,370.,400.,430.,
//  500.//,//1500.,
//};
//const int n_anabins_pt_gen = sizeof(anabins_pt_gen)/sizeof(double)-1;
//
//// ------ RECO SPECTRA BINS IN PROGRESS------
//const double anabins_pt_reco[] = {//analysis pt bins
//  60.,70.,80.,90.,100.,110.,
//  120.,130.,140.,150.,160.,170.,
//  180.,190.,200.,210.,220.,230.,
//  240.,250.,260.,270.,280.,290.,
//  310.,330.,350.,370.,400.,430.,
//  500.//,//1500.,
//};
//const int n_anabins_pt_reco = sizeof(anabins_pt_reco)/sizeof(double)-1;





//// ------ GEN SPECTRA BINS, IN PROGRESS------
const double anabins_pt_gen[] = {//analysis pt bins
  //1., 5., 6., 8., 10., 12., 15., 
  //18., //21., //24., 
  //28., //32., //37., 
  //43., //garbage bins
  //49., //typical gen pt cut
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
  686.//, 
  //737., 
  //790., 
  //846.,  
  //905.,
  //1000.////967.//,  
  //1032.//, 
  //1101.//, 
  //1172.//, 
  //1248., 
  //1327.,
  //1410.//, 
  //1497.//,//1588.//, 1684. //, 1784., 1890., 2000., 2116.//,2238.  
};
const int n_anabins_pt_gen = sizeof(anabins_pt_gen)/sizeof(double)-1;

// ------ RECO SPECTRA BINS IN PROGRESS------
const double anabins_pt_reco[] = {//analysis pt bins
  //1., 5., 6., 8., 10., 12., 15., 
  //18., //21., //24., 28., //32., //37., 
  //43., //garbage bins
  //49., //typical gen pt cut
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
  686.//, 
  //  737., 
  //790., 
  //846.,  
  //905.,   
  //967.//,  
  //1000.//
  //1032.//,
  //1101.//, 
  //1172.//, 
  //1248., 
//  1327.,
//  1410.//, 
  //1497.//, 1588. //, 1684. //, 1784., 1890., 2000., 2116.//,2238.  
};
const int n_anabins_pt_reco = sizeof(anabins_pt_reco)/sizeof(double)-1;




































/*
//// These analysis binnings (based on SMP binning) are excellent for Bayesian Unfolding of the following samples + settings
source runData.sh "0.0eta2.0_ptHi967_ptLoRec49"  "0.0eta2.0_ptHi967_ptLoRec49_ptLoGen49" "05" "09" "18" 1 1 0

const bool doToyErrs          =false; // error options: kCovToy, kCovariance, kErrors, kNoError
const bool doMCIntegralScaling=false;
const bool fillRespHists      =true;
//const bool useTH2ProjRespHist= false;
const bool doOverUnderflows   =false;//leave false almost always
const bool clearOverUnderflows=true; 

// settings that don't really get used
const bool zeroBins=false; //leave false almost always
//const bool normalizedMCMatrix=false;   

// ------ GEN SPECTRA BINS, IN PROGRESS------
const double anabins_pt_gen[] = {//analysis pt bins
  //1., 
  //5., 6., 8., 10., 12., 15., 
  //  18., 
  //  21., 
  //  24., 
  //  28., 
  //  32., 
  //  37., 
  //  43., //garbage bins
  //49., //typical gen pt cut
  //56., //typical recopt cut
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
  //395., 
  430., 
  //468., 
  //507., 
  //548., 
  592., 
  //638., 
  //686., 
  //737., 
  //790., 
  //846.,  
  //905.,   
  967.//,
  //1032.//, 
  //  1101.//, 1172. , 
  //  
  //  1248., 
  //  1327.,
  //  1410., 
  //  1497.//, 1588. //, 1684. //, 1784., 1890., 2000., 2116.//,2238.  
};
const int n_anabins_pt_gen = sizeof(anabins_pt_gen)/sizeof(double)-1;

// ------ RECO SPECTRA BINS IN PROGRESS------
const double anabins_pt_reco[] = {//analysis pt bins
  //1., 
  //5., 6., 8., 10., 12., 15., 
  //  18., 
  //  21., 
  //  24., 28.,
  //  32., 
  //  37., 
  //  43., //garbage bins
  //  49., //typical gen pt cut
  //56., //typical recopt cut
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
  //395., 
  430., 
  //468., 
  //507., 
  //548., 
  592., 
  //638., 
  //686., 
  //737., 
  //790., 
  //846.,  
  //905.,   
  967.//,
  //  1032.//, 
  //1032.//, 
  //1101.//, 1172. , 
  //1248., 
  //1327.,
  //1410., 
  //1497.//, 1588. //, 1684. //, 1784., 1890., 2000., 2116.//,2238.  
};
const int n_anabins_pt_reco = sizeof(anabins_pt_reco)/sizeof(double)-1;





 */
