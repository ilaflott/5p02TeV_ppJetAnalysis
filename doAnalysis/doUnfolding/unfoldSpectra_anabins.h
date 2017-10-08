

//// ------ GEN SPECTRA BINS ------
//const double anabins_pt_gen[] = {//analysis pt bins
//  //  0, 
//  //  3, 4, 5, 7, 9, 12, 15, 18, 21, 
//  28, // 24, 
//  30, // 28, 
//  32, 
//  37, 
//  43, 49, 56, 64, 74, 84, 97, 114, 
//  133, 153, 174, 196, 220, 245, 272, 300, 
//  330, 362, 395, 430, 468, 507, 548, 592,
//  638, 686 , 
//  800, 1000
//};
//const int n_anabins_pt_gen = sizeof(anabins_pt_gen)/sizeof(double)-1;
//
//// ------ RECO SPECTRA BINS ------
//const double anabins_pt_reco[] = {//analysis pt bins
//  //  0, 
//  //  3, 4, 5, 7, 9, 12, 15, 18, 21, 
//  24, 
//  28,
//  32,
//  37, 
//  43, 49, 56, 64, 74, 84, 97, 114, 
//  133, 153, 174, 196, 220, 245, 272, 300, 
//  330, 362, 395, 430, 468, 507, 548, 592,
//  638, 686 , 
//  800, 1000
//};
//const int n_anabins_pt_reco = sizeof(anabins_pt_reco)/sizeof(double)-1;



// for readForests jobs w/ 30/30 rec/gen pt cuts
// same as old binning... shifted up by 2 to respect the cuts at readforests-level
// ------ GEN SPECTRA BINS ------
const double anabins_pt_gen[] = {//analysis pt bins
  //  32.,
  37.,
  43., 
  49.,
  56.,
  64.,
  74.,
  84., 97., 114., 133., 153., 174., 196., 220., 245., 272., 300., 330., 362., 395., 430., 468.,
  507., 548., 592., 638.,
  686., 1000.,
  1300.
  //  1500.//  , //junk bin after this
  //1050.
};
const int n_anabins_pt_gen = sizeof(anabins_pt_gen)/sizeof(double)-1;

// ------ RECO SPECTRA BINS ------
const double anabins_pt_reco[] = {//analysis pt bins
  //  32.,
  //  37.,
  //  43., //junk bins above
  49.,
  56.,
  64.,
  74.,
  84., 97., 114., 133., 153., 174., 196., 220., 245., 272., 300., 330., 362., 395., 430., 468.,
  507., 548., 592., 638.,
  686., 1000.
  //  , //junk bin after this
  //1050.
};
const int n_anabins_pt_reco = sizeof(anabins_pt_reco)/sizeof(double)-1;






//// for readForests jobs w/ 40/40 rec/gen pt cuts
//// same as old binning... shifted up by 2 to respect the cuts at readforests-level
//// ------ GEN SPECTRA BINS ------
//const double anabins_pt_gen[] = {//analysis pt bins
//  //  24, 
//  //30,  34,  
//  40,   46, 
//  52, 59, 67, 77, 87, 
//  100, 117, 136, 156, 
//  177, 199, 223, 248, 
//  275, 303, 333, 365, 
//  398, 433, 471, 510, 
//  551, 595, 641, 689 , 
//  739, 811, 900
//  //802,     1000
//};
//const int n_anabins_pt_gen = sizeof(anabins_pt_gen)/sizeof(double)-1;
//
//// ------ RECO SPECTRA BINS ------
//const double anabins_pt_reco[] = {//analysis pt bins
//  //  24, 
//  //30,  34,  
//  40,   46, 
//  52, 59, 67, 77, 87, 
//  100, 117, 136, 156, 
//  177, 199, 223, 248, 
//  275, 303, 333, 365, 
//  398, 433, 471, 510, 
//  551, 595, 641, 689 , 
//  739, 811, 900
//  //802,     1000
//};
//const int n_anabins_pt_reco = sizeof(anabins_pt_reco)/sizeof(double)-1;

