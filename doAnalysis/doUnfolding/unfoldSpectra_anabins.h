






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



//same as old binning... shifted up by 2 to respect the cuts at readforests-level
// ------ GEN SPECTRA BINS ------
const double anabins_pt_gen[] = {//analysis pt bins
//  24, 
  30,  34,  39,   45, 
  51, 58, 66, 76, 86, 
  99, 116, 135, 155, 
  176, 198, 222, 247, 
  274, 302, 332, 364, 
  397, 432, 470, 509, 
  550, 594, 640, 688 , 
  738, 810, 900
  //802,     1000
};
const int n_anabins_pt_gen = sizeof(anabins_pt_gen)/sizeof(double)-1;

// ------ RECO SPECTRA BINS ------
const double anabins_pt_reco[] = {//analysis pt bins
//  24, 
  30,  34,  39,   45, 
  51, 58, 66, 76, 86, 
  99, 116, 135, 155, 
  176, 198, 222, 247, 
  274, 302, 332, 364, 
  397, 432, 470, 509, 
  550, 594, 640, 688 , 
  738, 810, 900
  //802,     1000
};
const int n_anabins_pt_reco = sizeof(anabins_pt_reco)/sizeof(double)-1;

