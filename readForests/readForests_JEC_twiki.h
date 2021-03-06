//these JEC below from 
//// https://twiki.cern.ch/twiki/pub/CMS/JECHI2015/jec_pp502_ppreco.tar.gz
const float JEC_etabins[]={
    -3, 
  -2.5, 
-2.043, 
 -1.74, 
-1.392, 
-1.131, 
-0.879, 
-0.609, 
-0.435, 
-0.261, 
-0.087, 
     0, 
 0.087, 
 0.261, 
 0.435, 
 0.609, 
 0.879, 
 1.131, 
 1.392, 
  1.74, 
 2.043, 
   2.5, 
   3.0
};
const int n_JEC_etabins=sizeof(JEC_etabins)/sizeof(float)-1;

std::string L2RelativeJEC_str="([0]+[1]/((log10(x)^2)+[2])+[3]*exp(-[4]*(log10(x)-[5])*(log10(x)-[5])))";

const float L2RelJEC_par0_arr[]=  {
0.152933,  
0.622836,  
0.623008,  
0.517384,  
0.550841,  
0.505359,  
0.489862,  
0.488863,  
0.489742,  
0.491527,  
0.505095,  
0.505014,  
0.485949,  
0.488406,  
0.496361,  
0.484733,  
0.511506,  
0.550964,  
0.524713,  
0.629732,  
0.615143,  
0.205773   
};
const float L2RelJEC_par1_arr[]=  {
  0.324127,   
  -0.315834,   
  -0.375219,   
  0.230233,   
  -0.0852162,   
  0.114691,   
  0.154937,   
  0.161449,   
  0.163213,   
  0.147434,   
  0.0414896,   
  0.0308973,   
  0.171986,   
  0.159175,   
  0.0933105,   
  0.19155,   
  0.0853956,   
  -0.0365914,   
  0.188975,   
  -0.464362,   
  -0.27785,   
  0.209332   
};
const float L2RelJEC_par2_arr[]=  {
    2.3506, 
   1.19754, 
   2.30358, 
   2.86129, 
   3.51549, 
   2.83332, 
   2.64737, 
   2.64227, 
   2.64392, 
   2.67644, 
   2.92395, 
   2.92523, 
   2.59705, 
   2.63732, 
   2.78026, 
    2.5591, 
   2.93091, 
   3.47219, 
   2.98709, 
   2.44627, 
  0.919933, 
   2.97844 

};
const float L2RelJEC_par3_arr[]=  {
  0.893557, 
  0.635574, 
  0.600538, 
  0.491214, 
  0.522515, 
  0.477504, 
  0.462788, 
  0.461645, 
  0.462613, 
  0.463973, 
  0.473599, 
   0.47346, 
  0.459349, 
  0.461304, 
  0.466895, 
  0.458908, 
  0.482391, 
  0.521164, 
  0.497621, 
  0.602711, 
  0.641698, 
  0.874848 
};
const float L2RelJEC_par4_arr[]=  {
    0.181443, 
   0.0606863, 
   0.0483874, 
  0.00346462, 
   0.0111259, 
 -0.00416623, 
 -0.00895999, 
 -0.00905401, 
 -0.00869922, 
 -0.00888437, 
 -0.00490226, 
 -0.00456554, 
 -0.00975919, 
 -0.00865905, 
 -0.00669459, 
 -0.00995563, 
 -0.00292809, 
   0.0119433, 
  0.00461541, 
   0.0470842, 
   0.0608597, 
    0.157687 
};

const float L2RelJEC_par5_arr[]=  {
    1.43978,
  -0.154422,
  -0.209178,
 -0.0446651,
  -0.154945,
 -0.0100519,
  0.0191157,
  0.0210632,
    0.01893,
  0.0243322,
  0.0405298,
  0.0403433,
  0.0199142,
  0.0185436,
  0.0316187,
  0.0117113,
-0.00565626,
  -0.120012,
 -0.0528985,
  -0.169094,
  -0.240562,
    1.27834
};

const float L2RelJEC_ptlo_arr[]=  {
  19.5918    ,
  20.3998    ,
  20.5817    ,
  20.9155    ,
  21.0888    ,
  21.6636    ,
  21.9848    ,
  21.9868    ,
  22.0411    ,
  21.9639    ,
  22.1868    ,
  22.2029    ,
  22.1095    ,
  22.0217    ,
  22.2924    ,
  21.9559    ,
  21.5945    ,
  20.9972    ,
  20.9149    ,
  20.7592    ,
  20.5834    ,
  19.5957    ,
};

const float L2RelJEC_pthi_arr[]=  {
 300.886 ,
 396.168 ,
 546.199 ,
 751.421 ,
 777.297 ,
 793.744 ,
 806.436 ,
 807.991 ,
 809.351 ,
 808.827 ,
 811.224 ,
 813.948 ,
 812.464 ,
 811.933 ,
 811.096 ,
 807.369 ,
  794.06 ,
 779.359 ,
 761.455 ,
 534.256 ,
 403.165 ,
 303.718 

};



std::string L3AbsoluteJEC_str="[0]+[1]/(pow(log10(x)-[2],2)+[2])+[3]*exp((-[4]*(log10(x)-[5])*(log10(x)-[5]))+([6]*(log10(x)-[5])))";
float L3AbsJEC_ptlo=4., L3AbsJEC_pthi=1000.;
float L3Abs_par0=0.846919  ;
float L3Abs_par1=0.557384  ;
float L3Abs_par2=2.44612   ;
float L3Abs_par3=0.00160516;
float L3Abs_par4=3.57438   ;
float L3Abs_par5=2.23248   ;
float L3Abs_par6=-7.90008  ;  
