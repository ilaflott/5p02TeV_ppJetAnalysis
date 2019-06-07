float jecunc_etabins[]={-5.4,-5,-4.4,-4,-3.5,-3,-2.8,-2.6,-2.4,-2.2,-2,-1.8,-1.6,-1.4,-1.2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.5,4,4.4,5, 5.4 };//final 5.4 is an extra eta bin boundary, makes 5 the final lower limit for finding the bin
const int n_jecunc_etabins=sizeof(jecunc_etabins)/sizeof(float)-1;
float jecunc_ptbins[]={9,11,13.5,16.5,19.5,22.5,26,30,34.5,40,46,52.5,60,69,79,90.5,105.5,123.5,143,163.5,185,208,232.5,258.5,286,331,396,468.5,549.5,639,738,847.5,968.5,1102,1249.5,1412,9999};//final 9999 is an extra pt bin boundary, make 1412 the final lower limit for finding the bin
const int n_jecunc_ptbins=sizeof(jecunc_ptbins)/sizeof(float)-1;




double jecunc[n_jecunc_etabins][n_jecunc_ptbins]={
  {0.0344,0.0332,0.0322,0.0313,0.0306,0.0301,0.0297,0.0292,0.0280,0.0270,0.0260,0.0253,0.0241,0.0228,0.0220,0.0216,0.0211,0.0208,0.0205,0.0203,0.0201,0.0201,0.0200,0.0200,0.0200,0.0200,0.0200,0.0200,0.0202,0.0203,0.0204,0.0206,0.0207,0.0209,0.0211,0.0213},
  {0.0344,0.0332,0.0322,0.0313,0.0306,0.0301,0.0297,0.0292,0.0280,0.0269,0.0260,0.0253,0.0241,0.0228,0.0216,0.0207,0.0200,0.0197,0.0194,0.0191,0.0190,0.0189,0.0188,0.0187,0.0187,0.0187,0.0187,0.0188,0.0189,0.0191,0.0192,0.0194,0.0196,0.0197,0.0199,0.0201},
  {0.0344,0.0332,0.0322,0.0313,0.0306,0.0301,0.0297,0.0292,0.0280,0.0269,0.0260,0.0253,0.0241,0.0227,0.0215,0.0206,0.0199,0.0197,0.0198,0.0203,0.0210,0.0213,0.0213,0.0212,0.0213,0.0213,0.0213,0.0214,0.0215,0.0216,0.0217,0.0219,0.0220,0.0222,0.0224,0.0225},
  {0.0344,0.0332,0.0321,0.0313,0.0306,0.0301,0.0296,0.0292,0.0279,0.0268,0.0259,0.0252,0.0241,0.0227,0.0215,0.0206,0.0199,0.0196,0.0197,0.0203,0.0210,0.0219,0.0230,0.0241,0.0253,0.0261,0.0261,0.0262,0.0264,0.0265,0.0266,0.0267,0.0268,0.0269,0.0271,0.0272},
  {0.0344,0.0332,0.0321,0.0312,0.0305,0.0301,0.0296,0.0292,0.0279,0.0268,0.0259,0.0252,0.0241,0.0226,0.0214,0.0205,0.0198,0.0196,0.0197,0.0202,0.0210,0.0219,0.0229,0.0241,0.0252,0.0272,0.0297,0.0323,0.0336,0.0336,0.0337,0.0338,0.0339,0.0340,0.0341,0.0343},
  {0.0418,0.0408,0.0399,0.0392,0.0387,0.0383,0.0379,0.0375,0.0366,0.0357,0.0351,0.0346,0.0340,0.0319,0.0302,0.0288,0.0278,0.0273,0.0275,0.0281,0.0291,0.0304,0.0320,0.0338,0.0355,0.0384,0.0422,0.0460,0.0501,0.0524,0.0524,0.0525,0.0526,0.0526,0.0528,0.0529},
  {0.0300,0.0286,0.0273,0.0262,0.0254,0.0248,0.0242,0.0237,0.0221,0.0207,0.0195,0.0186,0.0176,0.0162,0.0150,0.0137,0.0127,0.0117,0.0110,0.0107,0.0106,0.0107,0.0110,0.0114,0.0119,0.0128,0.0141,0.0155,0.0169,0.0182,0.0196,0.0198,0.0200,0.0202,0.0204,0.0206},
  {0.0285,0.0270,0.0257,0.0245,0.0236,0.0230,0.0224,0.0218,0.0201,0.0185,0.0171,0.0161,0.0150,0.0140,0.0130,0.0120,0.0111,0.0102,0.0094,0.0089,0.0085,0.0082,0.0080,0.0078,0.0079,0.0080,0.0082,0.0086,0.0091,0.0096,0.0101,0.0106,0.0111,0.0114,0.0118,0.0121},
  {0.0280,0.0265,0.0252,0.0240,0.0231,0.0225,0.0218,0.0212,0.0195,0.0178,0.0164,0.0153,0.0143,0.0133,0.0123,0.0114,0.0104,0.0095,0.0087,0.0081,0.0077,0.0074,0.0071,0.0070,0.0069,0.0070,0.0072,0.0074,0.0079,0.0083,0.0089,0.0093,0.0098,0.0103,0.0107,0.0110},
  {0.0283,0.0269,0.0255,0.0244,0.0235,0.0229,0.0222,0.0216,0.0199,0.0183,0.0170,0.0159,0.0149,0.0138,0.0128,0.0118,0.0108,0.0098,0.0089,0.0083,0.0078,0.0074,0.0071,0.0069,0.0069,0.0069,0.0070,0.0073,0.0077,0.0082,0.0090,0.0098,0.0103,0.0109,0.0114,0.0117},
  {0.0281,0.0266,0.0252,0.0241,0.0232,0.0225,0.0219,0.0213,0.0195,0.0179,0.0165,0.0154,0.0144,0.0134,0.0124,0.0114,0.0105,0.0096,0.0088,0.0082,0.0078,0.0075,0.0072,0.0070,0.0070,0.0070,0.0071,0.0073,0.0077,0.0081,0.0086,0.0090,0.0095,0.0100,0.0105,0.0109},
  {0.0281,0.0266,0.0252,0.0241,0.0232,0.0225,0.0218,0.0212,0.0195,0.0179,0.0165,0.0153,0.0144,0.0133,0.0124,0.0114,0.0104,0.0096,0.0088,0.0082,0.0077,0.0074,0.0072,0.0070,0.0069,0.0069,0.0071,0.0073,0.0077,0.0080,0.0086,0.0090,0.0094,0.0100,0.0105,0.0109},
  {0.0277,0.0262,0.0249,0.0237,0.0228,0.0221,0.0214,0.0208,0.0190,0.0174,0.0159,0.0148,0.0138,0.0128,0.0119,0.0110,0.0100,0.0092,0.0085,0.0079,0.0075,0.0072,0.0069,0.0068,0.0068,0.0068,0.0069,0.0071,0.0075,0.0078,0.0083,0.0087,0.0091,0.0096,0.0101,0.0105},
  {0.0275,0.0260,0.0246,0.0234,0.0225,0.0218,0.0211,0.0205,0.0187,0.0170,0.0155,0.0143,0.0133,0.0123,0.0114,0.0105,0.0096,0.0087,0.0080,0.0073,0.0069,0.0066,0.0063,0.0061,0.0061,0.0061,0.0062,0.0064,0.0068,0.0071,0.0076,0.0080,0.0084,0.0088,0.0092,0.0097},
  {0.0275,0.0260,0.0246,0.0234,0.0225,0.0218,0.0211,0.0205,0.0187,0.0170,0.0155,0.0143,0.0133,0.0123,0.0114,0.0105,0.0096,0.0087,0.0080,0.0073,0.0069,0.0066,0.0063,0.0061,0.0061,0.0061,0.0062,0.0064,0.0068,0.0071,0.0076,0.0080,0.0084,0.0088,0.0092,0.0097},
  {0.0275,0.0259,0.0246,0.0234,0.0224,0.0218,0.0211,0.0205,0.0186,0.0169,0.0154,0.0142,0.0133,0.0123,0.0113,0.0104,0.0095,0.0086,0.0079,0.0073,0.0069,0.0066,0.0063,0.0061,0.0061,0.0061,0.0062,0.0064,0.0068,0.0071,0.0075,0.0079,0.0083,0.0087,0.0092,0.0096},
  {0.0275,0.0259,0.0245,0.0234,0.0224,0.0217,0.0211,0.0204,0.0186,0.0169,0.0154,0.0142,0.0132,0.0122,0.0113,0.0104,0.0095,0.0086,0.0079,0.0073,0.0069,0.0066,0.0063,0.0061,0.0061,0.0061,0.0062,0.0064,0.0068,0.0071,0.0075,0.0079,0.0083,0.0088,0.0092,0.0096},
  {0.0275,0.0260,0.0246,0.0234,0.0225,0.0218,0.0211,0.0205,0.0187,0.0170,0.0155,0.0143,0.0133,0.0123,0.0114,0.0105,0.0095,0.0086,0.0079,0.0073,0.0069,0.0066,0.0063,0.0061,0.0061,0.0061,0.0062,0.0064,0.0068,0.0071,0.0076,0.0080,0.0085,0.0089,0.0093,0.0098},
  {0.0275,0.0260,0.0246,0.0234,0.0225,0.0218,0.0211,0.0205,0.0187,0.0170,0.0155,0.0143,0.0133,0.0123,0.0114,0.0105,0.0095,0.0086,0.0079,0.0073,0.0069,0.0066,0.0063,0.0061,0.0061,0.0061,0.0062,0.0064,0.0068,0.0072,0.0076,0.0080,0.0085,0.0089,0.0093,0.0098},
  {0.0280,0.0265,0.0252,0.0240,0.0231,0.0224,0.0218,0.0212,0.0194,0.0178,0.0164,0.0153,0.0142,0.0132,0.0122,0.0112,0.0101,0.0092,0.0083,0.0077,0.0071,0.0068,0.0064,0.0062,0.0061,0.0061,0.0062,0.0065,0.0069,0.0073,0.0079,0.0083,0.0089,0.0094,0.0099,0.0105},
  {0.0280,0.0265,0.0252,0.0240,0.0231,0.0224,0.0218,0.0212,0.0194,0.0178,0.0164,0.0153,0.0142,0.0132,0.0122,0.0112,0.0101,0.0092,0.0083,0.0077,0.0071,0.0068,0.0064,0.0062,0.0061,0.0061,0.0062,0.0065,0.0069,0.0073,0.0079,0.0083,0.0089,0.0094,0.0099,0.0105},
  {0.0275,0.0260,0.0246,0.0234,0.0225,0.0218,0.0211,0.0205,0.0187,0.0170,0.0155,0.0143,0.0133,0.0123,0.0114,0.0105,0.0095,0.0086,0.0079,0.0073,0.0069,0.0066,0.0063,0.0061,0.0061,0.0061,0.0062,0.0064,0.0068,0.0072,0.0076,0.0080,0.0085,0.0089,0.0093,0.0098},
  {0.0275,0.0260,0.0246,0.0234,0.0225,0.0218,0.0211,0.0205,0.0187,0.0170,0.0155,0.0143,0.0133,0.0123,0.0114,0.0105,0.0095,0.0086,0.0079,0.0073,0.0069,0.0066,0.0063,0.0061,0.0061,0.0061,0.0062,0.0064,0.0068,0.0071,0.0076,0.0080,0.0085,0.0089,0.0093,0.0098},
  {0.0275,0.0259,0.0245,0.0234,0.0224,0.0217,0.0211,0.0204,0.0186,0.0169,0.0154,0.0142,0.0132,0.0122,0.0113,0.0104,0.0095,0.0086,0.0079,0.0073,0.0069,0.0066,0.0063,0.0061,0.0061,0.0061,0.0062,0.0064,0.0068,0.0071,0.0075,0.0079,0.0083,0.0088,0.0092,0.0096},
  {0.0275,0.0259,0.0246,0.0234,0.0224,0.0218,0.0211,0.0205,0.0186,0.0169,0.0154,0.0142,0.0133,0.0123,0.0113,0.0104,0.0095,0.0086,0.0079,0.0073,0.0069,0.0066,0.0063,0.0061,0.0061,0.0061,0.0062,0.0064,0.0068,0.0071,0.0075,0.0079,0.0083,0.0087,0.0092,0.0096},
  {0.0275,0.0260,0.0246,0.0234,0.0225,0.0218,0.0211,0.0205,0.0187,0.0170,0.0155,0.0143,0.0133,0.0123,0.0114,0.0105,0.0096,0.0087,0.0080,0.0073,0.0069,0.0066,0.0063,0.0061,0.0061,0.0061,0.0062,0.0064,0.0068,0.0071,0.0076,0.0080,0.0084,0.0088,0.0092,0.0097},
  {0.0275,0.0260,0.0246,0.0234,0.0225,0.0218,0.0211,0.0205,0.0187,0.0170,0.0155,0.0143,0.0133,0.0123,0.0114,0.0105,0.0096,0.0087,0.0080,0.0073,0.0069,0.0066,0.0063,0.0061,0.0061,0.0061,0.0062,0.0064,0.0068,0.0071,0.0076,0.0080,0.0084,0.0088,0.0092,0.0097},
  {0.0277,0.0262,0.0249,0.0237,0.0228,0.0221,0.0214,0.0208,0.0190,0.0174,0.0159,0.0148,0.0138,0.0128,0.0119,0.0110,0.0100,0.0092,0.0085,0.0079,0.0075,0.0072,0.0069,0.0068,0.0068,0.0068,0.0069,0.0071,0.0075,0.0079,0.0083,0.0088,0.0092,0.0096,0.0101,0.0105},
  {0.0281,0.0266,0.0252,0.0241,0.0232,0.0225,0.0218,0.0212,0.0195,0.0179,0.0165,0.0153,0.0143,0.0133,0.0124,0.0114,0.0104,0.0096,0.0088,0.0082,0.0077,0.0074,0.0072,0.0070,0.0070,0.0070,0.0071,0.0073,0.0077,0.0081,0.0086,0.0090,0.0095,0.0100,0.0105,0.0110},
  {0.0281,0.0266,0.0252,0.0241,0.0232,0.0225,0.0218,0.0212,0.0195,0.0179,0.0165,0.0153,0.0144,0.0133,0.0124,0.0114,0.0105,0.0096,0.0088,0.0082,0.0077,0.0074,0.0072,0.0070,0.0070,0.0070,0.0071,0.0073,0.0077,0.0081,0.0086,0.0090,0.0095,0.0100,0.0105,0.0110},
  {0.0283,0.0269,0.0255,0.0244,0.0235,0.0229,0.0222,0.0216,0.0199,0.0183,0.0170,0.0159,0.0148,0.0138,0.0127,0.0117,0.0107,0.0097,0.0089,0.0082,0.0077,0.0074,0.0071,0.0069,0.0068,0.0068,0.0070,0.0072,0.0077,0.0082,0.0088,0.0095,0.0100,0.0105,0.0111,0.0115},
  {0.0281,0.0266,0.0252,0.0241,0.0232,0.0225,0.0218,0.0212,0.0195,0.0179,0.0165,0.0153,0.0143,0.0133,0.0123,0.0114,0.0104,0.0094,0.0087,0.0081,0.0077,0.0074,0.0071,0.0069,0.0069,0.0069,0.0071,0.0073,0.0078,0.0082,0.0088,0.0092,0.0098,0.0102,0.0106,0.0109},
  {0.0299,0.0285,0.0273,0.0262,0.0254,0.0248,0.0242,0.0236,0.0221,0.0207,0.0195,0.0185,0.0174,0.0161,0.0148,0.0135,0.0123,0.0111,0.0101,0.0093,0.0088,0.0085,0.0083,0.0082,0.0084,0.0088,0.0094,0.0102,0.0112,0.0122,0.0132,0.0140,0.0142,0.0144,0.0147,0.0150},
  {0.0300,0.0286,0.0273,0.0262,0.0254,0.0248,0.0242,0.0237,0.0221,0.0207,0.0195,0.0186,0.0176,0.0162,0.0150,0.0138,0.0127,0.0117,0.0111,0.0108,0.0106,0.0107,0.0110,0.0114,0.0119,0.0127,0.0140,0.0154,0.0168,0.0182,0.0196,0.0198,0.0200,0.0202,0.0204,0.0206},
  {0.0419,0.0409,0.0400,0.0393,0.0387,0.0383,0.0380,0.0376,0.0367,0.0358,0.0352,0.0346,0.0341,0.0320,0.0304,0.0290,0.0280,0.0275,0.0276,0.0283,0.0293,0.0306,0.0321,0.0339,0.0356,0.0385,0.0423,0.0461,0.0502,0.0525,0.0525,0.0526,0.0527,0.0527,0.0528,0.0529},
  {0.0344,0.0332,0.0321,0.0312,0.0305,0.0301,0.0296,0.0292,0.0279,0.0268,0.0259,0.0252,0.0241,0.0226,0.0214,0.0205,0.0198,0.0196,0.0197,0.0202,0.0210,0.0219,0.0229,0.0241,0.0252,0.0272,0.0297,0.0323,0.0336,0.0336,0.0337,0.0338,0.0339,0.0340,0.0341,0.0343},
  {0.0344,0.0332,0.0321,0.0313,0.0306,0.0301,0.0296,0.0292,0.0279,0.0268,0.0259,0.0252,0.0241,0.0227,0.0215,0.0206,0.0199,0.0196,0.0197,0.0203,0.0210,0.0219,0.0230,0.0241,0.0253,0.0261,0.0261,0.0262,0.0264,0.0265,0.0266,0.0267,0.0268,0.0269,0.0271,0.0272},
  {0.0344,0.0332,0.0322,0.0313,0.0306,0.0301,0.0297,0.0292,0.0280,0.0269,0.0260,0.0253,0.0241,0.0227,0.0215,0.0206,0.0199,0.0197,0.0198,0.0203,0.0210,0.0213,0.0213,0.0212,0.0213,0.0213,0.0213,0.0214,0.0215,0.0216,0.0217,0.0219,0.0220,0.0222,0.0224,0.0225},
  {0.0344,0.0332,0.0322,0.0313,0.0306,0.0301,0.0297,0.0292,0.0280,0.0269,0.0260,0.0253,0.0241,0.0228,0.0216,0.0207,0.0200,0.0197,0.0194,0.0191,0.0190,0.0189,0.0188,0.0187,0.0187,0.0187,0.0187,0.0188,0.0189,0.0191,0.0192,0.0194,0.0196,0.0197,0.0199,0.0201},
  {0.0344,0.0332,0.0322,0.0313,0.0306,0.0301,0.0297,0.0292,0.0280,0.0270,0.0260,0.0253,0.0241,0.0228,0.0220,0.0216,0.0211,0.0208,0.0205,0.0203,0.0201,0.0201,0.0200,0.0200,0.0200,0.0200,0.0200,0.0200,0.0202,0.0203,0.0204,0.0206,0.0207,0.0209,0.0211,0.0213}
};
