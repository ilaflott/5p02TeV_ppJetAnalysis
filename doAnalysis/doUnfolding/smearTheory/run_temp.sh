#!/bin/bash


rootcompile temp.C

parambinstart=5
nbins=20
parambin=$parambinstart
while [ $parambin -le $nbins ]
do
    ./temp.exe $parambin
    
    parambin=$(($parambin + 1))
done

#rm         logFits_NNPDF_NNLO_POWPY8_NPs/*
#mv *param* logFits_NNPDF_NNLO_POWPY8_NPs/.
#ls -1      logFits_NNPDF_NNLO_POWPY8_NPs/
#scp2Serin  logFits_NNPDF_NNLO_POWPY8_NPs/

rm         logFits_NNPDF_NNLO/*
mv *param* logFits_NNPDF_NNLO/.
ls -1      logFits_NNPDF_NNLO/
scp2Serin  logFits_NNPDF_NNLO/

return