#!/bin/bash

rootcompile quickPrintNLOSpectra.C
echo "quickPrinting NLOSpectra"

#FILE=fNLOJetsSpectra/R04/fnl5020_LO2_R04Jets_modify_CT10nlo_HS.root
#./quickPrintNLOSpectra.exe ${FILE} 1
#./quickPrintNLOSpectra.exe ${FILE} 2
#./quickPrintNLOSpectra.exe ${FILE} 3
#./quickPrintNLOSpectra.exe ${FILE} 4
#
#FILE=fNLOJetsSpectra/R04/fnl5020_LO2_R04Jets_modify_CT14nlo_HS.root
#./quickPrintNLOSpectra.exe ${FILE} 1
#./quickPrintNLOSpectra.exe ${FILE} 2
#./quickPrintNLOSpectra.exe ${FILE} 3
#./quickPrintNLOSpectra.exe ${FILE} 4
#
#
#FILE=fNLOJetsSpectra/R04/fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS.root
#./quickPrintNLOSpectra.exe ${FILE} 1
#./quickPrintNLOSpectra.exe ${FILE} 2
#./quickPrintNLOSpectra.exe ${FILE} 3
#./quickPrintNLOSpectra.exe ${FILE} 4
#
#
#FILE=fNLOJetsSpectra/R04/fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS.root
#./quickPrintNLOSpectra.exe ${FILE} 1
#./quickPrintNLOSpectra.exe ${FILE} 2
#./quickPrintNLOSpectra.exe ${FILE} 3
#./quickPrintNLOSpectra.exe ${FILE} 4

FILE=fNLOJetsSpectra/R04/fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC.root
./quickPrintNLOSpectra.exe ${FILE} 1
./quickPrintNLOSpectra.exe ${FILE} 2
./quickPrintNLOSpectra.exe ${FILE} 3
./quickPrintNLOSpectra.exe ${FILE} 4

scp2Serin fNLOJetsSpectra/R04/\*.pdf
