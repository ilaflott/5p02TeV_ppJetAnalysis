#!/bin/bash

##CONST
#jtID=1
#useNPCorr=1
##CONST
if [ $# -eq 6 ]
then
    R=$1
    condorDir=$2
    date_output=$3
    etaBin=$4
    etaBinOut=$5
    useSimpleBinning=$6
else
    echo "usage:"
    echo "source run_bayesUnfoldNLOMCSpectra.sh [R=3,4] [condorDir] [date_output] [etaBin] [etaBinOut] [useSimpleBinning]"
    return
fi




#### unfoldMCSpectra
echo ""
echo "running unfoldMCSpectra"
echo ""

#./bayesUnfoldNLOMCSpectra.exe ${condorDir}/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak${R}PFJets_${date_output}_JERS_${etaBin}  Py8_closureTest_${etaBinOut} ${jtID} ${useSimpleBinning}
#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_00eta20.root MC_CT10PDF_closureTest_${etaBinOut} ${useNPCorr} ${useSimpleBinning}

#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_00eta20.root MC_CT10nlo_closureTest_${etaBinOut} 0 ${useSimpleBinning}
#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_CT14nlo_HS_gaussSmear_00eta20.root MC_CT14nlo_closureTest_${etaBinOut} 0 ${useSimpleBinning}
#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS_gaussSmear_00eta20.root  MC_HERAPDF_closureTest_${etaBinOut} 0 ${useSimpleBinning}
#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS_gaussSmear_00eta20.root	   MC_MMHTnlo_closureTest_${etaBinOut} 0 ${useSimpleBinning}
./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_gaussSmear_00eta20.root MC_NNPDFnnlo_closureTest_${etaBinOut} 0 ${useSimpleBinning}




#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_00eta20.root MC_CT10nlo_closureTest_${etaBinOut} 1 ${useSimpleBinning}
#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_CT14nlo_HS_gaussSmear_00eta20.root MC_CT14nlo_closureTest_${etaBinOut} 1 ${useSimpleBinning}
#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS_gaussSmear_00eta20.root  MC_HERAPDF_closureTest_${etaBinOut} 1 ${useSimpleBinning}
#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS_gaussSmear_00eta20.root	   MC_MMHTnlo_closureTest_${etaBinOut} 1 ${useSimpleBinning}
./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_gaussSmear_00eta20.root MC_NNPDFnnlo_closureTest_${etaBinOut} 1 ${useSimpleBinning}













return