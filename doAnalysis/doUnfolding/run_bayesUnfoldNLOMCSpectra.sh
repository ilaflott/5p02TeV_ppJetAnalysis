#!/bin/bash

##CONST
#jtID=1
#useNPCorr=1
##CONST
if [ $# -eq 7 ]
then
    R=$1
    condorDir=$2
    date_output=$3
    etaBin=$4
    etaBinOut=$5
    useSimpleBinning=$6
    BayeskIter=$7
else
    echo "usage:"
    echo "source run_bayesUnfoldNLOMCSpectra.sh [R=3,4] [condorDir] [date_output] [etaBin] [etaBinOut] [useSimpleBinning] [BayeskIter]"
    return
fi




#### unfoldMCSpectra
echo ""
echo "running unfoldMCSpectra"
echo ""


./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_ptmax1000_gaussSmear_00eta20.root MC_NNPDFnnlo_closureTest_fitwgts_ptmax1000_${etaBinOut} 1 ${useSimpleBinning} ${BayeskIter}
#return
#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_spl3wgts_ptmax1000_extv2.2_gaussSmear_00eta20.root MC_NNPDFnnlo_closureTest_spl3wgts_ptmax1000_extv2.2_${etaBinOut} 1 ${useSimpleBinning} ${BayeskIter}
#return
#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_ptmax686_gaussSmear_00eta20.root MC_NNPDFnnlo_closureTest_fitwgts_ptmax686_${etaBinOut} 1 ${useSimpleBinning} ${BayeskIter}
#return
./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_spl3wgts_ptmax686_gaussSmear_00eta20.root MC_NNPDFnnlo_closureTest_spl3wgts_ptmax686_${etaBinOut} 1 ${useSimpleBinning} ${BayeskIter}
return


#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_ptmax686_gaussSmear_00eta20.root MC_NNPDFnnlo_closureTest_${etaBinOut} 0 ${useSimpleBinning}
#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_ptmax1000_gaussSmear_00eta20.root MC_NNPDFnnlo_closureTest_${etaBinOut} 0 ${useSimpleBinning}
#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_spl3wgts_ptmax686_gaussSmear_00eta20.root MC_NNPDFnnlo_closureTest_${etaBinOut} 0 ${useSimpleBinning}
return





#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_00eta20.root MC_CT10nlo_closureTest_${etaBinOut} 0 ${useSimpleBinning}
#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_CT14nlo_HS_gaussSmear_00eta20.root MC_CT14nlo_closureTest_${etaBinOut} 0 ${useSimpleBinning}
#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS_gaussSmear_00eta20.root  MC_HERAPDF_closureTest_${etaBinOut} 0 ${useSimpleBinning}
#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS_gaussSmear_00eta20.root	   MC_MMHTnlo_closureTest_${etaBinOut} 0 ${useSimpleBinning}

#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_00eta20.root MC_CT10nlo_closureTest_${etaBinOut} 1 ${useSimpleBinning}
#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_CT14nlo_HS_gaussSmear_00eta20.root MC_CT14nlo_closureTest_${etaBinOut} 1 ${useSimpleBinning}
#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS_gaussSmear_00eta20.root  MC_HERAPDF_closureTest_${etaBinOut} 1 ${useSimpleBinning}
#./bayesUnfoldNLOMCSpectra.exe /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory  /fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS_gaussSmear_00eta20.root	   MC_MMHTnlo_closureTest_${etaBinOut} 1 ${useSimpleBinning}
