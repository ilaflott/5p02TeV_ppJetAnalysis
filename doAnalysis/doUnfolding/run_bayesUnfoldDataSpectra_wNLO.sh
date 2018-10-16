#!/bin/bash


##CONST
jtID=1
##CONST

if [ $# -eq 10 ]
then
    R=$1
    data_condorDir=$2
    data_date=$3
    dataetaBin=$4
    MC_condorDir=$5
    MC_date=$6
    MCetaBin=$7
    etaBinOut=$8
    useSimpleBinning=$9
    BayeskIter=${10}
else
    echo "usage:"
    echo "source run_bayesUnfoldDataSpectra_wNLO.sh [R=3,4] [dataCondorDir] [data_date] [MCCondorDir] [MC_date] [MCetaBin] [etaBinOut] [useSimpleBinning] [BayeskIter]"
    return
fi

#### unfoldDataSpectra
echo ""
echo "running bayesUnfoldDataSpectra_wNLO"
echo ""

#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory /fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_00eta20.root              unf_HPtJetTrig_${etaBinOut}_CT10nlo      1 ${jtID} ${useSimpleBinning}
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory /fnl5020_LO2_R04Jets_modify_CT14nlo_HS_gaussSmear_00eta20.root              unf_HPtJetTrig_${etaBinOut}_CT14nlo      1 ${jtID} ${useSimpleBinning}
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory /fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS_gaussSmear_00eta20.root  unf_HPtJetTrig_${etaBinOut}_HERAPDF15NLO 1 ${jtID} ${useSimpleBinning}
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory /fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS_gaussSmear_00eta20.root      unf_HPtJetTrig_${etaBinOut}_MMHT2014nlo  1 ${jtID} ${useSimpleBinning}

#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_ptmax686_gaussSmear_00eta20.root unf_HPtJetTrig_${etaBinOut}_NNPDF30nnlo  1 ${jtID} ${useSimpleBinning}
./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_ptmax1000_gaussSmear_00eta20.root unf_HPtJetTrig_fitwgts_ptmax1000_${etaBinOut}_NNPDF30nnlo  1 ${jtID} ${useSimpleBinning} ${BayeskIter}
#return
./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_ptmax686_gaussSmear_00eta20.root unf_HPtJetTrig_fitwgts_ptmax686_${etaBinOut}_NNPDF30nnlo  1 ${jtID} ${useSimpleBinning} ${BayeskIter}
#return
./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_spl3wgts_ptmax686_gaussSmear_00eta20.root unf_HPtJetTrig_spl3wgts_ptmax686_${etaBinOut}_NNPDF30nnlo  1 ${jtID} ${useSimpleBinning} ${BayeskIter}
return



#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory /fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_00eta20.root              unf_HPtJetTrig_${etaBinOut}_CT10nlo      0 ${jtID} ${useSimpleBinning}
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory /fnl5020_LO2_R04Jets_modify_CT14nlo_HS_gaussSmear_00eta20.root              unf_HPtJetTrig_${etaBinOut}_CT14nlo      0 ${jtID} ${useSimpleBinning}
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory /fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS_gaussSmear_00eta20.root  unf_HPtJetTrig_${etaBinOut}_HERAPDF15NLO 0 ${jtID} ${useSimpleBinning}
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory /fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS_gaussSmear_00eta20.root      unf_HPtJetTrig_${etaBinOut}_MMHT2014nlo  0 ${jtID} ${useSimpleBinning}


#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_ptmax686_gaussSmear_00eta20.root unf_HPtJetTrig_${etaBinOut}_NNPDF30nnlo  0 ${jtID} ${useSimpleBinning}
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_ptmax1000_gaussSmear_00eta20.root unf_HPtJetTrig_${etaBinOut}_NNPDF30nnlo  0 ${jtID} ${useSimpleBinning}
./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_spl3wgts_ptmax686_gaussSmear_00eta20.root unf_HPtJetTrig_spl3wgts_ptmax686_${etaBinOut}_NNPDF30nnlo  0 ${jtID} ${useSimpleBinning}










return