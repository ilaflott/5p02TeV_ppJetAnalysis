#!/bin/bash


##CONST
jtID=1
NLOFileDir="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory"
##CONST

if [ $# -eq 8 ]
then
    R=$1
    data_condorDir=$2
    data_date=$3
    dataetaBin=$4
    etaBinOut=$5 
    useSimpleBinning=$6 
    BayeskIter=$7 
    useNPCorr=$8
else
    echo "usage:"
    echo "source run_bayesUnfoldDataSpectra_wNLO.sh [R=3,4] [dataCondorDir] [data_date] [datadirTag] [etaBinOut] [useSimpleBinning] [BayeskIter] [useNPCorr]"
    return
fi

#### unfoldDataSpectra
echo ""
echo "running bayesUnfoldDataSpectra_wNLO"
echo ""
#fnameFirstHalf="fnl5020_L02_R04Jets_modify_"
#fnameSecondHalf="_fitwgts_ptmax1000_gausSmear_00eta20.root"

#### ./bayesUnfoldDataSpectra_wNLO.exe  <targDataDir> <targMCDir> <targMCName> <baseOutputName> <applyNPCorrs> <doJetID> <useSimpleBins> <kIterInput>


##### BASE NLO FILES
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_00eta20.root              unf_HPtJetTrig_${etaBinOut}_CT10nlo      ${useNPCorr} ${jtID} ${useSimpleBinning}
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_CT14nlo_HS_gaussSmear_00eta20.root              unf_HPtJetTrig_${etaBinOut}_CT14nlo      ${useNPCorr} ${jtID} ${useSimpleBinning}
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS_gaussSmear_00eta20.root  unf_HPtJetTrig_${etaBinOut}_HERAPDF15NLO ${useNPCorr} ${jtID} ${useSimpleBinning}
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS_gaussSmear_00eta20.root      unf_HPtJetTrig_${etaBinOut}_MMHT2014nlo  ${useNPCorr} ${jtID} ${useSimpleBinning}
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_gaussSmear_00eta20.root      unf_HPtJetTrig_${etaBinOut}_NNPDF30nnlo  ${useNPCorr} ${jtID} ${useSimpleBinning}


##### NNPDF VARIATIONS
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_ptmax686_gaussSmear_00eta20.root unf_HPtJetTrig_${etaBinOut}_NNPDF30nnlo  ${useNPCorr} ${jtID} ${useSimpleBinning}
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_ptmax1000_gaussSmear_00eta20.root unf_HPtJetTrig_fitwgts_ptmax1000_${etaBinOut}_NNPDF30nnlo  ${useNPCorr} ${jtID} ${useSimpleBinning} ${BayeskIter}
#return
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_spl3wgts_ptmax1000_extv2.2_gaussSmear_00eta20.root unf_HPtJetTrig_spl3wgts_ptmax1000_extv2.2_${etaBinOut}_NNPDF30nnlo  ${useNPCorr} ${jtID} ${useSimpleBinning} ${BayeskIter}
#return
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_ptmax686_gaussSmear_00eta20.root unf_HPtJetTrig_fitwgts_ptmax686_${etaBinOut}_NNPDF30nnlo  ${useNPCorr} ${jtID} ${useSimpleBinning} ${BayeskIter}
#return



#### ./bayesUnfoldDataSpectra_wNLO.exe  <targDataDir> <targMCDir> <targMCName> <baseOutputName> <applyNPCorrs> <doJetID> <useSimpleBins> <kIterInput>
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_spl3wgts_ptmax686_gaussSmear_00eta20.root unf_HPtJetTrig_spl3wgts_ptmax686_${etaBinOut}_NNPDF30nnlo  ${useNPCorr} ${jtID} ${useSimpleBinning} ${BayeskIter}
./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_spl3wgts_ptmax686_gaussSmear_00eta20.root unf_HPtJetTrig_spl3wgts_ptmax686_${etaBinOut}_NNPDF30nnlo  ${useNPCorr} ${jtID} ${useSimpleBinning} ${BayeskIter}

return















#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_00eta20.root              unf_HPtJetTrig_${etaBinOut}_CT10nlo      ${useNPCorr} ${jtID} ${useSimpleBinning}
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_CT14nlo_HS_gaussSmear_00eta20.root              unf_HPtJetTrig_${etaBinOut}_CT14nlo      ${useNPCorr} ${jtID} ${useSimpleBinning}
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS_gaussSmear_00eta20.root  unf_HPtJetTrig_${etaBinOut}_HERAPDF15NLO ${useNPCorr} ${jtID} ${useSimpleBinning}
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS_gaussSmear_00eta20.root      unf_HPtJetTrig_${etaBinOut}_MMHT2014nlo  ${useNPCorr} ${jtID} ${useSimpleBinning}


#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_ptmax686_gaussSmear_00eta20.root unf_HPtJetTrig_${etaBinOut}_NNPDF30nnlo  ${useNPCorr} ${jtID} ${useSimpleBinning}
#./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_ptmax1000_gaussSmear_00eta20.root unf_HPtJetTrig_${etaBinOut}_NNPDF30nnlo  ${useNPCorr} ${jtID} ${useSimpleBinning}
./bayesUnfoldDataSpectra_wNLO.exe ${data_condorDir}/ppData_HighPtJetTrig_ak${R}PFJets_${data_date}_jetPlots_${dataetaBin} ${NLOFileDir} /fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_spl3wgts_ptmax686_gaussSmear_00eta20.root unf_HPtJetTrig_spl3wgts_ptmax686_${etaBinOut}_NNPDF30nnlo  ${useNPCorr} ${jtID} ${useSimpleBinning}










return