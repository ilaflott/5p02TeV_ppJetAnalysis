#!/bin/bash



if [[ "$#" -eq 0 ]]
then
    
    
    ##DEBUG
    #rootcompile smearTheorySpectra_gaussCoreJER.C    
    
    echo ""
    echo "RUN NLO smearing"
    echo ""

    #./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_CT10nlo_HS
    #./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_CT14nlo_HS
    #./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS 
    #./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS
    #./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC 	    
    
    #./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC 0    #fitweights
    ./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC 1    #splineweights
    
    #rm NNPDF_NNLO_00eta20_fitwgts_gausSmThy_plots/*
    rm NNPDF_NNLO_00eta20_spl3wgts_gausSmThy_plots/*
    
    root -l -b -q quickCheck_gausSmThy.C++
    
    echo ""
    echo "DONE NLO smearing."
    echo "DONE NLO smearing. scp2Serin"
    echo ""
    
    #scp2Serin fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_00eta20.root
    ##scp2Serin fnl5020_LO2_R04Jets_modify_CT14nlo_HS_gaussSmear_00eta20.root
    ##scp2Serin fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS_gaussSmear_00eta20.root
    ##scp2Serin fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS_gaussSmear_00eta20.root
    ##scp2Serin fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_gaussSmear_00eta20.root    
    #scp2Serin "fnl5020_LO2_R04Jets_modify_*_gaussSmear_00eta20.root"
    
    #cp fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_gaussSmear_00eta20.root NNPDF_NNLO_00eta20_gausSmThy_plots/.
    #cp fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_fitwgts_gaussSmear_00eta20.root NNPDF_NNLO_00eta20_fitwgts_gausSmThy_plots/.
    cp fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_spl3wgts_gaussSmear_00eta20.root NNPDF_NNLO_00eta20_spl3wgts_gausSmThy_plots/.
    
    scp2Serin NNPDF_NNLO_00eta20_spl3wgts_gausSmThy_plots    
    
    #scp2Serin NNPDF_NNLO_00eta20_fitwgts_gausSmThy_plots    
    #scp2Serin "NNPDF_NNLO_00eta20_*_gausSmThy_plots"
    


elif [[ "$#" -eq 1 ]]
then
    etaBin=$1
    echo ""
    echo "compiling JER_gaussCore_smearTheoryCurves_singleYbin"
    echo ""
    
    rooUnfoldCompile smearTheorySpectra_gaussCoreJER_singleYbin.C
    
    echo ""
    echo "RUN NLO smearing"
    echo ""
    
    ./smearTheorySpectra_gaussCoreJER_singleYbin.exe fnl5020_LO2_R04Jets_modify_CT10nlo_HS $etaBin
#./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_CT14nlo_HS $etaBin
#./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS $etaBin
#./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS $etaBin
#./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC $etaBin
    
    echo ""
    echo "DONE NLO smearing."
    echo "DONE NLO smearing. scp2Serin"
    echo ""
    
    outfile="fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_etaBin${etaBin}.root"
    scp2Serin $outfile
elif [[ "$#" -eq 2 ]]
then
    
    echo ""
    echo "compiling JER_gaussCore_smearTheoryCurves_singleYbin"
    echo ""
    
    rooUnfoldCompile smearTheorySpectra_gaussCoreJER_singleYbin.C
    
    echo ""
    echo "RUN NLO smearing"
    echo ""
    
    ./smearTheorySpectra_gaussCoreJER_singleYbin.exe fnl5020_LO2_R04Jets_modify_CT10nlo_HS 1
    ./smearTheorySpectra_gaussCoreJER_singleYbin.exe fnl5020_LO2_R04Jets_modify_CT10nlo_HS 2
    ./smearTheorySpectra_gaussCoreJER_singleYbin.exe fnl5020_LO2_R04Jets_modify_CT10nlo_HS 3
    ./smearTheorySpectra_gaussCoreJER_singleYbin.exe fnl5020_LO2_R04Jets_modify_CT10nlo_HS 4
#./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_CT14nlo_HS $etaBin
#./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS $etaBin
#./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS $etaBin
#./smearTheorySpectra_gaussCoreJER.exe fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC $etaBin
    
    echo ""
    echo "DONE NLO smearing."
    echo "DONE NLO smearing. scp2Serin"
    echo ""
    
    scp2Serin fnl5020_LO2_R04Jets_modify_CT10nlo_HS_gaussSmear_etaBin\?.root
else 
    echo "bad input. see macro's insides."
fi

echo ""
echo "done w/ run_smearTheory.sh"
echo ""

return

