#!/bin/bash




#!/bin/bash

Nasks=$1
T=$2
Ncount=0

#secondCount=0
while [[ $Ncount -lt $Nasks ]]
do





    
    
    
    
    
######################################################################################
    
    echo ""
    echo "num of jet80 files copied over so far is..."
    echo ""
    
    ls -d /cms/se/store/user/ilaflott/jobOutput2/2015pp_HighPtJet80/
    ls -lt /cms/se/store/user/ilaflott/jobOutput2/2015pp_HighPtJet80/*.root | wc -l
    echo "out of...."
    cat 5p02TeV_HighPtJet80_forests.txt | wc -l
    
######################################################################################

    echo ""
    echo "num of lower jets files copied over so far is..."
    echo ""
    
    ls -d /cms/se/store/user/ilaflott/jobOutput2/2015pp_HighPtLowerJets/
    ls -lt /cms/se/store/user/ilaflott/jobOutput2/2015pp_HighPtLowerJets/*.root | wc -l
    echo "out of...."
    cat 5p02TeV_HighPtLowerJets_forests.txt | wc -l
    
######################################################################################
    

    echo ""
    echo "num of Py8 files copied over so far is..."
    echo ""
    
    
    ls -d /cms/se/store/user/ilaflott/jobOutput2/Pythia8_Dijet_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/dijet*/
    ls -lt /cms/se/store/user/ilaflott/jobOutput2/Pythia8_Dijet_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/dijet*/*.root | wc -l
    echo "out of...."
    cat 5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt | wc -l
    
######################################################################################
    
    
    
    echo ""
    echo "num of merged Py8 files copied over so far is..."
    echo ""
    
    
    ls -d /cms/se/store/user/ilaflott/jobOutput2/Pythia8_Dijet_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/merged/
    ls -lt /cms/se/store/user/ilaflott/jobOutput2/Pythia8_Dijet_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/merged/*.root | wc -l
    echo "out of...."
    ls -l ../readForests/filelists/5p02TeV_Py8_CUETP8M1_QCDjet_Official_forests_individualBins/*.txt | wc -l
    
######################################################################################
    
    Ncount=$(( $Ncount + 1 ))
    sleep ${T}s
done

echo ""
echo "done asking"
echo ""

return
