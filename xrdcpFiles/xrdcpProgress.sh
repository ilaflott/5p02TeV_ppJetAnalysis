#!/bin/bash



######################################################################################

echo ""
echo "num of jet80 files copied over so far is..."
echo ""

ls -lt /cms/se/store/user/ilaflott/jobOutput2/2015pp_HighPtJet80/*.root | wc -l


echo ""
echo "out of...."
echo ""


cat 5p02TeV_HighPtJet80_forests.txt | wc -l

######################################################################################

echo ""
echo "num of lower jets files copied over so far is..."
echo ""

ls -lt /cms/se/store/user/ilaflott/jobOutput2/2015pp_HighPtLowerJets/*.root | wc -l


echo ""
echo "out of...."
echo ""


cat 5p02TeV_HighPtLowerJets_forests.txt | wc -l

######################################################################################

echo ""
echo "done"
echo ""

return