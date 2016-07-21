#!/bin/bash

# env setup for cgate/hidsk/T2_US_MIT cluster

echo 'sourcing cms software'
source /osg/app/cmssoft/cms/cmsset_default.sh
echo 'cmsenv'
cmsenv
echo 'voms initializing...'
voms-proxy-init -voms cms

# specific to this CMSSW area
export CMSSW_BASE=${PWD/%"/src"/}
#CMSSW_BASE=${CMSSW_BASE/%"/src"/}
export ROOUNFOLD_PATH=/net/hisrv0001/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/RooUnfold-svnTrunk_6.27.16/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOUNFOLD_PATH
source /net/hisrv0001/home/ilaflott/.bashrc # have to do this *again* for cdsrc alias to work
                                            # this gives me the willies, a bit
                                            # if odd behavior later, this could be a potential source of the problem