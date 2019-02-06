#!/bin/bash


#echo "sourcing condor+OSG software"
source /condor/HTCondor/current/condor.sh
source /osg/current/setup.sh

echo "setting environment variables"
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
export SCRAM_ARCH=slc6_amd64_gcc491

echo "sourcing cmssetup defaults"
source $VO_CMS_SW_DIR/cmsset_default.sh

##other env variables, no need to output to screen
#export CMSSW_GIT_REFERENCE=$VO_CMS_SW_DIR/cmssw.git
#export CVMFS_758=${VO_CMS_SW_DIR}/${SCRAM_ARCH}/cms/cmssw/CMSSW_7_5_8/src/

echo 'cmsenv'
cmsenv

##setup new git version
#export PATH=/cms/base/git-2.9.5/bin:${PATH}

## don't setup CRAB3 unless running grid jobs.
## Mixins.py error pops up if using cmsRun locally with CRAB3 envsetup because crab needs python 2.6, while CMSSW uses 2.7
## if CRAB3 + CMSSW must both be run in same session with one env setup, lookup the "bootstrap script" in CRAB3 twiki literature
#echo 'setting up CRAB3'
#source $VO_CMS_SW_DIR/crab3/crab.sh

#echo 'obtaining voms proxy'
#voms-proxy-init -voms cms

# convenient env variables
export CMSSW_BASE=${PWD/%"/src"/}
alias cdsrc="cd ${CMSSW_BASE}/src"

#export ROOUNFOLD_PATH=${CMSSW_BASE}/src/RooUnfold-svnTrunk
export ROOUNFOLD_PATH=${CMSSW_BASE}/src/RooUnfold-svnTrunk_master
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOUNFOLD_PATH



####OLD
##SCRAM_ARCH for specific CMSSW_XYZ area
#export SCRAM_ARCH=slc6_amd64_gcc491

#source /condor/HTCondor/current/condor.sh
#source /osg/current/setup.sh
###deprecated, instead do above
##source /osg/osg3.2/osg-client/setup.sh

#echo 'sourcing cms software'
#source $VO_CMS_SW_DIR/cmsset_default.sh
