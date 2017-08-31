##SCRAM_ARCH for specific CMSSW_XYZ area
export SCRAM_ARCH="slc6_amd64_gcc491" 

echo "sourcing condor+OSG software"
source /condor/HTCondor/current/condor.sh
source /osg/osg3.2/osg-client/setup.sh

echo 'sourcing cms software'
source "$VO_CMS_SW_DIR/cmsset_default.sh"

echo 'cmsenv'
cmsenv

## don't setup CRAB3 unless running grid jobs.
## Mixins.py error pops up if using cmsRun locally with CRAB3 envsetup because crab needs python 2.6, while CMSSW uses 2.7
## if CRAB3 + CMSSW must both be run in same session with one env setup, lookup the "bootstrap script" in CRAB3 twiki literature

echo 'setting up CRAB3'
source "$VO_CMS_SW_DIR/crab3/crab.sh"
export CRAB3_OUTPUT_DIR=/cms/se/store/user/ilaflott

echo 'obtaining voms proxy'
voms-proxy-init -voms cms

# convenient env variables
export CVMFS_758=/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_5_8/src/

export CMSSW_BASE=${PWD/%"/src"/}
alias cdsrc="cd ${CMSSW_BASE}/src"

export ROOUNFOLD_PATH="${CMSSW_BASE}/src/RooUnfold-svnTrunk"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOUNFOLD_PATH

