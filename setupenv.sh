# source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.csh
echo 'sourcing cms software'
source /osg/app/cmssoft/cms/cmsset_default.sh
echo 'cmsenv'
cmsenv
# source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.csh
echo 'voms initializing...'
voms-proxy-init -voms cms
