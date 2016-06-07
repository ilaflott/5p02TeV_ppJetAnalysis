#!/bin/bash

tar -xvf jetRAA_run_PbPb_MC.tar > /dev/null

export SCRAM_ARCH=slc6_amd64_gcc491
source /osg/app/cmssoft/cms/cmsset_default.sh

echo ""
echo "----------------------------------------------------"
echo "Job started on `date` at WN: `hostname` "
echo "Job is running on `uname -a`"

df -h

gcc --version

startfile=$1
endfile=$2
radius=$3
algo=$4
jetType=$5
outfile=$6
echo "Processing..."

root -b -l <<EOF
.x RAA_read_mc_pbpb.C+($startfile,$endfile,$radius,"$algo","$jetType","$outfile")
.q
EOF

# root -b -q RAA_read_mc_pbpb.C\+\($startfile,$endfile\)
mv $outfile /mnt/hadoop/cms/store/user/rkunnawa/rootfiles/JetRAA_5p02TeV/Nov24/.

echo "Done!"

echo "Copied output files to hadoop" 
