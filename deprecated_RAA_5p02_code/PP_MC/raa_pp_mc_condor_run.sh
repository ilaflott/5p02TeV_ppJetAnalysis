#!/bin/bash

tar -xvf jetRAA_run_pp_MC.tar > /dev/null

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
jetType=$4
outfile=$5
algo=$6
echo "Processing..."

root -b -l <<EOF
.x RAA_read_mc_pp.C+($startfile,$endfile,$radius, "$algo", "$jetType","$outfile")
.q
EOF

# root -b -q RAA_read_mc_pbpb.C\+\($startfile,$endfile\)
mv $outfile /mnt/hadoop/cms/store/user/rkunnawa/rootfiles/JetRAA_5p02TeV/Jan2016/.

echo "Done!"

echo "Copied output files to hadoop rootfiles/JetRAA_5p02TeV/Nov24/" 
