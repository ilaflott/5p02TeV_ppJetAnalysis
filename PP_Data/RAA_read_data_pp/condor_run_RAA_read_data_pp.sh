#!/bin/bash

# untar contents of tarball submitted
tar -xvf run_RAA_read_pp_data.tar > /dev/null

# env variables
export SCRAM_ARCH=slc6_amd64_gcc491
source /osg/app/cmssoft/cms/cmsset_default.sh

# debug
echo ""
echo "----------------------------------------------------"
echo "Job started on `date` at WN: `hostname` "
echo "Job is running on `uname -a`"

echo "root directory is..."
echo $ROOTSYS

echo "free disk space is..."
df -h

echo "gcc verion is..."
gcc --version

echo "working directory is..."
pwd

echo "contents before job processing..."
ls -al

# run script inputs, defined and set in corresponding submit script
startfile=$1
endfile=$2
radius=$3
jetType=$4
outfile=$5
destination=$6

echo "Processing..."
root -b -l <<EOF
.x RAA_read_data_pp.C+($startfile,$endfile,$radius,"$jetType","$outfile")
.q
EOF

#destination="/export/d00/scratch/ilaflott/5p02TeV_ppJetAnalysis/condor_output/PP_Data/RAA_read_data_pp"
echo "Copying file ${outfile} to ${destination}..."
mv ${outfile} ${destination}

echo "Done!"