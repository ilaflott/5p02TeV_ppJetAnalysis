#!/bin/bash

# run this with $ source raa_condor_submit.sh NJobs NFilesPerJob

counter=0
incrementer=1

destination=/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/JetRAA/June29/
filelist=jetRAA_PbPb_mc_forests.txt

nFiles=`wc -l < $filelist`
tardir=`pwd`
radius=4
algo="Vs"
jetType="PF"

echo "nFiles in list: $nFiles"
while [ $counter -lt $1 ]
do
    echo $counter >> Submitted
    
    Error="ak$algo$radius$jetType-PbPbMC-$endfile.err"
    Output="ak$algo$radius$jetType-PbPbMC-$endfile.out"
    Log="ak$algo$radius$jetType-PbPbMC-$endfile.log"
    
    startfile=$(( $counter * $2 ))
    endfile=$(( ($counter + 1) * $2 ))
    if [ $endfile -gt $nFiles ]; then
        let endfile=$nFiles
        let counter=$1
    fi

    outfile="PbPb_MC_5p02TeV_ak${algo}${radius}${jetType}_20_eta_20_${endfile}.root"
    
    # Condor submit file
    cat > subfile <<EOF

Universe       = vanilla
Environment = "HOSTNAME=$HOSTNAME"
# files will be copied back to this dir
# Initialdir     = .
#tell condor where my grid certificate it
#x509userproxy=/tmp/x509up_u2142
# run my script
Executable     = raa_pbpb_mc_condor_run.sh
+AccountingGroup = "group_cmshi.rkunnawa"
#+IsMadgraph = 1
Arguments      = $startfile $endfile $radius $algo $jetType $outfile
# input files. in this case, there are none.
Input          = /dev/null
# log files
Error          = LOG/$Error
Output         = LOG/$Output
Log            = LOG/$Log
# get the environment (path, etc.)
Getenv         = True
# prefer to run on fast computers
Rank           = kflops
# only run on 64 bit computers
Requirements   = Arch == "X86_64"
transfer_input_files = $tardir/jetRAA_run_PbPb_MC.tar
# should write all output & logs to a local directory
# and then transfer it back to Initialdir on completion
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
# specify any extra input files (for example, an orcarc file)
Queue
EOF

    # submit the job
    echo "submitting condor_run.sh $startfile $endfile to $destination ..." 
    condor_submit subfile
    counter=$(($counter + 1))
done
