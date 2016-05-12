#!/bin/bash
# usage: source condor_makeTarAndSubmit.sh NJobs NFilesPerJob

##
##COMPILE CODE AND CREATE TARBALL OF INPUT FILES##
##

# compile the script
root -l <<EOF
.L RAA_read_data_pp.C+
.q
EOF

#remove old tar ball if there and replace with new one
rm run_RAA_read_pp_data.tar

#create tar from filelist and the compiled macro
tar -zcvf run_RAA_read_pp_data.tar *pp_data_forests.txt RAA_read_data_pp*.*
rm *.d *.so *.pcm

echo "code compiled and tarball created"

##
##CREATE SUBMIT FILE AND SUBMIT TO CONDOR
##

#inputs to make the run script
NJobs=$1
NFilesPerJob=$2
radius=4
jetType="PF"
destination="/export/d00/scratch/ilaflott/5p02TeV_ppJetAnalysis/condor_output/PP_Data"

#input tarball location
tardir=`pwd`

#load the filelist to determine job splitting, should match filelist in tarball
#filelist=HighPtJet80_5p02TeV_pp_data_forests.txt
filelist=testFile_2015_5p02TeV_pp_data_forests.txt

#one condor job submit per NFilesPerJob until we submit NJobs
nFiles=`wc -l < $filelist`
echo "nFiles in list: $nFiles"

counter=0
NJobsSubmitted=0
while [ $counter -lt $NJobs ]
do

    # job splitting
    startfile=$(( $counter * $NFilesPerJob ))
    endfile=$(( ($counter + 1) * $NFilesPerJob ))

    # check for end of filelist
    if [ $endfile -gt $nFiles ]; then
        let endfile=$nFiles
        let counter=$NJobs
    fi

    # define output names for job submission
    Error="ak$radius-ppData-${startfile}_to_${endfile}.err"
    Output="ak$radius-ppData-${startfile}_to_${endfile}.out"
    Log="ak$radius-ppData-${startfile}_to_${endfile}.log"
    outfile="5p02TeV_pp_data_ak${radius}${jetType}_20_eta_20_${startfile}_to_${endfile}.root"
    
    # create the condor submit file
    cat > subfile <<EOF

Universe       = vanilla
Environment = "HOSTNAME=$HOSTNAME"
#tell condor where grid certificate is
#x509userproxy=/tmp/x509up_u2142
# tell it what to run
Executable     = condor_run_RAA_read_data_pp.sh
+AccountingGroup = "group_cmshi.ilaflott"
#arguments for the run script
Arguments      = $startfile $endfile $radius $jetType $outfile $destination
# input files to run script, none in this case
Input          = /dev/null
# log files
Error          = $destination/logs/$Error
Output         = $destination/logs/$Output
Log            = $destination/logs/$Log
# get the environment (path, etc.)
Getenv         = True
# prefer to run on fast computers
Rank           = kflops
# only run on 64 bit computers
Requirements   = Arch == "X86_64"
transfer_input_files = $tardir/run_RAA_read_pp_data.tar
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
# specify any extra input files (for example, an orcarc file)
Queue
EOF

    # submit the job defined in the above submit file
    echo "submitting RAA_read_data_pp job #${counter}..." 
    condor_submit subfile
    counter=$(($counter + 1))
    echo "submitted job # ${counter} of ${NJobs}"
    rm subfile
done
