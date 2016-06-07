#!bin/bash

root -l <<EOF
.L RAA_read_data_pbpb.C++
.q
EOF

rm jetRAA_run_PbPb_Data.tar
tar -zcvf jetRAA_run_PbPb_Data.tar jetRAA_PbPb_*.txt RAA_read_*.* 
