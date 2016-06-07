#!bin/bash

root -l <<EOF
.L RAA_read_data_pp.C++
.q
EOF

rm jetRAA_run_pp_Data.tar
tar -zcvf jetRAA_run_pp_Data.tar pp_*.txt RAA_read_*.*
