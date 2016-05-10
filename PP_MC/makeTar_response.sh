#!bin/bash

root -l <<EOF
.L deriveDijetResponse.C++
.q
EOF

rm run_pp_data_response.tar
tar -zcvf run_pp_data_response.tar pp_*.txt deriveDijet*.*
