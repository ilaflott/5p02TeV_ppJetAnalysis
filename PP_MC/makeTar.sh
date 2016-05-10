#!bin/bash

root -l <<EOF
.L RAA_read_mc_pp.C++
.q
EOF


rm jetRAA_run_pp_MC.tar
tar -zcvf jetRAA_run_pp_MC.tar *.txt RAA_read_*.* weights_*.root
