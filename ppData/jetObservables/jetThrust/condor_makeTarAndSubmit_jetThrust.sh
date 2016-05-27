
root -l <<EOF
.L thrust_HiForest.C++
.q
EOF


rm thrust_pp_MC.tar
tar -zcvf thrust_pp_MC.tar jetRAA_pp_*.txt thrust_HiFore*.* weights_pp*.root
