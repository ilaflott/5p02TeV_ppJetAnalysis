#!/bin/bash

root -b  <<EOF
.L utilities.h
.L unfold.C
bool t = true
unfold(4, t, t, t, 21, "Data")
.q
EOF
#unfold(4, t, t, t, 21, "MC")
