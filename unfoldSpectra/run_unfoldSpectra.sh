#!/bin/bash

root -b  <<EOF
.L unfoldSpectra.C
bool t = true
unfoldSpectra(4, t, t, t, 21, "Data")
.q
EOF
#.L utilities.h
#unfold(4, t, t, t, 21, "MC")
