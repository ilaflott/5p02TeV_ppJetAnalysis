#!/bin/bash

rootcompile vzEvtWeights.C
./vzEvtWeights.exe

scp2Serin \*png
scp2Serin Py8_vzEvtWeights_hVz_baryshift_neg_04_11_18.root