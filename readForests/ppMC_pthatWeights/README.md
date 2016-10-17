pthatEvtWeights is stuck into readFiles/readForests_ppMC instead of computeWeights/ intentionally.

this is because pthatEvtWeights runs on the forests in the same way as the readForests routine, 
and requires only the MC dataset to be correctly computed. Someone else did the comparison
with real data before us, for us, and gave us cross-section weights to be used in this later
computation

this is unlike computeWeights/vzEvtWeights.C in this sense; vzEvtWeights runs on histogram
output from both readForests_ppData and ppMC, and these weights are then meant to be applied
to the events in the same way the pthatEvtWeights are.

long story short; the two weight computations are not analogous at this point in time, 
because one depends on output from readFiles/Forests like unfoldSpectra or printQA plots, 
and the other goes over the forests like readFiles/Forests.

so future ian, just think about this before you decide to change folder structure again!