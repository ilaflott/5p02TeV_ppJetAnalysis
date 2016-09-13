9/6/2016, a bug is preventing me from opening files locally on hidsk, the file are stored on hadoop

noticed that condor jobs would run, and could still access the files as normal. So i wrote 
printRootFileContent.C to figure out what collections were available, and used it as an opportunity
to motivate coding up a more "generic" condor script, that takes in a code, compiles it, then 
submits the executable to condor for running.