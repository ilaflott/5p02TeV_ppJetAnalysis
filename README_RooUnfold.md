ian w laflotte, started on 6/15/2016, last updated on 7/21/16

RooUnfold is a tad tricky to implement at cgate/MIT. So i've noted down the entire process here. I tested it in CMSSW_7_5_8/ROOT6 on hidsk

**********
// GENERAL ROOUNFOLD IMPLEMENTATION --------------------------------------------------------
**********

*** Step 1. 

setupenv and download trunk

In a SCRAM-based CMSSW-area's src folder; download the RooUnfold svn trunk with

     source /osg/app/cmssoft/cms/cmsset_default.sh
     cmsenv
     voms-proxy-init -voms cms
     svn co https://svnsrv.desy.de/public/unfolding/RooUnfold/trunk  #downloads the directory "trunk"
     mv trunk RooUnfold-svnTrunk  #rename for convenience

*** Step 2. 

build/compile RooUnfold library/binaries

self-explanatory; do
         
     cd RooUnfold-svnTrunk
     make
     make bin

*** Step 3. export environment variables 

bash won't know where to look for relevant RooUnfold binaries at run-time when code is run unless we alter LD_LIBRARY_PATH manually. its also useful to tell bash exactly where RooUnfold-svnTrunk is for compiling code later; do

     export ROOUNFOLD_PATH=<your RooUnfold directory in a SCRAM-based CMSSW are>
     export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:ROOUNFOLD_PATH

*** Step 4. Compiling+running code that uses RooUnfold

to use RooUnfold in code, make sure to #include the proper RooUnfold headers (e.g. #include "RooUnfoldBayes.h"). 
To compile code including these libraries/headers/functions, one must point g++ to not only the locations of the headers; but also the location of the RooUnfold binary, and to use the RooUnfold library/binary flag
the relvant g++ compilation therefore is

     # NAME is a bash variable that holds the name of the script, e.g. unfoldSpectra.C
     g++ -L"${ROOUNFOLD_PATH}" -I"${ROOUNFOLD_PATH}" -lRooUnfold -I"${ROOUNFOLD_PATH}/src" $NAME $(root-config --cflags --libs) -Werror -Wall -O2 -o "${NAME/%.C/}.exe"

here, the part that says 

     -L"${ROOUNFOLD_PATH}" -I"${ROOUNFOLD_PATH}" -lRooUnfold

and 

     -I"${ROOUNFOLD_PATH}/src"

points g++ to the lcoation of the RooUnfold libraries/binaries/headers for linking at compile time. I make this a bash function called "rooUnfoldCompile", and it's very similar to my rootcompile function; it differs only by the aforementioned two pieces

**********
// DETAILS REGARDING MY OWN SETUP --------------------------------------------------------
**********

in my ${CMSSW_BASE}/src do 

     svn co https://svnsrv.desy.de/public/unfolding/RooUnfold/trunk  #downloads the directory "trunk"
     mv trunk RooUnfold-svnTrunk  #rename for convenience
     cd RooUnfold-svnTrunk
     make
     make bin

in .bashrc -> rooUnfoldCompile() function as described above, defined at login

in setupenv.sh -> LD_LIBRARY_PATH and ROOUNFOLD_PATH defined, do

     source setupenv.sh  # no need to enter grid pass

now I should be able to compile code that uses RooUnfold with
   
     rooUnfoldCompile myUnfoldCode.C  # should yield myUnfoldCode.exe
