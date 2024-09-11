

g++ -O3 -ffast-math -I "$(dirname "$PWD")" preprocessing_tests.cpp -o preprocessing_tests 
./preprocessing_tests data/2pel_chainA.pdb 5 ground_truths/2pel_chainA_closeK5
./preprocessing_tests data/2pel_chainA.pdb 2 ground_truths/2pel_chainA
./preprocessing_tests data/3cna_chainA.pdb 5 ground_truths/3cna_chainA_closeK5

#g++ -O3 -ffast-math -I /home/russ/Scripts/git/pylelab_USalign preprocessing_tests.cpp -o preprocessing_tests && ./preprocessing_tests data/2pel_chainA.pdb 5 pyle_updated && rm preprocessing_tests
#
#g++ -O3 -ffast-math -I /home/russ/Apps/USalign preprocessing_tests.cpp -o preprocessing_tests && ./preprocessing_tests data/2pel_chainA.pdb 5 pyle_old && rm preprocessing_tests

## this is the compile line that I was using for the preprocessing_tests and 
## similar to the pybind11 compilation line; but it differs from USalign's 
## default Make script's compilation arguments. For sanity sake, going to revert
## back to those commands.
##g++ -O3 -ffast-math -Wall -std=c++20 -fvisibility=hidden -g -I "$(dirname "$PWD")" preprocessing_tests.cpp -o preprocessing_tests #-static# -lm

