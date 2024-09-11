

g++ -O3 -ffast-math -I /home/russ/Scripts/git/USalign preprocessing_tests.cpp -o preprocessing_tests

## this is the compile line that I was using for the preprocessing_tests and 
## similar to the pybind11 compilation line; but it differs from USalign's 
## default Make script's compilation arguments. For sanity sake, going to revert
## back to those commands.
##g++ -O3 -ffast-math -Wall -std=c++20 -fvisibility=hidden -g -I /home/russ/Scripts/git/USalign preprocessing_tests.cpp -o preprocessing_tests #-static# -lm

