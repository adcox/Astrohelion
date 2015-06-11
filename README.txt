# Trajectory Propagation and Analysis Toolkit #

**Dependencies**
You will need the following libraries and associated headers installed before 
you can compile and run TPAT scripts
- Gnu Scientific Library
- MatIO

**To Install**
1. Extract all files from the compressed folder
2. Run 'make' to compile the source code into libraries for your system
3. Run 'sudo make install' to install the TPAT on your system

**To Run**
Scripts should include TPAT headers using the notation 
#include <tpat/tpat_cr3bp_traj.hpp>. To compile a script, use c++11 as 
the standard:

g++ -std=c++11 -ltpat -lgsl -lgslcblas -lmatio myScript.cpp -o MyScript