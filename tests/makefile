# This file will compile various tests
#
# Syntax Notes
# 	$^ gives the RHS of the :
#	$@ gives the LHS of the :
#	$< gives the first item in the dependencies list, i.e. RHS of :
#
#	This code: $(patsubst %.o,$(OBJ)/%.o, $^)
#	replaces all the object %.o files in the dependency list with a relative path using the $(OBJ) macro

# Instructions for creating libraries:
# <http://www.iram.fr/~roche/code/c++/AddNumbers.html>

############################################################
# Macros for compiling
############################################################

UNAME_S := $(shell uname -s)

# My headers
INC := ../include
# System headers; these are included to not throw warnings
INC_SYS  := ../include_extern
# Source files
SRC := ../src
# Destination directory for compiled objects; use one for optimized, second for debug versions
# OBJ := ../obj
OBJ := ../obj_debug
# Directory for compiled binaries
BIN := ../bin
# Location of library dependencies
LIB := ../lib

# Compiler specification and flags
ifeq ($(UNAME_S), Darwin)
	CXX := g++-6

	# This should make it easier to read symbols in stack trace (hypothetically)
	LDFLAGS += -Wl,-export_dynamic
else
	
	ifeq ($(UNAME_S), Linux)
		LDFLAGS += -L/usr/local/lib
	endif

	CXX := g++
endif

CXX += -std=c++11 -fopenmp
# CXX := g++-5 -std=c++11
CFLAGS += -ggdb -Wall -Wextra -Weffc++ -Wdisabled-optimization -Wold-style-cast -Wimport -Wmissing-declarations -Wmissing-field-initializers -pedantic
# CFLAGS += -O3 -Wall -Wextra -Weffc++ -Wdisabled-optimization -Wold-style-cast -Wimport -Wmissing-declarations -Wmissing-field-initializers -pedantic
COMP := $(CXX) $(CFLAGS)

# Library names and locations
LIBS = gsl gslcblas matio cspice boost_filesystem boost_system boost_unit_test_framework boost_test_exec_monitor
LDFLAGS += $(foreach lib, $(LIBS),-l$(lib))

# Get JUST the filenames, no filepaths, of the source files
SRC_FILES := $(notdir $(shell find $(SRC)/*.cpp))

# Get a list of all CPP files in SRC directory
SOURCES := $(addprefix $(SRC)/,$(SRC_FILES))

# Get list of all object files by copying source file names and 
# putting the OBJ path before the name
OBJECTS := $(patsubst %.cpp,$(OBJ)/%.o, $(SRC_FILES))

# Header files that don't have associated objects; we need the compiler to
# know that objects are dependent on these and to update if changes are made to them
IMPORTANT_HEADERS := AsciiOutput.hpp Common.hpp Utilities.hpp Engine.hpp

HEADER_DEPS := $(addprefix $(INC)/,$(IMPORTANT_HEADERS))

############################################################
## Other Macros
############################################################

MKDIR_P = mkdir -p

############################################################
.PHONY: directories

all: directories unit run

directories: $(OBJ)

unit: arcDataTest

run:
	$(BIN)/arcDataTest

############################################################
## TESTS - All executable %.out files go in the BIN directory
############################################################
# Order is important here:
# 1) Compiler and compiler flags
# 2) Includes
# 3) Compiled objects
# 4) Source / Executable
# 5) Linked libraries
# 6) Output

arcDataTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^  unittest_arcDataOps.cpp $(LDFLAGS) -o $(BIN)/$@

# arcDataTest: $(OBJECTS)
# 	$(COMP) -I $(INC) -isystem $(INC_SYS) $^  test_arcDataOps.cpp $(LDFLAGS) -o $(BIN)/$@

bcSysTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_bcr4bprSysData.cpp $(LDFLAGS) -o $(BIN)/$@

bcSysDefTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ unittest_bc4bp_def.cpp $(LDFLAGS) -o $(BIN)/$@

boostTest: 
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_boost.cpp $(LDFLAGS) -o $(BIN)/$@

calcTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_calc.cpp $(LDFLAGS) -o $(BIN)/$@

correctorTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_corrector.cpp $(LDFLAGS) -o $(BIN)/$@

conicTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_2bp.cpp $(LDFLAGS) -o $(BIN)/$@

famGen: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_famGeneration.cpp $(LDFLAGS) -o $(BIN)/$@

familyOps: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_cr3bpFamily.cpp $(LDFLAGS) -o $(BIN)/$@

eigenTest: 
	$(COMP) -I $(INC) -isystem $(INC_SYS) test_eigen.cpp $(LDFLAGS) -o $(BIN)/$@

lambertTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_lambert.cpp $(LDFLAGS) -o $(BIN)/$@

linMotionTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_linMotion.cpp $(LDFLAGS) -o $(BIN)/$@

lowThrustTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_cr3bp_lowThrust.cpp $(LDFLAGS) -o $(BIN)/$@

manifolds: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_manifolds.cpp $(LDFLAGS) -o $(BIN)/$@

matioTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_matio.cpp -lmatio -o $(BIN)/$@

manyArcTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_makeManyArcs.cpp $(LDFLAGS) -o $(BIN)/$@

multShootCons: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_multShootCons.cpp $(LDFLAGS) -o $(BIN)/$@

nodesetTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_nodeset.cpp $(LDFLAGS) -o $(BIN)/$@

qhManifoldTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_qhManifolds.cpp $(LDFLAGS) -o $(BIN)/$@

simEngineTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_simulationEngine.cpp $(LDFLAGS) -o $(BIN)/$@

sortEigTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_sortEigs.cpp $(LDFLAGS) -o $(BIN)/$@

sysSwitchTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_sysSwitch.cpp $(LDFLAGS) -o $(BIN)/$@	

trajTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_traj.cpp $(LDFLAGS) -o $(BIN)/$@

utilityTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_utilities.cpp $(LDFLAGS) -o $(BIN)/$@

temp: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ temp_testECIConversion.cpp $(LDFLAGS) -o $(BIN)/$@	

qh19: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ test_qh19.cpp $(LDFLAGS) -o $(BIN)/$@	
############################################################
## OBJECTS - All the %.o files go in the OBJ directory
############################################################

$(OBJ):
	$(MKDIR_P) $(OBJ)

$(OBJ)/%.o: $(SRC)/%.cpp $(HEADER_DEPS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) -c $< -o $@

############################################################
## UTILITY
############################################################
clean:
	@- $(RM) $(OBJ)/*.o

cleandist: clean

nuke:
	@- $(RM) $(OBJ)/*.o
	@- $(RM) $(BIN)/*.out

printVars:
	$(info $(SRC_FILES))
	$(info $(SOURCES))
	$(info $(OBJECTS))

list:
	@- echo "You may choose from the following options:"
	@- echo "  arcDataTest - test basic linkable structure, adding and deleting nodes, etc."
	@- echo "  bcSysTest - test BC4BP system data, saving, loading, etc."
	@- echo "  calcTest - test calculation functions"
	@- echo "  correctorTest - test multiple shooting corrector"
	@- echo "  familyOps - test the ability to locate orbits within a family"
	@- echo "  eigenTest - Test the Eigen library"
	@- echo "  lambertTest - Test the Lambert Arc Engine"
	@- echo "  linMotionTest - Test the linear motion generator"
	@- echo "  lowThrustTest - Test the low-thrust CR3BP EOM propagation"
	@- echo "  manifolds - Test manifold generation"
	@- echo "  matioTest - Test the MatIO library functions"
	@- echo "  manyArcTest - Test the capability of the simulation engine to generate many arcs in sequence"
	@- echo "  multShootCons - Test constraints in multiple shooting algorithms"
	@- echo "  nodesetTest - Test the nodeset object and functions"
	@- echo "  qhmanifoldTest - Run a script that propagates quasi-halo manifolds to test the full functionality of a real application"
	@- echo "  simEngineTest - Test the simulation engine"
	@- echo "  sortEigTest - Test the eigenvalue sorting algorithm"
	@- echo "  sysSwitchTest - Test functions that convert trajectories between systems and models"
	@- echo "  trajTest - Test the trajectory objects"
	@- echo "  utilityTest - Test the utility function library"
	@- echo "  temp - test some temporary file"



