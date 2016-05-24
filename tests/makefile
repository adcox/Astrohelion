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
# CXX := clang++ -std=c++11
CXX := g++-5 -std=c++11 -fopenmp
# CXX := g++-5 -std=c++11
CFLAGS += -ggdb -W -Wall -Wextra -Weffc++ -pedantic
# CFLAGS += -O3 -W -Wall -Wextra -Weffc++ -pedantic
COMP := $(CXX) $(CFLAGS)

# Library names and locations
LIBS = gsl gslcblas matio cspice boost_filesystem boost_system
LDFLAGS += $(foreach lib, $(LIBS),-l$(lib))

# Options that are platform dependent
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Linux)
	LDFLAGS += -L /usr/local/lib
endif

# This should make it easier to read symbols in stack trace (hypothetically)
LDFLAGS += -Wl,-export_dynamic

# Get JUST the filenames, no filepaths, of the source files
SRC_FILES := $(notdir $(shell find $(SRC)/*.cpp))

# Get a list of all CPP files in SRC directory
SOURCES := $(addprefix $(SRC)/,$(SRC_FILES))

# Get list of all object files by copying source file names and 
# putting the OBJ path before the name
OBJECTS := $(patsubst %.cpp,$(OBJ)/%.o, $(SRC_FILES))

# Header files that don't have associated objects; we need the compiler to
# know that objects are dependent on these and to update if changes are made to them
IMPORTANT_HEADERS := tpat_ascii_output.hpp tpat_constants.hpp tpat_utilities.hpp

HEADER_DEPS := $(addprefix $(INC)/,$(IMPORTANT_HEADERS))
############################################################
## Macros to make dependency lists shorter (don't have to put)
## OBJ in front of all the files, use a macro to do the substitution
############################################################


############################################################
.PHONY: printVars

all: arcDataTest bcSysTest calcTest correctorTest familyOps eigenTest linMotionTest matioTest multShootCons nodesetTest simEngineTest sysSwitchTest trajTest utilityTest
	

############################################################
## TESTS - All executable %.out files go in the BIN directory
############################################################

arcDataTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) test_arcDataOps.cpp -o $(BIN)/$@
	
bcSysTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) test_bcr4bprSysData.cpp -o $(BIN)/$@

boostTest:
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) test_boost.cpp -o $(BIN)/$@

calcTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) test_calc.cpp -o $(BIN)/$@

correctorTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) test_corrector.cpp -o $(BIN)/$@

familyOps: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) test_cr3bpFamily.cpp -o $(BIN)/$@

eigenTest:
	$(COMP) -I $(INC) -isystem $(INC_SYS) test_eigen.cpp -o $(BIN)/$@

linMotionTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) test_linMotion.cpp -o $(BIN)/$@

lowThrustTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) test_cr3bp_lowThrust.cpp -o $(BIN)/$@

manifolds: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) test_manifolds.cpp -o $(BIN)/$@

matioTest: $(OBJ)/tpat_utilities.o
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) -lmatio test_matio.cpp -o $(BIN)/$@

manyArcTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) test_makeManyArcs.cpp -o $(BIN)/$@

multShootCons: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) test_multShootCons.cpp -o $(BIN)/$@

nodesetTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) test_nodeset.cpp -o $(BIN)/$@

simEngineTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) test_simulationEngine.cpp -o $(BIN)/$@

sortEigTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) test_sortEigs.cpp -o $(BIN)/$@

sysSwitchTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) test_sysSwitch.cpp -o $(BIN)/$@	

trajTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) test_traj.cpp -o $(BIN)/$@

utilityTest: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) test_utilities.cpp -o $(BIN)/$@

temp: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) temp_prop.cpp -o $(BIN)/$@	

############################################################
## OBJECTS - All the %.o files go in the OBJ directory
############################################################

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
	@- echo "  linMotionTest - Test the linear motion generator"
	@- echo "  lowThrustTest - Test the low-thrust CR3BP EOM propagation"
	@- echo "  manifolds - Test manifold generation"
	@- echo "  matioTest - Test the MatIO library functions"
	@- echo "  manyArcTest - Test the capability of the simulation engine to generate many arcs in sequence"
	@- echo "  multShootCons - Test constraints in multiple shooting algorithms"
	@- echo "  nodesetTest - Test the nodeset object and functions"
	@- echo "  simEngineTest - Test the simulation engine"
	@- echo "  sortEigTest - Test the eigenvalue sorting algorithm"
	@- echo "  sysSwitchTest - Test functions that convert trajectories between systems and models"
	@- echo "  trajTest - Test the trajectory objects"
	@- echo "  utilityTest - Test the utility function library"
	@- echo "  temp - test some temporary file"



