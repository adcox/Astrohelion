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
INSTALL_DIR := /opt/local

## My headers
INC := ../include
## System headers; these are included to not throw warnings
INC_SYS  := -isystem $(INSTALL_DIR)/include
## Source files
SRC := ../src
## Destination directory for compiled objects; use one for optimized, second for debug versions
# OBJ := ../obj
OBJ := ../obj_debug
## Directory for compiled binaries
BIN := ../bin
## Location of library dependencies
LIB := $(INSTALL_DIR)/lib

# Compiler specification and flags
ifeq ($(UNAME_S), Darwin)
	CXX := g++-6
else
	CXX := g++
endif

CXX += -std=c++11 -fopenmp
CFLAGS += -ggdb -Wall -Wextra -Weffc++ -Wdisabled-optimization -Wold-style-cast -Wimport -Wmissing-declarations -Wmissing-field-initializers -pedantic
# CFLAGS += -O3 -Wall -Wextra -Weffc++ -Wdisabled-optimization -Wold-style-cast -Wimport -Wmissing-declarations -Wmissing-field-initializers -pedantic
COMP := $(CXX) $(CFLAGS)

## Library names and locations
LIBS = gsl gslcblas matio cspice boost_filesystem boost_system boost_unit_test_framework boost_test_exec_monitor
LDFLAGS += $(foreach lib, $(LIBS),-l$(lib))
LDFLAGS += -L$(LIB)

## Get JUST the filenames, no filepaths, of the source files
SRC_FILES := $(notdir $(shell find $(SRC)/*.cpp))

## Get a list of all CPP files in SRC directory
SOURCES := $(addprefix $(SRC)/,$(SRC_FILES))

## Get list of all object files by copying source file names and 
## putting the OBJ path before the name
OBJECTS := $(patsubst %.cpp,$(OBJ)/%.o, $(SRC_FILES))

## Header files that don't have associated objects; we need the compiler to
## know that objects are dependent on these and to update if changes are made to them
IMPORTANT_HEADERS := AsciiOutput.hpp Common.hpp Utilities.hpp Engine.hpp

HEADER_DEPS := $(addprefix $(INC)/,$(IMPORTANT_HEADERS))

############################################################
## Other Macros
############################################################

MKDIR_P = mkdir -p

############################################################
.PHONY: directories

all: directories unit run

directories:
	$(MKDIR_P) $(OBJ)
	$(MKDIR_P) $(BIN)/unit

unit: arcData arcset bc4bpDef calc cr3bpFam ltFam lambertArc linMotion\
 	manifoldEngine multShootCons_bc4bp multShootCons_cr3bp\
 	multShootCons_cr3bp_lt multShootEngine simEngine utilities event

run:
	$(BIN)/unit/arcData
	$(BIN)/unit/arcset
	$(BIN)/unit/bc4bpDef
	$(BIN)/unit/calc
	$(BIN)/unit/cr3bpFam
	$(BIN)/unit/ltFam
	#$(BIN)/unit/lambertArc	# Fails right now
	$(BIN)/unit/linMotion
	$(BIN)/unit/manifoldEngine
	$(BIN)/unit/multShootCons_bc4bp
	$(BIN)/unit/multShootCons_cr3bp
	$(BIN)/unit/multShootCons_cr3bp_lt
	$(BIN)/unit/multShootEngine
	$(BIN)/unit/simEngine
	$(BIN)/unit/utilities
	$(BIN)/unit/event

############################################################
## UNIT TESTS - All executable %.out files go in the BIN directory
############################################################
# Order is important here:
# 1) Compiler and compiler flags
# 2) Includes
# 3) Compiled objects
# 4) Source / Executable
# 5) Linked libraries
# 6) Output

arcData: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^  unittest_arcDataOps.cpp $(LDFLAGS) -o $(BIN)/unit/$@

arcset: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^  unittest_arcset.cpp $(LDFLAGS) -o $(BIN)/unit/$@

bc4bpDef: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^  unittest_bc4bp_def.cpp $(LDFLAGS) -o $(BIN)/unit/$@

calc: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^  unittest_calc.cpp $(LDFLAGS) -o $(BIN)/unit/$@

cr3bpFam: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^  unittest_cr3bpFamily.cpp $(LDFLAGS) -o $(BIN)/unit/$@	

event: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^  unittest_event.cpp $(LDFLAGS) -o $(BIN)/unit/$@	

lambertArc: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^  unittest_lambertArc.cpp $(LDFLAGS) -o $(BIN)/unit/$@	

linMotion: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^  unittest_linMotion.cpp $(LDFLAGS) -o $(BIN)/unit/$@	

ltFam: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^  unittest_cr3bp-lt_family.cpp $(LDFLAGS) -o $(BIN)/unit/$@	

manifoldEngine: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^  unittest_manifold.cpp $(LDFLAGS) -o $(BIN)/unit/$@	

multShootCons_bc4bp: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^  unittest_multShootCons_bc4bp.cpp $(LDFLAGS) -o $(BIN)/unit/$@

multShootCons_cr3bp: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^  unittest_multShootCons_cr3bp.cpp $(LDFLAGS) -o $(BIN)/unit/$@

multShootCons_cr3bp_lt: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^  unittest_multShootCons_cr3bp_lt.cpp $(LDFLAGS) -o $(BIN)/unit/$@

multShootEngine: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^  unittest_multShootEngine.cpp $(LDFLAGS) -o $(BIN)/unit/$@

simEngine: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^  unittest_simEngine.cpp $(LDFLAGS) -o $(BIN)/unit/$@

utilities: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^  unittest_utilities.cpp $(LDFLAGS) -o $(BIN)/unit/$@

############################################################
## OTHER TESTS - All executable %.out files go in the BIN directory
############################################################


# arcDataTest: $(OBJECTS)
# 	$(COMP) -I $(INC) $(INC_SYS) $^  test_arcDataOps.cpp $(LDFLAGS) -o $(BIN)/$@

bcSysTest: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^ test_bcr4bprSysData.cpp $(LDFLAGS) -o $(BIN)/$@

benchmark: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^ test_benchmark.cpp $(LDFLAGS) -o $(BIN)/$@

boostTest: 
	$(COMP) -I $(INC) $(INC_SYS) $^ test_boost.cpp $(LDFLAGS) -o $(BIN)/$@

calcTest: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^ test_calc.cpp $(LDFLAGS) -o $(BIN)/$@

conicTest: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^ test_2bp.cpp $(LDFLAGS) -o $(BIN)/$@

convertLT: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^  convertLTFiles.cpp $(LDFLAGS) -o $(BIN)/$@	

famGen: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^ test_famGeneration.cpp $(LDFLAGS) -o $(BIN)/$@

eigenTest: 
	$(COMP) -I $(INC) $(INC_SYS) test_eigen.cpp $(LDFLAGS) -o $(BIN)/$@

lowThrustTest: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^ test_cr3bp_lowThrust.cpp $(LDFLAGS) -o $(BIN)/$@

matioTest: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^ test_matio.cpp $(LDFLAGS) -o $(BIN)/$@

manyArcTest: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^ test_makeManyArcs.cpp $(LDFLAGS) -o $(BIN)/$@

qhManifoldTest: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^ test_qhManifolds.cpp $(LDFLAGS) -o $(BIN)/$@

sortEigTest: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^ test_sortEigs.cpp $(LDFLAGS) -o $(BIN)/$@

sysSwitchTest: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^ test_sysSwitch.cpp $(LDFLAGS) -o $(BIN)/$@	

temp: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^ temp1.cpp $(LDFLAGS) -o $(BIN)/$@	

temp2: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^ temp2.cpp $(LDFLAGS) -o $(BIN)/$@

temp3: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^ temp3.cpp $(LDFLAGS) -o $(BIN)/$@
############################################################
## OBJECTS - All the %.o files go in the OBJ directory
############################################################

$(OBJ):
	$(MKDIR_P) $(OBJ)

$(OBJ)/%.o: $(SRC)/%.cpp $(HEADER_DEPS)
	$(COMP) -I $(INC) $(INC_SYS) -c $< -o $@

############################################################
## UTILITY
############################################################
clean:
	@- $(RM) $(OBJ)/*.o

cleandist: clean

nuke:
	@- $(RM) $(OBJ)/*.o

printVars:
	$(info $(SRC_FILES))
	$(info $(SOURCES))
	$(info $(OBJECTS))

list:
	@- echo "Unit Tests:"
	@- echo "  arcData"
	@- echo "  arcset"
	@- echo "  bc4bpDef"
	@- echo "  cr3bpFam"
	@- echo "  event"
	@- echo "  lambertArc"
	@- echo "  linMotion"
	@- echo "  ltFam"
	@- echo "  manifoldEngine"
	@- echo "  multShootCons_bc4bp"
	@- echo "  multShootCons_cr3bp"
	@- echo "  multShootCons_cr3bp_lt"
	@- echo "  multShootEngine"
	@- echo "  simEngine"
	@- echo "  utilities"
	@- echo ""
	@- echo "Other Tests:"
	@- echo "  bcSysTest - test BC4BP system data, saving, loading, etc."
	@- echo "  calcTest - test calculation functions"
	@- echo "  eigenTest - Test the Eigen library"
	@- echo "  lowThrustTest - Test the low-thrust CR3BP EOM propagation"
	@- echo "  matioTest - Test the MatIO library functions"
	@- echo "  manyArcTest - Test the capability of the simulation engine to generate many arcs in sequence"
	@- echo "  qhmanifoldTest - Run a script that propagates quasi-halo manifolds to test the full functionality of a real application"
	@- echo "  sortEigTest - Test the eigenvalue sorting algorithm"
	@- echo "  sysSwitchTest - Test functions that convert trajectories between systems and models"
	@- echo "  temp - test some temporary file"



