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

# Paths for files
INC := ../include
SRC := ../src
OBJ := ../obj_debug
BIN := ../bin
LIB := ../lib

# Compiler specification and flags
CXX := clang++ -std=c++11
CFLAGS += -ggdb -W -Wall -Wextra -Weffc++ -pedantic 
COMP := $(CXX) $(CFLAGS)

# Library names and locations
LIBS = gsl gslcblas matio cspice
LDFLAGS += $(foreach lib, $(LIBS),-l$(lib))

# Options that are platform dependent
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Linux)
	LDFLAGS += -L /usr/local/lib
endif

# Get JUST the filenames, no filepaths, of the source files
SRC_FILES := $(notdir $(shell find $(SRC)/*.cpp))

# Get a list of all CPP files in SRC directory
SOURCES := $(addprefix $(SRC)/,$(SRC_FILES))

# Get list of all object files by copying source file names and 
# putting the OBJ path before the name
OBJECTS := $(patsubst %.cpp,$(OBJ)/%.o, $(SRC_FILES))

# Header files that don't have associated objects; we need the compiler to
# know that objects are dependent on these and to update if changes are made to them
IMPORTANT_HEADERS := tpat.hpp tpat_ascii_output.hpp tpat_constants.hpp tpat_exceptions.hpp tpat_utilities.hpp

HEADER_DEPS := $(addprefix $(INC)/,$(IMPORTANT_HEADERS))
############################################################
## Macros to make dependency lists shorter (don't have to put)
## OBJ in front of all the files, use a macro to do the substitution
############################################################


############################################################
.PHONY: printVars

all: prop
	

############################################################
## TESTS - All executable %.out files go in the BIN directory
############################################################

famStability: $(OBJECTS)
	$(COMP) -I $(INC) $^ $(LDFLAGS) Fam_PeriStability.cpp -o $(BIN)/$@

generateLiss: $(OBJECTS)
	$(COMP) -I $(INC) $^ $(LDFLAGS) GenerateLiss.cpp -o $(BIN)/$@

prop: $(OBJECTS)
	$(COMP) -I $(INC) $^ $(LDFLAGS) PropFromScienceOrbits.cpp -o $(BIN)/$@
	
manPeri: $(OBJECTS)
	$(COMP) -I $(INC) $^ $(LDFLAGS) EML2Halo_UM_Man_Periapses.cpp -o $(BIN)/$@
	
manTubeGridding: $(OBJECTS)
	$(COMP) -I $(INC) $^ $(LDFLAGS) IceCube_ManTubeForwardProp.cpp -o $(BIN)/$@

spEffects: $(OBJECTS)
	$(COMP) -I $(INC) $^ $(LDFLAGS) LPF_LocateSP_ManyBodies.cpp -o $(BIN)/$@

############################################################
## OBJECTS - All the %.o files go in the OBJ directory
############################################################

$(OBJ)/%.o: $(SRC)/%.cpp $(HEADER_DEPS)
	$(COMP) -I $(INC) -c $< -o $@

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



