# This file will compile various tests
#
# Syntax Notes
# 	$^ gives the RHS of the :
#	$@ gives the LHS of the :
#	$< gives the first item in the dependencies list, i.e. RHS of :
#
#	This code: $(patsubst %.o,$(OBJ)/%.o, $^)
#	replaces all the object %.o files in the dependency list with a relative path using the $(OBJ) macro

############################################################
# Macros for compiling
############################################################

# Paths for files
INC := ../include
SRC := ../src
OBJ := ../obj
BIN := ../bin

# Compiler specification and flags
CXX := clang++ -std=c++11
CFLAGS += -g -W -Wall -Wextra -pedantic 
COMP := $(CXX) $(CFLAGS)

# Library names and locations
LIBS = gsl gslcblas matio
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

############################################################
## Macros to make dependency lists shorter (don't have to put)
## OBJ in front of all the files, use a macro to do the substitution
############################################################


############################################################
.PHONY: printVars

all: simEngineTest
	

############################################################
## TESTS - All executable %.out files go in the BIN directory
############################################################

calcTest: $(OBJECTS)
	$(COMP) -I $(INC) $^ $(LDFLAGS) TestCalc.cpp -o $(BIN)/$@

linMotionTest: $(OBJECTS)
	$(COMP) -I $(INC) $^ $(LDFLAGS) LinMotionTest.cpp -o $(BIN)/$@

simEngineTest: $(OBJECTS)
	$(COMP) -I $(INC) $^ $(LDFLAGS) SimulationEngine_test.cpp -o $(BIN)/$@

############################################################
## OBJECTS - All the %.o files go in the OBJ directory
############################################################

$(OBJ)/%.o: $(SRC)/%.cpp
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



