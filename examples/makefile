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
LIBS = gsl gslcblas matio cspice boost_filesystem boost_system
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

all: directories examples run

directories:
	$(MKDIR_P) $(OBJ)
	$(MKDIR_P) $(BIN)

examples: finiteDiff_checkAMatrix

run:
	$(BIN)/finiteDiff_checkAMatrix

############################################################
## EXAMPLE SCRIPTS
############################################################

finiteDiff_checkAMatrix: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_SYS) $^ finiteDiff_checkAMatrix.cpp $(LDFLAGS) -o $(BIN)/$@

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
	@- echo "Example Scripts:"
	@- echo "  finiteDiff_checkAMatrix - example of finite differencing to check A matrix accuracy"



