# Makefile for the entire library (not for individual tests)
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
INC := include
# System headers; these are included to not throw warnings
INC_SYS  := include_extern
# Source files
SRC := src
# Destination directory for compiled objects; use one for optimized, second for debug versions
# OBJ := obj
OBJ := obj_debug
# Directory for compiled binaries
BIN := bin
# Location of library dependencies
LIB := lib

# Compiler specification and flags
# CXX := clang++ -std=c++11
CXX := g++-6 -std=c++11 -fopenmp
CFLAGS += -ggdb -W -Wall -Wextra -Weffc++ -pedantic 
# CFLAGS += -O3 -W -Wall -Wextra -Weffc++ -pedantic
COMP := $(CXX) $(CFLAGS)

# Library names and locations
LIBS = gsl gslcblas matio cspice boost_filesystem boost_system
LDFLAGS += $(foreach lib, $(LIBS),-l$(lib))

SYS_INC_DIR := /usr/local/include/tpat

# Options that are platform dependent
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Linux)
	LDFLAGS += -L /usr/local/lib
else ifeq ($(UNAME_S), Darwin)
	# Do special things for OS X
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

all:
ifeq ($(UNAME_S), Linux)
	@echo Making Linux libraries
	@make libtpat.a
#	@make libtpat.so
else ifeq ($(UNAME_S), Darwin)
	@echo Making OS X libraries
	@make libtpat.a
	@make libtpat.dylib
endif

check:
	make -C tests/
	bin/calcTest
	bin/linMotionTest
	bin/matrixTest
	bin/utilityTest
	bin/sysSwitchTest
	bin/nodesetTest
	bin/simEngineTest

install:
ifeq ($(UNAME_S), Linux)
	@echo Installing Linux libraries and headers
	# Do stuff to install linux libraries
else ifeq ($(UNAME_S), Darwin)
	@echo Installing OS X libraries and headers
	@if [ ! -d $(SYS_INC_DIR) ]; then mkdir $(SYS_INC_DIR); fi
	cp -r $(INC)/* $(SYS_INC_DIR)
	cp $(LIB)/libtpat.a /usr/local/lib/
	cp $(LIB)/libtpat.dylib /usr/local/lib/
endif

libtpat.a: $(OBJECTS)
	ar rcs $(LIB)/$@ $^

libtpat.so: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) -shared -o $(LIB)/$@

libtpat.dylib: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_SYS) $^ $(LDFLAGS) -shared -o $(LIB)/$@

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
	@- $(RM) $(LIB)/libtpat.*

printVars:
	$(info $(SRC_FILES))
	$(info $(SOURCES))
	$(info $(OBJECTS))




