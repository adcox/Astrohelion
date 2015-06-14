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

# Paths for files
INC := include
SRC := src
OBJ := obj
BIN := bin
LIB := lib

# Compiler specification and flags
CXX := clang++ -std=c++11
CFLAGS += -W -Wall -Wextra -pedantic -O3
COMP := $(CXX) $(CFLAGS)

# Library names and locations
LIBS = gsl gslcblas matio cspice
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
IMPORTANT_HEADERS := tpat.hpp tpat_ascii_output.hpp tpat_constants.hpp tpat_exceptions.hpp

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

install:
ifeq ($(UNAME_S), Linux)
	@echo Installing Linux libraries and headers
	# Do stuff to install linux libraries
else ifeq ($(UNAME_S), Darwin)
	@echo Installing OS X libraries and headers
	@if [ ! -d $(SYS_INC_DIR) ]; then mkdir $(SYS_INC_DIR); fi
	cp $(INC)/*.hpp $(SYS_INC_DIR)
	cp $(LIB)/libtpat.a /usr/local/lib/
	cp $(LIB)/libtpat.dylib /usr/local/lib/
endif

libtpat.a: $(OBJECTS)
	ar rcs $(LIB)/$@ $^

libtpat.so: $(OBJECTS)
	$(COMP) -I $(INC) $^ $(LDFLAGS) -shared -o $(LIB)/$@

libtpat.dylib: $(OBJECTS)
	$(COMP) -I $(INC) $^ $(LDFLAGS) -shared -o $(LIB)/$@

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
	@- $(RM) $(LIB)/libtpat.*

printVars:
	$(info $(SRC_FILES))
	$(info $(SOURCES))
	$(info $(OBJECTS))




