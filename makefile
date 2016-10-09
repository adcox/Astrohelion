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

UNAME_S := $(shell uname -s)

# My headers
INC := include
# System headers; these are included to not throw warnings
INC_EXTERN := include_extern
# Source files
SRC := src
# Destination directory for compiled objects; use one for optimized, second for debug versions
OBJ := obj
# OBJ := obj_debug
# Directory for compiled binaries
BIN := bin
# Location of library dependencies
LIB := lib

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

# Library names and locations
LIBS = gsl gslcblas matio cspice boost_filesystem boost_system
LDFLAGS += $(foreach lib, $(LIBS),-l$(lib))

SYS_INC_DIR := /usr/local/include/astrohelion
SYS_INC_EXTERN_DIR := /usr/local/include/astrohelion_extern

# Options that are platform dependent

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
IMPORTANT_HEADERS := Core.hpp AsciiOutput.hpp Common.hpp Exceptions.hpp Utilities.hpp

HEADER_DEPS := $(addprefix $(INC)/,$(IMPORTANT_HEADERS))

############################################################
## Other Macros
############################################################

MKDIR_P = mkdir -p

############################################################
.PHONY: directories

directories: $(OBJ)

$(OBJ):
	$(MKDIR_P) $(OBJ)

all: directories
ifeq ($(UNAME_S), Linux)
	@echo Making Linux libraries
	@make libastrohelion.a
#	@make libastrohelion.so
else ifeq ($(UNAME_S), Darwin)
	@echo Making OS X libraries
	@make libastrohelion.a
	@make libastrohelion.dylib
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
	cp -r $(INC_EXTERN) $(SYS_INC_EXTERN_DIR)
	cp $(LIB)/libastrohelion.a /usr/local/lib/
	cp $(LIB)/libastrohelion.dylib /usr/local/lib/
endif

docs:
	doxygen doxy/dox_config

libastrohelion.a: $(OBJECTS)
	ar rcs $(LIB)/$@ $^

libastrohelion.so: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_EXTERN) $^ $(LDFLAGS) -shared -o $(LIB)/$@

libastrohelion.dylib: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_EXTERN) $^ $(LDFLAGS) -shared -o $(LIB)/$@

############################################################
## OBJECTS - All the %.o files go in the OBJ directory
############################################################

$(OBJ)/%.o: $(SRC)/%.cpp $(HEADER_DEPS)
	$(COMP) -I $(INC) -isystem $(INC_EXTERN) -c $< -o $@

############################################################
## UTILITY
############################################################
clean:
	@- $(RM) $(OBJ)/*.o

cleandist: clean

nuke:
	@- $(RM) $(OBJ)/*.o
	@- $(RM) $(LIB)/libastrohelion.*

printVars:
	$(info $(SRC_FILES))
	$(info $(SOURCES))
	$(info $(OBJECTS))




