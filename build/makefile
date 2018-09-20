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
# <http://www.yolinux.com/TUTORIALS/LibraryArchives-StaticAndDynamic.html>
# 
############################################################
# Macros for compiling
############################################################

UNAME_S := $(shell uname -s)
INSTALL_DIR := /opt/local

# My headers
INC := include
# System headers; these are included to not throw warnings
INC_EXTERN := $(INSTALL_DIR)/include
# Source files
SRC := src
# Destination directory for compiled objects; use one for optimized, second for debug versions
OBJ := obj
# OBJ := obj_debug
# Directory for compiled binaries
BIN := bin/astrohelion
# Location of library dependencies
LIB := $(INSTALL_DIR)/lib

# Compiler specification and flags
ifeq ($(UNAME_S), Darwin)
	CXX := g++-6
else
	CXX := g++
endif

CXX += -std=c++11 -fopenmp
# CFLAGS += -ggdb -Wall -Wextra -Weffc++ -Wdisabled-optimization -Wold-style-cast -Wimport -Wmissing-declarations -Wmissing-field-initializers -pedantic
CFLAGS += -O3 -Wall -Wextra -Weffc++ -Wdisabled-optimization -Wold-style-cast -Wimport -Wmissing-declarations -Wmissing-field-initializers -pedantic -Wno-deprecated-declarations
COMP := $(CXX) $(CFLAGS)
# SO_FLAGS are set to blank by default and only used during linux shared library creation
SO_FLAGS := 

# Library names and locations
LIBS = gsl gslcblas matio cspice
LDFLAGS += $(foreach lib, $(LIBS),-l$(lib))
LDFLAGS += -L$(LIB)

SYS_INC_DIR := $(INSTALL_DIR)/include/astrohelion
SYS_INC_EXTERN_DIR := $(INSTALL_DIR)/include

# Get JUST the filenames, no filepaths, of the source files
SRC_FILES := $(notdir $(shell find $(SRC)/*.cpp))

# Get a list of all CPP files in SRC directory
SOURCES := $(addprefix $(SRC)/,$(SRC_FILES))

# Get list of all object files by copying source file names and 
# putting the OBJ path before the name
OBJECTS := $(patsubst %.cpp,$(OBJ)/%.o, $(SRC_FILES))

# Header files that don't have associated objects; we need the compiler to
# know that objects are dependent on these and to update if changes are made to them
IMPORTANT_HEADERS := Core.hpp AsciiOutput.hpp Common.hpp Utilities.hpp
# IMPORTANT_HEADERS := AsciiOutput.hpp Common.hpp Utilities.hpp

HEADER_DEPS := $(addprefix $(INC)/,$(IMPORTANT_HEADERS))

############################################################
## Other Macros
############################################################

MKDIR_P = mkdir -p

############################################################
.PHONY: directories

all: directories static clean dynamic

directories: $(OBJ)

static:
ifeq ($(UNAME_S), Linux)
	@echo Making Linux static libraries
	@make libastrohelion.a
else ifeq ($(UNAME_S), Darwin)
	@echo Making Apple OS static libraries
	@make libastrohelion.a
endif

dynamic:
ifeq ($(UNAME_S), Linux)
	@echo Making Linux dynamic libraries
	@make libastrohelion.so.1
else ifeq ($(UNAME_S), Darwin)
	@echo Making Apple OS dynamic libraries
	@make libastrohelion.dylib
endif

check:
	make -C tests/
	$(BIN)/calcTest
	$(BIN)/linMotionTest
	$(BIN)/matrixTest
	$(BIN)/utilityTest
	$(BIN)/sysSwitchTest
	$(BIN)/nodesetTest
	$(BIN)/simEngineTest

# Install Astrohelion in the INSTALL_DIR
#	Step 1: Create the include directory if it doesnot exist
#	Step 2: Copy the include files to the install include directory
#	Step 3: (Linux Only) Create symbolic links to the shared library.
#		Two links are created:
#			- lib*.so.1 -> lib*.so.1.0 allows run time binding to work
#			- lib*.so -> lib*.so.1.0 allows naming convention for compile
#				flag (-lastrohelion) to work
install:
	@if [ ! -d $(SYS_INC_DIR) ]; then mkdir -p $(SYS_INC_DIR); fi
	cp -r $(INC)/* $(SYS_INC_DIR)
ifeq ($(UNAME_S), Linux)
	ln -sf $(INSTALL_DIR)/lib/libastrohelion.so.1.0 $(INSTALL_DIR)/lib/libastrohelion.so.1
	ln -sf $(INSTALL_DIR)/lib/libastrohelion.so.1.0 $(INSTALL_DIR)/lib/libastrohelion.so
endif

docs:
	doxygen doxy/dox_config

# Create the static library (Linux and Apple OS)
libastrohelion.a: $(OBJECTS)
	ar rcs $(LIB)/$@ $^

# Create the dynamic library (Linux)
# 	Step 1: Create the object code using the SO_FLAGS set to -fPIC.
# 		This tells the compiler to output position independent code,
# 		which is required for shared libraries
# 	Step 2: Create the shared library. We append .0 to the library
# 		name, so it will be libastrohelion.so.1.0. Other, similar, 
# 		names are used for symlinks in the `install` step
libastrohelion.so.1: SO_FLAGS = -fPIC
libastrohelion.so.1: $(OBJECTS)
	$(COMP) -shared -Wl,-soname,$@ -o $(LIB)/$@.0 $^ $(LDFLAGS)

# Create the dynamic library (Apple OS)
libastrohelion.dylib: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_EXTERN) $^ $(LDFLAGS) -shared -o $(LIB)/$@

############################################################
## OBJECTS - All the %.o files go in the OBJ directory
############################################################

$(OBJ):
	$(MKDIR_P) $(OBJ)

$(OBJ)/%.o: $(SRC)/%.cpp $(HEADER_DEPS)
	$(COMP) $(SO_FLAGS) -I $(INC) -isystem $(INC_EXTERN) -c $< -o $@

############################################################
## UTILITY
############################################################
clean:
	@- $(RM) $(OBJ)/*.o

nuke:
	@- $(RM) $(OBJ)/*.o
	@- $(RM) $(LIB)/libastrohelion.*

printVars:
	$(info $(SRC_FILES))
	$(info $(SOURCES))
	$(info $(OBJECTS))




