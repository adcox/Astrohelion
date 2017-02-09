UNAME_S := $(shell uname -s)

# My headers
INC := ../../include
# System headers; these are included to not throw warnings
INC_EXT  := ../../include_extern
# Source files
SRC := ../../src
# Directory for object files
OBJ := ../../obj_debug

# Compiler specification and flags
ifeq ($(UNAME_S), Darwin)
	CXX := g++-6
else
	
	ifeq ($(UNAME_S), Linux)
		LDFLAGS += -L/usr/local/lib
	endif

	CXX := g++-5
endif

CXX += -std=c++11 -fopenmp
# CXX := g++-5 -std=c++11
CFLAGS += -Wall -Wextra -Weffc++ -Wdisabled-optimization -Wold-style-cast -Wimport -Wmissing-declarations -Wmissing-field-initializers -pedantic
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

MKDIR_P = mkdir -p

############################################################
## MACRO builds
############################################################

.PHONY: directories

all: directories unit run

directories: $(OBJ)

unit: arcDataOps bcSysDef simEngine

run:
	./arcDataOps
	./bcSysDef
	./simEngine
	
############################################################
## UNIT TESTS
############################################################

arcDataOps: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_EXT) $^ ../unittest_arcDataOps.cpp $(LDFLAGS) -o $@

bcSysDef: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_EXT) $^ ../unittest_bc4bp_def.cpp $(LDFLAGS) -o $@

simEngine: $(OBJECTS)
	$(COMP) -I $(INC) -isystem $(INC_EXT) $^ ../unittest_simEngine.cpp $(LDFLAGS) -o $@

############################################################
## OBJECTS - All the %.o files go in the OBJ directory
############################################################

$(OBJ):
	$(MKDIR_P) $(OBJ)

$(OBJ)/%.o: $(SRC)/%.cpp $(HEADER_DEPS)
	$(COMP) -I $(INC) -isystem $(INC_EXT) -c $< -o $@