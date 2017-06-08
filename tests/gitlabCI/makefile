UNAME_S := $(shell uname -s)

# My headers
INC := ../../include
# System headers; these are included to not throw warnings
INC_EXT  := -isystem/opt/local/include
# Source files
SRC := ../../src
# Directory for object files
OBJ := ../../obj_debug

# Compiler specification and flags
ifeq ($(UNAME_S), Darwin)
	CXX := g++-6
else
	
	ifeq ($(UNAME_S), Linux)
		LDFLAGS += -L/opt/local/lib
	endif

	CXX := g++
endif

CXX += -std=c++11 -fopenmp
# CXX := g++-5 -std=c++11
CFLAGS += -Wall -Wextra -Weffc++ -Wdisabled-optimization -Wold-style-cast -Wimport -Wmissing-declarations -Wmissing-field-initializers -pedantic
# CFLAGS += -O3 -Wall -Wextra -Weffc++ -Wdisabled-optimization -Wold-style-cast -Wimport -Wmissing-declarations -Wmissing-field-initializers -pedantic
COMP := $(CXX) $(CFLAGS)

# Library names and locations
LDFLAGS += -Wl,-rpath-link,/opt/local/lib -L/opt/local/lib
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

unit: arcData bc4bpDef cr3bpFam lambertArc linMotion manifoldEngine\
	multShootCons_bc4bp multShootCons_cr3bp multShootCons_cr3bp_lt multShootEngine\
	arcset simEngine utilities

run:
	./arcData.bin
	./bc4bpDef.bin
	#./cr3bpFam.bin		# Requires family data...
	#./lambertArc.bin	# Fails right now
	./linMotion.bin
	./manifoldEngine.bin
	./multShootCons_bc4bp.bin
	./multShootCons_cr3bp.bin
	./multShootCons_cr3bp_lt.bin
	./multShootEngine.bin
	./arcset.bin
	./simEngine.bin
	./utilities.bin
	
############################################################
## UNIT TESTS
############################################################

arcData: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_EXT) $^  ../unittest_arcDataOps.cpp $(LDFLAGS) -o $@.bin

arcset: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_EXT) $^  ../unittest_arcset.cpp $(LDFLAGS) -o $@.bin
	
bc4bpDef: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_EXT) $^  ../unittest_bc4bp_def.cpp $(LDFLAGS) -o $@.bin

cr3bpFam: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_EXT) $^  ../unittest_cr3bpFamily.cpp $(LDFLAGS) -o $@.bin	

lambertArc: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_EXT) $^  ../unittest_lambertArc.cpp $(LDFLAGS) -o $@.bin	

linMotion: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_EXT) $^  ../unittest_linMotion.cpp $(LDFLAGS) -o $@.bin	

manifoldEngine: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_EXT) $^  ../unittest_manifold.cpp $(LDFLAGS) -o $@.bin	

multShootCons_bc4bp: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_EXT) $^  ../unittest_multShootCons_bc4bp.cpp $(LDFLAGS) -o $@.bin

multShootCons_cr3bp: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_EXT) $^  ../unittest_multShootCons_cr3bp.cpp $(LDFLAGS) -o $@.bin

multShootCons_cr3bp_lt: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_EXT) $^  ../unittest_multShootCons_cr3bp_lt.cpp $(LDFLAGS) -o $@.bin

multShootEngine: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_EXT) $^  ../unittest_multShootEngine.cpp $(LDFLAGS) -o $@.bin

simEngine: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_EXT) $^  ../unittest_simEngine.cpp $(LDFLAGS) -o $@.bin

utilities: $(OBJECTS)
	$(COMP) -I $(INC) $(INC_EXT) $^  ../unittest_utilities.cpp $(LDFLAGS) -o $@.bin

############################################################
## OBJECTS - All the %.o files go in the OBJ directory
############################################################

directories:
	$(MKDIR_P) $(OBJ)
	$(MKDIR_P) data

$(OBJ)/%.o: $(SRC)/%.cpp $(HEADER_DEPS)
	$(COMP) -I $(INC) $(INC_EXT) -c $< -o $@