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
INC = ../include
SRC = ../src
OBJ = ../obj
BIN = ../bin

# COMP = clang++ -W -Wall -Wextra -pedantic -std=c++0x
COMP = clang++ -W -Wall -Wextra -pedantic -std=c++11 -g
GSL_FLAGS = -lgsl -lgslcblas
MATIO_FLAGS = -lmatio

# Options that are platform dependent
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Linux)
	MATIO_FLAGS += -L /usr/local/lib
endif
############################################################
## Macros to make dependency lists shorter (don't have to put)
## OBJ in front of all the files, use a macro to do the substitution
############################################################

# Base objects and stand-alone classes
BASE_OBJ = bodyData.o calculations.o constraint.o event.o matrix.o nodeset.o sysData.o\
trajectory.o utilities.o

# Objects derived from base objects
DERIV_OBJ = bcr4bpNodeset.o bcr4bpSysData.o bcr4bpTraj.o cr3bpNodeset.o cr3bpSysData.o\
cr3bpTraj.o

# Heavy lifters
ENGINES = corrector.o simEngine.o

# List of all possible dependencies
ALL = $(BASE_OBJ) $(DERIV_OBJ) $(ENGINES)

# dependencies for simEngineTest
SIM_ENGINE_DEP = $(ALL)

# dependencies for sysDataTest
SYS_DATA_DEP = sysData.o cr3bpSysData.o bodyData.o

# dependencies for bcr4bprSysDataTest
BCR_SYS_DATA_DEP = sysData.o bcr4bpSysData.o bodyData.o

# dependencies for trajTest
TRAJ_TEST_DEP = bodyData.o matrix.o sysData.o cr3bpSysData.o cr3bpTraj.o trajectory.o

# dependencies for nodesetTest
NODESET_TEST_DEP = $(ALL)

# Where MATLAB is located on this system
MAT_ROOT = /Applications/MATLAB_R2014a.app

############################################################

all: bodyDataTest

############################################################
## TESTS - All executable %.out files go in the BIN directory
############################################################

bcr4bprSysDataTest: $(patsubst %.o,$(OBJ)/%.o, $(BCR_SYS_DATA_DEP))
	$(COMP) -I $(INC) $^ BCR4BPRSysData_test.cpp -o $(BIN)/$@
	
bodyDataTest: $(OBJ)/bodyData.o
	$(COMP) -I $(INC) $(SRC)/adtk_body_data.cpp BodyData_test.cpp -o $(BIN)/$@

cr3bpTest: $(OBJ)/bodyData.o $(OBJ)/calculations.o $(OBJ)/matrix.o
	$(COMP) -I $(INC) $^ $(GSL_FLAGS) CR3BP_test.cpp -o $(BIN)/$@

matrixTest: $(OBJ)/matrix.o
	$(COMP) -I $(INC) $^ $(GSL_FLAGS) Matrix_test.cpp -o $(BIN)/$@

nodesetTest: $(patsubst %.o,$(OBJ)/%.o, $(NODESET_TEST_DEP))
	$(COMP) -I $(INC) $^ $(GSL_FLAGS) $(MATIO_FLAGS) Nodeset_test.cpp -o $(BIN)/$@

precisionTest:
	$(COMP) NumericalPrecision_test.cpp -o $(BIN)/$@

simEngineTest: $(patsubst %.o,$(OBJ)/%.o, $(SIM_ENGINE_DEP))
	$(COMP) -I $(INC) $^ $(GSL_FLAGS) $(MATIO_FLAGS) SimulationEngine_test.cpp -o $(BIN)/$@

sysDataTest: $(patsubst %.o,$(OBJ)/%.o, $(SYS_DATA_DEP))
	$(COMP) -I $(INC) $^ CR3BPSysData_Test.cpp -o $(BIN)/$@

trajTest: $(patsubst %.o,$(OBJ)/%.o, $(TRAJ_TEST_DEP))
	$(COMP) -I $(INC) $^ $(GSL_FLAGS) $(MATIO_FLAGS) Trajectory_test.cpp -o $(BIN)/$@

matIOTest:
	$(COMP) $(MATIO_FLAGS) MatIO_test.cpp -o $(BIN)/$@

utilityTest: $(OBJ)/utilities.o
	$(COMP) -I $(INC) $^ UtilitiesTest.cpp -o $(BIN)/$@

############################################################
## OBJECTS - All the %.o files go in the OBJ directory
############################################################

$(OBJ)/bcr4bpNodeset.o: $(SRC)/adtk_bcr4bpr_nodeset.cpp $(INC)/adtk_bcr4bpr_nodeset.hpp
	$(COMP) -I $(INC) -c $< -o $@

$(OBJ)/bcr4bpSysData.o: $(SRC)/adtk_bcr4bpr_sys_data.cpp $(INC)/adtk_bcr4bpr_sys_data.hpp
	$(COMP) -I $(INC) -c $< -o $@

$(OBJ)/bcr4bpTraj.o: $(SRC)/adtk_bcr4bpr_traj.cpp $(INC)/adtk_bcr4bpr_traj.hpp
	$(COMP) -I $(INC) -c $< -o $@

$(OBJ)/bodyData.o: $(SRC)/adtk_body_data.cpp $(INC)/adtk_body_data.hpp
	$(COMP) -I $(INC) -c $< -o $@

$(OBJ)/calculations.o: $(SRC)/adtk_calculations.cpp $(INC)/adtk_calculations.hpp
	$(COMP) -I $(INC) -c $< -o $@

$(OBJ)/corrector.o: $(SRC)/adtk_correction_engine.cpp $(INC)/adtk_correction_engine.hpp
	$(COMP) -I $(INC) -c $< -o $@

$(OBJ)/constraint.o: $(SRC)/adtk_constraint.cpp $(INC)/adtk_constraint.hpp
	$(COMP) -I $(INC) -c $< -o $@

$(OBJ)/cr3bpNodeset.o: $(SRC)/adtk_cr3bp_nodeset.cpp $(INC)/adtk_cr3bp_nodeset.hpp
	$(COMP) -I $(INC) -c $< -o $@

$(OBJ)/cr3bpSysData.o: $(SRC)/adtk_cr3bp_sys_data.cpp $(INC)/adtk_cr3bp_sys_data.hpp
	$(COMP) -I $(INC) -c $< -o $@

$(OBJ)/cr3bpTraj.o: $(SRC)/adtk_cr3bp_traj.cpp $(INC)/adtk_cr3bp_traj.hpp
	$(COMP) -I $(INC) -c $< -o $@

$(OBJ)/event.o: $(SRC)/adtk_event.cpp $(INC)/adtk_event.hpp
	$(COMP) -I $(INC) -c $< -o $@

$(OBJ)/matrix.o: $(SRC)/adtk_matrix.cpp $(INC)/adtk_matrix.hpp
	$(COMP) -I $(INC) -c $< -o $@

$(OBJ)/nodeset.o: $(SRC)/adtk_nodeset.cpp $(INC)/adtk_nodeset.hpp
	$(COMP) -I $(INC) -c $< -o $@

$(OBJ)/simEngine.o: $(SRC)/adtk_simulation_engine.cpp $(INC)/adtk_simulation_engine.hpp
	$(COMP) -I $(INC) -c $< -o $@

$(OBJ)/sysData.o: $(SRC)/adtk_sys_data.cpp $(INC)/adtk_sys_data.hpp
	$(COMP) -I $(INC) -c $< -o $@

$(OBJ)/trajectory.o: $(SRC)/adtk_trajectory.cpp $(INC)/adtk_trajectory.hpp
	$(COMP) -I $(INC) -c $< -o $@

$(OBJ)/utilities.o: $(SRC)/adtk_utilities.cpp $(INC)/adtk_utilities.hpp
	$(COMP) -I $(INC) -c $< -o $@

## UTILITY
clean:
	rm $(OBJ)/*.o

nuke:
	rm $(OBJ)/*.o
	rm $(BIN)/*.out
