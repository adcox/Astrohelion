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

COMP = g++ -W -Wall -Wextra -pedantic -std=c++0x
GSL_FLAGS = -lgsl -lgslcblas

############################################################
## Macros to make dependency lists shorter (don't have to put)
## OBJ in front of all the files, use a macro to do the substitution
############################################################

# dependencies for simEngineTest
SIM_ENGINE_DEP = bodyData.o calculations.o cr3bpSysData.o matrix.o simEngine.o sysData.o trajectory.o cr3bpTraj.o
# dependencies for sysDataTest
SYS_DATA_DEP = sysData.o cr3bpSysData.o bodyData.o
# dependencies for trajTest
TRAJ_TEST_DEP = bodyData.o matrix.o sysData.o cr3bpSysData.o cr3bpTraj.o trajectory.o

############################################################

all: bodyDataTest

############################################################
## TESTS - All executable %.out files go in the BIN directory
############################################################

bodyDataTest: $(OBJ)/bodyData.o
	$(COMP) -I $(INC) $(SRC)/adtk_body_data.cpp BodyData_test.cpp -o $(BIN)/$@

cr3bpTest: $(OBJ)/bodyData.o $(OBJ)/calculations.o $(OBJ)/matrix.o
	$(COMP) -I $(INC) $^ $(GSL_FLAGS) CR3BP_test.cpp -o $(BIN)/$@

matrixTest: $(OBJ)/matrix.o
	$(COMP) -I $(INC) $^ $(GSL_FLAGS) Matrix_test.cpp -o $(BIN)/$@

precisionTest:
	$(COMP) NumericalPrecision_test.cpp -o $(BIN)/$@

simEngineTest: $(patsubst %.o,$(OBJ)/%.o, $(SIM_ENGINE_DEP))
	$(COMP) -I $(INC) $^ $(GSL_FLAGS) SimulationEngine_test.cpp -o $(BIN)/$@

sysDataTest: $(patsubst %.o,$(OBJ)/%.o, $(SYS_DATA_DEP))
	$(COMP) -I $(INC) $^ CR3BPSysData_Test.cpp -o $(BIN)/$@

trajTest: $(patsubst %.o,$(OBJ)/%.o, $(TRAJ_TEST_DEP))
	$(COMP) -I $(INC) $^ $(GSL_FLAGS) Trajectory_test.cpp -o $(BIN)/$@

############################################################
## OBJECTS - All the %.o files go in the OBJ directory
############################################################

$(OBJ)/bodyData.o: $(SRC)/adtk_body_data.cpp $(INC)/adtk_body_data.hpp
	$(COMP) -I $(INC) -c $< -o $(OBJ)/$@

$(OBJ)/cr3bpSysData.o: $(SRC)/adtk_cr3bp_sys_data.cpp $(INC)/adtk_cr3bp_sys_data.hpp
	$(COMP) -I $(INC) -c $< -o $(OBJ)/$@

$(OBJ)/cr3bpTraj.o: $(SRC)/adtk_cr3bp_traj.cpp $(INC)/adtk_cr3bp_traj.hpp
	$(COMP) -I $(INC) -c $< -o $(OBJ)/$@

$(OBJ)/calculations.o: $(SRC)/adtk_calculations.cpp $(INC)/adtk_calculations.hpp
	$(COMP) -I $(INC) -c $< -o $(OBJ)/$@ 

$(OBJ)/matrix.o: $(SRC)/adtk_matrix.cpp $(INC)/adtk_matrix.hpp
	$(COMP) -I $(INC) -c $< -o $@

$(OBJ)/simEngine.o: $(SRC)/adtk_simulation_engine.cpp $(INC)/adtk_simulation_engine.hpp
	$(COMP) -I $(INC) -c $< -o $(OBJ)/$@

$(OBJ)/sysData.o: $(SRC)/adtk_sys_data.cpp $(INC)/adtk_sys_data.hpp
	$(COMP) -I $(INC) -c $< -o $(OBJ)/$@

$(OBJ)/trajectory.o: $(SRC)/adtk_trajectory.cpp $(INC)/adtk_trajectory.hpp
	$(COMP) -I $(INC) -c $< -o $(OBJ)/$@

## UTILITY

clean:
	rm $(OBJ)/*.o

nuke:
	rm $(OBJ)/*.o
	rm $(BIN)/*.out
