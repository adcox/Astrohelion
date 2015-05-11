# This file will compile various tests

INC = ../include
SRC = ../src

COMP = g++

all: bodyDataTest

clean:
	rm *.out

bodyDataTest:
	$(COMP) -Wall -I $(INC) $(SRC)/adtk_body_data.cpp BodyData_test.cpp -o BodyDataTest.out

cr3bpTest:
	$(COMP) -Wall -I $(INC) $(SRC)/adtk_body_data.cpp $(SRC)/adtk_matrix.cpp $(SRC)/adtk_integrators.cpp CR3BP_test.cpp -lgsl -lgslcblas -o CR3BPTest.out

matrixTest:
	$(COMP) -Wall -I $(INC) $(SRC)/adtk_matrix.cpp Matrix_Test.cpp -lgsl -lgslcblas -o MatrixTest.out