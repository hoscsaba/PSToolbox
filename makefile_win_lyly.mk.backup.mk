CXX = g++
CFLAGS = -g -std=c++11 -pedantic
#TARGETS = my_tools Gas IdealGas FrozenMixtureLiquidGas Units LWP SCP Reservoir Valve Connector Valve_with_Absorber  CoolPropGas
TARGETS = my_tools SurgeChamber Gas IdealGas FrozenMixtureLiquidGas Units LWP SCP Reservoir PSToolboxRunner PSToolboxBaseEdge Valve Connector Valve_with_Absorber EpanetReader CheckValve Pump
INC = -I"C:\Program Files (x86)\Eigen\include\eigen3"
# LINK = -lmy_tools -lpython2.7
LINK = -lmy_tools\

SRC_DIR = src
INC+=-I$(SRC_DIR)

INC_CP = -I/Users/hoscsaba/program/CoolProp/include -I/Users/hoscsaba/program/CoolProp/externals/fmtlib
LINK_CP = -L/Users/hoscsaba/program/CoolProp/build1 -lCoolProp

all:$(TARGETS)
	del libmy_tools.a
	ar rvs libmy_tools.a my_tools.o
	ranlib libmy_tools.a
	del libPSToolbox.a
	ar rvs libPSToolbox.a PSToolBoxRunner.o PSToolboxBaseEdge.o Gas.o IdealGas.o FrozenMixtureLiquidGas.o Units.o LWP.o SCP.o Reservoir.o Valve.o Connector.o Valve_with_Absorber.o CheckValve.o Pump.o EpanetReader.o 

	ranlib libPSToolbox.a

#Compile rule for each source file
my_tools: my_tools.cpp 
	$(CXX) $(INC) $(CFLAGS) my_tools.cpp -c -o my_tools.o

SurgeChamber: $(SRC_DIR)/SurgeChamber.cpp 
	$(CXX) $(INC) $(CFLAGS) $(SRC_DIR)/SurgeChamber.cpp -c -o SurgeChamber.o

PSToolboxRunner: $(SRC_DIR)/PSToolboxRunner.cpp 
	$(CXX) $(INC) $(CFLAGS) $(SRC_DIR)/PSToolboxRunner.cpp -c -o PSToolboxRunner.o

PSToolboxBaseEdge: $(SRC_DIR)/PSToolboxBaseEdge.cpp 
	$(CXX) $(INC) $(CFLAGS) $(SRC_DIR)/PSToolboxBaseEdge.cpp -c -o PSToolboxBaseEdge.o

Gas: $(SRC_DIR)/Gas.cpp 
	$(CXX) $(CFLAGS) $(SRC_DIR)/Gas.cpp -c -o Gas.o

IdealGas: $(SRC_DIR)/IdealGas.cpp 
	$(CXX) $(CFLAGS) $(SRC_DIR)/IdealGas.cpp -c -o IdealGas.o

FrozenMixtureLiquidGas: $(SRC_DIR)/FrozenMixtureLiquidGas.cpp 
	$(CXX) $(CFLAGS) $(SRC_DIR)/FrozenMixtureLiquidGas.cpp -c -o FrozenMixtureLiquidGas.o

Units: $(SRC_DIR)/Units.cpp
	$(CXX) $(CFLAGS) $(SRC_DIR)/Units.cpp -c -o Units.o

LWP: $(SRC_DIR)/LWP.cpp
	$(CXX) $(INC) $(CFLAGS) $(SRC_DIR)/LWP.cpp -c -o LWP.o

SCP: $(SRC_DIR)/SCP.cpp
	$(CXX) $(INC) $(CFLAGS) $(SRC_DIR)/SCP.cpp -c -o SCP.o

Reservoir: $(SRC_DIR)/Reservoir.cpp
	$(CXX) $(INC) $(CFLAGS) $(SRC_DIR)/Reservoir.cpp -c -o Reservoir.o
# 	$(CXX) $(INC) $(LINK) $(CFLAGS) Reservoir.cpp -c -o Reservoir.o

Valve: $(SRC_DIR)/Valve.cpp
	$(CXX) $(INC) $(CFLAGS) $(SRC_DIR)/Valve.cpp -c -o Valve.o
# 	$(CXX) $(INC) $(LINK) $(CFLAGS) Valve.cpp -c -o Valve.o

Pump: $(SRC_DIR)/Pump.cpp
	$(CXX) $(INC) $(CFLAGS) $(SRC_DIR)/Pump.cpp -c -o Pump.o

EpanetReader: $(SRC_DIR)/EpanetReader.cpp
	$(CXX) $(INC) $(CFLAGS) $(SRC_DIR)/EpanetReader.cpp -c -o EpanetReader.o

CheckValve: $(SRC_DIR)/CheckValve.cpp
	$(CXX) $(INC) $(CFLAGS) $(SRC_DIR)/CheckValve.cpp -c -o CheckValve.o

Valve_with_Absorber: $(SRC_DIR)/Valve_with_absorber.cpp
	$(CXX) $(INC) $(CFLAGS) $(SRC_DIR)/Valve_with_Absorber.cpp -c -o Valve_with_Absorber.o

Connector: $(SRC_DIR)/Connector.cpp
	$(CXX) $(INC) $(CFLAGS) $(SRC_DIR)/Connector.cpp -c -o Connector.o

CoolPropGas: $(SRC_DIR)/CoolPropGas.cpp
	$(CXX) $(INC_CP) $(LINK_CP) $(CFLAGS) $(SRC_DIR)/CoolPropGas.cpp -c -o CoolPropGas.o

CoolPropMixture: $(SRC_DIR)/CoolPropMixture.cpp
	$(CXX) $(INC_CP) $(LINK_CP) $(CFLAGS) $(SRC_DIR)/CoolPropMixture.cpp -c -o CoolPropMixture.o

CoolPropHA: $(SRC_DIR)/CoolPropHA.cpp
	$(CXX) $(INC_CP) $(LINK_CP) $(CFLAGS) $(SRC_DIR)/CoolPropHA.cpp -c -o CoolPropHA.o


clean:
	$(RM) *.o *.a