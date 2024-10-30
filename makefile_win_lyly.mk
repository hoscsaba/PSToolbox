CXX = g++
CFLAGS = -g -std=c++11 -pedantic
#TARGETS = my_tools Gas IdealGas FrozenMixtureLiquidGas Units LWP SCP Reservoir Valve Connector Valve_with_Absorber  CoolPropGas
TARGETS = my_tools SurgeChamber Gas IdealGas FrozenMixtureLiquidGas Units LWP SCP Reservoir PSToolboxRunner PSToolboxBaseEdge Valve Connector Valve_with_Absorber EpanetReader CheckValve Pump
INC = -I"C:\Program Files (x86)\Eigen\include\eigen3"
# LINK = -lmy_tools -lpython2.7
LINK = -lmy_tools

INC_CP = -I/Users/hoscsaba/program/CoolProp/include -I/Users/hoscsaba/program/CoolProp/externals/fmtlib
LINK_CP = -L/Users/hoscsaba/program/CoolProp/build1 -lCoolProp

all:$(TARGETS)
	del libmy_tools.a
	ar rvs libmy_tools.a my_tools.o
	ranlib libmy_tools.a
	del libPSToolbox.a
	ar rvs libPSToolbox.a PSToolBoxRunner.o PSToolboxBaseEdge.o Gas.o IdealGas.o FrozenMixtureLiquidGas.o Units.o LWP.o SCP.o Reservoir.o Valve.o Connector.o Valve_with_Absorber.o CheckValve.o Pump.o EpanetReader.o SurgeChamber.o

	ranlib libPSToolbox.a

#Compile rule for each source file
my_tools: my_tools.cpp 
	$(CXX) $(INC) $(CFLAGS) my_tools.cpp -c -o my_tools.o

SurgeChamber: SurgeChamber.cpp 
	$(CXX) $(INC) $(CFLAGS) SurgeChamber.cpp -c -o SurgeChamber.o

PSToolboxRunner: PSToolboxRunner.cpp 
	$(CXX) $(INC) $(CFLAGS) PSToolboxRunner.cpp -c -o PSToolboxRunner.o

PSToolboxBaseEdge: PSToolboxBaseEdge.cpp 
	$(CXX) $(INC) $(CFLAGS) PSToolboxBaseEdge.cpp -c -o PSToolboxBaseEdge.o

Gas: Gas.cpp 
	$(CXX) $(CFLAGS) Gas.cpp -c -o Gas.o

IdealGas: IdealGas.cpp 
	$(CXX) $(CFLAGS) IdealGas.cpp -c -o IdealGas.o

FrozenMixtureLiquidGas: FrozenMixtureLiquidGas.cpp 
	$(CXX) $(CFLAGS) FrozenMixtureLiquidGas.cpp -c -o FrozenMixtureLiquidGas.o

Units: Units.cpp
	$(CXX) $(CFLAGS) Units.cpp -c -o Units.o

LWP: LWP.cpp
	$(CXX) $(INC) $(CFLAGS) LWP.cpp -c -o LWP.o

SCP: SCP.cpp
	$(CXX) $(INC) $(CFLAGS) SCP.cpp -c -o SCP.o

Reservoir: Reservoir.cpp
	$(CXX) $(INC) $(CFLAGS) Reservoir.cpp -c -o Reservoir.o
# 	$(CXX) $(INC) $(LINK) $(CFLAGS) Reservoir.cpp -c -o Reservoir.o

Valve: Valve.cpp
	$(CXX) $(INC) $(CFLAGS) Valve.cpp -c -o Valve.o
# 	$(CXX) $(INC) $(LINK) $(CFLAGS) Valve.cpp -c -o Valve.o

Pump: Pump.cpp
	$(CXX) $(INC) $(CFLAGS) Pump.cpp -c -o Pump.o

EpanetReader: EpanetReader.cpp
	$(CXX) $(INC) $(CFLAGS) EpanetReader.cpp -c -o EpanetReader.o

CheckValve: CheckValve.cpp
	$(CXX) $(INC) $(CFLAGS) CheckValve.cpp -c -o CheckValve.o

Valve_with_Absorber: Valve_with_absorber.cpp
	$(CXX) $(INC) $(CFLAGS) Valve_with_Absorber.cpp -c -o Valve_with_Absorber.o

Connector: Connector.cpp
	$(CXX) $(INC) $(CFLAGS) Connector.cpp -c -o Connector.o

CoolPropGas: CoolPropGas.cpp
	$(CXX) $(INC_CP) $(LINK_CP) $(CFLAGS) CoolPropGas.cpp -c -o CoolPropGas.o

CoolPropMixture: CoolPropMixture.cpp
	$(CXX) $(INC_CP) $(LINK_CP) $(CFLAGS) CoolPropMixture.cpp -c -o CoolPropMixture.o

CoolPropHA: CoolPropHA.cpp
	$(CXX) $(INC_CP) $(LINK_CP) $(CFLAGS) CoolPropHA.cpp -c -o CoolPropHA.o


clean:
	$(RM) *.o *.a