CXX = g++
CFLAGS = -g -std=c++11 -pedantic
#TARGETS = my_tools Gas IdealGas FrozenMixtureLiquidGas Units LWP SCP Reservoir Valve Connector Valve_with_Absorber  CoolPropGas
TARGETS = my_tools Gas IdealGas FrozenMixtureLiquidGas Units LWP SCP Reservoir Valve Connector Valve_with_Absorber 
INC = -IC:/ProgramData/chocolatey/lib/eigen/include
# LINK = -lmy_tools -lpython2.7
LINK = -lmy_tools

INC_CP = -I/Users/hoscsaba/program/CoolProp/include -I/Users/hoscsaba/program/CoolProp/externals/fmtlib
LINK_CP = -L/Users/hoscsaba/program/CoolProp/build1 -lCoolProp

all:$(TARGETS)
	del libmy_tools.a
	ar rvs libmy_tools.a my_tools.o
	ranlib libmy_tools.a
	del libPSToolbox.a
	ar rvs libPSToolbox.a Gas.o IdealGas.o FrozenMixtureLiquidGas.o Units.o LWP.o SCP.o Reservoir.o Valve.o Connector.o Valve_with_Absorber.o
	ranlib libPSToolbox.a

my_tools: my_tools.cpp 
	$(CXX) $(INC) $(CFLAGS) my_tools.cpp -c -o my_tools.o

Gas: Gas.cpp 
	$(CXX) $(CFLAGS) Gas.cpp -c -o Gas.o

IdealGas: IdealGas.cpp 
	$(CXX) $(CFLAGS)  -c -o IdealGas.o  IdealGas.cpp

FrozenMixtureLiquidGas: FrozenMixtureLiquidGas.cpp 
	$(CXX) $(CFLAGS) FrozenMixtureLiquidGas.cpp -c -o FrozenMixtureLiquidGas.o

Units: Units.cpp
	$(CXX) $(CFLAGS) Units.cpp -c -o Units.o

LWP: 
	$(CXX) $(INC) $(CFLAGS) LWP.cpp -c -o LWP.o

SCP: 
	$(CXX) $(INC) $(CFLAGS) SCP.cpp -c -o SCP.o

Reservoir: 
	$(CXX) $(INC) $(CFLAGS) Reservoir.cpp -c -o Reservoir.o
# 	$(CXX) $(INC) $(LINK) $(CFLAGS) Reservoir.cpp -c -o Reservoir.o

Valve:
	$(CXX) $(INC) $(CFLAGS) Valve.cpp -c -o Valve.o
# 	$(CXX) $(INC) $(LINK) $(CFLAGS) Valve.cpp -c -o Valve.o

Valve_with_Absorber:
	$(CXX) $(INC) $(CFLAGS) Valve_with_Absorber.cpp -c -o Valve_with_Absorber.o

Connector: Connector.cpp
	$(CXX) $(INC) $(CFLAGS) Connector.cpp -c -o Connector.o

#CoolPropGas: CoolPropGas.cpp
#	$(CXX) $(INC_CP) $(LINK_CP) $(CFLAGS) CoolPropGas.cpp -c -o CoolPropGas.o

clean:
	$(RM) *.o *.a
