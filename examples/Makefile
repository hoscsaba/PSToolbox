CXX = g++

TARGETS = \
test1_reservoir_with_ideal_gas \
test2_reservoir_valve_with_liquid \
test3_reservoir_pipe_with_ideal_gas\
test4_reservoir_valve_with_ideal_gas\
test5_2pipes_with_liquid \
test6_3pipes_with_liquid \
test7_2pipes_with_ideal_gas \
test8_frozen_mixture \
test11_runner \
test12_Pump \
test13_Pump_and_CheckValve.cpp \
test14_Pump_non-reflecting_BC \
test15_epanet.cpp \
test16_3pipes_with_ideal_gas

INC = -I/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 \
-I/Users/hoscsaba/Documents/GitHub/PSToolbox/src\
-L/Users/hoscsaba/Documents/GitHub/PSToolbox/lib\
-L/opt/homebrew/opt/python@3.12/Frameworks/Python.framework/Versions/3.12/lib/python3.12/config-3.12-darwin\
-I/Library/Developer/CommandLineTools/Library/Frameworks/Python3.framework/Versions/3.9/Headers\
-I/opt/homebrew/Cellar/numpy/2.1.1/lib/python3.12/site-packages/numpy/_core/include/

INC_CP = -I/Users/hoscsaba//Documents/GitHub/CoolProp/include -I/Users/hoscsaba//Documents/GitHub/CoolProp/externals/fmtlib
#LINK = -lPSToolbox -lmy_tools -lpython3.12
LINK = -lPSToolbox
LINK_CP = -L/Users/hoscsaba/program/CoolProp/build1 -lCoolProp
CFLAGS = -w -std=c++11 -pedantic -O3 -Wall -Wno-c++11-long-long

all:$(TARGETS)

test1_reservoir_with_ideal_gas: test1_reservoir_with_ideal_gas.cpp
	clear
	$(CXX) $(INC) $(LINK) $(CFLAGS) test1_reservoir_with_ideal_gas.cpp
	./a.out

test2_reservoir_valve_with_liquid: test2_reservoir_valve_with_liquid.cpp
	clear
	$(CXX) $(INC) $(LINK) $(CFLAGS) test2_reservoir_valve_with_liquid.cpp
	./a.out

test3_reservoir_pipe_with_ideal_gas: test3_reservoir_pipe_with_ideal_gas.cpp
	clear
	$(CXX) $(INC) $(LINK) $(CFLAGS) test3_reservoir_pipe_with_ideal_gas.cpp
	./a.out

test4_reservoir_valve_with_ideal_gas: test4_reservoir_valve_with_ideal_gas.cpp
	clear
	$(CXX) $(INC) $(LINK) $(CFLAGS) test4_reservoir_valve_with_ideal_gas.cpp
	./a.out

test5_2pipes_with_liquid: test5_2pipes_with_liquid.cpp
	clear
	$(CXX) $(INC) $(LINK) $(CFLAGS) test5_2pipes_with_liquid.cpp
	./a.out

test6_3pipes_with_liquid: test6_3pipes_with_liquid.cpp
	clear
	$(CXX) $(INC) $(LINK) $(CFLAGS) test6_3pipes_with_liquid.cpp
	./a.out

test7_2pipes_with_ideal_gas: test7_2pipes_with_ideal_gas.cpp
	clear
	$(CXX) $(INC) $(LINK) $(CFLAGS) test7_2pipes_with_ideal_gas.cpp
	./a.out

test8_frozen_mixture: test8_frozen_mixture.cpp
	clear
	$(CXX) $(INC) $(LINK) $(CFLAGS) test8_frozen_mixture.cpp
	./a.out

test9_CoolProp_base: test9_CoolProp_base.cpp
	clear
	$(CXX) $(INC_CP) $(LINK_CP) $(CFLAGS) test9_CoolProp_base.cpp
	./a.out

test10_reservoir_valve_with_hydrogen: test10_reservoir_valve_with_hydrogen.cpp
	clear
	$(CXX) $(INC) $(INC_CP) $(LINK) $(LINK_CP) $(CFLAGS) test10_reservoir_valve_with_hydrogen.cpp
	./a.out

test11_runner: test11_runner.cpp
	clear
	$(CXX) $(INC) $(LINK)  $(CFLAGS) test11_runner.cpp
	./a.out

test12_Pump: test12_Pump.cpp
	clear
	$(CXX) $(INC) $(INC_CP) $(LINK) $(LINK_CP) $(CFLAGS) test12_Pump.cpp
	./a.out

test13_Pump_and_CheckValve: test13_Pump_and_CheckValve.cpp
	clear
	$(CXX) $(INC) $(LINK) $(CFLAGS) test13_Pump_and_CheckValve.cpp
	./a.out

test14_Pump_non-reflecting_BC: test14_Pump_non-reflecting_BC.cpp
	clear
	$(CXX) $(INC) $(INC_CP) $(LINK) $(LINK_CP) $(CFLAGS) test14_Pump_non-reflecting_BC.cpp
	./a.out

test15_epanet: test15_epanet.cpp
	clear
	$(CXX) $(INC) $(LINK)	 $(CFLAGS) test15_epanet.cpp
	./a.out

test16_3pipes_with_ideal_gas: test16_3pipes_with_ideal_gas.cpp
	clear
	$(CXX) $(INC) $(LINK)	 $(CFLAGS) test16_3pipes_with_ideal_gas.cpp
	./a.out

clean:
	$(RM) *.o *.a
