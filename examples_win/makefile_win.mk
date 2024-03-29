CXX = g++

TARGETS = \
test1_reservoir_with_ideal_gas \
test2_reservoir_valve_with_liquid \
#test3_reservoir_pipe_with_ideal_gas\
#test4_reservoir_valve_with_ideal_gas\
#test5_2pipes_with_liquid \
#test6_3pipes_with_liquid \
#test7_2pipes_with_ideal_gas \
#test8_frozen_mixture

INC = -IC:/ProgramData/chocolatey/lib/eigen/include -ID:/BME-GPK_BSc/TDK/Szelep/Code/PSToolbox -LD:/BME-GPK_BSc/TDK/Szelep/Code/PSToolbox
LINK = -lPSToolbox -lmy_tools #-lpython2.7 -lCoolProp
CFLAGS = -std=c++11 -pedantic -O3 -Wall 

all:$(TARGETS)

test1_reservoir_with_ideal_gas: test1_reservoir_with_ideal_gas.cpp
	$(CXX) $(INC) test1_reservoir_with_ideal_gas.cpp $(LINK) $(CFLAGS)
	./a.exe

test2_reservoir_valve_with_liquid: test2_reservoir_valve_with_liquid.cpp
	$(CXX) $(INC) test2_reservoir_valve_with_liquid.cpp $(LINK)  $(CFLAGS) 
	./a.exe

test3_reservoir_pipe_with_ideal_gas: test3_reservoir_pipe_with_ideal_gas.cpp
	$(CXX) $(INC) test3_reservoir_pipe_with_ideal_gas.cpp $(LINK)   $(CFLAGS)
	./a.exe

test4_reservoir_valve_with_ideal_gas: test4_reservoir_valve_with_ideal_gas.cpp
	$(CXX) $(INC) test4_reservoir_valve_with_ideal_gas.cpp $(LINK)   $(CFLAGS)
	./a.exe

test5_2pipes_with_liquid: test5_2pipes_with_liquid.cpp
	$(CXX) $(INC) test5_2pipes_with_liquid.cpp $(LINK)   $(CFLAGS)
	./a.exe

test6_3pipes_with_liquid: test6_3pipes_with_liquid.cpp
	$(CXX) $(INC) test6_3pipes_with_liquid.cpp $(LINK)  $(CFLAGS) 
	./a.exe

test7_2pipes_with_ideal_gas: test7_2pipes_with_ideal_gas.cpp
	$(CXX) $(INC) test7_2pipes_with_ideal_gas.cpp $(LINK)   $(CFLAGS)
	./a.exe

test8_frozen_mixture: test8_frozen_mixture.cpp
	$(CXX) $(INC) test8_frozen_mixture.cpp $(LINK)   $(CFLAGS)
	./a.exe

test9_CoolProp_base: test9_CoolProp_base.cpp
 	clear
 	$(CXX) $(INC_CP) $(LINK_CP) $(CFLAGS) test9_CoolProp_base.cpp
 	./a.exe

# test10_reservoir_valve_with_hydrogen: test10_reservoir_valve_with_hydrogen.cpp
# 	clear
# 	$(CXX) $(INC) $(INC_CP) $(LINK) $(LINK_CP) $(CFLAGS) test10_reservoir_valve_with_hydrogen.cpp
# 	./a.out

clean:
	$(RM) *.o *.a
