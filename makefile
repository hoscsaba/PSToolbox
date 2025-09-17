# Compiler and flags
CXX = g++
# -w disables warnings.
CXXFLAGS = -w -O3 -Wall -Wextra -std=c++17 -g -I$(shell brew --prefix eigen)/include/eigen3

# Directories
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
LIB_DIR = lib

# Target
STATIC_LIB = $(LIB_DIR)/libPSToolbox.a

# Target
# TARGET = $(BIN_DIR)/main

# Source and object files
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
# Exclude specific cpp files
SRCS := $(filter-out $(SRC_DIR)/CoolPropGas.cpp $(SRC_DIR)/CoolPropHA.cpp $(SRC_DIR)/PSToolboxPlotter.cpp, $(SRCS))
OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))

# Default target
all: $(STATIC_LIB)

# Create static library from object files
$(STATIC_LIB): $(OBJS)
	@mkdir -p $(LIB_DIR)
#	libtool -static -o libmy_tools.a my_tools.o
#	libtool -static -o libPSToolbox.a RigidPipe.o PSToolBoxRunner.o PSToolboxBaseEdge.o Gas.o IdealGas.o FrozenMixtureLiquidGas.o Units.o LWP.o SCP.o Reservoir.o Valve.o Connector.o Valve_with_Absorber.o CheckValve.o Pump.o EpanetReader.o SurgeChamber.o Damper.o
	ar rcs $@ $^

# Compile source files into object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -rf $(OBJ_DIR) $(LIB_DIR)

.PHONY: all clean
