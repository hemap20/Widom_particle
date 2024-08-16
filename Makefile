# Compiler
CXX = g++
# Compiler flags
CXXFLAGS = -std=c++11 -Wall
# Target executable
TARGET = main

# Source files
SRCS = mc_main.cpp random.cpp generate_positions.cpp pe_i.cpp output_func.cpp mc_eq.cpp pe_total.cpp mc_move.cpp #dist.cpp forces.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Header files
HDRS = pe_i.h generate_positions.h random.h mc_eq.h output_func.h pe_total.h mc_move.h #dist.h forces.h

# Default target
all: $(TARGET)

# Rule for building the target executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Rule for building object files
%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJS) $(TARGET)

# Phony targets
.PHONY: all clean
