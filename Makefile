# Compiler
CXX = g++
# Compiler flags
CXXFLAGS = -std=c++11 -Wall
# GSL library flags
LIBS = -lgsl -lgslcblas

# Source files
SRCS = main.cpp input_func.cpp output_func.cpp pairwise_dist.cpp pot_energy.cpp forces.cpp insert.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable
TARGET = main

# Rule to link the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS)

# Rule to compile .cpp files to .o files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Phony target to clean objects and executable
.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET)
