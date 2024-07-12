# Compiler
CXX = g++
# Compiler flags
CXXFLAGS = -std=c++11 -Wall -Wextra
# CFLAGS = -Wall -I~/Codes/basics/Basics/gsl-2.8/gsl  # Adjust the path as necessary
# GSL library flags
LIBS = 

# Source files (adjust paths as per your directory structure)
SRCS = main.cpp input_func.cpp output_func.cpp pairwise_dist.cpp pot_energy.cpp forces.cpp insert.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable
TARGET = main

# Rule to link the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS)

# Rule to compile .C files to .o files
%.o: %.C
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Phony target to clean objects and executable
.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET)