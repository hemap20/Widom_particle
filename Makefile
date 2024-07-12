# Compiler
CXX = g++
# Compiler flags
CXXFLAGS = -std=c++11 -Wall -I~/Codes/basics/Basics/gsl-2.8/gsl
# GSL library flags
LIBS = -lgsl -lgslcblas

# Source files (adjust paths as per your directory structure)
SRCS = basics/Basics/main.C basics/Basics/input_func.C basics/Basics/output_func.C basics/Basics/pairwise_dist.C basics/Basics/pot_energy.C basics/Basics/forces.C basics/Basics/insert.C

# Object files
OBJS = $(SRCS:.C=.o)

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
