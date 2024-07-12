# Compiler
CXX = g++
# Compiler flags
CFLAGS = -O3 -Wall -L/home/hema/lib -I/home/hema/gsl/include   # Adjust the path as necessary
# CFLAGS = -Wall -I~/Codes/basics/Basics/gsl-2.8/gsl  # Adjust the path as necessary
# GSL library flags
LIBS = -L/usr/local/lib -lgsl -lgslcblas -lm

# Source files (adjust paths as per your directory structure)
SRCS = main.C input_func.C output_func.C pairwise_dist.C pot_energy.C forces.C insert.C

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
