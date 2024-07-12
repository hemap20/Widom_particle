# Compiler
CXX = g++
# Compiler flags
CXXFLAGS = -std=c++11 -Wall -I/usr/include/gsl
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
