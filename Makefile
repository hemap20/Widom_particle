# Compiler
CXX = g++
# Compiler flags
CXXFLAGS = -std=c++11 -Wall -Wextra
# Target executable
TARGET = main

# Source files
SRCS = main.cpp rand_pos.cpp positions.cpp pe_i.cpp
# Object files
OBJS = $(SRCS:.cpp=.o)

# Header files
HDRS = pe_i.h positions.h rand_pos.h

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
