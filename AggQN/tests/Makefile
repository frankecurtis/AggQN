# Copyright (C) 2020 Albert Berahas, Frank E. Curtis, Baoyu Zhou
#
# This code is published under the ??? License.
#
# Author(s) : Albert Berahas, Frank E. Curtis, Baoyu Zhou

# C++ compiler
CXX = g++

# C++ compiler flags
CXXFLAGS = -g -Wall -std=c++11

# Set sources, etc.
headers = $(wildcard *.hpp)
sources = $(wildcard *.cpp)
objects = $(sources:.cpp=.o)
depends = $(sources:.cpp=.d)

# Executable(s)
EXES = $(sources:.cpp=)

# Libraries
AggQNLIB = $(AGGQNDIR)/AggQN/src/libAggQN.a

# Includes
INCLUDES = -I $(AGGQNDIR)/AggQN/src

# Rule for all
all: $(EXES)

# Rule for executable
$(EXES): % : %.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(AggQNLIB) -L $(LAPACKDIR) -ldl -lblas -llapack

# Dependencies for executable
$(EXES): $(AGGQNLIB)

# Rules for libraries
$(AGGQNLIB):
	$(MAKE) --directory=$(AGGQNDIR)/AggQN/src

# Dependencies for libraries
$(AggQNLIB): $(wildcard $(AGGQNDIR)/AggQN/src/*.hpp) $(wildcard $(AGGQNDIR)/AggQN/src/*.cpp)

# Rule for objects
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Objects depend on headers and sources
$(objects): $(headers)
$(objects): $(sources)

# Clean
clean:
	rm -f $(objects) $(depends)

# Very clean
veryclean: clean
	rm -f $(EXES)
