# the compiler: gcc for C program, define as g++ for C++
CXX = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CXXFLAGS  = -g -Wall -std=c++11 -pthread
RELFLAGS = -O3 -std=c++11 -pthread

# the build target executable:
SOURCES = simple.cc sim.cc test.cc
INCLUDES = sim.hh test.hh
TARGET = sim

all: $(TARGET)

$(TARGET): $(SOURCES) $(INCLUDES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES)

release: clean
	$(CXX) $(RELFLAGS) -o $(TARGET) $(SOURCES)

clean:
	$(RM) $(TARGET)
