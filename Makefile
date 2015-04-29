CXX = g++

COMMONFLAGS = -pthread -Wall
CXXFLAGS  = -g -std=c++11 $(COMMONFLAGS)
RELFLAGS = -O3 -std=c++11 $(COMMONFLAGS)


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
