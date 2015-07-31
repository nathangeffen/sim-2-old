CXX = g++
MODEL = dumb

COMMONFLAGS = -pthread -Wall
DEVFLAGS  = -g -rdynamic -std=c++11 $(COMMONFLAGS)
RELFLAGS = -O3 -std=c++11 $(COMMONFLAGS)


# the build target executable:
ADD_SOURCES =
ADD_INCLUDES =
SOURCES = sim.cc test.cc
INCLUDES = sim.hh test.hh debug.hh
TARGET = sim_$(MODEL)

all: $(TARGET)

$(TARGET): $(SOURCES) $(INCLUDES) $(MODEL).cc $(ADD_SOURCES) $(ADD_INCLUDES)
	$(CXX) $(DEVFLAGS) $(CXXFLAGS) -o $(TARGET) $(SOURCES) $(MODEL).cc $(ADD_SOURCES)

release: clean
	$(CXX) $(RELFLAGS) $(CXXFLAGS) -o $(TARGET) $(SOURCES) $(MODEL).cc $(ADD_SOURCES)

clean:
	$(RM) $(TARGET)
