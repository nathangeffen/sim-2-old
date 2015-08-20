CXX = g++
MODEL = dumb

CXXFLAGS = -pthread -Wall -std=c++11
DEVFLAGS  = -g -rdynamic
RELFLAGS = -O3
LDFLAGS =

# the build target executable:
SOURCES = sim.cc test.cc $(MODEL).cc
OBJECTS = $(SOURCES:.cc=.o)
DEPEND =  $(OBJECTS:%.o=%.d)
EXECUTABLE = sim-$(MODEL)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(DEVFLAGS) $(LDFLAGS) $(OBJECTS) -o $(EXECUTABLE)-dev

%.o: %.cc
	$(CXX) -c $(DEVFLAGS) $(CXXFLAGS)  -MD -MP -MF ${@:.o=.d} $< -o $@

release: clean
	$(CXX) $(RELFLAGS) $(CXXFLAGS) -o $(EXECUTABLE)-rel $(SOURCES)

clean:
	$(RM) $(EXECUTABLE)-release $(EXECUTABLE)-dev *.o

-include $(DEPEND)
