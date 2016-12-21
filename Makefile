CXX = g++

CSMINCDIRS = -I./linux64/include -I./linux64/include/csm
CSMLIBDIRS = -L./linux64/lib

CXXFLAGS = -Wall -ansi -fPIC
CSHAREDFLAGS = -shared

LDLIBS = -ldl -lcsmapi

OBJS = PinholePluginCSM302.o PinholeSensorModelCSM302.o
EXEC = pinhole.o

#LDFLAGS = -Wl,-E -Xlinker '-z' -Xlinker 'origin' -Xlinker "-rpath" -Xlinker ".:./linux64/lib"
LDFLAGS = -Wl,-E,-rpath="./linux64/lib"

.PHONY: all clean plugin 

all: $(EXEC)
	$(CXX) $(CSMLIBDIRS) $(CXXFLAGS) $(LDFLAGS) $(EXEC) -o pinhole $(LDLIBS)

$(EXEC): plugin

plugin: $(OBJS)
	mkdir -p plugins
	$(CXX) $(CSHAREDFLAGS) $(CXXFLAGS) $(OBJS) -o PinholePluginCSM302.so
	mv PinholePluginCSM302.so plugins/

%.o : %.cpp
	$(CXX) $(CSMINCDIRS) -c $(CXXFLAGS) $< -o $@

clean:
	rm -rf *.o
	rm -rf *.so
	rm -f pinhole
