# Various flags
CXX  = g++
LINK = $(CXX)
#CXXFLAGS = -I -Wall -g 
CXXFLAGS = -g -Wall -O3 -fopenmp #-I #-Wall -O3 -funroll-loops -pipe 
#CXXFLAGS = -g -Wall -fopenmp #-I #-Wall -O3 -funroll-loops -pipe 
LFLAGS =  -g -fopenmp -lm

TARGET  = ompRelaxmap

HEADER  = Node.h Module.h FileIO.h timing.h
FILES = OmpRelaxmap.cpp Node.cpp Module.cpp FileIO.cpp timing.cpp

OBJECTS = $(FILES:.cpp=.o)

$(TARGET): ${OBJECTS}
	$(LINK) $(LFLAGS) $^ -o $@

all: $(TARGET)

clean:
	rm -f $(OBJECTS)

distclean:
	rm -f $(OBJECTS) $(TARGET)

# Compile and dependency
$(OBJECTS): $(HEADER) Makefile




