CC=g++
SDIR=./src
ODIR=./objects
IDIR=./include
BOOSTLIB=/opt/boost_1_82_0/stage/lib
BOOSTINC=/opt/boost_1_82_0
CFLAGS=-Wall -I/usr/local/include -I$(BOOSTINC) -I$(IDIR) -std=gnu++14 -c -D_USE_MATH_DEFINES -O3 -fopenmp
LDFLAGS=-L/usr/local/lib -L$(BOOSTLIB) -lfftw3 -lfftw3f -lm -fopenmp
_SRCS=Params.cpp Pulse.cpp Noise.cpp generate_dataset.cpp
_HEADS=Constants.hpp DataOps.hpp Params.hpp Pulse.hpp Noise.hpp generate_dataset.hpp 

OBJECTS=$(patsubst %,$(ODIR)/%,$(_SRCS:.cpp=.o))
HEADERS=$(patsubst %,$(IDIR)/%,$(_HEADS))
SOURCES=$(patsubst %,$(SDIR)/%,$(_SRCS))
EXECUTABLE ?= generate_dataset

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

$(ODIR)/%.o: $(SDIR)/%.cpp $(HEADERS) 
	$(CC) $(CFLAGS) -c $< -o $@ 

.PHONY: clean

clean:
	rm $(ODIR)/*.o $(EXECUTABLE)

# DO NOT DELETE
