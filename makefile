CC=g++
SDIR=./src
ODIR=./objects
IDIR=./include
BOOSTLIB=/opt/boost_1_82_0/stage/lib
BOOSTINC=/opt/boost_1_82_0
HDF5LIB=-L/usr/local/lib -L/usr/local/hdf5/lib /usr/local/hdf5/lib/libhdf5_hl_cpp.a /usr/local/hdf5/lib/libhdf5_cpp.a /usr/local/hdf5/lib/libhdf5_hl.a /usr/local/hdf5/lib/libhdf5.a -lz -ldl -lm -Wl,-rpath -Wl,/usr/local/lib -ldl -lpthread -lrt
CFLAGS=-Wall -I/usr/local/include -I$(BOOSTINC) -I$(IDIR) -std=gnu++14 -c -D_USE_MATH_DEFINES -O3 -fopenmp
LDFLAGS=-L/usr/local/lib L$(BOOSTLIB) -L$(HDF5LIB) -L/usr/local/hdf5/lib /usr/local/hdf5/lib/libhdf5_hl_cpp.a /usr/local/hdf5/lib/libhdf5_cpp.a /usr/local/hdf5/lib/libhdf5_hl.a /usr/local/hdf5/lib/libhdf5.a -lz -ldl -lm -Wl,-rpath -Wl,/usr/local/lib -ldl -lpthread -lrt -lgomp -lfftw3 -lfftw3f -lm 
INCFLAGS=-I/usr/local/include -I/usr/local/hdf5/include -I/opt/hdf5-serial/src/H5FDsubfiling -I$(IDIR)

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
