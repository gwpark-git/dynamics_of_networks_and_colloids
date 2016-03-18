
CC=icpc
CFLAGS=-c -openmp -O2 -Wall -mkl -L/usr/local/include/ -L/usr/local/lib/ -lgsl
LDFLAGS=-openmp -O2 -Wall -mkl -L/usr/local/include/ -L/usr/local/lib/ -lgsl
SOURCES=src/Brownian_simulation.cpp lib/association.cpp lib/connectivity.cpp lib/geometry.cpp lib/handle_association.cpp lib/matrix.cpp lib/parallel.cpp lib/potential.cpp lib/random.cpp lib/read_file_condition.cpp lib/time_evolution.cpp lib/trajectory.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=stochastic_simulation

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
