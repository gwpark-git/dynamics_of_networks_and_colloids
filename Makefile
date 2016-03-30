
CC=icpc
OPT_LEVEL=-O2
CFLAGS=-c -openmp $(OPT_LEVEL) -Wall -mkl -L/usr/local/include/ -L/usr/local/lib/ -lgsl
LDFLAGS=-openmp $(OPT_LEVEL) -Wall -mkl -L/usr/local/include/ -L/usr/local/lib/ -lgsl
SOURCES=src/Brownian_simulation.cpp lib/association.cpp lib/connectivity.cpp lib/geometry.cpp lib/handle_association.cpp lib/matrix.cpp lib/parallel.cpp lib/potential.cpp lib/random.cpp lib/read_file_condition.cpp lib/time_evolution.cpp lib/trajectory.cpp lib/cell_list.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=stochastic_simulation

all: $(SOURCES) $(EXECUTABLE)

debug: CFLAGS += -g -no_pie -Wl
debug: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f stochastic_simulation lib/*.o src/*.o

