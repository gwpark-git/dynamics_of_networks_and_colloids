
CC=icpc
GIT_VERSION := $(shell git describe --dirty --always --tags)
CFLAGS=-c -I/usr/local/include -Wall -DVERSION=\"$(GIT_VERSION)\" 
# LDFLAGS=-Wall -lmkl_intel_lp64 -lmkl_core  -lgsl
LDFLAGS=-Wall -L/usr/local/lib -lmkl_intel_lp64 -lmkl_core  -lgsl -Wl,-rpath,${MKLROOT}/lib
ifdef SDKROOT
        CFLAGS += -isysroot ${SDKROOT}
        LDFLAGS += -L${SDKROOT}/usr/lib
endif
# CFLAGS=-std=c++11 -c -Wall -mkl -I/usr/local/include/ -DVERSION=\"$(GIT_VERSION)\" 
# LDFLAGS=-Wall -mkl -L/usr/local/lib/ -lgsl 
SCOPE_FLAGS_I=-I/opt/exp_soft/unina.it/gsl-2.1/include/
SCOPE_FLAGS_L=-L/opt/exp_soft/unina.it/gsl-2.1/lib/
OPTFLAGS=-O2 -lmkl_intel_thread -liomp5 -lpthread -qopenmp
# OPTFLAGS=-openmp -O2
DEBUGFLAGS=-no_pie -g
SOURCES=src/stochastic_simulation.cpp src/brownian.cpp src/dumbbell_model.cpp src/repulsive_brownian.cpp src/stochastic_HEUR_flowers.cpp lib/association.cpp lib/connectivity.cpp lib/geometry.cpp lib/handle_association.cpp lib/matrix.cpp lib/parallel.cpp lib/potential.cpp lib/random.cpp lib/file_IO.cpp lib/time_evolution.cpp lib/trajectory.cpp lib/cell_list.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=stochastic_simulation

all: $(SOURCES) $(EXECUTABLE)

scope: CFLAGS += $(OPTFLAGS) $(SCOPE_FLAGS_I)
scope: LDFLAGS += $(OPTFLAGS) $(SCOPE_FLAGS_L)
scope: $(SOURCES) $(EXECUTABLE)

opt: CFLAGS += $(OPTFLAGS)
opt: LDFLAGS += $(OPTFLAGS)
opt: $(SOURCES) $(EXECUTABLE)

debug: CFLAGS += $(DEBUGFLAGS)
debug: LDFLAGS += $(DEBUGFLAGS)
debug: $(SOURCES) $(EXECUTABLE)

debug_scope_opt: CFLAGS += $(DEBUGFLAGS) $(SCOPE_FLAGS_I) $(OPTFLAGS)
debug_scope_opt: LDFLAGS += $(DEBUGFLAGS) $(SCOPE_FLAGS_L) $(OPTFLAGS)
debug_scope_opt: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f stochastic_simulation lib/*.o src/*.o

