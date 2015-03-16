include Makefile.config
# use make USE_GUROBI=1 to compile for Gurobi
CC=g++
ifeq ($(DEBUG),1)
 OGDF_LIBRARY_DIR=$(OGDF_LIBRARY_DIR_DEBUG)
else
 OGDF_LIBRARY_DIR=$(OGDF_LIBRARY_DIR_RELEASE)
endif
INCLUDES=-I$(OGDF_INCLUDE_DIR)
CXXFLAGS=-std=c++11 -W -Wall -Wextra -Wno-unused-parameter
LDFLAGS=-pthread -lm -L$(OGDF_LIBRARY_DIR)
ifeq ($(DEBUG),1)
 CXXFLAGS += -O0 -g -ggdb -W -DOGDF_DEBUG
else
 ifeq ($(DEBUG),2)
  # profiler mode but release
  CXXFLAGS += -O3 -g -ggdb -pg
  LDFLAGS += -pg
 else
  CXXFLAGS += -O3
 endif
endif
LDLIBS=-lOGDF -lCOIN
ifdef USE_GUROBI
 BINARY_SUFFIX = _gurobi
 INCLUDES += -I$(GUROBI_INCLUDE_DIR)
 CXXFLAGS += -DLP_SOLVER_GUROBI
 LDFLAGS += -L$(GUROBI_LIBRARY_DIR)
 LDLIBS += -lgurobi$(GUROBI_VERSION) -lgurobi_c++
else
 BINARY_SUFFIX = _cplex
 INCLUDES += -I$(CPLEX_INCLUDE_DIR) -I$(CONCERT_INCLUDE_DIR)
 CXXFLAGS += -DIL_STD
 LDFLAGS += -L$(CPLEX_LIBRARY_DIR) -L$(CONCERT_LIBRARY_DIR)
 LDLIBS += -lilocplex -lcplex -lconcert
endif
CXXFLAGS += $(INCLUDES)
BINARIES=solve$(BINARY_SUFFIX) convert-individual delnode partitioner

all: $(BINARIES)

solve$(BINARY_SUFFIX): solve$(BINARY_SUFFIX).o MincutPushRelabel.o
convert-individual: convert-individual.o
delnode: delnode.o
partitioner: partitioner.o

solve$(BINARY_SUFFIX).o: WeightedNodeColoredGraph.h
convert-individual.o: WeightedNodeColoredGraph.h
delnode.o: WeightedNodeColoredGraph.h
partitioner.o: WeightedNodeColoredGraph.h


solve$(BINARY_SUFFIX).cc: solve.cc
	# This is bad style somehow, but I don't know how to effectively tell the Makefile that it should use solve.cc
	cp "$<" "$@"

%o: %cc

../Makefile.config:
	$(error Need to know OGDF, CPLEX and CONCERT paths. Please rename Makefile.config.default to Makefile.config and edit appropriately)

.PHONY: clean distclean
clean:
	$(RM) \
		$(BINARIES) \
		$(BINARIES:%=%.o)
