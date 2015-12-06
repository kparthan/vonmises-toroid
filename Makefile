# kernel-compatible version
#CFLAGS=-std=c++98 -c -O3 -I/home/parthan/external_libs/ -fopenmp
#LDFLAGS=-static -lboost_program_options -lboost_filesystem -fopenmp -lnlopt -lgsl -lgslcblas -lm

CFLAGS=-std=c++0x -c -O3 -fopenmp
#CFLAGS=-std=c++0x -g -c -fopenmp
#CFLAGS=-std=c++0x -pg -c -fopenmp
#LDFLAGS=-pg -fopenmp -lnlopt -lm -lboost_program_options -lboost_system -lboost_filesystem -lgsl -lgslcblas 
LDFLAGS=-fopenmp -lnlopt -lm -lboost_program_options -lboost_system -lboost_filesystem -lgsl -lgslcblas 

OBJECTS = main.o \
  Support.o \
  Structure.o \
  Normal.o  \
  vMC.o \
  Mixture_vMC.o \
  BVM_Ind.o \
  BVM_Sine.o \
  MarginalDensitySine.o \
  BVM_Cosine.o \
  MarginalDensityCosine.o \
  KappaSolver.o \
  OptimizeSine.o \
  OptimizeInd.o \
  OptimizeCosine.o \
  Mixture_Ind.o \
  Mixture_Sine.o \
  Test.o \
  Experiments.o

all: main 

main: $(OBJECTS)
	g++ $(OBJECTS) -o $@ $(LDFLAGS) 

main.o: main.cpp Header.h Support.h
	g++ $(CFLAGS) $< -o $@

Support.o: Support.cpp Support.h Header.h 
	g++ $(CFLAGS) $< -o $@

Structure.o: Structure.cpp Structure.h Header.h 
	g++ $(CFLAGS) $< -o $@

Normal.o: Normal.cpp Normal.h 
	g++ $(CFLAGS) $< -o $@

vMC.o: vMC.cpp vMC.h Header.h
	g++ $(CFLAGS) $< -o $@

Mixture_vMC.o: Mixture_vMC.cpp Mixture_vMC.h Header.h
	g++ $(CFLAGS) $< -o $@

BVM_Ind.o: BVM_Ind.cpp BVM_Ind.h Header.h
	g++ $(CFLAGS) $< -o $@

BVM_Sine.o: BVM_Sine.cpp BVM_Sine.h Header.h
	g++ $(CFLAGS) $< -o $@

MarginalDensitySine.o: MarginalDensitySine.cpp MarginalDensitySine.h Header.h
	g++ $(CFLAGS) $< -o $@

BVM_Cosine.o: BVM_Cosine.cpp BVM_Cosine.h Header.h
	g++ $(CFLAGS) $< -o $@

MarginalDensityCosine.o: MarginalDensityCosine.cpp MarginalDensityCosine.h Header.h
	g++ $(CFLAGS) $< -o $@

KappaSolver.o: KappaSolver.cpp KappaSolver.h Header.h
	g++ $(CFLAGS) $< -o $@

OptimizeInd.o: OptimizeInd.cpp OptimizeInd.h Header.h
	g++ $(CFLAGS) $< -o $@

OptimizeSine.o: OptimizeSine.cpp OptimizeSine.h Header.h
	g++ $(CFLAGS) $< -o $@

OptimizeCosine.o: OptimizeCosine.cpp OptimizeCosine.h Header.h
	g++ $(CFLAGS) $< -o $@

Mixture_Ind.o: Mixture_Ind.cpp Mixture_Ind.h Header.h
	g++ $(CFLAGS) $< -o $@

Mixture_Sine.o: Mixture_Sine.cpp Mixture_Sine.h Header.h
	g++ $(CFLAGS) $< -o $@

Test.o: Test.cpp Test.h Header.h Support.h
	g++ $(CFLAGS) $< -o $@

Experiments.o: Experiments.cpp Experiments.h Header.h Support.h
	g++ $(CFLAGS) $< -o $@

clean:
	rm -f *.o *~ main gmon.out 

