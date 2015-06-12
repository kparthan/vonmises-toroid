# kernel-compatible version
#CFLAGS=-std=c++98 -c -O3 -I/home/parthan/external_libs/ -fopenmp
#LDFLAGS=-static -lboost_program_options -lboost_filesystem -fopenmp -lnlopt -lgsl -lgslcblas -lm

CFLAGS=-std=c++0x -c -O3 -fopenmp
#CFLAGS=-std=c++0x -g -c -fopenmp
LDFLAGS=-lboost_program_options -lboost_system -lboost_filesystem -fopenmp -lnlopt -lm 

OBJECTS = main.o \
  Support.o \
  Normal.o  \
  vMC.o \
  vMF.o \
  FB4.o \
  MultivariateNormal.o \
  ACG.o \
  Bingham.o \
  Kent.o \
  Kent_EccTrans.o \
  Kent_UnifTrans.o \
  FB6.o \
  Optimize.o \
  Optimize2.o \
  Structure.o \
  Mixture.o \
  Mixture_vMF.o \
  Test.o \
  Experiments.o

all: main 

main: $(OBJECTS)
	g++ $(OBJECTS) -o $@ $(LDFLAGS) 

main.o: main.cpp 
	g++ $(CFLAGS) $< -o $@

Support.o: Support.cpp Support.h Header.h UniformRandomNumberGenerator.h
	g++ $(CFLAGS) $< -o $@

Normal.o: Normal.cpp Normal.h 
	g++ $(CFLAGS) $< -o $@

vMF.o: vMF.cpp vMF.h Header.h
	g++ $(CFLAGS) $< -o $@

vMC.o: vMC.cpp vMC.h Header.h
	g++ $(CFLAGS) $< -o $@

FB4.o: FB4.cpp FB4.h Header.h
	g++ $(CFLAGS) $< -o $@

MultivariateNormal.o: MultivariateNormal.cpp MultivariateNormal.h Header.h
	g++ $(CFLAGS) $< -o $@

ACG.o: ACG.cpp ACG.h Header.h
	g++ $(CFLAGS) $< -o $@

Bingham.o: Bingham.cpp Bingham.h Header.h
	g++ $(CFLAGS) $< -o $@

Kent.o: Kent.cpp Kent.h Header.h
	g++ $(CFLAGS) $< -o $@

Kent_EccTrans.o: Kent_EccTrans.cpp Kent_EccTrans.h Header.h
	g++ $(CFLAGS) $< -o $@

Kent_UnifTrans.o: Kent_UnifTrans.cpp Kent_UnifTrans.h Header.h
	g++ $(CFLAGS) $< -o $@

FB6.o: FB6.cpp FB6.h Header.h
	g++ $(CFLAGS) $< -o $@

Structure.o: Structure.cpp Structure.h Header.h
	g++ $(CFLAGS) $< -o $@

Mixture.o: Mixture.cpp Mixture.h Header.h
	g++ $(CFLAGS) $< -o $@

Mixture_vMF.o: Mixture_vMF.cpp Mixture_vMF.h Header.h
	g++ $(CFLAGS) $< -o $@

Optimize.o: Optimize.cpp Optimize.h Header.h
	g++ $(CFLAGS) $< -o $@

Optimize2.o: Optimize2.cpp Optimize2.h Header.h
	g++ $(CFLAGS) $< -o $@

Test.o: Test.cpp Test.h Header.h
	g++ $(CFLAGS) $< -o $@

Experiments.o: Experiments.cpp Experiments.h Header.h
	g++ $(CFLAGS) $< -o $@

clean:
	rm -f *.o *~ main gmon.out 

