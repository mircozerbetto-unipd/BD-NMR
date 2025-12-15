include ../conf/make.COMPILERS.CPU
include ../conf/make.in

# Override for debugging
#COMPOPT=-g -pg -O0 -std=c++17 -DADD_ -DDEBUG

all: prunetrj trj2omega relax nmr

prunetrj: prune_trj.o
	$(CPP) $(COMPOPT) prune_trj.o $(LL) $(II) -o prunetrj

prune_trj.o: prune_trj.cpp
	$(CPP) $(COMPOPT) $(II) -c prune_trj.cpp

trj2omega: acfDriver.o nmr.o spline.o fit_cf.o extract_Omega.o
	$(CPP) $(COMPOPT) acfDriver.o nmr.o spline.o fit_cf.o extract_Omega.o $(LL) $(II) -lEulerQuaternion -lfftw3 -lcminpack -o getOmega

relax: acfDriver.o nmr.o spline.o fit_cf.o relax.o
	$(CPP) $(COMPOPT) acfDriver.o nmr.o spline.o fit_cf.o relax.o $(LL) $(II) -lEulerQuaternion -lfftw3 -lcminpack -o relax

nmr: acfDriver.o nmr.o spline.o fit_cf.o main.o
	$(CPP) $(COMPOPT) acfDriver.o nmr.o spline.o fit_cf.o main.o $(LL) $(II) -lEulerQuaternion -lfftw3 -lcminpack -o nmr

acfDriver.o: acfDriver.cpp acfDriver.h
	$(CPP) $(COMPOPT) $(II) -c acfDriver.cpp

nmr.o: nmr.cpp nmr.h
	$(CPP) $(COMPOPT) $(II) -c nmr.cpp

spline.o: spline.cpp spline.h
	$(CPP) $(COMPOPT) $(II) -c spline.cpp

fit_cf.o: fit_cf.cpp
	$(CPP) $(COMPOPT) $(II) -c fit_cf.cpp

main.o: main.cpp
	$(CPP) $(COMPOPT) $(II) -c main.cpp

extract_Omega.o: extract_Omega.cpp
	$(CPP) $(COMPOPT) $(II) -c extract_Omega.cpp

relax.o: relax.cpp
	$(CPP) $(COMPOPT) $(II) -c relax.cpp

clean:
	rm -rf ./nmr ./relax ./getOmega ./*.o

