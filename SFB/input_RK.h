#ifndef INPUT_H_
#define INPUT_H_

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <string>
#include <ctime>

#include "types.h"

#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "magma_v2.h"
#include "magma_lapack.h"

class input
{
	public:

		input();
		~input();

		void parseInputFile(char *);

		uint64_t getNSteps();
		uint64_t getNDump();

		int getNx();

		double getIntegrationDt();
		void getRCDiffusionTensor(double*);

		void getInitialState(Quaternion *, double *);

		uint64_t getSeed();

		double getFreqScale();
		double getTimeScale();
		
	private:

		void checkInputIntegrity();
		void readDiffusionTensor();
		void setInitialConfiguration();

		int nStepsFound, nDumpFound, nxFound, rcdFound, dtFound, seedFound, nErrors;

		uint64_t nSteps, nDump;
		int nx, nc;
		uint64_t seed;
		double dt, scale;
		double *D, *xi;
		Quaternion qi;

		std::string RCDiffFile;

};

#endif

