#ifndef ACFDRIVER_H_
#define ACFDRIVER_H_

#include <cstdio>
#include <cstdlib>

#include <fstream>
#include <sstream>
#include <iostream>

#include <cmath>
#include <complex>
#include <vector>

#include <string>

#include "fftw3.h"

#include "types.h"

class acfDriver
{
	public:

	 acfDriver();
	~acfDriver();

	void readDataFromFile(std::string);
	void setData(cdvector);
	void setCorrFun(dvector);
	void setTimeData(double, double);

	void calculateACF(int);
	void calculateRealSpectralDensity(void);

	dvector getACF();
	dvector getRealSpectralDensity();
	void saveACF(std::string);
	void saveRealSpectralDensity(std::string);

	dvector getFrequencies(void);

	// Methods for cross-correlations
	void setData(cdvector, cdvector);
	void setTimeSeries(cdvector, cdvector);
	void calculateCCF(int);
	cdvector getCCF(void);
	void setCrossCorrFun(cdvector);
	void calculateRealCrossSpectralDensity(void);

	private:

	int nRun, nCor;
	double t0, dt, w0, dw;

	dvector acf, spd, dataCF0;
	cdvector data0;
	double *corrfunIn;
	fftw_complex *dataIn, *dataTmp, *dataOut;

	// Data for cross-correlations
	cdvector dataA0, dataB0;
	fftw_complex *dataAIn, *dataBIn, *dataAOut, *dataBOut;
	cdvector ccf, dataCCF0;
};

#endif
