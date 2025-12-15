#include "acfDriver.h"

// Constructor
acfDriver::acfDriver()
{
	nRun = 0;
	nCor = 0;
	t0 = 0.0;
	dt = 1.0; // ns
}

// Deconstructor
acfDriver::~acfDriver(){}

// Methods to put time series in
void acfDriver::readDataFromFile(std::string f)
{
	double vr, vi;

	std::fstream inFile;
	inFile.open(f.c_str(), std::ios::in);
	while (!inFile.eof())
	{
		inFile >> vr >> vi;
		if (inFile.eof()) break;
		data0.push_back(dcomplex(vr,vi));
	}
	inFile.close();

	nRun = data0.size();

	return;
}

void acfDriver::setData(cdvector d)
{
	data0.clear();
	data0 = d;
	nRun = data0.size();
	return;
}

void acfDriver::setCorrFun(dvector d)
{
	acf.clear();
	acf = d;
	nRun = acf.size();
	return;
}

// Method to determine ''real'' time (in ns) and frequency (in ns^-1)
void acfDriver::setTimeData(double startT, double deltaT)
{
	t0 = startT;
	dt = deltaT;
	w0 = 0.0;
	return;
}

// Method to calculate the ACF
void acfDriver::calculateACF(int n)
{
	nCor = n;
	if (nCor > nRun)
	{
		std::cout << std::endl << "ERROR in " << __FILE__ << ", line " << __LINE__ << ": correlation time (nCor = " << nCor << ") must be less or equal to MD run time (nRun = " << nRun << ")." << std::endl << std::endl;
		exit(1);
	}

	double norm = 1.0 / double(nRun);
	fftw_plan p;
	dataIn = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nRun);
	dataTmp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nRun);
	dataOut = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nRun);

	// Step 1: FT transform of time serie
	p = fftw_plan_dft_1d(nRun, dataIn, dataOut, FFTW_FORWARD, FFTW_ESTIMATE);
	for (int i = 0; i < nRun; i++)
	{
		dataIn[i][0] = data0[i].real(); dataIn[i][1] = data0[i].imag();
	}
	fftw_execute(p);

	// Step 2: inverse FT transform of squared absolute value of transformed time serie
	for (int i = 0; i < nRun; i++)
		dataTmp[i][0] = dataOut[i][0] * dataOut[i][0] + dataOut[i][1] * dataOut[i][1];

	p = fftw_plan_dft_1d(nRun, dataIn, dataOut, FFTW_BACKWARD, FFTW_ESTIMATE);
	for (int i = 0; i < nRun; i++)
	{
		dataIn[i][0] = dataTmp[i][0];
		dataIn[i][1] = 0.0;
	}
	fftw_execute(p);

	acf.clear();
	for (int i = 0; i < nCor; i++) acf.push_back(dataOut[i][0] * norm * norm);

	fftw_destroy_plan(p);
	fftw_free(dataIn);
	fftw_free(dataTmp);
	fftw_free(dataOut);

	return;
}

// Method to calculate the Real Part of the Spectral Density
void acfDriver::calculateRealSpectralDensity(void)
{
	nCor = nRun;

	double norm = dt;
	fftw_plan p;
	corrfunIn = (double*)fftw_malloc(sizeof(double) * nRun);
	dataOut = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nRun);

	// Step 1: FT transform of time serie
	p = fftw_plan_dft_r2c_1d(nRun, corrfunIn, dataOut, FFTW_ESTIMATE);
	for (int i = 0; i < nRun; i++)
	{
		corrfunIn[i] = acf[i];
	}
	fftw_execute(p);

	spd.clear();
	for (int i = 0; i < nRun; i++) spd.push_back(dataOut[i][0] * norm);

	fftw_destroy_plan(p);
	fftw_free(corrfunIn);
	fftw_free(dataOut);

	return;
}

// Output methods
dvector acfDriver::getACF()
{
	return acf;
}

dvector acfDriver::getRealSpectralDensity()
{
	return spd;
}

void acfDriver::saveACF(std::string f)
{
	std::fstream outFile;
	outFile.open(f.c_str(), std::ios::out);
	for (int i = 0; i < acf.size(); ++i)
		outFile << (t0 + (double)i*dt) << "\t" << acf.at(i) << std::endl;
	outFile.close();
	return;
}

void acfDriver::saveRealSpectralDensity(std::string f)
{
	std::fstream outFile;
	outFile.open(f.c_str(), std::ios::out);
	for (int i = 0; i < nRun/2; ++i)
		outFile << 1.0e9 * 2.0 * M_PI * (double)i / (dt * (double)nRun) << "\t" << 1.0e-9 * spd.at(i) << std::endl;
	outFile.close();
	return;
}

dvector acfDriver::getFrequencies(void)
{
	dvector f(nRun, 0.0);
	for (int i = 0; i < nRun; ++i)
		f[i]= 1.0e9 * 2.0 * M_PI * (double)i / (dt * (double)nRun);
	return f;
}

////////////////////////////////////
// Methods for cross-correlations //
////////////////////////////////////

void acfDriver::setData(cdvector dA, cdvector dB)
{
	dataA0.clear();
	dataA0 = dA;

	dataB0.clear();
	dataB0 = dB;

	nRun = dataA0.size();

	return;
}

void acfDriver::calculateCCF(int n)
{
	nCor = n;
	if (nCor > nRun)
	{
		std::cout << std::endl << "ERROR in " << __FILE__ << ", line " << __LINE__ << ": correlation time (nCor = " << nCor << ") must be less or equal to MD run time (nRun = " << nRun << ")." << std::endl << std::endl;
		exit(1);
	}

	double norm = 1.0 / double(nRun);
	fftw_plan p;
	dataAIn =  (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nRun);
	dataBIn =  (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nRun);
	dataTmp =  (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nRun);
	dataAOut = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nRun);
	dataBOut = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nRun);

	// Step 1: FT transform of time series
	p = fftw_plan_dft_1d(nRun, dataAIn, dataAOut, FFTW_FORWARD, FFTW_ESTIMATE);
	for (int i = 0; i < nRun; i++)
	{
		dataAIn[i][0] = dataA0[i].real(); dataAIn[i][1] = dataA0[i].imag();
	}
	fftw_execute(p);

	p = fftw_plan_dft_1d(nRun, dataBIn, dataBOut, FFTW_FORWARD, FFTW_ESTIMATE);
	for (int i = 0; i < nRun; i++)
	{
		dataBIn[i][0] = dataB0[i].real(); dataBIn[i][1] = dataB0[i].imag();
	}
	fftw_execute(p);

	// Step 2: inverse FT transform of product (A* x B) of transformed time serie
	for (int i = 0; i < nRun; i++)
	{
		dataTmp[i][0] = dataAOut[i][0] * dataBOut[i][0] + dataAOut[i][1] * dataBOut[i][1];
		dataTmp[i][1] = dataAOut[i][0] * dataBOut[i][1] - dataAOut[i][1] * dataBOut[i][0];
	}

	p = fftw_plan_dft_1d(nRun, dataAIn, dataAOut, FFTW_BACKWARD, FFTW_ESTIMATE);
	for (int i = 0; i < nRun; i++)
	{
		dataAIn[i][0] = dataTmp[i][0];
		dataAIn[i][1] = dataTmp[i][1];
	}
	fftw_execute(p);

	ccf.clear();
	for (int i = 0; i < nCor; i++) ccf.push_back(std::complex<double>(dataAOut[i][0] * norm * norm, dataAOut[i][1] * norm * norm));

	fftw_destroy_plan(p);
	fftw_free(dataAIn);
	fftw_free(dataBIn);
	fftw_free(dataTmp);
	fftw_free(dataAOut);
	fftw_free(dataBOut);

	return;
}

void acfDriver::setCrossCorrFun(cdvector d)
{
	ccf.clear();
	ccf = d;
	nRun = ccf.size();
	return;
}

cdvector acfDriver::getCCF()
{
	return ccf;
}

void acfDriver::calculateRealCrossSpectralDensity(void)
{
	nCor = nRun;

	double norm = dt;
	fftw_plan p;
	dataIn  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nRun);
	dataOut = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nRun);

	// Step 1: FT transform of time serie
	p = fftw_plan_dft_1d(nRun, dataIn, dataOut, FFTW_FORWARD, FFTW_ESTIMATE);
	for (int i = 0; i < nRun; i++)
	{
		dataIn[i][0] = ccf[i].real();
		dataIn[i][1] = ccf[i].imag();
	}
	fftw_execute(p);

	spd.clear();
	for (int i = 0; i < nRun; i++) spd.push_back(dataOut[i][0] * norm);

	fftw_destroy_plan(p);
	fftw_free(dataIn);
	fftw_free(dataOut);

	return;
}

