#ifndef NMR_H
#define NMR_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>

#include "types.h"
#include "inputParser.h"

#define gyroH 2.675198e8    /* A m N^-1 s^-1 */  // Hydrogen magnetogyric ratio
#define hbar 1.054494e-34   /* N m  s        */  // Reduced Planck constant
#define MU0_OVER_4PI 1.0e-7 /* N A^-2        */  // Vacuum magnetic permittivity
#define TWO_OVER_FIFTHEEN  (2.0 / 15.0)
#define ONE_OVER_SIX (1.0 / 6.0)
#define SQRT_THREE_OVER_EIGHT sqrt(3.0 / 8.0)
#define PHASE(a) (a%2 == 0 ? 1.0 : -1.0)

class nmr {

public:
	nmr();
	virtual ~nmr();

	void setup(molecule, inputParserSrb);
	void setProbe(std::string);
	std::string getProbe(void);

	void setDeltaCSA(double);
	double getDeltaCSA(void);

	void setOmegaDC(double, double, double);
	void getOmegaDC(double*);

	void setOmegaDD(double, double, double);
	void getOmegaDD(double*);

	void setOmegaD(double, double, double);
	void getOmegaD(double*);

	int getNLarmorFreq(void);
	double getLarmorFreq(int);
	void setSmallJs(double *, int);
	void setSmallJCSA(int);
	void setSmallJDip2(int);

	void setAveragedRm3(double);
	double getAveragedRm3(void);

	double getRm3Overridden(void);	// AP: 29-11-20

	void setExpFile(std::string);
	std::string getExpFile(void);

	int getNFields(void);
	void setFields(double *, int);
	void getFields(double *);

	void calculateNMRdata(void);
	void calculateNMRdatainDF(void);
	dvector getNMRdata(void);

	std::string toString(void);

private:

	void makeLarmorFrequencies(void);
	void makeBigJs(void);
	double minFreq(double *, int);
	dcomplex D2MK(int, int, double, double, double);

	std::string probe;
	std::string expFile;
	double deltaCSA;
	double alphaDC, betaDC, gammaDC;
	double alphaDD, betaDD, gammaDD;
	double alphaD, betaD, gammaD;
	int nLarmorFreq; // Total number of Larmor frequencies
	double *larmorFreq, *freqH, *freqX;
	double *smalljs;
	double *smalljCSA;
	double *smalljDip2;
	double averageRm3;
	int nFields;
	double *fields;
	dvector nmrData;
	dvector expNmrData;

	double gyroX;
	double rm3_ov; // Sets a standard, constant r^3 factor

	double rex;

	int probeSet, fieldSet;
	int buildDip2; // For X-H2 probes // 0 = no, 1 = yes 

	dvector bigJdip, bigJcsa;

	bool isSmallJsSet;

};

#endif

