#ifndef INPUT_PARSER_SRB_H
#define INPUT_PARSER_SRB_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "zmat.h"

class inputParserSrb {

public:

	inputParserSrb();
	virtual ~inputParserSrb();

	const char* getPDBFile(void);
	void getRefAtoms(int *ids);
	double getReff(void);
	double getC(void);
	double getViscosity(void);
	double getTemperature(void);
	int getHydrodynamicInteractions(void);
	int getNfreq(void);
	double getDeltafreq(void);
	int getRank(void);

	// Rigid body related methods
	int getRigidBody(void);
	double getDxx(void);
	double getDyy(void);
	double getDzz(void);

	// NMR related methods
	int getNMR(void);
	double getDeltaCSA(void);
	std::string getProbe(void);
	std::string getExpFile(void);
	int getNFields(void);
	void getFields(double *);
	int getHatom(void);
	void setHatom(int);

	// K matrix related methods
	std::string getKMatrixMethod(void);
	std::string getHessianFile(void);
	std::string getGradientsFile(void);
	std::string getCovarFile(void);
	int ignoreNegativeEigs(void);

	// Order parameters
	int makeOrderParameters(void);
	vectorOfIntegers getHid(void);
	vectorOfIntegers getNid(void);
	vectorOfIntegers getCid(void);

	// Active atoms
	int getAllActive(void);

	// Main method
	void parseInputFile(char* inputFileName);

private:

	void checkInputConsistency();

	std::string pdbFileName;
	int *ids; // 1-based ID's of reference atoms
	int nfreq;
	int rank;
	double deltafreq; // fs^-1
	double Reff; // A
	double C;
	double viscosity; // Pa s
	int hydroint; // hydrodynamic interactions

	// Rigid body related quantities
	double temperature; // K
	double Dxx, Dyy, Dzz; // Hz

	// NMR related quantities
	std::string probe;
	double dcsa; // ppm
	std::string exp;
	int nFields;
	double *fields; // in input: MHz, internally converted to Hz
	int CH_Hid; // ID of the second H atom in CH2 probe. This will be 1-based, original PDB numbering

	// K Matrix
	int ignoreNegEigs;
	std::string kMethod, kFile1, kFile2;

	int pdbFileNameFound, idsFound, ReffFound, CFound, viscosityFound, temperatureFound, hydrointFound;
	int nfreqFound, deltafreqFound, rankFound;
	int rigidBodyFound, drrFound;
	int nmrFound, probeFound, dcsaFound, hatomFound, expFound, fieldFound;
	int kMethodFound;
	int ignoreNegEigsFound;

	// Order parameters
	int calculateOrderParameters;
	int nSites;
	vectorOfIntegers Hid, Nid, Cid;

	// Active atoms (coarse-graining)
	int allActiveFound;
	int allActive;


	int itmp;
};
#endif
