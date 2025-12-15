#ifndef DIFTEN_H
#define DIFTEN_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>

extern "C" {

#include "cblas.h"

static void
set_entry(double *A, int N, int i, int j, double val)
{
  A[j*N+i] = val;
}

static double
get_entry(const double *A, int N, int i, int j)
{
  return A[j*N+i];
}

static int
dsyevr(char JOBZ, char RANGE, char UPLO, int N,
       double *A, int LDA, double VL, double VU,
       int IL, int IU, double ABSTOL, int *M,
       double *W, double *Z, int LDZ, int *ISUPPZ,
       double *WORK, int LWORK, int *IWORK, int LIWORK)
{
  extern void dsyevr_(char *JOBZp, char *RANGEp, char *UPLOp, int *Np,
                      double *A, int *LDAp, double *VLp, double *VUp,
                      int *ILp, int *IUp, double *ABSTOLp, int *Mp,
                      double *W, double *Z, int *LDZp, int *ISUPPZ,
                      double *WORK, int *LWORKp, int *IWORK, int *LIWORKp,
                      int *INFOp);
  int INFO;
  dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU,
          &IL, &IU, &ABSTOL, M, W, Z, &LDZ, ISUPPZ,
          WORK, &LWORK, IWORK, &LIWORK, &INFO);
  return INFO;
}

static double
dlamch(char CMACH)
{
  extern double dlamch_(char *CMACHp);
  return dlamch_(&CMACH);
}
}

#include "zmat.h"

class diften {

public:

	diften();
	virtual ~diften();

	void setMolecule(molecule *);

	void setReff(double);
	double getReff(void);

	void setC(double);
	double getC(void);

	void setViscosity(double);
	double getViscosity(void);

	void setTemperature(double);
	double getTemperature(void);

	void setHydrodynamicInteractions(int);
	int getHydrodynamicInteractions(void);

	void buildFrictionInDF(int);
	void makeSqrtOfRotConfFriction(void);

	vectorOfDoubles getFrictionInDF(void);
	void getFrictionInDF(double*);
	vectorOfDoubles getRotoConformationalFrictionInDF(void);
	void getRotoConformationalFrictionInDF(double *);
	vectorOfDoubles getRotoConformationalBmatrixInDF(void);
	void getRotoConformationalBmatrixInDF(double *);
	vectorOfDoubles getDiffusionInDF(void);
	vectorOfDoubles getSqrtRotoConformationalFrictionInDF(void);
	void getRotationalFrictionInDF(double **);
	void getConformationalFrictionInDF(double **);
	void getRotoConformationalFrictionInDF(double **);
	void getFrictionTensor(double *);

	void outputFrictionInDF(std::string);
	void outputRotoConformationalFrictionInDF(std::string);
	void outputRotoConformationalBmatrixInDF(std::string, int);
	void outputDiffusionInDF(std::string);
	void outputSqrtRotoConformationalFrictionInDF(std::string);

private:

	void buildBmatrixInDF(void);
	void buildDiffusionInDF(void);

	int frictionBuilt;
	int hydroInt;
	double Reff, C, eta, T, x0;
	double *friction0, *friction, *diffusion, *bMatrix, *rotConfFriction;

	molecule *mol;

};

#endif
