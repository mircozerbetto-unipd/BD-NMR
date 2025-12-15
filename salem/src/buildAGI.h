#ifndef BUILD_AGI_
#define BUILD_AGI_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "zmat.h"

#ifdef __GPU__ // GPU uses magma/lapack
#define __MAGMA__
#include "cuda.h"
#include "magma_v2.h"
#include "testing.h"
extern "C"{
void dgebal_(const char *job, const magma_int_t *n, double *A, const magma_int_t *lda, magma_int_t *ilo, magma_int_t *ihi, double *scale, magma_int_t *info );
void dgebak_(const char *job, const char *side, const magma_int_t *n, const magma_int_t *ilo, const magma_int_t *ihi, const double *scale, const magma_int_t *m, double *V, const magma_int_t *ldv, magma_int_t *info);
}
#define TPB_1D 512 // Threads per block in 1D geometry
#define TPB_2D 22  // Threads per block in 2D geometry
#else // CPU uses LAPACK
#undef __MAGMA__
typedef struct { double r, i; } F77complex;
extern "C"{
void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info);
void dpotrf_(char *uplo,int *jb,double *A,int *lda,int *info);
void dpotri_(char *uplo,int *jb,double *A,int *lda,int *info);
void dgemm_(char* TRANSA, char* TRANSB, const int* M, const int* N, const int* K, double* alpha, double* A, const int* LDA, double* B, const int* LDB, double* beta, double* C, const int* LDC);
int zgeev_(char* JOBVL, char* JOBVR, int* N, F77complex* A, int* LDA, F77complex* W, F77complex* VL, int* LDVL, F77complex* VR, int* LDVR, F77complex* WORK, int* LWORK, double* RWORK, int* INFO);
}
#endif


class buildAGI {

public:

	buildAGI();
	buildAGI(molecule *);
	virtual ~buildAGI();

	void setMol(molecule *m);
	void buildMatrices(void);

	vectorOfDoubles getAmat(void);
	void getAmat(double *);
	vectorOfDoubles getGmat(void);
	void getGmat(double *);
	vectorOfDoubles getIten(void);
	void getIten(double *);
	vectorOfDoubles getinvIten(void);
	void getinvIten(double *);
	vectorOfDoubles getCoordinatesInCF(void);
	vectorOfDoubles getFirstDerivativesInCF(void);
	
private:

	molecule *mol;
	vectorOfDoubles Amat, Gmat, Iten, InvIten, ca, dca;
};

#endif
