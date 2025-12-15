////////////////////////////////////////////////////////////////////////////////////////////////////////////
// IMPORTANT: cuBLAS ''READS'' MATRICES IN COLUMN-MAJOR FORMAT, WHILE buildAGI.cpp RETURNS
// A, g, AND I MATRICES IN ROW-MAJOR FORMAT. THUS, WHEN USING THE MAGMA LIBRARY, THE
// GPU COPY OF THE MATRICES IS THE TRANSPOSE OF THE MATRIX. THIS FACT IS USED IN:
// 1) CPU-->GPU MATRIX COPY (magma_dsetmatrix): ROW AND COLUMN SIZES ARE INVERTED
// 2) GPU-->CPU MATRIX COPY (magma_dgetmatrix): THE TRANSPOSE OF THE MATRIX IS COLLECTED
// 3) MATRIX-MATRIX MULTIPLICATION (magma_dgemm): WHAT IS DONE IS THE TRANSPOSE OF THE MULTIPLICATION
//    OF THE TRANSPOSED.
// 
// THUS, IF ONE WANTS TO CALCULATE C = alpha*A*B + beta*C, with A(n,m), B(m,k), and C(n,k) THEN
// magma_dsetmatrix(m, n, A, m, d_Atr, m)
// magma_dsetmatrix(k, m, B, k, d_Btr, k)
// magma_dgemm(MagmaNoTrans, MagmaNoTrans, k, n, m, alpha, d_Btr, k, d_Atr, m, beta, d_Ctr, k)
// magma_dgetmatrix(k, n, d_Ctr, k, C, k);
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// WARNING: DO NOT USE MAGMA_DPRINT TO OUTPUT MATRICES SINCE IT OUTPUTS THE TRANSPOSED MATRIX
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include <chrono>

#include <exception>

#include "inputParser.h"
#include "buildAGI.h"
#include "diften.h"
#include "nmr.h"
#include "orderParameter.h"
#include "varcovar.h"
#include "getNH.h"
extern "C" {
#include "cquadpak.h"
}

#define __MAGMA__
#include "cuda.h"
#include "magma_v2.h"
#include "testing.h"
extern "C"{
void dgebal_(const char *job, const magma_int_t *n, double *A, const magma_int_t *lda, magma_int_t *ilo, magma_int_t *ihi, double *scale, magma_int_t *info );
void dgebak_(const char *job, const char *side, const magma_int_t *n, const magma_int_t *ilo, const magma_int_t *ihi, const double *scale, const magma_int_t *m, double *V, const magma_int_t *ldv, magma_int_t *info);
void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info);
void dgeev_(char* jobvl, char* jobvr, int* n, double* a, int* lda, double* wr, double* wi, double* vl, int* ldvl, double* vr, int* ldvr, double* work, int* lwork, int* info);
void dpotrf_(char *uplo,int *jb,double *A,int *lda,int *info);
void dpotri_(char *uplo,int *jb,double *A,int *lda,int *info);
}
#define TPB_1D 512 // Threads per block in 1D geometry
#define TPB_2D 22  // Threads per block in 2D geometry
#include "GPUkernels.cu"

#define CONVERT_C  // this uses automatic differentiation to convert variance-covariance matrix from Cartesian to internal coordinates

/*************************************/
/* Implements the integrand function */
/*************************************/

molecule srb_global_mol;
double *global_q0, *global_dq;
double *global_Kmat;
int *global_aid;
double function(double* x, int *nx)
{
	atom a;
	zmatline z;

	// Z-Matrix line 2
	a = srb_global_mol.getAtom(2);
	z = a.getZMatrixEntry();
	z.d = x[0];
	global_dq[0] = x[0] - global_q0[0];
	srb_global_mol.changeZMatrix(2, z);

	// Z-Matrix line 3
	a = srb_global_mol.getAtom(3);
	z = a.getZMatrixEntry();
	z.d = x[1];
	global_dq[1] = x[1] - global_q0[1];
	z.theta = x[2];
	global_dq[2] = x[2] - global_q0[2];
	srb_global_mol.changeZMatrix(3, z);

	// All the next lines
	int c = 3;
	for (int i = 4; i <= srb_global_mol.getNAtoms(); ++i)
	{
		a = srb_global_mol.getAtom(i);
		z = a.getZMatrixEntry();
		z.d = x[c];
		global_dq[c] = x[c] - global_q0[c];
		c++;
		z.theta = x[c];
		global_dq[c] = x[c] - global_q0[c];
		c++;
		z.phi = x[c];
		global_dq[c] = x[c] - global_q0[c];
		c++;
		srb_global_mol.changeZMatrix(i, z);
	}

	srb_global_mol.buildXYZfromZMatrix();

	// Calculate D200

	a = srb_global_mol.getAtomOld(global_aid[0]);
	vectorOfDoubles rH = srb_global_mol.getAtomPosition(a.id);
	a = srb_global_mol.getAtomOld(global_aid[1]);
	vectorOfDoubles rN = srb_global_mol.getAtomPosition(a.id);

	double dx = rH[0] - rN[0];
	double dy = rH[1] - rN[1];
	double dz = rH[2] - rN[2];
	double rNH = sqrt(dx * dx + dy * dy + dz * dz);
	dz /= rNH;
	dx = 0.5 * (3.0 * dz * dz - 1.0);

	// Calculate Peq

	double peq = 0.0;
	for (int i = 0; i < nx[0]; ++i)
	{
		for (int j = 0; j < nx[0]; ++j)
			peq += global_dq[i] * global_dq[j] * global_Kmat[i * nx[0] + j];
	}
std::cout << peq << std::endl;
	peq = exp(-peq);

	// Return function value
	return(dx * peq);
}

/***********************************************/
/* BUILD HESSIAN MATRIX IN INTERNAL COORDINATES*/
/***********************************************/

void read_hessian(molecule *mol, std::string hessianFile, /*std::string gradientsFile,*/ double *hessInt) 
{
	int natoms = mol->getNAtoms();
	int nk = 3 * natoms;

	int nActiveAtoms = mol->getNActiveAtoms();
	int nq = 3 * nActiveAtoms - 6;

	atom *activeAtoms = new atom[nActiveAtoms];
	mol->getActiveAtoms(activeAtoms);

	double *chessian = new double [nk * nk];
	double *hessian  = new double [nk * nk];

	// Load Force Constants

	std::fstream fileh;
	fileh.open(hessianFile.c_str(), std::ios::in);

	atom atom1, atom2;

	int itmp;
	fileh >> itmp;

	for (int i = 0; i < (nk * nk); ++i)
	{
		fileh >> chessian[i];
	}
	fileh.close();

	// Ensure symmetry (errors can arise from low precision in GROMACS)
	for (int i = 0; i < nk; ++i)
	{
		for (int j = i + 1; j < nk; ++j)
		{
			chessian[i * nk + j] = 0.5 * (chessian[i * nk + j] + chessian[j * nk + i]);
			chessian[j * nk + i] = chessian[i * nk + j];
		}
	}

#ifdef DEBUG
	std::fstream test_chess;
	test_chess.open("test_chessian.dat", std::ios::out);
	
	for (int i = 0; i < nk; ++i)
	{
		for (int j = 0; j < nk; ++j)
			test_chess << std::scientific << chessian[i * nk + j] << "\t";
		test_chess << std::endl;
	}
	test_chess.close();
#endif // DEBUG




////	// Find Eigenvalues and Eigenvectors. They will correspond to the
////	// first eigenv's of small k matrix (which is block diagonal)
////	magma_int_t infoK;
////	double *h_workK;
////	double *lambdaK;
////	magma_int_t lworkK;
////	magma_int_t *iworkK;
////	magma_int_t liworkK;
////	TESTING_CHECK(magma_dmalloc_cpu(&lambdaK, nk));
////
////	double *dKmat;
////	TESTING_CHECK(magma_dmalloc(&dKmat, nk * nk));
////	cudaMemcpy(dKmat, chessian, nk * nk * sizeof(double), cudaMemcpyHostToDevice);
////
////	// Query for workspace sizes
////	double aux_workK[1];
////	magma_int_t aux_iworkK[1];
////	magma_dsyevd_gpu (MagmaVec, MagmaLower, nk, dKmat, nk, lambdaK, chessian, nk, aux_workK, -1, aux_iworkK, -1, &infoK);
////	lworkK = (magma_int_t) aux_workK[0];
////	liworkK = aux_iworkK[0];
////	iworkK = (magma_int_t *) malloc (liworkK * sizeof (magma_int_t));
////	TESTING_CHECK(magma_dmalloc_cpu (&h_workK, lworkK)); // host mem. for workspace
////
////	// Compute the eigenvalues and eigenvectors for a symmetric, real matrix
////	//////////////////////////////////////////////////////////////////////////////////////////////////////////
////	// IMPORTANT: since this is a FORTRAN-based routine, the problem here solved is: Kmat = W^tr Lambda W
////	//            The matrix of eigenvectors is V = W^tr
////	//////////////////////////////////////////////////////////////////////////////////////////////////////////
////	magma_dsyevd_gpu(MagmaVec, MagmaLower, nk, dKmat, nk, lambdaK, chessian, nk, h_workK, lworkK, iworkK, liworkK, &infoK);
////	TESTING_CHECK(infoK);
////
////	cudaMemcpy(chessian, dKmat, nk * nk * sizeof(double), cudaMemcpyDeviceToHost);
////
////	// Free memory
////	magma_free(dKmat);
////	free(iworkK);
////	magma_free_cpu(h_workK);
////
////	std::cout << "* Hessian - Eigenvalues" << std::endl;
////	for (int i = 0; i < nk; ++i)
////		std::cout << i+1 << "  " << lambdaK[i] << std::endl;
////	std::cout << "*****************" << std::endl;
////
////	std::cout << "* Hessian - Eigenvectors" << std::endl;
////	for (int i = 0; i < nk; ++i)
////	{
////		for (int j = 0; j < nk; ++j)
////			std::cout << chessian[i * nk + j] << "\t";
////		std::cout << std::endl;
////	}
////	std::cout << "*****************" << std::endl;

	// Convert to new numbering

	int indi, indj;

	for (int i = 0; i < nk; i += 3)
	{
		atom1 = mol->getAtomOld(i/3 + 1);
		indi = atom1.getID() - 1;

		for (int j = 0; j < natoms; ++j)
		{
			atom2 = mol->getAtomOld(j + 1);
			indj = atom2.getID() - 1;

			hessian[3 * indi * nk + indj * 3]           = chessian[i * nk + j * 3];
			hessian[3 * indi * nk + indj * 3 + 1]       = chessian[i * nk + j * 3 + 1];
			hessian[3 * indi * nk + indj * 3 + 2]       = chessian[i * nk + j * 3 + 2];

			hessian[(3 * indi + 1) * nk + indj * 3]     = chessian[(i + 1) * nk + j * 3];
			hessian[(3 * indi + 1) * nk + indj * 3 + 1] = chessian[(i + 1) * nk + j * 3 + 1];
			hessian[(3 * indi + 1) * nk + indj * 3 + 2] = chessian[(i + 1) * nk + j * 3 + 2];

			hessian[(3 * indi + 2) * nk + indj * 3]     = chessian[(i + 2) * nk + j * 3];
			hessian[(3 * indi + 2) * nk + indj * 3 + 1] = chessian[(i + 2) * nk + j * 3 + 1];
			hessian[(3 * indi + 2) * nk + indj * 3 + 2] = chessian[(i + 2) * nk + j * 3 + 2];
		}
	}

	std::cout << std::endl << std::endl;

#ifdef DEBUG
	std::cout << "* Cartesian Hessian matrix:" << std::endl;
	for (int i = 0; i < nk; ++i)
	{
		for (int j = 0; j < nk; ++j)
			std::cout << std::scientific << hessian[i * nk + j] << "\t";
		std::cout << std::endl;
	}

	std::fstream test_hess;
	test_hess.open("test_hessian.dat", std::ios::out);
	
	for (int i = 0; i < nk; ++i)
	{
		for (int j = 0; j < nk; ++j)
			test_hess << std::scientific << hessian[i * nk + j] << "\t";
		test_hess << std::endl;
	}
	test_hess.close();
#endif // DEBUG

	printf("Now transforms to the internal coordinate reference system....\n ");

	int idxi;
	int lbi;

	double* fdcoi  = new double [nk * nq];
	double* tmpmtx = new double [nk * nq];
	double tmpres;

	vectorOfDoubles dummy(3);

	idxi = 0;

	for (int nai = 2; nai <= natoms; nai++)
	{
		if (mol->getAtom(nai).isActive()) // Include only internal coordinates of active atoms
		{
			lbi = (nai == 2 ? 1 : (nai == 3 ? 2 : 3));
			for (int qi = 1; qi <= lbi; qi++)
			{
				for(int naj = 1; naj <= natoms; naj++)
				{
					dummy = mol->getFirstDerivative(naj, qi, nai);
	
					fdcoi[((naj - 1) * 3 + 0) * nq + idxi] = dummy[0];
					fdcoi[((naj - 1) * 3 + 1) * nq + idxi] = dummy[1];
					fdcoi[((naj - 1) * 3 + 2) * nq + idxi] = dummy[2];
				}
				idxi++;
			}
		}
	}

#ifdef DEBUG
	std::fstream deriv_matrix;
	deriv_matrix.open("Derivative_matrix.dat", std::ios::out);
	
	for (int i = 0; i < nk; ++i)
	{
		for (int j = 0; j < nq; ++j)
		{
			deriv_matrix << std::scientific << fdcoi[i * nq + j] << "\t";
		}
		deriv_matrix << std::endl;
	}
	deriv_matrix.close();
#endif // DEBUG

	for(int i = 0; i < nk; i++)
	{
		for(int j = 0; j < nq; j++)
		{
			tmpres = 0.0;
			for(int k = 0; k < nk; k++)
			{
				tmpres += hessian[i * nk + k] * fdcoi[k * nq + j];
			}
			tmpmtx[i * nq + j] = tmpres;
		}
	}

	for(int i = 0; i < nq; i++)
	{
		for(int j = 0; j < nq;j++)
		{
			tmpres = 0.0;
			for(int k = 0; k < nk; k++)
			{
				tmpres += fdcoi[k * nq + i] * tmpmtx[k * nq + j];
			}
			hessInt[i * nq + j] = tmpres;
		}
	}

#ifdef DEBUG
	std::fstream hessInt_matrix;
	hessInt_matrix.open("HessInt_matrix.dat", std::ios::out);
	
	for (int i = 0; i < nq; ++i)
	{
		for (int j = 0; j < nq; ++j)
		{
			hessInt_matrix << std::scientific << hessInt[i * nq + j] << "\t";
		}
		hessInt_matrix << std::endl;
	}
	hessInt_matrix.close();
#endif // DEBUG

	delete[]activeAtoms;
	delete[]chessian;
	delete[]hessian;
	delete[]fdcoi;
	delete[]tmpmtx;
	return; 
}

double* multiply_matrices(double *A, double *B, magma_trans_t trA, magma_trans_t trB, int dimension)
{
        // NB: DGEMM reads matrices in column-major format. Therefore, instead
        //     of computing C = AB, the calculation C^tr = B^tr A^tr is carried
        //     out just passing to dgemm B and A with NO transposition
        
        magma_queue_t queue;
        magma_queue_create(0, &queue); // Create a MAGMA queue

        double *C = new double[dimension * dimension];

        double *d_A, *d_B, *d_C;

        magma_dmalloc(&d_A, dimension * dimension);
        magma_dmalloc(&d_B, dimension * dimension);
        magma_dmalloc(&d_C, dimension * dimension);

        // Copy A and B to device
        magma_dsetmatrix(dimension, dimension, A, dimension, d_A, dimension, queue);
        magma_dsetmatrix(dimension, dimension, B, dimension, d_B, dimension, queue);

        // Matrix multiplication
        double alpha = 1.0;
        double beta = 0.0;

        magma_dgemm(trB, trA, dimension, dimension, dimension, alpha, d_B, dimension, d_A, dimension, beta, d_C, dimension, queue);

        // Copy result back to std::vector
        magma_dgetmatrix(dimension, dimension, d_C, dimension, C, dimension, queue);

        magma_free(d_A);
        magma_free(d_B);
        magma_free(d_C);
        magma_queue_destroy(queue);

	return C;
}

void covar2hessian(molecule *mol, double *C, double *hessInt, double kBT) 
{
	int natoms = mol->getNAtoms();
	int nk = 3 * natoms;

	int nActiveAtoms = mol->getNActiveAtoms();
	int nq = 3 * nActiveAtoms - 6;

	atom *activeAtoms = new atom[nActiveAtoms];
	mol->getActiveAtoms(activeAtoms);

	double *hessian  = new double [nk * nk];

	// Calculate K = kB T inv(C)

	double *d_hessian;
	magma_dmalloc(&d_hessian, nk * nk);
	cudaMemcpy(d_hessian, C, nk * nk * sizeof(double), cudaMemcpyHostToDevice);

	magma_int_t infoK;
	double *h_workK;
	magma_int_t lworkK;
	magma_int_t *iworkK;
	magma_int_t liworkK;

	double *lambdaC;
	TESTING_CHECK(magma_dmalloc_cpu(&lambdaC, nk));

	// Query for workspace sizes
	double aux_workK[1];
	magma_int_t aux_iworkK[1];
	magma_dsyevd_gpu (MagmaVec, MagmaLower, nk, d_hessian, nk, lambdaC, hessian, nk, aux_workK, -1, aux_iworkK, -1, &infoK);
	lworkK = (magma_int_t) aux_workK[0];
	liworkK = aux_iworkK[0];
	iworkK = (magma_int_t *) malloc (liworkK * sizeof (magma_int_t));
	TESTING_CHECK(magma_dmalloc_cpu (&h_workK, lworkK)); // host mem. for workspace
	
	// Compute the eigenvalues and eigenvectors for a symmetric, real matrix
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// IMPORTANT: since this is a FORTRAN-based routine, the problem here solved is: C = W^tr Lambda W
	//            The matrix of eigenvectors is V = W^tr
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	magma_dsyevd_gpu(MagmaVec, MagmaLower, nk, d_hessian, nk, lambdaC, hessian, nq, h_workK, lworkK, iworkK, liworkK, &infoK);
	TESTING_CHECK(infoK);

	cudaMemcpy(hessian, d_hessian, nk * nk * sizeof(double), cudaMemcpyDeviceToHost);

	// Invert the eigenvalues of C to get those of K
	// Also, multiply by kBT to get correct units
	double *lambdaH = new double[nk * nk];
		for (int i = 0; i < nk; ++i)
	{
		for (int j = 0; j < nk; ++j)
			lambdaH[i * nk + j] = 0.0;
		if (fabs(lambdaC[i]) < 1.0e-10) // This handles null eigenvalues because of roto-translational symmetries
			lambdaH[i * nk + i] = 0.0;
		else
			lambdaH[i * nk + i] = kBT / lambdaC[i];
	}

	// Compute Hessian by pseudo-inversio of C accounting for roto-translational symmetries (see calculation of lambdaH)
	double *tmpM = new double[nk * nk];
	tmpM = multiply_matrices(lambdaH, hessian, MagmaNoTrans, MagmaNoTrans, nk);
	hessian = multiply_matrices(hessian, tmpM, MagmaTrans, MagmaNoTrans, nk);
	delete[]tmpM;
	delete[]lambdaH;
	
//#ifdef DEBUG
	std::fstream KcartFile;
	KcartFile.open("cartesian_K.dat", std::ios::out);
	for (int i = 0; i < nk; ++i)
	{
		for (int j = 0; j < nk; ++j)
			KcartFile << std::scientific << std::setprecision(16) << hessian[i * nk + j] << "\t";
		KcartFile << std::endl;
	}
	KcartFile.close();
//#endif

	printf("Now transforms to the internal coordinate reference system....\n ");

	int idxi;
	int lbi;

	double* fdcoi  = new double [nk * nq];
	double* tmpmtx = new double [nk * nq];
	double tmpres;

	vectorOfDoubles dummy(3);

	idxi = 0;

	for (int nai = 2; nai <= natoms; nai++)
	{
		if (mol->getAtom(nai).isActive()) // Include only internal coordinates of active atoms
		{
			lbi = (nai == 2 ? 1 : (nai == 3 ? 2 : 3));
			for (int qi = 1; qi <= lbi; qi++)
			{
				for(int naj = 1; naj <= natoms; naj++)
				{
					dummy = mol->getFirstDerivative(naj, qi, nai);
	
					fdcoi[((naj - 1) * 3 + 0) * nq + idxi] = dummy[0];
					fdcoi[((naj - 1) * 3 + 1) * nq + idxi] = dummy[1];
					fdcoi[((naj - 1) * 3 + 2) * nq + idxi] = dummy[2];
				}
				idxi++;
			}
		}
	}

//#ifdef DEBUG
	std::fstream deriv_matrix;
	deriv_matrix.open("Derivative_matrix.dat", std::ios::out);
	
	for (int i = 0; i < nk; ++i)
	{
		for (int j = 0; j < nq; ++j)
		{
			deriv_matrix << std::scientific << fdcoi[i * nq + j] << "\t";
		}
		deriv_matrix << std::endl;
	}
	deriv_matrix.close();
//#endif // DEBUG

	for(int i = 0; i < nk; i++)
	{
		for(int j = 0; j < nq; j++)
		{
			tmpres = 0.0;
			for(int k = 0; k < nk; k++)
			{
				tmpres += hessian[i * nk + k] * fdcoi[k * nq + j];
			}
			tmpmtx[i * nq + j] = tmpres;
		}
	}

	for(int i = 0; i < nq; i++)
	{
		for(int j = 0; j < nq; j++)
		{
			tmpres = 0.0;
			for(int k = 0; k < nk; k++)
			{
				tmpres += fdcoi[k * nq + i] * tmpmtx[k * nq + j];
			}
			hessInt[i * nq + j] = tmpres;
		}
	}

//#ifdef DEBUG
	std::fstream hessInt_matrix;
	hessInt_matrix.open("HessInt_matrix.dat", std::ios::out);
	
	for (int i = 0; i < nq; ++i)
	{
		for (int j = 0; j < nq; ++j)
		{
			hessInt_matrix << std::scientific << std::setprecision(16)  << hessInt[i * nq + j] << "\t";
		}
		hessInt_matrix << std::endl;
	}
	hessInt_matrix.close();
//#endif // DEBUG

	delete[]activeAtoms;
	delete[]hessian;
	delete[]fdcoi;
	delete[]tmpmtx;
	return; 
}

/******************/
/* SMALL u MATRIX */
/******************/
// u = U^tr / sqrt(kB T) = sqrt(lambda) * smallk / sqrt(kB T)
void buildSmallu(int nk,  double sqkBT, double *lambda, double *smallk, double *smallu)
{
	dim3 threadsPerBlock, numBlocks;
	threadsPerBlock.x = threadsPerBlock.y = TPB_2D;
	numBlocks.x = 1 + nk / threadsPerBlock.x;
	numBlocks.y = 1 + nk / threadsPerBlock.y;

	double *d_sqrtLambda, *d_smallk, *d_smallu;

	// Build the sqrt of the eigenvalues

	cudaMalloc(&d_sqrtLambda, nk * sizeof(double));
	cudaMemcpy(d_sqrtLambda, lambda, nk * sizeof(double), cudaMemcpyHostToDevice);

	gpu_sqrt_v<<<1 + nk / TPB_1D, TPB_1D>>>(nk, d_sqrtLambda);

	// calculate small u 

	cudaMalloc(&d_smallk, nk * nk * sizeof(double));
	cudaMemcpy(d_smallk, smallk, nk * nk * sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc(&d_smallu, nk * nk * sizeof(double));
	double m = 1. / sqkBT;
	gpu_diagA_times_B<<<numBlocks, threadsPerBlock>>>(nk, m, d_sqrtLambda, d_smallk, d_smallu);

	cudaMemcpy(smallu, d_smallu, nk * nk * sizeof(double), cudaMemcpyDeviceToHost);

	// free memory

	cudaFree(d_sqrtLambda);
	cudaFree(d_smallk);
	cudaFree(d_smallu);

	return;
}

/***********************************************************/
/* Integral function for <rNH^-3> calculation - 1D version */
/***********************************************************/
double rNH_eq, kNH, normNH;
double integrandRm3(double r)
{
	double f = 1.0 / (r * r * r);
	f *= normNH * exp(-0.5 * kNH * (r - rNH_eq) * (r - rNH_eq));
       	return f;	
}

/******************************************************/
/* Integral function for order parameters calculation */
/******************************************************/
#define kToll 500.0
double beta0, k_beta;
double integrandS(double q)
{
	double f = 1.0;
	f *= 0.5 * (3.0 * q * q - 1.0);
	f *= exp(0.5 * k_beta * (3.0 * q * q - 1.0));
      	return f;
}
double integrandNormS(double q)
{
	double f = 1.0;
	f *= exp(0.5 * k_beta * (3.0 * q * q - 1.0));
      	return f;
}

/********/
/* MAIN */
/********/
int main (int argc, char *argv[])
{

	// Setup chrono
	using std::chrono::high_resolution_clock;
	using std::chrono::duration_cast;
	using std::chrono::duration;
	using std::chrono::milliseconds;
	std::ofstream timesFileStr;
	timesFileStr.open("times.dat", std::ios::out);

	auto tTotal1 = high_resolution_clock::now(); // TIMER0: TOTAL EXECUTION TIME

	auto t1 = high_resolution_clock::now(); // TIMER1: INPUT READING AND SETUP
	// Check command line options
	if (argc < 2)
	{
		std::cout << std::endl << "ERROR: an input file must be specified when calling SALEM" << std::endl << std::endl;
		exit(1);
	}
	
	inputParserSrb ip;
	ip.parseInputFile(argv[1]);

	// This controls SALEM to express Cartesian coordinates in MF, dump the Z-matrix and XYZ, and quit
	int stopAtZmatrix = 0;
	if (argc == 3)
		sscanf(argv[2], "%d", &stopAtZmatrix);

	// Configure SeStO
	std::string str;
	char* zmtlibHome;
	zmtlibHome = getenv ("ZMATLIB_HOME");
	if (zmtlibHome == NULL)
	{
		std::cout << std::endl << "ERROR: the ZMATLIB_HOME envinronment variable is not set. Please set the path or use the script salem" << std::endl << std::endl;
		return 1;
	}
	config conf;
	str.assign(zmtlibHome); str.append("/sesto_1.0/sesto");
	conf.setSesto(str);
	str.assign(zmtlibHome); str.append("/data/vdw.dat");
	conf.setVdwFile(str);
	str.assign(zmtlibHome); str.append("/data/aminoacids.dat");
	conf.setAminoAcidsFile(str);
	str.assign(zmtlibHome); str.append("/data/atoms_masses.dat");
	conf.setMassFile(str);

	dim3 threadsPerBlock, numBlocks;

	TESTING_CHECK(magma_init());

	// create queue
	magma_device_t device;
	magma_getdevice(&device);
	magma_queue_t queue;
	magma_queue_create(device, &queue);
	
	// Set the thermal energy
	double kBT = 1.38064852e-23 * ip.getTemperature() * (1.0e20 * 1.0e-30 / 1.660538921e-27);  // kBT in amu A^2 / fs^2

	// Load the molecule
	molecule mol(&conf);
	mol.setIOFileName(ip.getPDBFile());
	mol.setAllActiveAtoms(ip.getAllActive());
	mol.loadMolecule();
	int natoms = mol.getNAtoms();
	int nActiveAtoms = mol.getNActiveAtoms();
	int nx  = 3 * natoms;
	int nxA = 3 * nActiveAtoms;
	int nq  = 3 * nActiveAtoms - 6;
	int rcdDim = nxA - 3;

	// Set the reference atoms
	int *refAtoms = new int[4];
	ip.getRefAtoms(refAtoms);
	mol.setMainDihedralAngle(refAtoms[0], refAtoms[1], refAtoms[2], refAtoms[3]);

	// Build the Z-Matrix
	mol.buildZMatrix();

	// Get the list of active atoms
	atom *activeAtoms = new atom[nActiveAtoms];
	mol.getActiveAtoms(activeAtoms);

	// Express XYZ in AF
	//
	// IMPORTANT: KEEP MOLECULE ATOMS EXPRESSED IN AF!
	//            From this point, no more changes in atoms positions
	//       
	//            If required, coordinates and derivatives in MF
	//	      can be accessed using the proper functions
	//            in the molecule object, where they are calculated on-the-fly
	mol.buildXYZfromZMatrix();

	// Dump the XYZ coordinates in AF
	mol.setIOFileName("rebuild.xyz");
	mol.dumpMolecule();

	if (stopAtZmatrix)
	{
		mol.expressXYZinMF();
		mol.setIOFileName("rebuild_MF.xyz");
		mol.dumpMolecule();
		exit(1);
	}
	auto t2 = high_resolution_clock::now(); // TIMER1
	duration<double, std::milli> ms_double = t2 - t1;
	timesFileStr << "Input parsing and calculation setup: " << std::scientific << ms_double.count() * 1.0e-3 << " s\n";

	auto t3 = high_resolution_clock::now();; // Add additional timer

	// Calculate K matrix
	double *Kmat;
	double *lambdaK; //AP:21-04-20 Line moved from its original place

	TESTING_CHECK(magma_dmalloc_pinned(&Kmat, nq * nq));

	std::string method = ip.getKMatrixMethod();
	if (!method.compare("hessian"))
	{
		t1 = high_resolution_clock::now(); // TIMER2a: PREPATATION AND INVERSION OF HESSIAN IN INTERNAL COORDINATES

		// Load Hessian Matrix
		read_hessian(&mol, ip.getHessianFile(), /*ip.getGradientsFile(),*/ Kmat);
#ifdef DEBUG	
		std::cout << "* Upper block of small k matrix:" << std::endl;
		std::ofstream f_Kmat;
		f_Kmat.open("Kmat.dat", std::ios::out);
		for (int i = 0; i < nq; ++i)
		{
			for (int j = 0; j < nq; ++j)
			{
				std::cout << std::scientific << std::setprecision(5)  << Kmat[i * nq + j] << "\t";
				f_Kmat    << std::scientific << std::setprecision(12) << Kmat[i * nq + j] << "\t";
			}
			std::cout << std::endl;
			f_Kmat    << std::endl;
		}
		std::cout << "*****************" << std::endl;
		f_Kmat.close();
#endif

		t2 = high_resolution_clock::now(); // TIMER2a

		// Find Eigenvalues and Eigenvectors. They will correspond to the
		// first eigenv's of small k matrix (which is block diagonal)
		magma_int_t infoK;
		double *h_workK;
		magma_int_t lworkK;
		magma_int_t *iworkK;
		magma_int_t liworkK;
		TESTING_CHECK(magma_dmalloc_cpu(&lambdaK, nq));

		double *dKmat;
		TESTING_CHECK(magma_dmalloc(&dKmat, nq * nq));
		cudaMemcpy(dKmat, Kmat, nq * nq * sizeof(double), cudaMemcpyHostToDevice);

		// Query for workspace sizes
		double aux_workK[1];
		magma_int_t aux_iworkK[1];
		magma_dsyevd_gpu (MagmaVec, MagmaLower, nq, dKmat, nq, lambdaK, Kmat, nq, aux_workK, -1, aux_iworkK, -1, &infoK);
		lworkK = (magma_int_t) aux_workK[0];
		liworkK = aux_iworkK[0];
		iworkK = (magma_int_t *) malloc (liworkK * sizeof (magma_int_t));
		TESTING_CHECK(magma_dmalloc_cpu (&h_workK, lworkK)); // host mem. for workspace
	
		// Compute the eigenvalues and eigenvectors for a symmetric, real matrix
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		// IMPORTANT: since this is a FORTRAN-based routine, the problem here solved is: Kmat = W^tr Lambda W
		//            The matrix of eigenvectors is V = W^tr
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		magma_dsyevd_gpu(MagmaVec, MagmaLower, nq, dKmat, nq, lambdaK, Kmat, nq, h_workK, lworkK, iworkK, liworkK, &infoK);
		TESTING_CHECK(infoK);

		cudaMemcpy(Kmat, dKmat, nq * nq * sizeof(double), cudaMemcpyDeviceToHost);

		// Free memory
		magma_free(dKmat);
		free(iworkK);
		magma_free_cpu(h_workK);

		// IGNORE THE SIGN OF THE EIGENVALUES, IF REQUESTED //
		//
		//        *** USE THIS OPTION WITH CAUTION ***      //
		if (ip.ignoreNegativeEigs())
		{
			for (int i = nq - 1; i >=0; --i)
				lambdaK[i] = fabs(lambdaK[i]);
		}

//#ifdef DEBUG
		std::cout << "* bigK - Eigenvalues" << std::endl;
		for (int i = 0; i < nq; ++i)
			std::cout << i+1 << "  " << lambdaK[i] << std::endl;
		std::cout << "*****************" << std::endl;
		std::cout << "* bigK - Eigenvectors" << std::endl;
		for (int i = 0; i < nq; ++i)
		{
			for (int j = 0; j < nq; ++j)
				std::cout << Kmat[i * nq + j] << "\t";
			std::cout << std::endl;
		}
		std::cout << "*****************" << std::endl;
//#endif // DEBUG

		t3 = high_resolution_clock::now(); // TIMER2a
		ms_double = t2 - t1;
		timesFileStr << "Preparation of the Hessian matrix in internal coordinates: " << std::scientific << ms_double.count() * 1.0e-3 << " s\n";
		ms_double = t3 - t2;
		timesFileStr << "Calculation of eigvenvalues and eigenvectors of the Hessian matrix: " << std::scientific << ms_double.count() *1.0e-3 << " s\n";

	}
	else // from covariance matrix
	{
		t1 = high_resolution_clock::now(); // TIMER2b: PREPATATION AND INVERSION OF VARIANCE-COVARIANCE MATRIX IN INTERNAL COORDINATES

		// Load the covariance matrix in Cartesian coordinates
		double *cartesianC = new double [nx * nx];
		std::fstream fC;

		fC.open(ip.getCovarFile().c_str(), std::ios::in);
#ifdef CONVERT_C
		for (int i = 0; i < nx; ++i)
		{
			for (int j = 0; j < nx; ++j)
				fC >> cartesianC[i * nx + j];
		}
#else
		double *cartesianC_initial_numbering = new double [nx * nx];
		for (int i = 0; i < nx; ++i)
		{
			for (int j = 0; j < nx; ++j)
				fC >> cartesianC_initial_numbering[i * nx + j];
		}

		// Swap rows and columns according to new atom numbering
		int indi, indj;
		atom atom1, atom2;
		for (int i = 0; i < nx; i += 3)
		{
			atom1 = mol.getAtomOld(i/3 + 1);
			indi = atom1.getID() - 1;

			for (int j = 0; j < nx; j += 3)
			{
				atom2 = mol.getAtomOld(j/3 + 1);
				indj = atom2.getID() - 1;

				cartesianC[3 * indi * nx + 3 * indj]     = cartesianC_initial_numbering[i * nx + j];
				cartesianC[3 * indi * nx + 3 * indj + 1] = cartesianC_initial_numbering[i * nx + j + 1];
				cartesianC[3 * indi * nx + 3 * indj + 2] = cartesianC_initial_numbering[i * nx + j + 2];
	                                                                                                             
				cartesianC[(3 * indi + 1) * nx + 3 * indj]     = cartesianC_initial_numbering[(i + 1) * nx + j];
				cartesianC[(3 * indi + 1) * nx + 3 * indj + 1] = cartesianC_initial_numbering[(i + 1) * nx + j + 1];
				cartesianC[(3 * indi + 1) * nx + 3 * indj + 2] = cartesianC_initial_numbering[(i + 1) * nx + j + 2];

				cartesianC[(3 * indi + 2) * nx + 3 * indj]     = cartesianC_initial_numbering[(i + 2) * nx + j];
				cartesianC[(3 * indi + 2) * nx + 3 * indj + 1] = cartesianC_initial_numbering[(i + 2) * nx + j + 1];
				cartesianC[(3 * indi + 2) * nx + 3 * indj + 2] = cartesianC_initial_numbering[(i + 2) * nx + j + 2];
			}
		}
		delete[]cartesianC_initial_numbering;
#endif
		fC.close();

//#ifdef DEBUG
		fC.open("swapped_C.dat", std::ios::out);
		for (int i = 0; i < nx; ++i)
		{
			for (int j = 0; j < nx; ++j)
				fC << std::scientific << std::setprecision(16) << cartesianC[i * nx + j] << " ";
			fC << std::endl;
		}
		fC.close();
//#endif // DEBUG

#ifdef CONVERT_C
		// Convert to internal coordinates
		double *internalC;
		TESTING_CHECK(magma_dmalloc_pinned(&internalC, nq * nq));
	
		////////////////////////////////////rotateVarcovarInAF(mol.getAtom(1), mol.getAtom(2), mol.getAtom(3), cartesianC, nx); I THINK THAT THE ROTATION IS NOT NEEDED
		varcovar(&mol, activeAtoms, cartesianC, internalC);
	
		delete[] cartesianC; // Free memory
	
//#ifdef DEBUG
		fC.open("internalCovarianceInAF.dat", std::ios::out);
		for (int i = 0; i < nq; ++i)
		{
			for (int j = 0; j < nq; ++j)
				fC << std::scientific << std::setprecision(16) << internalC[i * nq + j] << " ";
			fC << std::endl;
		}
		fC.close();
//#endif // DEBUG

		t2 = high_resolution_clock::now(); // TIMER2b
#else

		// compute Kmat = kBT * inv(cartesianC)

		covar2hessian(&mol, cartesianC, Kmat, kBT);
#endif
		// Find Eigenvalues and Eigenvectors. They will correspond to the
		// first eigenv's of small k matrix (which is block diagonal)

		// Allocate and copy matrix in GPU memory
		double *dKmat;
		TESTING_CHECK(magma_dmalloc(&dKmat, nq * nq));
#ifdef CONVERT_C
		cudaMemcpy(dKmat, internalC, nq * nq * sizeof(double), cudaMemcpyHostToDevice);
#else
		cudaMemcpy(dKmat, Kmat, nq * nq * sizeof(double), cudaMemcpyHostToDevice);
#endif

		magma_int_t infoK;
		double *h_workK;
		magma_int_t lworkK;
		magma_int_t *iworkK;
		magma_int_t liworkK;
		TESTING_CHECK(magma_dmalloc_cpu(&lambdaK, nq));

		// Query for workspace sizes
		double aux_workK[1];
		magma_int_t aux_iworkK[1];
		magma_dsyevd_gpu (MagmaVec, MagmaLower, nq, dKmat, nq, lambdaK, Kmat, nq, aux_workK, -1, aux_iworkK, -1, &infoK);
		lworkK = (magma_int_t) aux_workK[0];
		liworkK = aux_iworkK[0];
		iworkK = (magma_int_t *) malloc (liworkK * sizeof (magma_int_t));
		TESTING_CHECK(magma_dmalloc_cpu (&h_workK, lworkK)); // host mem. for workspace
	
		// Compute the eigenvalues and eigenvectors for a symmetric, real matrix
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		// IMPORTANT: since this is a FORTRAN-based routine, the problem here solved is: Kmat = W^tr Lambda W
		//            The matrix of eigenvectors is V = W^tr
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		magma_dsyevd_gpu(MagmaVec, MagmaLower, nq, dKmat, nq, lambdaK, Kmat, nq, h_workK, lworkK, iworkK, liworkK, &infoK);
		TESTING_CHECK(infoK);

		cudaMemcpy(Kmat, dKmat, nq * nq * sizeof(double), cudaMemcpyDeviceToHost);

#ifdef CONVERT_C		
		// Invert the eigenvalues of C to get those of K
		// Also, multiply by kBT to get correct units
		for (int i = 0; i < nq; ++i)
			lambdaK[i] = kBT / lambdaK[i];
#endif

		// Free memory
		TESTING_CHECK(magma_free(dKmat));
		free(iworkK);
		TESTING_CHECK(magma_free_cpu(h_workK));

		t3 = high_resolution_clock::now(); // TIMER2b
		ms_double = t2 - t1;
		timesFileStr << "Preparation of the var-covar matrix in internal coordinates: " << std::scientific << ms_double.count() * 1.0e-3 << " s\n";
		ms_double = t3 - t2;
		timesFileStr << "Calculation of eigvenvalues and eigenvectors of the var-covar matrix: " << std::scientific << ms_double.count() *1.0e-3 << " s\n";
	} // end K matrix calculation and diagonalization

	////////////////////////////////////////////////////////////////////////////
	// IMPORTANT NOTE: THE mol OBJECT CONTAINS CARTESIAN COORDINATES IN AF    //
	//                 AND HAS METHODS TO CALCULATE DERIVATIVES IN AF *ONLY*  //
	//                                                                        //
	//                 SINCE WE SHALL WORK IN DF, COORDINATES AND DERIVATIVES //
	//                 IN DF ARE BUILT ON THE FLY IN agi OBJECT WHEN CALLING  //
	//                 agi.buildMatrices(), USING COORDINATES AND DERIVATIVES //
	//                 IN AF FROM THE MOLECULE mol OBJECT                     //
	////////////////////////////////////////////////////////////////////////////

	t1 = high_resolution_clock::now(); // TIMER4: FRICTION TENSOR

	// FRICTION TENSOR AND INERTIA MATRIX
	std::cout << "* Calculating friction tensor..." << std::endl;
	diften dite;
	dite.setMolecule(&mol);
	dite.setReff(ip.getReff());
	dite.setC(ip.getC());
	dite.setViscosity(ip.getViscosity() * (1.0e-10 * 1.0e-15 / 1.660538921e-27));  // change to amu / A fs
	dite.setTemperature(ip.getTemperature());
	dite.setHydrodynamicInteractions(ip.getHydrodynamicInteractions());
	double sqkBT = sqrt(kBT);

	dite.buildFrictionInDF(0);
	int nf = 3 + nq; // dim of friction tensor
	double *rcFriction;
	double *extFriction;

	TESTING_CHECK(magma_dmalloc_pinned(&rcFriction, nf * nf));
	TESTING_CHECK(magma_dmalloc_pinned(&extFriction, nxA * nxA));

	dite.getRotoConformationalFrictionInDF(rcFriction);
	dite.getFrictionTensor(extFriction);

#ifdef DEBUG
	std::cout << "******************************************************" << std::endl;
	std::cout << " * Roto-conformational friction tensor " << std::endl;
	std::cout << "******************************************************" << std::endl;
	std::ofstream f_rcf;
	f_rcf.open("RC_friction.dat", std::ios::out);
	for (int i = 0; i < nf; ++i)
	{
		for (int j = 0; j < nf; ++j)
		{
			std::cout  << rcFriction[i * nf + j] * 1.0e15 * 1.660538921e-27 << "\t";
			f_rcf << std::scientific << std::setprecision(12) << rcFriction[i * nf + j] << "\t";
		}
		std::cout  << std::endl;
		f_rcf << std::endl;
	}
	std::cout << "******************************************************" << std::endl;
	f_rcf.close();
#endif

	t2 = high_resolution_clock::now(); // TIMER4

	// AP: 08-09-20 *** BEGIN PRINT DIFFUSION TENSOR *** //
	int INFO;
	double *Diff = new double[nxA * nxA];
	double *RCDiff = new double[rcdDim * rcdDim];
	for (int i = 0; i < nxA * nxA; ++i) Diff[i] = extFriction[i]; // Diff is in Salem scaled units
	char UPLO = 'L';
	dpotrf_(&UPLO, &nxA, Diff, &nxA, &INFO);
	dpotri_(&UPLO, &nxA, Diff, &nxA, &INFO);

	for (int i = 0; i < nxA; ++i)
	{
		for (int j = i + 1; j < nxA; ++j)
		{
			Diff[j * nxA + i] = Diff[i * nxA + j];
		}
	}

	for (int i = 0; i < nf; ++i)
	{
		for (int j = 0; j < nf; ++j)
			RCDiff[i * rcdDim + j] = Diff[(i + 3) * nx + (j + 3)] * kBT; // roto-conformational diffusion tensor stored in Salem scaled units 
	}

//#ifdef DEBUG	
	std::cout << "******************************************************" << std::endl;
	std::cout << " * Roto-conformational diffusion tensor " << std::endl;
	std::cout << "******************************************************" << std::endl;
	std::ofstream f_rcd;
	double kBT_factor = 1.38064852e-23 * ip.getTemperature();
	f_rcd.open("RC_diffusion.dat", std::ios::out);
	for (int i = 0; i < nf; ++i)
	{
		for (int j = 0; j < nf; ++j)
		{
			std::cout  << Diff[(i + 3) * nx + (j + 3)] * kBT_factor  / (1.0e15 * 1.660538921e-27) << "\t"; // Output diffusion tensor in SI units
			f_rcd << std::scientific << std::setprecision(12) << Diff[(i + 3) * nxA + (j + 3)] << "\t";
		}
		std::cout  << std::endl;
		f_rcd << std::endl;
	}
	std::cout << "******************************************************" << std::endl;
	f_rcd.close();
//#endif

	delete[] Diff;

	t3 = high_resolution_clock::now(); // TIMER4
	ms_double = t2 - t1;
	timesFileStr << "Computation of friction tensor: " << std::scientific << ms_double.count() * 1.0e-3 << " s\n";
	ms_double = t3 - t2;
	timesFileStr << "Computation of diffusion tensor (DO WE NEED THIS OR IT WAS ONLY FOR DEBUGGING BY ANDREA?): " << std::scientific << ms_double.count() *1.0e-3 << " s\n";

	// Copy parts of the diffusion tensor in matrices

	double *RDiff = new double[3 * 3];
	RDiff[0 * 3 + 0] = RCDiff[0 * rcdDim + 0];	RDiff[0 * 3 + 1] = RCDiff[0 * rcdDim + 1];	RDiff[0 * 3 + 2] = RCDiff[0 * rcdDim + 2];
	RDiff[1 * 3 + 0] = RCDiff[1 * rcdDim + 0];	RDiff[1 * 3 + 1] = RCDiff[1 * rcdDim + 1];	RDiff[1 * 3 + 2] = RCDiff[1 * rcdDim + 2];
	RDiff[2 * 3 + 0] = RCDiff[2 * rcdDim + 0];	RDiff[2 * 3 + 1] = RCDiff[2 * rcdDim + 1];	RDiff[2 * 3 + 2] = RCDiff[2 * rcdDim + 2];

	int cdDim = rcdDim - 3;
	double *CDiff = new double[cdDim * cdDim];
	double *RCBlockDiff = new double [3 * cdDim];

	for (int i = 0; i < cdDim; ++i)
	{
		for (int j = 0; j < cdDim; ++j)
			CDiff[i * cdDim + j] = RCDiff[(i + 3) * rcdDim + (j + 3)];
	}

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < cdDim; ++j)
			RCBlockDiff[i * cdDim + j] = RCDiff[i * rcdDim + (j + 3)];
	}

	delete[] RCDiff;

	// Diagonalize the rotational part of the diffusion tensor

	magma_int_t infoK;
	magma_int_t infoRDiff;
	double *h_workRDiff;
	magma_int_t lworkRDiff;
	magma_int_t *iworkRDiff;
	magma_int_t liworkRDiff;
	double *lambdaRDiff;
	TESTING_CHECK(magma_dmalloc_cpu(&lambdaRDiff, 3));

	double *dRDiff;
	TESTING_CHECK(magma_dmalloc(&dRDiff, 3 * 3));
	cudaMemcpy(dRDiff, RDiff, 3 * 3 * sizeof(double), cudaMemcpyHostToDevice);

	// Query for workspace sizes
	double aux_workRDiff[1];
	magma_int_t aux_iworkRDiff[1];
	magma_dsyevd_gpu (MagmaVec, MagmaLower, 3, dRDiff, 3, lambdaRDiff, RDiff, 3, aux_workRDiff, -1, aux_iworkRDiff, -1, &infoK);
	lworkRDiff = (magma_int_t) aux_workRDiff[0];
	liworkRDiff = aux_iworkRDiff[0];
	iworkRDiff = (magma_int_t *) malloc (liworkRDiff * sizeof (magma_int_t));
	TESTING_CHECK(magma_dmalloc_cpu (&h_workRDiff, lworkRDiff)); // host mem. for workspace
	
	// Compute the eigenvalues and eigenvectors for a symmetric, real matrix
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// IMPORTANT: since this is a FORTRAN-based routine, the problem here solved is: RDiff = W^tr Lambda W
	//            The matrix of eigenvectors is E = W^tr
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	magma_dsyevd_gpu(MagmaVec, MagmaLower, 3, dRDiff, 3, lambdaRDiff, RDiff, 3, h_workRDiff, lworkRDiff, iworkRDiff, liworkRDiff, &infoRDiff);
	TESTING_CHECK(infoRDiff);

	cudaMemcpy(RDiff, dRDiff, 3 * 3 * sizeof(double), cudaMemcpyDeviceToHost);

	// Free memory
	free(iworkRDiff);
	magma_free_cpu(h_workRDiff);

	// Compute the square root of K matrix
	
	double *dKmat, *dSqrtLambdaK, *dSqLVt, *dSqKmat;
	TESTING_CHECK(magma_dmalloc(&dKmat, nq * nq));    // device memory for Kmat (transposed of eigenvectors of K)
	TESTING_CHECK(magma_dmalloc(&dSqLVt, nq * nq));   // device memory for lambdaK^1/2 V^t 
	TESTING_CHECK(magma_dmalloc(&dSqKmat, nq * nq));  // device memory for Kmat^1/2
	TESTING_CHECK(magma_dmalloc(&dSqrtLambdaK, nq));  // device memory for square root eigenvalues of K

	cudaMemcpy(dKmat, Kmat, nq * nq * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dSqrtLambdaK, lambdaK, nq * sizeof(double), cudaMemcpyHostToDevice);
	
	gpu_sqrt_v<<<1 + nq / TPB_1D, TPB_1D>>>(nq, dSqrtLambdaK);

	threadsPerBlock.x = threadsPerBlock.y = TPB_2D;
	numBlocks.x = 1 + nq / threadsPerBlock.x;
	numBlocks.y = 1 + nq / threadsPerBlock.y;
	double m = 1. / sqkBT;
	gpu_diagA_times_B<<<numBlocks, threadsPerBlock>>>(nq, m, dSqrtLambdaK, dKmat, dSqLVt);

	magma_dgemm(MagmaNoTrans, MagmaTrans, nq, nq, nq, 1.0, dSqLVt, nq, dKmat, nq, 0.0, dSqKmat, nq, queue);

	double *sqKmat = new double[nq * nq];
	cudaMemcpy(sqKmat, dSqKmat, nq * nq * sizeof(double), cudaMemcpyDeviceToHost);

	TESTING_CHECK(magma_free (dKmat));
	TESTING_CHECK(magma_free (dSqrtLambdaK));
	TESTING_CHECK(magma_free (dSqLVt));

	std::ofstream sqKmatFile;
	sqKmatFile.open("sqKmat.dat", std::ios::out);
	for (int i = 0; i < nq; ++i)
	{
		for (int j = 0; j < nq; ++j)
			sqKmatFile << std::scientific << std::setprecision(12) <<  sqKmat[i * nq + j] << "\t";
		sqKmatFile << std::endl;
	}
	sqKmatFile.close();

	// Build the new C,C and R,C blocks of the diffusion tensor

	double *dRCDiff, *dCDiff, *dRCDiff1, *dCDiff1;
	TESTING_CHECK(magma_dmalloc(&dRCDiff,  3 * cdDim)); 
	TESTING_CHECK(magma_dmalloc(&dRCDiff1, 3 * cdDim)); 
	TESTING_CHECK(magma_dmalloc(&dCDiff,  cdDim * cdDim));
	TESTING_CHECK(magma_dmalloc(&dCDiff1, cdDim * cdDim));

	cudaMemcpy(dRCDiff, RCBlockDiff, 3     * cdDim * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dCDiff,  CDiff,       cdDim * cdDim * sizeof(double), cudaMemcpyHostToDevice);
	
	magma_dgemm(MagmaNoTrans, MagmaNoTrans, cdDim, 3, cdDim, 1.0, dSqKmat, cdDim, dRCDiff, cdDim, 0.0, dRCDiff1, cdDim, queue);
	magma_dgemm(MagmaNoTrans, MagmaNoTrans, cdDim, cdDim, cdDim, 1.0, dSqKmat, cdDim, dCDiff, cdDim, 0.0, dCDiff1, cdDim, queue);
	magma_dgemm(MagmaNoTrans, MagmaTrans, cdDim, cdDim, cdDim, 1.0, dCDiff1, cdDim, dSqKmat, cdDim, 0.0, dCDiff, cdDim, queue);

	
	cudaMemcpy(RCBlockDiff, dRCDiff1, 3     * cdDim * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(CDiff,       dCDiff,   cdDim * cdDim * sizeof(double), cudaMemcpyDeviceToHost);

	TESTING_CHECK(magma_free (dRCDiff));
	TESTING_CHECK(magma_free (dRCDiff1));
	TESTING_CHECK(magma_free (dCDiff));
	TESTING_CHECK(magma_free (dCDiff1));
	TESTING_CHECK(magma_free (dSqKmat));

#ifdef DEBUG
	std::ofstream RCDpFile;
	RCDpFile.open("rcdiff_prime.dat", std::ios::out);
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < cdDim; ++j)
			RCDpFile << std::scientific << std::setprecision(12) <<  RCBlockDiff[i * cdDim + j] << "\t";
		RCDpFile << std::endl;
	}
	RCDpFile.close();

	std::ofstream CDpFile;
	CDpFile.open("cdiff_prime.dat", std::ios::out);
	for (int i = 0; i < cdDim; ++i)
	{
		for (int j = 0; j < cdDim; ++j)
			CDpFile << std::scientific << std::setprecision(12) <<  CDiff[i * cdDim + j] << "\t";
		CDpFile << std::endl;
	}
	CDpFile.close();
#endif

	// Diagonalize the conformational part of the diffusion tensor

	magma_int_t infoCDiff;
	double *h_workCDiff;
	magma_int_t lworkCDiff;
	magma_int_t *iworkCDiff;
	magma_int_t liworkCDiff;
	double *lambdaCDiff;
	TESTING_CHECK(magma_dmalloc_cpu(&lambdaCDiff, cdDim));

	TESTING_CHECK(magma_dmalloc(&dCDiff, cdDim * cdDim));
	cudaMemcpy(dCDiff, CDiff, cdDim * cdDim * sizeof(double), cudaMemcpyHostToDevice);

	// Query for workspace sizes
	double aux_workCDiff[1];
	magma_int_t aux_iworkCDiff[1];
	magma_dsyevd_gpu (MagmaVec, MagmaLower, cdDim, dCDiff, cdDim, lambdaCDiff, CDiff, cdDim, aux_workCDiff, -1, aux_iworkCDiff, -1, &infoK);
	lworkCDiff = (magma_int_t) aux_workCDiff[0];
	liworkCDiff = aux_iworkCDiff[0];
	iworkCDiff = (magma_int_t *) malloc (liworkCDiff * sizeof (magma_int_t));
	TESTING_CHECK(magma_dmalloc_cpu (&h_workCDiff, lworkCDiff)); // host mem. for workspace
	
	// Compute the eigenvalues and eigenvectors for a symmetric, real matrix
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// IMPORTANT: since this is a FORTRAN-based routine, the problem here solved is: CDiff = W^tr Lambda W
	//            The matrix of eigenvectors is T = W^tr
	//////////////////////////////////////////////////////////////////////////////////////////////////////////

	magma_dsyevd_gpu(MagmaVec, MagmaLower, cdDim, dCDiff, cdDim, lambdaCDiff, CDiff, cdDim, h_workCDiff, lworkCDiff, iworkCDiff, liworkCDiff, &infoCDiff);
	TESTING_CHECK(infoCDiff);

	cudaMemcpy(CDiff, dCDiff, cdDim * cdDim * sizeof(double), cudaMemcpyDeviceToHost);

	magma_free(dCDiff);
	free(iworkCDiff);
	magma_free_cpu(h_workCDiff);

	// Calculate the new RC block of the diffusion tensor

	TESTING_CHECK(magma_dmalloc(&dRCDiff, 3 * cdDim));
	TESTING_CHECK(magma_dmalloc(&dRCDiff1, 3 * cdDim));
	TESTING_CHECK(magma_dmalloc(&dCDiff, cdDim * cdDim));
	
	cudaMemcpy(dRDiff, RDiff, 3 * 3 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dRCDiff, RCBlockDiff, 3 * cdDim * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dCDiff, CDiff, cdDim * cdDim * sizeof(double), cudaMemcpyHostToDevice);

	magma_dgemm(MagmaTrans, MagmaNoTrans, cdDim, cdDim, cdDim, 1.0, dCDiff, cdDim, dRCDiff, cdDim, 0.0, dRCDiff1, cdDim, queue);


#ifdef DEBUG
	cudaMemcpy(RCBlockDiff, dRCDiff1, 3 * cdDim * sizeof(double), cudaMemcpyDeviceToHost);
	RCDpFile.open("rcdiff_tmp.dat", std::ios::out);
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < cdDim; ++j)
			RCDpFile << std::scientific << std::setprecision(12) <<  RCBlockDiff[i * cdDim + j] << "\t";
		RCDpFile << std::endl;
	}
	RCDpFile.close();
#endif

	magma_dgemm(MagmaNoTrans, MagmaNoTrans, cdDim, 3, 3, 1.0, dRCDiff1, cdDim, dRDiff, 3, 0.0, dRCDiff, cdDim, queue);
	
	cudaMemcpy(RCBlockDiff, dRCDiff, 3 * cdDim * sizeof(double), cudaMemcpyDeviceToHost);

	magma_free(&dRDiff);
	magma_free(&dRCDiff);
	magma_free(&dRCDiff1);
	magma_free(&dCDiff);

	// Output the diffusion tensor
	RCDiff = new double[rcdDim * rcdDim];
	for (int i = rcdDim * rcdDim - 1; i >= 0; --i)
		RCDiff[i] = 0.0;

	RCDiff[0 * rcdDim + 0] = lambdaRDiff[0];
	RCDiff[1 * rcdDim + 1] = lambdaRDiff[1];
	RCDiff[2 * rcdDim + 2] = lambdaRDiff[2];

	for (int i = 0; i < cdDim; ++i)
		RCDiff[(i + 3) * rcdDim + (i + 3)] = lambdaCDiff[i];

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < cdDim; ++j)
		{
			RCDiff[i * rcdDim + (j + 3)] = RCBlockDiff[i * cdDim + j];
			RCDiff[(j + 3) * rcdDim + i] = RCBlockDiff[i * cdDim + j];
		}
	}

	// 1. Diffusion tensor
	std::ofstream rcDiffFile;
	rcDiffFile.open("rcdiff_z.dat", std::ios::out);
	for (int i = 0; i < rcdDim; ++i)
	{
		for (int j = 0; j < rcdDim; ++j)
			rcDiffFile << std::scientific << std::setprecision(12) <<  RCDiff[i * rcdDim + j] << "\t";
		rcDiffFile << std::endl;
	}
	rcDiffFile.close();

	// 2. Matrix to diagonalize the rotational part
	std::ofstream EFile;
	EFile.open("E.dat", std::ios::out);
	for (int i = 0; i < 3; ++ i)
	{
		for (int j = 0; j < 3; ++j)
			EFile << std::scientific << std::setprecision(12) << RDiff[i * 3 + j] << "\t";
		EFile << std::endl;
	}
	EFile.close();

	// 3. Matrix to diagonalize the conformational part
	std::ofstream TFile;
	TFile.open("T.dat", std::ios::out);
	for (int i = 0; i < cdDim; ++ i)
	{
		for (int j = 0; j < cdDim; ++j)
			TFile << std::scientific << std::setprecision(12) << CDiff[i * cdDim + j] << "\t";
		TFile << std::endl;
	}
	TFile.close();

	magma_int_t infoD;
	double *h_workD;
	magma_int_t lworkD;
	magma_int_t *iworkD;
	magma_int_t liworkD;
	double *lambdaD;
	TESTING_CHECK(magma_dmalloc_cpu(&lambdaD, rcdDim));

	TESTING_CHECK(magma_dmalloc(&dRCDiff, rcdDim * rcdDim));
	cudaMemcpy(dRCDiff, RCDiff, rcdDim * rcdDim * sizeof(double), cudaMemcpyHostToDevice);

	// Query for workspace sizes
	double aux_workD[1];
	magma_int_t aux_iworkD[1];
	magma_dsyevd_gpu (MagmaVec, MagmaLower, rcdDim, dRCDiff, rcdDim, lambdaD, RCDiff, rcdDim, aux_workD, -1, aux_iworkD, -1, &infoD);
	lworkD = (magma_int_t) aux_workD[0];
	liworkD = aux_iworkD[0];
	iworkD = (magma_int_t *) malloc (liworkD * sizeof (magma_int_t));
	TESTING_CHECK(magma_dmalloc_cpu (&h_workD, lworkD)); // host mem. for workspace
	
	// Compute the eigenvalues and eigenvectors for a symmetric, real matrix
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// IMPORTANT: since this is a FORTRAN-based routine, the problem here solved is: CDiff = W^tr Lambda W
	//            The matrix of eigenvectors is T = W^tr
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	magma_dsyevd_gpu(MagmaVec, MagmaLower, rcdDim, dRCDiff, rcdDim, lambdaD, RCDiff, rcdDim, h_workD, lworkD, iworkD, liworkD, &infoD);
	TESTING_CHECK(infoD);

	std::cout << std::endl << "EIGENVALUES OF FINAL RC DIFFUSION TENSOR" << std::endl << std::endl;
	for (int i = 0; i < rcdDim; ++i)
		std::cout << lambdaD[i] << std::endl;

	magma_free(dRCDiff);
	free(iworkD);
	magma_free_cpu(h_workD);

	// End of timers
	auto tTotal2 = high_resolution_clock::now(); // TIMER0: TOTAL EXECUTION TIME
	ms_double = tTotal2 - tTotal1;
	timesFileStr << "\n\nTOTAL EXECUTION TIME: " << std::scientific << ms_double.count() * 1.0e-3 / 60. << " min\n\n";
	timesFileStr.close();

	// BYE!	
	return 0;
}

