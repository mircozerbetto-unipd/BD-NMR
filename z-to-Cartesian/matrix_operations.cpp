#include "matrix_operations.h"

std::vector<double> invert_matrix(const std::vector<double>& matrix, int dimension)
{
	std::vector<double>invMatrix(dimension * dimension);
	for (int i = dimension * dimension - 1; i >= 0; --i)
		invMatrix[i] = matrix[i];

	double *d_invMatrix;
	magma_dmalloc(&d_invMatrix, dimension * dimension);
	cudaMemcpy(d_invMatrix, invMatrix.data(), dimension * dimension * sizeof(double), cudaMemcpyHostToDevice);

	int *ipiv = new int[dimension];
	magma_int_t info;
	magma_dgetrf_gpu(dimension, dimension, d_invMatrix, dimension, ipiv, &info);
	if (info != 0)
	{
		std::cerr << "Matrix inversion (magma_dgetrf_gpu) failed with code" << info << "\n";
	        exit(-1);
	}

	magma_int_t lwork = magma_get_dgetri_nb(dimension) * dimension;
	double *dwork;
	magma_dmalloc(&dwork, lwork);

	magma_dgetri_gpu(dimension, d_invMatrix, dimension, ipiv, dwork, lwork, &info);
	if (info != 0)
	{
		std::cerr << "Matrix inversion (magma_dgetri_gpu) failed with code" << info << "\n";
        	exit(-1);
	}

	cudaMemcpy(invMatrix.data(), d_invMatrix, dimension * dimension * sizeof(double), cudaMemcpyDeviceToHost);

	magma_free(d_invMatrix);

	return invMatrix;
}

std::vector<double> multiply_matrices(const std::vector<double>& A, const std::vector<double>& B, int dimension, magma_queue_t queue)
{
	// NB: DGEMM reads matrices in column-major format. Therefore, instead
	//     of computing C = AB, the calculation C^tr = B^tr A^tr is carried
	//     out just passing to dgemm B and A with NO transposition
	
	std::vector<double> C(dimension * dimension, 0.0);

	double *d_A, *d_B, *d_C;

	magma_dmalloc(&d_A, dimension * dimension);
	magma_dmalloc(&d_B, dimension * dimension);
	magma_dmalloc(&d_C, dimension * dimension);

	// Copy A and B to device
	magma_dsetmatrix(dimension, dimension, A.data(), dimension, d_A, dimension, queue);
	magma_dsetmatrix(dimension, dimension, B.data(), dimension, d_B, dimension, queue);

	// Matrix multiplication
	double alpha = 1.0;
	double beta = 0.0;

	magma_dgemm(MagmaNoTrans, MagmaNoTrans, dimension, dimension, dimension, alpha, d_B, dimension, d_A, dimension, beta, d_C, dimension, queue);

	// Copy result back to std::vector
	magma_dgetmatrix(dimension, dimension, d_C, dimension, C.data(), dimension, queue);

	magma_free(d_A);
	magma_free(d_B);
	magma_free(d_C);


    return C;
}

void multiply_matrix_vector(double *d_A, double *d_B, double *d_C, int dimension, magma_queue_t queue)
{
	// NB: DGEMM reads matrices in column-major format. Therefore, instead
	//     of computing C = AB, the calculation C^tr = B^tr A^tr is carried
	//     out just passing to dgemm B and A with NO transposition
	
	std::vector<double> C(dimension, 0.0);

	// Matrix multiplication
	double alpha = 1.0;
	double beta = 0.0;

	magma_dgemm(MagmaNoTrans, MagmaNoTrans, 1, dimension, dimension, alpha, d_B, 1, d_A, dimension, beta, d_C, 1, queue);

	// Copy result back to std::vector
        cudaMemcpy(C.data(), d_C, dimension * sizeof(double), cudaMemcpyDeviceToHost);

	return;
}

