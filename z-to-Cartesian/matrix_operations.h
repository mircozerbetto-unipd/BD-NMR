#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include <vector>
#include <cuda.h>
#include <cuda_runtime.h>
#include <magma_v2.h>
#include <iostream>

// Function to invert a matrix
std::vector<double> invert_matrix(const std::vector<double>& matrix, int dimension);

// Function to multiply two matrices
std::vector<double> multiply_matrices(const std::vector<double>& A, const std::vector<double>& B, int dimension, magma_queue_t queue);

// Function to multiply matrix by vecotr
void multiply_matrix_vector(double *d_A, double *d_B, double *d_C, int dimension, magma_queue_t queue);

#endif // MATRIX_OPERATIONS_H
