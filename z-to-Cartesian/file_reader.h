#ifndef FILE_READER_H
#define FILE_READER_H

#include <vector>
#include <string>
#include "zmat.h"

// Function to read a matrix from a file
std::vector<double> read_matrix_from_file(const std::string& filename, int dimension);

// Function to q0
std::vector<double> read_PDB_file(const std::string& filename, int *refAtoms, molecule *mol);

#endif // FILE_READER_H
