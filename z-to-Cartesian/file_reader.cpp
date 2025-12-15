#include "file_reader.h"
#include <fstream>
#include <iostream>
#include <sstream>

#include "zmat.h"

std::vector<double> read_matrix_from_file(const std::string& filename, int dimension)
{
	std::vector<double> A(dimension * dimension, 0.0);

	std::ifstream file(filename);

	if (!file)
	{
		std::cerr << "Unable to open file: " << filename << std::endl;
		exit(-1);
	}
    
	for (int i = 0; i < dimension; ++i)
	{
		for (int j = 0; j < dimension; ++j)
				file >> A[i * dimension + j];
	}

	file.close();

	return A;
}

std::vector<double> read_PDB_file(const std::string& filename, int *refAtoms, molecule *mol)
{

	mol->setIOFileName(filename);
        mol->setAllActiveAtoms(1);
        mol->loadMolecule();
	mol->setMainDihedralAngle(refAtoms[0], refAtoms[1], refAtoms[2], refAtoms[3]);

        int nAtoms = mol->getNAtoms();
	int dimension = 3 * nAtoms - 6;
	std::vector<double> q0(dimension, 0.0);


        // Build the Z-Matrix
        mol->buildZMatrix();
 
	// Assign q0 values
	atom at;

	at = mol->getAtom(2);
	q0[0] = at.getZMatValue(1);

	at = mol->getAtom(3);
	q0[1] = at.getZMatValue(1);
	q0[2] = at.getZMatValue(2);

   	int pos = 3; 

	for (int i = 4; i <= nAtoms; ++i)
	{
		at = mol->getAtom(i);
		q0[pos] = at.getZMatValue(1);
		q0[pos + 1] = at.getZMatValue(2);
		q0[pos + 2] = at.getZMatValue(3);
		pos += 3;
	}

	return q0;
}

