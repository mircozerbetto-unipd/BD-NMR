#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "zmat.h"
#include "searchNH.h"
#include "searchCH.h"
#include "searchCH2.h"

int main(int argv, char* argc[])
{

	if (argv < 6)
	{
		std::cout << std::endl << "ERROR: input arguments must be: file.pdb refAtomID1 refAtomID2 refAtomID3 refAtomID4" << std::endl << std::endl;
		return -1;
	}
	
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

	// Load the molecule
	molecule mol(&conf);
	mol.setIOFileName(argc[1]);
	mol.loadMolecule();
	int natoms = mol.getNAtoms();

	// Set the reference atoms
	int *refAtoms = new int[4];
	sscanf(argc[2], "%d", &refAtoms[0]);
	sscanf(argc[3], "%d", &refAtoms[1]);
	sscanf(argc[4], "%d", &refAtoms[2]);
	sscanf(argc[5], "%d", &refAtoms[3]);
	mol.setMainDihedralAngle(refAtoms[0], refAtoms[1], refAtoms[2], refAtoms[3]);
	mol.buildZMatrix();
	mol.buildXYZfromZMatrix();

	// N-H
	NHprobeVec nhp = searchNH(mol);
	std::cout << "Number of NH probes found = " << nhp.size() << std::endl;

	// CHx
	CHprobeVec chp = searchCH(mol);
	std::cout << "Number of CH probes found = " << chp.size() << std::endl;

	CH2probeVec ch2p = searchCH2(mol);
	std::cout << "Number of CH2 probes found = " << ch2p.size() << std::endl;

	return 0;
}

