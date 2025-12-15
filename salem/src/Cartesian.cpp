// HERE ARE THE ROUTINES REQUIRED FOR THE AUTOMATIC
// DIFFERENTIATION OF FUNCTIONS OF THE COORDINATES
// OF THE ATOMS WITH RESPECT TO INTERNAL COORDINATES
//
// NOTE: COORDINATES OF THE ATOMS ARE EXPRESSED IN MF
//

#include "ADEL/ADEL.h"

// coord: 0-based coordinate ID in range 0 <= coord < 3Natoms
// coordinates: x = 0, y = 1, z = 2
// example: coordinate y of atomID 8: coord = (atomID-1) * 3 * Natoms + 1
// result is independent on q1 and q2 since it is just stored in mol
double getCartesianCoordinate(molec* mol, int coord, real q1, real q2)
{
	int atomID = coord / 3;
	int x      = coord % 3;
	
	//vectorOfDoubles r = mol->getAtomPositionInMF(atomID + 1);
	vectorOfDoubles r = mol->getAtomPositionInCustomFrame(atomID + 1);
	//std::cout << std::endl << std::endl << "From-Cartesian: " << atomID << " " << r[0] << " " << r[1] << " " << r[2] << std::endl; //AP:19-04-20//
	switch (x)
	{
		case 0:
			return r[0];
		case 1:
			return r[1];
		case 2:
			return r[2];
		default:
		{
			std::cout << std::endl << std::endl << "ERROR in file " << __FILE__ << ", line " << __LINE__ << ": coordinate id " << coord << " is in wrong range. Should be between 0 and " << 3 * mol->getNAtoms() - 1 << std::endl << std::endl;
			exit(1);
			return 0.0;
		}
	}
}

// coord: 0-based coordinate ID in range 0 <= coord < 3Natoms
// coordinates: x = 0, y = 1, z = 2
// example: coordinate y of atomID 8: coord = (atomID-1) * 3 * Natoms + 1
//
// qid: 0-based internal coordinate ID in range 0 <= coord <= 3Natoms
// internal coordinates: d = 0, theta = 1, phi = 2
// example: coordinate phi of atomID 5: qid = (atomID - 1) * 3 * Natoms + 2
double getCartesianFirstDerivative(molec* mol, int coord, real q1, real q2)
{
	int atomID1 = coord / 3 + 1;
	int x       = coord % 3;

	int atomID2 = (int)q1 / 3 + 1;
	int q       = (int)q1 % 3 + 1;

	//vectorOfDoubles dxdq = mol->getFirstDerivativeInMF(atomID1, q, atomID2);
	vectorOfDoubles dxdq = mol->getFirstDerivativeInCustomFrame(atomID1, q, atomID2);

	switch (x)
	{
		case 0:
			return dxdq[0];
		case 1:
			return dxdq[1];
		case 2:
			return dxdq[2];
		default:
		{
			std::cout << std::endl << std::endl << "ERROR in file " << __FILE__ << ", line " << __LINE__ << ": coordinate id " << coord << " is in wrong range. Should be between 0 and " << 3 * mol->getNAtoms() - 1 << std::endl << std::endl;
			exit(1);
			return 0.0;
		}
	}
}

// coord: 0-based coordinate ID in range 0 <= coord < 3Natoms
// coordinates: x = 0, y = 1, z = 2
// example: coordinate y of atomID 8: coord = (atomID-1) * 3 * Natoms + 1
//
// qid: 0-based internal coordinate ID in range 0 <= coord <= 3Natoms
// internal coordinates: d = 0, theta = 1, phi = 2
// example: coordinate phi of atomID 5: qid = (atomID - 1) * 3 * Natoms + 2
double getCartesianSecondDerivative(molec* mol, int coord, real q1, real q2)
{
	int atomID1 = coord / 3 + 1;
	int x       = coord % 3;

	int atomID2 = (int)q1 / 3 + 1;
	int qa      = (int)q1 % 3 + 1;

	int atomID3 = (int)q2 / 3 + 1;
	int qb      = (int)q2 % 3 + 1;

	//vectorOfDoubles d2xdq2 = mol->getSecondDerivativeInMF(atomID1, qa, atomID2, qb, atomID3);
	vectorOfDoubles d2xdq2 = mol->getSecondDerivativeInCustomFrame(atomID1, qa, atomID2, qb, atomID3);

	switch (x)
	{
		case 0:
			return d2xdq2[0];
		case 1:
			return d2xdq2[1];
		case 2:
			return d2xdq2[2];
		default:
		{
			std::cout << std::endl << std::endl << "ERROR in file " << __FILE__ << ", line " << __LINE__ << ": coordinate id " << coord << " is in wrong range. Should be between 0 and " << 3 * mol->getNAtoms() - 1 << std::endl << std::endl;
			exit(1);
			return 0.0;
		}
	}
}
