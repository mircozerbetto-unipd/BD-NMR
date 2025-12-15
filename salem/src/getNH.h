#ifndef GETNH_H
#define GETNH_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "zmat.h"

class NHsites {

public:
	NHsites();
	virtual ~NHsites();

	void locateNHsites(molecule *mol);

	int getNumberOfSites(void);

	// NOTE:
	// i is the site number, 0-based;
	// adID will be the 1-based ID in the **OLD** numeration (i.e., same as input PDB)
	void getAtomsOfSite(int i, int *adID); 

	// NOTE: 1-based indexes in **OLD** numeration (i.e., same as input PDB)
	vectorOfIntegers getHydrogens(void);
	vectorOfIntegers getNitrogens(void);
	vectorOfIntegers getCarbons(void);
	vectorOfIntegers getProlines(void);

private:

	int nSites;

	vectorOfIntegers hydrogens;
	vectorOfIntegers nitrogens;
	vectorOfIntegers carbons;

};

#endif
