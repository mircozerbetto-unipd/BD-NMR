#ifndef SEARCH_NH_H
#define SEARCH_NH_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "zmat.h"
#include "rotomat.h"

struct NHprobe {
	int Hid, Nid, Cid, Oid;	// Reference atoms ID's 1, 2, 3, 4
	vectorOfDoubles T;	// Translation to a2 (N)
	vectorOfDoubles R;	// Rotation to AF
	int resNumber;		// Residue number of N atom
	std::string resName;	// Residue name of N atom
};

typedef std::vector<NHprobe> NHprobeVec;

NHprobeVec searchNH(molecule mol);

#endif
