#ifndef SEARCH_CH2_H
#define SEARCH_CH2_H

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

struct CH2probe {
	int Hid, Cid, Xid, Yid;	// Reference atoms ID's 1, 2, 3, 4
	int H2id;		// ID of the second H atom
	vectorOfDoubles R;	// Translation to a2 (C)
	vectorOfDoubles T;	// Rotation to AF
	int resNumber;		// Residue number of C atom
	std::string resName;	// Residue name of C atom
};

typedef std::vector<CH2probe> CH2probeVec;

CH2probeVec searchCH2(molecule mol);

#endif
