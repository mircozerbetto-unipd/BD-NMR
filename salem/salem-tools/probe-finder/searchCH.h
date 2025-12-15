#ifndef SEARCH_CH_H
#define SEARCH_CH_H

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

struct CHprobe {
	int Hid, Cid, Xid, Yid;	// Reference atoms ID's 1, 2, 3, 4
	vectorOfDoubles T;	// Translation to a2 (C)
	vectorOfDoubles R;	// Rotation to AF
	int resNumber;		// Residue number of C atom
	std::string resName;	// Residue name of C atom
};

typedef std::vector<CHprobe> CHprobeVec;

CHprobeVec searchCH(molecule mol);

#endif
