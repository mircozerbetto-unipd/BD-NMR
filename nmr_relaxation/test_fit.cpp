#include <cstdlib>
#include <cstdio>

#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

#include <string>

#include "types.h"

#include "fit_cf.h"

int main(void)
{
	int gsize = 50000;
	dvector tData = dvector(gsize, 0.0);
	dvector gDip1 = dvector(gsize, 0.0);

	std::ifstream ifile("./6_Dip1_acf_pre-fit_probe-2.dat");
	for (int i = 0; i < gsize; ++i)
		ifile >> tData[i] >> gDip1[i];
	ifile.close();
	
	gDip1 = fit_cf(tData, gDip1, 4);

	std::ofstream ofile("./6_Dip1_acf_probe-2.dat");
	for (int i = 0; i < gsize; ++i)
		ofile << tData[i] << "\t" << gDip1[i] << "\n";
	ofile.close();

	return 0;
}
