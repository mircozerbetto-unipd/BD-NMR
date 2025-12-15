// This object contains methods to calculate
// D^2_{K,0} Wigner matrices and all their
// first and second derivatives with respect
// to the angles alpha and beta
#ifndef D2K0_H
#define D2K0_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "types.h"

#define ALPHA   1
#define BETA    2
#define COSBETA 3

class d2k0 {

public:
	d2k0();
	virtual ~d2k0();

	void setAlpha(double);
	void setBeta(double);
	void setAngle(int, double);
	void setAngles(double, double);

	void setK(int);

	dcomplex getD2K0(void);
	dcomplex getFirstDerivative(int);
	dcomplex getSecondDerivative(int, int);

	dcomplex getSmalld2K0(void);
	dcomplex getFirstDerivativeSmalld2k0(void);
	dcomplex getSecondDerivativeSmalld2k0(void);
	dcomplex getFirstDerivativeSmalld2k0CosBeta(void);
	dcomplex getSecondDerivativeSmalld2k0CosBeta(void);

private:

	int k;
	double dk;
	double a, sa, ca;
	double b, sb, cb;
	dcomplex ea;
	dcomplex cee;

	int isSetAlpha, isSetBeta, isSetK;

	void checkAngle(int);
	void checkK(int);
	void checkData(void);
};

#endif // D2K0_H
