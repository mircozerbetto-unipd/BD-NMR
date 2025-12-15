// This object contains methods to calculate
// D^2_{K,0} Wigner matrices and all their
// first and second derivatives with respect
// to the angles alpha and beta
#include "d2k0.h"

/***************/
/* Constructor */
/***************/
d2k0::d2k0()
{
	isSetAlpha = 0;
	isSetBeta  = 0;
	isSetK     = 0;

	a  = 0.0;
	b  = 0.0;
	dk = 0.0;

	cee = std::complex<double>(0.0, 1.0);

	return;
}

/*****************/
/* Deconstructor */
/*****************/
d2k0::~d2k0()
{
	return;
}

/***********/
/* METHODS */
/***********/
void d2k0::setAlpha(double x)
{
	a = x;
	if (isSetK)
	{
		sa = sin(dk * a);
		ca = cos(dk * a);
	}
	isSetAlpha = 1;
	return;
}

void d2k0::setBeta(double x)
{
	b = x;
	sb = sin(b);
	cb = cos(b);
	isSetBeta = 1;
	return;
}

void d2k0::setAngle(int i, double x)
{
	checkAngle(i);
	switch (i)
	{
		case ALPHA:
		{
			a = x;
			if (isSetK)
			{
				sa = sin(dk * a);
				ca = cos(dk * a);
				ea = std::complex<double>(ca, -sa);
			}
			isSetAlpha = 1;
			break;
		}
		case BETA:
		{
			b = x;
			sb = sin(b);
			cb = cos(b);
			isSetBeta = 1;
			break;
		}
	}
	return;
}

void d2k0::setAngles(double x, double y)
{
	a = x;
	if (isSetK)
	{
		sa = sin(dk * a);
		ca = cos(dk * a);
		ea = std::complex<double>(ca, -sa);
	}
	b = y;
	sb = sin(b);
	cb = cos(b);
	isSetAlpha = 1;
	isSetBeta = 1;
	return;
}

void d2k0::setK(int i)
{
	checkK(i);
	k = i;
	dk = (double)k;
	if (isSetAlpha)
	{
		sa = sin(dk * a);
		ca = cos(dk * a);
		ea = std::complex<double>(ca, -sa);
	}
	isSetK = 1;
	return;
}

dcomplex d2k0::getD2K0(void)
{
	checkData();
	return(ea * getSmalld2K0());
}

dcomplex d2k0::getFirstDerivative(int i)
{
	checkAngle(i);
	checkData();
	dcomplex fd;
	switch (i)
	{
		case ALPHA:
		{
			fd = -dk * cee * getD2K0();
			break;
		}
		case BETA:
		{
			fd = ea * getFirstDerivativeSmalld2k0();
			break;
		}
		case COSBETA:
		{
			fd = ea * getFirstDerivativeSmalld2k0CosBeta();
		}
	}
	return fd;
}

dcomplex d2k0::getSecondDerivative(int i, int j)
{
	checkAngle(i);
	checkAngle(j);
	checkData();
	dcomplex sd;
	switch (i)
	{
		case ALPHA:
		{
			switch(j)
			{
				case ALPHA:
				{
					sd = - dk * dk * getD2K0();
					break;
				}
				case BETA:
				{
					sd = - dk * cee * ea * getFirstDerivativeSmalld2k0();
					break;
				}
				case COSBETA:
				{
					sd = - dk * cee * ea * getFirstDerivativeSmalld2k0CosBeta();
					break;
				}
			}
			break;
		}
		case BETA:
		{
			switch(j)
			{
				case ALPHA:
				{
					sd = - dk * cee * ea * getFirstDerivativeSmalld2k0();
					break;
				}
				case BETA:
				{
					sd = ea * getSecondDerivativeSmalld2k0();
					break;
				}
				case COSBETA:
				{
					std::cout << "ERROR: mixed beta / cosbeta derivative not implemented in " << __FILE__ << " (line " << __LINE__ << ")" << std::endl;
					exit(1);
					break;
				}
			}
			break;
		}
		case COSBETA:
		{
			switch(j)
			{
				case ALPHA:
				{
					sd = - dk * cee * ea * getFirstDerivativeSmalld2k0CosBeta();
					break;
				}
				case BETA:
				{
					std::cout << "ERROR: mixed beta / cosbeta derivative not implemented in " << __FILE__ << " (line " << __LINE__ << ")" << std::endl;
					exit(1);
					break;
				}
				case COSBETA:
				{
					sd = ea * getSecondDerivativeSmalld2k0CosBeta();
					break;
				}
			}
			break;
		}
	}
	return sd;
}

#define SQRT_THREE_OVER_EIGHT (sqrt(3.0 / 8.0))
dcomplex d2k0::getSmalld2K0(void)
{
	checkData();
	double smalld;
	switch (k)
	{
		case -2:
		case  2:
		{
			smalld = SQRT_THREE_OVER_EIGHT * sb * sb;
			break;
		}
		case -1:
		case  1:
		{
			smalld = - dk * 2.0 * SQRT_THREE_OVER_EIGHT * sb * cb;
			break;
		}
		case  0:
		{
			smalld = 0.5 * (3.0 * cb * cb - 1.0);
			break;
		}
	}
	return smalld;
}

dcomplex d2k0::getFirstDerivativeSmalld2k0(void) // Derivative is only with respect to beta
{
	checkData();
	double smalldFD;
	switch (k)
	{
		case -2:
		case  2:
		{
			smalldFD = 2.0 * SQRT_THREE_OVER_EIGHT * sb * cb;
			break;
		}
		case -1:
		case  1:
		{
			smalldFD = - dk * 2.0 * SQRT_THREE_OVER_EIGHT * (2.0 * cb * cb - 1.0);
			break;
		}
		case  0:
		{
			smalldFD = -3.0 * cb * sb;
			break;
		}
	}
	return smalldFD;
}

dcomplex d2k0::getSecondDerivativeSmalld2k0(void) // Derivative is only with respect to beta
{
	checkData();
	double smalldSD;
	switch (k)
	{
		case -2:
		case  2:
		{
			smalldSD = 2.0 * SQRT_THREE_OVER_EIGHT * (2.0 * cb * cb - 1.0);
			break;
		}
		case -1:
		case  1:
		{
			smalldSD = dk * 8.0 * SQRT_THREE_OVER_EIGHT * sb * cb;
			break;
		}
		case  0:
		{
			smalldSD = -3.0 * (2.0 * cb * cb - 1.0);
			break;
		}
	}
	return smalldSD;
}

dcomplex d2k0::getFirstDerivativeSmalld2k0CosBeta(void) // Derivative is only with respect to beta
{
	checkData();
	double smalldFD;
	switch (k)
	{
		case -2:
		case  2:
		{
			smalldFD = - 2.0 * SQRT_THREE_OVER_EIGHT * cb;
			break;
		}
		case -1:
		case  1:
		{
			if (fabs(cb - 1.0) > 1.0e-13)
			{
				double r = cb / sb;
				smalldFD = - dk * 2.0 * SQRT_THREE_OVER_EIGHT * (sb - cb * r);
			}
			else
				smalldFD = 0.0;
			break;
		}
		case  0:
		{
			smalldFD = 3.0 * cb;
			break;
		}
	}
	return smalldFD;
}

dcomplex d2k0::getSecondDerivativeSmalld2k0CosBeta(void) // Derivative is only with respect to beta
{
	checkData();
	double smalldSD;
	switch (k)
	{
		case -2:
		case  2:
		{
			smalldSD = -2.0 * SQRT_THREE_OVER_EIGHT;
			break;
		}
		case -1:
		case  1:
		{
			if (fabs(cb - 1.0) > 1.0e-13)
			{
				double r = cb / sb;
				smalldSD = dk * 8.0 * SQRT_THREE_OVER_EIGHT * r * (r * r + 3.0);
			}
			else
				smalldSD = 0.0;
			break;
		}
		case  0:
		{
			smalldSD = 3.0;
			break;
		}
	}
	return smalldSD;
}
/**********/
/* Checks */
/**********/
void d2k0::checkAngle(int i)
{
	if (i != ALPHA && i != BETA && i != COSBETA)
	{
		std::cout << std::endl << "ERROR in file " << __FILE__ << ", line " << __LINE__ << ": wrong angle index " << i << ". Please select betwee 1 (alpha), 2 (beta), or 3 (cos(beta))." << std::endl << std::endl;
		exit(1);
	}
	return;
}

void d2k0::checkK(int i)
{
	if (i < -2 || i > 2)
	{
		std::cout << std::endl << "ERROR in file " << __FILE__ << ", line " << __LINE__ << ": the k index must be in the range [-2, 2], while its actual value, " << i << ", is out of range." << std::endl << std::endl;
		exit(1);
	}
}

void d2k0::checkData(void)
{
	int nerr = 0;
	if (!isSetAlpha)
	{
		std::cout << "ERROR in d2k0 object (" << __FILE__ << ", " << __LINE__ << ": alpha angle has not been set." << std::endl;
		nerr++;
	}
	if (!isSetBeta)
	{
		std::cout << "ERROR in d2k0 object (" << __FILE__ << ", " << __LINE__ << ": beta angle has not been set." << std::endl;
		nerr++;
	}
	if (!isSetK)
	{
		std::cout << "ERROR in d2k0 object (" << __FILE__ << ", " << __LINE__ << ": K index has not been set." << std::endl;
		nerr++;
	}
	if (nerr) exit(1);
	return;
}

