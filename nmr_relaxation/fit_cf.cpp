#include "fit_cf.h"
#include <cminpack.h>
#define real __cminpack_real__


/////////////////////////////////////////////////
// Structure to pass correlation-function data //
/////////////////////////////////////////////////

typedef struct  {
	int nData;
	double *x, *y;
} fcndata_t;

//////////////////////////////////////////////
// Calculate the multi-exponential function //
//////////////////////////////////////////////

double mef(int nFit, const double *params, double x)
{
	double f = 0.0;
	double norm = 0.0;

	for (int i = 0; i < nFit/2; ++i)
	{
		f += params[i * 2 + 0] * params[i * 2 + 0] * exp(- x / params[i * 2 + 1]);
		norm += params[i * 2 + 0] * params[i * 2 + 0];
	}
	
	f /= norm;

	return f;
}

//////////////////////////////////////////
// Calculate Jacobian of the deviations //
//////////////////////////////////////////

void mef_jac(int nFit, int nData, const double *x, const double *xdata, int ldfjac, double *fjac)
{
	// Pre-calculate useful data
	double div = 0.0;
	for (int i = 0; i < nFit; i += 2)
	{
		div += x[i] * x[i];
	}
	div = 1.0 / div;

	// Jacobian of fvec (this is why of the minus sign)
	double f;
	for (int i = 0; i < nData; ++i)
	{
		f = mef(nFit, x, xdata[i]);
		for (int j = 0; j < nFit/2; ++j)
		{
			// pre-exponential factor of j-th exponential
			fjac[i + ldfjac * (j * 2 + 0)] = -2.0 * x[j * 2 + 0] * div * (std::exp(-xdata[i] / x[j * 2 + 1]) - f);// / (10.0 * (xdata[i] + 1.0) * (xdata[i] + 1.0));

			// characteristic time of j-th exponential
			fjac[i + ldfjac * (j * 2 + 1)] = -std::exp(-xdata[i] / x[j * 2 + 1]) * xdata[i] * div * x[j * 2 + 0] * x[j * 2 + 0] / (x[j * 2 + 1] * x[j * 2 + 1]) / (10.0 * (xdata[i] + 1.0) * (xdata[i] + 1.0));
			if (xdata[i] < 10.0)
			{
				fjac[i + ldfjac * (j * 2 + 0)] *= 1.0e4/ (xdata[i] + 1.0);
				fjac[i + ldfjac * (j * 2 + 1)] *= 1.0e4/ (xdata[i] + 1.0);
			}

		}
	}
	return;
}

//////////////////////////////////////////
// Function to calculate the deviations //
//////////////////////////////////////////

int fcn(void *p, int nData, int nFit, const double *x, double *fvec, double *fjac, 
	 int ldfjac, int iflag)
{      

	const double *ydata = ((fcndata_t*)p)->y;
	const double *xdata = ((fcndata_t*)p)->x;

	if (iflag == 0)
	{
		/* insert print statements here when nprint is positive. */
		/* if the nprint parameter to lmder is positive, the function is
		   called every nprint iterations with iflag=0, so that the
		   function may perform special operations, such as printing
		   residuals. */
		return 0;
	}

	if (iflag != 2)
	{
		/* compute residuals */
		for (int i = 0; i < nData; ++i)
		{
			fvec[i] = (ydata[i] - mef(nFit, x, xdata[i]));;
			if (xdata[i] < 10.0)
				fvec[i] *= 1.0e4/ (xdata[i] + 1.0);
		}

	}
	else
		/* compute Jacobian */
		mef_jac(nFit, nData, x, xdata, ldfjac, fjac);

	return 0;
}

/////////////////
// Fit routine //
/////////////////

dvector fit_cf(dvector x_data, dvector y_data, int nExp)
{
  	const int nData = x_data.size();
	const int nFit = 2 * nExp;

 	int ldfjac, lwa, info;
	int *ipvt;
	double tol, fnorm;
	double *x, *fvec, *fjac, *wa;

	ldfjac = nData;
	lwa = 5 * nFit + nData;
	wa = new double[lwa];
	ipvt = new int[nFit];
	x = new double[nFit];
	fvec = new double[nData];
	fjac = new double[nData * nFit];

	//////////////////////////////////////////////////
	// Initial parameter guesses (adjust as needed) //
	//////////////////////////////////////////////////

	if (nExp == 1)
	{
		x[0] = 1.0;
		x[1] = x_data.back() * 2.0 / 100.0; // 1 exponential
	}
	else if (nExp == 2)
	{
		x[0] = 0.9;
		x[1] =  x_data.back() * 2.0 / 100.0;
		x[2] = 0.1;
		x[3] = x_data.back() * 2.0 / 1000.0; // 2 exponentials
	}
	else if (nExp == 3)
	{
		x[0] = 0.8;
		x[1] = x_data.back() * 2.0 / 100.0;
		x[2] = 0.1;
		x[3] = x_data.back() * 2.0 / 1000.0;
		x[4] = 0.1;
		x[5] = x_data.back() * 2.0 / 10000.0; // 3 exponentials
	}
	else if (nExp == 4)
	{
		x[0] = 0.6;
		x[1] = x_data.back() * 2.0 / 100.0;
		x[2] = 0.2;
		x[3] = x_data.back() * 2.0 / 1000.0;
		x[4] = 0.1;
		x[5] = x_data.back() * 2.0 / 10000.0;
		x[6] = 0.1;
		x[7] = x_data.back() * 2.0 / 100000.0; // 4 exponentials
	}
	else
	{
		x[0] = 0.7;
		x[1] = x_data.back() * 2.0 / 100.0;
		x[2] = 0.3;
		x[3] = x_data.back() * 2.0 / 1000.0;
		x[4] = 0.2;
		x[5] = x_data.back() * 2.0 / 10000.0;
		x[6] = 0.2;
		x[7] = x_data.back() * 2.0 / 100000.0;
		x[8] = 0.2;
		x[9] = x_data.back() * 2.0 / 1000000.0; // 5 exponentials
	}

	//////////////////////////////////
	// Fit the correlation function //
	//////////////////////////////////

	fcndata_t data;
	data.nData = nData;
	data.x = x_data.data();
	data.y = y_data.data();

	tol = sqrt(__cminpack_func__(dpmpar)(1));
    
	int nDatamin = 2300; //nData;
	info = __cminpack_func__(lmder1)(fcn, &data, nDatamin, nFit, x,
                                        fvec, fjac, ldfjac, tol,
	                                ipvt, wa, lwa);

	fnorm = __cminpack_func__(enorm)(nDatamin, fvec);
	printf("      \n\nfinal l2 norm of the residuals%15.7g\n", (double)fnorm);
	printf("      exit parameter                %10i\n\n", info);
  	tol = __cminpack_func__(dpmpar)(1);

	///////////////////////
	// Print the results //
	///////////////////////

	double norm = 0.0;
	for (int i = 0; i < nFit; i += 2)
		norm += x[i] * x[i];

	std::cout << "Optimized function:" << std::endl << std::endl << "C(t / ns) = ";
	for (int i = 0; i < nFit/2; ++i)
	{
		std::cout << x[2 * i] * x[2 * i] / norm << " * exp(-t/" << x[2 * i + 1] << ")";
		if (i < nExp - 1)
			std::cout << " + ";
	}

	std::cout << std::endl << std::endl;

	dvector fcf = dvector(x_data.size(), 0.0);
	for (int i = 0; i < fcf.size(); ++i)
		fcf[i] = mef(nFit, x, x_data[i]);

    return fcf;
}

