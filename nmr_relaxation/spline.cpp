/***********************************************************************************
 * C++OPPS 1.0 - Interpretation of NMR relaxation in proteins                      *
 * Copyright (C) 2008  Mirco Zerbetto                                              * 
 *                                                                                 *
 * This program is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU General Public License                     *
 * as published by the Free Software Foundation; either version 2                  *
 * of the License, or any later version.                                           *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
 * GNU General Public License for more details.                                    *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License               *
 * along with this program; if not, write to the Free Software                     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. *
 ***********************************************************************************
 * Author: Mirco Zerbetto                                                          *
 * Dipartimento di Scienze Chimiche - Universita' di Padova - Italy                *
 * E-mail: mirco.zerbetto@unipd.it                                                 *
 ***********************************************************************************/

/*
 ============================================================================
 Name        : spline.cpp
 Author      : Mirco Zerbetto
 Version     : 1.0
 Copyright   : 2008 Mirco Zerbetto
 Description : Class to evaluate spline approximation to a function
 ============================================================================
 */

#include <cstdlib>

#include "spline.h"

/***************
 * Constructor *
 ***************/
spline::spline()
{
}

/**************
 * Destructor *
 **************/

spline::~spline()
{
}

/**************
 * Initiators *
 **************/

void spline::init(void)
{
	rank = 0;
	npsp = 0;
	return;
}

void spline::init(int r)
{
	init();
	rank = r;
	return;
}

/****************************************************
 * Set the experimental data to be fitted by spline *
 ****************************************************/

void spline::setExperimentalPoints(dvector x, dvector y)
{
	if (x.size() != y.size())
	{
		std::cout << "ERROR in spline::setExperimentalPoints. Size of x and y must be the same (found size x = " << x.size() << ", size y = " << y.size() << ")" << std::endl;
		exit(0);
	}
	xExp = x;
	yExp = y;
	npsp = x.size();
	yp2 = dvector(npsp,0.0);
	yp1 = ypn = (y[1] - y[0]) / (x[1] - x[0]);
	return;
}

/******************************************************************
 * Retrive all experimental data in an array of v = 2*npsp length *
 * x = 0 --> v.size()/2-1, y = v.size()/2 --> v.size()-1          *
 ******************************************************************/

dvector spline::getExperimentalPoints(void)
{
	if (!npsp)
	{
		std::cout << "ERROR in spline::getExperimentalPoints(void). Data has not been put yet."<<std::endl;
		exit(0);
	}
	dvector v = xExp;
	for (int i = 0; i < npsp; i++) v.push_back(yExp[i]);
	return v;
}

/************************************
 * Retrive x if i = 0 or y if i = 1 *
 ************************************/

dvector spline::getExperimentalPoints(int i)
{
	switch (i)
	{
		case 0:
			return xExp;
		case 1:
			return yExp;
		default:
		{
			std::cout << "ERROR in spline::getExpreimentalPoints(int i). Integer i must be 0 (get x) or 1 (get y)" << std::endl;
			exit(0);
		}
	}
}

/**********************************************
 * Return the yp2 vector obtained from spline *
 **********************************************/

dvector spline::getSplineCoeff(void)
{
	return yp2;
}

/********************************************
 * Return the number of experimental points *
 ********************************************/

int spline::getNPSP(void)
{
	return npsp;
}

/*************************************
 * Calculate the spline coefficients *
 *************************************/

void spline::makeSpline(void)
{
	int i,k;
	double p, qn, sig, un;
	dvector u(npsp,0.0);
	
	if (yp1 > INFTY)
		yp2[0] = 0.0;
	else
	{
		yp2[0] = -0.5;
		u[0] = (3.0 / (xExp[1] - xExp[0])) * ((yExp[1] - yExp[0])/(xExp[1] - xExp[0]) - yp1);
	}
	for (i = 1; i < npsp-1; i++)
	{
		sig = (xExp[i] - xExp[i-1]) / (xExp[i+1] - xExp[i-1]);
		p = 1.0 / (sig * yp2[i-1] + 2.0);
		yp2[i] = (sig - 1.0) * p;
		u[i] = (yExp[i+1] - yExp[i]) / (xExp[i+1] - xExp[i]) - (yExp[i] - yExp[i-1]) / (xExp[i] - xExp[i-1]);
		u[i] = (6.0 * u[i] / (xExp[i+1] - xExp[i-1]) - sig*u[i-1]) * p;
	}
	if (ypn > INFTY)
		qn = un = 0.0;
	else
	{
		qn = 0.5;
		un = (3.0 / (xExp[npsp-1] - xExp[npsp-2])) * (ypn - (yExp[npsp-1] - yExp[npsp-2])/(xExp[npsp-1] - xExp[npsp-2]));
	}
	yp2[npsp-1] = (un - qn * u[npsp-2]) / (qn * yp2[npsp-2] + 1.0);
	for (k = npsp-2; k>=0; k--)
	{
		yp2[k] *= yp2[k+1];
		yp2[k] += u[k];
	}
	return;
}

dvector spline::makeSpline(dvector x, dvector y, double y1, double yn)
{
	xExp = x;
	yExp = y;
	npsp = x.size();
	yp1 = y1;
	ypn = yn;
	yp2 = dvector(npsp,0.0);
	makeSpline();
	return yp2;
}

/******************************************
 * Evaluate spline function for a given x *
 ******************************************/

double spline::splint(double x)
{
	int klo = 0, khi = npsp-1, k;
	double h, invh, b,a;
	while(khi-klo>1)
	{
		k = (khi+klo)>>1;
		if (xExp[k] > x) khi = k;
		else klo = k;
	}
	h = xExp[khi] - xExp[klo];
	invh = 1.0 / h;
	if (fabs(h) < ZERO)
	{
		std::cout << "ERROR in spline::splint(double x, double* y). x = " << x << " is out of experimental range" << std::endl;
		exit(0);
	}
	a = (xExp[khi] - x) * invh;
	b = (x - xExp[klo]) * invh;
	return (a*yExp[klo] + b*yExp[khi] + ((a*a*a-a)*yp2[klo]+(b*b*b-b)*yp2[khi])*h*h*ONE_OVER_SIX);
}

double spline::splint(double x, int klo)
{
	int khi = klo++;
	double h, invh, b, a;
	h = xExp[khi] - xExp[klo];
	invh = 1.0 / h;
	a = (xExp[khi] - x) * invh;
	b = (x - xExp[klo]) * invh;
	return (a*yExp[klo] + b*yExp[khi] + ((a*a*a-a)*yp2[klo]+(b*b*b-b)*yp2[khi])*h*h*ONE_OVER_SIX);
}

double spline::splint(dvector xe, dvector ye, dvector y2, double x)
{
	xExp = xe;
	yExp = ye;
	npsp = xe.size();
	yp2 = y2;
	return splint(x);
}

double spline::splint(dvector xe, dvector ye, dvector y2, double x, int klo)
{
	xExp = xe;
	yExp = ye;
	npsp = xe.size();
	yp2 = y2;
	return splint(x,klo);
}
