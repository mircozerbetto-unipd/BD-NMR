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

#ifndef SPLINE_H_
#define SPLINE_H_

#include <iostream>
#include <math.h>
#include <vector>

#include "types.h"

#define INFTY 0.99e30
#define ZERO 1.0e-15
#define ONE_OVER_SIX 1.0/6.0

class spline
{
	public:
		spline();
		virtual ~spline();
		
		void init();
		void init(int);
		void setExperimentalPoints(dvector,dvector);
		dvector getExperimentalPoints(void);
		dvector getExperimentalPoints(int);
		dvector getSplineCoeff(void);
		int getNPSP(void);
		void makeSpline(void);
		dvector makeSpline(dvector,dvector,double,double);
		double splint(double,int);
		double splint(double);
		double splint(dvector,dvector,dvector,double);
		double splint(dvector,dvector,dvector,double,int);
		
	private:
		
		int rank;
		int npsp;
		double yp1,ypn;
		dvector xExp,yExp,yp2;
		
};

#endif /*SPLINE_H_*/
