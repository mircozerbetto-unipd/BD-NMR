/***********************************************************************************
 * C++OPPS 2.2 - Interpretation of NMR relaxation in proteins                      *
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

#ifndef BDNMR_TYPES_H_
#define BDNMR_TYPES_H_

#include <array>
#include <vector>
#include <complex>

// Real
typedef std::vector <int>            ivector;
typedef std::vector <double>         dvector;
typedef std::array<std::array<double, 3>, 3> dsqmatrix3;

// Complex
typedef std::complex<double>               dcomplex;
typedef std::complex <ivector>             civector;
typedef std::vector <std::complex<double>> cdvector;

// Structures
struct EulerAngles {
	double alpha;
	double beta;
	double gamma;
};

struct Quaternion {
	double q0;
	double q1;
	double q2;
	double q3;
};

#endif
