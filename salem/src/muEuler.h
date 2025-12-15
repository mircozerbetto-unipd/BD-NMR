//////////////////////////////////////////////////
// This driver calculates MF->muF Euler angles  //
// using AD to get derivatives with respect to  //
// the internal coordinates                     //
//////////////////////////////////////////////////
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cmath>

#include "zmat.h"

/////////////////////////////////////////////////////////////////////////////////
// INPUT:                                                                      //
// idH = 1-based index of atom H to build DF                                   //
// idN = 1-based index of atom N to build DF                                   //
// idC = 1-based index of atom C to build DF                                   //
// q1  = index of first internal coordinate, 0 <= q1 < (3Natoms-6)             //
// q2  = index of second internal coordinate, 0 <= q2 < (3Natoms-6)            //
// mol = molecule bearing both information on Z-Matrix and reference structure //
// OUTPUT:                                                                     //
// alpha = value of angle                                                      //
// fd  = array of first derivatives: fd = [dalpha/dq1  dalpha/dq2]             //
// sd  = array of second derivatives                                           //
//            | d^2alpha/dq1^2    d^2alpha/dq1dq2 |                            //
//       sd = |                                   |                            //
//            | d^2alpha/dq2dq1   d^2alpha/dq2^2  |                            //
//                                                                             //
// Analogous IN/OUT description for the other 2 angles                         //
/////////////////////////////////////////////////////////////////////////////////
void eulerD(int idH, int idN, int idC, int q1, int q2, molecule *mol, double *a, double *afd, double *asd, double *b, double *bfd, double *bsd);
void eulerC(int idH, int idN, int idC, int q1, int q2, molecule *mol, double bDC, double *a, double *afd, double *asd, double *b, double *bfd, double *bsd);
