///////////////////////////////////////////////////////////
// This driver transforms the Variance-Covariance matrix //
// in Cartesian coordinates to internal coordinates      //
// based on the Z-Matrix associated to the molecule      //
// and a reference structure of the molecule             //
///////////////////////////////////////////////////////////
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cmath>

#include "zmat.h"
#include "rotomat.h"

/////////////////////////////////////////////////////////////////////////////////
// The conversion routine                                                      //
// mol = molecule bearing both information on Z-Matrix and reference structure //
// M   = variance-covariance matrix in Cartesian coordinates 3nA x 3nA         //
// C   = variance-covariance matrix in internal coordinates (3nA-6) x (3nA-6)  //
/////////////////////////////////////////////////////////////////////////////////
void varcovar(molecule *mol, atom *activeAtoms, double *M, double *C);
double cov_betabeta(atom At1, atom At2, atom Bt1, atom Bt2, int dimM, double *M, double *beta0);
void rotateVarcovarInAF(atom a1, atom a2, atom a3, double *C, int nq);
