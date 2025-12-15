#ifndef ORDERPARAMETER_H
#define ORDERPARAMETER_H

#include "d2k0Gradients.h"

#define ONE_OVER_SQRT_TWO (1.0 / sqrt(2.0))

double S20(int Hid, int Nid, int Hcid, molecule* mol, double *invsmallu, int uDim, int nq);

void PISE_TEST(int Hid, int Nid, int Hcid, molecule* mol, int uDim, int nq, int kappa);
#endif
