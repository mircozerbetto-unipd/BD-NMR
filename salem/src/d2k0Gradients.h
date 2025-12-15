// Calculate first and second derivatives of
// Wigner matrices with respect to internal
// coordinates of the molecule
#ifndef D2K0_GRADIENTS_H
#define D2K0_GRADIENTS_H

#include "zmat.h"
#include "muEuler.h"
#include "d2k0.h"
#include "types.h"

#define DIP 1
#define CSA 2

class d2k0Grad {

public:
	d2k0Grad();
	virtual ~d2k0Grad();

	void setMolecule(molecule *);
	void setHid(int);
	void setNid(int);
	void setCid(int);
	void setK(int);
	void setBetaDC(double);

	void buildTaylorExpansionCoefficients(void);

	void getZeroOrder(int mu, dcomplex* D2K0);
	void getFirstDerivatives(int mu, dcomplex *g1);
	void getSecondDerivatives(int mu, dcomplex *g2);

	int getArraysMainDim(void);
	void getQIdx(int *);

private:

	molecule *mol;      // object pointing to molecule (which should bear q0 information)
	int Hid, Nid, Cid;  // 1-based ID's of H, N, C atoms
	int K;              // projection of angular momentum on Z-MF
	double betaDC;      // dipolar --> CSA frame beta tilt

	d2k0 wig;
	double aDq0, bDq0, aCq0, bCq0;
	void buildEulerAngles(vectorOfDoubles rH, vectorOfDoubles rN, vectorOfDoubles rC);

	// Dipolar interaction
	dcomplex D2K0q0D;    // value of D2K0 at q0
	dcomplex *g1D, *g2D; // arrays with non-zero first and second derivatives

	// CSA interaction
	dcomplex D2K0q0C;    // value of D2K0 at q0
	dcomplex *g1C, *g2C; // arrays with non-zero first and second derivatives

	int *iG;             // arrays of qi, and (qi,qj) indexes of internal coordinates of non-zero derivatives
	int dim;             // dimensions of arrays: g1(dim), g2(dim * dim), iG1(dim), iG2(2 * dim * dim)
};

#endif // D2K0_GRADIENTS_H
