// Calculate first and second derivatives of
// Wigner matrices with respect to internal
// coordinates of the molecule
#include "d2k0Gradients.h"

/***************/
/* Constructor */
/***************/
d2k0Grad::d2k0Grad()
{

	return;
}

/*****************/
/* Deconstructor */
/*****************/
d2k0Grad::~d2k0Grad()
{
	return;
}


/***********/
/* METHODS */
/***********/
void d2k0Grad::setMolecule(molecule *m)
{
	mol = m;
	return;
}

void d2k0Grad::setHid(int i)
{
	Hid = i;
	return;
}

void d2k0Grad::setNid(int i)
{
	Nid = i;
	return;
}

void d2k0Grad::setCid(int i)
{
	Cid = i;
	return;
}

void d2k0Grad::setK(int i)
{
	K = i;
	return;
}

void d2k0Grad::setBetaDC(double b)
{
	betaDC = b;
	return;
}

void d2k0Grad::buildEulerAngles(vectorOfDoubles rH, vectorOfDoubles rN, vectorOfDoubles rC)
{
	double n;

	double zD1 = rH[0] - rN[0];
	double zD2 = rH[1] - rN[1];
	double zD3 = rH[2] - rN[2];
	n   = 1.0 / sqrt(zD1 * zD1 + zD2 * zD2 + zD3 * zD3);
	zD1 *= n;
	zD2 *= n;
	zD3 *= n;

	double v1 = rC[0]- rN[0];
	double v2 = rC[1]- rN[1];
	double v3 = rC[2]- rN[2];

	double yD1 = v2 * zD3 - v3 * zD2;
	double yD2 = v3 * zD1 - v1 * zD3;
	double yD3 = v1 * zD2 - v2 * zD1;
	n = 1.0 / sqrt(yD1 * yD1 + yD2 * yD2 + yD3 * yD3);
	yD1 *= n;
	yD2 *= n;
	yD3 *= n;

	double xD1 = yD2 * zD3 - yD3 * zD2;
	double xD2 = yD3 * zD1 - yD1 * zD3;
	double xD3 = yD1 * zD2 - yD2 * zD1;

	// Apply a small beta-tilt to avoid singularities (gimbal lock at beta = 0)
	double deltaB = 1.0e-10;
	double tZ1 = (1.0 - 0.5 * deltaB * deltaB) * zD1 - deltaB * zD3;
	double tZ2 = zD2 - deltaB * zD3;
	double tZ3 = deltaB * zD1 + (1.0 - 0.5 * deltaB * deltaB) * zD3;
	aDq0 = atan2(zD2, zD1);
	bDq0 = acos(zD3);

	double tY3 = deltaB * yD1 + (1.0 - 0.5 * deltaB * deltaB) * yD3;
	double tX3 = deltaB * xD1 + (1.0 - 0.5 * deltaB * deltaB) * xD3;

	double gDq0 = atan2(tY3, -tX3);

	double caD = cos(aDq0), saD = sin(aDq0);
	double cbD = cos(bDq0), sbD = sin(bDq0);
	double cgD = cos(gDq0), sgD = sin(gDq0);

	double Exx = caD * cbD * cgD - saD * sgD;
	double Exy = saD * cbD * cgD + caD * sgD;
	double Exz = -sbD * cgD;
	double Eyx = -caD * cbD * sgD - saD * cgD;
	double Eyy = -saD * cbD * sgD + caD * cgD;
	double Eyz = sbD * sgD;
	double Ezx = caD * sbD;
	double Ezy = saD * sbD;
	double Ezz = cbD;

	v1 = zD1 - rN[0];
	v2 = zD2 - rN[1];
	v3 = zD3 - rN[2];
	xD1 = Exx * v1 + Exy * v2 + Exz * v3;
	xD2 = Eyx * v1 + Eyy * v2 + Eyz * v3;
	xD3 = Ezx * v1 + Ezy * v2 + Ezz * v3;

	cbD = cos(betaDC); sbD = sin(betaDC);

	Exx = cbD;
	Exy = 0.0;
	Exz = -sbD;
	Eyx = 0.0;
	Eyy = 1.0;
	Eyz = 0.0;
	Ezx = sbD;
	Ezy = 0.0;
	Ezz = cbD;

	zD1 = Exx * xD1 + Exy * xD2 + Exz * xD3 + rN[0];
	zD2 = Eyx * xD1 + Eyy * xD2 + Eyz * xD3 + rN[1];
	zD3 = Ezx * xD1 + Ezy * xD2 + Ezz * xD3 + rN[2];

	aCq0 = atan2(zD2, zD1);
	bCq0 = acos(zD3);

	return;
}

void d2k0Grad::buildTaylorExpansionCoefficients(void)
{
	atom aH = mol->getAtom(Hid);
	atom aN = mol->getAtom(Nid);
	atom aC = mol->getAtom(Cid);

	vectorOfIntegers Hchain = aH.getChain();
	vectorOfIntegers Cchain = aC.getChain();

	int nHc = Hchain.size();


/*     TO CANCEL 
	dim = 3 * (nHc - 2) + (Hchain[nHc - 2] > 2 ? 3 : (Hchain[nHc - 2] == 2 ? 2 : (Hchain[nHc - 2] == 0 ? 1 : 0)));
	if (Hchain.size() == Cchain.size())
		dim += Cchain[0] > 2 ? 3 : (Cchain[0] == 2 ? 2 : (Cchain[0] == 0 ? 1 : 0));
*/

	dim = 0;
	for (int i = 0; i < nHc - 2; ++i)
	{
		if (mol->getAtom(Hchain[i] + 1).isActive())
			dim += 3;
	}
	if (mol->getAtom(Hchain[nHc - 2] + 1).isActive())
	{
		if (Hchain[nHc - 2] > 2)
			dim += 3;
		else if (Hchain[nHc - 2] == 2)
			dim += 2;
		else if (Hchain[nHc - 2] == 0)
			dim += 1;
	}
	if (Hchain.size() == Cchain.size() && mol->getAtom(Cchain[0] + 1).isActive())
	{
		if (Cchain[0] > 2)
			dim += 3;
		else if (Cchain[0] == 2)
			dim += 2;
		else if (Cchain[0] == 0)
			dim += 1;
	}

	int dim2 = dim * dim;

	iG = new int[dim];

	g1D = new dcomplex[dim];
	g1C = new dcomplex[dim];

	g2D = new dcomplex[dim2];
	g2C = new dcomplex[dim2];

	////vectorOfDoubles rH = mol->getAtomPositionInMF(aH.id);
	////vectorOfDoubles rN = mol->getAtomPositionInMF(aN.id);
	////vectorOfDoubles rC = mol->getAtomPositionInMF(aC.id);
	vectorOfDoubles rH = mol->getAtomPositionInCustomFrame(aH.id);
	vectorOfDoubles rN = mol->getAtomPositionInCustomFrame(aN.id);
	vectorOfDoubles rC = mol->getAtomPositionInCustomFrame(aC.id);

	double al, be;
	double *alFD, *beFD;
	double *alSD, *beSD;

	alFD = new double[2];
	beFD = new double[2];

	alSD = new double[4];
	beSD = new double[4];

	// Zero order
	wig.setK(K);
	buildEulerAngles(rH, rN, rC);

	wig.setAngles(aDq0, bDq0);
	D2K0q0D = wig.getD2K0();

	wig.setAngles(aCq0, bCq0);
	D2K0q0C = wig.getD2K0();

	// First and second derivatives
	/* Fill iG */
	int qc = 0;
	for (int i = 0; i < nHc; ++i)
	{
		if (i < nHc - 2 && mol->getAtom(Hchain[i] + 1).isActive())
		{
			iG[qc] = (mol->getAtom(Hchain[i] + 1).activeID - 1) * 3 + 0;
			qc++;
			iG[qc] = (mol->getAtom(Hchain[i] + 1).activeID - 1) * 3 + 1;
			qc++;
			iG[qc] = (mol->getAtom(Hchain[i] + 1).activeID - 1) * 3 + 2;
			qc++;
		}
		else if (i == nHc - 2 && mol->getAtom(Hchain[i] + 1).isActive())
		{
			if (Hchain[i] > 2)
			{
				iG[qc] = (mol->getAtom(Hchain[i] + 1).activeID - 1) * 3 + 0;
				qc++;
				iG[qc] = (mol->getAtom(Hchain[i] + 1).activeID - 1) * 3 + 1;
				qc++;
				iG[qc] = (mol->getAtom(Hchain[i] + 1).activeID - 1) * 3 + 2;
				qc++;
			}
			else if (Hchain[i] == 2) // Atom 3 is active by definition
			{
				iG[qc] = Hchain[i] * 3 + 0;
				qc++;
				iG[qc] = Hchain[i] * 3 + 1;
				qc++;
			}
			else if (Hchain[i] == 0) // Atom 1 is active by definition
			{
				iG[qc] = (2 - 1) * 3 + 0; // distance 2-1 is d2 ==> qid = 3
				qc++;
			}
		}
	}
	if (Hchain.size() == Cchain.size() && mol->getAtom(Cchain[0] + 1).isActive())
	{
		if (Cchain[0] > 2)
		{
			iG[dim - 3] = (mol->getAtom(Cchain[0] + 1).activeID - 1) * 3 + 0;
			iG[dim - 2] = (mol->getAtom(Cchain[0] + 1).activeID - 1) * 3 + 1;
			iG[dim - 1] = (mol->getAtom(Cchain[0] + 1).activeID - 1) * 3 + 2;
		}
		else if (Cchain[0] == 2) // Atom 3 is active by definition
		{
			iG[dim - 2] = Cchain[0] * 3 + 0;
			iG[dim - 1] = Cchain[0] * 3 + 1;
		}
		else if (Cchain[0] == 0) // Atom 1 is active by definition
			iG[dim - 1] = (2 - 1) * 3 + 0;
	}

	/* Build derivatives */

	dcomplex Da, Db, Daa, Dbb, Dab;

	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			/* 1. Dipolar */
			eulerD(aH.id, aN.id, aC.id, iG[i], iG[j], mol, &al, alFD, alSD, &be, beFD, beSD);

			wig.setAngles(al, acos(be));
			Da  = wig.getFirstDerivative(ALPHA);
			Db  = wig.getFirstDerivative(BETA);
			Daa = wig.getSecondDerivative(ALPHA, ALPHA);
			Dbb = wig.getSecondDerivative(BETA, BETA);
			Dab = wig.getSecondDerivative(ALPHA, BETA);

			double mult = fabs(fabs(be) - 1.0) < 1.0e-13 ? 0.0 : -1.0 / sqrt(1.0 - be * be);
			beFD[0] *= mult;
			beFD[1] *= mult;
			beSD[1] *= mult * mult;

			if (i == j) g1D[i] = Da * alFD[0] + Db * beFD[0];
			g2D[i * dim + j]  = Daa * alFD[0] * alFD[1] + Da * alSD[1];
			g2D[i * dim + j] += Dbb * beFD[0] * beFD[1] + Db * beSD[1];
			g2D[i * dim + j] += Dab * (alFD[0] * beFD[1] + alFD[1] * beFD[0]);
			g2D[i * dim + j] *= 0.5;

			/* 2. CSA */
/*
			eulerC(aH.id, aN.id, aC.id, iG[i], iG[j], mol, betaDC, &al, alFD, alSD, &be, beFD, beSD);

			wig.setAngles(al, acos(be));
			Da  = wig.getFirstDerivative(ALPHA);
			Db  = wig.getFirstDerivative(COSBETA);
			Daa = wig.getSecondDerivative(ALPHA, ALPHA);
			Dbb = wig.getSecondDerivative(COSBETA, COSBETA);
			Dab = wig.getSecondDerivative(ALPHA, COSBETA);

			if (i == j) g1C[i] = Da * alFD[0] + Db * beFD[0];
			g2C[i * dim + j]  = Daa * alFD[0] * alFD[1] + Da * alSD[1];
			g2C[i * dim + j] += Dbb * beFD[0] * beFD[1] + Db * beSD[1];
			g2C[i * dim + j] += Dab * (alFD[0] * beFD[1] + alFD[1] * beFD[0]);
			g2C[i * dim + j] *= 0.5;
*/
		}
	}
	return;
}

void d2k0Grad::getZeroOrder(int mu, dcomplex* D2K0)
{
	switch (mu)
	{
		case DIP:
		{
			D2K0[0] = D2K0q0D;
			break;
		}
		case CSA:
		{
			D2K0[0] = D2K0q0C;
			break;
		}
		default:
		{
			std::cout << "Error in " << __FILE__ << ", line " << __LINE__ << ": magnetic interation " << mu << " not implemented. Please choose between DIP (or 1), or CSA (or 2)" << std::endl;
			exit(1);
		}
	}
	return;
}

void d2k0Grad::getFirstDerivatives(int mu, dcomplex *g1)
{
	switch (mu)
	{
		case DIP:
		{
			for (int i = 0; i < dim; ++i)
				g1[i]  = g1D[i];
			break;
		}
		case CSA:
		{
			for (int i = 0; i < dim; ++i)
				g1[i]  = g1C[i];
			break;
		}
		default:
		{
			std::cout << "Error in " << __FILE__ << ", line " << __LINE__ << ": magnetic interation " << mu << " not implemented. Please choose between DIP (or 1), or CSA (or 2)" << std::endl;
			exit(1);
		}
	}
	return;
}

void d2k0Grad::getSecondDerivatives(int mu, dcomplex *g2)
{
	switch (mu)
	{
		case DIP:
		{
			for (int i = 0; i < dim * dim; ++i)
				g2[i]  = g2D[i];
			break;
		}
		case CSA:
		{
			for (int i = 0; i < dim * dim; ++i)
				g2[i]  = g2C[i];
			break;
		}
		default:
		{
			std::cout << "Error in " << __FILE__ << ", line " << __LINE__ << ": magnetic interation " << mu << " not implemented. Please choose between DIP (or 1), or CSA (or 2)" << std::endl;
			exit(1);
		}
	}
	return;
}

int d2k0Grad::getArraysMainDim(void)
{
	return dim;
}

void d2k0Grad::getQIdx(int *idx)
{
	for (int i = 0; i < dim; ++i)
		idx[i] = iG[i];
	return;
}

