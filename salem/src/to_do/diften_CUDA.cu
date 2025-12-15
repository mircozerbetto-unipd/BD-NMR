#include "diften.h"

__global__ void gpu_fill_matrix(int n, double *M)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < n && c < n)
	{
		if (r > c)
			M[r * n + c] = M[c * n + r];
	}
	return;
}

__global__ void gpu_build_unconstrained_diffusion(int na, int nd, double C, double Reff, double *cm, double *Ra, double *dUD)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	if (i < na && j < na)
	{
		if (i == j)
		{
			/* Diagonal tensor */

			double mult = cm[i]; /* R0 / atoms[i].R0 */ // Commented part of the line are to be used (and adapted) if atoms have different hydrodynamic radii

			dUD[(i * 3 + 0) * nd + (i * 3 + 0)] = mult;
			dUD[(i * 3 + 0) * nd + (i * 3 + 1)] = 0.0;
			dUD[(i * 3 + 0) * nd + (i * 3 + 2)] = 0.0;
	
			dUD[(i * 3 + 1) * nd + (i * 3 + 0)] = 0.0;
			dUD[(i * 3 + 1) * nd + (i * 3 + 1)] = mult;
			dUD[(i * 3 + 1) * nd + (i * 3 + 2)] = 0.0;
	
			dUD[(i * 3 + 2) * nd + (i * 3 + 0)] = 0.0;
			dUD[(i * 3 + 2) * nd + (i * 3 + 1)] = 0.0;
			dUD[(i * 3 + 2) * nd + (i * 3 + 2)] = mult;
		}
		else
		{
	    		/* Out of diagonal tensors */

			double RX = Ra[i * 3 + 0] - Ra[j * 3 + 0];
			double RY = Ra[i * 3 + 1] - Ra[j * 3 + 1];
			double RZ = Ra[i * 3 + 2] - Ra[j * 3 + 2];
			double R = sqrt(RX * RX + RY * RY + RZ * RZ);

			double sumR0  = Reff + Reff;  /*atoms[i].R0 + atoms[j].R0;*/ 
			double sumR02 = sumR0 * sumR0 - 2.0 * Reff * Reff; /*atoms[i].R0 * atoms[j].R0;*/
			double avgR0 = 0.5 * sumR0;
	 
			double mult, iso, aniso;   
			if (R >= sumR0)
			{
				mult = C * Reff / (8.0 * R);
				iso = 1.0 + sumR02 / (3.0 * R * R);
				aniso = ( 1.0 - sumR02 / (R * R) ) / (R * R);
			}
			else
			{
				mult = C * Reff / (6.0 * avgR0);
				iso= 1.0 - 9.0 * R / (32.0 * avgR0);
				aniso = 3.0 / (32.0 * R * avgR0);
			}

			mult *= 0.5 * (cm[i] + cm[j]);

			dUD[(i * 3 + 0) * nd + (j * 3 + 0)] = mult * (iso + aniso * RX * RX);
			dUD[(i * 3 + 0) * nd + (j * 3 + 1)] = mult * aniso * RX * RY;
			dUD[(i * 3 + 0) * nd + (j * 3 + 2)] = mult * aniso * RX * RZ;

			dUD[(i * 3 + 1) * nd + (j * 3 + 0)] = mult * aniso * RY * RX;
			dUD[(i * 3 + 1) * nd + (j * 3 + 1)] = mult * (iso + aniso * RY * RY);
			dUD[(i * 3 + 1) * nd + (j * 3 + 2)] = mult * aniso * RY * RZ;

			dUD[(i * 3 + 2) * nd + (j * 3 + 0)] = mult * aniso * RZ * RX;
			dUD[(i * 3 + 2) * nd + (j * 3 + 1)] = mult * aniso * RZ * RY;
			dUD[(i * 3 + 2) * nd + (j * 3 + 2)] = mult * (iso + aniso * RZ * RZ);
	  	}

	}
	return;
}

/***************/
/* Constructor */
/***************/

diften::diften()
{
	frictionBuilt = 0;
	useHydrodynamicsInteractions = 1; // On by default
}


/*****************/
/* Deconstructor */
/*****************/

diften::~diften()
{
	return;
}

/********************/
/* Add the molecule */
/********************/

void diften::setMolecule(molecule *m)
{
	mol = m;

	// Init arrays (6 + nq) x (6 + nq)

	int n  = 3 * mol->getNAtoms();
	int na = 3 * mol->getNActiveAtoms();

	friction  = new double[na * na];
	diffusion = new double[na * na];
	bMatrix   = new double[n  * na]; // B is, in general, rectangular

	return;
}

/********************/
/* Fluid parameters */
/********************/

// R / A
void diften::setReff(double r)
{
	Reff = r;
	return;
}

double diften::getReff(void)
{
	return Reff;
}

// Boundary conditions
void diften::setC(double b)
{
	C = b;
	return;
}

double diften::getC(void)
{
	return C;
}

// viscosity / Pa s
void diften::setViscosity(double n)
{	
	eta = n; // Internally use A as unit for lengths
	return;
}

double diften::getViscosity(void)
{
	return (eta ); // Gives back the viscosity in Pa s
}

// T / K
void diften::setTemperature(double a)
{
	T = a;
	return;
}

double diften::getTemperature(void)
{
	return T;
}

// Friction correction multipliers
void diften::setCsiMult(double *m, int n)
{
	csiMult.clear();
	csiMult = vectorOfDoubles(n, 0.0);
	for (int i = 0; i < n; ++i)
		csiMult[i] = m[i];
	return;
}

/*******************************/
/* Method to build the Bmatrix */
/*******************************/

void diften::buildBmatrixInMF(void)
{
	int idj;
	int n  = mol->getNAtoms();
	int na = mol->getNActiveAtoms();
	int m  = 3 * na;
	int nq = m - 6;
	int aid, a2id;
	double M = 0.0, mi;
	atom at;

	vectorOfDoubles d(3), r(3);

	// Get the position of the center of mass

	mol->calculateCenterOfMass();
	vectorOfDoubles rCM = mol->getCenterOfMass();

	// Calculate derivatives of CM

	atom *activeAtoms = new atom[na];
	mol->getActiveAtoms(activeAtoms);

	vectorOfDoubles dRCM(nq * 3);

	// q = d2
	dRCM[0 * 3 + 0] = 0.0;
	dRCM[0 * 3 + 1] = 0.0;
	dRCM[0 * 3 + 2] = 0.0;
	for (int i = 1; i <= mol->getNAtoms(); ++i)
	{ 
		at = mol->getAtom(i);
		mi = at.getMass();
		M += mi;
		d = mol->getFirstDerivative(i, 1, 2);
		dRCM[0 * 3 + 0] += mi * d[0];
		dRCM[0 * 3 + 1] += mi * d[1];
		dRCM[0 * 3 + 2] += mi * d[2];
	}
	M = 1.0 / M;
	dRCM[0 * 3 + 0] *= M;
	dRCM[0 * 3 + 1] *= M;
	dRCM[0 * 3 + 2] *= M;

	// q = d3
	dRCM[1 * 3 + 0] = 0.0;
	dRCM[1 * 3 + 1] = 0.0;
	dRCM[1 * 3 + 2] = 0.0;
	for (int i = 1; i <= mol->getNAtoms(); ++i)
	{ 
		at = mol->getAtom(i);
		mi = at.getMass();
		d = mol->getFirstDerivative(i, 1, 3);
		dRCM[1 * 3 + 0] += mi * d[0];
		dRCM[1 * 3 + 1] += mi * d[1];
		dRCM[1 * 3 + 2] += mi * d[2];
	}
	dRCM[1 * 3 + 0] *= M;
	dRCM[1 * 3 + 1] *= M;
	dRCM[1 * 3 + 2] *= M;

	// q = theta3
	dRCM[2 * 3 + 0] = 0.0;
	dRCM[2 * 3 + 1] = 0.0;
	dRCM[2 * 3 + 2] = 0.0;
	for (int i = 1; i <= mol->getNAtoms(); ++i)
	{ 
		at = mol->getAtom(i);
		mi = at.getMass();
		d = mol->getFirstDerivative(i, 2, 3);
		dRCM[2 * 3 + 0] += mi * d[0];
		dRCM[2 * 3 + 1] += mi * d[1];
		dRCM[2 * 3 + 2] += mi * d[2];
	}
	dRCM[2 * 3 + 0] *= M;
	dRCM[2 * 3 + 1] *= M;
	dRCM[2 * 3 + 2] *= M;

	// All other internal coordinates

	for (int j = 3; j < nq; j += 3)
	{

		idj = activeAtoms[3 + j / 3 - 1].getID();

		dRCM[j * 3 + 0] = 0.0;
		dRCM[j * 3 + 1] = 0.0;
		dRCM[j * 3 + 2] = 0.0;

		dRCM[(j + 1) * 3 + 0] = 0.0;
		dRCM[(j + 1) * 3 + 1] = 0.0;
		dRCM[(j + 1) * 3 + 2] = 0.0;

		dRCM[(j + 2) * 3 + 0] = 0.0;
		dRCM[(j + 2) * 3 + 1] = 0.0;
		dRCM[(j + 2) * 3 + 2] = 0.0;

		for (int i = 1; i <= mol->getNAtoms(); ++i)
		{ 
			at = mol->getAtom(i);
			mi = at.getMass();

			d = mol->getFirstDerivative(i, 1, idj);
			dRCM[j * 3 + 0] += mi * d[0];
			dRCM[j * 3 + 1] += mi * d[1];
			dRCM[j * 3 + 2] += mi * d[2];

			d = mol->getFirstDerivative(i, 2, idj);
			dRCM[(j + 1) * 3 + 0] += mi * d[0];
			dRCM[(j + 1) * 3 + 1] += mi * d[1];
			dRCM[(j + 1) * 3 + 2] += mi * d[2];

			d = mol->getFirstDerivative(i, 3, idj);
			dRCM[(j + 2) * 3 + 0] += mi * d[0];
			dRCM[(j + 2) * 3 + 1] += mi * d[1];
			dRCM[(j + 2) * 3 + 2] += mi * d[2];
		}

		dRCM[j * 3 + 0] *= M;
		dRCM[j * 3 + 1] *= M;
		dRCM[j * 3 + 2] *= M;

		dRCM[(j + 1) * 3 + 0] *= M;
		dRCM[(j + 1) * 3 + 1] *= M;
		dRCM[(j + 1) * 3 + 2] *= M;

		dRCM[(j + 2) * 3 + 0] *= M;
		dRCM[(j + 2) * 3 + 1] *= M;
		dRCM[(j + 2) * 3 + 2] *= M;
	}

	//
	// Build the B matrix
	//
	//
	// Apply rotation to express the tensor in the frame that diagonalizes
	// the inertia tensor assuming the latter INDIPENDENT ON INTERNAL COORDINATES.
	// If the last condition is not fulfilled, the buildBmatrixInMF() code
	// must be modified

	mol->calculateInertiaTensorInCM();
	matrix3D Itens = mol->getInertiaTensor();

	jacobi diag;
 	diag.setMaxStep(500);
	diag.setMatrix(Itens);
	diag.diagonalize();
	diag.reorder();
	matrix3D Emat = diag.getEigenVectors3D();

	double CAx[9], CA[3], EC[3], dCA[3];
	double EI[9];
	EI[0] = Emat.xx; EI[1] = Emat.xy; EI[2] = Emat.xz;
	EI[3] = Emat.yx; EI[4] = Emat.yy; EI[5] = Emat.yz;
	EI[6] = Emat.zx; EI[7] = Emat.zy; EI[8] = Emat.zz;

	for (int i = 0; i < 3 * n; i += 3)
	{
		// Translational part

		bMatrix[(i + 0) * m + 0] = 1.0;
		bMatrix[(i + 0) * m + 1] = 0.0;
		bMatrix[(i + 0) * m + 2] = 0.0;
		bMatrix[(i + 1) * m + 0] = 0.0;
		bMatrix[(i + 1) * m + 1] = 1.0;
		bMatrix[(i + 1) * m + 2] = 0.0;
		bMatrix[(i + 2) * m + 0] = 0.0;
		bMatrix[(i + 2) * m + 1] = 0.0;
		bMatrix[(i + 2) * m + 2] = 1.0;

		// Rotational part (Omega = 0, 0, 0)

		aid = 1 + i / 3;

		r = mol->getAtomPosition(aid);
	
		CA[0] = r[0] - rCM[0];
		CA[1] = r[1] - rCM[1];
		CA[2] = r[2] - rCM[2];

		EC[0] = EI[0] * CA[0] + EI[1] * CA[1] + EI[2] * CA[2];
		EC[1] = EI[3] * CA[0] + EI[4] * CA[1] + EI[5] * CA[2];
		EC[2] = EI[6] * CA[0] + EI[7] * CA[1] + EI[8] * CA[2];

		CAx[0] =  0.0;
		CAx[1] =  EC[2];
		CAx[2] = -EC[1];
		CAx[3] = -EC[2];
		CAx[4] =  0.0;
		CAx[5] =  EC[0];
		CAx[6] =  EC[1];
		CAx[7] = -EC[0];
		CAx[8] =  0.0;

		bMatrix[(i + 0) * m + 3] = CAx[0];
		bMatrix[(i + 0) * m + 4] = CAx[1];
		bMatrix[(i + 0) * m + 5] = CAx[2];
		bMatrix[(i + 1) * m + 3] = CAx[3];
		bMatrix[(i + 1) * m + 4] = CAx[4];
		bMatrix[(i + 1) * m + 5] = CAx[5];
		bMatrix[(i + 2) * m + 3] = CAx[6];
		bMatrix[(i + 2) * m + 4] = CAx[7];
		bMatrix[(i + 2) * m + 5] = CAx[8];
		
		// Configurational part

		// q = d2
		d = mol->getFirstDerivative(aid, 1, 2);

		dCA[0] = d[0] - dRCM[0 * 3 + 0];
		dCA[1] = d[1] - dRCM[0 * 3 + 1];
		dCA[2] = d[2] - dRCM[0 * 3 + 2];

		EC[0] = EI[0] * dCA[0] + EI[1] * dCA[1] + EI[2] * dCA[2];
		EC[1] = EI[3] * dCA[0] + EI[4] * dCA[1] + EI[5] * dCA[2];
		EC[2] = EI[6] * dCA[0] + EI[7] * dCA[1] + EI[8] * dCA[2];

		bMatrix[(i + 0) * m + 6] = EC[0];
		bMatrix[(i + 1) * m + 6] = EC[1];
		bMatrix[(i + 2) * m + 6] = EC[2];

		// q = d3
		d = mol->getFirstDerivative(aid, 1, 3);

		dCA[0] = d[0] - dRCM[1 * 3 + 0];
		dCA[1] = d[1] - dRCM[1 * 3 + 1];
		dCA[2] = d[2] - dRCM[1 * 3 + 2];

		EC[0] = EI[0] * dCA[0] + EI[1] * dCA[1] + EI[2] * dCA[2];
		EC[1] = EI[3] * dCA[0] + EI[4] * dCA[1] + EI[5] * dCA[2];
		EC[2] = EI[6] * dCA[0] + EI[7] * dCA[1] + EI[8] * dCA[2];

		bMatrix[(i + 0) * m + 7] = EC[0];
		bMatrix[(i + 1) * m + 7] = EC[1];
		bMatrix[(i + 2) * m + 7] = EC[2];

		// q = theta3
		d = mol->getFirstDerivative(aid, 2, 3);

		dCA[0] = d[0] - dRCM[2 * 3 + 0];
		dCA[1] = d[1] - dRCM[2 * 3 + 1];
		dCA[2] = d[2] - dRCM[2 * 3 + 2];

		EC[0] = EI[0] * dCA[0] + EI[1] * dCA[1] + EI[2] * dCA[2];
		EC[1] = EI[3] * dCA[0] + EI[4] * dCA[1] + EI[5] * dCA[2];
		EC[2] = EI[6] * dCA[0] + EI[7] * dCA[1] + EI[8] * dCA[2];

		bMatrix[(i + 0) * m + 8] = EC[0];
		bMatrix[(i + 1) * m + 8] = EC[1];
		bMatrix[(i + 2) * m + 8] = EC[2];

		for (int j = 9; j < m; j += 3)
		{
			a2id = activeAtoms[j / 3].getID();

			// d
			d = mol->getFirstDerivative(aid, 1, a2id);

			dCA[0] = d[0] - dRCM[(j - 6 + 0) * 3 + 0];
			dCA[1] = d[1] - dRCM[(j - 6 + 0) * 3 + 1];
			dCA[2] = d[2] - dRCM[(j - 6 + 0) * 3 + 2];

			EC[0] = EI[0] * dCA[0] + EI[1] * dCA[1] + EI[2] * dCA[2];
			EC[1] = EI[3] * dCA[0] + EI[4] * dCA[1] + EI[5] * dCA[2];
			EC[2] = EI[6] * dCA[0] + EI[7] * dCA[1] + EI[8] * dCA[2];

			bMatrix[(i + 0) * m + (j + 0)] = EC[0];
			bMatrix[(i + 1) * m + (j + 0)] = EC[1];
			bMatrix[(i + 2) * m + (j + 0)] = EC[2];

			// theta
			d = mol->getFirstDerivative(aid, 2, a2id);

			dCA[0] = d[0] - dRCM[(j - 6 + 1) * 3 + 0];
			dCA[1] = d[1] - dRCM[(j - 6 + 1) * 3 + 1];
			dCA[2] = d[2] - dRCM[(j - 6 + 1) * 3 + 2];

			EC[0] = EI[0] * dCA[0] + EI[1] * dCA[1] + EI[2] * dCA[2];
			EC[1] = EI[3] * dCA[0] + EI[4] * dCA[1] + EI[5] * dCA[2];
			EC[2] = EI[6] * dCA[0] + EI[7] * dCA[1] + EI[8] * dCA[2];

			bMatrix[(i + 0) * m + (j + 1)] = EC[0];
			bMatrix[(i + 1) * m + (j + 1)] = EC[1];
			bMatrix[(i + 2) * m + (j + 1)] = EC[2];

			// phi
			d = mol->getFirstDerivative(aid, 3, a2id);

			dCA[0] = d[0] - dRCM[(j - 6 + 2) * 3 + 0];
			dCA[1] = d[1] - dRCM[(j - 6 + 2) * 3 + 1];
			dCA[2] = d[2] - dRCM[(j - 6 + 2) * 3 + 2];

			EC[0] = EI[0] * dCA[0] + EI[1] * dCA[1] + EI[2] * dCA[2];
			EC[1] = EI[3] * dCA[0] + EI[4] * dCA[1] + EI[5] * dCA[2];
			EC[2] = EI[6] * dCA[0] + EI[7] * dCA[1] + EI[8] * dCA[2];

			bMatrix[(i + 0) * m + (j + 2)] = EC[0];
			bMatrix[(i + 1) * m + (j + 2)] = EC[1];
			bMatrix[(i + 2) * m + (j + 2)] = EC[2];
		}
	}

	return;
}

/*********************************/
/* Methods to build the Friction */
/*********************************/

void diften::setUseHydrodynamicsInteractions(int h)
{
	useHydrodynamicsInteractions = h;
	return;
}

int diften::getUseHydrodynamicsInteractions(void)
{
	return useHydrodynamicsInteractions;
}

void diften::buildUnconstrainedDiffusion(void)
{

//#ifdef __GPU__ // TODO: check what's wrong
//
//	int na = mol->getNAtoms();
//	int nat3 = 3 * na;
//	unconstrainedDiffusion = new double [nat3 * nat3];
//	vectorOfDoubles Ri;
//	double *hRA, *dRA;
//	double *hCM, *dCM;
//
//	hRA = new double[nat3];
//	hCM = new double[na];
//
//	TESTING_CHECK(cudaMalloc(&dRA, nat3));
//	TESTING_CHECK(cudaMalloc(&dCM, na));
//
//	for (int i = 0; i < mol->getNAtoms(); ++i)
//	{
//		Ri = mol->getAtomPosition(i + 1);
//		hRA[i * 3 + 0] = Ri[0];
//		hRA[i * 3 + 1] = Ri[1];
//		hRA[i * 3 + 2] = Ri[2];
//		hCM[i]         = csiMult[mol->getAtom(i + 1).getResidueNumber() - 1];
//	}
//
//	cudaMemcpy(dRA, hRA, nat3 * sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(dCM, hCM, na * sizeof(double), cudaMemcpyHostToDevice);
//
//	delete[] hRA;
//	delete[] hCM;
//
//	double *dUD;
//	TESTING_CHECK(cudaMalloc(&dUD, nat3 * nat3 * sizeof(double)));
//	
//	dim3 threadsPerBlock, numBlocks;
//	threadsPerBlock.x = threadsPerBlock.y = TPB_2D;
//	numBlocks.x = 1 + na / threadsPerBlock.x;
//	numBlocks.y = 1 + na / threadsPerBlock.y;
//	gpu_build_unconstrained_diffusion<<<numBlocks, threadsPerBlock>>>(na, nat3, C, Reff, dCM, dRA, dUD);
//
//	cudaMemcpy(unconstrainedDiffusion, dUD, nat3 * nat3 * sizeof(double), cudaMemcpyDeviceToHost);
//
//	TESTING_CHECK(cudaFree(dRA));
//	TESTING_CHECK(cudaFree(dCM));
//	TESTING_CHECK(cudaFree(dUD));
//
//	for (int i = 0; i < nat3; ++i)
//	{
//		for (int j = 0; j < nat3; ++j)
//			std::cout << i << ", " << j << ") " << unconstrainedDiffusion[i * nat3 + j] << std::endl;
//	}
//	exit(1);
//
//#else // END GPU CODE

	int na = mol->getNAtoms();
	int nat3 = 3 * na;
	double RX, RY, RZ, R;
	double sumR0, sumR02, avgR0;
	double mult, iso, aniso, cmI, cmJ;
	unconstrainedDiffusion = new double [nat3 * nat3];
	vectorOfDoubles Ri, Rj;

	for (int i = 0; i < na; ++i)
	{
		// Diagonal tensor

		cmI = csiMult[mol->getAtom(i + 1).getResidueNumber() - 1];

		mult = 1.0 / cmI; // R0 / atoms[i].R0 // Commented part of the line are to be used (and adapted) if atoms have different hydrodynamic radii

		unconstrainedDiffusion[(i * 3 + 0) * nat3 + (i * 3 + 0)] = mult;
		unconstrainedDiffusion[(i * 3 + 0) * nat3 + (i * 3 + 1)] = 0.0;
		unconstrainedDiffusion[(i * 3 + 0) * nat3 + (i * 3 + 2)] = 0.0;
	
		unconstrainedDiffusion[(i * 3 + 1) * nat3 + (i * 3 + 0)] = 0.0;
		unconstrainedDiffusion[(i * 3 + 1) * nat3 + (i * 3 + 1)] = mult;
		unconstrainedDiffusion[(i * 3 + 1) * nat3 + (i * 3 + 2)] = 0.0;
	
		unconstrainedDiffusion[(i * 3 + 2) * nat3 + (i * 3 + 0)] = 0.0;
		unconstrainedDiffusion[(i * 3 + 2) * nat3 + (i * 3 + 1)] = 0.0;
		unconstrainedDiffusion[(i * 3 + 2) * nat3 + (i * 3 + 2)] = mult;

	    	Ri = mol->getAtomPosition(i + 1);

		for (int j= i + 1; j < na; ++j)
		{
	    		// Out of diagonal tensors

			cmJ = csiMult[mol->getAtom(j + 1).getResidueNumber() - 1];
	
		    	Rj = mol->getAtomPosition(j + 1);
	    
			RX = Ri[0] - Rj[0];
			RY = Ri[1] - Rj[1];
			RZ = Ri[2] - Rj[2];
			R = sqrt(RX * RX + RY * RY + RZ * RZ);

			sumR0  = cmI * Reff + cmJ * Reff;  //atoms[i].R0 + atoms[j].R0; 
			sumR02 = sumR0 * sumR0 - 2.0 * cmI * Reff * cmJ * Reff; //atoms[i].R0 * atoms[j].R0;
			avgR0 = 0.5 * sumR0;
	    
			if (R >= sumR0)
			{
				mult = C * Reff / (8.0 * R);
				iso = 1.0 + sumR02 / (3.0 * R * R);
				aniso = ( 1.0 - sumR02 / (R * R) ) / (R * R);
			}
			else
			{
				mult = C * Reff / (6.0 * avgR0);
				iso= 1.0 - 9.0 * R / (32.0 * avgR0);
				aniso = 3.0 / (32.0 * R * avgR0);
			}

			//mult *= 0.5 * (cmI + cmJ);

			unconstrainedDiffusion[(i * 3 + 0) * nat3 + (j * 3 + 0)] = unconstrainedDiffusion[(j * 3 + 0) * nat3 + (i * 3 + 0)] = mult * (iso + aniso * RX * RX);
			unconstrainedDiffusion[(i * 3 + 0) * nat3 + (j * 3 + 1)] = unconstrainedDiffusion[(j * 3 + 1) * nat3 + (i * 3 + 0)] = mult * aniso * RX * RY;
			unconstrainedDiffusion[(i * 3 + 0) * nat3 + (j * 3 + 2)] = unconstrainedDiffusion[(j * 3 + 2) * nat3 + (i * 3 + 0)] = mult * aniso * RX * RZ;

			unconstrainedDiffusion[(i * 3 + 1) * nat3 + (j * 3 + 0)] = unconstrainedDiffusion[(j * 3 + 0) * nat3 + (i * 3 + 1)] = mult * aniso * RY * RX;
			unconstrainedDiffusion[(i * 3 + 1) * nat3 + (j * 3 + 1)] = unconstrainedDiffusion[(j * 3 + 1) * nat3 + (i * 3 + 1)] = mult * (iso + aniso * RY * RY);
			unconstrainedDiffusion[(i * 3 + 1) * nat3 + (j * 3 + 2)] = unconstrainedDiffusion[(j * 3 + 2) * nat3 + (i * 3 + 1)] = mult * aniso * RY * RZ;

			unconstrainedDiffusion[(i * 3 + 2) * nat3 + (j * 3 + 0)] = unconstrainedDiffusion[(j * 3 + 0) * nat3 + (i * 3 + 2)] = mult * aniso * RZ * RX;
			unconstrainedDiffusion[(i * 3 + 2) * nat3 + (j * 3 + 1)] = unconstrainedDiffusion[(j * 3 + 1) * nat3 + (i * 3 + 2)] = mult * aniso * RZ * RY;
			unconstrainedDiffusion[(i * 3 + 2) * nat3 + (j * 3 + 2)] = unconstrainedDiffusion[(j * 3 + 2) * nat3 + (i * 3 + 2)] = mult * (iso + aniso * RZ * RZ);
	  	}
	}
//#endif // END CPU CODE
	return;
}

void diften::buildFrictionInMF(int buildDiff)
{

	x0 = C * Reff * eta * M_PI;

	buildBmatrixInMF();

	int n  = 3 * mol->getNAtoms();
	int na = 3 * mol->getNActiveAtoms();

	if (!useHydrodynamicsInteractions)
	{
		// Friction tensor is calculated without hydrodynamic interactions
		// F = x0 * Btr B

		for (int i = 0; i < na; ++i)
		{
			for (int j = 0; j < na; ++j)
			{
				friction[i * na + j] = 0.0;
				for (int k = 0; k < n; ++k)
					friction[i * na + j] += bMatrix[k * na + i] * bMatrix[k * na + j];
				friction[i * na + j] *= x0;
			}
		}
	}
	else
	{
		buildUnconstrainedDiffusion();
		int INFO;
#ifdef __GPU__
		TESTING_CHECK(magma_init());
		magma_device_t device;
		magma_getdevice(&device);
		magma_queue_t queue;
		magma_queue_create(device, &queue);

		/* PERFORMS INVERSION TO OBTAIN UNCONSTRAINED FRICTION */
		double *dUncFric;
		TESTING_CHECK(magma_dmalloc(&dUncFric, n * n));
		cudaMemcpy(dUncFric, unconstrainedDiffusion, n * n * sizeof(double), cudaMemcpyHostToDevice);

		TESTING_CHECK(magma_dpotrf_gpu(MagmaLower, n, dUncFric, n, &INFO));
		TESTING_CHECK(magma_dpotri_gpu(MagmaLower, n, dUncFric, n, &INFO));

		dim3 threadsPerBlock, numBlocks;
		threadsPerBlock.x = threadsPerBlock.y = TPB_2D;
		numBlocks.x = 1 + n / threadsPerBlock.x;
		numBlocks.y = 1 + n / threadsPerBlock.y;
		gpu_fill_matrix<<<numBlocks, threadsPerBlock>>>(n, dUncFric);

		/* MULTIPLY B^tr csi B TO OBTAIN THE CONSTRAINED FRICTION */
		double *dB;
		TESTING_CHECK(magma_dmalloc(&dB, n * na));
		magma_dsetmatrix(na, n, bMatrix, na, dB, na, queue);

		double *dUFB;
		TESTING_CHECK(magma_dmalloc(&dUFB, n * na));
		magma_dgemm(MagmaNoTrans, MagmaNoTrans, na, n, n, x0, dB, na, dUncFric, n, 0.0, dUFB, na, queue);
		TESTING_CHECK(magma_free(dUncFric));

		double *dF;
		TESTING_CHECK(magma_dmalloc(&dF, na * na));
		magma_dgemm(MagmaNoTrans, MagmaTrans, na, na, n, 1.0, dUFB, na, dB, na, 0.0, dF, na, queue);
		TESTING_CHECK(magma_free(dUFB));
		TESTING_CHECK(magma_free(dB));

		cudaMemcpy(friction, dF, na * na * sizeof(double), cudaMemcpyDeviceToHost);

		TESTING_CHECK(magma_free(dF));

		TESTING_CHECK(magma_finalize());
#else // END GPU CODE
		char UPLO='L'; 
		dpotrf_(&UPLO, &n, unconstrainedDiffusion, &n, &INFO);
		dpotri_(&UPLO, &n, unconstrainedDiffusion, &n, &INFO);

		for (int r = 0; r < n; ++r)
		{
			for (int c = r + 1; c < n; ++c)
				unconstrainedDiffusion[c * n + r] = unconstrainedDiffusion[r * n + c]; // TODO: check if correct or has to be inverted
		}
		
		char noTr= 'N';
		char Tr= 'T';
		double alpha = 1.0;
		double beta = 0.0;

		double *dUFB = new double[n * na];

		dgemm_(&noTr, &noTr, &na, &n, &n, &x0, bMatrix, &na, unconstrainedDiffusion, &n, &beta, UFB, &na);
		dgemm_(&noTr, &Tr, &na, &na, &n, &alpha, UFB, &na, bMatrix, &na, &beta, friction, &na);

		delete[] UFB;
#endif // END CPU CODE

		delete[] unconstrainedDiffusion;

	}

	// Build diffusion tensor (if needed)

	if (buildDiff)
		buildDiffusionInMF();

	frictionBuilt = 1;
	return;
}

/**********************************************************************/
/* Method to take square root of Roto-Conformational part of friction */
/**********************************************************************/

void diften::makeSqrtOfRotConfFriction(void)
{
	int n = mol->getNActiveAtoms();
	int m = 3 * n;
	int dim = m - 3;
	rotConfFriction = new double[dim * dim];
	
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
			rotConfFriction[i * dim + j] = friction[(i+3) * m + (j+3)];
	}

	// Calculate the square root

	double *B, *W, *Z, *WORK;
	int *ISUPPZ, *IWORK;
	int  i, j;
	int  M;

	/* allocate space for the output parameters and workspace arrays */
	W      = (double *)malloc(dim * sizeof(double));
	Z      = (double *)malloc(dim * dim * sizeof(double));
	ISUPPZ = (int *)malloc(2 * dim *sizeof(int));
	WORK   = (double *)malloc(26 * dim * sizeof(double));
	IWORK  = (int *)malloc(10 * dim * sizeof(int));

	 /* get the eigenvalues and eigenvectors */
	dsyevr('V', 'A', 'L', dim, rotConfFriction, dim, 0, 0, 0, 0, dlamch('S'), &M, W, Z, dim, ISUPPZ, WORK, 26*dim, IWORK, 10*dim);

	/* allocate and initialise a new matrix B=Z*D */
	double lambda;
	B = (double * )malloc(dim * dim * sizeof(double));
	for (j = 0; j < dim; ++j)
	{
		lambda = sqrt(W[j]);
		for (i = 0; i < dim; ++i)
			set_entry(B, dim, i, j, get_entry(Z, dim, i, j) * lambda);
	}

	/* calculate the square root A = B * Z^T */
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, dim, dim, dim, 1, B, dim, Z, dim, 0, rotConfFriction, dim);

	return;
}

/*********************************/
/* Method to build the Diffusion */
/*********************************/

void diften::buildDiffusionInMF(void)
{
	// TODO: inversion of friction tensor and multiplication by kBT
	// NOTE: kB must be espressed in A, for what concerns the length
}
	
/******************/
/* Return methods */
/******************/

vectorOfDoubles diften::getFrictionInMF(void)
{
	int n = mol->getNActiveAtoms();
	n *= n;
	vectorOfDoubles f(n);
	for (int i = 0; i < n; ++i) f[i] = friction[i];
		
	return f;
}


void diften::getFrictionInMF(double *fr)
{
	int n = 3 * mol->getNActiveAtoms();
	for (int i = 0; i < n * n; ++i)
		fr[i] = friction[i];
	return;
}

vectorOfDoubles diften::getRotoConformationalBmatrixInMF(void)
{
	int nr = 3 * mol->getNActiveAtoms();
	int nc = nr - 3;
	vectorOfDoubles m(nr * nc);

	for (int i = 0; i = nr; ++i)
	{
		for (int j = 0; j <  nc; ++j)
			m[i * nc + j] = bMatrix[i * nr + (j + 3)];
	}
	return m;
}

vectorOfDoubles diften::getDiffusionInMF(void)
{
	// TODO
	vectorOfDoubles d;
	return d;
}

vectorOfDoubles diften::getRotoConformationalFrictionInMF(void)
{
	int nr = 3 * mol->getNActiveAtoms();
	int nc = nr - 3;
	vectorOfDoubles rcf(nc * nc);
	for (int i = 0; i < nc; ++i)
	{
		for (int j = 0; j < nc; ++j)
			rcf[i * nc + j] = friction[(i + 3) * nr + (j + 3)];
	}
	return rcf;
}

void diften::getRotoConformationalFrictionInMF(double *rcf)
{
	int nr = 3 * mol->getNActiveAtoms();
	int nc = nr - 3;
	for (int i = 0; i < nc; ++i)
	{
		for (int j = 0; j < nc; ++j)
			rcf[i * nc + j] = friction[(i + 3) * nr + (j + 3)];
	}
	return;
}

vectorOfDoubles diften::getSqrtRotoConformationalFrictionInMF(void)
{
	int nr = 3 * mol->getNActiveAtoms();
	int nc = nr - 3;
	vectorOfDoubles rcf(nc * nc);
	for (int i = 0; i < nc * nc; ++i) rcf[i] = rotConfFriction[i];
	return rcf;
}

void diften::getRotationalFrictionInMF(double **FRR)
{
	int n = 3 * mol->getNActiveAtoms();

	FRR[0][0] = friction[3 * n + 3];
	FRR[0][1] = friction[3 * n + 4];
	FRR[0][2] = friction[3 * n + 5];

	FRR[1][0] = friction[4 * n + 3];
	FRR[1][1] = friction[4 * n + 4];
	FRR[1][2] = friction[4 * n + 5];

	FRR[2][0] = friction[5 * n + 3];
	FRR[2][1] = friction[5 * n + 4];
	FRR[2][2] = friction[5 * n + 5];

	return;
}

void diften::getConformationalFrictionInMF(double **FSS)
{
	int n  = 3 * mol->getNActiveAtoms();
	int nq = n - 6;
	for (int i = 0; i < nq; ++i)
	{
		for (int j = 0; j < nq; ++j)
			FSS[i][j] = friction[(i + 6) * n + (j + 6)];
	}

	return;
}

void diften::getRotoConformationalFrictionInMF(double **FRS)
{
	int n  = 3 * mol->getNActiveAtoms();
	int nq = n - 6;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < nq; ++j)
			FRS[i][j] = friction[(i + 3) * n + (j + 6)];
	}
	
	return;
}

/******************/
/* Output methods */
/******************/

void diften::outputFrictionInMF(std::string s)
{
	int n = 3 * mol->getNActiveAtoms();
	// Standard output
	if (!s.compare("screen"))
	{
		std::cout << std::endl << "Friction tensor in MF:" << std::endl << std::endl;
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
				std::cout << friction[i * n + j] << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	// Output on file
	else
	{
		std::fstream f;
		f.open(s.c_str(), std::ios::out);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
				f << std::scientific << std::setprecision(12) << friction[i * n + j] << "\t";
			f << std::endl;
		}
		f.close();
	}
	return;
}

void diften::outputRotoConformationalBmatrixInMF(std::string s, int multiplyBySqrtX0)
{
	int nr = 3 * mol->getNActiveAtoms();
	int nc = nr - 3;

	double mult = multiplyBySqrtX0 ? sqrt(x0) : 1.0;

	// Standard output
	if (!s.compare("screen"))
	{
		std::cout << std::endl << "Roto-Conformational B matrix in MF:" << std::endl << std::endl;
		for (int i = 0; i < nr; ++i)
		{
			for (int j = 0; j < nc; ++j)
				std::cout << mult * bMatrix[i * nr + ( j + 3)] << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	// Output on file
	else
	{
		std::fstream f;
		f.open(s.c_str(), std::ios::out);
		for (int i = 0; i < nr; ++i)
		{
			for (int j = 0; j < nc; ++j)
				f << mult * bMatrix[i * nr + ( j + 3)] << std::endl;
		}
		f.close();
	}
	return;
}

void diften::outputDiffusionInMF(std::string s)
{
	int n = 3 * mol->getNActiveAtoms();
	// Standard output
	if (!s.compare("screen"))
	{
		std::cout << std::endl << "Diffusion tensor in MF:" << std::endl << std::endl;
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
				std::cout << diffusion[i * n + j] << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	// Output on file
	else
	{
		std::fstream f;
		f.open(s.c_str(), std::ios::out);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
				f << diffusion[i * n + j] << std::endl;
		}
		f.close();
	}
	return;
}

void diften::outputRotoConformationalFrictionInMF(std::string s)
{
	int n = 3 * mol->getNActiveAtoms() - 3;
	int m = n + 3;
	// Standard output
	if (!s.compare("screen"))
	{
		std::cout << std::endl << "Sqrt roto-conformational friction tensor in MF:" << std::endl << std::endl;
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
				std::cout << friction[(i+3) * m + (j+3)] << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	// Output on file
	else
	{
		std::fstream f;
		f.open(s.c_str(), std::ios::out);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
				f << friction[(i+3) * m + (j+3)] << std::endl;
		}
		f.close();
	}
	

}

void diften::outputSqrtRotoConformationalFrictionInMF(std::string s)
{
	int n = 3 * mol->getNActiveAtoms() - 3;
	// Standard output
	if (!s.compare("screen"))
	{
		std::cout << std::endl << "Sqrt roto-conformational friction tensor in MF:" << std::endl << std::endl;
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
				std::cout << rotConfFriction[i * n + j] << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	// Output on file
	else
	{
		std::fstream f;
		f.open(s.c_str(), std::ios::out);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
				f << rotConfFriction[i * n + j] << std::endl;
		}
		f.close();
	}
	

}

