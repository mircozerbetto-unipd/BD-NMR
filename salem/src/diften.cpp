#include "diften.h"

// N.B.: the molecule object has atoms coordinates in AF, i.e. centered in atom 2 and
//       oriented with the Thomson scheme used to build the Z-Matrix.
//
//       The friction tensor (as well as all the other tensors in builgAgI.cpp) will be
//       computed in DF. This frame is obtained by a -90 deg rotation about the AF 
//       y axis.
//
//       Thus, in the buildMatrices() method, first Cartesian coordinates, as
//       well as derivatives, are expressed in DF. Then, the tensors are computed
//       in such a frame since the transformed (by simple rotation) coordinates
//       are used.
//

extern "C"{
void dpotrf_(char *uplo,int *jb,double *A,int *lda,int *info);
void dpotri_(char *uplo,int *jb,double *A,int *lda,int *info);
}

/***************/
/* Constructor */
/***************/

diften::diften()
{
	frictionBuilt = 0;
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

	friction0 = new double[n * n];
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

// Hydrodynamic Interactions (0: off | 1: on [Rotne-Prage model])
void diften::setHydrodynamicInteractions(int hi)
{
	hydroInt = hi;
	return;
}

int diften::getHydrodynamicInteractions(void)
{
	return hydroInt;
}

/*******************************/
/* Method to build the Bmatrix */
/*******************************/

void diften::buildBmatrixInDF(void)
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

	// Get the active atoms

	atom *activeAtoms = new atom[na];
	mol->getActiveAtoms(activeAtoms);

	//
	// Build the B matrix
	//

	// Apply rotation to express the tensor in the frame in DF

	matrix3D Emat = matrix3D(0.0,0.0,1.0,0.0,1.0,0.0,-1.0,0.0,0.0); // AP: 25-06-2020 // -90° // MZ: This is the active rotation (i.e. transforms the coordinates, not the reference system

	double CAx[9], CA[3], EC[3], dCA[3];
	double EI[9];
	EI[0] = Emat.xx; EI[1] = Emat.xy; EI[2] = Emat.xz;
	EI[3] = Emat.yx; EI[4] = Emat.yy; EI[5] = Emat.yz;
	EI[6] = Emat.zx; EI[7] = Emat.zy; EI[8] = Emat.zz;

	int count;

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

		CA[0] = r[0];
		CA[1] = r[1];
		CA[2] = r[2];

		EC[0] = EI[0] * CA[0] + EI[1] * CA[1] + EI[2] * CA[2];
		EC[1] = EI[3] * CA[0] + EI[4] * CA[1] + EI[5] * CA[2];
		EC[2] = EI[6] * CA[0] + EI[7] * CA[1] + EI[8] * CA[2];

		CAx[0] =  0.0;
		CAx[1] = -EC[2];
		CAx[2] =  EC[1];
		CAx[3] =  EC[2];
		CAx[4] =  0.0;
		CAx[5] = -EC[0];
		CAx[6] = -EC[1];
		CAx[7] =  EC[0];
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

		dCA[0] = d[0];
		dCA[1] = d[1];
		dCA[2] = d[2];

		EC[0] = EI[0] * dCA[0] + EI[1] * dCA[1] + EI[2] * dCA[2];
		EC[1] = EI[3] * dCA[0] + EI[4] * dCA[1] + EI[5] * dCA[2];
		EC[2] = EI[6] * dCA[0] + EI[7] * dCA[1] + EI[8] * dCA[2];

		bMatrix[(i + 0) * m + 6] = EC[0];
		bMatrix[(i + 1) * m + 6] = EC[1];
		bMatrix[(i + 2) * m + 6] = EC[2];

		// q = d3
		d = mol->getFirstDerivative(aid, 1, 3);

		dCA[0] = d[0];
		dCA[1] = d[1];
		dCA[2] = d[2];

		EC[0] = EI[0] * dCA[0] + EI[1] * dCA[1] + EI[2] * dCA[2];
		EC[1] = EI[3] * dCA[0] + EI[4] * dCA[1] + EI[5] * dCA[2];
		EC[2] = EI[6] * dCA[0] + EI[7] * dCA[1] + EI[8] * dCA[2];

		bMatrix[(i + 0) * m + 7] = EC[0];
		bMatrix[(i + 1) * m + 7] = EC[1];
		bMatrix[(i + 2) * m + 7] = EC[2];

		// q = theta3
		d = mol->getFirstDerivative(aid, 2, 3);

		dCA[0] = d[0];
		dCA[1] = d[1];
		dCA[2] = d[2];

		EC[0] = EI[0] * dCA[0] + EI[1] * dCA[1] + EI[2] * dCA[2];
		EC[1] = EI[3] * dCA[0] + EI[4] * dCA[1] + EI[5] * dCA[2];
		EC[2] = EI[6] * dCA[0] + EI[7] * dCA[1] + EI[8] * dCA[2];

		bMatrix[(i + 0) * m + 8] = EC[0];
		bMatrix[(i + 1) * m + 8] = EC[1];
		bMatrix[(i + 2) * m + 8] = EC[2];

		count = 9;

		for (int a2id = 4; a2id <= mol->getNAtoms(); ++a2id)
		{
			if (mol->getAtom(a2id).isActive())
			{
				// d
				d = mol->getFirstDerivative(aid, 1, a2id);

				dCA[0] = d[0];
				dCA[1] = d[1];
				dCA[2] = d[2];

				EC[0] = EI[0] * dCA[0] + EI[1] * dCA[1] + EI[2] * dCA[2];
				EC[1] = EI[3] * dCA[0] + EI[4] * dCA[1] + EI[5] * dCA[2];
				EC[2] = EI[6] * dCA[0] + EI[7] * dCA[1] + EI[8] * dCA[2];

				bMatrix[(i + 0) * m + count] = EC[0];
				bMatrix[(i + 1) * m + count] = EC[1];
				bMatrix[(i + 2) * m + count] = EC[2];

				count++;

				// theta
				d = mol->getFirstDerivative(aid, 2, a2id);

				dCA[0] = d[0];
				dCA[1] = d[1];
				dCA[2] = d[2];

				EC[0] = EI[0] * dCA[0] + EI[1] * dCA[1] + EI[2] * dCA[2];
				EC[1] = EI[3] * dCA[0] + EI[4] * dCA[1] + EI[5] * dCA[2];
				EC[2] = EI[6] * dCA[0] + EI[7] * dCA[1] + EI[8] * dCA[2];

				bMatrix[(i + 0) * m + count] = EC[0];
				bMatrix[(i + 1) * m + count] = EC[1];
				bMatrix[(i + 2) * m + count] = EC[2];

				count++;

				// phi
				d = mol->getFirstDerivative(aid, 3, a2id);

				dCA[0] = d[0];
				dCA[1] = d[1];
				dCA[2] = d[2];

				EC[0] = EI[0] * dCA[0] + EI[1] * dCA[1] + EI[2] * dCA[2];
				EC[1] = EI[3] * dCA[0] + EI[4] * dCA[1] + EI[5] * dCA[2];
				EC[2] = EI[6] * dCA[0] + EI[7] * dCA[1] + EI[8] * dCA[2];

				bMatrix[(i + 0) * m + count] = EC[0];
				bMatrix[(i + 1) * m + count] = EC[1];
				bMatrix[(i + 2) * m + count] = EC[2];

				count++;
			}
		}
	}

#ifdef DEBUG	
		std::cout << "* Bmatrix:" << std::endl;
		std::ofstream f_Bmat;
		f_Bmat.open("Bmat.dat", std::ios::out);
		for (int i = 0; i < m; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				std::cout << std::scientific << std::setprecision(5)  << bMatrix[i * m + j] << "\t";
				f_Bmat    << std::scientific << std::setprecision(12) << bMatrix[i * m + j] << "\t";
			}
			std::cout << std::endl;
			f_Bmat    << std::endl;
		}
		std::cout << "*****************" << std::endl;
		f_Bmat.close();
#endif

	return;
}

/********************************/
/* Method to build the Friction */
/********************************/

// Built in DF

void diften::buildFrictionInDF(int buildDiff)
{

	x0 = C * Reff * eta * M_PI;

	buildBmatrixInDF();

	// Friction tensor is calculated without hydrodynamic interactions
	// F = x0 * Btr B

	int n  = 3 * mol->getNAtoms();
	int na = 3 * mol->getNActiveAtoms();

	if (!hydroInt)
	{
		for (int i = 0; i < na; ++i)
		{
			for (int j = 0; j < na; ++j)
			{
				friction[i * na + j] = 0.0;
				for (int k = 0; k < n; ++k) friction[i * na + j] += bMatrix[k * na + i] * bMatrix[k * na + j];
				friction[i * na + j] *= x0;
			}
		}
	}
	else
	{
		// Get the active atoms
		atom *activeAtoms = new atom[na / 3];
		mol->getActiveAtoms(activeAtoms);


		atom at;
		vectorOfDoubles r1(3), r2(3),ri(3), rj(3);
		double RX, RY, RZ, R, mult, iso, aniso;
		char UPLO = 'L';
		int info;
		matrix3D Emat = matrix3D(0.0,0.0,1.0,0.0,1.0,0.0,-1.0,0.0,0.0);
		double EI[9];
		EI[0] = Emat.xx; EI[1] = Emat.xy; EI[2] = Emat.xz;
		EI[3] = Emat.yx; EI[4] = Emat.yy; EI[5] = Emat.yz;
		EI[6] = Emat.zx; EI[7] = Emat.zy; EI[8] = Emat.zz;


		for (int i = 0; i < mol->getNAtoms(); i++)
		{

			/* Diagonal tensor */

			friction0[(i * 3 + 0) * n + (i * 3 + 0)] = 1.0;
			friction0[(i * 3 + 0) * n + (i * 3 + 1)] = 0.0;
			friction0[(i * 3 + 0) * n + (i * 3 + 2)] = 0.0;

			friction0[(i * 3 + 1) * n + (i * 3 + 0)] = 0.0;
			friction0[(i * 3 + 1) * n + (i * 3 + 1)] = 1.0;
			friction0[(i * 3 + 1) * n + (i * 3 + 2)] = 0.0;

			friction0[(i * 3 + 2) * n + (i * 3 + 0)] = 0.0;
			friction0[(i * 3 + 2) * n + (i * 3 + 1)] = 0.0;
			friction0[(i * 3 + 2) * n + (i * 3 + 2)] = 1.0;

			r1 = mol->getAtomPosition(i + 1);
			ri[0] = EI[0] * r1[0] + EI[1] * r1[1] + EI[2] * r1[2];
			ri[1] = EI[3] * r1[0] + EI[4] * r1[1] + EI[5] * r1[2];
			ri[2] = EI[6] * r1[0] + EI[7] * r1[1] + EI[8] * r1[2];

			for (int j = i + 1; j < mol->getNAtoms(); j++)
			{

				/* Out of diagonal tensors */

				r2 = mol->getAtomPosition(j + 1);
				rj[0] = EI[0] * r2[0] + EI[1] * r2[1] + EI[2] * r2[2];
				rj[1] = EI[3] * r2[0] + EI[4] * r2[1] + EI[5] * r2[2];
				rj[2] = EI[6] * r2[0] + EI[7] * r2[1] + EI[8] * r2[2];
	    
				RX = ri[0] - rj[0]; RY = ri[1] - rj[1]; RZ = ri[2] - rj[2];
				R = sqrt(RX * RX + RY * RY + RZ * RZ);

				if (R >= 2.0 * Reff)
				{
					mult = 3.0 * Reff / (4.0 * R);
					iso = 1.0 + 2.0 * Reff * Reff / (3.0 * R * R);
					aniso = (1.0 - 2.0 * Reff * Reff / (R * R)) / (R * R);
				}
				else
				{
					mult = 1.0;
					iso = 1.0 - 9.0 * R / (32.0 * Reff);
					aniso = 3.0 / (32.0 * R * Reff);
				}

				friction0[(i * 3 + 0) * n + (j * 3 + 0)] = friction0[(j * 3 + 0) * n + (i * 3 + 0)] = mult * (iso + aniso * RX * RX);
				friction0[(i * 3 + 0) * n + (j * 3 + 1)] = friction0[(j * 3 + 1) * n + (i * 3 + 0)] = mult * aniso * RX * RY;
				friction0[(i * 3 + 0) * n + (j * 3 + 2)] = friction0[(j * 3 + 2) * n + (i * 3 + 0)] = mult * aniso * RX * RZ;
                                                                                                              
				friction0[(i * 3 + 1) * n + (j * 3 + 0)] = friction0[(j * 3 + 0) * n + (i * 3 + 1)] = mult * aniso * RY * RX;
				friction0[(i * 3 + 1) * n + (j * 3 + 1)] = friction0[(j * 3 + 1) * n + (i * 3 + 1)] = mult * (iso + aniso * RY * RY);
				friction0[(i * 3 + 1) * n + (j * 3 + 2)] = friction0[(j * 3 + 2) * n + (i * 3 + 1)] = mult * aniso * RY * RZ;
                                                                                                              
				friction0[(i * 3 + 2) * n + (j * 3 + 0)] = friction0[(j * 3 + 0) * n + (i * 3 + 2)] = mult * aniso * RZ * RX;
				friction0[(i * 3 + 2) * n + (j * 3 + 1)] = friction0[(j * 3 + 1) * n + (i * 3 + 2)] = mult * aniso * RZ * RY;
				friction0[(i * 3 + 2) * n + (j * 3 + 2)] = friction0[(j * 3 + 2) * n + (i * 3 + 2)] = mult * (iso + aniso * RZ * RZ);
			}
		}

#ifdef DEBUG
		std::ofstream f_dhi;
		f_dhi.open("DEBUG_hydroint.dat", std::ios::out);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				f_dhi << std::scientific << std::setprecision(12) << friction0[i * n + j] << "\t";
			}
			f_dhi << std::endl;
		}
		std::cout << "******************************************************" << std::endl;
		f_dhi.close();
#endif

		// Invert the 3Natoms x 3Natoms translational diffusion to obtain the friction

		dpotrf_(&UPLO,&n,friction0,&n,&info);
		dpotri_(&UPLO,&n,friction0,&n,&info);

		if(UPLO=='L')
		{
			for (int i = 0; i < n; ++i)
			{
				for (int j = i + 1; j < n; ++j)
				{
					friction0[j * n + i] = friction0[i * n + j];
				}
			}
		
		}
		else
		{
		}

#ifdef DEBUG
		std::ofstream f_dhi2;
		f_dhi2.open("DEBUG_hydroint_2.dat", std::ios::out);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				f_dhi2 << std::scientific << std::setprecision(12) << friction0[i * n + j] << "\t";
			}
			f_dhi2 << std::endl;
		}
		std::cout << "******************************************************" << std::endl;
		f_dhi2.close();
#endif

		// T-R-C friction = B^tr friction0 B

		double *dummy = new double[n * na];
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < na; ++j)
			{
				dummy[i * na + j] = 0.0;
				for (int k = 0; k < n; ++k)
					dummy[i * na + j] += friction0[i * n + k] * bMatrix[k * na + j];
			}
		}

		for (int i = 0; i < na; ++i)
		{
			for (int j = 0; j < na; ++j)
			{
				friction[i * na + j] = 0.0;
				for (int k = 0; k < n; ++k)
					friction[i * na + j] += bMatrix[k * na + i] * dummy[k * na + j];
				friction[i * na + j] *= x0;
			}
		}
	}

	// Build diffusion tensor (if needed)

	if (buildDiff)
		buildDiffusionInDF();

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
	int  M, I1, I2, INFO;

	char Vch, Ach, Lch, Sch;

	/* allocate space for the output parameters and workspace arrays */
	W      = (double *)malloc(dim * sizeof(double));
	Z      = (double *)malloc(dim * dim * sizeof(double));
	ISUPPZ = (int *)malloc(2 * dim *sizeof(int));
	WORK   = (double *)malloc(26 * dim * sizeof(double));
	IWORK  = (int *)malloc(10 * dim * sizeof(int));

	Vch = 'V';
	Ach = 'A';
	Lch = 'L';
	I1 = 26 * dim;
	I2 = 10 * dim;

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

void diften::buildDiffusionInDF(void)
{
	// TODO: inversion of friction tensor and multiplication by kBT
	// NOTE: kB must be espressed in A, for what concerns the length
}
	
/******************/
/* Return methods */
/******************/

vectorOfDoubles diften::getFrictionInDF(void)
{
	int n = 3 * mol->getNActiveAtoms();
	n *= n;
	vectorOfDoubles f(n);
	for (int i = 0; i < n; ++i) f[i] = friction[i];
		
	return f;
}


void diften::getFrictionInDF(double *fr)
{
	fr = friction;
	return;
}

vectorOfDoubles diften::getRotoConformationalBmatrixInDF(void)
{
	int nr = 3 * mol->getNActiveAtoms();
	int nc = nr - 3;
	vectorOfDoubles m(nc * nc);

	for (int i = 0; i = nc; ++i)
	{
		for (int j = 0; j <  nc; ++j)
			m[i * nc + j] = bMatrix[(i + 3) * nc + (j + 3)];
	}
	return m;
}

void diften::getRotoConformationalBmatrixInDF(double *m)
{
	int nr = 3 * mol->getNActiveAtoms();
	int nc = nr - 3;
	for (int i = 0; i = nc; ++i)
	{
		for (int j = 0; j <  nc; ++j)
			m[i * nc + j] = bMatrix[(i + 3) * nc + (j + 3)];
	}
	return;
}

void diften::getFrictionTensor(double *extf)
{
	int nr = 3 * mol->getNActiveAtoms();
	for (int i = 0; i < nr; ++i)
	{
		for (int j = 0; j <  nr; ++j)
			extf[i * nr + j] = friction[i * nr + j];
	}
	return;
}

vectorOfDoubles diften::getDiffusionInDF(void)
{
	// TODO
	vectorOfDoubles d;
	return d;
}

vectorOfDoubles diften::getRotoConformationalFrictionInDF(void)
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

void diften::getRotoConformationalFrictionInDF(double *rcf)
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

vectorOfDoubles diften::getSqrtRotoConformationalFrictionInDF(void)
{
	int nr = 3 * mol->getNActiveAtoms();
	int nc = nr - 3;
	vectorOfDoubles rcf(nc * nc);
	for (int i = 0; i < nc * nc; ++i) rcf[i] = rotConfFriction[i];
	return rcf;
}

void diften::getRotationalFrictionInDF(double **FRR)
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

void diften::getConformationalFrictionInDF(double **FSS)
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

void diften::getRotoConformationalFrictionInDF(double **FRS)
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

void diften::outputFrictionInDF(std::string s)
{
	int n = 3 * mol->getNActiveAtoms();
	// Standard output
	if (!s.compare("screen"))
	{
		std::cout << std::endl << "Friction tensor in DF:" << std::endl << std::endl;
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

void diften::outputRotoConformationalBmatrixInDF(std::string s, int multiplyBySqrtX0)
{
	int nr = 3 * mol->getNActiveAtoms();
	int nc = nr - 3;

	double mult = multiplyBySqrtX0 ? sqrt(x0) : 1.0;

	// Standard output
	if (!s.compare("screen"))
	{
		std::cout << std::endl << "Roto-Conformational B matrix in DF:" << std::endl << std::endl;
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

void diften::outputDiffusionInDF(std::string s)
{
	int n = 3 * mol->getNActiveAtoms();
	// Standard output
	if (!s.compare("screen"))
	{
		std::cout << std::endl << "Diffusion tensor in DF:" << std::endl << std::endl;
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

void diften::outputRotoConformationalFrictionInDF(std::string s)
{
	int n = 3 * mol->getNActiveAtoms() - 3;
	int m = n + 3;
	// Standard output
	if (!s.compare("screen"))
	{
		std::cout << std::endl << "Sqrt roto-conformational friction tensor in DF:" << std::endl << std::endl;
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

void diften::outputSqrtRotoConformationalFrictionInDF(std::string s)
{
	int n = 3 * mol->getNActiveAtoms() - 3;
	// Standard output
	if (!s.compare("screen"))
	{
		std::cout << std::endl << "Sqrt roto-conformational friction tensor in DF:" << std::endl << std::endl;
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

