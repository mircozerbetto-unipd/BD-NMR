#include "buildAGI.h"

// N.B.: the molecule object has atoms coordinates in AF, i.e. centered in atom 2 and
//       oriented with the Thomson scheme used to build the Z-Matrix.
//
//       All of the tensors (as well as the friction tensor in diften.cpp) will be
//       computed in DF. This frame is obtained by a -90 deg rotation about the AF 
//       y axis.
//
//       Thus, in the buildMatrices() method, first Cartesian coordinates, as
//       well as derivatives, are expressed in DF. Then, the tensors are computed
//       in such a frame since the transformed (by simple rotation) coordinates
//       are used.
//


/***************/
/* Constructor */
/***************/

buildAGI::buildAGI()
{
	return;
}

buildAGI::buildAGI(molecule *m)
{
	mol = m;

	int nq = 3 * mol->getNActiveAtoms() - 6;

	Amat = vectorOfDoubles(nq * 3);
	Gmat = vectorOfDoubles(nq * nq);
	Iten = vectorOfDoubles(3 * 3);
	InvIten = vectorOfDoubles(3 * 3);
	ca   = vectorOfDoubles(3 * mol->getNAtoms());
	dca  = vectorOfDoubles(3 * mol->getNAtoms() * nq);
	return;
}

/*****************/
/* Deconstructor */
/*****************/

buildAGI::~buildAGI()
{
	return;
}

void buildAGI::setMol(molecule *m)
{
	mol = m;

	int nq = 3 * mol->getNActiveAtoms() - 6;

	Amat = vectorOfDoubles(nq * 3);
	Gmat = vectorOfDoubles(nq * nq);
	Iten = vectorOfDoubles(3 * 3);
	InvIten = vectorOfDoubles(3 * 3);
	ca   = vectorOfDoubles(3 * mol->getNAtoms());
	dca  = vectorOfDoubles(3 * mol->getNAtoms() * nq);

	return;
}

/********************************/
/* Method to buld all the stuff */
/********************************/

// N.B.: all of the tensors will be expressed in DF

void buildAGI::buildMatrices(void)
{
	/**********************************************/
	/* CHECK IS EVERYTHNG IS OK WITH THE MOLECULE */
	/**********************************************/

	int check = mol->areRefAtomsSet();
	if (!check)
	{
		std::cout << "\n\nERROR: in " << __FILE__ << ", line " << __LINE__ <<": reference atoms on the molecule have not been set. Please, use the setMainDihedralAngle method of molecule.cpp class to set the atoms before calling this routine.\n\n";
		exit(1);
	}

	check = mol->hasZMatrix();
	if (!check)
	{
		std::cout << "\n\nERROR: in " << __FILE__ << ", line " << __LINE__ <<": no Z-Matrix built for the molecule. Please, call the buildZMatrix() method of the molecule.cpp class before calling this routine.\n\n";
		exit(1);
	}

	/*******************************/
	/* COORDINATES AND DERIVATIVES */
	/*******************************/

	int nActiveAtoms = mol->getNActiveAtoms();
	atom *activeAtoms = new atom[nActiveAtoms];
	mol->getActiveAtoms(activeAtoms);
	
	int nq = 3 * nActiveAtoms - 6;

	matrix3D Emat = matrix3D(0.0,0.0,1.0,0.0,1.0,0.0,-1.0,0.0,0.0); // AP: 25-06-2020 // -90° // MZ: This is the active rotation (i.e. transforms the coordinates, not the reference system

	// 1. coordinates

	vectorOfDoubles r(3), m(mol->getNAtoms());
#ifdef DEBUG
	std::ofstream f;
	f.open("CoordXYZ.dat", std::ios::out);
#endif
	for (int i = 0; i < mol->getNAtoms(); ++i)
	{
		r = mol->getAtomPosition(i+1);
		ca[i * 3 + 0] = Emat.xx * r[0] + Emat.xy * r[1] + Emat.xz * r[2];
		ca[i * 3 + 1] = Emat.yx * r[0] + Emat.yy * r[1] + Emat.yz * r[2];
		ca[i * 3 + 2] = Emat.zx * r[0] + Emat.zy * r[1] + Emat.zz * r[2];
#ifdef DEBUG
#ifdef OUT_COL_MAT
		f << ca[i * 3 + 0] << std::endl << ca[i * 3 + 1] << std::endl << ca[i * 3 + 2] << std::endl;
#else
		f << ca[i * 3 + 0] << "\t" << ca[i * 3 + 1] << "\t" << ca[i * 3 + 2] << std::endl;
#endif
#endif
	}
#ifdef DEBUG
	f.close();
#endif

	// 2. first derivatives
	vectorOfDoubles d(3);

#ifdef DEBUG
	f.open("Deriv_ca.dat", std::ios::out);
#endif

	int q, count;
	for (int i = 1; i <= mol->getNAtoms(); ++i)
	{
		// q1 = d2
		d = mol->getFirstDerivative(i, 1, 2);
		dca[(3 * (i-1) + 0) * nq + 0] = Emat.xx * d[0] + Emat.xy * d[1] + Emat.xz * d[2];
		dca[(3 * (i-1) + 1) * nq + 0] = Emat.yx * d[0] + Emat.yy * d[1] + Emat.yz * d[2];
		dca[(3 * (i-1) + 2) * nq + 0] = Emat.zx * d[0] + Emat.zy * d[1] + Emat.zz * d[2];
#ifdef DEBUG
#ifdef OUT_COL_MAT
		f << dca[(3 * (i-1) + 0) * nq + 0] << std::endl << dca[(3 * (i-1) + 1) * nq + 0] << std::endl << dca[(3 * (i-1) + 2) * nq + 0] << std::endl;
#else
		f << dca[(3 * (i-1) + 0) * nq + 0] << "\t" << dca[(3 * (i-1) + 1) * nq + 0] << "\t" << dca[(3 * (i-1) + 2) * nq + 0] << "\t";
#endif
#endif

		// q2 = d3
		d = mol->getFirstDerivative(i, 1, 3);
		dca[(3 * (i - 1) + 0) * nq + 1] = Emat.xx * d[0] + Emat.xy * d[1] + Emat.xz * d[2];
		dca[(3 * (i - 1) + 1) * nq + 1] = Emat.yx * d[0] + Emat.yy * d[1] + Emat.yz * d[2];
		dca[(3 * (i - 1) + 2) * nq + 1] = Emat.zx * d[0] + Emat.zy * d[1] + Emat.zz * d[2];
#ifdef DEBUG
#ifdef OUT_COL_MAT
		f << dca[(3 * (i - 1) + 0) * nq + 1] << std::endl << dca[(3 * (i - 1) + 1) * nq + 1] << std::endl << dca[(3 * (i - 1) + 2) * nq + 1] << std::endl;
#else
		f << dca[(3 * (i - 1) + 0) * nq + 1] << "\t" << dca[(3 * (i - 1) + 1) * nq + 1] << "\t" << dca[(3 * (i - 1) + 2) * nq + 1] << "\t";
#endif
#endif

		// q3 = theta3
		d = mol->getFirstDerivative(i, 2, 3);
		dca[(3 * (i - 1) + 0) * nq + 2] = Emat.xx * d[0] + Emat.xy * d[1] + Emat.xz * d[2];
		dca[(3 * (i - 1) + 1) * nq + 2] = Emat.yx * d[0] + Emat.yy * d[1] + Emat.yz * d[2];
		dca[(3 * (i - 1) + 2) * nq + 2] = Emat.zx * d[0] + Emat.zy * d[1] + Emat.zz * d[2];
#ifdef DEBUG
#ifdef OUT_COL_MAT
		f << dca[(3 * (i - 1) + 0) * nq + 2] << std::endl << dca[(3 * (i - 1) + 1) * nq + 2] << std::endl << dca[(3 * (i - 1) + 2) * nq + 2] << std::endl;
#else
		f << dca[(3 * (i - 1) + 0) * nq + 2] << "\t" << dca[(3 * (i - 1) + 1) * nq + 2] << "\t" << dca[(3 * (i - 1) + 2) * nq + 2] << "\t";
#endif
#endif

		// all other internal coordinates

		count = 3;
		for (int q = 4; q <= mol->getNAtoms(); ++q)
		{

			if (mol->getAtom(q).isActive())
			{
				// d_q
				d = mol->getFirstDerivative(i, 1, q);
				dca[(3 * (i - 1) + 0) * nq + count] = Emat.xx * d[0] + Emat.xy * d[1] + Emat.xz * d[2];
				dca[(3 * (i - 1) + 1) * nq + count] = Emat.yx * d[0] + Emat.yy * d[1] + Emat.yz * d[2];
				dca[(3 * (i - 1) + 2) * nq + count] = Emat.zx * d[0] + Emat.zy * d[1] + Emat.zz * d[2];
	#ifdef DEBUG
	#ifdef OUT_COL_MAT
				f << dca[(3 * (i - 1) + 0) * nq + count] << std::endl << dca[(3 * (i - 1) + 1) * nq + count] << std::endl << dca[(3 * (i - 1) + 2) * nq + count] << std::endl;
	#else
				f << dca[(3 * (i - 1) + 0) * nq + count] << "\t" << dca[(3 * (i - 1) + 1) * nq + count] << "\t" << dca[(3 * (i - 1) + 2) * nq + count] << "\t";
	#endif
	#endif
				count++;

				// theta_q
				d = mol->getFirstDerivative(i, 2, q);
				dca[(3 * (i - 1) + 0) * nq + count] = Emat.xx * d[0] + Emat.xy * d[1] + Emat.xz * d[2];
				dca[(3 * (i - 1) + 1) * nq + count] = Emat.yx * d[0] + Emat.yy * d[1] + Emat.yz * d[2];
				dca[(3 * (i - 1) + 2) * nq + count] = Emat.zx * d[0] + Emat.zy * d[1] + Emat.zz * d[2];
	#ifdef DEBUG
	#ifdef OUT_COL_MAT
				f << dca[(3 * (i - 1) + 0) * nq + count] << std::endl << dca[(3 * (i - 1) + 1) * nq + count] << std::endl << dca[(3 * (i - 1) + 2) * nq + count] << std::endl;
	#else
				f << dca[(3 * (i - 1) + 0) * nq + count] << "\t" << dca[(3 * (i - 1) + 1) * nq + count] << "\t" << dca[(3 * (i - 1) + 2) * nq + count] << "\t";
	#endif
	#endif
				count++;

				// phi_q
				d = mol->getFirstDerivative(i, 3, q);
				dca[(3 * (i - 1) + 0) * nq + count] = Emat.xx * d[0] + Emat.xy * d[1] + Emat.xz * d[2];
				dca[(3 * (i - 1) + 1) * nq + count] = Emat.yx * d[0] + Emat.yy * d[1] + Emat.yz * d[2];
				dca[(3 * (i - 1) + 2) * nq + count] = Emat.zx * d[0] + Emat.zy * d[1] + Emat.zz * d[2];
	#ifdef DEBUG
	#ifdef OUT_COL_MAT
				f << dca[(3 * (i - 1) + 0) * nq + count] << std::endl << dca[(3 * (i - 1) + 1) * nq + count] << std::endl << dca[(3 * (i - 1) + 2) * nq + count] << std::endl;
	#else
	#endif
				f << dca[(3 * (i - 1) + 0) * nq + count] << "\t" << dca[(3 * (i - 1) + 1) * nq + count] << "\t" << dca[(3 * (i - 1) + 2) * nq + count] << "\t";
	#endif
				count++;
			}
		}

#ifdef DEBUG
#ifndef OUT_COL_MAT
		f << std::endl;
#endif
#endif
	}
#ifdef DEBUG
	f.close();
#endif

	// 3. masses

#ifdef DEBUG
	f.open("AtomMass.dat", std::ios::out);
#endif
	atom at;
	for (int i = 0; i < m.size(); ++i)
	{
		at = mol->getAtom(i+1);
		m[i] = at.getMass();
#ifdef DEBUG
		f << m[i] << std::endl;
#endif
	}
#ifdef DEBUG
	f.close();
#endif

	/******************/
	/* INERTIA TENSOR */
	/******************/

	// Computed in DF //

	double Ixx = 0.0, Iyy = 0.0, Izz = 0.0, Ixy = 0.0, Ixz = 0.0, Iyz = 0.0;
	for (int i = 0; i < mol->getNAtoms(); ++i)
	{
		Ixx += m[i]  * (ca[i * 3 + 1] * ca[i * 3 + 1] + ca[i * 3 + 2] * ca[i * 3 + 2]);
		Iyy += m[i]  * (ca[i * 3 + 0] * ca[i * 3 + 0] + ca[i * 3 + 2] * ca[i * 3 + 2]);
		Izz += m[i]  * (ca[i * 3 + 0] * ca[i * 3 + 0] + ca[i * 3 + 1] * ca[i * 3 + 1]);
		Ixy -= m[i] * ca[i * 3 + 0] * ca[i * 3 + 1];
		Ixz -= m[i] * ca[i * 3 + 0] * ca[i * 3 + 2];
		Iyz -= m[i] * ca[i * 3 + 1] * ca[i * 3 + 2];
	}
#ifdef DEBUG
	std::cout << "Debug info: Check if I is diagonal\n";
	std::cout << "Ixx = " << Ixx << std::endl << "Iyy = " << Iyy << std::endl << "Izz = " << Izz << std::endl << "Ixy = " << Ixy << std::endl << "Ixz = " << Ixz << std::endl << "Iyz = " << Iyz << std::endl;
	std::fstream f_inertia;
	f_inertia.open("inertia.dat", std::ios::out);
	f_inertia << std::scientific << std::setprecision(16) <<Ixx << "\t" << Ixy << "\t" << Ixz << std::endl;
	f_inertia << std::scientific << std::setprecision(16) <<Ixy << "\t" << Iyy << "\t" << Iyz << std::endl;
	f_inertia << std::scientific << std::setprecision(16) <<Ixz << "\t" << Iyz << "\t" << Izz << std::endl;
	f_inertia.close();
	std::cout << std::endl;
#endif

	Iten[0] = Ixx;      Iten[1] = Ixy;      Iten[2] = Ixz;
	Iten[3] = Ixy;      Iten[4] = Iyy;      Iten[5] = Iyz;
	Iten[6] = Ixz;      Iten[7] = Iyz;      Iten[8] = Izz;

	// Compute the inverse of I //

	double DetIten = 0.0;

	DetIten = Iten[0] * ((Iten[4]*Iten[8])-(Iten[5]*Iten[7]));
	DetIten = DetIten - (Iten[1] * ((Iten[3]*Iten[8])-(Iten[5]*Iten[6])));
	DetIten = DetIten + (Iten[2] * ((Iten[3]*Iten[7])-(Iten[4]*Iten[6])));

	InvIten[0] = Iten[0];    InvIten[1] = Iten[3];    InvIten[2] = Iten[6];
	InvIten[3] = Iten[1];    InvIten[4] = Iten[4];    InvIten[5] = Iten[7];
	InvIten[6] = Iten[2];    InvIten[7] = Iten[5];    InvIten[8] = Iten[8];

	double AdjIten[9];
	AdjIten[0] = (InvIten[4]*InvIten[8]) - (InvIten[5]*InvIten[7]);
	AdjIten[1] = -((InvIten[3]*InvIten[8]) - (InvIten[5]*InvIten[6]));
	AdjIten[2] = (InvIten[3]*InvIten[7]) - (InvIten[4]*InvIten[6]);
	AdjIten[3] = -((InvIten[1]*InvIten[8]) - (InvIten[2]*InvIten[7]));
	AdjIten[4] = (InvIten[0]*InvIten[8]) - (InvIten[2]*InvIten[6]);
	AdjIten[5] = -((InvIten[0]*InvIten[7]) - (InvIten[1]*InvIten[6]));
	AdjIten[6] = (InvIten[1]*InvIten[5]) - (InvIten[2]*InvIten[4]);
	AdjIten[7] = -((InvIten[0]*InvIten[5]) - (InvIten[2]*InvIten[3]));
	AdjIten[8] = (InvIten[0]*InvIten[4]) - (InvIten[1]*InvIten[3]);

	for(int i=0; i<9; i++) InvIten[i] = AdjIten[i]/DetIten;

	/************/
	/* A MATRIX */
	/************/

	// Computed in DF //

#ifdef DEBUG
	f.open("A_Matrix.dat", std::ios::out);
#endif

	double dummy[3];
	for (int i = 0; i < nq; ++i)
	{
		Amat[i * 3 + 0] = 0.0;
		Amat[i * 3 + 1] = 0.0;
		Amat[i * 3 + 2] = 0.0;
		dummy[0] = 0.0;
		dummy[1] = 0.0;
		dummy[2] = 0.0;
		for (int j = 0; j < mol->getNAtoms(); ++j)
		{
			Amat[i * 3 + 0] += m[j] * (ca[j * 3 + 1] * dca[(j * 3 + 2) * nq + i] - ca[j * 3 + 2] * dca[(j * 3 + 1) * nq + i]);
			Amat[i * 3 + 1] += m[j] * (ca[j * 3 + 2] * dca[(j * 3 + 0) * nq + i] - ca[j * 3 + 0] * dca[(j * 3 + 2) * nq + i]);
			Amat[i * 3 + 2] += m[j] * (ca[j * 3 + 0] * dca[(j * 3 + 1) * nq + i] - ca[j * 3 + 1] * dca[(j * 3 + 0) * nq + i]);
		}
		dummy[0] = (InvIten[0] * Amat[i*3 + 0]) + (InvIten[1] * Amat[i*3 + 1]) + (InvIten[2] * Amat[i*3 + 2]);
		dummy[1] = (InvIten[3] * Amat[i*3 + 0]) + (InvIten[4] * Amat[i*3 + 1]) + (InvIten[5] * Amat[i*3 + 2]);
		dummy[2] = (InvIten[6] * Amat[i*3 + 0]) + (InvIten[7] * Amat[i*3 + 1]) + (InvIten[8] * Amat[i*3 + 2]);
		Amat[i * 3 + 0] = dummy[0];
                Amat[i * 3 + 1] = dummy[1];
                Amat[i * 3 + 2] = dummy[2];
#ifdef DEBUG
#ifdef OUT_COL_MAT
		f << Amat[i * 3 + 0] << std::endl << Amat[i * 3 + 1] << std::endl << Amat[i * 3 + 2] << std::endl;
#else
		f << std::scientific << std::setprecision(16) << Amat[i * 3 + 0] << "\t" << Amat[i * 3 + 1] << "\t" << Amat[i * 3 + 2] << std::endl;
#endif
#endif
	}
#ifdef DEBUG
	f.close();
#endif

	/***********************/
	/* A1 MATRIX = Atr I A */
	/***********************/
	
	vectorOfDoubles A1(nq * nq);

	for (int i = 0; i < nq; ++i)
	{
		for (int j = i; j < nq; ++j)
		{
			A1[i * nq + j] = 0.0;
			for (int k = 0; k < 3; ++k)
			{
				dummy[k] = 0.0;
				for (int l=0; l<3; ++l)                                            // AP:29-06-20 //
				{
					dummy[k] += Iten[k * 3 + l]*Amat[j * 3 + l];               // AP:29-06-20 //
				}
				A1[i * nq + j] += Amat[i * 3 + k] * dummy[k];                      // AP:29-06-20 //
			}
			A1[j * nq + i] = A1[i * nq + j];
		}
	}

	/************/
	/* g MATRIX */
	/************/

	// Computed in DF //

	// 1. Calculate the covariant metric tensor
	for (int i = 0; i < nq; ++i)
	{
		for (int j = i; j < nq; ++j)
		{
			Gmat[i * nq + j] = 0.0;
			for (int k = 0; k < mol->getNAtoms(); ++k)
			{
				Gmat[i * nq + j] += m[k] * dca[(3 * k + 0) * nq + i] * dca[(3 * k + 0) * nq + j];
				Gmat[i * nq + j] += m[k] * dca[(3 * k + 1) * nq + i] * dca[(3 * k + 1) * nq + j];
				Gmat[i * nq + j] += m[k] * dca[(3 * k + 2) * nq + i] * dca[(3 * k + 2) * nq + j];
			}
			Gmat[i * nq + j] -= A1[i * nq + j];
			Gmat[j * nq + i] = Gmat[i * nq + j];
		}
	}
#ifdef DEBUG
	f.open("g_covariant_Matrix.dat", std::ios::out);
	for (int i = 0; i < nq; ++i)
	{
		for (int j = i; j < nq; ++j)
			f << i+1 << "\t" << j+1 << "\t" << std::scientific << std::setprecision(12) << Gmat[i * nq + j] << std::endl;
	}
	f.close();
#endif

	// 2. Invert the covariant metric tensor to obtain the contravariant metric tensor
	int INFO;
	double *invGmat = new double[nq * nq];
	for (int i = 0; i < nq; ++i)
	{
		for (int j = i; j < nq; ++j)
			invGmat[i * nq + j] = invGmat[j * nq + i] = Gmat[i * nq + j];
	}
#ifdef __MAGMA__
	TESTING_CHECK(magma_init());
	TESTING_CHECK(magma_dpotrf(MagmaLower, nq, invGmat, nq, &INFO));
	TESTING_CHECK(magma_dpotri(MagmaLower, nq, invGmat, nq, &INFO));
	TESTING_CHECK(magma_finalize());
#else
	char UPLO = 'L';
	dpotrf_(&UPLO, &nq, invGmat, &nq, &INFO);
	dpotri_(&UPLO, &nq, invGmat, &nq, &INFO);
#endif
	for (int i = 0; i < nq; ++i)
	{
		for (int j = i; j < nq; ++j)
			Gmat[i * nq + j] = Gmat[j * nq + i] = invGmat[i * nq + j];
	}
	delete(invGmat);

#ifdef DEBUG
	f.open("g_contravariant_Matrix.dat", std::ios::out);
	for (int i = 0; i < nq; ++i)
	{
		for (int j = i; j < nq; ++j)
			f << i+1 << "\t" << j+1 << "\t" << std::scientific << std::setprecision(12) << Gmat[i * nq + j] << std::endl;
	}
	f.close();
#endif

	/********/
	/* Exit */
	/********/
	
	return;
}

/***************************/
/* Methods to return stuff */
/***************************/

vectorOfDoubles buildAGI::getAmat(void)
{
	return Amat;
}

void buildAGI::getAmat(double *a)
{
	for (int i = 0; i < Amat.size(); ++i)
		a[i] = Amat[i];
	return;
}

vectorOfDoubles buildAGI::getGmat(void)
{
	return Gmat;
}

void buildAGI::getGmat(double *g)
{
	for (int i = 0; i < Gmat.size(); ++i)
		g[i] = Gmat[i];
	return;
}

vectorOfDoubles buildAGI::getIten(void)
{
	return Iten;
}

void buildAGI::getIten(double *it)
{
	for (int i = 0; i < Iten.size(); ++i)
		it[i] = Iten[i];
	return;
}

vectorOfDoubles buildAGI::getinvIten(void)
{
	return Iten;
}

void buildAGI::getinvIten(double *iit)
{
	for (int i = 0; i < InvIten.size(); ++i)
		iit[i] = InvIten[i];
	return;
}

vectorOfDoubles buildAGI::getCoordinatesInCF(void)
{
	return ca;
}

vectorOfDoubles buildAGI::getFirstDerivativesInCF(void)
{
	return dca;
}

