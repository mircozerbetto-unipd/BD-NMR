///////////////////////////////////////////////////////////
// This driver transforms the Variance-Covariance matrix //
// in Cartesian coordinates to internal coordinates      //
// based on the Z-Matrix associated to the molecule      //
// and a reference structure of the molecule             //
// ////////////////////////////////////////////////////////
#include "muEuler.h"
#include "ADEL/ADEL.h"

// Order of the derivatives to be calculated
const int Order_mu = 2;
// Size of the variable q
const int qCount = 2;
// The molecule and the atom coordinate ID must be global
molecule *ggMol;
int g_rH[3], g_rN[3], g_rC[3];
double betaDC;
// Define types for the active datatypes for the 3 kinds of internal coordinates
typedef ADVariable<Order_mu, qCount> ADReverse_v;
// The algorithms
template<typename T>T Algorithm_alphaD(const T (&q)[qCount])
{
	T xH = cartesian(ggMol, g_rH[0], q[0], q[1]);
	T yH = cartesian(ggMol, g_rH[1], q[0], q[1]);
	T zH = cartesian(ggMol, g_rH[2], q[0], q[1]);
	T xN = cartesian(ggMol, g_rN[0], q[0], q[1]);
	T yN = cartesian(ggMol, g_rN[1], q[0], q[1]);
	T zN = cartesian(ggMol, g_rN[2], q[0], q[1]);
	T d = xH - xN;
	T n = d * d;
	d = yH - yN;
	n = n + d * d;
	d = zH - zN;
	n = n + d * d;
	n = sqrt(n);
	xH = xH - xN;
	yH = yH - yN;
	zH = zH - zN;
	xH = xH / n;
	yH = yH / n;
	zH = zH / n;
	// Apply a small beta-tilt to avoid singularities (gimbal lock at beta = 0)
	double deltaB = 1.0e-10;
	T tY = yH - deltaB * zH;
	T tX = (1.0 - 0.5 * deltaB * deltaB) * xH - deltaB * zH;
	T a;
	a = atan2(tY, tX); // atan2 gives angle in -pi:pi interval, while alpha should be in 0:2pi
	a = a + M_PI;
	return a;
}
template<typename T>T Algorithm_CosBetaD(const T (&q)[qCount])
{
	T xH = cartesian(ggMol, g_rH[0], q[0], q[1]);
	T yH = cartesian(ggMol, g_rH[1], q[0], q[1]);
	T zH = cartesian(ggMol, g_rH[2], q[0], q[1]);
	T xN = cartesian(ggMol, g_rN[0], q[0], q[1]);
	T yN = cartesian(ggMol, g_rN[1], q[0], q[1]);
	T zN = cartesian(ggMol, g_rN[2], q[0], q[1]);
	T d = xH - xN;
	T n = d * d;
	d = yH - yN;
	n = n + d * d;
	d = zH - zN;
	n = n + d * d;
	n = sqrt(n);
	xH = xH - xN;
	yH = yH - yN;
	zH = zH - zN;
	xH = xH / n;
	yH = yH / n;
	zH = zH / n;

	// Apply a small beta-tilt to avoid singularities (gimbal lock at beta = 0)
	double deltaB = 1.0e-10;
	T tZ = deltaB * xH + (1.0 - 0.5 * deltaB * deltaB) * zH;

	// return cos(b) to handle the problematic case when zH = 1.0
	return tZ;
}
template<typename T>T Algorithm_alphaC(const T (&q)[qCount])
{
	T zD1 = cartesian(ggMol, g_rH[0], q[0], q[1]) - cartesian(ggMol, g_rN[0], q[0], q[1]);
	T zD2 = cartesian(ggMol, g_rH[1], q[0], q[1]) - cartesian(ggMol, g_rN[1], q[0], q[1]);
	T zD3 = cartesian(ggMol, g_rH[2], q[0], q[1]) - cartesian(ggMol, g_rN[2], q[0], q[1]);
	T n   = zD1 * zD1;
	n = n + zD2 * zD2;
	n = n + zD3 * zD3;
	n = sqrt(n);
	zD1 = zD1 / n;
	zD2 = zD2 / n;
	zD3 = zD3 / n;

	T v1 = cartesian(ggMol, g_rC[0], q[0], q[1]) - cartesian(ggMol, g_rN[0], q[0], q[1]);
	T v2 = cartesian(ggMol, g_rC[1], q[0], q[1]) - cartesian(ggMol, g_rN[1], q[0], q[1]);
	T v3 = cartesian(ggMol, g_rC[2], q[0], q[1]) - cartesian(ggMol, g_rN[2], q[0], q[1]);

	T yD1 = v2 * zD3 - v3 * zD2;
	T yD2 = v3 * zD1 - v1 * zD3;
	T yD3 = v1 * zD2 - v2 * zD1;
	n = yD1 * yD1;
	n = n + yD2 * yD2;
	n = n + yD3 * yD3;
	n = sqrt(n);
	yD1 = yD1 / n;
	yD2 = yD2 / n;
	yD3 = yD3 / n;

	T xD1 = yD2 * zD3 - yD3 * zD2;
	T xD2 = yD3 * zD1 - yD1 * zD3;
	T xD3 = yD1 * zD2 - yD2 * zD1;

	T Exx = xD1;
	T Exy = xD2;
	T Exz = xD3;
	T Eyx = yD1;
	T Eyy = yD2;
	T Eyz = yD3;
	T Ezx = zD1;
	T Ezy = zD2;
	T Ezz = zD3;

	T rN1 = cartesian(ggMol, g_rN[0], q[0], q[1]);
	T rN2 = cartesian(ggMol, g_rN[1], q[0], q[1]);
	T rN3 = cartesian(ggMol, g_rN[2], q[0], q[1]);

	v1 = zD1 - rN1; v2 = zD2 - rN2; v3 = zD3 - rN3;
	xD1 = Exx * v1 + Exy * v2 + Exz * v3;
	xD2 = Eyx * v1 + Eyy * v2 + Eyz * v3;
	xD3 = Ezx * v1 + Ezy * v2 + Ezz * v3;

	T cbD = cos(betaDC), sbD = sin(betaDC);

	Exx = cbD;
	Exy = 0.0;
	Exz = -sbD;
	Eyx = 0.0;
	Eyy = 1.0;
	Eyz = 0.0;
	Ezx = sbD;
	Ezy = 0.0;
	Ezz = cbD;

	zD1 = Exx * xD1 + Exy * xD2 + Exz * xD3 + rN1;
	zD2 = Eyx * xD1 + Eyy * xD2 + Eyz * xD3 + rN2;
	zD3 = Ezx * xD1 + Ezy * xD2 + Ezz * xD3 + rN3;

	// Apply a small beta-tilt to avoid singularities (gimbal lock at beta = 0)
	double deltaB = 1.0e-10;
	T tY = zD2 - deltaB * zD3;
	T tX = (1.0 - 0.5 * deltaB * deltaB) * zD1 - deltaB * zD3;

	T a;
	a = atan2(tY, tX); // atan2 gives angle in -pi:pi interval, while alpha should be in 0:2pi
	a = a + M_PI;
	return a;
}
template<typename T>T Algorithm_CosBetaC(const T (&q)[qCount])
{
	T zD1 = cartesian(ggMol, g_rH[0], q[0], q[1]) - cartesian(ggMol, g_rN[0], q[0], q[1]);
	T zD2 = cartesian(ggMol, g_rH[1], q[0], q[1]) - cartesian(ggMol, g_rN[1], q[0], q[1]);
	T zD3 = cartesian(ggMol, g_rH[2], q[0], q[1]) - cartesian(ggMol, g_rN[2], q[0], q[1]);
	T n   = zD1 * zD1;
	n = n + zD2 * zD2;
	n = n + zD3 * zD3;
	n = sqrt(n);
	zD1 = zD1 / n;
	zD2 = zD2 / n;
	zD3 = zD3 / n;

	T v1 = cartesian(ggMol, g_rC[0], q[0], q[1]) - cartesian(ggMol, g_rN[0], q[0], q[1]);
	T v2 = cartesian(ggMol, g_rC[1], q[0], q[1]) - cartesian(ggMol, g_rN[1], q[0], q[1]);
	T v3 = cartesian(ggMol, g_rC[2], q[0], q[1]) - cartesian(ggMol, g_rN[2], q[0], q[1]);

	T yD1 = v2 * zD3 - v3 * zD2;
	T yD2 = v3 * zD1 - v1 * zD3;
	T yD3 = v1 * zD2 - v2 * zD1;
	n = yD1 * yD1;
	n = n + yD2 * yD2;
	n = n + yD3 * yD3;
	n = sqrt(n);
	yD1 = yD1 / n;
	yD2 = yD2 / n;
	yD3 = yD3 / n;

	T xD1 = yD2 * zD3 - yD3 * zD2;
	T xD2 = yD3 * zD1 - yD1 * zD3;
	T xD3 = yD1 * zD2 - yD2 * zD1;

	T Exx = xD1;
	T Exy = xD2;
	T Exz = xD3;
	T Eyx = yD1;
	T Eyy = yD2;
	T Eyz = yD3;
	T Ezx = zD1;
	T Ezy = zD2;
	T Ezz = zD3;

	T rN1 = cartesian(ggMol, g_rN[0], q[0], q[1]);
	T rN2 = cartesian(ggMol, g_rN[1], q[0], q[1]);
	T rN3 = cartesian(ggMol, g_rN[2], q[0], q[1]);

	v1 = zD1 - rN1; v2 = zD2 - rN2; v3 = zD3 - rN3;
	xD1 = Exx * v1 + Exy * v2 + Exz * v3;
	xD2 = Eyx * v1 + Eyy * v2 + Eyz * v3;
	xD3 = Ezx * v1 + Ezy * v2 + Ezz * v3;

	T cbD = cos(betaDC), sbD = sin(betaDC);

	Exx = cbD;
	Exy = 0.0;
	Exz = -sbD;
	Eyx = 0.0;
	Eyy = 1.0;
	Eyz = 0.0;
	Ezx = sbD;
	Ezy = 0.0;
	Ezz = cbD;

	zD1 = Exx * xD1 + Exy * xD2 + Exz * xD3 + rN1;
	zD2 = Eyx * xD1 + Eyy * xD2 + Eyz * xD3 + rN2;
	zD3 = Ezx * xD1 + Ezy * xD2 + Ezz * xD3 + rN3;

	// Apply a small beta-tilt to avoid singularities (gimbal lock at beta = 0)
	double deltaB = 1.0e-10;
	T tZ = deltaB * zD1 + (1.0 - 0.5 * deltaB * deltaB) * zD3;

	// return cos(b) to handle the problematic case when zH = 1.0
	return tZ;
}
// Get eulerD angles
void eulerD(int idH, int idN, int idC, int q1, int q2, molecule *mol, double *a, double *afd, double *asd, double *b, double *bfd, double *bsd)
{
	ggMol = mol;
	int na = ggMol->getNAtoms();
	int nx = 3 * na;

	g_rH[0] = (idH - 1) * 3 + 0; g_rH[1] = (idH - 1) * 3 + 1; g_rH[2] = (idH - 1) * 3 + 2;
	g_rN[0] = (idN - 1) * 3 + 0; g_rN[1] = (idN - 1) * 3 + 1; g_rN[2] = (idN - 1) * 3 + 2;
	g_rC[0] = (idC - 1) * 3 + 0; g_rC[1] = (idC - 1) * 3 + 1; g_rC[2] = (idC - 1) * 3 + 2;
	
	ADReverse_v Q[qCount];
	ADReverse_v AD_ALPHA, AD_COSBETA;
	
	Q[0].Value = q1;
	Q[0].ID    = 0;
	Q[0].Base  = true;

	Q[1].Value = q2;
	Q[1].ID    = 1;
	Q[1].Base  = true;
	AD_ALPHA  = Algorithm_alphaD<ADReverse_v>(Q);
	a[0]   = AD_ALPHA.Value;
	afd[0] = AD_ALPHA.Data[0];
	afd[1] = AD_ALPHA.Data[1];
	asd[0] = AD_ALPHA.Data[2];
	asd[1] = AD_ALPHA.Data[3];
	asd[2] = AD_ALPHA.Data[4];
	asd[3] = AD_ALPHA.Data[5];

	AD_COSBETA   = Algorithm_CosBetaD<ADReverse_v>(Q);
	b[0]   = AD_COSBETA.Value;
	bfd[0] = AD_COSBETA.Data[0];
	bfd[1] = AD_COSBETA.Data[1];
	bsd[0] = AD_COSBETA.Data[2];
	bsd[1] = AD_COSBETA.Data[3];
	bsd[2] = AD_COSBETA.Data[4];
	bsd[3] = AD_COSBETA.Data[5];

	/*std::cout << a[0] << ", " << afd[0] << ", " << afd[1] << ", " << asd[0] << ", " << asd[1] << ", " << asd[2] << ", " << asd[3] << std::endl;
	std::cout << b[0] << ", " << bfd[0] << ", " << bfd[1] << ", " << bsd[0] << ", " << bsd[1] << ", " << bsd[2] << ", " << bsd[3] << std::endl;*/
	return;
}
// Get euleg_rC angles
void eulerC(int idH, int idN, int idC, int q1, int q2, molecule *mol, double bDC, double *a, double *afd, double *asd, double *b, double *bfd, double *bsd)
{
	ggMol = mol;
	betaDC = bDC;
	int na = ggMol->getNAtoms();
	int nx = 3 * na;
	
	g_rH[0] = (idH - 1) * 3 + 0; g_rH[1] = (idH - 1) * 3 + 1; g_rH[2] = (idH - 1) * 3 + 2;
	g_rN[0] = (idN - 1) * 3 + 0; g_rN[1] = (idN - 1) * 3 + 1; g_rN[2] = (idN - 1) * 3 + 2;
	g_rC[0] = (idC - 1) * 3 + 0; g_rC[1] = (idC - 1) * 3 + 1; g_rC[2] = (idC - 1) * 3 + 2;
	
	ADReverse_v Q[qCount];
	ADReverse_v AD_ALPHA, AD_COSBETA;
	
	Q[0].Value = q1;
	Q[0].ID   = 0;
	Q[0].Base = true;

	Q[1].Value = q2;
	Q[1].ID   = 1;
	Q[1].Base = true;

	AD_ALPHA  = Algorithm_alphaC<ADReverse_v>(Q);
	a[0]   = AD_ALPHA.Value;
	afd[0] = AD_ALPHA.Data[0];
	afd[1] = AD_ALPHA.Data[1];
	asd[0] = AD_ALPHA.Data[2];
	asd[1] = AD_ALPHA.Data[3];
	asd[2] = AD_ALPHA.Data[4];
	asd[3] = AD_ALPHA.Data[5];
	
	AD_COSBETA   = Algorithm_CosBetaC<ADReverse_v>(Q);
	b[0]   = AD_COSBETA.Value;
	bfd[0] = AD_COSBETA.Data[0];
	bfd[1] = AD_COSBETA.Data[1];
	bsd[0] = AD_COSBETA.Data[2];
	bsd[1] = AD_COSBETA.Data[3];
	bsd[2] = AD_COSBETA.Data[4];
	bsd[3] = AD_COSBETA.Data[5];
	
	return;
}

