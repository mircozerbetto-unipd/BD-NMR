///////////////////////////////////////////////////////////
// This driver transforms the Variance-Covariance matrix //
// in Cartesian coordinates to internal coordinates      //
// based on the Z-Matrix associated to the molecule      //
// and a reference structure of the molecule             //
// ////////////////////////////////////////////////////////
#include "ADEL/ADEL.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cmath>

#include "zmat.h"

// Order of the derivatives to be calculated
const int Order = 1;
// Size of the variables on which distance (d), angle (a), and torsion (t) depend on
const int dCount = 6;
const int aCount = 9;
const int tCount = 12;
// Define types for the active datatypes for the 3 kinds of internal coordinates
typedef ADVariable<Order, dCount> ADReverse_d;
typedef ADVariable<Order, aCount> ADReverse_a;
typedef ADVariable<Order, tCount> ADReverse_t;
// These are the templates algorithms that calculate d, a, or t from Cartesian coordinates
/*********************************/
/* Bond distance A1-A2           */
/* d = | r2 - r1 |               */
/* r1 = [x[0] x[1] x[2]]^tr      */
/* r2 = [x[3] x[4] x[5]]^tr      */
/*********************************/
template<typename T>T Algorithm_d(const T (&x)[dCount])
{
	T d = (x[3] - x[0]) * (x[3] - x[0]);
	  d = d + (x[4] - x[1]) * (x[4] - x[1]);
	  d = d + (x[5] - x[2]) * (x[5] - x[2]);
	return sqrt(d);
}
/*********************************/
/* Bond angle A1-A2-A3           */
/* a = acos(b . c / |b||c|)      */
/* b = r1 - r2                   */
/* c = r3 - r2                   */
/* r1 = [x[0] x[1] x[2]]^tr      */
/* r2 = [x[3] x[4] x[5]]^tr      */
/* r3 = [x[6] x[7] x[8]]^tr      */
/*********************************/
template<typename T>T Algorithm_a(const T (&x)[aCount])
{
	T b[3], c[3];
	b[0] = x[0] - x[3]; b[1] = x[1] - x[4]; b[2] = x[2] - x[5];
	c[0] = x[6] - x[3]; c[1] = x[7] - x[4]; c[2] = x[8] - x[5];
	T nb = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
	T nc = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
	T a = (b[0] * c[0] + b[1] * c[1] + b[2] * c[2]) / (nb * nc);
	return acos(a);
}
/********************************/
/* Torsion angle A1-A2-A3-A4    */
/* t = acos(qa . qb / |qa||qb|) */
/* qa = b x a                   */
/* qb = c x b                   */
/* a = r2 - r1                  */
/* b = r3 - r2                  */
/* c = r4 - r3                  */
/* r1 = [x[0] x[ 1] x[ 2]]^tr   */
/* r2 = [x[3] x[ 4] x[ 5]]^tr   */
/* r3 = [x[6] x[ 7] x[ 8]]^tr   */
/* r3 = [x[9] x[10] x[11]]^tr   */
/*********************************/
template<typename T>T Algorithm_t(const T (&x)[tCount])
{
	T a[3], b[3], c[3], qa[3], qb[3];

	a[0] = x[3] - x[0]; a[1] = x[4]  - x[1]; a[2] = x[5]  - x[2];
	b[0] = x[6] - x[3]; b[1] = x[7]  - x[4]; b[2] = x[8]  - x[5];
	c[0] = x[9] - x[6]; c[1] = x[10] - x[7]; c[2] = x[11] - x[8];

	qa[0] = b[1] * a[2] - b[2] * a[1];
	qa[1] = b[2] * a[0] - b[0] * a[2];
	qa[2] = b[0] * a[1] - b[1] * a[0];

	qb[0] = c[1] * b[2] - c[2] * b[1];
	qb[1] = c[2] * b[0] - c[0] * b[2];
	qb[2] = c[0] * b[1] - c[1] * b[0];

	T nqa = sqrt(qa[0] * qa[0] + qa[1] * qa[1] + qa[2] * qa[2]);
	T nqb = sqrt(qb[0] * qb[0] + qb[1] * qb[1] + qb[2] * qb[2]);

	T t = (qa[0] * qb[0] + qa[1] * qb[1] + qa[2] * qb[2]) / (nqa * nqb);
	return acos(t);
}
//////////////////////////////////
// distance-distance covariance //
//////////////////////////////////
double cov_dd(atom At1, atom At2, atom Bt1, atom Bt2, int dimM, double *M)
{
	ADReverse_d Xd1[dCount], Xd2[dCount];
    	ADReverse_d D1, D2;

	Xd1[0].Value = At1.x(); Xd1[0].ID = 0; Xd1[0].Base = true;
	Xd1[1].Value = At1.y(); Xd1[1].ID = 1; Xd1[1].Base = true;
	Xd1[2].Value = At1.z(); Xd1[2].ID = 2; Xd1[2].Base = true;
	Xd1[3].Value = At2.x(); Xd1[3].ID = 3; Xd1[3].Base = true;
	Xd1[4].Value = At2.y(); Xd1[4].ID = 4; Xd1[4].Base = true;
	Xd1[5].Value = At2.z(); Xd1[5].ID = 5; Xd1[5].Base = true;
    	D1 = Algorithm_d<ADReverse_d>(Xd1);

	Xd2[0].Value = Bt1.x(); Xd2[0].ID = 0; Xd2[0].Base = true;
	Xd2[1].Value = Bt1.y(); Xd2[1].ID = 1; Xd2[1].Base = true;
	Xd2[2].Value = Bt1.z(); Xd2[2].ID = 2; Xd2[2].Base = true;
	Xd2[3].Value = Bt2.x(); Xd2[3].ID = 3; Xd2[3].Base = true;
	Xd2[4].Value = Bt2.y(); Xd2[4].ID = 4; Xd2[4].Base = true;
	Xd2[5].Value = Bt2.z(); Xd2[5].ID = 5; Xd2[5].Base = true;
    	D2 = Algorithm_d<ADReverse_d>(Xd2);

	double Cdd = 0.0;
	double Cst[9];

	// A1--B1
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[0] * (Cst[0] * D2.Data[0] + Cst[1] * D2.Data[1] + Cst[2] * D2.Data[2]);
	Cdd += D1.Data[1] * (Cst[3] * D2.Data[0] + Cst[4] * D2.Data[1] + Cst[5] * D2.Data[2]);
	Cdd += D1.Data[2] * (Cst[6] * D2.Data[0] + Cst[7] * D2.Data[1] + Cst[8] * D2.Data[2]);

	// A1--B2
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[0] * (Cst[0] * D2.Data[3] + Cst[1] * D2.Data[4] + Cst[2] * D2.Data[5]);
	Cdd += D1.Data[1] * (Cst[3] * D2.Data[3] + Cst[4] * D2.Data[4] + Cst[5] * D2.Data[5]);
	Cdd += D1.Data[2] * (Cst[6] * D2.Data[3] + Cst[7] * D2.Data[4] + Cst[8] * D2.Data[5]);

	// A2--B1
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[3] * (Cst[0] * D2.Data[0] + Cst[1] * D2.Data[1] + Cst[2] * D2.Data[2]);
	Cdd += D1.Data[4] * (Cst[3] * D2.Data[0] + Cst[4] * D2.Data[1] + Cst[5] * D2.Data[2]);
	Cdd += D1.Data[5] * (Cst[6] * D2.Data[0] + Cst[7] * D2.Data[1] + Cst[8] * D2.Data[2]);

	// A2--B2
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[3] * (Cst[0] * D2.Data[3] + Cst[1] * D2.Data[4] + Cst[2] * D2.Data[5]);
	Cdd += D1.Data[4] * (Cst[3] * D2.Data[3] + Cst[4] * D2.Data[4] + Cst[5] * D2.Data[5]);
	Cdd += D1.Data[5] * (Cst[6] * D2.Data[3] + Cst[7] * D2.Data[4] + Cst[8] * D2.Data[5]);

	return Cdd;
}
///////////////////////////////
// distance-angle covariance //
///////////////////////////////
double cov_da(atom At1, atom At2, atom Bt1, atom Bt2, atom Bt3, int dimM, double *M)
{
	ADReverse_d Xd1[dCount], D1; 
	ADReverse_a Xa2[aCount], A2;

	Xd1[0].Value = At1.x(); Xd1[0].ID = 0; Xd1[0].Base = true;
	Xd1[1].Value = At1.y(); Xd1[1].ID = 1; Xd1[1].Base = true;
	Xd1[2].Value = At1.z(); Xd1[2].ID = 2; Xd1[2].Base = true;
	Xd1[3].Value = At2.x(); Xd1[3].ID = 3; Xd1[3].Base = true;
	Xd1[4].Value = At2.y(); Xd1[4].ID = 4; Xd1[4].Base = true;
	Xd1[5].Value = At2.z(); Xd1[5].ID = 5; Xd1[5].Base = true;
    	D1 = Algorithm_d<ADReverse_d>(Xd1);

	Xa2[0].Value = Bt1.x(); Xa2[0].ID = 0; Xa2[0].Base = true;
	Xa2[1].Value = Bt1.y(); Xa2[1].ID = 1; Xa2[1].Base = true;
	Xa2[2].Value = Bt1.z(); Xa2[2].ID = 2; Xa2[2].Base = true;
	Xa2[3].Value = Bt2.x(); Xa2[3].ID = 3; Xa2[3].Base = true;
	Xa2[4].Value = Bt2.y(); Xa2[4].ID = 4; Xa2[4].Base = true;
	Xa2[5].Value = Bt2.z(); Xa2[5].ID = 5; Xa2[5].Base = true;
	Xa2[6].Value = Bt3.x(); Xa2[6].ID = 6; Xa2[6].Base = true;
	Xa2[7].Value = Bt3.y(); Xa2[7].ID = 7; Xa2[7].Base = true;
	Xa2[8].Value = Bt3.z(); Xa2[8].ID = 8; Xa2[8].Base = true;
    	A2 = Algorithm_a<ADReverse_a>(Xa2);

	double Cdd = 0.0;
	double Cst[9];

	// A1--B1
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[0] * (Cst[0] * A2.Data[0] + Cst[1] * A2.Data[1] + Cst[2] * A2.Data[2]);
	Cdd += D1.Data[1] * (Cst[3] * A2.Data[0] + Cst[4] * A2.Data[1] + Cst[5] * A2.Data[2]);
	Cdd += D1.Data[2] * (Cst[6] * A2.Data[0] + Cst[7] * A2.Data[1] + Cst[8] * A2.Data[2]);

	// A1--B2
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[0] * (Cst[0] * A2.Data[3] + Cst[1] * A2.Data[4] + Cst[2] * A2.Data[5]);
	Cdd += D1.Data[1] * (Cst[3] * A2.Data[3] + Cst[4] * A2.Data[4] + Cst[5] * A2.Data[5]);
	Cdd += D1.Data[2] * (Cst[6] * A2.Data[3] + Cst[7] * A2.Data[4] + Cst[8] * A2.Data[5]);

	// A1--B3
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[0] * (Cst[0] * A2.Data[6] + Cst[1] * A2.Data[7] + Cst[2] * A2.Data[8]);
	Cdd += D1.Data[1] * (Cst[3] * A2.Data[6] + Cst[4] * A2.Data[7] + Cst[5] * A2.Data[8]);
	Cdd += D1.Data[2] * (Cst[6] * A2.Data[6] + Cst[7] * A2.Data[7] + Cst[8] * A2.Data[8]);

	// A2--B1
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[3] * (Cst[0] * A2.Data[0] + Cst[1] * A2.Data[1] + Cst[2] * A2.Data[2]);
	Cdd += D1.Data[4] * (Cst[3] * A2.Data[0] + Cst[4] * A2.Data[1] + Cst[5] * A2.Data[2]);
	Cdd += D1.Data[5] * (Cst[6] * A2.Data[0] + Cst[7] * A2.Data[1] + Cst[8] * A2.Data[2]);

	// A2--B2
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[3] * (Cst[0] * A2.Data[3] + Cst[1] * A2.Data[4] + Cst[2] * A2.Data[5]);
	Cdd += D1.Data[4] * (Cst[3] * A2.Data[3] + Cst[4] * A2.Data[4] + Cst[5] * A2.Data[5]);
	Cdd += D1.Data[5] * (Cst[6] * A2.Data[3] + Cst[7] * A2.Data[4] + Cst[8] * A2.Data[5]);

	// A2--B3
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[3] * (Cst[0] * A2.Data[6] + Cst[1] * A2.Data[7] + Cst[2] * A2.Data[8]);
	Cdd += D1.Data[4] * (Cst[3] * A2.Data[6] + Cst[4] * A2.Data[7] + Cst[5] * A2.Data[8]);
	Cdd += D1.Data[5] * (Cst[6] * A2.Data[6] + Cst[7] * A2.Data[7] + Cst[8] * A2.Data[8]);

	return Cdd;
}
/////////////////////////////////
// distance-torsion covariance //
/////////////////////////////////
double cov_dt(atom At1, atom At2, atom Bt1, atom Bt2, atom Bt3, atom Bt4, int dimM, double *M)
{
	ADReverse_d Xd1[dCount], D1;
    	ADReverse_t Xt2[tCount], T2;

	Xd1[0].Value = At1.x(); Xd1[0].ID = 0; Xd1[0].Base = true;
	Xd1[1].Value = At1.y(); Xd1[1].ID = 1; Xd1[1].Base = true;
	Xd1[2].Value = At1.z(); Xd1[2].ID = 2; Xd1[2].Base = true;
	Xd1[3].Value = At2.x(); Xd1[3].ID = 3; Xd1[3].Base = true;
	Xd1[4].Value = At2.y(); Xd1[4].ID = 4; Xd1[4].Base = true;
	Xd1[5].Value = At2.z(); Xd1[5].ID = 5; Xd1[5].Base = true;
    	D1 = Algorithm_d<ADReverse_d>(Xd1);

	Xt2[ 0].Value = Bt1.x(); Xt2[ 0].ID =  0; Xt2[ 0].Base = true;
	Xt2[ 1].Value = Bt1.y(); Xt2[ 1].ID =  1; Xt2[ 1].Base = true;
	Xt2[ 2].Value = Bt1.z(); Xt2[ 2].ID =  2; Xt2[ 2].Base = true;
	Xt2[ 3].Value = Bt2.x(); Xt2[ 3].ID =  3; Xt2[ 3].Base = true;
	Xt2[ 4].Value = Bt2.y(); Xt2[ 4].ID =  4; Xt2[ 4].Base = true;
	Xt2[ 5].Value = Bt2.z(); Xt2[ 5].ID =  5; Xt2[ 5].Base = true;
	Xt2[ 6].Value = Bt3.x(); Xt2[ 6].ID =  6; Xt2[ 6].Base = true;
	Xt2[ 7].Value = Bt3.y(); Xt2[ 7].ID =  7; Xt2[ 7].Base = true;
	Xt2[ 8].Value = Bt3.z(); Xt2[ 8].ID =  8; Xt2[ 8].Base = true;
	Xt2[ 9].Value = Bt4.x(); Xt2[ 9].ID =  9; Xt2[ 9].Base = true;
	Xt2[10].Value = Bt4.y(); Xt2[10].ID = 10; Xt2[10].Base = true;
	Xt2[11].Value = Bt4.z(); Xt2[11].ID = 11; Xt2[11].Base = true;
    	T2 = Algorithm_t<ADReverse_t>(Xt2);

	double Cdd = 0.0;
	double Cst[9];

	// A1--B1
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[0] * (Cst[0] * T2.Data[0] + Cst[1] * T2.Data[1] + Cst[2] * T2.Data[2]);
	Cdd += D1.Data[1] * (Cst[3] * T2.Data[0] + Cst[4] * T2.Data[1] + Cst[5] * T2.Data[2]);
	Cdd += D1.Data[2] * (Cst[6] * T2.Data[0] + Cst[7] * T2.Data[1] + Cst[8] * T2.Data[2]);

	// A1--B2
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[0] * (Cst[0] * T2.Data[3] + Cst[1] * T2.Data[4] + Cst[2] * T2.Data[5]);
	Cdd += D1.Data[1] * (Cst[3] * T2.Data[3] + Cst[4] * T2.Data[4] + Cst[5] * T2.Data[5]);
	Cdd += D1.Data[2] * (Cst[6] * T2.Data[3] + Cst[7] * T2.Data[4] + Cst[8] * T2.Data[5]);

	// A1--B3
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[0] * (Cst[0] * T2.Data[6] + Cst[1] * T2.Data[7] + Cst[2] * T2.Data[8]);
	Cdd += D1.Data[1] * (Cst[3] * T2.Data[6] + Cst[4] * T2.Data[7] + Cst[5] * T2.Data[8]);
	Cdd += D1.Data[2] * (Cst[6] * T2.Data[6] + Cst[7] * T2.Data[7] + Cst[8] * T2.Data[8]);

	// A1--B4
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[0] * (Cst[0] * T2.Data[9] + Cst[1] * T2.Data[10] + Cst[2] * T2.Data[11]);
	Cdd += D1.Data[1] * (Cst[3] * T2.Data[9] + Cst[4] * T2.Data[10] + Cst[5] * T2.Data[11]);
	Cdd += D1.Data[2] * (Cst[6] * T2.Data[9] + Cst[7] * T2.Data[10] + Cst[8] * T2.Data[11]);

	// A2--B1
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[3] * (Cst[0] * T2.Data[0] + Cst[1] * T2.Data[1] + Cst[2] * T2.Data[2]);
	Cdd += D1.Data[4] * (Cst[3] * T2.Data[0] + Cst[4] * T2.Data[1] + Cst[5] * T2.Data[2]);
	Cdd += D1.Data[5] * (Cst[6] * T2.Data[0] + Cst[7] * T2.Data[1] + Cst[8] * T2.Data[2]);

	// A2--B2
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[3] * (Cst[0] * T2.Data[3] + Cst[1] * T2.Data[4] + Cst[2] * T2.Data[5]);
	Cdd += D1.Data[4] * (Cst[3] * T2.Data[3] + Cst[4] * T2.Data[4] + Cst[5] * T2.Data[5]);
	Cdd += D1.Data[5] * (Cst[6] * T2.Data[3] + Cst[7] * T2.Data[4] + Cst[8] * T2.Data[5]);

	// A2--B3
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[3] * (Cst[0] * T2.Data[6] + Cst[1] * T2.Data[7] + Cst[2] * T2.Data[8]);
	Cdd += D1.Data[4] * (Cst[3] * T2.Data[6] + Cst[4] * T2.Data[7] + Cst[5] * T2.Data[8]);
	Cdd += D1.Data[5] * (Cst[6] * T2.Data[6] + Cst[7] * T2.Data[7] + Cst[8] * T2.Data[8]);

	// A2--B4
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cdd += D1.Data[3] * (Cst[0] * T2.Data[9] + Cst[1] * T2.Data[10] + Cst[2] * T2.Data[11]);
	Cdd += D1.Data[4] * (Cst[3] * T2.Data[9] + Cst[4] * T2.Data[10] + Cst[5] * T2.Data[11]);
	Cdd += D1.Data[5] * (Cst[6] * T2.Data[9] + Cst[7] * T2.Data[10] + Cst[8] * T2.Data[11]);

	return Cdd;
}
////////////////////////////
// angle-angle covariance //
////////////////////////////
double cov_aa(atom At1, atom At2, atom At3, atom Bt1, atom Bt2, atom Bt3, int dimM, double *M)
{
	ADReverse_a Xa1[aCount], Xa2[aCount];
    	ADReverse_a A1, A2;

	Xa1[0].Value = At1.x(); Xa1[0].ID = 0; Xa1[0].Base = true;
	Xa1[1].Value = At1.y(); Xa1[1].ID = 1; Xa1[1].Base = true;
	Xa1[2].Value = At1.z(); Xa1[2].ID = 2; Xa1[2].Base = true;
	Xa1[3].Value = At2.x(); Xa1[3].ID = 3; Xa1[3].Base = true;
	Xa1[4].Value = At2.y(); Xa1[4].ID = 4; Xa1[4].Base = true;
	Xa1[5].Value = At2.z(); Xa1[5].ID = 5; Xa1[5].Base = true;
	Xa1[6].Value = At3.x(); Xa1[6].ID = 6; Xa1[6].Base = true;
	Xa1[7].Value = At3.y(); Xa1[7].ID = 7; Xa1[7].Base = true;
	Xa1[8].Value = At3.z(); Xa1[8].ID = 8; Xa1[8].Base = true;
    	A1 = Algorithm_a<ADReverse_a>(Xa1);

	Xa2[0].Value = Bt1.x(); Xa2[0].ID = 0; Xa2[0].Base = true;
	Xa2[1].Value = Bt1.y(); Xa2[1].ID = 1; Xa2[1].Base = true;
	Xa2[2].Value = Bt1.z(); Xa2[2].ID = 2; Xa2[2].Base = true;
	Xa2[3].Value = Bt2.x(); Xa2[3].ID = 3; Xa2[3].Base = true;
	Xa2[4].Value = Bt2.y(); Xa2[4].ID = 4; Xa2[4].Base = true;
	Xa2[5].Value = Bt2.z(); Xa2[5].ID = 5; Xa2[5].Base = true;
	Xa2[6].Value = Bt3.x(); Xa2[6].ID = 6; Xa2[6].Base = true;
	Xa2[7].Value = Bt3.y(); Xa2[7].ID = 7; Xa2[7].Base = true;
	Xa2[8].Value = Bt3.z(); Xa2[8].ID = 8; Xa2[8].Base = true;
    	A2 = Algorithm_a<ADReverse_a>(Xa2);

	double Cdd = 0.0;
	double Cst[9];

	// A1--B1
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[0] * (Cst[0] * A2.Data[0] + Cst[1] * A2.Data[1] + Cst[2] * A2.Data[2]);
	Cdd += A1.Data[1] * (Cst[3] * A2.Data[0] + Cst[4] * A2.Data[1] + Cst[5] * A2.Data[2]);
	Cdd += A1.Data[2] * (Cst[6] * A2.Data[0] + Cst[7] * A2.Data[1] + Cst[8] * A2.Data[2]);

	// A1--B2
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[0] * (Cst[0] * A2.Data[3] + Cst[1] * A2.Data[4] + Cst[2] * A2.Data[5]);
	Cdd += A1.Data[1] * (Cst[3] * A2.Data[3] + Cst[4] * A2.Data[4] + Cst[5] * A2.Data[5]);
	Cdd += A1.Data[2] * (Cst[6] * A2.Data[3] + Cst[7] * A2.Data[4] + Cst[8] * A2.Data[5]);

	// A1--B3
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[0] * (Cst[0] * A2.Data[6] + Cst[1] * A2.Data[7] + Cst[2] * A2.Data[8]);
	Cdd += A1.Data[1] * (Cst[3] * A2.Data[6] + Cst[4] * A2.Data[7] + Cst[5] * A2.Data[8]);
	Cdd += A1.Data[2] * (Cst[6] * A2.Data[6] + Cst[7] * A2.Data[7] + Cst[8] * A2.Data[8]);

	// A2--B1
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[3] * (Cst[0] * A2.Data[0] + Cst[1] * A2.Data[1] + Cst[2] * A2.Data[2]);
	Cdd += A1.Data[4] * (Cst[3] * A2.Data[0] + Cst[4] * A2.Data[1] + Cst[5] * A2.Data[2]);
	Cdd += A1.Data[5] * (Cst[6] * A2.Data[0] + Cst[7] * A2.Data[1] + Cst[8] * A2.Data[2]);

	// A2--B2
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[3] * (Cst[0] * A2.Data[3] + Cst[1] * A2.Data[4] + Cst[2] * A2.Data[5]);
	Cdd += A1.Data[4] * (Cst[3] * A2.Data[3] + Cst[4] * A2.Data[4] + Cst[5] * A2.Data[5]);
	Cdd += A1.Data[5] * (Cst[6] * A2.Data[3] + Cst[7] * A2.Data[4] + Cst[8] * A2.Data[5]);

	// A2--B3
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[3] * (Cst[0] * A2.Data[6] + Cst[1] * A2.Data[7] + Cst[2] * A2.Data[8]);
	Cdd += A1.Data[4] * (Cst[3] * A2.Data[6] + Cst[4] * A2.Data[7] + Cst[5] * A2.Data[8]);
	Cdd += A1.Data[5] * (Cst[6] * A2.Data[6] + Cst[7] * A2.Data[7] + Cst[8] * A2.Data[8]);

	// A3--B1
	Cst[0] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[1] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[2] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[4] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[5] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[7] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[8] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[6] * (Cst[0] * A2.Data[0] + Cst[1] * A2.Data[1] + Cst[2] * A2.Data[2]);
	Cdd += A1.Data[7] * (Cst[3] * A2.Data[0] + Cst[4] * A2.Data[1] + Cst[5] * A2.Data[2]);
	Cdd += A1.Data[8] * (Cst[6] * A2.Data[0] + Cst[7] * A2.Data[1] + Cst[8] * A2.Data[2]);

	// A3--B2
	Cst[0] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[1] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[2] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[4] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[5] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[7] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[8] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[6] * (Cst[0] * A2.Data[3] + Cst[1] * A2.Data[4] + Cst[2] * A2.Data[5]);
	Cdd += A1.Data[7] * (Cst[3] * A2.Data[3] + Cst[4] * A2.Data[4] + Cst[5] * A2.Data[5]);
	Cdd += A1.Data[8] * (Cst[6] * A2.Data[3] + Cst[7] * A2.Data[4] + Cst[8] * A2.Data[5]);

	// A3--B3
	Cst[0] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[1] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[2] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[4] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[5] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[7] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[8] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[6] * (Cst[0] * A2.Data[6] + Cst[1] * A2.Data[7] + Cst[2] * A2.Data[8]);
	Cdd += A1.Data[7] * (Cst[3] * A2.Data[6] + Cst[4] * A2.Data[7] + Cst[5] * A2.Data[8]);
	Cdd += A1.Data[8] * (Cst[6] * A2.Data[6] + Cst[7] * A2.Data[7] + Cst[8] * A2.Data[8]);

	return Cdd;
}
//////////////////////////////
// angle-torsion covariance //
//////////////////////////////
double cov_at(atom At1, atom At2, atom At3, atom Bt1, atom Bt2, atom Bt3, atom Bt4, int dimM, double *M)
{
	ADReverse_a Xa1[aCount], A1;
    	ADReverse_t Xt2[tCount], T2;

	Xa1[0].Value = At1.x(); Xa1[0].ID = 0; Xa1[0].Base = true;
	Xa1[1].Value = At1.y(); Xa1[1].ID = 1; Xa1[1].Base = true;
	Xa1[2].Value = At1.z(); Xa1[2].ID = 2; Xa1[2].Base = true;
	Xa1[3].Value = At2.x(); Xa1[3].ID = 3; Xa1[3].Base = true;
	Xa1[4].Value = At2.y(); Xa1[4].ID = 4; Xa1[4].Base = true;
	Xa1[5].Value = At2.z(); Xa1[5].ID = 5; Xa1[5].Base = true;
	Xa1[6].Value = At3.x(); Xa1[6].ID = 6; Xa1[6].Base = true;
	Xa1[7].Value = At3.y(); Xa1[7].ID = 7; Xa1[7].Base = true;
	Xa1[8].Value = At3.z(); Xa1[8].ID = 8; Xa1[8].Base = true;
    	A1 = Algorithm_a<ADReverse_a>(Xa1);

	Xt2[ 0].Value = Bt1.x(); Xt2[ 0].ID =  0; Xt2[ 0].Base = true;
	Xt2[ 1].Value = Bt1.y(); Xt2[ 1].ID =  1; Xt2[ 1].Base = true;
	Xt2[ 2].Value = Bt1.z(); Xt2[ 2].ID =  2; Xt2[ 2].Base = true;
	Xt2[ 3].Value = Bt2.x(); Xt2[ 3].ID =  3; Xt2[ 3].Base = true;
	Xt2[ 4].Value = Bt2.y(); Xt2[ 4].ID =  4; Xt2[ 4].Base = true;
	Xt2[ 5].Value = Bt2.z(); Xt2[ 5].ID =  5; Xt2[ 5].Base = true;
	Xt2[ 6].Value = Bt3.x(); Xt2[ 6].ID =  6; Xt2[ 6].Base = true;
	Xt2[ 7].Value = Bt3.y(); Xt2[ 7].ID =  7; Xt2[ 7].Base = true;
	Xt2[ 8].Value = Bt3.z(); Xt2[ 8].ID =  8; Xt2[ 8].Base = true;
	Xt2[ 9].Value = Bt4.x(); Xt2[ 9].ID =  9; Xt2[ 9].Base = true;
	Xt2[10].Value = Bt4.y(); Xt2[10].ID = 10; Xt2[10].Base = true;
	Xt2[11].Value = Bt4.z(); Xt2[11].ID = 11; Xt2[11].Base = true;
    	T2 = Algorithm_t<ADReverse_t>(Xt2);

	double Cdd = 0.0;
	double Cst[9];

	// A1--B1
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[0] * (Cst[0] * T2.Data[0] + Cst[1] * T2.Data[1] + Cst[2] * T2.Data[2]);
	Cdd += A1.Data[1] * (Cst[3] * T2.Data[0] + Cst[4] * T2.Data[1] + Cst[5] * T2.Data[2]);
	Cdd += A1.Data[2] * (Cst[6] * T2.Data[0] + Cst[7] * T2.Data[1] + Cst[8] * T2.Data[2]);

	// A1--B2
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[0] * (Cst[0] * T2.Data[3] + Cst[1] * T2.Data[4] + Cst[2] * T2.Data[5]);
	Cdd += A1.Data[1] * (Cst[3] * T2.Data[3] + Cst[4] * T2.Data[4] + Cst[5] * T2.Data[5]);
	Cdd += A1.Data[2] * (Cst[6] * T2.Data[3] + Cst[7] * T2.Data[4] + Cst[8] * T2.Data[5]);

	// A1--B3
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[0] * (Cst[0] * T2.Data[6] + Cst[1] * T2.Data[7] + Cst[2] * T2.Data[8]);
	Cdd += A1.Data[1] * (Cst[3] * T2.Data[6] + Cst[4] * T2.Data[7] + Cst[5] * T2.Data[8]);
	Cdd += A1.Data[2] * (Cst[6] * T2.Data[6] + Cst[7] * T2.Data[7] + Cst[8] * T2.Data[8]);

	// A1--B4
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[0] * (Cst[0] * T2.Data[9] + Cst[1] * T2.Data[10] + Cst[2] * T2.Data[11]);
	Cdd += A1.Data[1] * (Cst[3] * T2.Data[9] + Cst[4] * T2.Data[10] + Cst[5] * T2.Data[11]);
	Cdd += A1.Data[2] * (Cst[6] * T2.Data[9] + Cst[7] * T2.Data[10] + Cst[8] * T2.Data[11]);

	// A2--B1
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[3] * (Cst[0] * T2.Data[0] + Cst[1] * T2.Data[1] + Cst[2] * T2.Data[2]);
	Cdd += A1.Data[4] * (Cst[3] * T2.Data[0] + Cst[4] * T2.Data[1] + Cst[5] * T2.Data[2]);
	Cdd += A1.Data[5] * (Cst[6] * T2.Data[0] + Cst[7] * T2.Data[1] + Cst[8] * T2.Data[2]);

	// A2--B2
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[3] * (Cst[0] * T2.Data[3] + Cst[1] * T2.Data[4] + Cst[2] * T2.Data[5]);
	Cdd += A1.Data[4] * (Cst[3] * T2.Data[3] + Cst[4] * T2.Data[4] + Cst[5] * T2.Data[5]);
	Cdd += A1.Data[5] * (Cst[6] * T2.Data[3] + Cst[7] * T2.Data[4] + Cst[8] * T2.Data[5]);

	// A2--B3
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[3] * (Cst[0] * T2.Data[6] + Cst[1] * T2.Data[7] + Cst[2] * T2.Data[8]);
	Cdd += A1.Data[4] * (Cst[3] * T2.Data[6] + Cst[4] * T2.Data[7] + Cst[5] * T2.Data[8]);
	Cdd += A1.Data[5] * (Cst[6] * T2.Data[6] + Cst[7] * T2.Data[7] + Cst[8] * T2.Data[8]);

	// A2--B4
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[3] * (Cst[0] * T2.Data[9] + Cst[1] * T2.Data[10] + Cst[2] * T2.Data[11]);
	Cdd += A1.Data[4] * (Cst[3] * T2.Data[9] + Cst[4] * T2.Data[10] + Cst[5] * T2.Data[11]);
	Cdd += A1.Data[5] * (Cst[6] * T2.Data[9] + Cst[7] * T2.Data[10] + Cst[8] * T2.Data[11]);

	// A3--B1
	Cst[0] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[1] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[2] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[4] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[5] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[7] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[8] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[6] * (Cst[0] * T2.Data[0] + Cst[1] * T2.Data[1] + Cst[2] * T2.Data[2]);
	Cdd += A1.Data[7] * (Cst[3] * T2.Data[0] + Cst[4] * T2.Data[1] + Cst[5] * T2.Data[2]);
	Cdd += A1.Data[8] * (Cst[6] * T2.Data[0] + Cst[7] * T2.Data[1] + Cst[8] * T2.Data[2]);

	// A3--B2
	Cst[0] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[1] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[2] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[4] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[5] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[7] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[8] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[6] * (Cst[0] * T2.Data[3] + Cst[1] * T2.Data[4] + Cst[2] * T2.Data[5]);
	Cdd += A1.Data[7] * (Cst[3] * T2.Data[3] + Cst[4] * T2.Data[4] + Cst[5] * T2.Data[5]);
	Cdd += A1.Data[8] * (Cst[6] * T2.Data[3] + Cst[7] * T2.Data[4] + Cst[8] * T2.Data[5]);

	// A3--B3
	Cst[0] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[1] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[2] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[4] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[5] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[7] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[8] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[6] * (Cst[0] * T2.Data[6] + Cst[1] * T2.Data[7] + Cst[2] * T2.Data[8]);
	Cdd += A1.Data[7] * (Cst[3] * T2.Data[6] + Cst[4] * T2.Data[7] + Cst[5] * T2.Data[8]);
	Cdd += A1.Data[8] * (Cst[6] * T2.Data[6] + Cst[7] * T2.Data[7] + Cst[8] * T2.Data[8]);

	// A3--B4
	Cst[0] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[1] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[2] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[4] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[5] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[7] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[8] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cdd += A1.Data[6] * (Cst[0] * T2.Data[9] + Cst[1] * T2.Data[10] + Cst[2] * T2.Data[11]);
	Cdd += A1.Data[7] * (Cst[3] * T2.Data[9] + Cst[4] * T2.Data[10] + Cst[5] * T2.Data[11]);
	Cdd += A1.Data[8] * (Cst[6] * T2.Data[9] + Cst[7] * T2.Data[10] + Cst[8] * T2.Data[11]);

	return Cdd;
}
////////////////////////////////
// torsion-torsion covariance //
////////////////////////////////
double cov_tt(atom At1, atom At2, atom At3, atom At4, atom Bt1, atom Bt2, atom Bt3, atom Bt4, int dimM, double *M)
{
	ADReverse_t Xt1[tCount], Xt2[tCount];
    	ADReverse_t T1, T2;

	Xt1[ 0].Value = At1.x(); Xt1[ 0].ID =  0; Xt1[ 0].Base = true;
	Xt1[ 1].Value = At1.y(); Xt1[ 1].ID =  1; Xt1[ 1].Base = true;
	Xt1[ 2].Value = At1.z(); Xt1[ 2].ID =  2; Xt1[ 2].Base = true;
	Xt1[ 3].Value = At2.x(); Xt1[ 3].ID =  3; Xt1[ 3].Base = true;
	Xt1[ 4].Value = At2.y(); Xt1[ 4].ID =  4; Xt1[ 4].Base = true;
	Xt1[ 5].Value = At2.z(); Xt1[ 5].ID =  5; Xt1[ 5].Base = true;
	Xt1[ 6].Value = At3.x(); Xt1[ 6].ID =  6; Xt1[ 6].Base = true;
	Xt1[ 7].Value = At3.y(); Xt1[ 7].ID =  7; Xt1[ 7].Base = true;
	Xt1[ 8].Value = At3.z(); Xt1[ 8].ID =  8; Xt1[ 8].Base = true;
	Xt1[ 9].Value = At4.x(); Xt1[ 9].ID =  9; Xt1[ 9].Base = true;
	Xt1[10].Value = At4.y(); Xt1[10].ID = 10; Xt1[10].Base = true;
	Xt1[11].Value = At4.z(); Xt1[11].ID = 11; Xt1[11].Base = true;
    	T1 = Algorithm_t<ADReverse_t>(Xt1);

	Xt2[ 0].Value = Bt1.x(); Xt2[ 0].ID =  0; Xt2[ 0].Base = true;
	Xt2[ 1].Value = Bt1.y(); Xt2[ 1].ID =  1; Xt2[ 1].Base = true;
	Xt2[ 2].Value = Bt1.z(); Xt2[ 2].ID =  2; Xt2[ 2].Base = true;
	Xt2[ 3].Value = Bt2.x(); Xt2[ 3].ID =  3; Xt2[ 3].Base = true;
	Xt2[ 4].Value = Bt2.y(); Xt2[ 4].ID =  4; Xt2[ 4].Base = true;
	Xt2[ 5].Value = Bt2.z(); Xt2[ 5].ID =  5; Xt2[ 5].Base = true;
	Xt2[ 6].Value = Bt3.x(); Xt2[ 6].ID =  6; Xt2[ 6].Base = true;
	Xt2[ 7].Value = Bt3.y(); Xt2[ 7].ID =  7; Xt2[ 7].Base = true;
	Xt2[ 8].Value = Bt3.z(); Xt2[ 8].ID =  8; Xt2[ 8].Base = true;
	Xt2[ 9].Value = Bt4.x(); Xt2[ 9].ID =  9; Xt2[ 9].Base = true;
	Xt2[10].Value = Bt4.y(); Xt2[10].ID = 10; Xt2[10].Base = true;
	Xt2[11].Value = Bt4.z(); Xt2[11].ID = 11; Xt2[11].Base = true;
    	T2 = Algorithm_t<ADReverse_t>(Xt2);

	double Cdd = 0.0;
	double Cst[9];

	// A1--B1
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cdd += T1.Data[0] * (Cst[0] * T2.Data[0] + Cst[1] * T2.Data[1] + Cst[2] * T2.Data[2]);
	Cdd += T1.Data[1] * (Cst[3] * T2.Data[0] + Cst[4] * T2.Data[1] + Cst[5] * T2.Data[2]);
	Cdd += T1.Data[2] * (Cst[6] * T2.Data[0] + Cst[7] * T2.Data[1] + Cst[8] * T2.Data[2]);

	// A1--B2
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cdd += T1.Data[0] * (Cst[0] * T2.Data[3] + Cst[1] * T2.Data[4] + Cst[2] * T2.Data[5]);
	Cdd += T1.Data[1] * (Cst[3] * T2.Data[3] + Cst[4] * T2.Data[4] + Cst[5] * T2.Data[5]);
	Cdd += T1.Data[2] * (Cst[6] * T2.Data[3] + Cst[7] * T2.Data[4] + Cst[8] * T2.Data[5]);

	// A1--B3
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cdd += T1.Data[0] * (Cst[0] * T2.Data[6] + Cst[1] * T2.Data[7] + Cst[2] * T2.Data[8]);
	Cdd += T1.Data[1] * (Cst[3] * T2.Data[6] + Cst[4] * T2.Data[7] + Cst[5] * T2.Data[8]);
	Cdd += T1.Data[2] * (Cst[6] * T2.Data[6] + Cst[7] * T2.Data[7] + Cst[8] * T2.Data[8]);

	// A1--B4
	Cst[0] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[1] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[2] = M[((At1.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[4] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[5] = M[((At1.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[7] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[8] = M[((At1.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cdd += T1.Data[0] * (Cst[0] * T2.Data[9] + Cst[1] * T2.Data[10] + Cst[2] * T2.Data[11]);
	Cdd += T1.Data[1] * (Cst[3] * T2.Data[9] + Cst[4] * T2.Data[10] + Cst[5] * T2.Data[11]);
	Cdd += T1.Data[2] * (Cst[6] * T2.Data[9] + Cst[7] * T2.Data[10] + Cst[8] * T2.Data[11]);

	// A2--B1
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cdd += T1.Data[3] * (Cst[0] * T2.Data[0] + Cst[1] * T2.Data[1] + Cst[2] * T2.Data[2]);
	Cdd += T1.Data[4] * (Cst[3] * T2.Data[0] + Cst[4] * T2.Data[1] + Cst[5] * T2.Data[2]);
	Cdd += T1.Data[5] * (Cst[6] * T2.Data[0] + Cst[7] * T2.Data[1] + Cst[8] * T2.Data[2]);

	// A2--B2
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cdd += T1.Data[3] * (Cst[0] * T2.Data[3] + Cst[1] * T2.Data[4] + Cst[2] * T2.Data[5]);
	Cdd += T1.Data[4] * (Cst[3] * T2.Data[3] + Cst[4] * T2.Data[4] + Cst[5] * T2.Data[5]);
	Cdd += T1.Data[5] * (Cst[6] * T2.Data[3] + Cst[7] * T2.Data[4] + Cst[8] * T2.Data[5]);

	// A2--B3
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cdd += T1.Data[3] * (Cst[0] * T2.Data[6] + Cst[1] * T2.Data[7] + Cst[2] * T2.Data[8]);
	Cdd += T1.Data[4] * (Cst[3] * T2.Data[6] + Cst[4] * T2.Data[7] + Cst[5] * T2.Data[8]);
	Cdd += T1.Data[5] * (Cst[6] * T2.Data[6] + Cst[7] * T2.Data[7] + Cst[8] * T2.Data[8]);

	// A2--B4
	Cst[0] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[1] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[2] = M[((At2.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[4] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[5] = M[((At2.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[7] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[8] = M[((At2.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cdd += T1.Data[3] * (Cst[0] * T2.Data[9] + Cst[1] * T2.Data[10] + Cst[2] * T2.Data[11]);
	Cdd += T1.Data[4] * (Cst[3] * T2.Data[9] + Cst[4] * T2.Data[10] + Cst[5] * T2.Data[11]);
	Cdd += T1.Data[5] * (Cst[6] * T2.Data[9] + Cst[7] * T2.Data[10] + Cst[8] * T2.Data[11]);

	// A3--B1
	Cst[0] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[1] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[2] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[4] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[5] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[7] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[8] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cdd += T1.Data[6] * (Cst[0] * T2.Data[0] + Cst[1] * T2.Data[1] + Cst[2] * T2.Data[2]);
	Cdd += T1.Data[7] * (Cst[3] * T2.Data[0] + Cst[4] * T2.Data[1] + Cst[5] * T2.Data[2]);
	Cdd += T1.Data[8] * (Cst[6] * T2.Data[0] + Cst[7] * T2.Data[1] + Cst[8] * T2.Data[2]);

	// A3--B2
	Cst[0] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[1] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[2] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[4] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[5] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[7] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[8] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cdd += T1.Data[6] * (Cst[0] * T2.Data[3] + Cst[1] * T2.Data[4] + Cst[2] * T2.Data[5]);
	Cdd += T1.Data[7] * (Cst[3] * T2.Data[3] + Cst[4] * T2.Data[4] + Cst[5] * T2.Data[5]);
	Cdd += T1.Data[8] * (Cst[6] * T2.Data[3] + Cst[7] * T2.Data[4] + Cst[8] * T2.Data[5]);

	// A3--B3
	Cst[0] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[1] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[2] = M[((At3.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[4] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[5] = M[((At3.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[7] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[8] = M[((At3.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cdd += T1.Data[6] * (Cst[0] * T2.Data[6] + Cst[1] * T2.Data[7] + Cst[2] * T2.Data[8]);
	Cdd += T1.Data[7] * (Cst[3] * T2.Data[6] + Cst[4] * T2.Data[7] + Cst[5] * T2.Data[8]);
	Cdd += T1.Data[8] * (Cst[6] * T2.Data[6] + Cst[7] * T2.Data[7] + Cst[8] * T2.Data[8]);

	// A3--B4
	Cst[0] = M[((At4.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[1] = M[((At4.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[2] = M[((At4.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At4.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[4] = M[((At4.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[5] = M[((At4.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At4.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[7] = M[((At4.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[8] = M[((At4.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cdd += T1.Data[6] * (Cst[0] * T2.Data[9] + Cst[1] * T2.Data[10] + Cst[2] * T2.Data[11]);
	Cdd += T1.Data[7] * (Cst[3] * T2.Data[9] + Cst[4] * T2.Data[10] + Cst[5] * T2.Data[11]);
	Cdd += T1.Data[8] * (Cst[6] * T2.Data[9] + Cst[7] * T2.Data[10] + Cst[8] * T2.Data[11]);

	// A4--B1
	Cst[0] = M[((At4.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[1] = M[((At4.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[2] = M[((At4.id0 - 1) * 3 + 0) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At4.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[4] = M[((At4.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[5] = M[((At4.id0 - 1) * 3 + 1) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At4.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 0)]; Cst[7] = M[((At4.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 1)]; Cst[8] = M[((At4.id0 - 1) * 3 + 2) * dimM + ((Bt1.id0 - 1) * 3 + 2)];
	Cdd += T1.Data[ 9] * (Cst[0] * T2.Data[0] + Cst[1] * T2.Data[1] + Cst[2] * T2.Data[2]);
	Cdd += T1.Data[10] * (Cst[3] * T2.Data[0] + Cst[4] * T2.Data[1] + Cst[5] * T2.Data[2]);
	Cdd += T1.Data[11] * (Cst[6] * T2.Data[0] + Cst[7] * T2.Data[1] + Cst[8] * T2.Data[2]);

	// A4--B2
	Cst[0] = M[((At4.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[1] = M[((At4.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[2] = M[((At4.id0 - 1) * 3 + 0) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At4.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[4] = M[((At4.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[5] = M[((At4.id0 - 1) * 3 + 1) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At4.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 0)]; Cst[7] = M[((At4.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 1)]; Cst[8] = M[((At4.id0 - 1) * 3 + 2) * dimM + ((Bt2.id0 - 1) * 3 + 2)];
	Cdd += T1.Data[ 9] * (Cst[0] * T2.Data[3] + Cst[1] * T2.Data[4] + Cst[2] * T2.Data[5]);
	Cdd += T1.Data[10] * (Cst[3] * T2.Data[3] + Cst[4] * T2.Data[4] + Cst[5] * T2.Data[5]);
	Cdd += T1.Data[11] * (Cst[6] * T2.Data[3] + Cst[7] * T2.Data[4] + Cst[8] * T2.Data[5]);

	// A4--B3
	Cst[0] = M[((At4.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[1] = M[((At4.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[2] = M[((At4.id0 - 1) * 3 + 0) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At4.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[4] = M[((At4.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[5] = M[((At4.id0 - 1) * 3 + 1) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At4.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 0)]; Cst[7] = M[((At4.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 1)]; Cst[8] = M[((At4.id0 - 1) * 3 + 2) * dimM + ((Bt3.id0 - 1) * 3 + 2)];
	Cdd += T1.Data[ 9] * (Cst[0] * T2.Data[6] + Cst[1] * T2.Data[7] + Cst[2] * T2.Data[8]);
	Cdd += T1.Data[10] * (Cst[3] * T2.Data[6] + Cst[4] * T2.Data[7] + Cst[5] * T2.Data[8]);
	Cdd += T1.Data[11] * (Cst[6] * T2.Data[6] + Cst[7] * T2.Data[7] + Cst[8] * T2.Data[8]);

	// A4--B4
	Cst[0] = M[((At4.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[1] = M[((At4.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[2] = M[((At4.id0 - 1) * 3 + 0) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[3] = M[((At4.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[4] = M[((At4.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[5] = M[((At4.id0 - 1) * 3 + 1) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cst[6] = M[((At4.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 0)]; Cst[7] = M[((At4.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 1)]; Cst[8] = M[((At4.id0 - 1) * 3 + 2) * dimM + ((Bt4.id0 - 1) * 3 + 2)];
	Cdd += T1.Data[ 9] * (Cst[0] * T2.Data[9] + Cst[1] * T2.Data[10] + Cst[2] * T2.Data[11]);
	Cdd += T1.Data[10] * (Cst[3] * T2.Data[9] + Cst[4] * T2.Data[10] + Cst[5] * T2.Data[11]);
	Cdd += T1.Data[11] * (Cst[6] * T2.Data[9] + Cst[7] * T2.Data[10] + Cst[8] * T2.Data[11]);

	return Cdd;
}
/////////////////////////////////////////////////////////////////////////////////
// The conversion routine                                                      //
// nA  = number of atoms                                                       //
// mol = molecule bearing both information on Z-Matrix and reference structure //
// M   = variance-covariance matrix in Cartesian coordinates 3nA x 3nA         //
// C   = variance-covariance matrix in internal coordinates (3nA-6) x (3nA-6)  //
/////////////////////////////////////////////////////////////////////////////////
void varcovar(int nA, molecule *mol, double *M, double *C)
{
	int qrow, qcol, nq, nx;
	atom a1, a2;
	atom a1A, a1B, a1C, a2A, a2B, a2C;
	zmatline z1, z2;

	nx = 3 * nA;
	nq = 3 * nA - 6;

	qrow = 0;
	for (int ia = 2; ia <= nA; ++ia)
	{
		a1 = mol->getAtom(ia);
		z1 = a1.getZMatrixEntry();

		qcol = qrow;

		for (int ja = ia; ja <= nA; ++ja)
		{
			a2 = mol->getAtom(ja);
			z2 = a2.getZMatrixEntry();

			if (z1.A)
			{
				/* cov(d1-d2) */
				if (z2.A)
				{
					a1A = mol->getAtom(z1.A);
					a2A = mol->getAtom(z2.A);
					C[qrow * nq + qcol] = cov_dd(a1, a1A, a2, a2A, nx, M);
					C[qcol * nq + qrow] = C[qrow * nq + qcol];
				}
				/* cov(d1-a2) */
				if (z2.B)
				{
					a1A = mol->getAtom(z1.A);
					a2A = mol->getAtom(z2.A);
					a2B = mol->getAtom(z2.B);
					C[qrow * nq + (qcol + 1)] = cov_da(a1, a1A, a2, a2A, a2B, nx, M);
					C[(qcol + 1) * nq + qrow] = C[qrow * nq + (qcol + 1)];
				}
				/* cov(d1-t2) */
				if (z2.C)
				{
					a1A = mol->getAtom(z1.A);
					a2A = mol->getAtom(z2.A);
					a2B = mol->getAtom(z2.B);
					a2C = mol->getAtom(z2.C);
					C[qrow * nq + (qcol + 2)] = cov_dt(a1, a1A, a2, a2A, a2B, a2C, nx, M);
					C[(qcol + 2) * nq + qrow] = C[qrow * nq + (qcol + 2)];
				}
				qrow++;
			}
			if (z1.B)
			{
				/* cov(a1-d2) */
				if (z2.A)
				{
					a1A = mol->getAtom(z1.A);
					a1B = mol->getAtom(z1.B);
					a2A = mol->getAtom(z2.A);
					C[qrow * nq + qcol] = cov_da(a2, a2A, a1, a1A, a1B, nx, M);
					C[qcol * nq + qrow] = C[qrow * nq + qcol];
				}
				/* cov(a1-a2) */
				if (z2.B)
				{
					a1A = mol->getAtom(z1.A);
					a1B = mol->getAtom(z1.B);
					a2A = mol->getAtom(z2.A);
					a2B = mol->getAtom(z2.B);
					C[qrow * nq + (qcol + 1)] = cov_aa(a1, a1A, a1B, a2, a2A, a2B, nx, M);
					C[(qcol + 1) * nq + qrow] = C[qrow * nq + (qcol + 1)];
				}
				/* cov(a1-t2) */
				if (z2.C)
				{
					a1A = mol->getAtom(z1.A);
					a1B = mol->getAtom(z1.B);
					a2A = mol->getAtom(z2.A);
					a2B = mol->getAtom(z2.B);
					a2C = mol->getAtom(z2.C);
					C[qrow * nq + (qcol + 2)] = cov_at(a1, a1A, a1B, a2, a2A, a2B, a2C, nx, M);
					C[(qcol + 2) * nq + qrow] = C[qrow * nq + (qcol + 2)];
				}
				qrow++;
			}
			if (z1.C)
			{
				/* cov(t1-d2) */
				if (z2.A)
				{
					a1A = mol->getAtom(z1.A);
					a1B = mol->getAtom(z1.B);
					a1C = mol->getAtom(z1.C);
					a2A = mol->getAtom(z2.A);
					C[qrow * nq + qcol] = cov_dt(a2, a2A, a1, a1A, a1B, a1C, nx, M);
					C[qcol * nq + qrow] = C[qrow * nq + qcol];
				}
				/* cov(t1-a2) */
				if (z2.B)
				{
					a1A = mol->getAtom(z1.A);
					a1B = mol->getAtom(z1.B);
					a1C = mol->getAtom(z1.C);
					a2A = mol->getAtom(z2.A);
					a2B = mol->getAtom(z2.B);
					C[qrow * nq + (qcol + 1)] = cov_at(a2, a2A, a2B, a1, a1A, a1B, a1C, nx, M);
					C[(qcol + 1) * nq + qrow] = C[qrow * nq + (qcol + 1)];
				}
				/* cov(t1-t2) */
				if (z2.C)
				{
					a1A = mol->getAtom(z1.A);
					a1B = mol->getAtom(z1.B);
					a1C = mol->getAtom(z1.C);
					a2A = mol->getAtom(z2.A);
					a2B = mol->getAtom(z2.B);
					a2C = mol->getAtom(z2.C);
					C[qrow * nq + (qcol + 2)] = cov_tt(a1, a1A, a1B, a1C, a2, a2A, a2B, a2C, nx, M);
					C[(qcol + 2) * nq + qrow] = C[(qcol + 2) * nq + qcol];
				}
				qrow++;
			}
			qcol += (z2.C ? 3 : (z2.B ? 2 : z2.A ? 1 : 0));
			qrow -= (z1.C ? 3 : (z1.B ? 2 : z1.A ? 1 : 0));
		}
		qrow += (z1.C ? 3 : (z1.B ? 2 : z1.A ? 1 : 0));
	}

/*
	for (qrow = 0; qrow < nq; ++qrow)
	{
		for (qcol = qrow + 1; qcol < nq; ++qcol)
			C[qcol * nq + qrow] = 0.0;//C[qrow * nq + qcol];
	}
*/

	return;
}

int main(void)
{
	// Configure SeStO
	string str;
	char* zmtlibHome;
	zmtlibHome = getenv ("ZMATLIB_HOME");
	if (zmtlibHome == NULL)
	{
		cout << endl << "ERROR: the ZMATLIB_HOME envinronment variable is not set. Please set the path or use the script salem" << endl << endl;
		return 1;
	}
	config conf;
	str.assign(zmtlibHome); str.append("/sesto_1.0/sesto");
	conf.setSesto(str);
	str.assign(zmtlibHome); str.append("/data/vdw.dat");
	conf.setVdwFile(str);
	str.assign(zmtlibHome); str.append("/data/aminoacids.dat");
	conf.setAminoAcidsFile(str);
	str.assign(zmtlibHome); str.append("/data/atoms_masses.dat");
	conf.setMassFile(str);
	
	// Load the molecule
	molecule mol(&conf);
	mol.setIOFileName("H2O_1.pdb");
	mol.loadMolecule();
	int natoms = mol.getNAtoms();
	int nx = 3 * natoms;
	int nq = 3 * natoms - 6;

	zmatline *z = new zmatline[3];
	z[0].A = 0; z[0].B = 0; z[0].C = 0; z[0].d = 0.0; z[0].theta = 0.0;    z[0].phi = 0.0;
	z[1].A = 1; z[1].B = 0; z[1].C = 0; z[1].d = 1.0; z[1].theta = 0.0;    z[1].phi = 0.0;
	z[2].A = 1; z[2].B = 2; z[2].C = 0; z[2].d = 1.0; z[2].theta = 1.7627; z[2].phi = 0.0;
	mol.setZMatrix(z);
	
	double *M = new double[nx * nx];
	
	fstream f;
	f.open("mcovar.dat", ios::in);
	for (int i = 0; i < nx; ++i)
	{
		for (int j = 0; j < nx; ++j)
			f >> M[i * nx + j];
	}
	f.close();

	double *C = new double [nq * nq];

	varcovar(natoms, &mol, M, C);

	f.open("mcovar_q.dat", ios::out);
	for (int i = 0; i < nq; ++i)
	{
		for (int j = 0; j < nq; ++j)
			f << C[i * nq + j] << " ";
		f << endl;
	}
	f.close();
	
}
