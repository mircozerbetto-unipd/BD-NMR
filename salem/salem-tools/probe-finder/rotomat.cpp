#include "rotomat.h"

vectorOfDoubles rotomat(atom a1, atom a2, atom a3)
{
	double norm;
	vectorOfDoubles rA1 = a1.getPosition0();
	vectorOfDoubles rA2 = a2.getPosition0();
	vectorOfDoubles rA3 = a3.getPosition0();

	vectorOfDoubles R(9);

	// x
	R[0] = rA2[0] - rA1[0];
	R[1] = rA2[1] - rA1[1];
	R[2] = rA2[2] - rA1[2];
	norm = 1.0 / sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
	R[0] *= norm;
	R[1] *= norm;
	R[2] *= norm;

	// supporting v
	R[3] = rA3[0] - rA2[0];
	R[4] = rA3[1] - rA2[1];
	R[5] = rA3[2] - rA2[2];
	norm = 1.0 / sqrt(R[3] * R[3] + R[4] * R[4] + R[5] * R[5]);
	R[3] *= norm;
	R[4] *= norm;
	R[5] *= norm;

	// z
	R[6] = R[1] * R[5] - R[2] * R[4];
	R[7] = R[2] * R[3] - R[0] * R[5];
	R[8] = R[0] * R[4] - R[1] * R[3];
	norm = 1.0 / sqrt(R[6] * R[6] + R[7] * R[7] + R[8] * R[8]);
	R[6] *= norm;
	R[7] *= norm;
	R[8] *= norm;

	// y
	R[3] = R[7] * R[2] - R[8] * R[1];
	R[4] = R[8] * R[0] - R[6] * R[2];
	R[5] = R[6] * R[1] - R[7] * R[0];
	norm = 1.0 / sqrt(R[3] * R[3] + R[4] * R[4] + R[5] * R[5]);
	R[3] *= norm;
	R[4] *= norm;
	R[5] *= norm;
	
	return R;
}

