#include "orderParameter.h"

inline int qidx(int qid)
{
	switch(qid)
	{
		case 3:
			return 0;
		case 6:
			return 1;
		case 7:
			return 2;
		default:
			return (qid - 6);
	}
}

double S20(int Hid, int Nid, int Cid, molecule* mol, double *invsmallu, int uDim, int nq)
{
	d2k0Grad W;

	W.setMolecule(mol);

	W.setHid(Hid);
	W.setNid(Nid);
	W.setCid(Cid);

	W.setK(0);

	W.setBetaDC(0.0); // Tilt not needed to calculate N-H order parameters

	W.buildTaylorExpansionCoefficients();

	dcomplex czero = std::complex<double>(0.0, 0.0);
	dcomplex D200;

	int dim = W.getArraysMainDim();
	int *ig = new int[dim];
	dcomplex *g1 = new dcomplex[dim];
	dcomplex *g2 = new dcomplex[dim * dim];

	W.getQIdx(ig);

	W.getZeroOrder(DIP, &D200);
	W.getFirstDerivatives(DIP, g1);
	W.getSecondDerivatives(DIP, g2);

	double S = D200.real();

	std::cout << " - Smax = " << S << std::endl;/////////////////////////////////////////////////////////////////////////////////////

	dcomplex tmp;
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			tmp = czero;
			for (int m = 0; m < nq; ++m)
			{
				tmp.real(tmp.real() + invsmallu[qidx(ig[i]) * uDim + m] * invsmallu[qidx(ig[j]) * uDim + m]);
			}
			tmp *= g2[i * dim + j] * 0.5;
			S += tmp.real();
		}
	}
	return S;
}
