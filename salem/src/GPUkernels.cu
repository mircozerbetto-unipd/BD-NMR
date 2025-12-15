// Zero out k
__global__ void gpu_zero_smallk(int nk, double *smallk)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < nk && c < nk)
		smallk[r * nk + c] = 0.0;
	return;
}

// Copy eigenvectors of K into k
__global__ void gpu_copy_K(int nq, int nk, double *smallk, double *K)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < nq && c < nq)
	{
		// TO DEL // if (r <= c)
			smallk[r * nk + c] = K[r * nq + c];
		// TO DEL // else
		// TO DEL // 	smallk[r * nk + c] = K[c * nq + r];
	}
	return;
}

// Copy eigenvectors of k_LB into k
// nq: dimension of K (3N - 6)
// nkLB: dimension of k_LB (3N-3)
// nk: dimension of k (6N-9)
__global__ void gpu_copy_g(int nq, int nkLB, int nk, double *smallk, double *kLB)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < nkLB && c < nkLB)
	{
		// TO DEL //if (r <= c)
		// TO DEL //	smallk[(r + 3 + shift) * nk + (c + 3 + shift)] = G[r * nq + c];
		// TO DEL //else
		// TO DEL //	smallk[(r + 3 + shift) * nk + (c + 3 + shift)] = G[c * nq + r];
		smallk[(r + nq) * nk + (c + nq)] = kLB[r * nkLB + c];
	}
	return;
}

// Copy gA and (A^tr g) into k
__global__ void gpu_copy_gA(int shift, int nq, int nk, double *smallk, double *GA)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < nq && c < 3)
	{
		smallk[(r + 3 + shift) * nk + (c + shift)] = -GA[r * 3 + c];
		smallk[(c + shift) * nk + (r + 3 + shift)] = -GA[r * 3 + c];
	}
	return;
}

// Copy [inv(I) + A^tr gA] into k
__global__ void gpu_copy_invIAgA(int shift, int nk, double *smallk, double *invI, double *AgA)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < 3 && c < 3)
	{
		if (r <= c)
			smallk[(r + shift) * nk + (c + shift)] = invI[r * 3 + c] + AgA[r * 3 + c];
		else
			smallk[(r + shift) * nk + (c + shift)] = invI[c * 3 + r] + AgA[c * 3 + r];
	}
	return;
}

// Copy kLB into k
// shift: location of the sub-block of smallk where to put smallkLB
// nkLB: dimension of smallkLB
// nk: dimension of smallk
__global__ void gpu_copy_smallkLB(int shift, int nkLB, int nk, double *smallk, double *smallkLB)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < nkLB && c < nkLB)
		smallk[(r + shift) * nk + (c + shift)] = smallkLB[r * nkLB + c];
	return;
}

// Sqrt of the elements of an array (NB: destroyes the original array)
__global__ void gpu_sqrt_v(int dim, double *v)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < dim)
		v[i] = sqrt(v[i]);
	return;
}

// Sqrt of the inverse of the elements of an array (NB: destroyes the original array)
__global__ void gpu_sqrt_inv_v(int dim, double *v)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < dim)
		v[i] = sqrt(1.0 / v[i]);
	return;
}

// Mult A * B where: A = diagonal elements of a diagonal matrix; B = matrix (square)
__global__ void gpu_diagA_times_B(int dim, double mult, double *A, double *B, double *C)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < dim && c < dim)
		C[r * dim + c] = mult * (A[r] * B[r * dim + c]);

	return;
}

// Mult A * B^tr where: A = diagonal elements of a diagonal matrix; B = matrix (square)
__global__ void gpu_diagA_times_Btr(int dim, double mult, double *A, double *B, double *C)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < dim && c < dim)
		C[r * dim + c] = mult * (A[r] * B[c * dim + r]);
	return;
}

// Mult B^tr * A where: A = diagonal elements of a diagonal matrix; B = matrix (square)
__global__ void gpu_Btr_times_diagA(int dim, double mult, double *B, double *A, double *C)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < dim && c < dim)
		C[r * dim + c] = mult * (B[c * dim + r] * A[c]);
	return;
}

// Build the B matrix (reduced inertia tensor)
__global__ void gpu_build_Bmat_tr(int nq, int nk, int nf, double kBT, double *friction, double *Bmat)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < nk && c < nk)
	{
		if (r >= nq && c >= nq)
			Bmat[r * nk + c] = kBT * friction[(r - nq) * nf + (c - nq)];
		else if ((r < nq && c < (3 + nq)) || (r < (3 + nq) && c < nq))
			Bmat[r * nk + c] = 0.0;
		else if (r < nq && c >= (3 + nq))
			Bmat[r * nk + c] = (r == (c - 3 - nq)) ?  kBT : 0.0;
		else if (r >= (3 + nq) && c < nq)
			Bmat[r * nk + c] = ((r - 3 - nq) == c) ? -kBT : 0.0;
		else
			; // do nothing (anyway, there should be no excluded cases)
	}
	return;
}

// Method to buld w = s + wX
//
// WARNING: wx is stored by row, while wxs will be stored column-wise to meet
//          the F77-like implementation of LAPACK libratries in MAGMA
__global__ void gpu_make_wxs(int nk, double s, double *wx, magmaDoubleComplex *wxs)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < nk && c < nk)
		wxs[r + c * nk] = r == c ? MAGMA_Z_MAKE(wx[r * nk + c], s) : MAGMA_Z_MAKE(wx[r * nk + c], 0.0);
	return;
}

// Store wint column-wise
__global__ void gpu_makeColumnWise(int nk, double *dWintTmp, magmaDoubleComplex *dWint)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < 3 && c < nk)
		dWint[r + c * 3] = MAGMA_Z_MAKE(dWintTmp[r * nk + c], 0.0);
	return;
}

__global__ void gpu_makeMbg(int nk, int b, int g, magmaDoubleComplex *wi, magmaDoubleComplex *M)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < nk && c < nk)
	{
		//M[r + c * nk] = MAGMA_Z_ADD(wi[g + r * 3], wi[g + c * 3]);
		//M[r + c * nk] = MAGMA_Z_MUL(M[r + c * nk], wi[b + r * 3]);
		//if (r == c)
		//	M[r + c * nk] = MAGMA_Z_ADD(wi[g + r * 3], wi[g + c * 3]);
		//else
			M[r + c * nk] = wi[g + c * 3];
		M[r + c * nk] = MAGMA_Z_MUL(M[r + c * nk], wi[b + r * 3]);
	}
	return;
}

#define CP(j,k) (sqrt((double)(j * (j + 1)) - (double)(k * (k + 1))))
#define CM(j,k) (sqrt((double)(j * (j + 1)) - (double)(k * (k - 1))))
// M matrices are stored in column-wise
__global__ void gpu_buildMmatrices(int JJ, magmaDoubleComplex *Mx, magmaDoubleComplex *My, magmaDoubleComplex *Mz)
{
	int k1 = threadIdx.x + blockIdx.x * blockDim.x - JJ;
	int k2 = threadIdx.y + blockIdx.y * blockDim.y - JJ;
	int dimM = 2 * JJ + 1;
	if ((k1 + JJ) < dimM && (k2 + JJ) < dimM)
	{
		Mx[(k1 + JJ) + (k2 + JJ) * dimM] = MAGMA_Z_ZERO;
		My[(k1 + JJ) + (k2 + JJ) * dimM] = MAGMA_Z_ZERO;
		Mz[(k1 + JJ) + (k2 + JJ) * dimM] = MAGMA_Z_ZERO;
		switch (k1 - k2)
		{
			case 0:
			{
				Mz[(k1 + JJ) + (k2 + JJ) * dimM] = MAGMA_Z_MAKE(0.0, (double)k1);
				break;
			}
			case -1:
			{
				double coeff = CM(JJ, k2);
				Mx[(k1 + JJ) + (k2 + JJ) * dimM] = MAGMA_Z_MAKE(0.0, 0.5 * coeff);
				My[(k1 + JJ) + (k2 + JJ) * dimM] = MAGMA_Z_MAKE(-0.5 * coeff, 0.0);
				break;
			}
			case 1:
			{
				double coeff = CP(JJ, k2);
				Mx[(k1 + JJ) + (k2 + JJ) * dimM] = MAGMA_Z_MAKE(0.0, 0.5 * coeff);
				My[(k1 + JJ) + (k2 + JJ) * dimM] = MAGMA_Z_MAKE(0.5 * coeff, 0.0);
				break;
			}
		}
	}
	return;
}

__global__ void gpu_zeroOutArray(int n, magmaDoubleComplex *a)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < n)
		a[i] = MAGMA_Z_ZERO;
	return;
}

// Make a matrix of inverse of Rs eigenvalues
__global__ void gpu_invRsEig(int dimM, magmaDoubleComplex *eig, magmaDoubleComplex *invEig)
{
	int k1 = threadIdx.x + blockIdx.x * blockDim.x;
	int k2 = threadIdx.y + blockIdx.y * blockDim.y;
	if (k1 < dimM && k2 < dimM)
	{
		if (k1 == k2)
			invEig[k1 + k2 * dimM] = MAGMA_Z_DIV(MAGMA_Z_ONE, eig[k1]);
		else
			invEig[k1 + k2 * dimM] = MAGMA_Z_ZERO;
	}
	return;
}

// Compute Rs
__global__ void gpu_makeRs(int dimM, magmaDoubleComplex s, magmaDoubleComplex *R2, magmaDoubleComplex *R4, magmaDoubleComplex *Rs, int srb)
{
	int k1 = threadIdx.x + blockIdx.x * blockDim.x;
	int k2 = threadIdx.y + blockIdx.y * blockDim.y;
	if (k1 < dimM && k2 < dimM)
	{
		if (srb)
		{
			Rs[k1 + k2 * dimM] = MAGMA_Z_ADD(R2[k1 + k2 * dimM], R4[k1 + k2 * dimM]);
			// Pade approximant
//			magmaDoubleComplex mR4 = MAGMA_Z_MAKE(-MAGMA_Z_REAL(R4[k1 + k2 * dimM]), -MAGMA_Z_IMAG(R4[k1 + k2 * dimM]));
//			Rs[k1 + k2 * dimM] = MAGMA_Z_MUL(R2[k1 + k2 * dimM], R2[k1 + k2 * dimM]);
//			Rs[k1 + k2 * dimM] = MAGMA_Z_DIV(Rs[k1 + k2 * dimM], MAGMA_Z_ADD(R2[k1 + k2 * dimM], mR4));
		}
		else
		{
			Rs[k1 + k2 * dimM] = R2[k1 + k2 * dimM];
		}

		if (k1 == k2)
			Rs[k1 + k2 * dimM] = MAGMA_Z_ADD(Rs[k1 + k2 * dimM], s);
	}
	return;
}

// Branch 2 methods

__global__ void gpu_buildGammaMatrix(int nb, int nk, int nKK, magmaDoubleComplex *Mx, magmaDoubleComplex *My, magmaDoubleComplex *Mz, magmaDoubleComplex *wx, magmaDoubleComplex *wi, magmaDoubleComplex *gamma)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < nb && c < nb)
	{
		if (r < nKK && c < nKK)
			gamma[r + c * nb] = MAGMA_Z_ZERO;
		else if (r < nKK && c >= nKK)
		{
			int Kr = r;
			int sc = c / nKK;
			int Kc = c - nKK * sc;
			sc--;
			gamma[r + c * nb] = MAGMA_Z_ZERO;
			gamma[r + c * nb] = MAGMA_Z_ADD(gamma[r + c * nb], MAGMA_Z_MUL(wi[0 + sc * 3], Mx[Kr + Kc * nKK]));
			gamma[r + c * nb] = MAGMA_Z_ADD(gamma[r + c * nb], MAGMA_Z_MUL(wi[1 + sc * 3], My[Kr + Kc * nKK]));
			gamma[r + c * nb] = MAGMA_Z_ADD(gamma[r + c * nb], MAGMA_Z_MUL(wi[2 + sc * 3], Mz[Kr + Kc * nKK]));
		}
		else if (r >= nKK && c < nKK)
		{
			int sr = r / nKK;
			int Kr = r - nKK * sr;
			sr--;
			int Kc = c;
			gamma[r + c * nb] = MAGMA_Z_ZERO;
			gamma[r + c * nb] = MAGMA_Z_ADD(gamma[r + c * nb], MAGMA_Z_MUL(wi[0 + sr * 3], Mx[Kr + Kc * nKK]));
			gamma[r + c * nb] = MAGMA_Z_ADD(gamma[r + c * nb], MAGMA_Z_MUL(wi[1 + sr * 3], My[Kr + Kc * nKK]));
			gamma[r + c * nb] = MAGMA_Z_ADD(gamma[r + c * nb], MAGMA_Z_MUL(wi[2 + sr * 3], Mz[Kr + Kc * nKK]));
		}
		else
		{
			int sr = r / nKK;
			int Kr = r - nKK * sr;
			sr--;
			int sc = c / nKK;
			int Kc = c - nKK * sc;
			sc--;
			if (Kr == Kc)
				gamma[r + c * nb] = wx[sr + sc * nk];
			else
				gamma[r + c * nb] = MAGMA_Z_ZERO;
		}
	}
	return;
}

__global__ void gpu_ShiftGamma(int nb, double s, magmaDoubleComplex *gamma, magmaDoubleComplex *gammaS)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < nb && c < nb)
		gammaS[r + c * nb] = (r == c) ? MAGMA_Z_ADD(MAGMA_Z_MAKE(0.0, s), gamma[r + c * nb]) : gamma[r + c * nb];
	return;
}


__global__ void gpu_buildStartingVector(int nb, int nKK, magmaDoubleComplex *stvec)
{
		int r = threadIdx.x + blockIdx.x * blockDim.x;
		int c = threadIdx.y + blockIdx.y * blockDim.y;
		if (r < nb && c < nKK)
			stvec[r + c * nb] = (r == c) ? MAGMA_Z_ONE : MAGMA_Z_ZERO;
		return;
}

// Build F1(s)
__global__ void gpu_build_F1(int nk, magmaDoubleComplex s, magmaDoubleComplex *w, magmaDoubleComplex *F1)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	if (r < nk)
		F1[r] = MAGMA_Z_DIV(MAGMA_Z_ONE, MAGMA_Z_ADD(s, w[r]));
	return;
}

// Build F2(s)
__global__ void gpu_build_F2(int nk, magmaDoubleComplex s, magmaDoubleComplex *w, magmaDoubleComplex *F2)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	if (r < nk)
		F2[r] = MAGMA_Z_DIV(MAGMA_Z_MAKE(2.0, 0.0), MAGMA_Z_ADD(s, MAGMA_Z_ADD(w[r], w[r])));
	return;
}

// Build F1,1(s)
__global__ void gpu_build_F11(int nk, int n, magmaDoubleComplex cfreq, magmaDoubleComplex *w, magmaDoubleComplex *F11)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;


	// Each block loads its elements into shared memory

	if (i < n)
	{
		int j;
		int s = nk - 1;
		int p,q;
		for (j = 2; j <= nk; ++j)
		{
			if (i < s)
			{
				p = j - 2;
				break;
			}
			s += nk - j;
		}
		s -= (nk - j) + 1;
		q = i - s + p + 1;
		F11[i] = MAGMA_Z_DIV(MAGMA_Z_ONE, MAGMA_Z_ADD(cfreq, MAGMA_Z_ADD(w[p], w[q])));
	}
	return;		
}

// Returns H = (F1 * wp^tr)^tr where: 
//  F1 = diagonal elements of F1 (nk x nk)
//  wp = matrix 3 x nk STORED COLUM WISE
//  H  = the result, 3 x nk, STORED COLUMN WISE
__global__ void gpu_F1_times_wp_tr(int dim, magmaDoubleComplex *F1, magmaDoubleComplex *wp, magmaDoubleComplex *H)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < 3 && c < dim)
		H[r + c * 3] = MAGMA_Z_MUL(F1[c], wp[r + c * 3]);
	return;
}

// Methods to make Cpq matrix

// z- and z+
// Remember: wp and wm are stored column wise
__device__ magmaDoubleComplex gpu_zMP(int mu, int nu, int albe, int pq, magmaDoubleComplex *dWmp)
{
	magmaDoubleComplex z = MAGMA_Z_ZERO;
	if (albe == mu)
		z = MAGMA_Z_ADD(z, dWmp[pq + nu * 3]);
	if (albe == nu)
		z = MAGMA_Z_ADD(z, dWmp[pq + mu * 3]);
	return z;
}

// NOTE: Cpq is stored COLUM WISE!
// Remember: wp and wm are stored column wise
__global__ void gpu_build_c1_11(int nk, int pidx, int qidx, magmaDoubleComplex freq, magmaDoubleComplex *dWm, magmaDoubleComplex *dWp, magmaDoubleComplex *dWeig, magmaDoubleComplex *dCpq)
{
	int alpha = threadIdx.x + blockIdx.x * blockDim.x;
	int beta  = threadIdx.y + blockIdx.y * blockDim.y;
	if (alpha < nk && beta < nk)
	{
		dCpq[alpha + beta * nk] = MAGMA_Z_ZERO;
		magmaDoubleComplex cab;
		for (unsigned int mu = 0; mu < nk; ++mu)
		{
			for (unsigned int nu = mu + 1; nu < nk; ++nu)
			{
				cab = MAGMA_Z_MUL(gpu_zMP(mu, nu, alpha, pidx, dWm), gpu_zMP(mu, nu, beta, qidx, dWp));
				cab = MAGMA_Z_DIV(cab, MAGMA_Z_ADD(freq, MAGMA_Z_ADD(dWeig[mu], dWeig[nu])));
				dCpq[alpha + beta * nk] = MAGMA_Z_ADD(dCpq[alpha + beta * nk], cab);
			}	
		}
	}
	return;
}

// Builds Mpq *** STORED COLUMN WISE
// WARNING: Cpq is destroyed, overwritten by Mpq
// Remember: wp and wm are stored column wise
__global__ void gpu_build_Mpq(int nk, int pidx, int qidx, magmaDoubleComplex *dWm, magmaDoubleComplex *dWp, magmaDoubleComplex *dF2, magmaDoubleComplex *dCpq)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < nk && c < nk && r==c)
		dCpq[r + c * nk] = MAGMA_Z_ADD(dCpq[r + c * nk], MAGMA_Z_MUL(dWm[pidx + r * 3], MAGMA_Z_MUL(dWp[qidx + r * 3], dF2[r])));
	return;
}

// Builds F1 x Mpq x F1
// WARNING: Mpq is destroyed, overwritten by the new natrix
// Remember: Mpq, and thus the resulting matrix, are STORED COLUMN WISE
__global__ void gpu_build_F1_Mpq_F1(int nk, int pidx, int qidx, magmaDoubleComplex *dF1, magmaDoubleComplex *dCpq)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c  = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < nk && c < nk)
		dCpq[r + c * nk] = MAGMA_Z_MUL(dF1[r], MAGMA_Z_MUL(dCpq[r + c * nk], dF1[c]));
	return;
}

// Convert Real matrix to Complex matrix
__global__ void gpu_RealToComplex(int n, double *Rmat, magmaDoubleComplex *Cmat)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c  = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < n && c < n)
		Cmat[r + c * n] = MAGMA_Z_MAKE(Rmat[r * n + c], 0.0);
	return;
}

// Takes the complex conjugate of a matrix
__global__ void gpu_matrixConj(int n, magmaDoubleComplex *Cmat)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c  = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < n && c < n)
		Cmat[r * n + c] = MAGMA_Z_CONJ(Cmat[r * n + c]);
	return;
}

// Copy a complex matrix into another
__global__ void gpu_copyA2B(int n, magmaDoubleComplex *A, magmaDoubleComplex *B)
{
	int r = threadIdx.x + blockIdx.x * blockDim.x;
	int c  = threadIdx.y + blockIdx.y * blockDim.y;
	if (r < n && c < n)
		B[r * n + c] = A[r * n + c];
	return;
}

// Reduction
__global__ void gpu_Sum_A(int a, int b, int g, int d, int n, magmaDoubleComplex *wm, magmaDoubleComplex *wp, magmaDoubleComplex *F1, magmaDoubleComplex * F2, magmaDoubleComplex *total)
{
	int tid = threadIdx.x;
	int i = blockIdx.x * blockDim.x + threadIdx.x;


	// Each block loads its elements into shared memory

	extern __shared__ magmaDoubleComplex x[];
	if (i < n)
	{
		x[tid] = MAGMA_Z_MUL(wm[a + i * 3], MAGMA_Z_MUL(wm[b + i * 3], MAGMA_Z_MUL(wp[g + i * 3], wp[d + i * 3])));
		x[tid] = MAGMA_Z_MUL(x[tid], MAGMA_Z_MUL(F2[i], MAGMA_Z_MUL(F1[i], F1[i]))); 
	}
	else
		x[tid] = MAGMA_Z_ZERO;
	__syncthreads();


	// Build summation tree over elements.

	for (int offset = blockDim.x / 2; offset > 0; offset >>= 1)
	{
		if(tid < offset) x[tid] = MAGMA_Z_ADD(x[tid], x[tid + offset]);
		__syncthreads();
	}


	// Thread 0 adds the partial sum to the total sum

	if( threadIdx.x == 0 )
		total[blockIdx.x] = x[0];

	return;
}

__global__ void gpu_Sum_B(int a, int b, int g, int d, int nk, int n, magmaDoubleComplex *wm, magmaDoubleComplex *wp, magmaDoubleComplex *F1, magmaDoubleComplex *F11, magmaDoubleComplex *total)
{
	int tid = threadIdx.x;
	int i = blockIdx.x * blockDim.x + threadIdx.x;


	// Each block loads its elements into shared memory

	extern __shared__ magmaDoubleComplex x[];
	if (i < n)
	{
		magmaDoubleComplex tmp;
		int j;
		int s = nk - 1;
		int p,q;
		for (j = 2; j <= nk; ++j)
		{
			if (i < s)
			{
				p = j - 2;
				break;
			}
			s += nk - j;
		}
		s -= (nk - j) + 1;
		q = i - s + p + 1;
		
/*
		x[tid] = MAGMA_Z_ADD(F1[p], F1[q]);
		x[tid] = MAGMA_Z_MUL(x[tid], x[tid]);
		x[tid] = MAGMA_Z_MUL(x[tid], MAGMA_Z_MUL(wm[a + p * 3], MAGMA_Z_MUL(wm[b + q * 3], MAGMA_Z_MUL(wp[g + q * 3], wp[d + p * 3]))));
		x[tid] = MAGMA_Z_MUL(x[tid], F11[i]);
*/

		tmp = MAGMA_Z_MUL(wm[a + q * 3], MAGMA_Z_MUL(wm[b + p * 3], MAGMA_Z_MUL(wp[g + p * 3], wp[d + q * 3])));
		tmp = MAGMA_Z_MUL(tmp, MAGMA_Z_MUL(F11[i],MAGMA_Z_MUL(F1[q], F1[q])));
		x[tid] = tmp;

		tmp = MAGMA_Z_MUL(wm[a + q * 3], MAGMA_Z_MUL(wm[b + p * 3], MAGMA_Z_MUL(wp[g + q * 3], wp[d + p * 3])));
		tmp = MAGMA_Z_MUL(tmp, MAGMA_Z_MUL(F11[i],MAGMA_Z_MUL(F1[q], F1[p])));
		x[tid] = MAGMA_Z_ADD(x[tid], tmp);

		tmp = MAGMA_Z_MUL(wm[a + p * 3], MAGMA_Z_MUL(wm[b + q * 3], MAGMA_Z_MUL(wp[g + p * 3], wp[d + q * 3])));
		tmp = MAGMA_Z_MUL(tmp, MAGMA_Z_MUL(F11[i],MAGMA_Z_MUL(F1[p], F1[q])));
		x[tid] = MAGMA_Z_ADD(x[tid], tmp);

		tmp = MAGMA_Z_MUL(wm[a + p * 3], MAGMA_Z_MUL(wm[b + q * 3], MAGMA_Z_MUL(wp[g + q * 3], wp[d + p * 3])));
		tmp = MAGMA_Z_MUL(tmp, MAGMA_Z_MUL(F11[i],MAGMA_Z_MUL(F1[p], F1[p])));
		x[tid] = MAGMA_Z_ADD(x[tid], tmp);
	}
	else
		x[tid] = MAGMA_Z_ZERO;
	__syncthreads();


	// Build summation tree over elements.

	for (int offset = blockDim.x / 2; offset > 0; offset >>=1)
	{
		if(tid < offset) x[tid] = MAGMA_Z_ADD(x[tid], x[tid + offset]);
		__syncthreads();
	}


	// Thread 0 adds the partial sum to the total sum

	if( threadIdx.x == 0 )
		total[blockIdx.x] = x[0];

	return;
}

__global__ void gpu_block_sum(magmaDoubleComplex *input, magmaDoubleComplex *results, size_t n)
{
	extern __shared__ magmaDoubleComplex sdata[];
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int tx = threadIdx.x;

	magmaDoubleComplex x = MAGMA_Z_ZERO;
	if (i < n)
		x = input[i];
	sdata[tx] = x;
	__syncthreads();

	for (int offset = blockDim.x / 2; offset > 0; offset >>= 1)
	{
		if (tx < offset)
			sdata[tx] = MAGMA_Z_ADD(sdata[tx], sdata[tx + offset]);
		__syncthreads();
	}

	if (threadIdx.x == 0)
		results[blockIdx.x] = sdata[0];

	return;
}
