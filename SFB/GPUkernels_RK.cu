/*******************************/
/* Zero out rectangular matrix */
/*******************************/
// To be optimized using 2D kernel
__global__ void gpuZeroOutRectMat(double *A, int nr, int nc)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < nr)
	{
		for (int j = 0; j < nc; ++j)
			A[i * nc + j] = 0.0;
	}

	return;
}

/*********************/
/* Fill the W matrix */
/*********************/
__global__ void gpuInitW(double *W, int nx)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < nx)
		W[(4 + i) * (3 + nx) + (3 + i)] = 1.0;

	return;
}

/************************************/
/* Generate the random displacement */
/************************************/
__device__ void buildSmallB(double* q, double* b)
{
	b[0] = -q[1] * 0.5;  b[1] = -q[2] * 0.5;  b[2] = -q[3] * 0.5; 
	b[3] =  q[0] * 0.5;  b[4] = -q[3] * 0.5;  b[5] =  q[2] * 0.5; 
	b[6] =  q[3] * 0.5;  b[7] =  q[0] * 0.5;  b[8] = -q[1] * 0.5; 
	b[9] = -q[2] * 0.5;  b[10] = q[1] * 0.5;  b[11] = q[0] * 0.5; 

	return;
}

__global__ void gpuPrepareRandomJump(double *states, int ndof, double *W, double *smallb)
{
        buildSmallB(states, smallb);
	for (int i = 0; i < 4; ++i)
	{
	        for (int j = 0; j < 3; ++j)
	            W[i * (ndof - 1) + j] = smallb[i * 3 + j];
	}

	return;
}

/********************************/
/* Calculate the harmonic force */
/********************************/
__global__ void gpuSFBHarmonicForce(double *states, double *F, int nX)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if (tid < nX)
		F[3 + tid] = -states[4 + tid];

	return;
}

/************************/
/* Compute states(t+dt) */
/************************/
// Virtual (same as Euler scheme propagation)
__global__ void gpuComputeNextState(double *states, double *drift, double *dt, double *rndJump, double *mult, double *oldStates, int ndof)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if (tid < ndof)
	{
		oldStates[tid] = states[tid];
		states[tid] += drift[tid] * dt[0] + rndJump[tid] * mult[0];
	}

	return;
}

// Runge-Kutta
__global__ void gpuComputeRKNextState(double *states, double *drift, double *drift_virtual, double *rndJump, double *rndJump_virtual, double *dt, double *mult, double *oldStates, int ndof)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if (tid < ndof)
	{
		states[tid] = oldStates[tid] + 0.5 * dt[0] * (drift[tid] + drift_virtual[tid]) + 0.5 * mult[0] * (rndJump[tid] + rndJump_virtual[tid]);
	}

	return;
}

__global__ void gpuNormalizeQuaternion(double *states, double lambdaQ)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if (tid < 4)
		states[tid] += lambdaQ * states[tid];

	return;
}

