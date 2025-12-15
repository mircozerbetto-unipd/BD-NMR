#include "sfb_RK.h"

#define TPB_1D 512 // Threads per block in 1D geometry
#define BPG    30  // Blocks per grid

/////////////////////////////////////////////
// Make CURAND verbose to understan errors //
/////////////////////////////////////////////

static const char *curandGetErrorString(curandStatus_t error)
{
    switch (error)
    {
        case CURAND_STATUS_SUCCESS:
            return "CURAND_STATUS_SUCCESS";

        case CURAND_STATUS_VERSION_MISMATCH:
            return "CURAND_STATUS_VERSION_MISMATCH";

        case CURAND_STATUS_NOT_INITIALIZED:
            return "CURAND_STATUS_NOT_INITIALIZED";

        case CURAND_STATUS_ALLOCATION_FAILED:
            return "CURAND_STATUS_ALLOCATION_FAILED";

        case CURAND_STATUS_TYPE_ERROR:
            return "CURAND_STATUS_TYPE_ERROR";

        case CURAND_STATUS_OUT_OF_RANGE:
            return "CURAND_STATUS_OUT_OF_RANGE";

        case CURAND_STATUS_LENGTH_NOT_MULTIPLE:
            return "CURAND_STATUS_LENGTH_NOT_MULTIPLE";

        case CURAND_STATUS_DOUBLE_PRECISION_REQUIRED:
            return "CURAND_STATUS_DOUBLE_PRECISION_REQUIRED";

        case CURAND_STATUS_LAUNCH_FAILURE:
            return "CURAND_STATUS_LAUNCH_FAILURE";

        case CURAND_STATUS_PREEXISTING_FAILURE:
            return "CURAND_STATUS_PREEXISTING_FAILURE";

        case CURAND_STATUS_INITIALIZATION_FAILED:
            return "CURAND_STATUS_INITIALIZATION_FAILED";

        case CURAND_STATUS_ARCH_MISMATCH:
            return "CURAND_STATUS_ARCH_MISMATCH";

        case CURAND_STATUS_INTERNAL_ERROR:
            return "CURAND_STATUS_INTERNAL_ERROR";
    }

    return "<unknown>";
}

//////////////////////////////////////////////////////////
// NB: since cuBLAS it is based on FORTRAN, matrices    //
//     are interpreted in column-major order.           //
//     Therefore AB = C is computed as                  //
//                                                      //
//                     B^tr A^tr = C^tr                 //
//                                                      //
//     this is achieved by passing B and A to cuBLAS    //
//     in row-major order, with no transposition.       //
//     In reading in the matrices, cuBLAS interprets    //
//     them as transposed. The output, which is C^tr    //
//     for cuBLAS, is simply C read in row-major order. //
//////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
	// Check command line options
	if (argc < 2)
	{
		std::cout << std::endl << "ERROR: an input file must be specified when calling SFB" << std::endl << std::endl;
		exit(1);
	}

	////////////////////////
	// Init MAGMA library //
	////////////////////////

	magma_init();
	magma_queue_t queue;
	magma_queue_create(0, &queue);

	///////////
	// Input //
	///////////

	input BDsetup;
	BDsetup.parseInputFile(argv[1]);

	///////////////////////////
	// Initial configuration //
	///////////////////////////
	
	int ndof = 4 + BDsetup.getNx();
	int nRC =  3 + BDsetup.getNx();
	int nRC_even = (nRC % 2 == 0) ? nRC : nRC + 1; // the normal rnd numbers generator requires an even array length
	double *xi, *states, *d_states;

	Quaternion qi = Quaternion();
	xi = new double[BDsetup.getNx()];
	states = new double[ndof];

	memset(states, 0, ndof * sizeof(double));

	BDsetup.getInitialState(&qi, xi);
	states[0] = qi.q0;
	states[1] = qi.q1;
	states[2] = qi.q2;
	states[3] = qi.q3;
	for (int i = 0; i < BDsetup.getNx(); ++i)
		states[4 + i] = xi[i];

	cudaMalloc((void**)&d_states, ndof * sizeof(double)); 
	cudaMemcpy(d_states, states, ndof * sizeof(double), cudaMemcpyHostToDevice);

	////////////////////////////////////////////////
	// Cholesky decomposition of diffusion tensor //
	////////////////////////////////////////////////

	double *D = new double [nRC * nRC];
    	double *L = new double [nRC * nRC]; // Cholesky decomposition, lower diagonal

	BDsetup.getRCDiffusionTensor(D);
	BDsetup.getRCDiffusionTensor(L);
#ifdef DEBUG
	std::cout << "* Input diffusion tensor: *" << std::endl;
	for (int i = 0; i < nRC; ++i)
	{
		for (int j = 0; j < nRC; ++j)
			std::cout << D[i * nRC + j] << "\t";
		std::cout << std::endl;
	}
	std::cout << "************************************" << std::endl << std::endl;
#endif

	int info;
    	magma_dpotrf(MagmaUpper, nRC, L, nRC, &info); // D = L * L^tr
	if (info != 0)
	{
		std::cerr << "Cholesky decomposition failed with error: " << info << std::endl;
		exit(-1);
	}
	for (int i = 0; i < nRC; i++)
	{
		for (int j = i + 1; j < nRC; j++)
			L[i * nRC + j] = 0.0;
	}

#ifdef DEBUG
	std::cout << "* Cholesky factor L *" << std::endl;
	for (int i = 0; i < nRC; i++)
	{
		for (int j = 0; j < nRC; j++)
			std::cout << L[i * nRC + j] << "\t";
		std::cout << std::endl;
	}
	std::cout << "************************************" << std::endl << std::endl;
#endif

	magma_queue_destroy(queue);
	magma_finalize();

	double *d_D, *d_L; // GPU copies of D and L

	cudaMalloc((void**)&d_D, nRC * nRC * sizeof(double));
	cudaMalloc((void**)&d_L, nRC * nRC * sizeof(double));

	cudaMemcpy(d_D, D, nRC * nRC * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_L, L, nRC * nRC * sizeof(double), cudaMemcpyHostToDevice);

	///////////////////////////////////////////////////////
	// Conversion matrix from quaternion to Euler angles //
	///////////////////////////////////////////////////////
	
	double *d_W;
	cudaMalloc((void**)&d_W, ndof * nRC * sizeof(double)); 
	gpuZeroOutRectMat<<<1 + ndof / TPB_1D, TPB_1D>>>(d_W, ndof, nRC);
	gpuInitW<<<1 + BDsetup.getNx() / TPB_1D, TPB_1D>>>(d_W, BDsetup.getNx());

	/////////////////////////////////////
	// Init force/torque on the device //
	/////////////////////////////////////
	
	double *d_F;
	cudaMalloc((void**)&d_F, nRC * sizeof(double)); 
	gpuZeroOutRectMat<<<1 + nRC / TPB_1D, TPB_1D>>>(d_F, nRC, 1);
#ifdef DEBUG
	double *F;
	F = new double[nRC];
#endif

	/////////////////////////////////
	// Setup random numbers array  //
	/////////////////////////////////
	
	double *d_rnd;
	cudaMalloc(&d_rnd, nRC_even * sizeof(double));
#ifdef DEBUG
	double *rnd;
	rnd = new double[nRC_even];
#endif

	/////////////////////////////////////
	// Allocate matrices on the device //
	/////////////////////////////////////

	double *d_smallb, *d_bigB;
	cudaMalloc((void**)&d_smallb, 4 * 3 * sizeof(double));
	cudaMalloc((void**)&d_bigB, 4 * 4 * sizeof(double));

	double *d_rndJump, *d_rndJump_virtual, *d_drift, *d_drift_virtual;
	cudaMalloc((void**)&d_rndJump, ndof * sizeof(double));
	cudaMalloc((void**)&d_rndJump_virtual, ndof * sizeof(double)); // Random jump from the virtual state in Runge-Kutta time integration
	cudaMalloc((void**)&d_drift, ndof * sizeof(double));
	cudaMalloc((void**)&d_drift_virtual, ndof * sizeof(double)); // Drift from the virtual state in Runge-Kutta time integration
	
	double *d_WL;
	cudaMalloc((void**)&d_WL, ndof * nRC * sizeof(double));

	double *d_oldStates;
	cudaMalloc((void**)&d_oldStates, ndof * sizeof(double));

	///////////////////////
	// Output trajectory //
	///////////////////////
	
	std::ofstream trjFile("trj.dat"); // To get from input

	////////////////////////////
	// Propagate the dynamics //
	////////////////////////////

	double mult, *d_mult;
	mult = std::sqrt(2.0 * BDsetup.getIntegrationDt());
	cudaMalloc((void**)&d_mult, sizeof(double));
	cudaMemcpy(d_mult, &mult, sizeof(double), cudaMemcpyHostToDevice);

	double dt, *d_dt;
	dt = BDsetup.getIntegrationDt();
	cudaMalloc((void**)&d_dt, sizeof(double));
	cudaMemcpy(d_dt, &dt, sizeof(double), cudaMemcpyHostToDevice);

	double time = 0.0, bNorm, cNorm, lambdaQ;
    	const double delta = 0.0;
	const double scale = 1.0;

	cublasHandle_t handle;
	cublasCreate(&handle);

	EulerQuaternion EQ;
	EulerAngles ea;

	////////////////////////////////////
	// Setup random numbers generator //
	////////////////////////////////////

	curandGenerator_t randGen;
	curandCreateGenerator(&randGen, CURAND_RNG_PSEUDO_DEFAULT);

	//////////////////////////////////
	// Output initial configuration //
	//////////////////////////////////

	cudaMemcpy(states, d_states, ndof * sizeof(double), cudaMemcpyDeviceToHost);
	qi.q0 = states[0]; qi.q1 = states[1]; qi.q2 = states[2]; qi.q3 = states[3];
	EQ.setQuaternion(qi);
	ea = EQ.getEulerAngles();

#ifdef DEBUG
	std::cout << "* Initial configuration *" << std::endl;
	for (int ir = 0; ir < ndof; ++ir)
		std::cout << states[ir] << std::endl;
	std::cout << "************************************" << std::endl << std::endl;
#endif

#ifdef DEBUG
	double *rndJump = new double[ndof];
	double *drift = new double[ndof];
	double *W = new double [ndof * nRC];
	double *WL = new double [ndof * nRC];
#endif

	std::cout << std::endl;
	std::cout << "___________________________________" << std::endl;
	std::cout << "                                   " << std::endl;
	std::cout << "        Start BD trajetory         " << std::endl;
	std::cout << "___________________________________" << std::endl;
	std::cout << std::endl;

#ifdef DEBUG
	for (uint64_t it = 1; it <= 3; ++it)
#else
	for (uint64_t it = 1; it <= BDsetup.getNSteps(); ++it)
#endif
	{
		time += BDsetup.getIntegrationDt(); // this is scaled by 1/Dzz

		/**************************************************************/	
		/* Prepare the sequence of random numbers to use in this step */
		/**************************************************************/

		// Generate nRC_even random normal numbers, rnd
		curandStatus_t status = curandGenerateNormalDouble(randGen, d_rnd, nRC_even, 0.0, 1.0);
		if (status != CURAND_STATUS_SUCCESS)
		{
			printf("CURAND normal double generation failed!\n");
			printf("CURAND error code: %s\n", curandGetErrorString(status));
			// Add more debugging information here, like the value of nRC_even
			printf("nRC_even value when error occurred: %zu\n", nRC_even);
			// Potentially exit if generation fails:
			return 1; // Or handle error as appropriate in your program
		}

#ifdef DEBUG
		cudaMemcpy(rnd, d_rnd, nRC_even * sizeof(double), cudaMemcpyDeviceToHost);
		// PS print only nRC (not nRC_even), since they are those effectively employed
		std::cout << "* Random numbers used in this step *" << std::endl;
		for (int ir = 0; ir < nRC; ++ir)
			std::cout << rnd[ir] << std::endl;
		std::cout << "************************************" << std::endl << std::endl;
#endif

		/* =============================== */
		/* === A) VIRTUAL DISPLACEMENT === */
		/* =============================== */

		/************************/	
		/* Random displacements */
		/************************/

		// Update W matrix
		gpuPrepareRandomJump<<<1, 1>>>(d_states, ndof, d_W, d_smallb);
#ifdef DEBUG
		cudaMemcpy(W, d_W, ndof * nRC * sizeof(double), cudaMemcpyDeviceToHost);
		std::cout << "* W at this step *" << std::endl;
		for (int ir = 0; ir < ndof; ++ir)
		{
			for (int ic = 0; ic < nRC; ++ic)
				std::cout << W[ir * nRC + ic] << "\t";
			std::cout<<std::endl;
		}
		std::cout << "************************************" << std::endl << std::endl;
#endif

		// WL = W * L
		cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nRC, ndof, nRC, &scale, d_L, nRC, d_W, nRC, &delta, d_WL, nRC);
#ifdef DEBUG
		cudaMemcpy(WL, d_WL, ndof * nRC * sizeof(double), cudaMemcpyDeviceToHost);
		std::cout << "* WL at this step *" << std::endl;
		for (int ir = 0; ir < ndof; ++ir)
		{
			for (int ic = 0; ic < nRC; ++ic)
				std::cout << WL[ir * nRC + ic] << "\t";
			std::cout<<std::endl;
		}
		std::cout << "************************************" << std::endl << std::endl;
#endif

		// rndJump = W * L * rnd = WL * rnd
		cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 1, ndof, nRC, &scale, d_rnd, 1, d_WL, nRC, &delta, d_rndJump, 1);

#ifdef DEBUG
		cudaMemcpy(rndJump, d_rndJump, ndof * sizeof(double), cudaMemcpyDeviceToHost);
		std::cout << "* Random displacement of this step *" << std::endl;
		for (int ir = 0; ir < ndof; ++ir)
			std::cout << rndJump[ir] << std::endl;
		std::cout << "************************************" << std::endl << std::endl;
#endif

		/*********/
		/* Drift */
		/*********/
		// Compute force
		gpuSFBHarmonicForce<<<1 + BDsetup.getNx() / TPB_1D, TPB_1D>>>(d_states, d_F, BDsetup.getNx());
#ifdef DEBUG
		cudaMemcpy(F, d_F, nRC * sizeof(double), cudaMemcpyDeviceToHost);
		std::cout << "* Force at this step *" << std::endl;
		for (int ir = 0; ir < nRC; ++ir)
			std::cout << F[ir] << std::endl;
		std::cout << "************************************" << std::endl << std::endl;
#endif
		// WL = W * D
		cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nRC, ndof, nRC, &scale, d_D, nRC, d_W, nRC, &delta, d_WL, nRC);
#ifdef DEBUG
		cudaMemcpy(WL, d_WL, ndof * nRC * sizeof(double), cudaMemcpyDeviceToHost);
		std::cout << "* WD at this step *" << std::endl;
		for (int ir = 0; ir < ndof; ++ir)
		{
			for (int ic = 0; ic < nRC; ++ic)
				std::cout << WL[ir * nRC + ic] << "\t";
			std::cout<<std::endl;
		}
		std::cout << "************************************" << std::endl << std::endl;
#endif
		// drift = WL * F
		cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 1, ndof, nRC, &scale, d_F, 1, d_WL, nRC, &delta, d_drift, 1);
#ifdef DEBUG
		cudaMemcpy(drift, d_drift, ndof * sizeof(double), cudaMemcpyDeviceToHost);
		std::cout << "* Drift at this step *" << std::endl;
		for (int ir = 0; ir < ndof; ++ir)
			std::cout << drift[ir] << std::endl;
		std::cout << "************************************" << std::endl << std::endl;
#endif

		/*************************************************/
		/* Compute next state by shift + q normalization */
		/*************************************************/
		// Shift
		// At this point: d_oldstates = x_k, d_drift = drift(x_k), d_rndJump = rndJump(x_k), and d_states = (virtual x_k+1 point)
		gpuComputeNextState<<<1 + ndof / TPB_1D, TPB_1D>>>(d_states, d_drift, d_dt, d_rndJump, d_mult, d_oldStates, ndof); 
#ifdef DEBUG
		cudaMemcpy(states, d_states, ndof * sizeof(double), cudaMemcpyDeviceToHost);
		std::cout << "* Next virtual state (not normalized) *" << std::endl;
		for (int ir = 0; ir < ndof; ++ir)
			std::cout << states[ir] << std::endl;
		std::cout << "************************************" << std::endl << std::endl;
#endif

		// Normalization using Laplace multiplier
		cublasDdot(handle, 4, d_states, 1, d_oldStates, 1, &bNorm);
		bNorm *= 2.0;

		cublasDdot(handle, 4, d_states, 1, d_states, 1, &cNorm);
		cNorm -= 1.0;

		lambdaQ = 0.5 * (-bNorm + std::sqrt(bNorm * bNorm  - 4.0 * cNorm));
#ifdef DEBUG
		std::cout << "* Lambda Q *" << std::endl;
		std::cout << lambdaQ << std::endl;
		std::cout << "************************************" << std::endl << std::endl;
#endif

		gpuNormalizeQuaternion<<<1 + 4 / TPB_1D, TPB_1D>>>(d_states, lambdaQ);
#ifdef DEBUG
		cudaMemcpy(states, d_states, ndof * sizeof(double), cudaMemcpyDeviceToHost);
		std::cout << "* Next virtual state (normalized) *" << std::endl;
		for (int ir = 0; ir < ndof; ++ir)
			std::cout << states[ir] << std::endl;
		std::cout << "Norm of quaternion: " << std::sqrt(states[0] * states[0] + states[1] * states[1] + states[2] * states[2] + states[3] * states[3]) << std::endl;
		std::cout << "************************************" << std::endl << std::endl;
#endif

		/* =================================== */
		/* === B) RUNGE-KUTTA DISPLACEMENT === */
		/* =================================== */

		/************************/	
		/* Random displacements */
		/************************/

		// Update W matrix
		gpuPrepareRandomJump<<<1, 1>>>(d_states, ndof, d_W, d_smallb);
#ifdef DEBUG
		cudaMemcpy(W, d_W, ndof * nRC * sizeof(double), cudaMemcpyDeviceToHost);
		std::cout << "* W at this step *" << std::endl;
		for (int ir = 0; ir < ndof; ++ir)
		{
			for (int ic = 0; ic < nRC; ++ic)
				std::cout << W[ir * nRC + ic] << "\t";
			std::cout<<std::endl;
		}
		std::cout << "************************************" << std::endl << std::endl;
#endif

		// WL = W * L
		cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nRC, ndof, nRC, &scale, d_L, nRC, d_W, nRC, &delta, d_WL, nRC);
#ifdef DEBUG
		cudaMemcpy(WL, d_WL, ndof * nRC * sizeof(double), cudaMemcpyDeviceToHost);
		std::cout << "* WL at this step *" << std::endl;
		for (int ir = 0; ir < ndof; ++ir)
		{
			for (int ic = 0; ic < nRC; ++ic)
				std::cout << WL[ir * nRC + ic] << "\t";
			std::cout<<std::endl;
		}
		std::cout << "************************************" << std::endl << std::endl;
#endif

		// rndJump_virtual = W * L * rnd = WL * rnd
		// This is the random jump from the virtual x_k+1 state (stored in d_states)
		cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 1, ndof, nRC, &scale, d_rnd, 1, d_WL, nRC, &delta, d_rndJump_virtual, 1);

#ifdef DEBUG
		cudaMemcpy(rndJump, d_rndJump_virtual, ndof * sizeof(double), cudaMemcpyDeviceToHost);
		std::cout << "* Random displacement from the virtual state of this step *" << std::endl;
		for (int ir = 0; ir < ndof; ++ir)
			std::cout << rndJump[ir] << std::endl;
		std::cout << "************************************" << std::endl << std::endl;
#endif

		/*********/
		/* Drift */
		/*********/
		// Compute force
		gpuSFBHarmonicForce<<<1 + BDsetup.getNx() / TPB_1D, TPB_1D>>>(d_states, d_F, BDsetup.getNx());
#ifdef DEBUG
		cudaMemcpy(F, d_F, nRC * sizeof(double), cudaMemcpyDeviceToHost);
		std::cout << "* Force at this step *" << std::endl;
		for (int ir = 0; ir < nRC; ++ir)
			std::cout << F[ir] << std::endl;
		std::cout << "************************************" << std::endl << std::endl;
#endif
		// WL = W * D
		cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nRC, ndof, nRC, &scale, d_D, nRC, d_W, nRC, &delta, d_WL, nRC);
#ifdef DEBUG
		cudaMemcpy(WL, d_WL, ndof * nRC * sizeof(double), cudaMemcpyDeviceToHost);
		std::cout << "* WD at this step *" << std::endl;
		for (int ir = 0; ir < ndof; ++ir)
		{
			for (int ic = 0; ic < nRC; ++ic)
				std::cout << WL[ir * nRC + ic] << "\t";
			std::cout<<std::endl;
		}
		std::cout << "************************************" << std::endl << std::endl;
#endif
		// drift_virtual = WL * F
		// This is the drift from the virtual x_k+1 state (stored in d_states)
		cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 1, ndof, nRC, &scale, d_F, 1, d_WL, nRC, &delta, d_drift_virtual, 1);
#ifdef DEBUG
		cudaMemcpy(drift, d_drift_virtual, ndof * sizeof(double), cudaMemcpyDeviceToHost);
		std::cout << "* Drift from the virtual state of this step *" << std::endl;
		for (int ir = 0; ir < ndof; ++ir)
			std::cout << drift[ir] << std::endl;
		std::cout << "************************************" << std::endl << std::endl;
#endif

		/*************************************************/
		/* Compute next state by shift + q normalization */
		/*************************************************/
		// Shift
		// states = oldStates + 0.5 * (drift + drift_virtual) + 0.5 * mult * (rndJump + rndJump_virtual)
		// at the end, oldStates = x_k and states = x_k+1 (the virtual state is destroyed)
		gpuComputeRKNextState<<<1 + ndof / TPB_1D, TPB_1D>>>(d_states, d_drift, d_drift_virtual, d_rndJump, d_rndJump_virtual, d_dt, d_mult, d_oldStates, ndof); 
#ifdef DEBUG
		cudaMemcpy(states, d_states, ndof * sizeof(double), cudaMemcpyDeviceToHost);
		std::cout << "* Next state (not normalized) *" << std::endl;
		for (int ir = 0; ir < ndof; ++ir)
			std::cout << states[ir] << std::endl;
		std::cout << "************************************" << std::endl << std::endl;
#endif

		// Normalization using Laplace multiplier
		cublasDdot(handle, 4, d_states, 1, d_oldStates, 1, &bNorm);
		bNorm *= 2.0;

		cublasDdot(handle, 4, d_states, 1, d_states, 1, &cNorm);
		cNorm -= 1.0;

		lambdaQ = 0.5 * (-bNorm + std::sqrt(bNorm * bNorm  - 4.0 * cNorm));
#ifdef DEBUG
		std::cout << "* Lambda Q *" << std::endl;
		std::cout << lambdaQ << std::endl;
		std::cout << "************************************" << std::endl << std::endl;
#endif

		gpuNormalizeQuaternion<<<1 + 4 / TPB_1D, TPB_1D>>>(d_states, lambdaQ);
#ifdef DEBUG
		cudaMemcpy(states, d_states, ndof * sizeof(double), cudaMemcpyDeviceToHost);
		std::cout << "* Next state (normalized) *" << std::endl;
		for (int ir = 0; ir < ndof; ++ir)
			std::cout << states[ir] << std::endl;
		std::cout << "Norm of quaternion: " << std::sqrt(states[0] * states[0] + states[1] * states[1] + states[2] * states[2] + states[3] * states[3]) << std::endl;
		std::cout << "************************************" << std::endl << std::endl;
#endif

		/**********************/
		/* Dump configuration */
		/**********************/
		if (!(it%BDsetup.getNDump()))
		{
			std::cout << "*******************************************" << std::endl;
			std::cout << "Step: " << it << ", time: " << time  * BDsetup.getTimeScale() * 0.001 << " ps" << std::endl;
			std::cout << "*******************************************" << std::endl;

			cudaMemcpy(states, d_states, ndof * sizeof(double), cudaMemcpyDeviceToHost);
			qi.q0 = states[0]; qi.q1 = states[1]; qi.q2 = states[2]; qi.q3 = states[3];
			EQ.setQuaternion(qi);
			ea = EQ.getEulerAngles();

			trjFile << it << "\t" << time * BDsetup.getTimeScale() * 0.001 << "\t"; // in ps
			trjFile << std::scientific << std::setprecision(8) << ea.alpha << "\t" << ea.beta << "\t" << ea.gamma;
			for (int j = 0; j < BDsetup.getNx(); ++j)
				trjFile << "\t" << std::scientific << std::setprecision(8) << states[4 + j];
			trjFile << std::endl;
		}
	}	

	trjFile.close();

	return 0;
}

