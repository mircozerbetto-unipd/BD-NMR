#include "input_RK.h"

/***************/
/* Constructor */
/***************/
input::input()
{
	nStepsFound =  0;
	nDumpFound  =  0;
	nxFound     =  0;
	rcdFound    =  0;
	dtFound     =  0;
	seedFound   =  0;
	
	nErrors     =  0;

	nSteps      = 0;

	return;
}

/*****************/
/* Deconstructor */
/*****************/
input::~input()
{
}

/********************/
/* Parse input file */
/********************/
void input::parseInputFile(char* inputFileName)
{
	std::fstream f;
	f.open(inputFileName, std::ios::in);
	std::string line, keyword;
	std::stringstream sline(std::ios_base::in);
	std::cout << "**********************" << std::endl;
	std::cout << "* Parsing input file *" << std::endl;
	std::cout << "**********************" << std::endl;
	while (getline(f, line))
	{
		sline.clear();
		sline.str(line);
		sline >> keyword;
		transform(keyword.begin(), keyword.end(), keyword.begin(), ::tolower);

		if (!keyword.compare("nsteps"))
		{
			sline >> nSteps;
			if (nSteps <= 0)
			{
				std::cout << std::endl << "- nSteps will be estimated" << std::endl;
				//nErrors++;
			}
			else
				std::cout << "- nSteps: " << nSteps <<std::endl;

			nStepsFound = 1;
		}
		else if (!keyword.compare("ndump"))
		{
			sline >> nDump;
			if (nDump <= 0)
			{
				std::cout << std::endl << "- nDump will be estimated" << std::endl;
				//std::cout << std::endl << "** ERROR ** Number nDump must be > 0" << std::endl;
				//nErrors++;
			}
			else if (nDump > nSteps && nSteps > 0)
			{
				std::cout << std::endl << "** ERROR ** nDump must be smaller or equal to nSteps" << std::endl;
				nErrors++;
			}
			else
				std::cout << "- nDump: " << nDump << std::endl;

			nDumpFound = 1;
		}
		else if (!keyword.compare("nintcoor"))
		{
			sline >> nx;
			if (nx <= 0)
			{
				std::cout << "** ERROR ** in this implementation nintcorr must be strictly > 0" << std::endl;
				nErrors++;
			}
			else
				std::cout << "- Number of internal coordinates (nintcoor): " << nx << std::endl;
			nxFound = 1;
		}
		else if (!keyword.compare("rcdiff"))
		{
			sline >> RCDiffFile;
			std::cout << "- RC diffusion tensor fle: " << RCDiffFile << std::endl;
			rcdFound = 1;
		}
		else if (!keyword.compare("dt"))
		{
			sline >> dt;
			if (dt <= 0.0)
			{
				std::cout << "- dt will be estimated" << std::endl;
			}
			else
				std::cout << "- Integration time step dt: " << dt << " fs" << std::endl;
			dtFound = 1;
		}
		else if (!keyword.compare("seed"))
		{
			sline >> seed;
			if (seed < 1)
			{
				std::cout << "** WARNING ** seed <= 0 in the input file; the time and date will be used to generate the seed" << std::endl;
				seed = (uint64_t)time(0);
			}
			else
				std::cout << "- Seed for random numbers generator: " << seed << std::endl;
			seedFound = 1;
		}
		else
		{
			std::cout << "** ERROR ** unknown keyword: " << keyword << std::endl;
			nErrors++;
		}
	}
	std::cout << "**********************" << std::endl;
	std::cout << "* End of input       *" << std::endl;
	std::cout << "**********************" << std::endl;

	checkInputIntegrity();
	readDiffusionTensor();
	setInitialConfiguration();

	return;
}

/********************************/
/* Check integrity of the input */
/********************************/
void input::checkInputIntegrity()
{
	if (!nStepsFound)
	{
		std::cout << "** ERROR ** nSteps not found" << std::endl;
		nErrors++;
	}
	if (!nDumpFound)
	{
		std::cout << "** ERROR ** nDump not found" << std::endl;
		nErrors++;
	}
	if (!nxFound)
	{
		std::cout << "** ERROR ** nIntCoor not found" << std::endl;
		nErrors++;
	}
	if (!rcdFound)
	{
		std::cout << "** ERROR ** rcdFile not keyword not found" << std::endl;
		nErrors++;
	}
	if (!dtFound)
	{
		std::cout << "** ERROR ** dt not found" << std::endl;

		nErrors++;
	}
	if (!seedFound)
	{
		std::cout << "** ERROR ** seed for random numbers not found" << std::endl;
		nErrors++;
	}
	if (nErrors > 0)
		exit(1);

	return;
}

/********************************/
/* Read in the diffusion tensor */
/********************************/
void input::readDiffusionTensor(void)
{
	nc     = 3 + nx;
	D  = new double[nc * nc];

	scale = -1.0;

	std::fstream fD;
	fD.open(RCDiffFile.c_str(), std::ios::in);
	for (int i = 0; i < nc; ++i)
	{
		for (int j = 0; j < nc; ++j)
		{
			fD >> D[i * nc + j];
		}
		if ((D[i * nc + i] > scale))
			scale = D[i * nc + i];
	}
	fD.close();


	// Determine dt, nSteps, nDump if requested
	if (((int)dt < 0) || (nSteps <= 0) || (nDump <= 0))
	{

		double wmax = 0.0, wmin = 999999999.9;

		std::cout << "***********************************" << std::endl;
		std::cout << "Automatic generation of trajectory " << std::endl;
		std::cout << "***********************************" << std::endl << std::endl;
		std::cout << "Finding eigenvalues of the diffusion tensor" << std::endl;

		magma_init();
		magma_queue_t queue;
		magma_queue_create(0, &queue);

		magma_int_t infoK;
		double *h_workK;
		magma_int_t lworkK;
		magma_int_t *iworkK;
		magma_int_t liworkK;

		double *dD;
		magma_dmalloc(&dD, nc * nc);
		cudaMemcpy(dD, D, nc * nc * sizeof(double), cudaMemcpyHostToDevice);

		double *lambdaD = new double[nc];
		double *wD = new double[nc * nc];

		// Query for workspace sizes
		double aux_workK[1];
		magma_int_t aux_iworkK[1];
		magma_dsyevd_gpu (MagmaNoVec, MagmaLower, nc, dD, nc, lambdaD, wD, nc, aux_workK, -1, aux_iworkK, -1, &infoK);
		lworkK = (magma_int_t) aux_workK[0];
		liworkK = aux_iworkK[0];
		iworkK = (magma_int_t *) malloc (liworkK * sizeof (magma_int_t));
		magma_dmalloc_cpu (&h_workK, lworkK); // host mem. for workspace

		// Compute eigenvalues	
		magma_dsyevd_gpu(MagmaNoVec, MagmaLower, nc, dD, nc, lambdaD, wD, nc, h_workK, lworkK, iworkK, liworkK, &infoK);

		// Free memory
		magma_free(dD);
		magma_free_cpu(h_workK);
		free(iworkK);
		free(wD);
		magma_finalize();

		// Output the eigenvalues
		for (int i = 0; i < nc; ++i)
		{
			std::cout << std::fixed << std::setprecision(0) << i+1 << ") " << std::scientific << std::setprecision(8) << lambdaD[i] << std::endl;
			wmin = wmin < lambdaD[i] ? wmin : lambdaD[i];
			wmax = wmax > lambdaD[i] ? wmax : lambdaD[i];
		}
		std::cout << "***********************************" << std::endl << std::endl;
	
		int resolution = 10000;
		int length = 10;
		int nRep = 20;

		dt = 1.0 / (8.0 * wmax);
		double tcorr = 1.0 / wmin;
		std::cout << "tcorr / ns = " << tcorr * 1.e-6 << std::endl;
		double tdump = tcorr / (double)resolution;

		nDump = (uint64_t)floor(tdump / dt);
		nSteps = (uint64_t)resolution * (uint64_t)length * (uint64_t)nRep * (uint64_t)nDump;

		std::cout << "dt / fs = " << dt << std::endl;
		std::cout << "nSteps = " << nSteps << std::endl;
		std::cout << "nDump = " << nDump << std::endl;
		std::cout << "nChunks = " << nRep << std::endl;

		std::ofstream trjData("trjEstimatedInput.dat");
		trjData << "dt " << dt << std::endl;
		trjData << "nSteps " << nSteps << std::endl;
		trjData << "nDump " << nDump << std::endl;
		trjData << "nChunks " << nRep << std::endl;
		trjData.close();

		scale = wmax;
	}

	std::cout << std::endl;
	std::cout << "*******************************************" << std::endl;
	std::cout << "Scale factor for frequencies: " << scale << " fs^-1" << std::endl;
	std::cout << "*******************************************" << std::endl;

	// Scale the diffution tensor
	double iscale = 1.0 / scale;
	for (int i = (nc * nc - 1); i >= 0; --i)
		D[i] *= iscale;

	// Scale dt
	dt *= scale;
	if (dt > 0.5)
		std::cout << std::endl << "WARNING: dt * max(D) > 0.5 The simulation is likely to become instable. It is suggested to reduce the time step of integration." << std::endl;

	return;
}

/*********************************/
/* Set the initial configuration */
/*********************************/
// It is assumed no orientation and that
// the internal coordinates are initially
// in the zero energy configuration
void input::setInitialConfiguration(void)
{

	// No rotation
	qi.q0 = 1.0;
	qi.q1 = 0.0;
	qi.q2 = 0.0;
	qi.q3 = 0.0;

	// Minimum energy conformation
	xi = new double[nx];
	for (int i = nx - 1; i >= 0; --i)
		xi[i] = 0.0;

	return;
}

/*********************/
/* Return input data */
/*********************/
uint64_t input::getNSteps()
{
	return nSteps;
}

uint64_t input::getNDump()
{
	return nDump;
}

int input::getNx()
{
	return nx;
}

double input::getIntegrationDt()
{
	return dt;
}

void input::getRCDiffusionTensor(double *Dout)
{
	for (int i = 0; i < nc; ++i)
	{
		for (int j = 0; j < nc; ++j)
			Dout[i * nc + j] = D[i * nc + j];
	}

	return;
}

void input::getInitialState(Quaternion *qiout, double *xiout)
{
	qiout[0].q0 = qi.q0;
	qiout[0].q1 = qi.q1;
	qiout[0].q2 = qi.q2;
	qiout[0].q3 = qi.q3;

	for (int i = 0; i < nx; ++i)
		xiout[i] = xi[i];
}

uint64_t input::getSeed()
{
	return seed;
}

double input::getFreqScale()
{
	return scale;
}

double input::getTimeScale()
{
	return (1.0 / scale);
}

