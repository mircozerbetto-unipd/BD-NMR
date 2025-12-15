#include "inputParser.h"

/***************/
/* Constructor */
/***************/
inputParserSrb::inputParserSrb()
{
	pdbFileNameFound   = 0;
	idsFound           = 0;
	ReffFound          = 0;
	CFound             = 0;
	viscosityFound     = 0;
	temperatureFound   = 0;
	hydrointFound      = 0;
	nfreqFound         = 0;
	deltafreqFound     = 0;
	rankFound          = 0;
	rigidBodyFound     = 0;
	drrFound           = 0;
	nmrFound           = 0;
	probeFound         = 0;
	dcsaFound          = 0;
	hatomFound         = 0;
	expFound           = 0; 
	fieldFound         = 0;
	kMethodFound       = 0;
	allActiveFound     = 0;
	ignoreNegEigsFound = 0;

	calculateOrderParameters = 0;

	C = 6.0;
	Reff = 2.0;
	temperature = 298.15;
	hydroint = 0;

	ids = new int[4];

	ignoreNegEigs = 1;

	return;
}

/*****************/
/* Deconstructor */
/*****************/
inputParserSrb::~inputParserSrb()
{
	return;
}

/***********/
/* METHODS */
/***********/
const char* inputParserSrb::getPDBFile(void)
{
	return pdbFileName.c_str();
}

void inputParserSrb::getRefAtoms(int *atID)
{
	atID[0] = ids[0];
	atID[1] = ids[1];
	atID[2] = ids[2];
	atID[3] = ids[3];
	return;
}

double inputParserSrb::getReff(void)
{
	return Reff; // A
}

double inputParserSrb::getC(void)
{
	return C;
}

double inputParserSrb::getViscosity(void)
{
	return viscosity; // Pa s
}

double inputParserSrb::getTemperature(void)
{
	return temperature; // K
}

int inputParserSrb::getHydrodynamicInteractions(void)
{
	return hydroint; // 0: off | 1: on [Rotne-Prage model]
}

int inputParserSrb::getNfreq(void)
{
	return nfreq;
}

double inputParserSrb::getDeltafreq(void)
{
	return deltafreq; // fs^-1
}

int inputParserSrb::getRank(void)
{
	return rank;
}

int inputParserSrb::getRigidBody(void)
{
	return rigidBodyFound;
}

double inputParserSrb::getDxx(void)
{
	return Dxx;
}

double inputParserSrb::getDyy(void)
{
	return Dyy;
}

double inputParserSrb::getDzz(void)
{
	return Dzz;
}

int inputParserSrb::getNMR(void)
{
	return nmrFound;
}

double inputParserSrb::getDeltaCSA(void)
{
	return dcsa;
}

std::string inputParserSrb::getProbe(void)
{
	return probe;
}

std::string inputParserSrb::getExpFile(void)
{
	if (expFound)
		return exp;
	else
		return "nofile";
}

int inputParserSrb::getNFields(void)
{
	return nFields;
}

void inputParserSrb::getFields(double *f)
{
	for (int i = 0; i < nFields; ++i)
		f[i] = fields[i];
	return;
}

int inputParserSrb::getAllActive(void)
{
	return allActive;
}

// Second H atom in CH2 probe.
// CH_Hid is 1-based, original PDB numbering, ID
void inputParserSrb::setHatom(int h)
{
	CH_Hid = h;
	return;
}
int inputParserSrb::getHatom(void)
{
	return CH_Hid;
}

std::string inputParserSrb::getKMatrixMethod(void)
{
	return kMethod;
}

std::string inputParserSrb::getHessianFile(void)
{
	return kFile1;
}

std::string inputParserSrb::getGradientsFile(void)
{
	return kFile2;
}

std::string inputParserSrb::getCovarFile(void)
{
	return kFile1;
}

// Order parameters
int inputParserSrb::makeOrderParameters(void)
{
	return calculateOrderParameters;
}

vectorOfIntegers inputParserSrb::getHid(void)
{
	return Hid;
}

vectorOfIntegers inputParserSrb::getNid(void)
{
	return Nid;
}

vectorOfIntegers inputParserSrb::getCid(void)
{
	return Cid;
}

// Parse input file
void inputParserSrb::parseInputFile(char* inputFileName)
{
	int tmpInt;
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

		if (!keyword.compare("pdb"))
		{
			sline >> pdbFileName;
			pdbFileNameFound = 1;
			std::cout << "- Input PDB file: " << pdbFileName << std::endl;
		}
		else if (!keyword.compare("refatoms"))
		{
			sline >> ids[0] >> ids[1] >> ids[2] >> ids[3];
			idsFound = 1;
			std::cout <<"- ID of reference atoms: " << ids[0] << " " << ids[1] << " " << ids[2] << " " << ids[3] << std::endl;
		}
		else if (!keyword.compare("reff"))
		{
			sline >> Reff;
			ReffFound = 1;
			std::cout << "- Effective radius for friction tensor: " << Reff << " Angstroms" << std::endl;
		}
		else if (!keyword.compare("c"))
		{
			sline >> C;
			CFound = 1;
			std::cout <<"- C (hydrodynamics boundary conditions): " << C << std::endl;
		}
		else if (!keyword.compare("viscosity"))
		{
			sline >> viscosity;
			viscosityFound = 1;
			std::cout << "- Medium viscosity: " << viscosity << " Pa s" << std::endl;
		}
		else if (!keyword.compare("temperature"))
		{
			sline >> temperature;
			temperatureFound = 1;
			std::cout << "- Temperature: " << temperature << " K" << std::endl;
		}
		else if (!keyword.compare("hydroint"))
		{
			sline >> hydroint;
			hydrointFound = 1;
			if (!hydroint) std::cout << "- Hydrodynamic Interactions: OFF" << std::endl;
			if (hydroint) std::cout << "- Hydrodynamic Interactions: ON [Rotne-Prage model]" << std::endl;
			
		}
		else if (!keyword.compare("nfreq"))
		{
			sline >> nfreq;
			nfreqFound = 1;
			std::cout << "- Number of spectral density evaluations: " << nfreq << std::endl;
		}
		else if (!keyword.compare("deltafreq"))
		{
			sline >> deltafreq;
			deltafreqFound = 1;
			std::cout << "- Frequency interval: " << deltafreq << "fs^-1" << std::endl;
		}
		else if (!keyword.compare("rank"))
		{
			sline >> rank;
			rankFound = 1;
			std::cout << "- Rank of physical observable: " << rank << std::endl;
		}
		else if (!keyword.compare("rigidbody"))
		{
			sline >> rigidBodyFound;
			std::cout << "- Calculating rigid body J(w): " << (rigidBodyFound ? "yes" : "no") << std::endl;
		}
		else if (!keyword.compare("drr"))
		{
			sline >> Dxx >> Dyy >> Dzz;
			// Convert to fs^-1
			Dxx *= 1.0e-15;
			Dyy *= 1.0e-15;
			Dzz *= 1.0e-15;
			drrFound = 1;
			std::cout << "- Rigid body rot. diff. tensor / fs: " << Dxx << ", " << Dyy << ", " << Dzz << std::endl;
		}
		else if (!keyword.compare("nmr"))
		{
			sline >> nmrFound;
			if (nmrFound < 0 || nmrFound > 1)
			{
				std::cout << std::endl << std::endl << "ERROR in input file: nmr keyword can assume only two values: 0 (no) or 1 (yes)." << std::endl << std::endl;
				exit(1);
			}	
			std::cout << "- Calculating NMR data: " << (nmrFound ? "yes" : "no") << std::endl;
		}
		else if (!keyword.compare("probe"))
		{
			sline >> probe;
			transform(probe.begin(), probe.end(), probe.begin(), ::tolower);
			if (probe.compare("nh") && probe.compare("ch") && probe.compare("ch2"))
			{
				std::cout << std::endl << std::endl << "ERROR in input file: only NH, CH and CH2 probes are available at present." << std::endl << std::endl;
				exit(1);
			}	
			probeFound = 1;
			std::cout << "- NMR probe: " << probe << std::endl;
		}
		else if (!keyword.compare("deltacsa"))
		{
			sline >> dcsa;
			dcsaFound = 1;
			std::cout << "- Delta CSA: " << dcsa << " ppm" << std::endl;
		}
		else if (!keyword.compare("hatom"))
		{
			sline >> CH_Hid;
			hatomFound = 1;
			std::cout << "- Second H id of CH2 probe: " << CH_Hid << std::endl;
		}
		else if (!keyword.compare("expfile"))
		{
			sline >> exp;
			expFound = 1;
			std::cout << "- File with experimental values: " << exp << std::endl;
		}
		else if (!keyword.compare("freq"))
		{
			sline >> nFields;
			fields = new double[nFields];
			std::cout << "- NMR data will be calculated at these frequencies (in MHz): ";
			for (int i = 0; i < nFields-1; ++i)
			{
				sline >> fields[i];
				std::cout << fields[i] << ", ";
				fields[i] *= 1.0e6;
			}
			sline >> fields[nFields - 1];
			std::cout << fields[nFields - 1] << std::endl;
			fields[nFields - 1] *= 1.0e6;
			fieldFound = 1;
		}
		else if (!keyword.compare("kmatrix"))
		{
			sline >> kMethod;
			transform(kMethod.begin(), kMethod.end(), kMethod.begin(), ::tolower);
			if (!kMethod.compare("hessian"))
			{
				sline >> kFile1;// >> kFile2;
				std::cout  << "- Using energy Hessian for force constants. Hessian from file: " << kFile1 << std::endl; //<< "; gradients from file " << kFile2 << std::endl;
			}
			else if (!kMethod.compare("covar"))
			{
				sline >> kFile1;
				std::cout << "- Using covariance matrix for force constants. Covariance matrix in Cartesian coordinates from file: " << kFile1 << std::endl;
			}
			else
			{
				std::cout << "ERROR: method " << kMethod << " for K matrix calculation not recognized. Available options are: hessian | covar" << std::endl;
				exit(1);
			}
			kMethodFound = 1;
		}
		else if (!keyword.compare("ignore_negative_hess_eigs"))
		{
			sline >> ignoreNegEigs;
			ignoreNegEigsFound = 1;
			if (ignoreNegEigs)
				std::cout << "- Sign of negative eigs of the Hessian will be ingnored. USE THIS OPTION WITH CAUTION AFTER HAVING ANALYZED THE HESSIAN THOROUGHLY." << std::endl;
			
		}
		else if (!keyword.compare("orderparam"))
		{
			calculateOrderParameters = 1;
			sline >> nSites;
			std::cout << "- Order parameters for " << nSites << " sites will be calculated:" << std::endl;
			Hid = vectorOfIntegers(nSites, 0);
			Nid = vectorOfIntegers(nSites, 0);
			Cid = vectorOfIntegers(nSites, 0);
			for (int i = 0; i < nSites; ++i)
			{
				sline >> Hid[i] >> Nid[i] >> Cid[i];
				std::cout << "   site " << i+1 << ": H(" << Hid[i] << "), N(" << Nid[i] << "), C(" << Cid[i] << ")" << std::endl;
			}
		}
		else if (!keyword.compare("activeatoms"))
		{
			sline >> itmp;
			if (itmp)
				allActive = 0;
			else
				allActive = 1;
			allActiveFound = 1;
			if (allActive)
				std::cout << "- All atoms will be included in the model" << std::endl;
			else
				std::cout << "- Only active atoms (with occupancy factor = 1) will be included in the model" << std::endl;
		}
		else
			std::cout << "WARNING: Ignoring unknown keyword: " << keyword << std::endl;
	}
	std::cout << "**********************" << std::endl;
	std::cout << "* End of input       *" << std::endl;
	std::cout << "**********************" << std::endl;

	checkInputConsistency();

	return;
}

void inputParserSrb::checkInputConsistency()
{
	int nerr = 0;
	if (!pdbFileNameFound)
	{
		std::cout << "ERROR: a PDB file must be specified via the ''pdb'' keyword in the input file." << std::endl;
		nerr ++;
	}
	if (!idsFound)
	{
		std::cout << "ERROR: the 4 reference atoms must be specified using the ''refAtoms ID1 ID2 ID3 ID4'' keyword." << std::endl;
		nerr ++;
	}
	if (!ReffFound)
	{
		std::cout << "WARNING: no Reff specified. Using default value of 2.0 Angstroms." << std::endl;
	}
	if (!CFound)
	{
		std::cout << "WARNING: no C specified. Uding default value of 6 (stick boundary conditions)." << std::endl;
	}
	if (!viscosityFound)
	{
		std::cout << "ERROR: viscosity (in Pa s) must be specified via the keyword ''viscosity''." << std::endl;
		nerr ++;
	}
	if (!temperatureFound)
	{
		std::cout << "WARNING: no temperature specified. Using default value of 298.15 K." << std::endl;
	}
	if (!hydrointFound)
	{
		std::cout << "WARNING: no hydrodynamics interactions specified. Using default value of 0; that is hydrodynamic interactions are not taken into account" << std::endl;
	}
//	if (!nfreqFound && !nmrFound)
//	{
//		std::cout << "ERROR: number of spectral density evaluations must be specified via the keyword ''nfreq''." << std::endl;
//		nerr++;
//	}
//	if (!deltafreqFound && !nmrFound)
//	{
//		std::cout << "ERROR: frequency interval (in fs^-1) must be specified for spectral density evaluations via the keyword ''deltafreq''." << std::endl;
//		nerr++;
//	}
//	if (!rankFound && !nmrFound)
//	{
//		std::cout << "ERROR: rank of the physical observable must be specified via the keyword ''rank''." << std::endl;
//		nerr++;
//	}
	if (rigidBodyFound && !drrFound)
	{
		std::cout << "ERROR: for rigid body calculation the principal values of the diffusion tensor must be given in input using the keyword ''drr''." << std::endl;
		nerr++;
	}
//	if (nmrFound && !dcsaFound)
//	{
//		std::cout << "ERROR: in NMR calculation, Delta CSA must be given in input using the keyword ''deltacsa''." << std::endl;
//		nerr++;
//	}
//	if (nmrFound && !probeFound)
//	{
//		std::cout << "ERROR: in NMR calculation, the probe must be given in input using the keyword ''probe''." << std::endl;
//		nerr++;
//	}
//	if (nmrFound && !probe.compare("ch2") && !hatomFound)
//	{
//		std::cout << "ERROR: in NMR calculation of CH2 probe, the ID of the second H atom must be set using the keyword ''Hatom''." << std::endl;
//		nerr++;
//	}
//	if (nmrFound && !expFound)
//	{
//		std::cout << "WARNING: no experimental file specified." << std::endl;
//	}
//	if (nmrFound && !fieldFound)
//	{
//		std::cout << "ERROR: in NMR calculation, at least 1 spectrometer frequency must be given in input using the keyword ''freq''." << std::endl;
//		nerr++;
//	}
	if (!kMethodFound)
	{
		std::cout << "ERROR: method for K matrix calculation must be specified using the keyword ''kmatrix''." << std::endl;
		nerr++;
	}
	if (!allActiveFound)
	{
		std::cout << "WARNING: ''activeAtoms'' keyword not found; including all the atoms in the model." << std::endl;
	}

	if (nerr)
	{
		std::cout << std::endl << std::endl << "Terminating SALEM because of errors." << std::endl << std::endl;
		exit(1);
	}
	return;
}

int inputParserSrb::ignoreNegativeEigs(void)
{
	return ignoreNegEigs;
}

