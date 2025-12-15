//////////////////////////////////////
//             NMR Module           //
//                                  //
// NB: ALL QUANTITIES IN S.I. UNITS //
//////////////////////////////////////
#include "nmr.h"

/***************/
/* Constructor */
/***************/
nmr::nmr()
{
	probeSet = 0;
	fieldSet = 0;
	expFile  = "nofile";
	buildDip2 = 0;

	// FOR THE MOMENT Rex IS HARD CODED TO 0
	rex = 0.0;

	return;
}

/*****************/
/* Deconstructor */
/*****************/
nmr::~nmr()
{
	return;
}

/***********/
/* METHODS */
/***********/
void nmr::setProbe(std::string p)
{
	probe = p;
	transform(probe.begin(), probe.end(), probe.begin(), ::tolower);
	if (!probe.compare("nh"))      // For 15-N NMR
	{
		gyroX = -2.7116e7;
		rm3_ov = 1.0 / (1.015 * 1.015 * 1.015);
		probeSet = 1;
	}
	else if (!probe.compare("ch"))      // For 13-C NMR
	{
		gyroX = 6.72828e7;
		rm3_ov = 1.0 / (1.13 * 1.13 * 1.13);
		probeSet = 1;
	}
	else if (!probe.compare("ch2"))      // For 13-C NMR
	{
		gyroX = 6.72828e7;
		rm3_ov = 1.0 / (1.13 * 1.13 * 1.13);
		probeSet = 1;
		buildDip2 = 1;
	}
	else
	{
		std::cout << std::endl << std::endl << "ERROR: probe " << p << " not implemented/recognized." << std::endl << std::endl;
		exit(1);
	}
	return;
}

// Setup the NMR calculation

void nmr::setup(molecule mol, inputParserSrb ip)
{
	/* 1. Fields */

	int nf = ip.getNFields();
	double *fields = new double[nf];
	ip.getFields(fields);


	/* 2. Setup probe */

	setProbe(ip.getProbe());
	setFields(fields, nf);

	/* 3. Spin data */

	setDeltaCSA(ip.getDeltaCSA());

	/* 4. Geometry */

	// 4.a r_XH^3
	setAveragedRm3(rm3_ov * 1.0e30); // nmr module works in SI // Fixed X-H length, hard-coded

	// 4.b Omega_Dip.->CSA
	setOmegaDC(0.0, -17.0 * M_PI / 180.0, 0.0); // This is hard-coded for proteins

	// 4.c Omega_Dip1->Dip2	double vnorm;

	if (buildDip2)
	{
		double vnorm;
		matrix3D E_AD2, E_AD1_tr, E_D1D2;

		atom H2 = mol.getAtomOld(ip.getHatom());
		vectorOfDoubles rH2 = H2.getPosition();

		vectorOfDoubles rH1 = mol.getAtomPosition(1);
		vectorOfDoubles rC  = mol.getAtomPosition(2);
		vectorOfDoubles rA  = mol.getAtomPosition(3);

		vectorOfDoubles xD2(3), yD2(3), zD2(3), v(3);

		v[0] = rA[0] - rC[0]; v[1] = rA[1] - rC[1]; v[2] = rA[2] - rC[2];

		zD2[0] = rH2[0] - rC[0]; zD2[1] = rH2[1] - rC[1]; zD2[2] = rH2[2] - rC[2];
		vnorm = 1.0 / sqrt(zD2[0] * zD2[0] + zD2[1] * zD2[1] + zD2[2] * zD2[2]);
		zD2[0] *= vnorm; zD2[1] *= vnorm; zD2[2] *= vnorm;

		xD2[0] = v[1] * zD2[2] - v[2] * zD2[1];
		xD2[1] = v[2] * zD2[0] - v[0] * zD2[2];
		xD2[2] = v[0] * zD2[1] - v[1] * zD2[0];
		vnorm = 1.0 / sqrt(xD2[0] * xD2[0] + xD2[1] * xD2[1] + xD2[2] * xD2[2]);
		xD2[0] *= vnorm; xD2[1] *= vnorm; xD2[2] *= vnorm;

		yD2[0] = zD2[1] * xD2[2] - zD2[2] * xD2[1];
		yD2[1] = zD2[2] * xD2[0] - zD2[0] * xD2[2];
		yD2[2] = zD2[0] * xD2[1] - zD2[1] * xD2[0];

		E_AD2.xx = xD2[0]; E_AD2.xy = xD2[1]; E_AD2.xz = xD2[2];
		E_AD2.yx = yD2[0]; E_AD2.yy = yD2[1]; E_AD2.yz = yD2[2];
		E_AD2.zx = zD2[0]; E_AD2.zy = zD2[1]; E_AD2.zz = zD2[2];

		E_AD1_tr.xx =  0.0; E_AD1_tr.xy =  0.0; E_AD1_tr.xz = 1.0;
		E_AD1_tr.yx =  0.0; E_AD1_tr.yy =  1.0; E_AD1_tr.yz = 0.0;
		E_AD1_tr.zx = -1.0; E_AD1_tr.zy =  0.0; E_AD1_tr.zz = 0.0;

		E_D1D2.multiply(E_AD2, E_AD1_tr);
		double alphaDD = atan2(E_D1D2.zy, E_D1D2.zx);
		double betaDD  = acos(E_D1D2.zz);
		double gammaDD = atan2(E_D1D2.yz, -E_D1D2.xz);
	}

	/* 5. File with experimental values */
	setExpFile(ip.getExpFile());

	return;
}

// GET/SET properties

std::string nmr::getProbe(void)
{
	return probe;
}

void nmr::setDeltaCSA(double d)
{
	deltaCSA = d;
	return;
}

double nmr::getDeltaCSA(void)
{
	return deltaCSA;
}

void nmr::setOmegaDC(double a, double b, double g)
{
	alphaDC = a;
	betaDC  = b;
	gammaDC = g;
	return;
}

void nmr::getOmegaDC(double* eul)
{
	eul[0] = alphaDC;
	eul[1] = betaDC;
	eul[2] = gammaDC;
	return;
}

void nmr::setOmegaDD(double a, double b, double g)
{
	alphaDD = a;
	betaDD  = b;
	gammaDD = g;
	return;
}

void nmr::getOmegaDD(double* eul)
{
	eul[0] = alphaDD;
	eul[1] = betaDD;
	eul[2] = gammaDD;
	return;
}

void nmr::setOmegaD(double a, double b, double g)
{
	alphaD = a;
	betaD  = b;
	gammaD = g;
	return;	
}

void nmr::getOmegaD(double *eul)
{
	eul[0] = alphaD;
	eul[1] = betaD;
	eul[2] = gammaD;
	return;
}

int nmr::getNLarmorFreq(void)
{
	return nLarmorFreq;
}

double nmr::getLarmorFreq(int n)
{
	if (n >= 0 && n < nLarmorFreq)
		return larmorFreq[n];
	else
	{
		std::cout << std::endl << std::endl << "ERROR in " << __FILE__ << ", line " << __LINE__ << ": index of Larmor freq. array out of boundary. Should be >= 0 and < " << nLarmorFreq << ", while the actual value is " << n << std::endl << std::endl;
		exit(1);
	}
	return 0.0;
}

// THE ARRAY OF SMALL j's HAS DIMENSION nLarmorFreq x 9
// FOR EACH FREQUENCY, THE 9 jKK' ARE ORDERED AS
// (0,0), (1,1), (2,2), (-2,2), (-1,1), (-1,2), (0,1), (0,2), (1,2)
void nmr::setSmallJs(double *js, int n)
{
	if (n >= 0 && n < nLarmorFreq)
	{
		for (int i = 0; i < 9; ++i)
		{
			smalljs[n * 9 + i] = js[i];
		}
		setSmallJCSA(n);
		if (!probe.compare("ch2"))
			setSmallJDip2(n);
	}
	else
	{
		std::cout << std::endl << std::endl << "ERROR in " << __FILE__ << ", line " << __LINE__ << ": index of Larmor freq. array out of boundary. Should be >= 0 and < " << nLarmorFreq << ", while the actual value is " << n << std::endl << std::endl;
		exit(1);
	}
	return;
}

void nmr::setSmallJCSA(int n)
{
	double tmpsmalljs[25];
	dcomplex dummy;
	dcomplex tmp;
	int cnt;
	 tmpsmalljs[0] = smalljs[n*9 + 2];  tmpsmalljs[1] = smalljs[n*9 + 8];  tmpsmalljs[2] = smalljs[n*9 + 7];  tmpsmalljs[3] = smalljs[n*9 + 5];  tmpsmalljs[4] = smalljs[n*9 + 3];
	 tmpsmalljs[5] =-smalljs[n*9 + 8];  tmpsmalljs[6] = smalljs[n*9 + 1];  tmpsmalljs[7] = smalljs[n*9 + 6];  tmpsmalljs[8] = smalljs[n*9 + 4];  tmpsmalljs[9] = smalljs[n*9 + 5];
	tmpsmalljs[10] = smalljs[n*9 + 7]; tmpsmalljs[11] =-smalljs[n*9 + 6]; tmpsmalljs[12] = smalljs[n*9 + 0]; tmpsmalljs[13] = smalljs[n*9 + 6]; tmpsmalljs[14] = smalljs[n*9 + 7];
	tmpsmalljs[15] =-smalljs[n*9 + 5]; tmpsmalljs[16] = smalljs[n*9 + 4]; tmpsmalljs[17] =-smalljs[n*9 + 6]; tmpsmalljs[18] = smalljs[n*9 + 1]; tmpsmalljs[19] = smalljs[n*9 + 8];
	tmpsmalljs[20] = smalljs[n*9 + 3]; tmpsmalljs[21] =-smalljs[n*9 + 5]; tmpsmalljs[22] = smalljs[n*9 + 7]; tmpsmalljs[23] =-smalljs[n*9 + 8]; tmpsmalljs[24] = smalljs[n*9 + 2];

	for (int k = 0; k <= 2; k++)
	{
		cnt = 0;
		tmp.real(0.0); tmp.imag(0.0);
		for (int i = -2; i <= 2; i++)
		{
			dummy = D2MK(i, k, alphaDC, betaDC, gammaDC);
			for (int j = -2; j <= 2; j++)
			{
				tmp += dummy * D2MK(j, k, alphaDC, betaDC, gammaDC) * tmpsmalljs[cnt];
				cnt++;
			}
		}
		smalljCSA[n * 3 + k] = tmp.real();
	}

	return;
}

// Build the spectral densities required for the 2nd C-H probe in CH2 moieties
void nmr::setSmallJDip2(int n)
{
	double tmpsmalljs[25];
	dcomplex dummy;
	dcomplex tmp;
	int cnt;
	 tmpsmalljs[0] = smalljs[n*9 + 2];  tmpsmalljs[1] = smalljs[n*9 + 8];  tmpsmalljs[2] = smalljs[n*9 + 7];  tmpsmalljs[3] = smalljs[n*9 + 5];  tmpsmalljs[4] = smalljs[n*9 + 3];
	 tmpsmalljs[5] =-smalljs[n*9 + 8];  tmpsmalljs[6] = smalljs[n*9 + 1];  tmpsmalljs[7] = smalljs[n*9 + 6];  tmpsmalljs[8] = smalljs[n*9 + 4];  tmpsmalljs[9] = smalljs[n*9 + 5];
	tmpsmalljs[10] = smalljs[n*9 + 7]; tmpsmalljs[11] =-smalljs[n*9 + 6]; tmpsmalljs[12] = smalljs[n*9 + 0]; tmpsmalljs[13] = smalljs[n*9 + 6]; tmpsmalljs[14] = smalljs[n*9 + 7];
	tmpsmalljs[15] =-smalljs[n*9 + 5]; tmpsmalljs[16] = smalljs[n*9 + 4]; tmpsmalljs[17] =-smalljs[n*9 + 6]; tmpsmalljs[18] = smalljs[n*9 + 1]; tmpsmalljs[19] = smalljs[n*9 + 8];
	tmpsmalljs[20] = smalljs[n*9 + 3]; tmpsmalljs[21] =-smalljs[n*9 + 5]; tmpsmalljs[22] = smalljs[n*9 + 7]; tmpsmalljs[23] =-smalljs[n*9 + 8]; tmpsmalljs[24] = smalljs[n*9 + 2];

	for (int k = 0; k <= 2; k++)
	{
		cnt = 0;
		tmp.real(0.0); tmp.imag(0.0);
		for (int i = -2; i <= 2; i++)
		{
			dummy = D2MK(i, k, alphaDD, betaDD, gammaDD);
			for (int j = -2; j <= 2; j++)
			{
				tmp += dummy * D2MK(j, k, alphaDD, betaDD, gammaDD) * tmpsmalljs[cnt];
				cnt++;
			}
		}
		smalljDip2[n * 3 + k] = tmp.real();
	}

	return;
}

void nmr::setAveragedRm3(double rm3)
{
	averageRm3 = rm3;
	return;
}

double nmr::getAveragedRm3(void)
{
	return averageRm3;
}

double nmr::getRm3Overridden(void)
{
	return rm3_ov;
}

void nmr::setExpFile(std::string ef)
{

	if (!ef.compare("nofile"))
		return;

	expFile = ef;
	std::fstream f;
	
	expNmrData.clear();
	int n;
	double tmp;

	f.open(expFile.c_str(), std::ios::in);
	f >> n;

	for (int i = 0; i < 3 * n; ++i)
	{
		f >> tmp;
		expNmrData.push_back(tmp);
	}

	f.close();

	return;
}

std::string nmr::getExpFile(void)
{
	return expFile;
}

int nmr::getNFields(void)
{
	return nFields;
}

void nmr::setFields(double *f, int nf)
{
	nFields = nf;
	fields = new double[nFields];
	for (int i = 0; i < nFields; ++i)
		fields[i] = f[i];
	fieldSet = 1;
	makeLarmorFrequencies();
	return;
}

void nmr:: getFields(double *f)
{
	for (int i = 0; i < nFields; ++i)
		f[i] = fields[i];
	return;
}

void nmr::makeLarmorFrequencies(void)
{
	if (probeSet && fieldSet)
	{
		// Fill arrays of Larmor frequencies
		nLarmorFreq = nFields * 5;
		larmorFreq = new double[nLarmorFreq];
		freqH = new double[nFields];
		freqX = new double[nFields];
		for (int f = 0; f < nFields; ++f)
		{
			freqH[f] = fields[f] * 2.0 * M_PI;
			freqX[f] = gyroX * freqH[f] / gyroH;
			larmorFreq[f * 5 + 0]  = 0.0;
			larmorFreq[f * 5 + 1] =  freqH[f];
			if (!probe.compare("nh")) larmorFreq[f * 5 + 2] = -freqX[f];
			if (!probe.compare("ch") || !probe.compare("ch2")) larmorFreq[f * 5 + 2] = freqX[f];
			larmorFreq[f * 5 + 3] =  freqH[f] - freqX[f];
			larmorFreq[f * 5 + 4] =  freqH[f] + freqX[f];
		}

		// Init the arrays of small j's
		smalljs = new double[nLarmorFreq * 9]; // 9 are the independent couples K,K'
		smalljCSA = new double[nLarmorFreq * 3];
		if(buildDip2)
			smalljDip2 = new double[nLarmorFreq * 3];
	}
	return;
}

void nmr::calculateNMRdatainDF(void)
{
	nmrData.clear();

	double dipFactor, csaFactor;
	
	dipFactor  =  MU0_OVER_4PI * gyroH * gyroX * hbar * averageRm3;
	dipFactor *= dipFactor;
	dipFactor *= 0.1;
	
	double freqHmin = minFreq(freqH, nFields);

	double exchangeCorr, T1, T2, NOE;

// AP: 13-07-20 : SITEMARE CICLO SU nLarmorFreq = nFields * 5 // <-- ???
	for (int ifield = 0; ifield < nFields; ifield++)
	{

		csaFactor = freqX[ifield] * freqX[ifield] * deltaCSA * deltaCSA * TWO_OVER_FIFTHEEN * 1.0e-12; // last multiplication since deltaCSA is in ppm
	
	        exchangeCorr = (freqH[ifield]/freqHmin);
	        exchangeCorr *= exchangeCorr;
	        exchangeCorr *= rex;

		T1 = dipFactor * (smalljs[3*9 + 0] + 3.0 * smalljs[2*9 + 1] + 6.0 * smalljs[4*9 + 2]) + csaFactor * smalljCSA[2*3 + 1];
		if (!probe.compare("ch2"))
			T1 += dipFactor * (smalljDip2[3*3 + 0] + 3.0 * smalljDip2[2*3 + 1] + 6.0 * smalljDip2[4*3 + 2]);
		T1 = 1.0/T1;
		nmrData.push_back(T1);

		T2 = 0.5 * dipFactor * (4.0 * smalljs[0*9 + 0] + smalljs[3*9 + 0]+ 3.0 * smalljs[2*9 + 1] + 6.0 * smalljs[1*9 + 1] + 6.0 * smalljs[4*9 + 2]) + ONE_OVER_SIX * csaFactor * (3.0 * smalljCSA[2*3 + 1] + 4.0 * smalljCSA[0*3 + 0]) + exchangeCorr;
		if (!probe.compare("ch2"))
			T2 += 0.5 * dipFactor * (4.0 * smalljDip2[0*3 + 0] + smalljDip2[3*3 + 0]+ 3.0 * smalljDip2[2*3 + 1] + 6.0 * smalljDip2[1*3 + 1] + 6.0 * smalljDip2[4*3 + 2]);
		T2 = 1.0/T2;
		nmrData.push_back(T2);

		NOE = 1.0 + (gyroH / gyroX) * dipFactor * (6.0 * smalljs[4*9 + 2] - smalljs[3*9 + 0]) * T1;
		if (!probe.compare("ch2"))
			NOE += (gyroH / gyroX) * dipFactor * (6.0 * smalljDip2[4*3 + 2] - smalljDip2[3*3 + 0]) * T1;
		nmrData.push_back(NOE);

		// Dipolar cross-relaxations
		// TODO: Build Dip1-Dip2 cross spectral densitied bigKdip
		//       before uncommenting the following lines
		/*
		if (!probe.compare("ch2"))
		{
			double dipFactor2 = 10.0*dipFactor/freqScale;
			double CCRRT1 = 0.6 * dipFactor2*bigKdip[2];
			relaxationTimes.push_back(CCRRT1);
			double CCRRT2 = 0.3 * dipFactor2*((4.0/3.0)*bigKdip[0] + bigKdip[2]);
			relaxationTimes.push_back(CCRRT2);
		}
		*/
	
	}
	return;
}

double nmr::minFreq(double *f, int n)
{
	double m = f[0];
	for (int i = 1; i < n; ++i)
		m = (f[i] < m) ? f[i] : m;
	return m;
}

dcomplex nmr::D2MK(int m, int k, double a, double b, double g)
{
	double cb, sb;
	double smalld;
	int h = (m + 2) * 5 + (k + 2);
	switch (h)
	{
		case 0:
		case 24:
		{
			cb = cos(b);
			smalld  = 1.0 + cb;
			smalld *= smalld * 0.25;
			break;
		}
		case 1:
		case 5:
		case 19:
		case 23:
		{
			sb = sin(b);
			cb = cos(b);
			smalld = -0.5 * sb * (1.0 + cb);
			smalld *= (h == 1 || h == 19) ? -1.0 : 1.0;
			break;
		}
		case 2:
		case 10:
		case 14:
		case 22:
		{
			sb = sin(b);
			smalld = SQRT_THREE_OVER_EIGHT * sb * sb;
			break;
		}
		case 3:
		case 9:
		case 15:
		case 21:
		{
			sb = sin(b);
			cb = cos(b);
			smalld = -0.5 * sb * (1.0 - cb);
			smalld *= (h == 3 || h == 9) ? -1.0 : 1.0;
			break;
		}
		case 4:
		case 20:
		{
			cb = cos(b);
			smalld = 1.0 - cb;
			smalld *= smalld * 0.25;
			break;
		}
		case 6:
		case 18:
		{
			cb = cos(b);
			smalld = 0.5 * (2.0 * cb * cb + cb - 1.0);
			break;
		}
		case 7:
		case 11:
		case 13:
		case 17:
		{
			sb = sin(2.0 * b);
			smalld = -SQRT_THREE_OVER_EIGHT * sb;
			smalld *= (h ==7 || h ==13) ? -1.0 : 1.0;
			break;
		}
		case 8:
		case 16:
		{
			cb = cos(b);
			smalld = 0.5 * (-2.0 * cb * cb + cb + 1.0);
			break;
		}
		case 12:
		{
			cb = cos(b);
			smalld = 0.5 * (3.0 * cb * cb - 1.0);
			break;
		}
		default:
		{
			std::cout << std::endl << std::endl << "ERROR in " << __FILE__ << ", line " << __LINE__ << ": m = " << m << "or k = " << k << " out of range in calling D2MK." << std::endl << std::endl;
			exit(1);
		}
	}
	
	double ca = cos((double)m * a), sa = sin((double)m * a);
	double cg = cos((double)k * g), sg = sin((double)k * g);

	dcomplex ea = std::complex<double>(ca, sa);
	dcomplex eg = std::complex<double>(cg, sg);
	
	return (ea * smalld * eg);
}

dvector nmr::getNMRdata(void)
{
	return nmrData;
}

std::string nmr::toString(void)
{
	std::stringstream s;
	for (int i = 0; i < nFields; i++)
	{
		s << fields[i] << " Hz" << std::endl;

		s << "T1 / ms: " << 1000.0 * nmrData[i * 3 + 0];
		if (expFile.compare("nofile"))
		{
			s << "\t" << 1000.0 * expNmrData[i * 3 + 0] << "\t" << 100.0 * (expNmrData[i * 3 + 0] - nmrData[i * 3 + 0]) / expNmrData[i * 3 + 0];
		}
		s << std::endl;
	
		s << "T2 / ms: " << 1000.0 * nmrData[i * 3 + 1];
		if (expFile.compare("nofile"))
		{
			s << "\t" << 1000.0 * expNmrData[i * 3 + 1] << "\t" << 100.0 * (expNmrData[i * 3 + 1] - nmrData[i * 3 + 1]) / expNmrData[i * 3 + 1];
		}
		s << std::endl;
	
		s << "NOE: " << nmrData[i * 3 + 2];
		if (expFile.compare("nofile"))
		{
			s << "\t" << expNmrData[i * 3 + 2] << "\t" << 100.0 * (expNmrData[i * 3 + 2] - nmrData[i * 3 + 2]) / expNmrData[i * 3 + 2];
		}
		s << std::endl;
	}
	return s.str();
}

