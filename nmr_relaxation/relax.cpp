#include <cstdlib>
#include <cstdio>

#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

#include <string>

#include "types.h"

#include "EulerQuaternion.h"
#include "acfDriver.h"
#include "nmr.h"
#include "fit_cf.h"

#undef BUILD_DIPOLAR_CSA_CCF

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		std::cout << std::endl;
		std::cout << "Usage: nmr data_file";
		std::cout << std::endl;
		std::cout << std::endl;
		exit(1);
	}

	nmr relax;
	relax.setDataFile((std::string)argv[1]);
	relax.readData();

	int probeCode, nTimes, nRep, nAtoms, nExp;
	double timeAlpha;
	double *omegaDip1, *omegaCsa1, *omegaDip2;
	cdvector *D200Dip1, *D200Csa1, *D200Dip2;
	std::string trjFile;

	////////////////
	// Trajectory //
	////////////////

	std::cout << "Number of atoms: ";
	std::cin  >> nAtoms;

	std::cout << "Number of trajectory snapshots: ";
	std::cin  >> nTimes;

	std::cout << "Number of chunks into which the trajectory will be cut: ";
	std::cin  >> nRep;

	std::cout << "Cartesian trajectory file: ";
	std::cin  >> trjFile;

	std::cout << "Number of exponential functions to fit correlation functions [1 to 5, suggested 4]: ";
	std::cin  >> nExp;

	std::cout << "Time stretching coefficient: ";
	std::cin  >> timeAlpha;

	////////////////////////////////////
	// Probe type and reference atoms //
	////////////////////////////////////

	std::cout << "Probe (1 - NH, 2 - CH, 3 - CH2): ";
	std::cin  >> probeCode;

	int nProbes;
	std::cout << "N probes = ";
	std::cin  >> nProbes;

	int *id1, *id2, *id3;
	id1 = new int[nProbes];
	id2 = new int[nProbes];
	id3 = new int[nProbes];

	D200Dip1 = new cdvector[nProbes];
	D200Csa1 = new cdvector[nProbes];
	D200Dip2 = new cdvector[nProbes];

	for (int np = 0; np < nProbes; ++np)
	{
		D200Dip1[np] = cdvector(nTimes, std::complex<double>(0.0, 0.0));
		D200Dip2[np] = cdvector(nTimes, std::complex<double>(0.0, 0.0));
		D200Csa1[np] = cdvector(nTimes, std::complex<double>(0.0, 0.0));

		switch (probeCode)
		{
			case 1:
			{
				std::cout << "ID's of atoms H, N, X (an atom bonded to N) of probe " << np+1 << ": ";
				break;
			}
			case 2:
			{
				std::cout << "ID's of atoms H, C, X (an atom bonded to C) of probe " << np+1 << ": ";
				break;
			}
			case 3:
			{
				std::cout << "ID's of atoms H, C, H of probe " << np+1 << ": ";
				break;
			}
			default:
			{
				std::cout << "Invalid probe code";
				return -1;
			}
		}	
	
		std::cin  >> id1[np] >> id2[np] >> id3[np];
	}

	///////////////////////////////////////////////
	// Read the time series of the Euler angles //
	///////////////////////////////////////////////

	int snapshot = 0;
	double dtmp;
	double *time = new double[nTimes];
	std::string lineDip1, lineCsa1, lineDip2;

	std::ifstream fileDip1("6_D200Dip1.dat");
	if (!fileDip1.is_open()) {
		std::cerr << "Error: Could not open the file 6_D200Dip1.dat" << std::endl;
		exit(-1);
	}

	std::ifstream fileCsa1("6_D200Csa1.dat");
	if (!fileCsa1.is_open()) {
		std::cerr << "Error: Could not open the file 6_D200Csa1.dat" << std::endl;
		exit(-1);
	}

	std::ifstream fileDip2("6_D200Dip2.dat");
	if (!fileDip2.is_open()) {
		std::cerr << "Error: Could not open the file 6_D200Dip2.dat" << std::endl;
		exit(-1);
	}

	std::cout << std::endl << std::endl << "Reading D200 time series..." <<std::endl;
	int grain = nTimes / 10;
	while (getline(fileDip1, lineDip1))
	{
		getline(fileCsa1, lineCsa1);
		getline(fileDip2, lineDip2);

		const char* p_dip1 = lineDip1.c_str();
	        const char* p_csa1 = lineCsa1.c_str();
        	const char* p_dip2 = lineDip2.c_str();

		int chars_read;

		sscanf(p_dip1, "%lf%n", &time[snapshot], &chars_read);
	        p_dip1 += chars_read;

		sscanf(p_csa1, "%lf%n", &dtmp, &chars_read);
	        p_csa1 += chars_read;

	        sscanf(p_dip2, "%lf%n", &dtmp, &chars_read);
	        p_dip2 += chars_read;

		for (int np = 0; np < nProbes; ++np)
		{
			sscanf(p_dip1, "%lf%n", &dtmp, &chars_read);
			D200Dip1[np][snapshot].real(dtmp);
			p_dip1 += chars_read;

			sscanf(p_csa1, "%lf%n", &dtmp, &chars_read);
			D200Csa1[np][snapshot].real(dtmp);
			p_csa1 += chars_read;

			sscanf(p_dip2, "%lf%n", &dtmp, &chars_read);
			D200Dip2[np][snapshot].real(dtmp);
			p_dip2 += chars_read;

		}
		snapshot++;
		if (!(snapshot%grain))
			std::cout << "..." << snapshot << " / " << nTimes << std::endl;
	}

	fileDip1.close();
	fileCsa1.close();
	fileDip2.close();

	////////////////////////
	// Spectral densities //
	////////////////////////

	int nData = nTimes / nRep;
	int ncorr = nData / 2;

	double norm, dataDt;
	dataDt = (time[1] - time[0]) * 1.0e-3; // dt in ns
	dataDt *= timeAlpha; // Apply time strecthing, which is equivalent to scale the diffusion tensor by a constant 1/timeAlpha
	acfDriver ad;
	ad.setTimeData(0.0, dataDt);

	dvector  gDip1Dip1;
	dvector  gCsa1Csa1;
	cdvector gDip1Csa1;

	dvector g = dvector(ncorr, 0.0);
	cdvector cg = cdvector(ncorr, std::complex<double>(0.0, 0.0));

	std::string results;

	std::string fileName;

	for (int np = 0; np < nProbes; ++np)
	{
		gDip1Dip1 = dvector(ncorr, 0.0);
		gCsa1Csa1 = dvector(ncorr, 0.0);
		gDip1Csa1 = cdvector(ncorr, std::complex<double>(0.0, 0.0));

		for (int p = 0; p < nRep; ++p)
		{
			////////////////////////////////////////////
			// Dip1 and Csa1 are common to all probes //
			////////////////////////////////////////////

			// Dip1-Dip1
			ad.setData({D200Dip1[np].begin() + p * nData, D200Dip1[np].begin() + (p + 1) * nData});
			ad.calculateACF(ncorr);
			g = ad.getACF();
			for (int i = 0; i < g.size(); ++i)
				gDip1Dip1[i] += g[i];
		
			// Csa1-Csa1
			ad.setData({D200Csa1[np].begin() + p * nData, D200Csa1[np].begin() + (p + 1) * nData});
			ad.calculateACF(ncorr);
			g = ad.getACF();
			for (int i = 0; i < g.size(); ++i)
				gCsa1Csa1[i] += g[i];
		
			// Dip1-Csa1
#ifdef BUILD_DIPOLAR_CSA_CCF
			ad.setData({D200Dip1[np].begin() + p * nData, D200Dip1[np].begin() + (p + 1) * nData}, {D200Csa1[np].begin() + p * nData, D200Csa1[np].begin() + (p + 1) * nData});
			ad.calculateCCF(ncorr);
			cg = ad.getCCF();
			for (int i = 0; i < cg.size(); ++i)
				gDip1Csa1[i] += cg[i];
#endif
		}

		// Fit correlation functions and compute Fourier transform //

		dvector tData = dvector(gDip1Dip1.size(), 0.0);

		// J Dip1-Dip1
		std::cout << std::endl << "** Calculating Jdip1-dip1 for probe " << np + 1 << std::endl << std::endl;

		double g0 =  1.0 / gDip1Dip1[0];
		for (int i = 0; i < gDip1Dip1.size(); ++i)
		{
			gDip1Dip1[i] *= g0;
			tData[i] = (double)i * dataDt;
		}
		ad.setCorrFun(gDip1Dip1);

		fileName = "6_Dip1_acf_pre-fit_probe-" + std::to_string(np + 1) + ".dat";
		ad.saveACF(fileName.c_str());
		std::cout << "- Fitting Dip1-Dip1 ACF with " << nExp << " exponential(s)..." << std::endl;
		gDip1Dip1 = fit_cf(tData, gDip1Dip1, nExp);
		ad.setCorrFun(gDip1Dip1);
		fileName = "6_Dip1_acf_probe-" + std::to_string(np + 1) + ".dat";
		ad.saveACF(fileName.c_str());
		ad.saveACF(fileName.c_str());
		for (int i = 0; i < 10240; ++i)
			gDip1Dip1.push_back(0.0);
		ad.setCorrFun(gDip1Dip1);
		ad.calculateRealSpectralDensity();
		dvector JDip1 = ad.getRealSpectralDensity();
		
		// J Csa1-Csa1
		std::cout << std::endl << "** Calculating Jcsa1-csa1 for probe " << np + 1 << std::endl << std::endl;

		g0 =  1.0 / gCsa1Csa1[0];
		for (int i = 0; i < gCsa1Csa1.size(); ++i)
			gCsa1Csa1[i] *= g0;
		ad.setCorrFun(gCsa1Csa1);
		fileName = "6_Csa1_acf_pre-fit_probe-" + std::to_string(np + 1) + ".dat";
		ad.saveACF(fileName.c_str());
		std::cout << "- Fitting Csa1-Csa1 ACF with " << nExp << " exponential(s)..." << std::endl;
		gCsa1Csa1 = fit_cf(tData, gCsa1Csa1, nExp);
		ad.setCorrFun(gCsa1Csa1);
		fileName = "6_Csa1_acf_probe-" + std::to_string(np + 1) + ".dat";
		ad.saveACF(fileName.c_str());
		for (int i = 0; i < 10240; ++i)
			gCsa1Csa1.push_back(0.0);
		ad.setCorrFun(gCsa1Csa1);
		ad.calculateRealSpectralDensity();
		dvector JCsa1 = ad.getRealSpectralDensity();

		// J Dip1-Csa1
#ifdef BUILD_DIPOLAR_CSA_CCF
		std::cout << std::endl << "** Calculating Jdip1-csa1 for probe " << np + 1 << std::endl << std::endl;

		dvector gDip1Csa1Real = dvector(nData, 0.0);
		g0 = 1.0 / gDip1Csa1[0].real();
		for (int i = 0; i < gDip1Csa1.size(); ++i)
			gDip1Csa1Real[i] = gDip1Csa1[i].real() * g0; // since we know that this is a real function
		ad.setCorrFun(gDip1Csa1Real);
		fileName = "6_Dip1-Csa1_ccf_pre-fit_probe-" + std::to_string(np + 1) + ".dat";
		ad.saveACF(fileName.c_str());
		std::cout << "- Fitting Dip1-Csa1 CCF (which is real) with " << nExp << " exponential(s)..." << std::endl;
		gDip1Csa1Real = fit_cf(tData, gDip1Csa1Real, nExp);
		ad.setCorrFun(gDip1Csa1Real);
		fileName = "6_Dip1-Csa1_ccf_probe-" + std::to_string(np + 1) + ".dat";
		ad.saveACF(fileName.c_str());
		for (int i = 0; i < gDip1Csa1.size(); ++i)
			gDip1Csa1[i] = std::complex<double>(gDip1Csa1Real[i], 0.0);
		for (int i = 0; i < 10240; ++i)
			gDip1Csa1.push_back(std::complex<double>(0.0, 0.0));
		ad.setCrossCorrFun(gDip1Csa1);
		ad.calculateRealCrossSpectralDensity();
		dvector JcrossD1C = ad.getRealSpectralDensity();
#else
		dvector JcrossD1C = dvector(ncorr, 0.0);
#endif

		// Write j.dat (freq in 1/s, J in s)
		fileName = "6_j_probe-" + std::to_string(np + 1) + ".dat";
		std::ofstream jFile(fileName.c_str());

		dvector f = ad.getFrequencies();

		switch (probeCode)
		{
			case 1:
			case 2:
			{
				for (int i = 0; i < JDip1.size()/2; ++i)
					jFile << std::scientific << std::setprecision(8) << f[i] << " " << JDip1[i] * 1.0e-9 << " " << JCsa1[i] * 1.0e-9 << " " << JcrossD1C[i] * 1.0e-9 << std::endl;
				jFile.close();
				break;
			}
			case 3:  // TODO: NOTMALIZATION OF CORRELATION FUNCTIONS
			{
				dvector  gDip2Dip2 = dvector(ncorr, 0.0);
				cdvector gDip2Csa1 = cdvector(ncorr, std::complex<double>(0.0, 0.0));
				cdvector gDip1Dip2 = cdvector(ncorr, std::complex<double>(0.0, 0.0));

				for (int p = 0; p < nRep; ++p)
				{
					////////////////////////////////////
					// Dip2, Dip1-Dip2, and Dip2-Csa1 //
					////////////////////////////////////

					// Dip2-Dip2
					ad.setData({D200Dip2[np].begin() + p * nData, D200Dip2[np].begin() + (p + 1) * nData});
					ad.calculateACF(ncorr);
					g = ad.getACF();
					for (int i = 0; i < g.size(); ++i)
						gDip2Dip2[i] += g[i];
				
					// Dip1-Dip2
					ad.setData({D200Dip1[np].begin() + p * nData, D200Dip1[np].begin() + (p + 1) * nData}, {D200Dip2[np].begin() + p * nData, D200Dip2[np].begin() + (p + 1) * nData});
					ad.calculateCCF(ncorr);
					cg = ad.getCCF();
					for (int i = 0; i < cg.size(); ++i)
						gDip1Dip2[i] += cg[i];

					// Dip2-Csa1
					ad.setData({D200Dip2[np].begin() + p * nData, D200Dip2[np].begin() + (p + 1) * nData}, {D200Csa1[np].begin() + p * nData, D200Csa1[np].begin() + (p + 1) * nData});
					ad.calculateCCF(ncorr);
					cg = ad.getCCF();
					for (int i = 0; i < cg.size(); ++i)
						gDip2Csa1[i] += cg[i];
				
				}

				// J Dip2-Dip2
				std::cout << std::endl << "** Calculating Jdip2-dip2 for probe " << np + 1 << std::endl << std::endl;

				for (int i = 0; i < gDip2Dip2.size(); ++i)
				{
					gDip2Dip2[i] /= (double)nRep;
					tData[i] = (double)i * dataDt;
				}
				ad.setCorrFun(gDip2Dip2);
				fileName = "6_Dip2_acf_pre-fit_probe-" + std::to_string(np + 1) + ".dat";
				ad.saveACF(fileName.c_str());
				std::cout << "- Fitting Dip2-Dip2 ACF with " << nExp << " exponential(s)..." << std::endl;
				gDip2Dip2 = fit_cf(tData, gDip2Dip2, nExp);
				ad.setCorrFun(gDip2Dip2);
				fileName = "6_Dip2_acf_probe-" + std::to_string(np + 1) + ".dat";
				ad.saveACF(fileName.c_str());
				for (int i = 0; i < 10240; ++i)
					gDip2Dip2.push_back(0.0);
				ad.setCorrFun(gDip2Dip2);
				ad.calculateRealSpectralDensity();
				dvector JDip2 = ad.getRealSpectralDensity();
		
				// J Dip1-Dip2
				std::cout << std::endl << "** Calculating Jdip1-dip2 for probe " << np + 1 << std::endl << std::endl;

				dvector gDip1Dip2Real = dvector(nData, 0.0);
				for (int i = 0; i < gDip1Dip2.size(); ++i)
				{
					gDip1Dip2[i] /= std::complex<double>((double)nRep, 0.0);
					gDip1Dip2Real[i] = gDip1Dip2[i].real(); // since we know that this is a real function
				}
				ad.setCorrFun(gDip1Dip2Real);
				fileName = "6_Dip1-Dip2_ccf_pre-fit_probe-" + std::to_string(np + 1) + ".dat";
				ad.saveACF(fileName.c_str());
				std::cout << "- Fitting Dip1-Dip2 CCF (which is real) with " << nExp << " exponential(s)..." << std::endl;
				gDip1Dip2Real = fit_cf(tData, gDip1Dip2Real, nExp);
				ad.setCorrFun(gDip1Dip2Real);
				fileName = "6_Dip1-Dip2_ccf_probe-" + std::to_string(np + 1) + ".dat";
				ad.saveACF(fileName.c_str());
				for (int i = 0; i < gDip1Dip2.size(); ++i)
					gDip1Dip2[i] = std::complex<double>(gDip1Dip2Real[i], 0.0);
				ad.setCrossCorrFun(gDip1Dip2);
				for (int i = 0; i < 10240; ++i)
					gDip1Dip2.push_back(std::complex<double>(0.0, 0.0));
				ad.setCrossCorrFun(gDip1Dip2);
				ad.calculateRealCrossSpectralDensity();
				dvector JcrossD1D2 = ad.getRealSpectralDensity();

				// J Dip2-Csa1
				std::cout << std::endl << "** Calculating Jdip2-dcsa1 for probe " << np + 1 << std::endl << std::endl;

				dvector gDip2Csa1Real = dvector(nData, 0.0);
				for (int i = 0; i < gDip2Csa1.size(); ++i)
				{
					gDip2Csa1[i] /= std::complex<double>((double)nRep, 0.0);
					gDip2Csa1Real[i] = gDip2Csa1[i].real(); // since we know that this is a real function
				}
				ad.setCorrFun(gDip2Csa1Real);
				fileName = "6_Dip2-Csa1_ccf_pre-fit_probe-" + std::to_string(np + 1) + ".dat";
				ad.saveACF(fileName.c_str());
				std::cout << "- Fitting Dip2-Csa1 CCF (which is real) with " << nExp << " exponential(s)..." << std::endl;
				gDip2Csa1Real = fit_cf(tData, gDip2Csa1Real, nExp);
				ad.setCorrFun(gDip2Csa1Real);
				fileName = "6_Dip2-Csa1_ccf_probe-" + std::to_string(np + 1) + ".dat";
				ad.saveACF(fileName.c_str());
				for (int i = 0; i < gDip2Csa1.size(); ++i)
					gDip2Csa1[i] = std::complex<double>(gDip2Csa1Real[i], 0.0);
				ad.setCrossCorrFun(gDip2Csa1);
				for (int i = 0; i < 10240; ++i)
					gDip2Csa1.push_back(std::complex<double>(0.0, 0.0));
				ad.setCrossCorrFun(gDip2Csa1);
				ad.calculateRealCrossSpectralDensity();
				dvector JcrossD2C = ad.getRealSpectralDensity();


				// Write j.dat file

				dvector f = ad.getFrequencies();

				fileName = "6_j_probe-" + std::to_string(np + 1) + ".dat";
				std::ofstream jFile(fileName.c_str());
				for (int i = 0; i < JDip1.size()/2; ++i)
				{
					jFile << std::scientific << std::setprecision(8) << f[i] << " ";
					jFile << JDip1[i] * 1.0e-9 << " " << JDip2[i] * 1.0e-9 << " " << JCsa1[i] * 1.0e-9 << " ";
					jFile << JcrossD1D2[i] * 1.0e-9 << " " << JcrossD1C[i] * 1.0e-9 << " " << JcrossD2C[i] * 1.0e-9 << std::endl;
				}
				jFile.close();

				break;
			}
		}

		fileName = "6_j_probe-" + std::to_string(np + 1) + ".dat";
		relax.setSpdFile(fileName.c_str());
		relax.readSpd();

//		double dCSA = 34., delta = 2.0;
//		for (int i = 0; i <= 20; ++i)
//		{
//			relax.setDeltaCsa(dCSA - (double)i * delta);
			relax.calcNMR();
			results = relax.toString();

			std::cout << std::endl << std::endl;
			std::cout << "****************************************" << std::endl;
			std::cout << "*** NMR RELAXATION RATES OF PROBE " << np + 1 << " ***" << std::endl;
			std::cout << "****************************************" << std::endl;
			std::cout << results << std::endl;
//		}

	}

	return 0;
}
