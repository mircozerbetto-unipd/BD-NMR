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


void calculateEulerAngles(int t, int probe, double *rA1, double *rA2, double *rA3, double *omegaDip1, double *omegaCsa1, double *omegaDip2)
{
	
	double norm;
	double R[9];
	EulerAngles Omega;
	EulerQuaternion eq;

	////////////////////////////////////////////////////////
	// First dipolar frame - same for all kinds of probes //	
	////////////////////////////////////////////////////////

	// zD1
	R[6] = rA2[0] - rA1[0];
	R[7] = rA2[1] - rA1[1];
	R[8] = rA2[2] - rA1[2];
	norm = 1.0 / sqrt(R[6] * R[6] + R[7] * R[7] + R[8] * R[8]);
	R[6] *= norm;
	R[7] *= norm;
	R[8] *= norm;

	// supporting v
	R[0] = rA3[0] - rA2[0];
	R[1] = rA3[1] - rA2[1];
	R[2] = rA3[2] - rA2[2];
	norm = 1.0 / sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
	R[0] *= norm;
	R[1] *= norm;
	R[2] *= norm;

	// yD1 = v x zD1
	R[3] = R[1] * R[8] - R[2] * R[7];
	R[4] = R[2] * R[6] - R[0] * R[8];
	R[5] = R[0] * R[7] - R[1] * R[6];
	norm = 1.0 / sqrt(R[3] * R[3] + R[4] * R[4] + R[5] * R[5]);
	R[3] *= norm;
	R[4] *= norm;
	R[5] *= norm;

	// xD1 = yD1 x zD1
	R[0] = R[4] * R[8] - R[5] * R[7];
	R[1] = R[5] * R[6] - R[3] * R[8];
	R[2] = R[3] * R[7] - R[4] * R[6];
	norm = 1.0 / sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
	R[0] *= norm;
	R[1] *= norm;
	R[2] *= norm;

	eq.setRotationMatrix(R);
	Omega = eq.getEulerAngles();
	omegaDip1[t * 3 + 0] = Omega.alpha;
	omegaDip1[t * 3 + 1] = Omega.beta;
	omegaDip1[t * 3 + 2] = Omega.gamma;

	////////////////////////
	// CSA for all probes //
	////////////////////////
	
	// Here, CSA frame is fixed at (0, 17 deg, 0) from dipolar frame for NH, while no tilt for CHn
	if (probe == 1)
	{

		eq.setEulerAngles(0.0, 17.0 * M_PI / 180.0, 0.0);
		double *RDC = new double[9];
		eq.getRotationMatrix(RDC);

		double RC[9];
		
		RC[0] = RDC[0] * R[0] + RDC[1] * R[3] + RDC[2] * R[6];
		RC[1] = RDC[0] * R[1] + RDC[1] * R[4] + RDC[2] * R[7];
		RC[2] = RDC[0] * R[2] + RDC[1] * R[5] + RDC[2] * R[8];
		
		RC[3] = RDC[3] * R[0] + RDC[4] * R[3] + RDC[5] * R[6];
		RC[4] = RDC[3] * R[1] + RDC[4] * R[4] + RDC[5] * R[7];
		RC[5] = RDC[3] * R[2] + RDC[4] * R[5] + RDC[5] * R[8];
		
		RC[6] = RDC[6] * R[0] + RDC[7] * R[3] + RDC[8] * R[6];
		RC[7] = RDC[6] * R[1] + RDC[7] * R[4] + RDC[8] * R[7];
		RC[8] = RDC[6] * R[2] + RDC[7] * R[5] + RDC[8] * R[8];

		eq.setRotationMatrix(RC);
		Omega = eq.getEulerAngles();
		omegaCsa1[t * 3 + 0] = Omega.alpha;
		omegaCsa1[t * 3 + 1] = Omega.beta;
		omegaCsa1[t * 3 + 2] = Omega.gamma;
	}
	else
	{
		omegaCsa1[t * 3 + 0] = omegaDip1[t * 3 + 0];
		omegaCsa1[t * 3 + 1] = omegaDip1[t * 3 + 1];
		omegaCsa1[t * 3 + 2] = omegaDip1[t * 3 + 2];
	}

	/////////////////////////
	// Dip2 for CH2 probes //
	/////////////////////////

	if (probe == 3)
	{
		// zD2
		R[6] = rA3[0] - rA1[0];
		R[7] = rA3[1] - rA1[1];
		R[8] = rA3[2] - rA1[2];
		norm = 1.0 / sqrt(R[6] * R[6] + R[7] * R[7] + R[8] * R[8]);
		R[6] *= norm;
		R[7] *= norm;
		R[8] *= norm;

		// supporting v
		R[0] = rA1[0] - rA2[0];
		R[1] = rA1[1] - rA2[1];
		R[2] = rA1[2] - rA2[2];
		norm = 1.0 / sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
		R[0] *= norm;
		R[1] *= norm;
		R[2] *= norm;

		// yD2 = -v x zD2
		R[3] = R[1] * R[8] - R[2] * R[7];
		R[4] = R[2] * R[6] - R[0] * R[8];
		R[5] = R[0] * R[7] - R[1] * R[6];
		norm = -1.0 / sqrt(R[3] * R[3] + R[4] * R[4] + R[5] * R[5]);
		R[3] *= norm;
		R[4] *= norm;
		R[5] *= norm;

		// xD2 = yD2 x zD2
		R[0] = R[4] * R[8] - R[5] * R[7];
		R[1] = R[5] * R[6] - R[3] * R[8];
		R[2] = R[3] * R[7] - R[4] * R[6];
		norm = 1.0 / sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
		R[0] *= norm;
		R[1] *= norm;
		R[2] *= norm;

		eq.setRotationMatrix(R);
		Omega = eq.getEulerAngles();
		omegaDip2[t * 3 + 0] = Omega.alpha;
		omegaDip2[t * 3 + 1] = Omega.beta;
		omegaDip2[t * 3 + 2] = Omega.gamma;
	}

	return;
}

void updateD200(int snapshot, int probeCode, double *omegaDip1, double *omegaDip2, double *omegaCsa1, dcomplex *D200Dip1, dcomplex *D200Csa1, dcomplex *D200Dip2)
{

	double cb;

	//////////////////////////////////////
	// D200Dip1 is common to all probes //
	//////////////////////////////////////

	cb = cos(omegaDip1[snapshot * 3 + 1]);
	D200Dip1[0].real(0.5 * (3.0 * cb * cb - 1.0)); 
	D200Dip1[0].imag(0.0); 

	////////////////////
	// D200Csa1 also //
	////////////////////

	cb = cos(omegaCsa1[snapshot * 3 + 1]);
	D200Csa1[0].real(0.5 * (3.0 * cb * cb - 1.0)); 
	D200Csa1[0].imag(0.0);

	//////////////////////
	// D200Dip2 for CH2 //
	//////////////////////

	if (probeCode == 3)
	{
		cb = cos(omegaDip2[snapshot * 3 + 1]);
		D200Dip2[0].real(0.5 * (3.0 * cb * cb - 1.0)); 
		D200Dip2[0].imag(0.0);
	}

	return;
}

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

///////////////////////////////////////////////////////////////////////
// BYPASS EVERYTHING TO FIR rNH, deltaCSA
///////////////////////////////////////////////////////////////////////
////relax.setSpdFile("6_j_probe-1.dat");
////relax.readSpd();
////
////relax.calcNMR();
////std::string results0 = relax.toString();
////
////std::cout << std::endl << std::endl;
////std::cout << "****************************" << std::endl;
////std::cout << "*** NMR RELAXATION RATES ***" << std::endl;
////std::cout << "****************************" << std::endl;
////std::cout << results0 << std::endl;
////return 0;
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

	int probeCode, nTimes, nRep, nAtoms, nExp;
	int id1, id2, id3;
	double timeAlpha;
	double *omegaDip1, *omegaCsa1, *omegaDip2;
	cdvector D200Dip1, D200Csa1, D200Dip2;
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

	switch (probeCode)
	{
		case 1:
		{
			omegaDip1 = new double[3 * nTimes];
			omegaCsa1 = new double[3 * nTimes];
			omegaDip2 = NULL;

			D200Dip1 = cdvector(nTimes);
			D200Csa1 = cdvector(nTimes);

			std::cout << "ID's of atoms H, N, X (an atom bonded to N): ";
			break;
		}
		case 2:
		{
			omegaDip1 = new double[3 * nTimes];
			omegaCsa1 = new double[3 * nTimes];
			omegaDip2 = NULL;

			D200Dip1 = cdvector(nTimes);
			D200Csa1 = cdvector(nTimes);

			std::cout << "ID's of atoms H, C, X (an atom bonded to C): ";
			break;
		}
		case 3:
		{
			omegaDip1 = new double[3 * nTimes];
			omegaDip2 = new double[3 * nTimes];
			omegaCsa1 = new double[3 * nTimes];

			D200Dip1 = cdvector(nTimes);
			D200Dip2 = cdvector(nTimes);
			D200Csa1 = cdvector(nTimes);

			std::cout << "ID's of atoms H, C, H: ";
			break;
		}
		default:
		{
			std::cout << "Invalid probe code";
			return -1;
		}
	}

	std::cin  >> id1 >> id2 >> id3;

	///////////////////////////////////////////////
	// Build the time series of the Euler angles //
	///////////////////////////////////////////////

	std::ifstream file(trjFile.c_str());
	if (!file.is_open()) {
		std::cerr << "Error: Could not open the file " << trjFile << std::endl;
		exit(-1);
	}

	std::string buffer;
	std::string line;
	size_t bytes_read = 0;

	int itmp;
	double dtmp;
	double *r1 = new double[3];
	double *r2 = new double[3];
	double *r3 = new double[3];
	std::string stmp;
	double *time = new double[nTimes];

	size_t chunk_size = 5 * 104857600; // 500 MB

	std::streampos fsize = file.tellg();	// The file pointer is currently at the beginning
	file.seekg(0, std::ios::end);	// Place the file pointer at the end of file
	fsize = file.tellg() - fsize;
	file.seekg(0, std::ios::beg);	// Place the file pointer at the beginning of file

	int nChunks = std::floor((double)fsize / (double)chunk_size);

	int snapshot = 0;
	int xyzline = 0;
	int atomsRead = 0;

	int count = 1;
	std::cout << "\n\n";

	while (true)
	{
		// Read a chunk of data
		char* chunk = new char[chunk_size + 1]; // +1 for the null terminator
		file.read(chunk, chunk_size);
		bytes_read = file.gcount(); // Actual number of bytes read

		if (bytes_read == 0)
		{
			delete[] chunk; // Clean up memory
			break; // End of file
		}

		chunk[bytes_read] = '\0'; // Null terminate the chunk
		buffer.append(chunk); // Append the chunk to the buffer
		delete[] chunk; // Clean up memory

		// Process each line in the current buffer
		size_t start = 0;

		while (true) {
			size_t end = buffer.find('\n', start);
			if (end == std::string::npos)
			{
				break; // No more complete lines in the buffer
			}
			line = buffer.substr(start, end - start); // Extract the line

			/*******************************/
			/* PROCESS LINES OF THE BUFFER */
			/*******************************/

			std::istringstream iss(line);

			// Read n atoms
			if (xyzline == 0)
			{
				iss >> itmp;
				xyzline++;
			}
			// Read snapshot
			else if (xyzline == 1)
			{
				iss >> time[snapshot] >> stmp;
				xyzline++;
			}
			// Read atoms
			else
			{
				atomsRead ++;

				if (atomsRead == id1)
					iss >> stmp >> r1[0] >> r1[1] >> r1[2];
				else if (atomsRead == id2)
					iss >> stmp >> r2[0] >> r2[1] >> r2[2];
				else if (atomsRead == id3)
					iss >> stmp >> r3[0] >> r3[1] >> r3[2];
				//else
				//	iss >> stmp >> dtmp >> dtmp >> dtmp;
			}

			if (atomsRead == nAtoms)
			{
				calculateEulerAngles(snapshot, probeCode, r1, r2, r3, omegaDip1, omegaCsa1, omegaDip2);
				updateD200(snapshot, probeCode, omegaDip1, omegaDip2, omegaCsa1, &D200Dip1[snapshot], &D200Csa1[snapshot], &D200Dip2[snapshot]);
				snapshot ++;
				xyzline = 0;
				atomsRead = 0;
			}

			/*******************************/

			start = end + 1; // Move to the next line
		}
		buffer.erase(0, start); // Remove processed lines from the buffer

		std::cout << "Processing trajectory file " << (double)(count * 100)/(double)nChunks << "%" << std::endl;
		count++;

	}

	// Process any remaining data in the buffer
	if (!buffer.empty())
	{
		std::cout << "Buffer is not empty and contains:" << std::endl;
		std::cout << buffer << std::endl;

		/*******************************/
		/* PROCESS LINES OF THE BUFFER */
		/*******************************/

		std::istringstream iss(line);

		// Read n atoms
		if (xyzline == 0)
		{
			iss >> itmp;
			xyzline++;
		}
		// Read snapshot
		else if (xyzline == 1)
		{
			iss >> time[snapshot] >> stmp;
			/// TO CANCEL /// std::cout << "Processing snapshot at time " << time[snapshot] << " ps." << std::endl;
			xyzline++;
		}
		// Read atoms
		else
		{
			atomsRead ++;

			if (atomsRead == id1)
				iss >> stmp >> r1[0] >> r1[1] >> r1[2];
			else if (atomsRead == id2)
				iss >> stmp >> r2[0] >> r2[1] >> r2[2];
			else if (atomsRead == id3)
				iss >> stmp >> r3[0] >> r3[1] >> r3[2];
			else
				iss >> stmp >> dtmp >> dtmp >> dtmp;
		}

		if (atomsRead == nAtoms)
		{
			calculateEulerAngles(snapshot, probeCode, r1, r2, r3, omegaDip1, omegaCsa1, omegaDip2);
			updateD200(snapshot, probeCode, omegaDip1, omegaDip2, omegaCsa1, &D200Dip1[snapshot], &D200Csa1[snapshot], &D200Dip2[snapshot]);

			snapshot++;
			xyzline = 0;
			atomsRead = 0;
		}

		/*******************************/
	}
	
	file.close(); // Close the file

	std::cout << std::endl << "Finished reading " << snapshot << " trajectory snapshots" << std::endl << std::endl;

	////////////////////////
	// Spectral densities //
	////////////////////////

	int nData = nTimes / nRep;
	int ncorr = nData / 2;

	double norm, dataDt;
	dvector g;
	cdvector cg;

	acfDriver ad;
	dataDt = (time[1] - time[0]) * 1.0e-3; // dt in ns
	dataDt *= timeAlpha; // Apply time strecthing, which is equivalent to scale the diffusion tensor by a constant 1/timeAlpha
	ad.setTimeData(0.0, dataDt);

	dvector  gDip1Dip1 = dvector(ncorr, 0.0);
	dvector  gCsa1Csa1 = dvector(ncorr, 0.0);
	cdvector gDip1Csa1 = cdvector(ncorr, std::complex<double>(0.0, 0.0));

	for (int p = 0; p < nRep; ++p)
	{
		std::cout << p + 1 << " / " << nRep << std::endl;

		////////////////////////////////////////////
		// Dip1 and Csa1 are common to all probes //
		////////////////////////////////////////////

		// Dip1-Dip1
		ad.setData({D200Dip1.begin() + p * nData, D200Dip1.begin() + (p + 1) * nData});
		ad.calculateACF(ncorr);
		g = ad.getACF();
		for (int i = 0; i < g.size(); ++i)
			gDip1Dip1[i] += g[i];
	
		// Csa1-Csa1
		ad.setData({D200Csa1.begin() + p * nData, D200Csa1.begin() + (p + 1) * nData});
		ad.calculateACF(ncorr);
		g = ad.getACF();
		for (int i = 0; i < g.size(); ++i)
			gCsa1Csa1[i] += g[i];
	
		// Dip1-Csa1
		ad.setData({D200Dip1.begin() + p * nData, D200Dip1.begin() + (p + 1) * nData}, {D200Csa1.begin() + p * nData, D200Csa1.begin() + (p + 1) * nData});
		ad.calculateCCF(ncorr);
		cg = ad.getCCF();
		for (int i = 0; i < cg.size(); ++i)
			gDip1Csa1[i] += cg[i];
	
	}
	std::cout << std::endl;

	dvector tData = dvector(gDip1Dip1.size(), 0.0);

	// J Dip1-Dip1
	double g0 =  1.0 / gDip1Dip1[0];
	for (int i = 0; i < gDip1Dip1.size(); ++i)
	{
		gDip1Dip1[i] *= g0;
		tData[i] = (double)i * dataDt;
	}
	ad.setCorrFun(gDip1Dip1);
	ad.saveACF("Dip1_acf_pre-fit.dat");
	std::cout << "- Fitting Dip1-Dip1 ACF with " << nExp << " exponential(s)..." << std::endl;
	gDip1Dip1 = fit_cf(tData, gDip1Dip1, nExp);
	ad.setCorrFun(gDip1Dip1);
	ad.saveACF("Dip1_acf.dat");
	for (int i = 0; i < 10240; ++i)
		gDip1Dip1.push_back(0.0);
        ad.setCorrFun(gDip1Dip1);
        ad.calculateRealSpectralDensity();
	dvector JDip1 = ad.getRealSpectralDensity();
	
	// J Csa1-Csa1
	g0 =  1.0 / gCsa1Csa1[0];
	for (int i = 0; i < gCsa1Csa1.size(); ++i)
		gCsa1Csa1[i] *= g0;
	ad.setCorrFun(gCsa1Csa1);
	ad.saveACF("Csa1_acf_pre-fit.dat");
	std::cout << "- Fitting Csa1-Csa1 ACF with " << nExp << " exponential(s)..." << std::endl;
	gCsa1Csa1 = fit_cf(tData, gCsa1Csa1, nExp);
	ad.setCorrFun(gCsa1Csa1);
	ad.saveACF("Csa1_acf.dat");
	for (int i = 0; i < 10240; ++i)
		gCsa1Csa1.push_back(0.0);
	ad.setCorrFun(gCsa1Csa1);
	ad.calculateRealSpectralDensity();
	dvector JCsa1 = ad.getRealSpectralDensity();

	// J Dip1-Csa1
	dvector gDip1Csa1Real = dvector(nData, 0.0);
	g0 = 1.0 / gDip1Csa1[0].real();
	for (int i = 0; i < gDip1Csa1.size(); ++i)
		gDip1Csa1Real[i] = gDip1Csa1[i].real() * g0; // since we know that this is a real function
	ad.setCorrFun(gDip1Csa1Real);
	ad.saveACF("Dip1-Csa1_ccf_pre-fit.dat");
	std::cout << "- Fitting Dip1-Csa1 CCF (which is real) with " << nExp << " exponential(s)..." << std::endl;
	gDip1Csa1Real = fit_cf(tData, gDip1Csa1Real, nExp);
	ad.setCorrFun(gDip1Csa1Real);
	ad.saveACF("Dip1-Csa1_ccf.dat");
	for (int i = 0; i < gDip1Csa1.size(); ++i)
		gDip1Csa1[i] = std::complex<double>(gDip1Csa1Real[i], 0.0);
	for (int i = 0; i < 10240; ++i)
		gDip1Csa1.push_back(std::complex<double>(0.0, 0.0));
	ad.setCrossCorrFun(gDip1Csa1);
	ad.calculateRealCrossSpectralDensity();
	dvector JcrossD1C = ad.getRealSpectralDensity();

	// Write j.dat (freq in 1/s, J in s)
	std::ofstream jFile("j.dat");

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
		case 3:
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
				ad.setData({D200Dip2.begin() + p * nData, D200Dip2.begin() + (p + 1) * nData});
				ad.calculateACF(ncorr);
				g = ad.getACF();
				for (int i = 0; i < g.size(); ++i)
					gDip2Dip2[i] += g[i];
			
				// Dip1-Dip2
				ad.setData({D200Dip1.begin() + p * nData, D200Dip1.begin() + (p + 1) * nData}, {D200Dip2.begin() + p * nData, D200Dip2.begin() + (p + 1) * nData});
				ad.calculateCCF(ncorr);
				cg = ad.getCCF();
				for (int i = 0; i < cg.size(); ++i)
					gDip1Dip2[i] += cg[i];

				// Dip2-Csa1
				ad.setData({D200Dip2.begin() + p * nData, D200Dip2.begin() + (p + 1) * nData}, {D200Csa1.begin() + p * nData, D200Csa1.begin() + (p + 1) * nData});
				ad.calculateCCF(ncorr);
				cg = ad.getCCF();
				for (int i = 0; i < cg.size(); ++i)
					gDip2Csa1[i] += cg[i];
			
			}

			// J Dip2-Dip2
			for (int i = 0; i < gDip2Dip2.size(); ++i)
			{
				gDip2Dip2[i] /= (double)nRep;
				tData[i] = (double)i * dataDt;
			}
			ad.setCorrFun(gDip2Dip2);
			ad.saveACF("Dip2_acf_pre-fit.dat");
			std::cout << "- Fitting Dip2-Dip2 ACF with " << nExp << " exponential(s)..." << std::endl;
			gDip2Dip2 = fit_cf(tData, gDip2Dip2, nExp);
			ad.setCorrFun(gDip2Dip2);
			ad.saveACF("Dip2_acf.dat");
			for (int i = 0; i < 10240; ++i)
				gDip2Dip2.push_back(0.0);
			ad.setCorrFun(gDip2Dip2);
			ad.calculateRealSpectralDensity();
			dvector JDip2 = ad.getRealSpectralDensity();
	
			// J Dip1-Dip2
			dvector gDip1Dip2Real = dvector(nData, 0.0);
			for (int i = 0; i < gDip1Dip2.size(); ++i)
			{
				gDip1Dip2[i] /= std::complex<double>((double)nRep, 0.0);
				gDip1Dip2Real[i] = gDip1Dip2[i].real(); // since we know that this is a real function
			}
			ad.setCorrFun(gDip1Dip2Real);
			ad.saveACF("Dip1-Dip2_ccf_pre-fit.dat");
			std::cout << "- Fitting Dip1-Dip2 CCF (which is real) with " << nExp << " exponential(s)..." << std::endl;
			gDip1Csa1Real = fit_cf(tData, gDip1Dip2Real, nExp);
			ad.setCorrFun(gDip1Dip2Real);
			ad.saveACF("Dip1-Dip2_ccf.dat");
			for (int i = 0; i < gDip1Dip2.size(); ++i)
				gDip1Dip2[i] = std::complex<double>(gDip1Dip2Real[i], 0.0);
			ad.setCrossCorrFun(gDip1Dip2);
			for (int i = 0; i < 10240; ++i)
				gDip1Dip2.push_back(std::complex<double>(0.0, 0.0));
			ad.setCrossCorrFun(gDip1Dip2);
			ad.calculateRealCrossSpectralDensity();
			dvector JcrossD1D2 = ad.getRealSpectralDensity();

			// J Dip2-Csa1
			dvector gDip2Csa1Real = dvector(nData, 0.0);
			for (int i = 0; i < gDip2Csa1.size(); ++i)
			{
				gDip2Csa1[i] /= std::complex<double>((double)nRep, 0.0);
				gDip2Csa1Real[i] = gDip2Csa1[i].real(); // since we know that this is a real function
			}
			ad.setCorrFun(gDip2Csa1Real);
			ad.saveACF("Dip2-Csa1_ccf_pre-fit.dat");
			std::cout << "- Fitting Dip2-Csa1 CCF (which is real) with " << nExp << " exponential(s)..." << std::endl;
			gDip2Csa1Real = fit_cf(tData, gDip2Csa1Real, nExp);
			ad.setCorrFun(gDip2Csa1Real);
			ad.saveACF("Dip2-Csa1_ccf.dat");
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

			std::ofstream jFile("j.dat");
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

	relax.setSpdFile("j.dat");
	relax.readSpd();

	relax.calcNMR();
	std::string results = relax.toString();

	std::cout << std::endl << std::endl;
	std::cout << "****************************" << std::endl;
	std::cout << "*** NMR RELAXATION RATES ***" << std::endl;
	std::cout << "****************************" << std::endl;
	std::cout << results << std::endl;

	return 0;
}
