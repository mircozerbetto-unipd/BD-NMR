#include <cstdlib>
#include <cstdio>

#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

#include <string>
#include <string_view> // For efficient string handling
#include <charconv>    // For std::from_chars (C++17+)

#include <algorithm>   // For std::min

#include "types.h"

#include "EulerQuaternion.h"
#include "acfDriver.h"
#include "nmr.h"
#include "fit_cf.h"

void skip_leading_whitespace(const char*& ptr, const char* end)
{
	while (ptr < end && (*ptr == ' ' || *ptr == '\t' || *ptr == '\r'))
		ptr++;
}

bool parse_double_from_sv(const std::string_view& sv, double& value)
{
	if (sv.empty()) return false;
	std::string temp_str(sv.data(), sv.length());
	char* end_ptr;
	errno = 0; // Clear errno before call
	value = std::strtod(temp_str.c_str(), &end_ptr);

	// Check for conversion errors:
	if (end_ptr == temp_str.c_str() || // No characters converted
		(errno == ERANGE && (value == HUGE_VAL || value == -HUGE_VAL)) // Overflow/underflow
	   )
	{
        	return false; // Conversion failed or out of range
	}

	return true;
}

void calculateEulerAngles(int probe, double *rA1, double *rA2, double *rA3, double *omegaDip1, double *omegaCsa1, double *omegaDip2, double betaDipCsa)
{
	
	double norm;
	double R[9];
	EulerAngles Omega;
	EulerQuaternion eq;

	////////////////////////////////////////////////////////
	// First dipolar frame - same for all kinds of probes //	
	////////////////////////////////////////////////////////

	// zD1
	R[6] = rA1[0] - rA2[0];
	R[7] = rA1[1] - rA2[1];
	R[8] = rA1[2] - rA2[2];
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

	eq.setRotationMatrix_noQ(R);
	Omega = eq.getEulerAngles();
	omegaDip1[0] = Omega.alpha;
	omegaDip1[1] = Omega.beta;
	omegaDip1[2] = Omega.gamma;

	////////////////////////
	// CSA for all probes //
	////////////////////////
	
	// Here, CSA frame is assumed to be at (0, betaDipCsa, 0) deg from dipolar frame for NH, while no tilt for CHn
	if (probe == 1 || probe == 2)
	{

		eq.setEulerAngles_noQ(0.0, betaDipCsa * M_PI / 180.0, 0.0);
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

		eq.setRotationMatrix_noQ(RC);
		Omega = eq.getEulerAngles();
		omegaCsa1[0] = Omega.alpha;
		omegaCsa1[1] = Omega.beta;
		omegaCsa1[2] = Omega.gamma;
	}

	/////////////////////////
	// Dip2 for CH2 probes //
	/////////////////////////

	if (probe == 3)
	{
		// zD2
		R[6] = rA3[0] - rA2[0];
		R[7] = rA3[1] - rA2[1];
		R[8] = rA3[2] - rA2[2];
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

		eq.setRotationMatrix_noQ(R);
		Omega = eq.getEulerAngles();
		omegaDip2[0] = Omega.alpha;
		omegaDip2[1] = Omega.beta;
		omegaDip2[2] = Omega.gamma;
	}

	return;
}

void updateD200(int probeCode, double *omegaDip1, double *omegaDip2, double *omegaCsa1, dcomplex *D200Dip1, dcomplex *D200Csa1, dcomplex *D200Dip2)
{

	double cb;

	//////////////////////////////////////
	// D200Dip1 is common to all probes //
	//////////////////////////////////////

	cb = cos(omegaDip1[1]);
	D200Dip1[0].real(0.5 * (3.0 * cb * cb - 1.0)); 
	D200Dip1[0].imag(0.0); 

	////////////////////
	// D200Csa1 also //
	////////////////////

	cb = cos(omegaCsa1[1]);
	D200Csa1[0].real(0.5 * (3.0 * cb * cb - 1.0)); 
	D200Csa1[0].imag(0.0);

	//////////////////////
	// D200Dip2 for CH2 //
	//////////////////////

	if (probeCode == 3)
	{
		cb = cos(omegaDip2[1]);
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

	int probeCode, nTimes, nRep, nAtoms, nExp;
	double timeAlpha;
	double *omegaDip1, *omegaCsa1, *omegaDip2, *betaDipCsa;
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

	int nProbes;
	std::cout << "N probes = ";
	std::cin  >> nProbes;

	int *id1, *id2, *id3;
	id1 = new int[nProbes];
	id2 = new int[nProbes];
	id3 = new int[nProbes];
	betaDipCsa = new double[nProbes];

	for (int np = 0; np < nProbes; ++np)
	{
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
	
		std::cin  >> id1[np] >> id2[np] >> id3[np] >> betaDipCsa[np];
	}

	omegaDip1 = new double[nProbes * 3];
	omegaCsa1 = new double[nProbes * 3];
	omegaDip2 = new double[nProbes * 3];
	
	D200Dip1 = cdvector(nProbes);
	D200Dip2 = cdvector(nProbes);
	D200Csa1 = cdvector(nProbes);

	///////////////////////////////////////////////
	// Build the time series of the Euler angles //
	///////////////////////////////////////////////

	std::ofstream dip1File("6_D200Dip1.dat");
	std::ofstream csa1File("6_D200Csa1.dat");
	std::ofstream dip2File("6_D200Dip2.dat");

	std::string buffer;
	std::string line;
	size_t bytes_read = 0;

	int itmp;
	double dtmp;
	double *r1 = new double[3 * nProbes];
	double *r2 = new double[3 * nProbes];
	double *r3 = new double[3 * nProbes];
	std::string stmp;
	double *time = new double[nTimes];

	size_t chunk_size = 1024 * 1024 * 1024; // 1 GB
	std::ifstream file(trjFile.c_str(), std::ios::binary); // Open in binary mode for raw data
	if (!file.is_open())
	{
		std::cerr << "Error opening file: " << trjFile.c_str() << std::endl;
		return -1;
	}

	file.seekg(0, std::ios::end);
	long long fsize = file.tellg();
	file.seekg(0, std::ios::beg);

	int nChunks = std::floor((double)fsize / (double)chunk_size);
	if (fsize % chunk_size != 0) // Account for the last partial chunk
	{
		nChunks++;
	}

	std::vector<char> chunk_buffer(chunk_size); // Use std::vector for automatic memory management
	std::string partial_line_accumulator; // To store partial lines spanning chunk boundaries

	int snapshot = 0;
	int xyzline = 0;
	int atomsRead = 0;

	int count = 0;
	int current_chunk_idx = 0;
	std::cout << "\n\n";

	while (true)
	{
		// Read a chunk of data
		file.read(chunk_buffer.data(), chunk_size);
		std::streamsize bytes_read = file.gcount();

		if (bytes_read == 0)
		{
			break; // End of file
		}

		partial_line_accumulator.append(chunk_buffer.data(), bytes_read);

		size_t current_pos_in_accumulator = 0;
		size_t newline_pos = 0;

		std::cout << "Processing trajectory file " << (double)(count * 100)/(double)nChunks << "%" << std::endl;

		while ((newline_pos = partial_line_accumulator.find('\n', current_pos_in_accumulator)) != std::string_view::npos)
		{
			std::string_view line_sv(partial_line_accumulator.data() + current_pos_in_accumulator, newline_pos - current_pos_in_accumulator);

			/*******************************/
			/* PROCESS LINES OF THE BUFFER */
			/*******************************/

			const char* current_ptr = line_sv.data();
			const char* end_ptr = line_sv.data() + line_sv.length();

			// Read n atoms
			if (xyzline == 0)
			{
				skip_leading_whitespace(current_ptr, end_ptr);
				auto [p, ec] = std::from_chars(current_ptr, end_ptr, itmp);
				xyzline++;
			}
			// Read snapshot
			else if (xyzline == 1)
			{
				skip_leading_whitespace(current_ptr, end_ptr);
				double f_val;
				size_t first_space_after_float = line_sv.find(' ', current_ptr - line_sv.data());
				std::string_view float_sv = line_sv.substr(current_ptr - line_sv.data(), first_space_after_float - (current_ptr - line_sv.data()));
				parse_double_from_sv(float_sv, f_val);
				time[snapshot] = f_val;	
				xyzline++;
			}
			// Read atoms
			else
			{
				atomsRead ++;

				size_t first_space_pos = line_sv.find('\t');
				std::string_view atom_name_sv = line_sv.substr(0, first_space_pos);

				current_ptr = line_sv.data() + first_space_pos + 1;

				double x, y, z;

				skip_leading_whitespace(current_ptr, end_ptr);
				const char* start_of_num = current_ptr;
				while (current_ptr < end_ptr && *current_ptr != '\t') current_ptr++;
				parse_double_from_sv(std::string_view(start_of_num, current_ptr - start_of_num), x);

				skip_leading_whitespace(current_ptr, end_ptr);
				start_of_num = current_ptr;
				while (current_ptr < end_ptr && *current_ptr != '\t') current_ptr++;
				parse_double_from_sv(std::string_view(start_of_num, current_ptr - start_of_num), y);
				

				skip_leading_whitespace(current_ptr, end_ptr);
				start_of_num = current_ptr;
				parse_double_from_sv(std::string_view(start_of_num, end_ptr - start_of_num), z);

				for (int np = 0; np < nProbes; ++np)
				{
					if (atomsRead == id1[np])
					{
						r1[np * 3 + 0] = x;
						r1[np * 3 + 1] = y;
						r1[np * 3 + 2] = z;
					}
					else if (atomsRead == id2[np])
					{
						r2[np * 3 + 0] = x;
						r2[np * 3 + 1] = y;
						r2[np * 3 + 2] = z;
					}
					else if (atomsRead == id3[np])
					{
						r3[np * 3 + 0] = x;
						r3[np * 3 + 1] = y;
						r3[np * 3 + 2] = z;
					}
				}
			}

			if (atomsRead == nAtoms)
			{
				dip1File << time[snapshot] << "\t";
				csa1File << time[snapshot] << "\t";
				dip2File << time[snapshot] << "\t";

				for (int np = 0; np < nProbes; ++np)
				{
					calculateEulerAngles(probeCode, &r1[np * 3], &r2[np * 3], &r3[np * 3], &omegaDip1[np * 3], &omegaCsa1[np * 3], &omegaDip2[np * 3], betaDipCsa[np]);
					updateD200(probeCode, &omegaDip1[np * 3], &omegaDip2[np * 3], &omegaCsa1[np * 3], &D200Dip1[np], &D200Csa1[np], &D200Dip2[np]);

					dip1File << std::scientific << std::setprecision(12) << D200Dip1[np].real() << "\t";
					csa1File << std::scientific << std::setprecision(12) << D200Csa1[np].real() << "\t";
					dip2File << std::scientific << std::setprecision(12) << D200Dip2[np].real() << "\t";
				}

				dip1File << "\n";
				csa1File << "\n";
				dip2File << "\n";

				snapshot ++;
				xyzline = 0;
				atomsRead = 0;
			}

			/*******************************/

			current_pos_in_accumulator = newline_pos + 1;

		}

		if (current_pos_in_accumulator > 0) {
			partial_line_accumulator = partial_line_accumulator.substr(current_pos_in_accumulator);
		}

		count++;

	}

	// Process any remaining data in the buffer
	if (!partial_line_accumulator.empty())
	{
		/* PROCESS CONTENT *************/

		const char* current_ptr = partial_line_accumulator.data();
		const char* end_ptr = partial_line_accumulator.data() + partial_line_accumulator.length();

		// Read n atoms
		if (xyzline == 0)
		{
			skip_leading_whitespace(current_ptr, end_ptr);
			auto [p, ec] = std::from_chars(current_ptr, end_ptr, itmp);
			xyzline++;
		}
		// Read snapshot
		else if (xyzline == 1)
		{
			skip_leading_whitespace(current_ptr, end_ptr);
			double f_val;
			size_t first_space_after_float = partial_line_accumulator.find(' ', current_ptr - partial_line_accumulator.data());
			std::string_view float_sv = partial_line_accumulator.substr(current_ptr - partial_line_accumulator.data(), first_space_after_float - (current_ptr - partial_line_accumulator.data()));
			parse_double_from_sv(float_sv, f_val);
			time[snapshot] = f_val;	
			xyzline++;
		}
		// Read atoms
		else
		{
			atomsRead ++;

			size_t first_space_pos = partial_line_accumulator.find('\t');
			std::string_view atom_name_sv = partial_line_accumulator.substr(0, first_space_pos);
			current_ptr = partial_line_accumulator.data() + first_space_pos + 1;
			double x, y, z;

			skip_leading_whitespace(current_ptr, end_ptr);
			const char* start_of_num = current_ptr;
			while (current_ptr < end_ptr && *current_ptr != '\t') current_ptr++;
			parse_double_from_sv(std::string_view(start_of_num, current_ptr - start_of_num), x);

			skip_leading_whitespace(current_ptr, end_ptr);
			start_of_num = current_ptr;
			while (current_ptr < end_ptr && *current_ptr != '\t') current_ptr++;
			parse_double_from_sv(std::string_view(start_of_num, current_ptr - start_of_num), y);
			

			skip_leading_whitespace(current_ptr, end_ptr);
			start_of_num = current_ptr;
			parse_double_from_sv(std::string_view(start_of_num, end_ptr - start_of_num), z);

			for (int np = 0; np < nProbes; ++np)
			{
				if (atomsRead == id1[np])
				{
					r1[np * 3 + 0] = x;
					r1[np * 3 + 1] = y;
					r1[np * 3 + 2] = z;
				}
				else if (atomsRead == id2[np])
				{
					r2[np * 3 + 0] = x;
					r2[np * 3 + 1] = y;
					r2[np * 3 + 2] = z;
				}
				else if (atomsRead == id3[np])
				{
					r3[np * 3 + 0] = x;
					r3[np * 3 + 1] = y;
					r3[np * 3 + 2] = z;
				}
			}
		}

		if (atomsRead == nAtoms)
		{
			dip1File << time[snapshot] << "\t";
			csa1File << time[snapshot] << "\t";
			dip2File << time[snapshot] << "\t";

			for (int np = 0; np < nProbes; ++np)
			{
				calculateEulerAngles(probeCode, &r1[np * 3], &r2[np * 3], &r3[np * 3], &omegaDip1[np * 3], &omegaCsa1[np * 3], &omegaDip2[np * 3], betaDipCsa[np]);
				updateD200(probeCode, &omegaDip1[np * 3], &omegaDip2[np * 3], &omegaCsa1[np * 3], &D200Dip1[np], &D200Csa1[np], &D200Dip2[np]);

				dip1File << std::scientific << std::setprecision(12) << D200Dip1[np].real() << "\t";
				csa1File << std::scientific << std::setprecision(12) << D200Csa1[np].real() << "\t";
				dip2File << std::scientific << std::setprecision(12) << D200Dip2[np].real() << "\t";
			}

			dip1File << "\n";
			csa1File << "\n";
			dip2File << "\n";

			snapshot ++;
			xyzline = 0;
			atomsRead = 0;
		}

	}
	
	file.close(); // Close the file

	dip1File.close();
	csa1File.close();
	dip2File.close();

	std::cout << std::endl << "Finished reading " << snapshot << " trajectory snapshots" << std::endl << std::endl;

	return 0;
}
