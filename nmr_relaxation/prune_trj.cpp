#include <cstdlib>
#include <cstdio>

#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

#include <string>

#include "types.h"


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

	int nProbes;
	std::cout << "N probes = ";
	std::cin  >> nProbes;

	int *id1, *id2, *id3;
	id1 = new int[nProbes];
	id2 = new int[nProbes];
	id3 = new int[nProbes];

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
	
		std::cin  >> id1[np] >> id2[np] >> id3[np];
	}

	///////////////////////////////////////////////
	// Build the time series of the Euler angles //
	///////////////////////////////////////////////

	std::ofstream prunedTrjFile("6_pruned_trj.dat");

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
	double *r1 = new double[3 * nProbes];
	double *r2 = new double[3 * nProbes];
	double *r3 = new double[3 * nProbes];
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

	int count = 0;
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

		std::cout << "Processing trajectory file " << (double)(count * 100)/(double)nChunks << "%" << std::endl;

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
				prunedTrjFile << 3 * nProbes << std::endl;
				xyzline++;
			}
			// Read snapshot
			else if (xyzline == 1)
			{
				iss >> time[snapshot] >> stmp;
				prunedTrjFile << time[snapshot] << " ps" << std::endl;
				xyzline++;
			}
			// Read atoms
			else
			{
				atomsRead ++;
				for (int np = 0; np < nProbes; ++np)
				{
					if (atomsRead == id1[np])
						iss >> stmp >> r1[np * 3 + 0] >> r1[np * 3 + 1] >> r1[np * 3 + 2];
					else if (atomsRead == id2[np])
						iss >> stmp >> r2[np * 3 + 0] >> r2[np * 3 + 1] >> r2[np * 3 + 2];
					else if (atomsRead == id3[np])
						iss >> stmp >> r3[np * 3 + 0] >> r3[np * 3 + 1] >> r3[np * 3 + 2];
				}
			}

			if (atomsRead == nAtoms)
			{
				for (int np = 0; np < nProbes; ++np)
				{
					prunedTrjFile << "A1P" << np + 1 << "\t" << std::scientific << std::setprecision(12) << r1[np * 3 + 0] << "\t" << r1[np * 3 + 1] << "\t" << r1[np * 3 + 2] << std::endl;
					prunedTrjFile << "A2P" << np + 1 << "\t" << std::scientific << std::setprecision(12) << r2[np * 3 + 0] << "\t" << r2[np * 3 + 1] << "\t" << r2[np * 3 + 2] << std::endl;
					prunedTrjFile << "A3P" << np + 1 << "\t" << std::scientific << std::setprecision(12) << r3[np * 3 + 0] << "\t" << r3[np * 3 + 1] << "\t" << r3[np * 3 + 2] << std::endl;
				}

				snapshot ++;
				xyzline = 0;
				atomsRead = 0;
			}

			/*******************************/

			start = end + 1; // Move to the next line
		}
		buffer.erase(0, start); // Remove processed lines from the buffer

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
			prunedTrjFile << 3 * nProbes << std::endl;
			xyzline++;
		}
		// Read snapshot
		else if (xyzline == 1)
		{
			iss >> time[snapshot] >> stmp;
			prunedTrjFile << time[snapshot] << " ps" << std::endl;
			xyzline++;
		}
		// Read atoms
		else
		{
			atomsRead ++;
			for (int np = 0; np < nProbes; ++np)
			{
				if (atomsRead == id1[np])
					iss >> stmp >> r1[np * 3 + 0] >> r1[np * 3 + 1] >> r1[np * 3 + 2];
				else if (atomsRead == id2[np])
					iss >> stmp >> r2[np * 3 + 0] >> r2[np * 3 + 1] >> r2[np * 3 + 2];
				else if (atomsRead == id3[np])
					iss >> stmp >> r3[np * 3 + 0] >> r3[np * 3 + 1] >> r3[np * 3 + 2];
			}
		}

		if (atomsRead == nAtoms)
		{
			for (int np = 0; np < nProbes; ++np)
			{
				prunedTrjFile << "A1P" << np + 1 << "\t" << std::scientific << std::setprecision(12) << r1[np * 3 + 0] << "\t" << r1[np * 3 + 1] << "\t" << r1[np * 3 + 2] << std::endl;
				prunedTrjFile << "A2P" << np + 1 << "\t" << std::scientific << std::setprecision(12) << r2[np * 3 + 0] << "\t" << r2[np * 3 + 1] << "\t" << r2[np * 3 + 2] << std::endl;
				prunedTrjFile << "A3P" << np + 1 << "\t" << std::scientific << std::setprecision(12) << r3[np * 3 + 0] << "\t" << r3[np * 3 + 1] << "\t" << r3[np * 3 + 2] << std::endl;
			}

			snapshot ++;
			xyzline = 0;
			atomsRead = 0;
		}

		/*******************************/
	}
	
	file.close(); // Close the file

	std::cout << std::endl << "Finished pruning " << snapshot << " trajectory snapshots" << std::endl << std::endl;

	return 0;
}
