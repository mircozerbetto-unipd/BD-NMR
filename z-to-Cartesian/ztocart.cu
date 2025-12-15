#include <iostream>
#include <sstream>
#include <cstdio>
#include "matrix_operations.h"
#include "file_reader.h"

#include "EulerQuaternion.h"

#include "zmat.h"

#undef DEBUG

// Update q
#define TPB_1D 512
__global__ void updateQ(double *q0, double *deltaQ, double *q, int nq)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < nq)
		q[i] = q0[i] + deltaQ[i];
	
	return;
}

// Update Z-Matrix
void setNewZMATRIX(molecule *mol, double *q)
{
	zmatline zml;

	// Atom 2
	zml = mol->getAtom(2).getZMatrixEntry();
	zml.d = q[0];
	mol->changeZMatrix(2, zml);

	// Atom 3
	zml = mol->getAtom(3).getZMatrixEntry();
	zml.d = q[1];
	zml.theta = q[2];
	mol->changeZMatrix(3, zml);
	
	int pos = 3; 
	
	for (int i = 4; i <= mol->getNAtoms(); ++i)
	{
		zml = mol->getAtom(i).getZMatrixEntry();
		zml.d = q[pos];
		zml.theta = q[pos + 1];
		zml.phi = q[pos + 2];
		mol->changeZMatrix(i, zml);
		pos += 3;
	}

	return;
}

// Apply rotation to molecule
void rotateXYZ(molecule *mol, EulerQuaternion ea)
{
	vectorOfDoubles xyz0 = vectorOfDoubles(3, 0.0);
	vectorOfDoubles xyz  = vectorOfDoubles(3, 0.0);

	double *R = new double[9];
	ea.getRotationMatrix(R);

	for (int i = 1; i <= mol->getNAtoms(); ++i)
	{
		xyz0 = mol->getAtom(i).getPosition();
		xyz[0] = R[0] * xyz0[0] + R[1] * xyz0[1] + R[2] * xyz0[2];
		xyz[1] = R[3] * xyz0[0] + R[4] * xyz0[1] + R[5] * xyz0[2];
		xyz[2] = R[6] * xyz0[0] + R[7] * xyz0[1] + R[8] * xyz0[2];
		mol->setAtomPosition(i, xyz);
	}

	return;
} 

int updateTrj_buffer(char* buffer, int current_pos, int max_size, molecule *mol, double time)
{
	int written_chars = 0;
	int nAtoms = mol->getNAtoms();

	// Use snprintf to write to the buffer, starting at buffer[current_pos]

	// 1. Write header line 1: NAtoms
	written_chars += snprintf(buffer + current_pos + written_chars, max_size - (current_pos + written_chars), " %d\n", nAtoms);

	// 2. Write header line 2: Time (use high precision)
	written_chars += snprintf(buffer + current_pos + written_chars, max_size - (current_pos + written_chars)," %.15f ps\n", time);

     	// 3. Write coordinate lines
	for (int i = 1; i <= nAtoms; i++)
		written_chars += snprintf(buffer + current_pos + written_chars, max_size - (current_pos + written_chars), "%s\t%.15f\t%.15f\t%.15f\n", mol->getAtom(i).getAtomType().c_str(), mol->getAtom(i).x(), mol->getAtom(i).y(), mol->getAtom(i).z());

	// Check if the output truncated due to buffer size (optional but recommended)
	if (written_chars >= max_size - current_pos) {
		std::cerr << "Warning: Buffer capacity exceeded!" << std::endl;
	}

    return written_chars;
}

int main()
{
	////////////////
	// Init MAGMA //
	////////////////

	magma_init();
	magma_queue_t queue;
	magma_queue_create(0, &queue); // Create a MAGMA queue


	//////////////////////////
	// Init molecule object //
	//////////////////////////

	std::string str;
	char* zmtlibHome;
	zmtlibHome = getenv ("ZMATLIB_HOME");
	if (zmtlibHome == NULL)
	{
		std::cout << std::endl << "ERROR: the ZMATLIB_HOME envinronment variable is not set. Please set the path or use the script salem" << std::endl << std::endl;
		exit(-1);
	}
	config conf;
	str.assign(zmtlibHome); str.append("/sesto_1.0/sesto");
	conf.setSesto(str);
	str.assign(zmtlibHome); str.append("/data/vdw.dat");
	conf.setVdwFile(str);
	str.assign(zmtlibHome); str.append("/data/aminoacids.dat");
	conf.setAminoAcidsFile(str);
	str.assign(zmtlibHome); str.append("/data/atoms_masses.dat");
	conf.setMassFile(str);

	molecule mol(&conf);

	////////////////////////////////
	// Input dimension and scalar //
	////////////////////////////////

	int nZ, nTimes;
	int *refAtoms = new int[4];

	std::cout << "Number of internal coordinates: ";
	std::cin >> nZ;
	std::cout << nZ;
	std::cout << "Number of trajectory snapshots: ";
	std::cin >> nTimes;
	std::cout << nTimes;

	////////////////
	// Read files //
	////////////////

	std::string PDBFile, TFile, sqKFile, trjFile;

	std::cout << "Reference PDB file: ";
	std::cin >> PDBFile;
	std::cout << PDBFile;
	std::cout << "List of reference atoms (4 integers, sepatated by space): ";
	std::cin >> refAtoms[0] >> refAtoms[1] >> refAtoms[2] >> refAtoms[3];
	std::cout << refAtoms[0] << ", " << refAtoms[1] << ", " << refAtoms[2] << ", " << refAtoms[3];
	std::vector<double> q0 = read_PDB_file(PDBFile, refAtoms, &mol);

	std::cout << "T matrix file: ";
	std::cin >> TFile;
	std::cout << TFile;
	std::vector<double> T = read_matrix_from_file(TFile, nZ);

	std::cout << "sqrt(K) matrix file: ";
	std::cin >> sqKFile;
	std::cout << sqKFile;
	std::vector<double> sqK = read_matrix_from_file(sqKFile, nZ);
	
	std::cout << "SFB Trajectory file: ";
	std::cin >> trjFile;

	/////////////////////////
	// Inverting T and sqK //
	/////////////////////////

	std::vector<double> invT = invert_matrix(T, nZ);
	std::vector<double> invSqK = invert_matrix(sqK, nZ);

	/////////////////////
	// Multiplications //
	/////////////////////

	std::vector<double> A = multiply_matrices(invSqK, invT, nZ, queue);
	double *d_A;
	magma_dmalloc(&d_A, nZ * nZ);
	cudaMemcpy(d_A, A.data(), nZ * nZ * sizeof(double), cudaMemcpyHostToDevice);

	//////////////////////
	// Recover q from z //
	//////////////////////

	double *q = nullptr;
       	double *zt = nullptr;
	cudaMallocHost((void**)&q, nZ * sizeof(double));
	cudaMallocHost((void**)&zt, nZ * sizeof(double));


	double *d_q0, *d_DeltaQ, *d_q;
	magma_dmalloc(&d_q0, nZ);
	cudaMemcpy(d_q0, q0.data(), nZ * sizeof(double), cudaMemcpyHostToDevice);

	magma_dmalloc(&d_DeltaQ, nZ);
	magma_dmalloc(&d_q, nZ);

	double *d_zt;
	magma_dmalloc(&d_zt, nZ);

	FILE *iFile = fopen(trjFile.c_str(), "r");

	int itmp;
	double time, al, be, ga;
	EulerQuaternion ea;

	/////////////////////////////////////
	// Output trajectory in XYZ format //
	/////////////////////////////////////

	std::string trjOutFileName = "trj.xyz";
	FILE *oFile = fopen(trjOutFileName.c_str(), "w");


#ifdef DEBUG
	std::ofstream qfile("q.dat");
#endif

	////////////////////
	// Buffered input //
	///////////////////

	// Define the chunk size for reading the input file (e.g., 100 MB)
	const size_t INPUT_CHUNK_SIZE = 100 * 1024 * 1024; // 100 MB
	char* input_buffer = new char[INPUT_CHUNK_SIZE + 1]; // +1 for null terminator

	char* current_read_pos = input_buffer; // Pointer for parsing
	char* buffer_end = input_buffer; // Marker for end of valid data
	size_t bytes_left = 0; // Bytes remaining in the buffer

	/////////////////////
	// Buffered output //
	/////////////////////
	// Define the size of the chunk (in number of snapshots) to write at once
	const int CHUNK_SIZE = 100; 

	// Approximate maximum line length (AtomType + 3 doubles + tabs/newline)
	// Assuming max 200 chars per line for safety
	const int MAX_LINE_LENGTH = 200; 
	int nAtoms = mol.getNAtoms();

	// Buffer size = (Header lines + Data lines) * Max line length * CHUNK_SIZE
	// Header is 2 lines per snapshot. Data is nAtoms lines per snapshot.
	const size_t BUFFER_CAPACITY = (2 + nAtoms) * MAX_LINE_LENGTH * CHUNK_SIZE;

	// Dynamically allocate the buffer
	char* buffer = new char[BUFFER_CAPACITY];
	int buffer_pos = 0; // Current position in the buffer
	int snapshot_count = 0; // Snapshots currently in the buffer

	////////
	// GO //
	////////

////////	int snapshots_processed = 0;
////////    
////////	while (snapshots_processed < nTimes) 
////////	{
////////		// -----------------------------------------------------
////////		// Input Reading/Parsing (Replaces fscanf calls)
////////		// -----------------------------------------------------
////////
////////		// Estimate bytes needed for one snapshot:
////////		// Header (itmp, time, angles) + nZ lines (zt)
////////		// Max (5 doubles * 1 line + nZ doubles * 1 line) * MAX_LINE_LENGTH
////////		const int BYTES_NEEDED_PER_SNAPSHOT = (5 + nZ) * 20; // A safe, generous estimate
////////	
////////		// Check if enough data is left in the buffer to read one snapshot
////////		if (buffer_end - current_read_pos < BYTES_NEEDED_PER_SNAPSHOT)
////////		{
////////			// Move any remaining partial data to the start of the buffer
////////			if (current_read_pos != input_buffer && bytes_left > 0)
////////			{
////////				memmove(input_buffer, current_read_pos, bytes_left);
////////			}
////////	    
////////			// Set the new starting point for the fread call
////////			size_t read_start_pos = bytes_left;
////////		    
////////			// Read next chunk from the file
////////			size_t bytes_read = fread(input_buffer + read_start_pos, 1, INPUT_CHUNK_SIZE - read_start_pos, iFile);
////////		    
////////			if (bytes_read == 0)
////////			{
////////				// End of file or error reading
////////				break;
////////			}
////////		    
////////			// Update buffer pointers
////////			bytes_left += bytes_read;
////////			buffer_end = input_buffer + bytes_left;
////////			*buffer_end = '\0'; // Null-terminate for sscanf safety
////////			current_read_pos = input_buffer;
////////		}
////////
////////		// --- Attempt to Parse the Snapshot Data ---
////////		
////////		int bytes_consumed = 0; // Bytes used by the sscanf block
////////		int n_scanned;
////////		
////////		// 1. Scan Header Line (itmp, time, al, be, ga)
////////		// Note: sscanf returns the number of successfully matched items
////////		n_scanned = sscanf(current_read_pos, "%d %lf %lf %lf %lf %n", &itmp, &time, &al, &be, &ga, &bytes_consumed);
////////		
////////		if (n_scanned != 5)
////////		{
////////			// Not enough data for the full header, break out of loop
////////			// This might be a partial line at the very end of the file
////////			break;
////////		}
////////		current_read_pos += bytes_consumed; // Advance past the header line
////////		
////////		// 2. Scan the nZ coordinate lines into zt
////////		for (int i = 0; i < nZ; ++i)
////////		{
////////			bytes_consumed = 0;
////////			n_scanned = sscanf(current_read_pos, "%lf %n", &zt[i], &bytes_consumed);
////////			
////////			if (n_scanned != 1)
////////			{
////////				std::cerr << "Error parsing zt[" << i << "] at snapshot " << snapshots_processed << std::endl;
////////				goto cleanup_and_exit; // Use a goto to safely exit the nested loop structure
////////			}
////////			current_read_pos += bytes_consumed; // Advance past the zt value
////////		}
////////		
////////		// Update the count
////////		snapshots_processed++;
////////		
////////		// -----------------------------------------------------
////////		// GPU Computation and Output Buffering (As you defined)
////////		// -----------------------------------------------------
////////		
////////		// Recover q from z
////////		cudaMemcpyAsync(d_zt, zt, nZ * sizeof(double), cudaMemcpyHostToDevice);
////////		multiply_matrix_vector(d_A, d_zt, d_DeltaQ, nZ, queue);
////////		updateQ<<<1 + nZ / TPB_1D, TPB_1D>>>(d_q0, d_DeltaQ, d_q, nZ);
////////		cudaMemcpyAsync(q, d_q, nZ * sizeof(double), cudaMemcpyDeviceToHost);
////////		magma_queue_sync(queue);
////////
////////EBUG
////////		for (int i = 0; i < nZ; ++i)
////////			qfile << q[i] << "\t";
////////		qfile << std::endl;
////////
////////
////////		// Recover Cartesian coordinates from q and Euler angles
////////		setNewZMATRIX(&mol, q);
////////		ea.setEulerAngles(al, be, ga);
////////		rotateXYZ(&mol, ea);
////////
////////		// Buffer the trajectory output
////////		int chars_added = updateTrj_buffer(buffer, buffer_pos, BUFFER_CAPACITY, &mol, time);
////////		buffer_pos += chars_added;
////////		snapshot_count++;
////////		
////////		// Write chunk if buffer is full
////////		if (snapshot_count >= CHUNK_SIZE)
////////		{
////////		    fwrite(buffer, 1, buffer_pos, oFile);
////////		    buffer_pos = 0;
////////		    snapshot_count = 0;
////////		}
////////	}
////////
////////cleanup_and_exit:
////////	// Process any remaining data in the buffer (Last partial output chunk)
////////	if (buffer_pos > 0)
////////	{
////////		std::cout << "Writing final partial output buffer to file." << std::endl;
////////		fwrite(buffer, 1, buffer_pos, oFile);
////////	}
////////
////////	// Cleanup the buffer
////////	delete[] buffer;
////////	delete[] input_buffer;

	for (int nt = 0; nt < nTimes; ++nt)
	{
		fscanf(iFile, "%d %lf %lf %lf %lf", &itmp, &time, &al, &be, &ga);
		for (int i = 0; i < nZ; ++i)
			fscanf(iFile, "%lf", &zt[i]);

		// Recover q from z
		cudaMemcpyAsync(d_zt, zt, nZ * sizeof(double), cudaMemcpyHostToDevice);
		multiply_matrix_vector(d_A, d_zt, d_DeltaQ, nZ, queue);
		updateQ<<<1 + nZ / TPB_1D, TPB_1D>>>(d_q0, d_DeltaQ, d_q, nZ);
		cudaMemcpyAsync(q, d_q, nZ * sizeof(double), cudaMemcpyDeviceToHost);
		magma_queue_sync(queue);

#ifdef DEBUG
		for (int i = 0; i < nZ; ++i)
			qfile << q[i] << "\t";
		qfile << std::endl;
#endif

		// Recover Cartesian coordinates from q and Euler angles
		setNewZMATRIX(&mol, q);
		ea.setEulerAngles(al, be, ga);
		rotateXYZ(&mol, ea);

		// Append the snapshot data to the buffer
		int chars_added = updateTrj_buffer(buffer, buffer_pos, BUFFER_CAPACITY, &mol, time);
		buffer_pos += chars_added;
		snapshot_count++;

    		// Check if the buffer is full (or if we reached the CHUNK_SIZE)
		if (snapshot_count >= CHUNK_SIZE)
		{
			// Write the entire buffer content to the file in one go
			size_t written = fwrite(buffer, 1, buffer_pos, oFile);

	                if (written != buffer_pos)
		       	{
				std::cerr << "Error: Failed to write full buffer to file." << std::endl;
			}

			// Reset the buffer state
			buffer_pos = 0;
			snapshot_count = 0;
		}

	}

	// Process any remaining data in the buffer (Last partial chunk)
	if (buffer_pos > 0)
	{
		std::cout << "Writing final partial buffer to file." << std::endl;
		fwrite(buffer, 1, buffer_pos, oFile);
	}

	// Cleanup the buffer
	delete[] buffer;

	fclose(iFile);
	fclose(oFile);

#ifdef DEBUG
	qfile.close();
#endif
	cudaFreeHost(q);
	cudaFreeHost(zt);

	magma_queue_destroy(queue);
	magma_finalize();

	return 0;
}
