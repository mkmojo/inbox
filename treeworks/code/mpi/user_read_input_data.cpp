/***
 *  Project: TreeWorks
 *
 *  File: user_read_input_data.cpp
 *  Created: Mar 17, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#include <fstream>

#include "user_read_input_data.hpp"


bool read_point_data(const char* input_file, const MPI_data& mpi_data, InputData& data)
{
	// Format of input file:
	// first line = number of datapoints
	// second line = nuber of levels
	// rest lines = the three dimensions of each point in a new line (tab separated)
	
	double read_data_time_s = MPI_Wtime();

	std::ifstream input(input_file);
	if(!input) return false;

	unsigned int n_global = 0, levels = 0;

	input >> n_global; 		// total number of data points
	input >> levels;	// total number of levels for the octree

	if(n_global < mpi_data.size() ) {
		// not enough data points
		return false;
	}
	if(levels > MAX_LEVELS) {
		// cannot fit the cell id into 64 bits, so exit
		return false;
	}

	unsigned int offset = 0, n_local = 0;

	if(n_global % mpi_data.size() == 0) {
		n_local = n_global / mpi_data.size();
		offset = mpi_data.rank() * n_local;
	} else {
		unsigned int rem = n_global % mpi_data.size();

		if(mpi_data.rank() < rem) {
			n_local = (n_global / mpi_data.size()) + 1;
			offset = mpi_data.rank() * n_local;
		} else {
			n_local = n_global / mpi_data.size();
			offset = mpi_data.rank() * n_local + rem;
		}
	}

	if(!data.input_data_alloc(n_local, n_global, levels)) return false;

	// IMPROVE: SEEK OR SOMETHING
	std::string junk;
	for(unsigned int i = 0; i <= offset; ++i) std::getline(input, junk);

	for(unsigned int i = 0; i < n_local; ++i)
	{
		if(input.eof()) return false; // file ends unexpectedly

		input >> data.point_array_[i].x_;
		input >> data.point_array_[i].y_;
		input >> data.point_array_[i].z_;
	}

	input.close();

	MPI_Barrier(mpi_data.comm());
	double read_data_time_e = MPI_Wtime();

	if(mpi_data.rank() == 0)
		std::cout << "Input reading time: "
			<< (read_data_time_e - read_data_time_s) * 1000 << "ms." << std::endl;

	return true;
} // read_point_data
