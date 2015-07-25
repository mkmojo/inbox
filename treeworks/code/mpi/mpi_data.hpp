/***
 *  Project: TreeWorks
 *
 *  File: mpi_data.hpp
 *  Created: Mar 17, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#ifndef MPI_DATA_HPP
#define MPI_DATA_HPP

#include <mpi.h>

class MPI_data {
	private:
		MPI_Comm comm_;	// the MPI communicator
		int rank_;		// rank of the processor
		int size_; 		// size of the comm

	public:
		MPI_data(MPI_Comm comm = MPI_COMM_WORLD) { 
			comm_ = comm;
			MPI_Comm_size(comm_, &size_);
			MPI_Comm_rank(comm_, &rank_);
		}

		MPI_Comm comm() const { return comm_; }
		int size() const { return size_; }
		int rank() const { return rank_; }

}; // class MPI_data

#endif // MPI_DATA_HPP
