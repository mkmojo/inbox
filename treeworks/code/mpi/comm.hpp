/***
 *  Project: TreeWorks
 *
 *  File: comm.hpp
 *  Created: Mar 19, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#ifndef COMM_HPP
#define COMM_HPP

#include "mpi_data.hpp"
#include "sample_sort.hpp"

namespace twlib {

	// MERGE WITH BEN'S ALL_GATHER
	// assume T is OctreeData<InputData::Point>::OctreePoint
	template <typename T>
	bool gather_cells(const T& cell, T* &array, const MPI_data& mpi_data) {

		array = new T[1];
		array[0] = cell;
		int n = 1;

		// just use ben's all gather
		all_gather(array, n, mpi_data.comm());

		return true;
	} // gather_cells()

} // namespace twlib

#endif // COMM_HPP
