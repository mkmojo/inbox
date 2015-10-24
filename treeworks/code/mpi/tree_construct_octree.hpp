/***
 *  Project: TreeWorks
 *
 *  File: tree_construct_octree.hpp
 *  Created: Mar 16, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#ifndef MPI_CONSTRUCT_TREE_OCTREE_HPP
#define MPI_CONSTRUCT_TREE_OCTREE_HPP

#include "mpi_data.hpp"
#include "input_data.hpp"
#include "sample_sort.hpp"

namespace tw
{
namespace mpi
{
	template <typename point_type, typename tree_type>
	bool tree_construct_compressed_octree(const point_type* input_points,
		const unsigned int& n_local, const int& levels, MPI_Comm comm, tree_type& t) {

		MPI_data mpi_data(comm);
		OctreeData octree_data;

		// create and initialize the octree data from input_points
		octree_data.init<point_type>(input_points, n_local, levels, mpi_data);
		octree_data.sample_sort(mpi_data);
		octree_data.count_unique();

		if(mpi_data.rank() == 0) octree_data.print_data();

		// construct the tree using the prepared data
		t.construct_c_octree(octree_data, mpi_data);

		return true;
	} // construct_tree_comp_octree()

} // namespace mpi
} // namespace tw

#endif // MPI_CONSTRUCT_TREE_OCTREE_HPP
