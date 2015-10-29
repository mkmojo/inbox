/***
 *  Project: TreeWorks
 *
 *  File: user_test.cpp
 *  Created: Mar 16, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#include <iostream>
#include <mpi.h>

// the library files
#include "base_tree.hpp"
#include "tree_compute.hpp"
#include "specialized_generate.hpp"

// user defined files
#include "mpi_data.hpp"
#include "user_input.hpp"
#include "user_generate_function.hpp"
#include "user_combine_function.hpp"
#include "user_value_type.hpp"
#include "user_read_input_data.hpp"

// generate and combine functons for KNN
#include "knn_generate.hpp"
#include "knn_combine.hpp"

// generate and combine functions for FMM
#include "adjacent_nodes_generate_function.hpp"
#include "interaction_list_generate_function.hpp"
#include "fmm_combine_functions.hpp"


int main(int argc, char** argv)
{
	if(argc != 2) {
		std::cerr << "Usage: a.out <filename>" << std::endl;
		return 1;
	}

	MPI_Init(&argc, &argv);

	MPI_data mpi_data(MPI_COMM_WORLD);

	// prepare input data
	InputData my_input_data;
	if(!read_point_data(argv[1], mpi_data, my_input_data)) {
		std::cerr << "Error in reading input data." << std::endl;
		MPI_Finalize();
		return 1;
	}

	typedef tw::mpi::BaseTree<OctreeNodeValueType, InputData::Point> octree_type;

	octree_type my_tree(my_input_data.point_array_, my_input_data.n_,
			my_input_data.levels_, mpi_data.comm());

	//UserCombineFunction dummy_combine;

	// test the functions
/*	tw::mpi::TreeCompute demo_compute_one(mpi_data.comm());

	// set the initial mass and velocity
	LocalComputationGenerate local_gen;
	LocalCombineFunction local_combine;
	demo_compute_one(my_tree, local_gen, local_combine);

	// perform an upward accumulation to compute total mass for each node
	UpwardAccumulationGenerate2 upward_acc_gen2;
	UpwardCombineFunction upward_combine;
	demo_compute_one(my_tree, upward_acc_gen2, upward_combine);

	// compute new velocities for each node in downward accumulation
	DownwardAccumulationGenerate2 downward_acc_gen2;
	DownwardCombineFunction downward_combine;
	demo_compute_one(my_tree, downward_acc_gen2, downward_combine);

	UpwardAccumulationGenerate upward_acc_gen;
	demo_compute_one(my_tree, upward_acc_gen, upward_combine);

	DownwardAccumulationGenerate downward_acc_gen;
	demo_compute_one(my_tree, downward_acc_gen, downward_combine);
*/

	/**
	 * KNN
	 */

/*	tw::mpi::TreeCompute knn_compute(mpi_data.comm());
	KNNGenerate knn_generate;
	KNNCombine knn_combine;
	demo_compute_one(my_tree, knn_generate, knn_combine);
*/

	/**
	 * FMM
	 */

	tw::mpi::TreeCompute fmm_compute(mpi_data.comm());

	// 1. local multipole expansions (local computations)
	LocalComputationGenerate local_gen;
	LocalMultipoleExpansionsCombineFunction local_multipole_combine;
	fmm_compute(my_tree, local_gen, local_multipole_combine);

	// 2. interpolations (upward accumulations)
	UpwardAccumulationGenerate2 upward_acc_gen2;
	InterpolationsCombineFunction interpolations_combine;
	fmm_compute(my_tree, upward_acc_gen2, interpolations_combine);

	// 3. translations (i-list computations)
	InteractionListGenerate ilist_generate;
	TranslationsCombineFunction translations_combine;
	fmm_compute(my_tree, ilist_generate, translations_combine);

	// 4. anterpolations (downward accumulations)
	DownwardAccumulationGenerate2 downward_acc_gen2;
	AnterpolationsCombineFunction anterpolations_combine;
	fmm_compute(my_tree, downward_acc_gen2, anterpolations_combine);

	// 5. nearfields (nearfield computations)
	AdjacentNodesGenerate nearfield_generate;
	NearfieldCombineFunction nearfield_combine;
	fmm_compute(my_tree, nearfield_generate, nearfield_combine);

	// 6. total fields (summations)
	LocalCombineFunction local_combine;
	fmm_compute(my_tree, local_gen, local_combine);

	// always call finalize to the tree before finalizing MPI
	my_tree.finalize();

	MPI_Finalize();
	
	return 0;
} // main()
