/***
 *  Project: TreeWorks
 *
 *  File: tree_compute.hpp
 *  Created: Mar 25, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#ifndef MPI_TREE_COMPUTE_HPP
#define MPI_TREE_COMPUTE_HPP

#include <vector>
#include <iterator>

#include "mpi_data.hpp"
#include "base_tree.hpp"
#include "specialized_generate.hpp"

#include "../jaz/sys_tools.hpp"

namespace tw
{
namespace mpi
{
	class TreeCompute {
		private:
			enum AlgorithmChoice {
				none,
				local_computation,
				no_dependency,
				upward_acc_rev,
				upward_acc_special,
				upward_acc_general,
				downward_acc_special,
				downward_acc_general,
				downward_acc_rev
			}; // enum AlgorithmChoice

			MPI_data mpi_data_;

		public:

			TreeCompute(MPI_Comm comm) { mpi_data_ = comm; }


			// the upward accumulation special version
			template<typename Tree>
			bool generate_interaction_set(const Tree& t, UpwardAccumulationGenerate generate,
					const typename Tree::tree_node &node,
					std::vector<std::vector<typename Tree::tree_node> >& interaction_sets,
					AlgorithmChoice& generate_type) {

				//std::cout << mpi_data_.rank() << ". This is the curious case "
				//	<< "of upward accumulation." << std::endl;

				typedef typename Tree::tree_node tree_node;

				std::vector<tree_node> temp_interaction_set;
				std::back_insert_iterator<std::vector<tree_node> >
						temp_interaction_set_iter(temp_interaction_set);

				generate_type = upward_acc_special;

				return generate(t, node, temp_interaction_set_iter);
			} // generate_interaction_set()


			// the downward accumulation special version
			template<typename Tree>
			bool generate_interaction_set(const Tree& t, DownwardAccumulationGenerate generate,
					const typename Tree::tree_node &node,
					std::vector<std::vector<typename Tree::tree_node> >& interaction_sets,
					AlgorithmChoice& generate_type) {

				//std::cout << mpi_data_.rank() << ". This is the curious case "
				//	<< "of downward accumulation." << std::endl;

				typedef typename Tree::tree_node tree_node;

				std::vector<tree_node> temp_interaction_set;
				std::back_insert_iterator<std::vector<tree_node> >
						temp_interaction_set_iter(temp_interaction_set);

				generate_type = downward_acc_special;

				return generate(t, node, temp_interaction_set_iter);
			} // generate_interaction_set()


			// the general version
			template<typename Tree, typename GenerateFunction>
			bool generate_interaction_set(const Tree& t, GenerateFunction generate,
					const typename Tree::tree_node &node,
					std::vector<std::vector<typename Tree::tree_node> >& interaction_sets,
					AlgorithmChoice& generate_type) {

				//std::cout << mpi_data_.rank() << ". This is the curious case of generality."
				//	<< std::endl;

				typedef typename Tree::tree_node tree_node;

				std::vector<tree_node> temp_interaction_set;
				std::back_insert_iterator<std::vector<tree_node> >
						temp_interaction_set_iter(temp_interaction_set);

				bool temp_flag = generate(t, node, temp_interaction_set_iter);
				interaction_sets.push_back(temp_interaction_set);
				
				generate_type = none;

				// test
				//std::cout << mpi_data_.rank() << "." << " "
				//	<< temp_interaction_set.size() << std::endl;
				/*std::cout << mpi_data_.rank() << ". I-set of " << node.small_cell()
					<< ", size = " << temp_interaction_set.size() << std::endl;
				for(int i = 0; i < temp_interaction_set.size(); ++ i) {
					std::cout << "  " << i << ". " << temp_interaction_set[i].is_leaf()
						<< " " << temp_interaction_set[i].level()
						<< " " << temp_interaction_set[i].small_cell()
						<< " " << temp_interaction_set[i].large_cell()
						<< " (" << temp_interaction_set[i].parent().proc_id_
						<< ", " << temp_interaction_set[i].parent().index_
						<< ") " << temp_interaction_set[i].num_children()
						<< " " << temp_interaction_set[i].num_points()
						<< std::endl;
				} // for */

				return temp_flag;
			} // generate_interaction_set()


			template <typename Tree, typename GenerateFunction, typename CombineFunction>
			bool operator()(Tree& t, GenerateFunction generate, CombineFunction combine) {

				typedef Tree tree_type;
				typedef typename Tree::tree_node tree_node;
				typedef typename Tree::index_type index_type;
				typedef typename Tree::iterator node_iterator;

				if(mpi_data_.rank() == 0) std::cout << "+ performing tree compute ... ";

				MPI_Barrier(mpi_data_.comm());

				double tree_compute_total_time_s = MPI_Wtime();
				double generate_time_s = MPI_Wtime();

				// there are t.size() nodes in the local tree
				// apply generate function to each of them
				// and obtain a list of tree_nodes for each
				std::vector<std::vector<tree_node> > interaction_sets;
				bool dependency_flag = false;
				AlgorithmChoice generate_type = none;

				// apply generate function on all nodes of the tree
				for(node_iterator iter = t.begin(); iter != t.end(); ++ iter) {
					//unsigned long int memory = jaz::mem_usage();
					//if(mpi_data_.rank() == 0)
					//	 std::cout << mpi_data_.rank() << "." << temp << " Memory used in bytes = "
					//		 << memory << std::endl;

					//if(mpi_data_.rank() == 0)
					//std::cout << mpi_data_.rank() << ". [" << (*iter).index().proc_id_ << ", "
					//	<< (*iter).index().index_ << "], parent = ("
					//	<< (*iter).parent().proc_id_ << ", "
					//	<< (*iter).parent().index_ << ")" << std::endl;

					AlgorithmChoice temp_generate_type = none;
					bool temp_flag = generate_interaction_set(t, generate, *iter, interaction_sets,
							temp_generate_type);

					// check that flag for each call to generate returns the same thing
					if(iter == t.begin()) {
						dependency_flag = temp_flag;
						generate_type = temp_generate_type;
					} else {
						if(temp_flag != dependency_flag) {
							std::cerr << "Grave error: DEPENDENCY_FLAG of all interaction sets "
								<< "are not the same! Aborting." << std::endl;
							return false;
						} // if
						if(temp_generate_type != generate_type) {
							std::cerr << "Grave error: Conflict in generate_type! Aborting."
								<< std::endl;
							return false;
						} // if
					} // if-else
				} // for

				MPI_Barrier(mpi_data_.comm());
				double generate_time_e = MPI_Wtime();
				double detection_time_s = MPI_Wtime();

				AlgorithmChoice combine_case = none;

				// identify the local computation case: each node has itself in its i-set
				// each node should have 1 node in its i-set and it shud be itself
				// (the dependency flag may be either true or false)
				if(generate_type == none) {
					int i = 0;
					for(node_iterator iter = t.begin(); iter != t.end(); ++ iter, ++ i) {
						// check if i-set sizes are == 1
						if(interaction_sets[i].size() != 1) break;
						// check if the node in i-set is itself
						if(!(interaction_sets[i][0].index() == (*iter).index())) break;
					} // for
					if(i == interaction_sets.size()) generate_type = local_computation;
				} // if

				if(dependency_flag == false) {
					// This is simple, just use each node in the interaction
					// set of each local node and perform the computations. Call this case no_dep.
					if(generate_type == local_computation) combine_case = local_computation;
					else combine_case = no_dependency;
				} else {
					// if(dependency_flag == true), then using the interaction set of each local node,
					// find a local consensus of the levels, the special cases of children/parent only,
					// and conditions for uniqueness of parent, and then find a global consensus.
					// If there is no global consensus, notify that dependency cannot be satisfied.
					// Else, the following cases occur for each local node u and remote node v:
					// 	* upward_acc_rev. each node in at most 1 i-set, and level(u) > level(v) 
					// 	* b. each node in at most 1 i-set, and level(u) < level(v)
					// 	* c. each node has at most 1 in i-set, and level(u) > level(v)
					// 	* downward_acc_rev. each node has at most 1 in i-set, and level(u) < level(v)
					// Cases upward_acc_rev and b result in upward accumulation (where in
					// upward_acc_rev, the tree is upside-down, and b is the original case of
					// upward tree accumulation when all v are u's children).
					// Cases c and downward_acc_rev result in downward accumulation (where in
					// downward_acc_rev, the tree is upside-down, and c is the original case of
					// downward tree accumulation).
					// Case b -> either one of the following:
					// 	* upward_acc_special: where each all (and only) children are in the i-set.
					// 	* upward_acc_general: the other cases.
					// Case c -> either one of the following:
					// 	* downward_acc_special: where only the parent is present in i-set.
					// 	* downward_acc_general: the other cases.

					if(generate_type == local_computation) combine_case = local_computation;
					else if(generate_type == upward_acc_special) combine_case = upward_acc_special;
					else if(generate_type == downward_acc_special) combine_case = downward_acc_special;
					else {
						// detect the cases for the special cases of upward and downward accumulations

						// check for special downward accumulation:
						if(combine_case == none) {
							//std::cout << "Checking for the special case of downward accumulation"
							//	<< std::endl;
							// for each node check if the I-set has only one node in it
							// and check about the levels of the nodes
							int i = 0;
							for(node_iterator iter = t.begin(); iter != t.end(); ++ iter, ++ i) {
								if(!(*iter).is_root()) {
									// check if i-set sizes are == 1
									if(interaction_sets[i].size() != 1) break;
									// check if the node in i-set of each node is its parent
									if(!interaction_sets[i][0].is_parent(*iter)) break;
								} // if
							} // for
							if(i == interaction_sets.size()) combine_case = downward_acc_special;
						} // if

						// check for special upward accumulation:
						if(combine_case == none) {
							//std::cout << "Checking for the special case of upward accumulation"
							//	<< std::endl;
							// for each node check if the I-set has same # of nodes as its # of children
							// and check for each node if it is its child
							int i = 0;
							for(node_iterator iter = t.begin(); iter != t.end(); ++ iter, ++ i) {
								if(!(*iter).is_leaf()) {
									// check if i-set sizes are == num of children
									if((*iter).num_children() != interaction_sets[i].size()) break;
									// check if the node in i-set of each node is its parent
									bool break_flag = false;
									for(int j = 0; j < interaction_sets[i].size(); ++j) {
										if(!interaction_sets[i][j].is_child(*iter)) {
											break_flag = true;
											break;
										} // if
									} // for
									if(break_flag) break;
								} // if
							} // for
							if(i == interaction_sets.size()) combine_case = upward_acc_special;
						} // if

						// implement algorithm detection for other general cases
						// ...

					} // if-else
				} // if-else

				MPI_Barrier(mpi_data_.comm());

				// find consensus combine_case among all procs
				AlgorithmChoice* consensus = new AlgorithmChoice[mpi_data_.size()];
				MPI_Allgather(&combine_case, 1, MPI_INT, consensus, 1, MPI_INT, mpi_data_.comm());
				for(int i = 0; i < mpi_data_.size(); ++ i) {
					if(consensus[i] != combine_case) {
						std::cerr << "Error in obtaining consensus for computations! Aborting."
							<< std::endl;
						return false;
					} // if
				} // for

				double detection_time_e = MPI_Wtime();

				// Currently only the special cases, upward_acc_special and downward_acc_special,
				// are implemented, were the nodes in the all i-sets are all children, or the parent,
				// respectively. In these cases, new dependency forest is not constructed.

				double compute_time_s = 0.0;
				double compute_time_e = 0.0;

				switch(combine_case) {
					case local_computation:
						// Local computation: apply the combine function to each node locally
						//std::cout << "Local computations case." << std::endl;
						compute_time_s = MPI_Wtime();
						t.local_compute(combine, mpi_data_);
						compute_time_e = MPI_Wtime();
						break;

					case no_dependency:
						// look into the paper dealing with this case and do the corresponding
						// implementation for special cases if needed.
						// All the nodes in the interaction set are already available at the
						// local processor, since they were required from the generate function.
						//std::cout << mpi_data_.rank() << ". No dependency case." << std::endl;
						compute_time_s = MPI_Wtime();
						t.no_dependency_compute(combine, interaction_sets, mpi_data_);
						compute_time_e = MPI_Wtime();
						break;

					case upward_acc_special:
						// The Upward Accumulation: i-set has all and only the children, for all nodes
						//std::cout << "Upward accumulation special case." << std::endl;
						compute_time_s = MPI_Wtime();
						t.upward_accumulate(combine, mpi_data_);
						compute_time_e = MPI_Wtime();
						break;

					case downward_acc_special:
						// The Downward Accumulation: i-set has only one node, and it is the
						// parent, for all nodes
						//std::cout << "Downward accumulation special case." << std::endl;
						compute_time_s = MPI_Wtime();
						t.downward_accumulate(combine, mpi_data_);
						compute_time_e = MPI_Wtime();
						break;

					case upward_acc_rev:
					case upward_acc_general:
					case downward_acc_general:
					case downward_acc_rev:
						std::cout << "Not yet implemented!" << std::endl;
						break;

					default:
						std::cerr << "Something went very wrong in algorithm detection." << std::endl;
						return false;
				} // switch

				MPI_Barrier(mpi_data_.comm());
				double tree_compute_total_time_e = MPI_Wtime();

				if(mpi_data_.rank() == 0) {
					std::cout << "done: "
						<< (tree_compute_total_time_e - tree_compute_total_time_s) * 1000 << "ms"
						<< " [g: " << (generate_time_e - generate_time_s) * 1000 << "ms"
						<< ", d: " << (detection_time_e - detection_time_s) * 1000 << "ms"
						<< ", c: " << (compute_time_e - compute_time_s) * 1000 << "ms]"
						<< std::endl;
				} // if

				// test
				/*if(mpi_data_.rank() == 0)
					t.print_tree(); */

				return true;
			} // operator()()

	}; // class TreeCompute

} // namespace mpi
} // namespace tw

#endif // MPI_TREE_COMPUTE_HPP
