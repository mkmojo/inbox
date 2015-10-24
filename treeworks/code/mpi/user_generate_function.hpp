/***
 *  Project: TreeWorks
 *
 *  File: user_generate_function.hpp
 *  Created: Mar 25, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#ifndef USER_GENERATE_FUNCTION_HPP
#define USER_GENERATE_FUNCTION_HPP

#include "base_tree.hpp"
#include "user_value_type.hpp"
#include "user_input.hpp"

class GenerateFunction {
	private:

	public:

		typedef InputData::Point point_type;
		typedef OctreeNodeValueType value_type;
		typedef tw::mpi::BaseTree<value_type, point_type> tree_type;
		typedef tree_type::index_type index_type;
		typedef tree_type::tree_node tree_node;

		// returns the dependency_flag
		bool operator()(const tree_type& t, const tree_node& local_node,
				std::back_insert_iterator<std::vector<tree_node> > &out) {

			// ...

			return false;
		} // operator()()

}; // class GenerateFunction


#endif // USER_GENERATE_FUNCTION_HPP
