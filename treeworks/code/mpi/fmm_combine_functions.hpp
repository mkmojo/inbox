/***
 *  Project:
 * 
 *  File: fmm_combine_functions.hpp
 *  Created: Sep 20, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#ifndef FMM_COMBINE_FUNCTIONS_HPP
#define FMM_COMBINE_FUNCTIONS_HPP

#include "base_tree.hpp"
#include "user_value_type.hpp"
#include "user_input.hpp"
#include <unistd.h>

class LocalMultipoleExpansionsCombineFunction {
	public:
		typedef OctreeNodeValueType value_type;
		typedef InputData::Point point_type;
		typedef tw::mpi::BaseTree<value_type, point_type>::tree_node tree_node;

		bool operator()(tree_node& u, const tree_node& v) {
			
			if(u.is_leaf()) {
				usleep(1000);	// sleep for 1ms
			} // if
			return true;
		} // operator()()
}; // class LocalMultipoleExpansionsCombineFunction

class InterpolationsCombineFunction {	// upward accumulation
	public:
		typedef OctreeNodeValueType value_type;
		typedef InputData::Point point_type;
		typedef tw::mpi::BaseTree<value_type, point_type>::tree_node tree_node;

		bool operator()(tree_node& u, const tree_node& v) {
			
			usleep(5500);	// sleep for 5.5ms

			return true;
		} // operator()()
}; // class InterpolationsCombineFunction

class TranslationsCombineFunction {		// i-list computations
	public:
		typedef OctreeNodeValueType value_type;
		typedef InputData::Point point_type;
		typedef tw::mpi::BaseTree<value_type, point_type>::tree_node tree_node;

		bool operator()(tree_node& u, const tree_node& v) {
			
			usleep(5800);	// sleep for 5.8ms

			return true;
		} // operator()()
}; // class TranslationsCombineFunction

class AnterpolationsCombineFunction {	// downward accumulation
	public:
		typedef OctreeNodeValueType value_type;
		typedef InputData::Point point_type;
		typedef tw::mpi::BaseTree<value_type, point_type>::tree_node tree_node;

		bool operator()(tree_node& u, const tree_node& v) {
			
			usleep(6000);	// sleep for 6ms

			return true;
		} // operator()()
}; // class AnterpolationsCombineFunction

class NearfieldCombineFunction {	// downward accumulation
	public:
		typedef OctreeNodeValueType value_type;
		typedef InputData::Point point_type;
		typedef tw::mpi::BaseTree<value_type, point_type>::tree_node tree_node;

		bool operator()(tree_node& u, const tree_node& v) {
			
			usleep(19000);	// sleep for 19ms

			return true;
		} // operator()()
}; // class AnterpolationsCombineFunction

#endif // FMM_COMBINE_FUNCTIONS_HPP
