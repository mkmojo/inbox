/***
 *  Project: TreeWorks
 *
 *  File: user_combine_function.hpp
 *  Created: Mar 25, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#ifndef USER_COMBINE_FUNCTION_HPP
#define USER_COMBINE_FUNCTION_HPP

#include "base_tree.hpp"
#include "user_value_type.hpp"
#include "user_input.hpp"

class UserCombineFunction {
	private:

	public:
		typedef OctreeNodeValueType value_type;
		typedef InputData::Point point_type;
		typedef tw::mpi::BaseTree<value_type, point_type>::tree_node tree_node;

		bool operator()(tree_node& u, const tree_node& v) {

			// dummy
			// ...

			return true;
		} // operator()()

}; // class UserCombineFunction


class LocalCombineFunction {
	private:

	public:
		typedef OctreeNodeValueType value_type;
		typedef InputData::Point point_type;
		typedef tw::mpi::BaseTree<value_type, point_type>::tree_node tree_node;

		bool operator()(tree_node& u, const tree_node& v) {

			if(u.is_leaf()) {
				value_type temp_value;
				temp_value.mass_ = u.num_points();
				temp_value.velocity_ = rand() % 100;
				u.set_value(temp_value);
			} // if

			return true;
		} // operator()()

}; // class LocalCombineFunction


class UpwardCombineFunction {
	private:

	public:
		typedef OctreeNodeValueType value_type;
		typedef InputData::Point point_type;
		typedef tw::mpi::BaseTree<value_type, point_type>::tree_node tree_node;

		bool operator()(tree_node& u, const tree_node& v) {

			value_type temp_value = u.value();
			temp_value.mass_ += v.value().mass_;
			if(temp_value.velocity_ < v.value().velocity_)
				temp_value.velocity_ = v.value().velocity_;
			u.set_value(temp_value);

			return true;
		} // operator()()

}; // class UpwardCombineFunction


class DownwardCombineFunction {
	private:

	public:
		typedef OctreeNodeValueType value_type;
		typedef InputData::Point point_type;
		typedef tw::mpi::BaseTree<value_type, point_type>::tree_node tree_node;

		bool operator()(tree_node& u, const tree_node& v) {

			value_type temp_value = u.value();
			temp_value.velocity_ = v.value().velocity_;
			u.set_value(temp_value);

			return true;
		} // operator()()

}; // class DownwardCombineFunction


#endif // USER_COMBINE_FUNCTION_HPP
