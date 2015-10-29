/***
 *  Project: TreeWorks
 *
 *  File: user_value_type.hpp
 *  Created: Mar 25, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#ifndef USER_VALUE_TYPE_HPP
#define USER_VALUE_TYPE_HPP

#include <climits>

#include "input_data.hpp"
#include "user_head.hpp"

struct OctreeNodeValueType {
	double mass_;
	double velocity_;

	OctreeData::Coord nn_coord[K];	// coordinates of NNs
	double nn_dist[K];				// distances to NNs
	double d_k;		// Kth smallest distance
	int count;		// how many NNs are present


	OctreeNodeValueType() { mass_ = 0.0; velocity_ = 0.0; d_k = INT_MAX; count = 0; }
	void print() {
		std::cout << "[" << mass_ << ", " << velocity_ << "]";
	} // print()
}; // struct OctreeValueType

#endif // USER_VALUE_TYPE_HPP
