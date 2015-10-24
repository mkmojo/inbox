/***
 *  Project:
 * 
 *  File: user_input.hpp
 *  Created: Apr 10, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#ifndef USER_INPUT_HPP
#define USER_INPUT_HPP

#include <iostream>

// This is specific to octree's user's input data (defines the point type)
// isolate this from this file into a user defined file, and also define the x(), y() and z() operations
struct InputData {

	struct Point {
		double x_;
		double y_;
		double z_;

		Point() { x_ = 0.0; y_ = 0.0; z_ = 0.0; }

		double x() const { return x_; }
		double y() const { return y_; }
		double z() const { return z_; }
	}; // struct Point

	Point *point_array_;    // array to store all the input data points
	unsigned int n_;        // local number of points
	unsigned int n_global_; // total number of data points
	int levels_;   // number of levels in the tree

	InputData() { point_array_ = 0; n_ = 0; n_global_ = 0; levels_ = 0; }
	~InputData() { delete[] point_array_; }


	// allocate memory for the input data
	bool input_data_alloc(unsigned int n_local, unsigned int n_global, unsigned int l) {
		delete[] point_array_;

		n_ = n_local;
		n_global_ = n_global;
		levels_ = l;
		point_array_ = new (std::nothrow) Point[n_];

		if(point_array_ == 0) return false;
		return true;
	} // input_data_alloc

	void print_data(void) {
		std::cout << "Total number of data points = " << n_global_ << std::endl;
		std::cout << "Local number of data points = " << n_ << std::endl;
		std::cout << "Number of levels required = " << levels_ << std::endl;

		for(int i = 0; i < n_; ++i) {
			std::cout << i+1 << ". (" << point_array_[i].x_ << ", " << point_array_[i].y_
					<< ", " << point_array_[i].z_ << ")" << std::endl;
		}
	} // print_data

}; // struct InputData

#endif // USER_INPUT_HPP
