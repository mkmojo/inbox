/***
 *  Project: TreeWorks
 * 
 *  File: utility.hpp
 *  Created: May 13, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <cstdint>
//#include <stdint.h>
#include <cmath>

template<typename T>
void swap(T& a, T& b) {
	// ...
} // swap()

// compare doubles
bool is_equal(double a, double b) {
	if(fabs(a - b) > 1e-9) return false;
	return true;
} // is_equal()

// returns the position of leftmost 1
int log_2(uint64_t a) {
	int k = 0;

	while(a) {
		k++;
		a = a >> 1;
	}
                
	return k;
} // log_2()

// function to find if a is contained in b (strictly)
bool is_contained(uint64_t a, uint64_t b) {
	// get the positions of the leftmost 1
	int k = log_2(a);
	int l = log_2(b);

	if(a == 0 || b == 0) return false;

	if(b == 1) return true; // b is the whole domain

	if(l >= k) return false;    // cell a cannot have larger or equal sized key as b

	uint64_t temp1 = 0, temp2 = 0, temp3 = 0;
	temp1 = a << (64 - k);
	temp2 = b << (64 - l);
	temp3 = temp1 ^ temp2;

	int m = log_2(temp3);

	if(m > 64 - l) return false;

	return true;
} // is_contained()

// compute the large cell of a node given its small cell 'cell'
// and its parent's small cell 'parent'.
// The large cell is the cell obtained by first division of parent's
// small cell, that contains the node's small cell.
uint64_t compute_large_cell(uint64_t parent, uint64_t cell) {
	int k = log_2(parent);
	int l = log_2(cell);

	uint64_t large_cell = 0;

	large_cell = cell >> (l - k - 3);

	return large_cell;
} // compute_large_cell()

// Compute the smallest cell which contains the two given cells a and b
uint64_t compute_small_cell(uint64_t a, uint64_t b) {

	if(is_contained(a, b)) return b;
	if(is_contained(b, a)) return a;

	int k = log_2(a);
	int l = log_2(b);

	uint64_t temp1 = 0, temp2 = 0, temp3 = 0, small_cell = 0;
	temp1 = a << (64 - k);
	temp2 = b << (64 - l);

	temp3 = temp1 ^ temp2;

	int m = log_2(temp3);
	int shift = 63 - ((63 - m) / 3) * 3;

	small_cell = temp1 >> shift;

	return small_cell;
} // compute_small_cell()

// compare two cells according to postordering
// true if a < b
bool compare(const uint64_t& a, const uint64_t& b) {
	// if either is contained in the other, return the order
	if( is_contained(a, b) ) return true;
	if( is_contained(b, a) ) return false;

	int k = log_2(a);
	int l = log_2(b);

	if(k > l) {
		uint64_t low = a >> (k - l);
		if(low < b) return true;
		else return false;
	} else {
		if(l > k) {
			uint64_t low = b >> (l - k);
			if(a < low) return true;
			else return false;
		} else {
			if(a < b) return true;
			else return false;
		} // if
	} // if
} // compare()

// given a cell and domain size, compute the coordinates of its corners
bool compute_coordinates(const uint64_t in_cell, const double domain,
		double& x1, double& x2, double& y1, double& y2, double& z1, double& z2) {

	int cell_bits = log_2(in_cell);

	uint64_t cell = in_cell << (64 - cell_bits + 1);
	int num_bits = (cell_bits - 1) / 3;

	unsigned int x = 0, y = 0, z = 0;
	uint64_t mask = 1;
	mask = mask << 63;

	for(int i = 0; i < num_bits; ++ i) {
		x = x << 1;
		if(cell & mask) x = x | 1;
		cell = cell << 1;

		y = y << 1;
		if(cell & mask) y = y | 1;
		cell = cell << 1;

		z = z << 1;
		if(cell & mask) z = z | 1;
		cell = cell << 1;
	} // for

	mask = 1 << num_bits;
	double sidelength = domain / mask;

	x1 = x * sidelength;
	y1 = y * sidelength;
	z1 = z * sidelength;
	x2 = (x + 1) * sidelength;
	y2 = (y + 1) * sidelength;
	z2 = (z + 1) * sidelength;

	return true;
} // compute_coordinates()

// compute eucledian distance between two given points
double distance(const double a_x, const double a_y, const double a_z,
		const double b_x, const double b_y, const double b_z) {

	return (sqrt(pow((a_x - b_x), 2) + pow((a_y - b_y), 2) + pow((a_z - b_z), 2)));

} // distance()

#endif // UTILITY_HPP
