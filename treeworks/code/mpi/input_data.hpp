/***
 *  Project: TreeWorks
 *
 *  File: input_data.hpp
 *  Created: Mar 17, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#ifndef INPUT_DATA_HPP
#define INPUT_DATA_HPP

#include <iostream>

#include "mpi_data.hpp"
#include "sample_sort.hpp"

inline int min(int a, int b) { return ((a < b) ? a : b); }

// Octree specific data (construction user)
class OctreeData {

	public:

		typedef struct Coord {
			double x_;
			double y_;
			double z_;
			
			// add a user defined value type for the points
			// ...

			Coord() { x_ = 0.0; y_ = 0.0; z_ = 0.0; }

			// add functions to obtain and set the value field
			// ...

		} coord_type;

		typedef struct OctreePoint {
			coord_type data_point_;
			unsigned int x_key_;	// the x, y and z keys
			unsigned int y_key_;
			unsigned int z_key_;
			uint64_t cell_key_; // the key (cell) to which the point belongs (64b)

			bool operator<(const OctreePoint& b) {
				uint64_t key_a = 0, key_b = 0;

				key_a = cell_key_;
				key_b = b.cell_key_;

				return (key_a < key_b);
			} // operator<

			bool operator>(const OctreePoint& b) {
				uint64_t key_a = 0, key_b = 0;

				key_a = cell_key_;
				key_b = b.cell_key_;

				return (key_a > key_b);
			} // operator>

			bool operator==(const OctreePoint& b) {
				uint64_t key_a = 0, key_b = 0;

				key_a = cell_key_;
				key_b = b.cell_key_;

				return (key_a == key_b);
			} // operator==

			bool operator!=(const OctreePoint& b) {
				uint64_t key_a = 0, key_b = 0;

				key_a = cell_key_;
				key_b = b.cell_key_;

				return (key_a != key_b);
			} // operator!=

			bool operator<=(const OctreePoint& b) {
				uint64_t key_a = 0, key_b = 0;

				key_a = cell_key_;
				key_b = b.cell_key_;

				return (key_a <= key_b);
			} // operator<=

			bool operator>=(const OctreePoint& b) {
				uint64_t key_a = 0, key_b = 0;

				key_a = cell_key_;
				key_b = b.cell_key_;

				return (key_a >= key_b);
			} // operator>=

		} octree_point; // struct OctreePoint

	private:

		OctreePoint *octree_point_array_;
		unsigned int n_; 		// local number of points
		unsigned int n_global_;	// total number of points
		unsigned int levels_;	// number of levels in the tree
		double domain_;			// the domain of the points
		unsigned int n_unique_; // local number of unique cells

	public:

		OctreeData() { octree_point_array_ = 0; n_ = 0; n_global_ = 0; levels_ = 0;
				domain_ = 0.0; n_unique_ = 0; }

		~OctreeData() { if(octree_point_array_ != 0) delete[] octree_point_array_; }

		OctreePoint operator[](unsigned int i) const {
			return octree_point_array_[i];
		} // operator[]

		unsigned int n() const { return n_; }
		unsigned int n_global() const { return n_global_; }
		unsigned int levels() const { return levels_; }
		unsigned int n_unique() const { return n_unique_; }
		double domain() const { return domain_; }

		// The input would be an array of the point_type (from the user) and n_local, levels
		template<typename point_type>
		bool init(const point_type* &input_points,
				const unsigned int& n_local, const int& levels,
				const MPI_data& mpi_data) {

			n_ = n_local;
			levels_ = levels;

			// compute the total number of points
			MPI_Allreduce(&n_, &n_global_, 1, MPI_UNSIGNED, MPI_SUM, mpi_data.comm());

			// array size is 3n to make sure there is enough space for sorted points
			octree_point_array_ = new /*(std::nothrow)*/ OctreePoint[3 * n_];
			if(octree_point_array_ == 0) {
				std::cerr << "Error in allocating memory for octree_point_array_" << std::endl;
				return false;
			} // if

			for(unsigned int i = 0; i < n_ ; ++i) {
				octree_point_array_[i].data_point_.x_ = input_points[i].x();
				octree_point_array_[i].data_point_.y_ = input_points[i].y();
				octree_point_array_[i].data_point_.z_ = input_points[i].z();

				// also put the user defined values associated with the points
				// ...
			}

			// translate each point to the first quadrant
			translate(mpi_data);

			// compute the cell each point belongs to and its corresponding key
			compute_cells(mpi_data);

			return true;
		} // init()


		bool count_unique() {
			n_unique_ = 1;
			for(unsigned int i = 0; i < n_ - 1; ++ i) {
				if(octree_point_array_[i].cell_key_ != octree_point_array_[i+1].cell_key_)
					++ n_unique_;
			} // for

			return true;
		} // count_unique()


		bool sample_sort(const MPI_data& mpi_data) {
			int n = n_;
			twlib::sample_sort<OctreePoint>(octree_point_array_, n, mpi_data.comm());
			n_ = n;

			return true;
		} // sample_sort()


		bool copy_octree_point(OctreePoint& dest, const OctreePoint& src) {

			dest.data_point_.x_ = src.data_point_.x_;
			dest.data_point_.y_ = src.data_point_.y_;
			dest.data_point_.z_ = src.data_point_.z_;
			dest.x_key_ = src.x_key_;
			dest.y_key_ = src.y_key_;
			dest.z_key_ = src.z_key_;
			dest.cell_key_ = src.cell_key_;

			return true;
		} // copy_octree_point()


		bool load_balance(const MPI_data& mpi_data) {

			MPI_Barrier(mpi_data.comm());

			// obtain number of points on all processors
			unsigned int * n_curr_array = new unsigned int[mpi_data.size()];
			MPI_Allgather(&n_, 1, MPI_UNSIGNED, n_curr_array, 1, MPI_UNSIGNED, mpi_data.comm());

			unsigned int n_local = 0;
			if(n_global_ % mpi_data.size() == 0) {
				n_local = n_global_ / mpi_data.size();
			} else {
				unsigned int rem = n_global_ % mpi_data.size();
				if(mpi_data.rank() < rem) n_local = (n_global_ / mpi_data.size()) + 1;
				else n_local = n_global_ / mpi_data.size();
			} // if-else

			// obtain supposed to be number of points on all processors (used to search for proc id)
			unsigned int * n_tobe_array = new unsigned int[mpi_data.size()];
			MPI_Allgather(&n_local, 1, MPI_UNSIGNED, n_tobe_array, 1, MPI_UNSIGNED, mpi_data.comm());

			for(int i = 1; i < mpi_data.size(); ++ i) n_tobe_array[i] += n_tobe_array[i - 1];
			// n_tobe_array now contains the last global index + 1 for each proc

			// find global current start and end indices for local n
			unsigned int curr_start = 0, curr_end = 0;
			for(int i = 0; i < mpi_data.rank(); ++ i) curr_start += n_curr_array[i];
			curr_end = curr_start + n_ - 1;

			// find the procs where start and end should go
			int proc_start = -1, proc_end = -1;
			for(int i = 0; i < mpi_data.size(); ++ i) {
				if(curr_start < n_tobe_array[i]) {
					proc_start = i;
					break;
				} // if
			} // for
			for(int i = 0; i < mpi_data.size(); ++ i) {
				if(curr_end < n_tobe_array[i]) {
					proc_end = i;
					break;
				} // if
			} // for

			// compute send counts
			int* send_count = new int[mpi_data.size()];
			for(int i = 0; i < mpi_data.size(); ++ i) send_count[i] = 0;
			for(int i = 0, p = proc_start, j = curr_start; i < n_;) {
				send_count[p] = (min(curr_end, n_tobe_array[p] - 1) - j + 1);
				j += send_count[p];
				i += send_count[p];
				++ p;
			} // for

			int* recv_count= 0;
			int recv_size;
			OctreePoint* recv_buff = 0;
			twlib::all_to_all(send_count, octree_point_array_, recv_count, recv_buff,
					recv_size, mpi_data.comm());

			// send last cell_key_ to next processor
			// recv last cell_key_ from prev processor
			// if my first cell_key_ == the recv cell_key_
			//     send all my points with same cell_key_ to prev proc

			uint64_t last_cell = recv_buff[recv_size - 1].cell_key_;
			uint64_t recv_last_cell = 0;
			MPI_Status status;
			if(mpi_data.rank() != mpi_data.size() - 1)
				MPI_Send(&last_cell, 8, MPI_BYTE, mpi_data.rank() + 1, 0, mpi_data.comm());
			if(mpi_data.rank() != 0)
				MPI_Recv(&recv_last_cell, 8, MPI_BYTE, mpi_data.rank() - 1, 0, mpi_data.comm(), &status);

			int count = 0;
			for(int i = 0; i < mpi_data.size(); ++ i) send_count[i] = 0;
			if(mpi_data.rank() != 0) {
				// count the number of points to send to prev proc
				for(int i = 0; recv_buff[i].cell_key_ == recv_last_cell && i < recv_size; ++ i) {
					++ send_count[mpi_data.rank() - 1];
					++ count;
				} // for
			} // if
			int size = (count == 0) ? 1 : count;
			OctreePoint* temp_send_buff = new OctreePoint[size];
			for(int i = 0; i < count; ++ i) temp_send_buff[i] = recv_buff[i];
			int* temp_recv_count= 0;
			int temp_recv_size;
			OctreePoint* temp_recv_buff = 0;
			twlib::all_to_all(send_count, temp_send_buff, temp_recv_count, temp_recv_buff,
					temp_recv_size, mpi_data.comm());

			MPI_Barrier(mpi_data.comm());

			OctreePoint* temp_point = octree_point_array_;
			octree_point_array_ = 0;

			// now construct the local octree_point_array_

			n_ = recv_size - count + temp_recv_size;
			octree_point_array_ = new OctreePoint[n_];

			int i = 0;
			for(; i < recv_size - count; ++ i)
				octree_point_array_[i] = recv_buff[count + i];
			for(int j = 0; i < n_; ++ i, ++ j)
				octree_point_array_[i] = temp_recv_buff[j];

			delete[] temp_recv_buff;
			delete[] temp_recv_count;
			delete[] temp_send_buff;
			delete[] recv_buff;
			delete[] recv_count;
			delete[] send_count;
			delete[] n_tobe_array;
			delete[] n_curr_array;

			delete[] temp_point;

			return true;
		} // load_balance()


		void print_data(void) {
			std::cout << "Total number of data points = " << n_global_ << std::endl;
			std::cout << "Local number of data points = " << n_ << std::endl;
			std::cout << "Number of levels required = " << levels_ << std::endl;
			std::cout << "Domain size = " << domain_ << std::endl;

			for(int i = 0; i < n_; ++i) {
				std::cout << i+1 << ". (" << octree_point_array_[i].data_point_.x_
						<< ", " << octree_point_array_[i].data_point_.y_ << ", "
						<< octree_point_array_[i].data_point_.z_ << ") -> ";
				std::cout << "(" << octree_point_array_[i].x_key_ << ", "
						<< octree_point_array_[i].y_key_ << ", " << octree_point_array_[i].z_key_
						<< ") -> " << octree_point_array_[i].cell_key_ << std::endl;
			}
		} // print_data()


	private:

		bool translate(const MPI_data& mpi_data) {
			coord_type min_coord = find_min_coords(mpi_data);

			for(unsigned int i = 0; i < n_; ++i) {
				octree_point_array_[i].data_point_.x_ -= min_coord.x_;
				octree_point_array_[i].data_point_.y_ -= min_coord.y_;
				octree_point_array_[i].data_point_.z_ -= min_coord.z_;
			}

			return true;
		} // translate()


		bool compute_cells(const MPI_data& mpi_data) {
			
			domain_ = compute_domain(mpi_data);
			
			unsigned int num_cells = 1 << (levels_ - 1);
			double cell_size = domain_ / num_cells;

			for(unsigned int i = 0; i < n_; ++i) {
				octree_point_array_[i].x_key_ =
					(unsigned int) (octree_point_array_[i].data_point_.x_ / cell_size);
				if(octree_point_array_[i].x_key_ >= num_cells) --octree_point_array_[i].x_key_;

				octree_point_array_[i].y_key_ =
					(unsigned int) (octree_point_array_[i].data_point_.y_ / cell_size);
				if(octree_point_array_[i].y_key_ >= num_cells) --octree_point_array_[i].y_key_;

				octree_point_array_[i].z_key_ =
					(unsigned int) (octree_point_array_[i].data_point_.z_ / cell_size);
				if(octree_point_array_[i].z_key_ >= num_cells) --octree_point_array_[i].z_key_;

				octree_point_array_[i].cell_key_ = compute_key(octree_point_array_[i].x_key_,
							octree_point_array_[i].y_key_, octree_point_array_[i].z_key_);
			} // for

			return true;
		} // compute_cells()


		uint64_t compute_key(unsigned int x, unsigned int y, unsigned int z) {
			// number of bits required for each x, y and z = height of the tree = levels_ - 1
			// therefore, height of more than 21 is not supported
			unsigned int height = levels_ - 1;

			uint64_t key = 1;
			uint64_t x_64 = (uint64_t) x << (64 - height);
			uint64_t y_64 = (uint64_t) y << (64 - height);
			uint64_t z_64 = (uint64_t) z << (64 - height);
			uint64_t mask = 1;
			mask = mask << 63; // leftmost bit is 1, rest 0

			for(unsigned int i = 0; i < height; ++i) {
				key = key << 1;
				if(x_64 & mask) key = key | 1;
				x_64 = x_64 << 1;

				key = key << 1;
				if(y_64 & mask) key = key | 1;
				y_64 = y_64 << 1;

				key = key << 1;
				if(z_64 & mask) key = key | 1;
				z_64 = z_64 << 1;
			} // for

			return key;
		} // compute_key()


		double compute_domain(const MPI_data& mpi_data) {
			coord_type max_coord = find_max_coords(mpi_data);

			double domain = 0.0;
			double temp_domain = (max_coord.x_ > max_coord.y_ && max_coord.x_ > max_coord.z_) ?
					max_coord.x_ :
					((max_coord.y_ > max_coord.z_) ? max_coord.y_ : max_coord.z_ );

			MPI_Allreduce(&temp_domain, &domain, 1, MPI_DOUBLE, MPI_MAX, mpi_data.comm());

			return domain;
		} // compute_domain()


		coord_type find_min_coords(const MPI_data& mpi_data) {
			coord_type temp_min_coord = octree_point_array_[0].data_point_;
			coord_type min_coord;

			for(unsigned int i = 1; i < n_; ++i) {
				if(octree_point_array_[i].data_point_.x_ < temp_min_coord.x_)
					temp_min_coord.x_ = octree_point_array_[i].data_point_.x_;
				if(octree_point_array_[i].data_point_.y_ < temp_min_coord.y_)
					temp_min_coord.y_ = octree_point_array_[i].data_point_.y_;
				if(octree_point_array_[i].data_point_.z_ < temp_min_coord.z_)
					temp_min_coord.z_ = octree_point_array_[i].data_point_.z_;
			}

			// again, assuming that the coordinates are doubles
			MPI_Allreduce(&temp_min_coord.x_, &min_coord.x_, 1, MPI_DOUBLE,
					MPI_MIN, mpi_data.comm());
			MPI_Allreduce(&temp_min_coord.y_, &min_coord.y_, 1, MPI_DOUBLE,
					MPI_MIN, mpi_data.comm());
			MPI_Allreduce(&temp_min_coord.z_, &min_coord.z_, 1, MPI_DOUBLE,
					MPI_MIN, mpi_data.comm());

			return min_coord;
		} // find_min_coords()


		coord_type find_max_coords(const MPI_data& mpi_data) {
			coord_type temp_max_coord = octree_point_array_[0].data_point_;
			coord_type max_coord;

			for(unsigned int i = 1; i < n_; ++i) {
				if(octree_point_array_[i].data_point_.x_ > temp_max_coord.x_)
					temp_max_coord.x_ = octree_point_array_[i].data_point_.x_;
				if(octree_point_array_[i].data_point_.y_ > temp_max_coord.y_)
					temp_max_coord.y_ = octree_point_array_[i].data_point_.y_;
				if(octree_point_array_[i].data_point_.z_ > temp_max_coord.z_)
					temp_max_coord.z_ = octree_point_array_[i].data_point_.z_;
			}

			// again, assuming that the coordinates are doubles
			MPI_Allreduce(&temp_max_coord.x_, &max_coord.x_, 1, MPI_DOUBLE,
					MPI_MAX, mpi_data.comm());
			MPI_Allreduce(&temp_max_coord.y_, &max_coord.y_, 1, MPI_DOUBLE,
					MPI_MAX, mpi_data.comm());
			MPI_Allreduce(&temp_max_coord.z_, &max_coord.z_, 1, MPI_DOUBLE,
					MPI_MAX, mpi_data.comm());

			return max_coord;
		} // find_max_coords()

}; // struct OctreeData

#endif // INPUT_DATA_HPP
