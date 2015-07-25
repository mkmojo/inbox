/***
 *  Project: TreeWorks
 *
 *  File: base_tree.hpp
 *  Created: Mar 15, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#ifndef MPI_BASE_TREE_HPP
#define MPI_BASE_TREE_HPP

#include <vector>
#include <map>
#include <iterator>
#include <utility>

#include "utility.hpp"
#include "comm.hpp"
#include "input_data.hpp"
#include "sample_sort.hpp"

// jaz library for memory usage tool
#include "../jaz/sys_tools.hpp"

namespace tw
{
namespace mpi
{
	const int MAX_CHILDREN = 8;
	const int MAX_POINTS = 5;

	struct IndexType {
		int proc_id_;
		int index_;

		IndexType() { proc_id_ = -1; index_ = -1; }
		bool operator==(const IndexType& index) {
			return (index.proc_id_ == proc_id_ && index.index_ == index_);
		}
	}; // struct IndexType

	// this is specific to compressed octrees
	struct TreeTypeData {
		uint64_t small_cell_;
		uint64_t large_cell_;
	}; // struct TreeTypeData

	// -----------------------------------------
	// a structure for communication and serialization purposes only,
	// to avoid vector usage and have a contiguous memory!!!
	// is there a better way to do communication and serialization of tree nodes???
	// -----------------------------------------
	template <typename ValueType>
	struct SerialTreeNode {
		typedef ValueType value_type;
		typedef IndexType index_type;
		typedef TreeTypeData tree_type_data;
		typedef OctreeData::OctreePoint octree_point;

		bool is_leaf_;
		bool is_root_;
		int level_;
		index_type index_;
		index_type parent_;
		int num_children_;						// number of children
		index_type children_[MAX_CHILDREN];		// HARD CODED FOR OCTREES !!!
		tree_type_data ttd_;
		value_type value_;			// FOR NOW THIS IS A POD, SO JUST COPY IT !!!
		int num_points_;			// number of points (if leaf)
		octree_point point_list_[MAX_POINTS];	// list of points
	}; // struct SerialTreeNode
	// -----------------------------------------

	struct ChildData {
		int parent_index_;
		int proc_id_;
		int index_;

		ChildData() { parent_index_ = -1; proc_id_ = 0; index_ = -1; }
	}; // struct ChildData

	// The base tree class
	template <typename ValueType, typename PointType>
	class BaseTree {

		public:

			typedef IndexType index_type;
			typedef ValueType value_type;
			typedef PointType point_type;	// the user input point for octree
			typedef TreeTypeData tree_type_data;
			typedef OctreeData::OctreePoint octree_point;
			typedef OctreeData::coord_type data_point;

			struct ltindex
			{
				bool operator()(const uint64_t a, const uint64_t b) const {
					return (a < b);
				}
			};

			// The tree node class
			class TreeNode {

				public:
					// Types:
					// children iterators
					typedef std::vector<index_type>::iterator children_iterator;
					typedef std::vector<index_type>::const_iterator const_children_iterator;
					// point list iterators
					typedef std::vector<octree_point>::iterator point_iterator;
					typedef std::vector<octree_point>::const_iterator const_point_iterator;

					// constructor
					TreeNode() {
						is_leaf_ = false; is_root_ = false; level_ = 1;
						ttd_.small_cell_ = 0; ttd_.large_cell_ = 0;
					} // TreeNode()

					// copy constructor
					TreeNode(const TreeNode& node) {

						typedef const_children_iterator children_iter;
						typedef const_point_iterator point_iter;

						is_leaf_ = node.is_leaf();
						is_root_ = node.is_root();
						level_ = node.level();
						value_ = node.value(); // value is POD !!!!
						index_ = node.index();
						parent_index_ = node.parent();
						ttd_.small_cell_ = node.small_cell();
						ttd_.large_cell_ = node.large_cell();
	
						if(node.num_children() != 0) {
							std::pair<children_iter, children_iter> children = node.children();
							children_index_.clear();
							for(children_iter iter = children.first; iter != children.second; ++iter)
								children_index_.push_back(*iter);
						} // if
	
						if(node.num_points() != 0) {
							std::pair<point_iter, point_iter> points = node.points(); 
							point_list_.clear();
							for(point_iter iter = points.first; iter != points.second; ++ iter)
								point_list_.push_back(*iter);
						} // if
					} // TreeNode()

					// children access iterators
					std::pair<const_children_iterator, const_children_iterator> children() const {
						const_children_iterator begin = children_index_.begin();
						const_children_iterator end = children_index_.end();
						return std::pair<const_children_iterator,
							   const_children_iterator>(begin, end);
					} // children()

					// other access functions
					value_type value() const { return value_; }
					index_type parent() const { return parent_index_; }
					int num_children() const { return children_index_.size(); }
					index_type index() const { return index_; }
					bool is_leaf() const { return is_leaf_; }
					bool is_root() const { return is_root_; }
					int level() const { return level_; }

					// check if this node is parent of a given node
					bool is_parent(const TreeNode& node) {
						return (index_ == node.parent());
					} // is_parent()
					// check if this node is a child of a given node
					// change this to parent checking after fixing the child issue in construction ...
					bool is_child(const TreeNode& node) {
						typedef const_children_iterator children_iter;
						std::pair<children_iter, children_iter> children = node.children();
						for(children_iter iter = children.first; iter != children.second; ++ iter) {
							if(index_ == (*iter)) return true;
						} // for
						return false;
					} // is_child()


					// functions to update data at the node
					void set_value(value_type val) { value_ = val; }	// value is POD
					void set_index(index_type i) { index_ = i; }
					void set_parent(index_type i) { parent_index_ = i; }
					void set_parent_proc(int p) { parent_index_.proc_id_ = p; }
					void set_parent_index(int i) { parent_index_.index_ = i; }

					// children functions
					void add_child(index_type i) { children_index_.push_back(i); }
					void clear_children() { children_index_.clear(); }
					void remove_last_children(int count) {
						for(int i = 0; i < count; ++ i) children_index_.pop_back();
					} // remove_last_children()
					void update_last_child(index_type i) {
						children_index_.pop_back(); 
						children_index_.push_back(i);
					} // update_last_child()
					void set_child_index(const int i, int new_index) {
						children_index_[i].index_ = new_index;
					} // set_child_index()
					int remove_remote_children(int rank) {
						int count = 0;
						children_iterator iter = children_index_.begin();
						while(iter != children_index_.end()) {
							if((*iter).proc_id_ != rank) {
								iter = children_index_.erase(iter);
								++ count;
							} else ++ iter;
						} // while
						return count;
					} // remove_remote_children()
					int remove_local_children(int rank) {
						int count = 0;
						children_iterator iter = children_index_.begin();
						while(iter != children_index_.end()) {
							if((*iter).proc_id_ == rank) {
								iter = children_index_.erase(iter);
								++ count;
							} else ++ iter;
						} // while
						return count;
					} // remove_local_children()

					// misc functions
					void make_leaf() { is_leaf_ = true; }
					void make_root() { is_root_ = true; }
					void set_level(int l) { level_ = l; }

					// the following assume it is a compressed octree with large and small cells
					uint64_t large_cell() const { return ttd_.large_cell_; }
					uint64_t small_cell() const { return ttd_.small_cell_; }
					void set_large_cell(uint64_t cell) { ttd_.large_cell_ = cell; }
					void set_small_cell(uint64_t cell) { ttd_.small_cell_ = cell; }

					// DEFINE A ITERATOR TO ITERATE THROUGH THE point_list_, BUT RETURNING
					// ONLY THE USER DEFINED Point STORED AT EACH OctreePoint
					// ...

					// data point functions (for octrees)
					int num_points() const { return point_list_.size(); }
					void clear_points() { point_list_.clear(); }
					std::pair<const_point_iterator, const_point_iterator> points() const {
						const_point_iterator begin = point_list_.begin();
						const_point_iterator end = point_list_.end();
						return std::pair<const_point_iterator, const_point_iterator>(begin, end);
					} // points()
					bool add_point(octree_point p) {
						if(is_leaf()) {
							point_list_.push_back(p);
							return true;
						}
						return false;
					} // add_point()

					// operators for comparison of tree nodes with other nodes and cells
					bool operator<(const TreeNode& b) const { return compare(b); }
					bool operator>=(const TreeNode& b) const { return (!compare(b)); }
					bool operator==(const TreeNode& b) const {
						return (ttd_.small_cell_ == b.small_cell()); }
					bool operator<(const uint64_t& cell) const { return compare(cell); }
					bool operator>(const uint64_t& cell) const {
						if(!compare(cell) && ttd_.small_cell_ != cell) return true;	}
					bool operator==(const uint64_t& cell) const {
						return (ttd_.small_cell_ == cell); }

					// other miscellaneous functions for the octree node
					// check if this node is an ancestor of the given 'node'
					bool is_ancestor(const TreeNode& node) {
						return is_contained(node.small_cell(), ttd_.small_cell_);
					} // is_ancestor()

				private:
					// compare the tree_node with a cell according to postordering
					// true if this < cell
					bool compare(const uint64_t& cell) const {

						// if either is contained in the other, return the order
						if( is_contained(ttd_.small_cell_, cell) ) return true;
						if( is_contained(cell, ttd_.small_cell_) ) return false;

						int k = log_2(ttd_.small_cell_);
						int l = log_2(cell);

						if(k > l) {
							uint64_t low = ttd_.small_cell_ >> (k - l);
							if(low < cell) return true;
							else return false;
						} else {
							if(l > k) {
								uint64_t low = cell >> (l - k);
								if(ttd_.small_cell_ < low) return true;
								else return false;
							} else {
								if(ttd_.small_cell_ < cell) return true;
								else return false;
							} // if
						} // if
					} // compare()

					// compare two tree_nodes according to postordering
					// true if this < b
					bool compare(const TreeNode& b) const {

						// if either is contained in the other, return the order
						if( is_contained(ttd_.small_cell_, b.small_cell()) ) return true;
						if( is_contained(b.small_cell(), ttd_.small_cell_) ) return false;

						int k = log_2(ttd_.small_cell_);
						int l = log_2(b.small_cell());

						if(k > l) {
							uint64_t low = ttd_.small_cell_ >> (k - l);
							if(low < b.small_cell()) return true;
							else return false;
						} else {
							if(l > k) {
								uint64_t low = b.small_cell() >> (l - k);
								if(ttd_.small_cell_ < low) return true;
								else return false;
							} else {
								if(ttd_.small_cell_ < b.small_cell()) return true;
								else return false;
							} // if
						} // if
					} // compare()


				private:
					bool is_leaf_;
					bool is_root_;
					int level_;				// the level of the node in the tree

					value_type value_;		// this is the app value to be stored at a node
					index_type index_;
					index_type parent_index_;
					std::vector<index_type> children_index_; // list of children indices

					tree_type_data ttd_;	// this is the specific tree structure data
					std::vector<octree_point> point_list_; 
											// list of points belonging to the leaf node

			}; // TreeNode

			typedef TreeNode tree_node;

			class LevelCombine {
				public:
					bool operator()(tree_node& u, const tree_node& v) {

						u.set_level(u.level() + v.level());
						return true;
					} // operator()()
			}; // LevelCombine

			// Types
			typedef LevelCombine level_combine;
			typedef typename std::vector<tree_node>::iterator iterator;
			typedef typename std::vector<tree_node>::const_iterator const_iterator;
			typedef unsigned long size_type;

			// for MPI communications - REDO !!!
			typedef SerialTreeNode<value_type> serial_tree_node;
			typedef ChildData child_data;

		public:

			/**
			 * Constructor -- constructs the tree given list of datapoints from the user
			 */
			BaseTree(const point_type* input_points, const unsigned int& n_local,
					const int& levels, MPI_Comm comm) { 

				MPI_data mpi_data(comm);
				mpi_data_ = mpi_data;

				OctreeData octree_data;

				// create and initialize the octree data from input_points
				octree_data.init<point_type>(input_points, n_local, levels, mpi_data_);
				octree_data.sample_sort(mpi_data_);
				octree_data.load_balance(mpi_data_);
				octree_data.count_unique();

				//if(mpi_data_.rank() == mpi_data_.size() - 1) octree_data.print_data();

				is_accumulation_ready_ = false;
				is_residual_tree_ = false;

				int data_size = octree_data.n();
				int *gather_size = new int[mpi_data_.size()];
				if(gather_size == NULL) {
					std::cerr << "Error allocating memory (in BaseTree())." << std::endl;
					return;
				} // if
				MPI_Allgather(&data_size, 1, MPI_INT, gather_size, 1, MPI_INT, mpi_data_.comm());

				if(mpi_data_.rank() == 0) {
					data_size = 0;
					for(int i = 0; i < mpi_data_.size(); ++ i) data_size += gather_size[i];
					std::cout << std::endl << "*** TREEWORKS v0.9a ***" << std::endl;
					std::cout << "Number of processes: " << mpi_data_.size();
					std::cout << ", Input size: " << data_size << std::endl;
				} // if

				delete[] gather_size;

				// construct the tree using the prepared octree data
				construct_c_octree(octree_data);
			} // BaseTree()


			~BaseTree() {
			} // ~BaseTree()


			/**
			 * Some basic functions
			 */

			size_type size() const { return node_list_.size(); }
			const_iterator begin() const { return node_list_.begin(); }
			const_iterator end() const { return node_list_.end(); }
			iterator begin() { return node_list_.begin(); }
			iterator end() { return node_list_.end(); }

			// function to finalize the tree - always has to be called before finalizing MPI	
			void finalize() {
				MPI_Barrier(mpi_data_.comm());
				// free the one-sided communication window
				MPI_Win_free(&node_list_win_);
				// free the searialized node list memory
				MPI_Free_mem(serialized_node_list_);
			} // finalize()


			// operator[] should be for global tree
			tree_node operator[](index_type n) const {

				tree_node empty_node;
				if(n.index_ == -1) return empty_node;
				if(n.proc_id_ < 0 || n.proc_id_ >= mpi_data_.size()) {
					std::cerr << "Error in index type! proc_id_ out of bounds. Aborting." << std::endl;
					return empty_node;
				} // if

/*				int print_proc = 2;
				if(mpi_data_.rank() == print_proc)
					std::cout << mpi_data_.rank() << " <-- ("
						<< n.proc_id_ << ", " << n.index_ << ")";
*/
				if(n.proc_id_ == mpi_data_.rank()) { // if the node is local to the processor

/*					if(mpi_data_.rank() == print_proc)
						std::cout << std::endl;;
*/
					return node_list_[n.index_];
				} else { // the node needs to be fetched from remote processor

					serial_tree_node* temp_serial_node = new serial_tree_node;
					if(temp_serial_node == 0) {
						std::cout << "Error in allocating memory for temp_serial_node! Aborting."
							<< std::endl;
						return empty_node;
					} // if

					// one-sided communication -- use MPI_Get() to obtain a remote node

/*					if(mpi_data_.rank() == print_proc)
						std::cout << " *****" << std::endl;;
*/
					int r1 = MPI_Win_lock(MPI_LOCK_SHARED, n.proc_id_, 0, node_list_win_);
/*					if(mpi_data_.rank() == print_proc)
						std::cout << " ------ LOCK = " << r1 << " "
							<< n.proc_id_ << " " << node_list_win_ << std::endl;;
*/
					int r2 = MPI_Get(temp_serial_node, sizeof(serial_tree_node), MPI_BYTE, n.proc_id_,
							n.index_, sizeof(serial_tree_node), MPI_BYTE, node_list_win_);
/*					if(mpi_data_.rank() == print_proc)
						std::cout << " ------ GET = " << r2 << " "
							<< n.proc_id_ << " " << n.index_ << " " << node_list_win_ << std::endl;;
*/
					int r3 = MPI_Win_unlock(n.proc_id_, node_list_win_);
/*					if(mpi_data_.rank() == print_proc)
						std::cout << " ------ UNLOCK = " << r3 << " "
							<< n.proc_id_ << " " << node_list_win_ << std::endl;;
*/
					tree_node temp_node;
					unpack_node(temp_serial_node[0], temp_node);

/*					if(mpi_data_.rank() == print_proc) {
						std::cout << temp_node.is_leaf()
							<< " " << temp_node.level()
							<< " " << temp_node.small_cell()
							<< " " << temp_node.large_cell()
							<< " ("	<< temp_node.parent().proc_id_
							<< ", " << temp_node.parent().index_
							<< ") " << temp_node.num_children()
							<< " " << temp_node.num_points()
							<< std::endl;
					} // if
*/
					delete temp_serial_node;

					return temp_node;
				} // if-else
			} // operator[]()


			// these obtain nodes need to be made one-sided communication functions!!!
			// ...
			// given a vector of cell keys (small cell), obtain the corresponding tree_nodes
			bool obtain_nodes(const std::vector<uint64_t>& cells,
					std::vector<tree_node> &out) const {

				// assume that cells does not have duplicates -- removal of duplicates should be
				// done prior to obtain_nodes

				std::cout << mpi_data_.rank() << ". Number of cells to obtain = "
					<< cells.size() << std::endl;

				int *send_count = new int[mpi_data_.size()];
				int *ptr = new int[mpi_data_.size()];
				for(int i = 0; i < mpi_data_.size(); ++ i) {
					send_count[i] = 0;
					ptr[i] = 0;
				} // for

				// test
				if(mpi_data_.rank() == 0) {
					std::cout << mpi_data_.rank() << ". Boundaries = ";
					for(int i = 0; i < mpi_data_.size() * 2; ++ i)
						std::cout << boundaries_[i].small_cell() << " ";
					std::cout << std::endl;
				} // if


				if(mpi_data_.rank() == 0) {
					std::cout << mpi_data_.rank() << ". Request cells (" << cells.size() << ") = ";
					for(int i = 0; i < cells.size(); ++ i)
						std::cout << cells[i] << " ";
					std::cout << std::endl;
				} // if
				// compute request send counts
				std::cout << mpi_data_.rank() << ". Processors = ";
				for(int i = 0; i < cells.size(); ++ i) {
					// use binary searches in boundaries_
					int proc_id = compute_proc(cells[i]);	// replace with a generic binary search !!!!
					if(proc_id != mpi_data_.rank()) ++ send_count[proc_id];

					std::cout << proc_id << " ";
				} // for
				std::cout << std::endl;


				// test
				std::cout << mpi_data_.rank() << " -- ";
				for(int i = 0; i < mpi_data_.size(); ++ i)
					std::cout << send_count[i] << " ";
				std::cout << std::endl;


				int send_total = 0;
				for(int i = 0; i < mpi_data_.size(); ++ i) {
				   ptr[i] = send_total;
				   send_total += send_count[i];
				} // for

				std::cout << mpi_data_.rank() << ". A" << std::endl;

				// prepare the send buffers to send the requests
				uint64_t *send_request = new uint64_t[send_total];
				for(int i = 0; i < cells.size(); ++ i) {
					// use binary searches in boundaries_
					int proc_id = compute_proc(cells[i]);	// replace with a generic binary search !!!!
					if(proc_id == mpi_data_.rank()) continue;
					send_request[ptr[proc_id]] = cells[i];
					++ ptr[proc_id];
				} // for

				std::cout << mpi_data_.rank() << ". B" << std::endl;

				// send and receive all requests
				uint64_t *recv_request = 0;
				int *recv_count = 0;
				int recv_size = 0;
				twlib::all_to_all(send_count, send_request, recv_count, recv_request, recv_size,
						mpi_data_.comm());

				std::cout << mpi_data_.rank() << ". C" << std::endl;

				// now that the requests have been received, prepare the buffers to send in response
				// to the requests.

				serial_tree_node *send_buff = new serial_tree_node[recv_size];
				tree_node empty_node;

				for(int i = 0; i < recv_size; ++ i) {
					int index = node_search(node_list_, recv_request[i]);	// replace with a generic
																			// binary search !!!!
					if(index == -1) { // if no such node exists
						// pack the empty node (small_cell and large_cell are 0)
						pack_node(send_buff[i], empty_node);
					} else {
						pack_node(send_buff[i], node_list_[index]);
					} // if-else
				} // for

				std::cout << mpi_data_.rank() << ". D" << std::endl;

				// perform the communication to send the requested nodes
				serial_tree_node *recv_buff = 0;
				int *temp_recv_count = 0;
				twlib::all_to_all(recv_count, send_buff, temp_recv_count, recv_buff, recv_size,
						mpi_data_.comm());

				// unpack and sort the received nodes
				std::vector<tree_node> recv_nodes;
				tree_node temp_node;
				for(int i = 0; i < recv_size; ++ i) {
					unpack_node(recv_buff[i], temp_node);
					if(temp_node.small_cell() != 0) recv_nodes.push_back(temp_node);
				} // for
				postorder_sort(0, recv_nodes.size() - 1, recv_nodes);

				std::cout << mpi_data_.rank() << ". E" << std::endl;

				// prepare the out list of the nodes
				for(int i = 0; i < cells.size(); ++ i) {
					int proc_id = compute_proc(cells[i]);
					if(proc_id == mpi_data_.rank()) { // if the node is local
						int index = node_search(node_list_, cells[i]);
						if(index != -1) out.push_back(node_list_[index]);
					} else {
						int index = node_search(recv_nodes, cells[i]);
						if(index != -1 && recv_nodes[index].small_cell() != 0)
							out.push_back(recv_nodes[index]);
					} // if-else
				} // for

				std::cout << mpi_data_.rank() << ". F" << std::endl;

				// remove duplicates
				postorder_sort(0, out.size() - 1, out);
				typename std::vector<tree_node>::iterator iter1 = out.begin();
				typename std::vector<tree_node>::iterator iter2 = out.begin();
				++ iter2;
				while(iter2 != out.end()) {
					if(*iter1 == *iter2) {
						out.erase(iter2);
						iter2 = iter1;
						++ iter2;
					} else {
						++ iter1;
						++ iter2;
					} // if-else
				} // while


				// delete all allocated memories
				delete[] recv_buff;
				delete[] temp_recv_count;
				delete[] send_buff;
				delete[] recv_request;
				delete[] recv_count;
				delete[] send_request;
				delete[] ptr;
				delete[] send_count;

				// NOTE:
				// Since the tree is stored in postorder, spatial locality is maintain to a large extent
				// and obtaining nodes will create multiple copies of many nodes which are local
				// to the processor. Try to avoid this by making a list of only those tree_nodes
				// which are obtained from other processors. But how ???

				std::cout << mpi_data_.rank() << ". Obtained nodes." << std::endl;

				return true;
			} // obtain_nodes()


			// given vector of index_type indices of the nodes, obtain the tree_nodes
			bool obtain_nodes(const std::back_insert_iterator<std::vector<index_type> >& indices,
					std::back_insert_iterator<std::vector<tree_node> > &out) {
				// ...
				return true;
			} // obtain_nodes()


			int compute_proc(const uint64_t &cell) const {
				// perform a binary search of cell on the boundaries_ array
				// and return the processor id (position = i, then proc_id = i/2, and
				// i HAS to be an odd number). otherwise return -1.
				
				for(int i = 0; i < boundaries_.size() ; ++ i) {
					if(boundaries_[i] < cell) continue;
					if(boundaries_[i] == cell || boundaries_[i] > cell) return i/2;
				} // for

				// FIX IT ...

/*				int left = 0, right = boundaries_.size() - 1;
				int mid = (left + right) >> 1;
				int proc_id = -1;

				while(left <= right) {
					// base cases
					if(boundaries_[left] == cell) return (left >> 1);
					if(boundaries_[right] == cell) return (right >> 1);
					if(right == left + 1 && boundaries_[left] < cell && boundaries_[right] > cell)
						return (left >> 1);

					// general
					if(boundaries_[mid] < cell) {
						if(mid % 2 == 0) left = mid;
						else left = mid + 1;
					} else {
						if(mid % 2 == 0) right = mid + 1;
						else right = mid;
					} // if-else

					mid = (left + right) >> 1;
				} // while
*/
				return -1;
			} // compute_proc()


			// this binary search should become generic and goto twlib file
			int node_search(const std::vector<tree_node>& node_list, const uint64_t& cell) const {

				int left = 0, right = node_list.size() - 1;
				int mid = (left + right) >> 1;

				while(left <= right) {
					if(node_list[left].is_leaf() && node_list[left] == cell) return left;
					if(node_list[right].is_leaf() && node_list[right] == cell) return right;
					if(node_list[mid].is_leaf() &&
							(is_contained(cell, node_list[mid].large_cell())
							 || node_list[mid].small_cell() == cell))
						return mid;

					if(node_list[mid] > cell)
						right = mid - 1;
					else
						left = mid + 1;
					mid = (left + right) >> 1;
				} // while

				return -1;
			} // node_search()


			/**
			 * Octree specific functions
			 */

			// compute the same sized cell keys lying in the nearfield of a given cell
			// (within the domain, and adjacent cells including self, hence, max of 27)
			bool compute_adjacent_cells(const uint64_t& cell,
					std::vector<uint64_t>& adjacent_cells) const {

				data_point center;
				double sidelength = 0.0;
				compute_center_and_sidelength(cell, center, sidelength);

				uint64_t temp_cell = 0;

				if(compute_cell(center.x_ + sidelength, center.y_, center.z_, sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_ + sidelength, center.y_, center.z_ + sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_ + sidelength, center.y_, center.z_ - sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);

				if(compute_cell(center.x_ + sidelength, center.y_ + sidelength, center.z_,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_ + sidelength, center.y_ + sidelength, center.z_ + sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_ + sidelength, center.y_ + sidelength, center.z_ - sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);

				if(compute_cell(center.x_ + sidelength, center.y_ - sidelength, center.z_,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_ + sidelength, center.y_ - sidelength, center.z_ + sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_ + sidelength, center.y_ - sidelength, center.z_ - sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);

				if(compute_cell(center.x_, center.y_, center.z_, sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_, center.y_, center.z_ + sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_, center.y_, center.z_ - sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);

				if(compute_cell(center.x_, center.y_ + sidelength, center.z_,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_, center.y_ + sidelength, center.z_ + sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_, center.y_ + sidelength, center.z_ - sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);

				if(compute_cell(center.x_, center.y_ - sidelength, center.z_,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_, center.y_ - sidelength, center.z_ + sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_, center.y_ - sidelength, center.z_ - sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);

				if(compute_cell(center.x_ - sidelength, center.y_, center.z_, sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_ - sidelength, center.y_, center.z_ + sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_ - sidelength, center.y_, center.z_ - sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);

				if(compute_cell(center.x_ - sidelength, center.y_ + sidelength, center.z_,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_ - sidelength, center.y_ + sidelength, center.z_ + sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_ - sidelength, center.y_ + sidelength, center.z_ - sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);

				if(compute_cell(center.x_ - sidelength, center.y_ - sidelength, center.z_,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_ - sidelength, center.y_ - sidelength, center.z_ + sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);
				if(compute_cell(center.x_ - sidelength, center.y_ - sidelength, center.z_ - sidelength,
							sidelength, temp_cell))
					adjacent_cells.push_back(temp_cell);

				return true;
			} // compute_adjacent_cells

			// given a cell key, compute the center and sidelength
			bool compute_center_and_sidelength(const uint64_t& in_cell,
					data_point& center, double& sidelength) const {

				int cell_bits = log_2(in_cell);

				uint64_t cell = in_cell << (64 - cell_bits + 1);
				int num_bits = (cell_bits - 1) / 3;

				unsigned int x = 0, y = 0, z = 0;
				uint64_t mask = 1;
				mask = mask << 63;
				
				for(int i = 0; i < num_bits; ++ i) {
					x = x << 1;
					if(cell & mask)	x = x | 1;
					cell = cell << 1;

					y = y << 1;
					if(cell & mask) y = y | 1;
					cell = cell << 1;

					z = z << 1;
					if(cell & mask) z = z | 1;
					cell = cell << 1;
				} // for

				mask = 1 << num_bits;
				sidelength = domain_ / mask;

				center.x_ = (x + 0.5) * sidelength;
				center.y_ = (y + 0.5) * sidelength;
				center.z_ = (z + 0.5) * sidelength;

				return true;
			} // compute_center_and_sidelength()

			// given center and sidelength, compute the key of the cell
			bool compute_cell(const double& center_x, const double& center_y, const double& center_z,
					const double& sidelength, uint64_t& cell) const {

				unsigned int x = 0, y = 0, z = 0;

				if(center_x > 0 && center_x < domain_ &&
						center_y > 0 && center_y < domain_ &&
						center_z > 0 && center_z < domain_) {

					x = (unsigned int) (center_x / sidelength);
					y = (unsigned int) (center_y / sidelength);
					z = (unsigned int) (center_z / sidelength);

					unsigned int height = log_2((uint64_t) (domain_ / sidelength)) - 1;

					cell = compute_key(x, y, z, height);

					return true;
				} // if

				return false;
			} // compute_cell()

			// given the x, y and z coordinates of the cell and the height, compute the cell key
			uint64_t compute_key(unsigned int x, unsigned int y, unsigned int z,
					unsigned int height) const {

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


			// given a cell, compute its immediate parent cell (remove last 3 bits)
			uint64_t compute_parent_cell(const uint64_t cell) {	return (cell >> 3);	}

			// compute the size of the small cell of the node
			double cell_size(const tree_node& node) const {

				int cell_bits = log_2(node.small_cell());
				int num_bits = (cell_bits - 1) / 3;

				uint64_t mask = 1 << num_bits;
				double sidelength = domain_ / mask;

				return sidelength;
			} // cell_size()


			// compute the cell size of the smallest cell (small cell is leaf)
			double min_cell_size() const {

				int i = 0;
				while((!node_list_[i].is_leaf()) && i < node_list_.size()) ++ i;
				if(i == node_list_.size()) i = 0;

				return cell_size(node_list_[i]);
			} // min_cell_size()


			// compute the coordinates of the small cell of node
			bool compute_node_coordinates(const TreeNode& node,
					double& x1, double& x2, double& y1, double& y2, double& z1, double& z2) const {

				compute_coordinates(node.small_cell(), domain_, x1, x2, y1, y2, z1, z2);

				return true;
			} // compute_coordinates()


			// check if the two nodes are adjacent
			bool adjacent(const TreeNode& a, const TreeNode& b) const {

				double xa1 = 0.0, xa2 = 0.0, ya1 = 0.0, ya2 = 0.0, za1 = 0.0, za2 = 0.0;
				double xb1 = 0.0, xb2 = 0.0, yb1 = 0.0, yb2 = 0.0, zb1 = 0.0, zb2 = 0.0;

				compute_node_coordinates(a, xa1, xa2, ya1, ya2, za1, za2);
				compute_node_coordinates(b, xb1, xb2, yb1, yb2, zb1, zb2);

				if(is_equal(xa1, xb1) || is_equal(xa1, xb2) || is_equal(xa2, xb1) ||
						is_equal(xa2, xb2) || is_equal(ya1, yb1) || is_equal(ya1, yb2) ||
						is_equal(ya2, yb1) || is_equal(ya2, yb2) ||	is_equal(za1, zb1) ||
						is_equal(za1, zb2) || is_equal(za2, zb1) || is_equal(za2, zb2))
					return true;

				return false;
			} // adjacent()


			// compute the largest distances between a and b (within small cells)
			double largest_distance(const TreeNode& a, const TreeNode& b) const {

				double ax1 = 0.0, ax2 = 0.0, ay1 = 0.0, ay2 = 0.0, az1 = 0.0, az2 = 0.0;
				double bx1 = 0.0, bx2 = 0.0, by1 = 0.0, by2 = 0.0, bz1 = 0.0, bz2 = 0.0;

				// ??1 is smaller than ??2
				compute_coordinates(a.small_cell(), domain_, ax1, ax2, ay1, ay2, az1, az2);
				compute_coordinates(b.small_cell(), domain_, bx1, bx2, by1, by2, bz1, bz2);

				double a_x = 0.0, b_x = 0.0;
				if(fabs(ax1 - bx2) > fabs(ax2 - bx1)) {
					a_x = ax1; b_x = bx2;
				} else {
					a_x = ax2; b_x = bx1;
				} // if-else
				double a_y = 0.0, b_y = 0.0;
				if(fabs(ay1 - by2) > fabs(ay2 - by1)) {
					a_y = ay1; b_y = by2;
				} else {
					a_y = ay2; b_y = by1;
				} // if-else
				double a_z = 0.0, b_z = 0.0;
				if(fabs(az1 - bz2) > fabs(az2 - bz1)) {
					a_z = az1; b_z = bz2;
				} else {
					a_z = az2; b_z = bz1;
				} // if-else

				return distance(a_x, a_y, a_z, b_x, b_y, b_z);

			} // largest_distnace()

			// compute the smallest distances between a and b (within small cells)
			double smallest_distance(const TreeNode& a, const TreeNode& b) const {

				if(is_contained(a.small_cell(), b.small_cell())
						|| is_contained(b.small_cell(), a.small_cell())
						|| a.small_cell() == b.small_cell())
					return 0.0;

				double ax1 = 0.0, ax2 = 0.0, ay1 = 0.0, ay2 = 0.0, az1 = 0.0, az2 = 0.0;
				double bx1 = 0.0, bx2 = 0.0, by1 = 0.0, by2 = 0.0, bz1 = 0.0, bz2 = 0.0;

				// ??1 is smaller than ??2
				compute_coordinates(a.small_cell(), domain_, ax1, ax2, ay1, ay2, az1, az2);
				compute_coordinates(b.small_cell(), domain_, bx1, bx2, by1, by2, bz1, bz2);

				double a_x = 0.0, b_x = 0.0;
				if(fabs(ax1 - bx2) < fabs(ax2 - bx1)) {
					a_x = ax1; b_x = bx2;
				} else {
					a_x = ax2; b_x = bx1;
				} // if-else
				double a_y = 0.0, b_y = 0.0;
				if(fabs(ay1 - by2) < fabs(ay2 - by1)) {
					a_y = ay1; b_y = by2;
				} else {
					a_y = ay2; b_y = by1;
				} // if-else
				double a_z = 0.0, b_z = 0.0;
				if(fabs(az1 - bz2) < fabs(az2 - bz1)) {
					a_z = az1; b_z = bz2;
				} else {
					a_z = az2; b_z = bz1;
				} // if-else

				return distance(a_x, a_y, a_z, b_x, b_y, b_z);

			} // smallest_distnace()


			/**
			 * Accumulation functions
			 */

			template<typename CombineFunction>
			bool upward_accumulate(CombineFunction& combine, const MPI_data& mpi_data) {

				if(!is_residual_tree_) construct_residual_tree(mpi_data);

				// improve the following to see if updates
				// to the contraction information can be avoided
				// ...

				start_end_.resize(node_list_.size());
				node_flag_.resize(node_list_.size());

				std::vector<int> upward_array(node_list_.size());

				// initialize the above vectors
				for(int i = 0; i < node_list_.size(); ++ i) {
					start_end_[i] = 'N';
					node_flag_[i] = 'A';
					upward_array[i] = node_list_[i].num_children();
				}

				// local accumulation, and mark the finished nodes as 'D'
				for(int i = 0; i < node_list_.size(); ++ i) {
					if(upward_array[i] == 0 && node_list_[i].parent().proc_id_ == mpi_data.rank()) {
						if(node_list_[i].parent().index_ != -1) {
							// combine is applied to the parent
							combine(node_list_[node_list_[i].parent().index_], node_list_[i]);
							-- upward_array[node_list_[i].parent().index_];
						}
						node_flag_[i] = 'D';
					} //else if(i > node_list_.size() - new_node_count_) break;
				} // for

				int iter = 0; 	// the number of iteration, used in parallel accumulation
				int* recv_upch = new int[mpi_data.size()];	// same as rcv_upch in Bhaanu's code
				for(int i = 0; i < mpi_data.size(); ++ i) recv_upch[i] = 0;

				// perform accumulation in parallel
				while(1) {
					int i = 0;
					while(i < mpi_data.size() - 1 && recv_upch[i] < 0) ++ i;
					if(i == mpi_data.size() - 1 && recv_upch[i] == 0) break;
					else {
						++ iter;
						MPI_Barrier(mpi_data.comm());
						parallel_upward_accumulation(iter, upward_array,
								recv_upch, start_end_, node_flag_, combine, mpi_data);
					} // if-else
				} // while

				delete[] recv_upch;

				return true;
			} // upward_accumulate


			template<typename CombineFunction>
			bool downward_accumulate(CombineFunction& combine, const MPI_data& mpi_data) {

				if(!is_residual_tree_) construct_residual_tree(mpi_data);
				if(!is_accumulation_ready_) construct_contraction(mpi_data);

				// make a copy of the initial node and use it for the default value of value
				tree_node default_node;
				int i = 0;
				while((node_list_[i].is_leaf() || node_list_[i].is_root()) &&
						i < node_list_.size())
					++ i;
				default_node.set_value(node_list_[i].value());

				// prepare upward_array
				std::vector<int> upward_array(node_list_.size());
				for(i = 0; i < node_list_.size(); ++ i)
					upward_array[i] = node_list_[i].num_children();

				for(i = 0; i < residual_array_.size(); ++ i) node_flag_[residual_array_[i]] = 'R';

				for(i = accumulation_iter_; i > 0; -- i)
					parallel_downward_accumulation(i, default_node, upward_array, combine, mpi_data);

				for(i = node_list_.size() - 1; i >= 0; -- i) {
					if(node_flag_[i] != 'R' && node_list_[i].parent().index_ != -1)
						combine(node_list_[i], node_list_[node_list_[i].parent().index_]);
				} // for

				return true;
			} // downward_accumulate


			/**
			 * Local Computations
			 */

			template<typename CombineFunction>
			bool local_compute(CombineFunction& combine, const MPI_data& mpi_data) {
				// for each node in the node_list_, apply combine on it
				for(int i = 0; i < node_list_.size(); ++ i) {
					combine(node_list_[i], node_list_[i]);
				} // for

				return true;
			} // local_compute()


			/**
			 * No dependency computations
			 */

			template<typename CombineFunction>
			bool no_dependency_compute(CombineFunction& combine,
					const std::vector<std::vector<tree_node> >& interaction_sets,
					const MPI_data& mpi_data) {
				// for each node in node_list_, apply combine to it for all nodes in its i-set
				for(int i = 0; i < node_list_.size(); ++ i) {
					for(int j = 0; j < interaction_sets[i].size(); ++ j) {
						combine(node_list_[i], interaction_sets[i][j]);
					} // for
				} // for

				return true;
			} // no_dependency_compute()


			void print_tree() {
				typedef typename tree_node::const_children_iterator children_iter;

				for(unsigned int i = 0; i < node_list_.size(); ++ i) {
					std::cout << i << ". ";
					node_list_[i].value().print();
					std::cout << " " << node_list_[i].is_leaf()
						<< " " << node_list_[i].level()
						<< " " << node_list_[i].small_cell()
						<< " " << node_list_[i].large_cell()
						<< " ("	<< node_list_[i].parent().proc_id_
						<< ", " << node_list_[i].parent().index_
						<< ") " << node_list_[i].num_children()
						<< " " << node_list_[i].num_points()
						<< " [";
					std::pair<children_iter, children_iter> children = node_list_[i].children();
					for(children_iter iter = children.first; iter != children.second;
							++ iter)
						std::cout << "(" << (*iter).proc_id_ << ", " << (*iter).index_ << ") ";
					std::cout << "]" << std::endl;
				} // for
			} // print_tree()

		private:

			/**
			 * Tree Construction
			 * octree_data is list input octree data points, sorted according to their cells
			 */
			bool construct_c_octree(const OctreeData& octree_data) {

				if(mpi_data_.rank() == 0) std::cout << "+ constructing compressed octree ... ";

				typedef typename tree_node::const_children_iterator children_iter;

				// create tree boundary array: all gather of first cell of each processor
				typedef OctreeData::OctreePoint octree_point;
				octree_point *boundary_array;
				octree_point first_cell = octree_data[0];
				twlib::gather_cells(first_cell, boundary_array, mpi_data_);

				MPI_Barrier(mpi_data_.comm());
				double construct_time_s = MPI_Wtime();

				// create temporary nodes
				std::vector<tree_node> temp_nodes(2 * octree_data.n_unique());
				if(temp_nodes.size() != 2 * octree_data.n_unique()) {
					std::cerr << "Error with memory allocation (temp_nodes in construct_c_octree())."
						<< std::endl;
					return false;
				} // if

				// store all leaf nodes in temp_nodes
				unsigned int j = octree_data.n_unique() - 1;
				for(unsigned int i = octree_data.n() - 1; i > 0; --i) {
					if(octree_data[i] != octree_data[i-1]) {
						temp_nodes[j].set_small_cell(octree_data[i].cell_key_);
						--j;
					} // if
				} // for
				temp_nodes[0].set_small_cell(octree_data[0].cell_key_);

				// insert the first cell
				temp_nodes[0].set_large_cell(1);	// set its large_cell as the whole domain
				index_type init_index = make_index(mpi_data_.rank(), -1);
				temp_nodes[0].set_parent(init_index);
				temp_nodes[0].make_leaf();
				
				j = octree_data.n_unique();

				// insert all leaves
				for(unsigned int i = 1; i < octree_data.n_unique(); ++i) {
					insert_leaf(temp_nodes[i], temp_nodes, i, j, octree_data.n_unique());
					temp_nodes[i].make_leaf();
				} // for

				// insert the leaf borrowed from the next processor
				if(mpi_data_.rank() != mpi_data_.size() - 1) {
					tree_node borrowed_leaf;
					borrowed_leaf.set_small_cell(boundary_array[mpi_data_.rank() + 1].cell_key_);
					borrowed_leaf.make_leaf();
					insert_leaf(borrowed_leaf, temp_nodes, octree_data.n_unique(), j,
							octree_data.n_unique());
				} // if

				// construct the postordering of the local tree
				int total_nodes = j;
				int count = 0;		// number of nodes
				int local_root = get_root(temp_nodes);
				construct_postorder(temp_nodes, local_root, count);

				// put the datapoints into the leaves
				insert_data_points(octree_data, total_nodes);

				MPI_Barrier(mpi_data_.comm());

				// identify 'out of order' nodes and put them on the correct processor
				// (last_small_cell is the small cell of the last leaf)
				uint64_t last_small_cell =
						temp_nodes[octree_data.n_unique() - 1].small_cell();
				globalize_tree(last_small_cell, octree_data.n_unique(), total_nodes,
						octree_data.levels(), boundary_array);

				// set the last node on last processor as root
				if(mpi_data_.rank() == mpi_data_.size() - 1)
					node_list_[node_list_.size() - 1].make_root();

				// set the index for each node
				for(int i = 0; i < node_list_.size(); ++ i)
					node_list_[i].set_index(make_index(mpi_data_.rank(), i));

				levels_ = octree_data.levels();
				// number of leaves = octree_data.n_unique() (if needed in future)
				domain_ = octree_data.domain();

				delete[] boundary_array;

				MPI_Barrier(mpi_data_.comm());
				double construct_time_e = MPI_Wtime();

				int num_nodes = node_list_.size();
				int* total_num_nodes = new int[mpi_data_.size()];
				if(total_num_nodes == NULL) {
					std::cerr << "Error in memory allocation (total_num_nodes in construct_c_octree())."
						<< std::endl;
					return false;
				} // if
				MPI_Allgather(&num_nodes, 1, MPI_INT, total_num_nodes, 1, MPI_INT, mpi_data_.comm());
				if(mpi_data_.rank() == 0) {
					std::cout << "done: "
						<< (construct_time_e - construct_time_s) * 1000 << "ms." << std::endl;
					num_nodes = 0;
					for(int i = 0; i < mpi_data_.size(); ++i) num_nodes += total_num_nodes[i];
					std::cout << "Number of nodes: " << num_nodes;
					std::cout << ", Levels: " << levels_ << std::endl;
				} // if

				delete[] total_num_nodes;

				// Tree construction is complete,
				// now perform post processing to make the tree ready

				double postprocess_time_s = MPI_Wtime();

				// make the tree ready for accumulation computations:
				// (probably merge residual tree construction and contraction together ... )
				construct_residual_tree(mpi_data_);
				construct_contraction(mpi_data_);
				compute_levels(mpi_data_);

				// Collect the leftmost and rightmost nodes cells from all processors.
				// This is useful in identifying the processor id where any given cell should reside.
				serial_tree_node *leftmost, *rightmost; // the local leftmost and rightmost nodes
				leftmost = new serial_tree_node;
				rightmost = new serial_tree_node;
				serial_tree_node *all_leftmost, *all_rightmost; // arrays for the results
				int recv_size = 0;

				pack_node(leftmost[0], node_list_[0]);
				pack_node(rightmost[0], node_list_[node_list_.size() - 1]);

				twlib::all_gather(leftmost, 1, all_leftmost, recv_size, mpi_data_.comm());
				twlib::all_gather(rightmost, 1, all_rightmost, recv_size, mpi_data_.comm());

				tree_node temp_node;
				for(int i = 0; i < mpi_data_.size(); ++ i) {
					unpack_node(all_leftmost[i], temp_node);
					boundaries_.push_back(temp_node);
					unpack_node(all_rightmost[i], temp_node);
					boundaries_.push_back(temp_node);
				} // for

				delete[] all_rightmost;
				delete[] all_leftmost;
				delete rightmost;
				delete leftmost;

				// Make the tree ready for One-sided communications
				//
				// Serialize the node_list_ and store it in a linear array
				// Each node in this linear array will be of exactly the same size
				// (use the serial_tree_node). This is not the general case since this
				// will be specialized only for octrees. Think of a better way to do this
				// later ...

				/*unsigned long int memory_a = jaz::mem_usage();
				if(mpi_data_.rank() == 0)
					std::cout << "Memory usage: " << (memory_a / 1024) << " KB" << std::endl; */

				int size = node_list_.size() * sizeof(serial_tree_node);
				int ret = MPI_Alloc_mem(size, MPI_INFO_NULL, &serialized_node_list_);
				if(ret == MPI_ERR_NO_MEM) {
					std::cerr << "Error in allocating serial MPI memory! Aborting." << std::endl;
					return false;
				} // if

				// serialize the tree
				serialize_tree(node_list_, serialized_node_list_);

				// Create a window on this linear array
				MPI_Win_create(serialized_node_list_, size, sizeof(serial_tree_node),
						            MPI_INFO_NULL, mpi_data_.comm(), &node_list_win_);

				MPI_Barrier(mpi_data_.comm());
				double postprocess_time_e = MPI_Wtime();

				/*std::cout << mpi_data_.rank() << ". Last Node = "
					<< serialized_node_list_[node_list_.size() - 1].is_leaf_
					<< " " << serialized_node_list_[node_list_.size() - 1].level_
					<< " " << serialized_node_list_[node_list_.size() - 1].ttd_.small_cell_
					<< " " << serialized_node_list_[node_list_.size() - 1].ttd_.large_cell_
					<< " ("	<< serialized_node_list_[node_list_.size() - 1].parent_.proc_id_
					<< ", " << serialized_node_list_[node_list_.size() - 1].parent_.index_
					<< ") " << serialized_node_list_[node_list_.size() - 1].num_children_
					<< " " << serialized_node_list_[node_list_.size() - 1].num_points_
					<< std::endl; */

				if(mpi_data_.rank() == 0) std::cout << "postprocessing done: "
					<< (postprocess_time_e - postprocess_time_s) * 1000 << "ms." << std::endl;

				unsigned long int memory = jaz::mem_usage();
				if(mpi_data_.rank() == 0)
					std::cout << "Memory usage: " << (memory / 1024) << " KB" << std::endl;

				//std::cout << mpi_data_.rank() << ". win = " << node_list_win_
				//	<< ", addres = " << &node_list_win_ << std::endl;

				/*if(mpi_data_.rank() == mpi_data_.size() - 1) {
					std::cout << "========== " << mpi_data_.rank() << " ==========" << std::endl;
					std::cout << "Points = " << octree_data.n() 
						<< ", Unique cells = " << octree_data.n_unique()
						<< ", Levels = " << levels_ << std::endl;
					std::cout << "Total nodes = " << total_nodes << " (size = " << node_list_.size()
						<< ")" << std::endl;
					int point_count = 0;
					for(unsigned int i = 0; i < node_list_.size(); ++ i) {
						std::cout << i << ". " << node_list_[i].is_leaf()
							<< " " << node_list_[i].level()
							<< " " << node_list_[i].small_cell()
							<< " " << node_list_[i].large_cell()
							<< " ("	<< node_list_[i].parent().proc_id_
							<< ", " << node_list_[i].parent().index_
							<< ") " << node_list_[i].num_children()
							<< " " << node_list_[i].num_points()
							<< " [";
						std::pair<children_iter, children_iter> children = node_list_[i].children();
						int k = 0;
						for(children_iter iter = children.first; iter != children.second;
								++ iter, ++ k)
							std::cout << "(" << (*iter).proc_id_ << ", " << (*iter).index_ << ") ";
						std::cout << "]";
						if(k > MAX_CHILDREN) std::cout << " <--------- ERROR!";
						std::cout << std::endl;
						point_count += node_list_[i].num_points();
					} // for
					std::cout << "Total number of points = " << point_count << std::endl;
				} // if */

				return true;
			} // construct_c_octree()


			// to serialize the node_list_ into array of serial_tree_node
			// for one-sided communication facilitation
			bool serialize_tree(const std::vector<tree_node>& node_list,
					serial_tree_node* &serialized_node_list) {

				if(mpi_data_.rank() == 0) std::cout << "+ serializing tree ... ";
				for(int i = 0; i < node_list.size(); ++ i)
					pack_node(serialized_node_list[i], node_list[i]);

				if(mpi_data_.rank() == 0) std::cout << "done." << std::endl;

				return true;
			} // serialize_tree


			bool compute_levels(const MPI_data& mpi_data) {
				// compute levels of the nodes of the tree (NOT cells) using accumulations
				if(mpi_data.rank() == 0) std::cout << "+ computing levels of nodes ... ";

				level_combine level_combine_function;
				downward_accumulate(level_combine_function, mpi_data);

				if(mpi_data.rank() == 0) std::cout << "done." << std::endl;

				return true;
			} // compute_levels()


			// merge the local trees from all procs into one global tree
			bool globalize_tree(const uint64_t last_small_cell, unsigned int m,
					int& total_nodes, unsigned int levels,
					const OctreeData::OctreePoint* boundary_array) {

				typedef typename tree_node::const_children_iterator children_iter;

				// count the number of 'out of order' (oo) nodes (nodes to right of borrowed leaf)
				int num_oo_nodes = 0;
				if(mpi_data_.rank() != mpi_data_.size() - 1) {
					int i = total_nodes - 1;
					while(node_list_[i].small_cell() != last_small_cell) {
						if(is_contained(boundary_array[mpi_data_.rank() + 1].cell_key_,
									node_list_[i].small_cell())) {
							++ num_oo_nodes;
							-- i;	
						} else {
							break;
						} // if-else
					} // while
				} // if

				// prepare the send buffers containing the oo nodes to be sent,
				// identify the processor id to where to send the oo nodes, and
				// perform the communications to send and receive the oo nodes

				int* send_count = new int[mpi_data_.size()];
				int* ptr = new int[mpi_data_.size()];

				// communication buffer cannot be of type tree_node due to presence of vectors!!!!
				serial_tree_node* send_buff = new serial_tree_node[num_oo_nodes];

				if(send_count == NULL || ptr == NULL || send_buff == NULL) {
					std::cerr << "Error in memory allocation (globalize_tree())." << std::endl;
					return false;
				} // if

				int* recv_count;
				serial_tree_node* recv_buff;

				// initialize counts
				for(int i = 0; i < mpi_data_.size(); ++i) {
					send_count[i] = 0;
					ptr[i] = 0;
				} // for
				// compute send counts for each processor
				if(num_oo_nodes != 0) {
					for(int i = total_nodes - 1; i >= total_nodes - num_oo_nodes; -- i) {
						int j = mpi_data_.size() - 1;
						while( !is_contained(boundary_array[j].cell_key_, 
									node_list_[i].small_cell()) ) {
							-- j;
						} // while
						++ send_count[j];
					} // for
				} // if
				int total_send_count = send_count[0];
				for(int i = 1; i < mpi_data_.size(); ++i) {
					ptr[i] = total_send_count;
					total_send_count += send_count[i];
				} // for

				int* proc_array = new int[num_oo_nodes];
					// to store the proc id where oo nodes are sent
				if(proc_array == NULL) {
					std::cerr << "Error in memory allocation (globalize_tree())." << std::endl;
					return false;
				} // if

				// construct the send buffer
				if(num_oo_nodes != 0) {
					for(int i = total_nodes - 1; i >= total_nodes - num_oo_nodes; -- i) {
						// find which processor this node should go to
						int j = mpi_data_.size() - 1;
						while( !is_contained(boundary_array[j].cell_key_, 
									node_list_[i].small_cell()) )
							-- j;

						// Count the number of children of this OO node that are also OO
						// The structure member num_children_ of a node indicates the number
						// of children local to the processor. That is why count must be subtracted
						// from node_list_[i].num_children_ since those OO child nodes will
						// no longer be local to the processor. This decrement is accounted
						// for at a later stage by an equal number of increments (when the OO
						// children are looking for their parent nodes.)
						int count = 0;
						int k = 0;
						std::pair<children_iter, children_iter> children = node_list_[i].children();
						for(children_iter iter = children.first; iter != children.second;
								++ iter, ++ k) {
							if((*iter).index_ >= total_nodes - num_oo_nodes)
								++ count;
							else // set the new parent processor location
								node_list_[(*iter).index_].set_parent_proc(j);
						} // for
						node_list_[i].remove_last_children(count);
						pack_node(send_buff[ptr[j]], node_list_[i]);
						node_list_.pop_back();

						// proc_array stores the proc id where OO is sent
						proc_array[ptr[j]] = j;
						++ ptr[j];
					} // for
				} // if

				// All processors now exchange the respective OO nodes.
				int recv_size = 0;
				twlib::all_to_all(send_count, send_buff, recv_count, recv_buff,
						recv_size, mpi_data_.comm());	// recv_size is total count

				// make a note of all local nodes whose parent index needs to be updated later
				std::vector<std::pair<int, int> > update_indices;
				for(int i = 0; i < num_oo_nodes; ++ i) {
					for(int j = 0; j < send_buff[i].num_children_; ++ j) {
						int index = send_buff[i].children_[j].index_;
						if(index < total_nodes - num_oo_nodes)
							update_indices.push_back(std::make_pair(index, i));
					} // for
				} // for

				// Merge the received oo nodes to local nodes and
				// let the other 'sending' procs know about the new location

				new_node_count_ = 0;
				int* notify_array = new int[recv_size];
				bool* present_array = new bool[recv_size];
				if(notify_array == NULL || present_array == NULL) {
					std::cerr << "Error in memory allocation (globalize_tree())." << std::endl;
					return false;
				} // if
				std::vector<tree_node> temp_recv_buff;

				// the node list is scanned as many times as recv_size - could be improved ???
				// (may be use binary search first in the local nodes, and then scan the temp_recv_buff)
				// check if a received node is already present in local node list
				// if it is, merge its children array, otherwise, store the new node in temp_recv_buff
				for(int i = 0; i < recv_size; ++ i) {
					present_array[i] = false;
					tree_node temp_node;
					unpack_node(recv_buff[i], temp_node);

					for(int j = 0; j < total_nodes - num_oo_nodes + new_node_count_; ++ j) {
						if(node_list_[j].small_cell() == temp_node.small_cell()) {
							// for every child in temp_node, add it to node_list_[j]
							std::pair<children_iter, children_iter> children = temp_node.children();
							for(children_iter iter = children.first; iter != children.second; ++ iter)
								node_list_[j].add_child(*iter);
							present_array[i] = true;

							// need to check and add the child to node in temp_recv_buff
							for(int k = 0; k < new_node_count_; ++ k) {
								if(temp_recv_buff[k].small_cell() == temp_node.small_cell()) {
									std::pair<children_iter, children_iter> children =
											temp_node.children();
									for(children_iter iter = children.first;
											iter != children.second; ++ iter)
										temp_recv_buff[k].add_child(*iter);
								} // if
							} // for

							break;
						} // if-else
					} // while

					// if it is not present, it is unpacked into temp_recv_buff
					if(!present_array[i]) {
						temp_recv_buff.push_back(temp_node);

						// also temporarily insert it at the end of node_list_ in order to
						// increase node_list_ and take care of duplicates
						node_list_.push_back(temp_node);

						++ new_node_count_;
					} // if
				} // for

				// sort the new received nodes in temp_recv_buff
				if(new_node_count_ > 1) postorder_sort(0, new_node_count_ - 1, temp_recv_buff);

				// remove the last new_node_count_ nodes from node_list_
				for(int i = 0; i < new_node_count_; ++ i) node_list_.pop_back();

				// find the 'supposed to be' positions of these new nodes
				// (again, may be improve by using binary search ??? )
				std::vector<int> positions;
				for(int i = 0; i < new_node_count_; ++ i) {
					int j = 0;
					while(j < node_list_.size()) {
						if(node_list_[j] < temp_recv_buff[i]) ++ j;
						else break;
					} // while
					positions.push_back(j);
				} // for

				// accordingly increase all the local indices for each position
				for(int i = new_node_count_ - 1; i >= 0; -- i) {
					for(int j = 0; j < node_list_.size(); ++ j) {
						index_type parent = node_list_[j].parent();
						if(parent.proc_id_ == mpi_data_.rank() && parent.index_ >= positions[i])
							node_list_[j].set_parent_index(parent.index_ + 1);
						// do for children also
						std::pair<children_iter, children_iter> children = node_list_[j].children();
						int k = 0;
						for(children_iter iter = children.first;
								iter != children.second; ++ iter, ++ k) {
							if((*iter).proc_id_ == mpi_data_.rank() &&
									(*iter).index_ >= positions[i])
								node_list_[j].set_child_index(k, (*iter).index_ + 1);
						} // for
					} // for
					for(int j = 0; j < total_send_count; ++ j)
						for(int k = 0; k < send_buff[j].num_children_; ++ k)
							if(send_buff[j].children_[k].proc_id_ == mpi_data_.rank() &&
								send_buff[j].children_[k].index_ >= positions[i])
									++ send_buff[j].children_[k].index_;
					for(int j = 0; j < update_indices.size(); ++ j)
						if(update_indices[j].first >= positions[i])
							++ update_indices[j].first;

					positions[i] += i;
				} // for

				// and then insert the un-present received nodes into node_list_
				for(int i = 0; i < new_node_count_; ++ i) {
					iterator pos = node_list_.begin() + positions[i];
					node_list_.insert(pos, temp_recv_buff[i]);
				} // for

				// and construct the new_node_index_ array
				new_node_index_.clear();
				for(int i = 0; i < new_node_count_; ++ i)
					new_node_index_.push_back(positions[i]);

				// construct the notify_array for the positions of the received nodes
				for(int i = 0; i < recv_size; ++ i) {
					tree_node temp_node;
					unpack_node(recv_buff[i], temp_node);
					for(int j = 0; j < node_list_.size(); ++ j)
						if(node_list_[j].small_cell() == temp_node.small_cell())
							notify_array[i] = j;
				} // for

				// communicate back the location of the newly added nodes
				int* to_be_discarded;
				int* recv_back;
				twlib::all_to_all(recv_count, notify_array, to_be_discarded, recv_back,
						recv_size, mpi_data_.comm());

				total_nodes = node_list_.size(); // total number of nodes

				// update parent index info using recv_back
				for(int i = 0; i < update_indices.size(); ++ i)
					node_list_[update_indices[i].first].set_parent_index(
														recv_back[update_indices[i].second]);

				// update the parent and child information of the each new node
				for(int i = 0; i < new_node_count_; ++ i) {
					int index = new_node_index_[i];
					int j = index - 1;

					// the previous node is its child

					node_list_[j].set_parent(make_index(mpi_data_.rank(), index));

					// compute set the large cell information
					uint64_t temp_large_cell = compute_large_cell(
							node_list_[index].small_cell(), node_list_[j].small_cell());
					node_list_[j].set_large_cell(temp_large_cell);

					node_list_[index].add_child(make_index(mpi_data_.rank(), j));

					if(index != total_nodes - 1) {
						uint64_t temp_small_cell = compute_small_cell(
								node_list_[index].small_cell(),	node_list_[index + 1].small_cell());
						j = index + 1; // the parent can only be on the right side
						while(j < total_nodes) {
							if(node_list_[j].small_cell() == temp_small_cell) {
								node_list_[index].set_parent(make_index(mpi_data_.rank(), j));
								temp_large_cell = compute_large_cell(node_list_[j].small_cell(),
										node_list_[index].small_cell());
								node_list_[index].set_large_cell(temp_large_cell);
								break;
							} else ++ j;
						} // while
						if(j == total_nodes) {
							int k = 0;
							while(k < num_oo_nodes) {
								if(send_buff[k].ttd_.small_cell_ == temp_small_cell) {
									node_list_[index].set_parent(
											make_index(proc_array[k], recv_back[k]));
									temp_large_cell = compute_large_cell(
											send_buff[k].ttd_.small_cell_,
											node_list_[index].small_cell());
									node_list_[index].set_large_cell(temp_large_cell);
									break;
								} else ++ k;
							} // while
						} // if
					} else if(mpi_data_.rank() != mpi_data_.size() - 1) {
						uint64_t temp_small_cell = compute_small_cell(
								node_list_[index].small_cell(),
								boundary_array[mpi_data_.rank() + 1].cell_key_);
						int k = 0;
						while(k < num_oo_nodes) {
							if(send_buff[k].ttd_.small_cell_ == temp_small_cell) {
								node_list_[index].set_parent(
										make_index(proc_array[k], recv_back[k]));
								temp_large_cell = compute_large_cell(
										send_buff[k].ttd_.small_cell_,
										node_list_[index].small_cell());
								node_list_[index].set_large_cell(temp_large_cell);
								break;
							} else ++ k;
						} // while
					} else {
						node_list_[index].set_parent(make_index(mpi_data_.rank(), -1));
					} // if-else
				} // for

				// rectify the index information for all remote children:

				// scan the local nodes, and collect info from all those nodes whose parent is remote
				// the info collected is the parent index,
				// and the node's index

				// initialize send counts
				for(int i = 0; i < mpi_data_.size(); ++i) {
					send_count[i] = 0;
					ptr[i] = 0;
				} // for
				// compute send count
				for(int i = 0; i < node_list_.size(); ++ i) {
					index_type parent = node_list_[i].parent();
					if(parent.proc_id_ != mpi_data_.rank())	++ send_count[parent.proc_id_];
				} // for
				total_send_count = send_count[0];
				for(int i = 1; i < mpi_data_.size(); ++ i) {
					ptr[i] = total_send_count;
					total_send_count += send_count[i];
				} // for
				// construct the send_child_buff
				child_data* send_child_buff = new child_data[total_send_count];
				if(send_child_buff == NULL) {
					std::cerr << "Error in memory allocation (globalize_tree())." << std::endl;
					return false;
				} // if
				for(int i = 0; i < node_list_.size(); ++ i) {
					index_type parent = node_list_[i].parent();
					if(parent.proc_id_ != mpi_data_.rank() && parent.index_ != -1) {
						send_child_buff[ptr[parent.proc_id_]].parent_index_ = parent.index_;
						send_child_buff[ptr[parent.proc_id_]].proc_id_ = mpi_data_.rank();
						send_child_buff[ptr[parent.proc_id_]].index_ = i;
						++ ptr[parent.proc_id_];
					} // if
				} // for

				recv_size = 0;
				int* recv_child_count;
				child_data* recv_child_buff;
				twlib::all_to_all(send_count, send_child_buff, recv_child_count, recv_child_buff,
						recv_size, mpi_data_.comm());	// recv_size is total count

				delete[] present_array;
				delete[] recv_buff;
				delete[] recv_count;
				delete[] to_be_discarded;
				delete[] recv_back;
				delete[] send_count;
				delete[] ptr;
				delete[] send_buff;
				delete[] proc_array;
				delete[] notify_array;

				// put the received child node info into multimap
				// key = cell, value = <processor, index>
				typedef std::multimap<const int, std::pair<int, int>, ltindex> child_map;
				typedef typename std::multimap<const int, std::pair<int, int>, ltindex>::iterator
					child_map_iter;

				child_map children_map;
				for(int i = 0; i < recv_size; ++ i) {
					children_map.insert(std::pair<const int, std::pair<int, int> >(
							recv_child_buff[i].parent_index_,
							std::make_pair(recv_child_buff[i].proc_id_, recv_child_buff[i].index_)));
				} // for

				delete[] recv_child_buff;
				delete[] recv_child_count;
				delete[] send_child_buff;

				// for each local node, search for the received cell in it and
				// do the needful (update children)
				int i = 0; // to scan through node_list_ only once
				for(child_map_iter iter = children_map.begin(); iter != children_map.end(); ) {
					// search for (*iter).first in node_list_
					int index = (*iter).first;
					node_list_[index].remove_remote_children(mpi_data_.rank());

					// insert all the new children
					while((*iter).first == index && iter != children_map.end()) {
						node_list_[index].add_child(make_index((*iter).second.first,
									(*iter).second.second));
						++ iter;
					} // while
				} // for

				// do the same for local nodes
				child_map local_children;
				for(int i = 0; i < node_list_.size(); ++ i) {
					index_type parent = node_list_[i].parent();
					if(parent.proc_id_ == mpi_data_.rank() && parent.index_ != -1) {
						local_children.insert(std::pair<const int, std::pair<int, int> >(
							parent.index_,
							std::make_pair(mpi_data_.rank(), i)));
					} // if
				} // for

				for(child_map_iter iter = local_children.begin(); iter != local_children.end(); ) {
					int index = (*iter).first;
					node_list_[index].remove_local_children(mpi_data_.rank());

					// insert all the new children
					while((*iter).first == index && iter != local_children.end()) {
						node_list_[index].add_child(make_index((*iter).second.first,
									(*iter).second.second));
						++ iter;
					} // while
				} // for

				return true;
			} // globalize_tree()


			bool construct_residual_tree(const MPI_data& mpi_data) {
				if(mpi_data.rank() == 0) std::cout << "+ constructing residual tree ... ";

				residual_array_.clear();
				std::vector<int> upward_array;

				// construct upward_array (stores the number of children accumulated)
				for(int i = 0; i < node_list_.size(); ++ i)
					upward_array.push_back(node_list_[i].num_children());

				// perform the local accumulation and construct local residual_array
				int j = 0;
				for(int i = 0; i < node_list_.size(); ++ i) {
					if(upward_array[i] == 0) {
						if(node_list_[i].parent().proc_id_ == mpi_data.rank()) {
							if(node_list_[i].parent().index_ != -1) // is not the root
								-- upward_array[node_list_[i].parent().index_];
						} else
							residual_array_.push_back(i);
					} else
						residual_array_.push_back(i);
				} // for

				is_residual_tree_ = true;
				if(mpi_data.rank() == 0) std::cout << "done." << std::endl;

				return true;
			} // construct_residual_tree()


			bool construct_contraction(const MPI_data& mpi_data) {

				if(!is_residual_tree_) construct_residual_tree(mpi_data);
				if(mpi_data.rank() == 0)
					std::cout << "+ constructing tree contraction structures ... ";

				unsigned int size = node_list_.size();

				// contraction structures to construct
				iteration_array_.resize(size);	
				start_end_.resize(size);
				node_flag_.resize(size);

				std::vector<int> upward_array(size);

				// check sizes
				if(iteration_array_.size() != size || start_end_.size() != size
						|| node_flag_.size() != size || upward_array.size() != size) {
					std::cerr << "Error in memory (construct_contraction())." << std::endl;
					return false;
				} // if

				// initialize the above vectors
				for(int i = 0; i < size; ++ i) {
					iteration_array_[i] = 0;
					start_end_[i] = 'N';
					node_flag_[i] = 'A';
					upward_array[i] = node_list_[i].num_children();
				} // for

				// local accumulation, and mark the finished nodes as 'D'
				for(int i = 0; i < size; ++ i) {
					if(upward_array[i] == 0 && node_list_[i].parent().proc_id_ == mpi_data.rank()) {
						if(node_list_[i].parent().index_ != -1) // is not the root
							-- upward_array[node_list_[i].parent().index_];
						node_flag_[i] = 'D';
						iteration_array_[i] = 0;
					} // if
				} // for

				accumulation_iter_ = 0;
				int* recv_upch = new int[mpi_data.size()];	// same as rcv_upch in Bhaanu's code
				for(int i = 0; i < mpi_data.size(); ++ i) recv_upch[i] = 0;

				// perform the basic accumulation in parallel
				while(1) {
					int i = 0;
					while(i < mpi_data.size() - 1 && recv_upch[i] < 0) ++ i;
					if(i == mpi_data.size() - 1 && recv_upch[i] == 0) break;
					else {
						//if(mpi_data.rank() == 0) std::cout << accumulation_iter_ << std::endl;
						++ accumulation_iter_;
						MPI_Barrier(mpi_data.comm());
						parallel_construct_contraction(upward_array, recv_upch, mpi_data);
					} // if-else
				} // while

				if(mpi_data.rank() == mpi_data.size() - 1 && residual_array_.size() != 0)
					iteration_array_[residual_array_[residual_array_.size() - 1]] =
																			accumulation_iter_ + 1;

				is_accumulation_ready_ = true;
				if(mpi_data.rank() == 0) std::cout << "done." << std::endl;

				delete[] recv_upch;

				return true;
			} // construct_contraction()


			bool parallel_construct_contraction(std::vector<int> &upward_array,
					int* &recv_upch, const MPI_data& mpi_data) {

				int* send_count = new int[mpi_data.size()];
				int* ptr = new int[mpi_data.size()];

				// initialize the arrays
				for(int i = 0; i < mpi_data.size(); ++ i) {
					send_count[i] = 0;
					ptr[i] = 0;
				} // for

				// accumulate leaf residual nodes, and compute send counts and prepare send_buff
				// compute the send counts
				for(int i = 0; i < residual_array_.size(); ++ i) {
					if(node_flag_[residual_array_[i]] != 'D' &&
							upward_array[residual_array_[i]] == 0) {
						index_type parent = node_list_[residual_array_[i]].parent();
						if(parent.proc_id_ != mpi_data.rank())
							++ send_count[parent.proc_id_];
					} // if
				} // for
				int total_count = send_count[0];
				for(int i = 0; i < mpi_data.size(); ++ i) {
					ptr[i] = total_count;
					total_count += send_count[i];
				} // for

				serial_tree_node *send_buff = new serial_tree_node[total_count];
				if(send_buff == NULL) {
					std::cerr << "Error allocating memory (parallel_construct_contraction())."
						<< std::endl;
					return false;
				} // if

				// construct send_buff
				for(int i = 0; i < residual_array_.size(); ++ i) {
					if(node_flag_[residual_array_[i]] != 'D' &&
							upward_array[residual_array_[i]] == 0) {
						index_type parent = node_list_[residual_array_[i]].parent();

						node_flag_[residual_array_[i]] = 'D';
						iteration_array_[residual_array_[i]] = accumulation_iter_;

						if(parent.proc_id_ != mpi_data.rank()) {
							pack_node(send_buff[ptr[parent.proc_id_]],
								node_list_[residual_array_[i]]); // value is also copied
							++ ptr[parent.proc_id_];
						} else {
							if(parent.index_ != -1)
								-- upward_array[parent.index_];
						} // if-else
					} // if
				} // for

				// perform the communication
				serial_tree_node *recv_buff;
				int *recv_count;
				int recv_size = 0;
				twlib::all_to_all(send_count, send_buff, recv_count, recv_buff,	recv_size,
						mpi_data.comm());

				// perform accumulation using received buffer
				for(int i = 0; i < recv_size; ++ i)
					if(recv_buff[i].parent_.index_ != -1)
						-- upward_array[recv_buff[i].parent_.index_];

				int i = 0, upch = 0;
				while(i < residual_array_.size()) {
					if(node_flag_[residual_array_[i]] == 'D') ++ i;
					else break;
				} // while
				if(i == residual_array_.size()) upch = -1;	// all nodes have completed accumulation
				else upch = upward_array[residual_array_[i]];

				// gather upch values from all processors
				MPI_Allgather(&upch, 1, MPI_INT, recv_upch, 1, MPI_INT, mpi_data.comm());

				// Compress leaf chains:

				for(int i = 0; i < residual_array_.size(); ++ i) {
					int j = i + 1;
					if(node_flag_[residual_array_[i]] != 'D') {
						if(upward_array[residual_array_[i]] == 0) {
							while(j < residual_array_.size() && node_flag_[residual_array_[j]] == 'D')
								++ j;
							if(j < residual_array_.size()) {
								if(upward_array[residual_array_[j]] == 1)
									node_flag_[residual_array_[i]] = 'S';
								else
									node_flag_[residual_array_[i]] = 'E';
							} else if(mpi_data.rank() != mpi_data.size() - 1) {
								int k = mpi_data.rank() + 1;
								while(k < mpi_data.size() && recv_upch[k] == -1) ++ k;
								if(k < mpi_data.size() && recv_upch[k] == 1)
									node_flag_[residual_array_[i]] = 'S';
								else
									node_flag_[residual_array_[i]] = 'E';
							} else
								node_flag_[residual_array_[i]] = 'E';
						} else if(upward_array[residual_array_[i]] == 1) { //////////////
							while(j < residual_array_.size() && node_flag_[residual_array_[j]] == 'D')
								++ j;
							if(j < residual_array_.size()) {
								if(upward_array[residual_array_[j]] == 1)
									node_flag_[residual_array_[i]] = 'U';
								else
									node_flag_[residual_array_[i]] = 'E';
							} else if(mpi_data.rank() != mpi_data.size() - 1) {
								int k = mpi_data.rank() + 1;
								while(k < mpi_data.size() && recv_upch[k] == -1) ++ k;
								if(k < mpi_data.size() && recv_upch[k] == 1)
									node_flag_[residual_array_[i]] = 'U';
								else
									node_flag_[residual_array_[i]] = 'E';
							} else
								node_flag_[residual_array_[i]] = 'E';
						} else {
							node_flag_[residual_array_[i]] = 'E'; //////////////////
						} // if-else
					} // if
				} // for

				// segemented parallel prefix for leaf chain compression
				spp_contraction_construction(upward_array, mpi_data);

				// check if accumulation is finished, or more to do
				i = 0;
				while(i < residual_array_.size()) {
					if(node_flag_[residual_array_[i]] == 'D') ++ i;
					else break;
				} // while
				if(i == residual_array_.size()) upch = -1;
				else upch = upward_array[residual_array_[i]];

				// collect recv_upch from all processors
				MPI_Allgather(&upch, 1, MPI_INT, recv_upch, 1, MPI_INT, mpi_data.comm());

				delete[] recv_count;
				delete[] recv_buff;
				delete[] send_buff;
				delete[] ptr;
				delete[] send_count;

				return true;
			} // parallel_construct_contraction()


			bool spp_contraction_construction(std::vector<int> &upward_array,
					const MPI_data& mpi_data) {

				// gather flag of last residual node with 'S' or 'E' on each processor
				char* recv_flag = new char[mpi_data.size()];
				int i = residual_array_.size() - 1;
				char flag = 'N';
				while(i >= 0) {
					if(node_flag_[residual_array_[i]] == 'S' ||
							node_flag_[residual_array_[i]] == 'E') {
						if(node_flag_[residual_array_[i]] == 'S') flag = 'S';
						else flag = 'E';
						break;
					} else
						-- i;
				} // while
				MPI_Allgather(&flag, 1, MPI_CHAR, recv_flag, 1, MPI_CHAR, mpi_data.comm());

				flag = 'N'; ////////////////
				i = mpi_data.rank() - 1;
				while(i >= 0) { ///////////////////
					if(recv_flag[i] == 'S' || recv_flag[i] == 'E') {
						if(recv_flag[i] == 'S') flag = 'N';
						else flag = 'E';
						break;
					} else
						-- i;
				} // while

				bool is_modified = false;
				lp_contraction_construction(upward_array, flag, mpi_data);
				pp_contraction_construction(flag, is_modified, mpi_data);
				if(is_modified) up_contraction_construction(upward_array, mpi_data);

				// mark the nodes in the chain as 'D'
				for(i = 0; i < residual_array_.size(); ++ i) {
					if(node_flag_[residual_array_[i]] != 'D') {
						if(node_flag_[residual_array_[i]] == 'S') {
							node_flag_[residual_array_[i]] = 'D';
							iteration_array_[residual_array_[i]] = accumulation_iter_;
							start_end_[residual_array_[i]] = 'S';
						} // if
						if(node_flag_[residual_array_[i]] == 'U') {
							if(upward_array[residual_array_[i]] == 0) {
								node_flag_[residual_array_[i]] = 'D';
								iteration_array_[residual_array_[i]] = accumulation_iter_;
								start_end_[residual_array_[i]] = 'U';
							} // if
						} // if
					} // if
				} // for

				delete[] recv_flag;

				return true;
			} // spp_contraction_construction()


			bool lp_contraction_construction(std::vector<int> &upward_array,
					char& flag, const MPI_data& mpi_data) {

				int i = 0;
				while(i < residual_array_.size()) {
					if(node_flag_[residual_array_[i]] == 'S') {
						flag = 'S';
						++ i;
						while(i < residual_array_.size() && node_flag_[residual_array_[i]] == 'D')
							++ i;
						if(i < residual_array_.size()) {
							-- upward_array[residual_array_[i]];
							start_end_[residual_array_[i]] = 'E';
						} //if
					} else if(node_flag_[residual_array_[i]] == 'U') {
						++ i;
						while(i < residual_array_.size() && node_flag_[residual_array_[i]] == 'D')
							++ i;
						if(i < residual_array_.size() && flag == 'S') {
							-- upward_array[residual_array_[i]];
							start_end_[residual_array_[i]] = 'E';
						} // if
					} else if(node_flag_[residual_array_[i]] == 'E') {
						flag = 'E';
					 	++ i;
					} else
						++ i;
				} // while

				return true;
			} // lp_contraction_construction()


			bool pp_contraction_construction(char& flag, bool& is_modified, const MPI_data& mpi_data) {

				char* recv_flags = new char[mpi_data.size()];
				MPI_Allgather(&flag, 1, MPI_CHAR, recv_flags, 1, MPI_CHAR, mpi_data.comm());

				is_modified = false;
				int i = 0, j = 0, k = 0;
				while(i < residual_array_.size()) {
					if(node_flag_[residual_array_[i]] == 'S') {
						is_modified = false;
						break;
					} else if(node_flag_[residual_array_[i]] == 'U'	||
							node_flag_[residual_array_[i]] == 'E') {
						j = mpi_data.rank() - 1;
						while(j >= 0) {
							if(recv_flags[j] == 'S') {
								is_modified = true;
								break;
							} else if(recv_flags[j] == 'E') {
								is_modified = false;
								break;
							} else 
								-- j;
						} // while
						break;
					} else
						++ i;
				} // while

				delete[] recv_flags;

				return true;
			} // pp_contraction_construction()


			bool up_contraction_construction(std::vector<int> &upward_array,
					const MPI_data& mpi_data) {

				int i = 0;
				while(i < residual_array_.size()) {
					if(node_flag_[residual_array_[i]] != 'D') {
						if(node_flag_[residual_array_[i]] == 'E') {
							-- upward_array[residual_array_[i]];
							start_end_[residual_array_[i]] = 'E';
							break;
						} else if(node_flag_[residual_array_[i]] == 'U') {
							-- upward_array[residual_array_[i]];
							start_end_[residual_array_[i]] = 'U';
							++ i;
						} else
							break;
					} else
						++ i;
				} // while

				return true;
			} // up_contraction_construction()


			template<typename CombineFunction>
			bool parallel_upward_accumulation(const int& iter, std::vector<int> &upward_array,
					int* &recv_upch, std::vector<char> &start_end, std::vector<char> &node_flag,
					CombineFunction& combine, const MPI_data& mpi_data) {

				int* send_count = new int[mpi_data.size()];
				int* ptr = new int[mpi_data.size()];

				for(int i = 0; i < mpi_data.size(); ++ i) {
					send_count[i] = 0;
					ptr[i] = 0;
				} // for

				// compute send counts
				for(int i = 0; i < residual_array_.size(); ++ i) {
					if(node_flag[residual_array_[i]] != 'D' && upward_array[residual_array_[i]] == 0) {
						index_type parent = node_list_[residual_array_[i]].parent();
						if(parent.proc_id_ != mpi_data.rank())
							++ send_count[parent.proc_id_];
					} // if
				} // for

				int total_count = send_count[0];
				for(int i = 0; i < mpi_data.size(); ++ i) {
					ptr[i] = total_count;
					total_count += send_count[i];
				} // for

				//serial_tree_node *send_buff = new serial_tree_node[18 * mpi_data.size() * levels_];
				serial_tree_node *send_buff = new serial_tree_node[total_count];
				if(send_buff == NULL) {
					std::cerr << "Error allocating memory (parallel_upward_accumulation())."
						<< std::endl;
					return false;
				} // if

				// accumulate leaf residual nodes and prepare send_buff
				for(int i = 0; i < residual_array_.size(); ++ i) {
					if(node_flag[residual_array_[i]] != 'D' && upward_array[residual_array_[i]] == 0) {
						index_type parent = node_list_[residual_array_[i]].parent();

						node_flag[residual_array_[i]] = 'D';

						if(parent.proc_id_ != mpi_data.rank()) {
							pack_node(send_buff[ptr[parent.proc_id_]],
									node_list_[residual_array_[i]]); // value is also copied
							++ ptr[parent.proc_id_];
						} else if(parent.index_ != -1) {
							combine(node_list_[parent.index_], node_list_[residual_array_[i]]);
							-- upward_array[parent.index_];
						} // if-else
					} // if
				} // for

				serial_tree_node *recv_buff;
				int *recv_count;
				int recv_size = 0;
				// perform the communication
				twlib::all_to_all(send_count, send_buff, recv_count, recv_buff,	recv_size,
						mpi_data.comm());
				
				for(int i = 0; i < recv_size; ++ i) {
					tree_node temp_node;
					unpack_node(recv_buff[i], temp_node);
					if(recv_buff[i].parent_.index_ != -1) {
						combine(node_list_[recv_buff[i].parent_.index_], temp_node);
						-- upward_array[recv_buff[i].parent_.index_];
					} // if
				} // for

				int i = 0, upch = 0;
				while(i < residual_array_.size()) {
					if(node_flag[residual_array_[i]] == 'D') ++ i;
					else break;
				} // while
				if(i == residual_array_.size()) upch = -1;
				else upch = upward_array[residual_array_[i]];

				// gather upch values from all processors
				MPI_Allgather(&upch, 1, MPI_INT, recv_upch, 1, MPI_INT, mpi_data.comm());

				// Compress leaf chains:

				for(int i = 0; i < residual_array_.size(); ++ i) {
					int j = i + 1;
					if(node_flag[residual_array_[i]] != 'D') {
						if(upward_array[residual_array_[i]] == 0) {
							while(j < residual_array_.size() && node_flag[residual_array_[j]] == 'D')
								++ j;
							if(j < residual_array_.size()) {
								if(upward_array[residual_array_[j]] == 1)
									node_flag[residual_array_[i]] = 'S';
								else
									node_flag[residual_array_[i]] = 'E';
							} else {
								if(mpi_data.rank() != mpi_data.size() - 1) {
									int k = mpi_data.rank() + 1;
									while(k < mpi_data.size() && recv_upch[k] == -1) ++ k;
									if(k < mpi_data.size() && recv_upch[k] == 1)
										node_flag[residual_array_[i]] = 'S';
									else
										node_flag[residual_array_[i]] = 'E';
								} else
									node_flag[residual_array_[i]] = 'E';
							} // if-else
						} else if(upward_array[residual_array_[i]] == 1) { ///////////
							while(j < residual_array_.size() && node_flag[residual_array_[j]] == 'D')
								++ j;
							if(j < residual_array_.size()) {
								if(upward_array[residual_array_[j]] == 1)
									node_flag[residual_array_[i]] = 'U';
								else
									node_flag[residual_array_[i]] = 'E';
							} else if(mpi_data.rank() != mpi_data.size() - 1) {
									int k = mpi_data.rank() + 1;
									while(k < mpi_data.rank() && recv_upch[k] == -1) ++ k;
									if(k < mpi_data.rank() && recv_upch[k] == 1)
										node_flag[residual_array_[i]] = 'U';
									else
										node_flag[residual_array_[i]] = 'E';
							} else
								node_flag[residual_array_[i]] = 'E';
						} else {
							node_flag_[residual_array_[i]] = 'E'; //////////
						} // if-else
					} // if
				} // for

				// compress leaf chains
				segmented_parallel_prefix(iter, upward_array, start_end, node_flag, combine, mpi_data);

				i = 0;
				while(i < residual_array_.size()) {
					if(node_flag[residual_array_[i]] == 'D') ++ i;
					else break;
				} // while
				if(i == residual_array_.size()) upch = -1;
				else upch = upward_array[residual_array_[i]];

				// collect recv_upch from all processors
				MPI_Allgather(&upch, 1, MPI_INT, recv_upch, 1, MPI_INT, mpi_data.comm());

				delete[] recv_count;
				delete[] recv_buff;
				delete[] send_buff;
				delete[] ptr;
				delete[] send_count;

				return true;
			} // parallel_upward_accumulation()


			template<typename CombineFunction>
			bool segmented_parallel_prefix(const int& iter, std::vector<int> &upward_array,
					std::vector<char> &start_end, std::vector<char> &node_flag,
					CombineFunction& combine, const MPI_data& mpi_data) {

				// gather flag of last node on each processor
				char* recv_flag = new char[mpi_data.size()];
				int i = residual_array_.size() - 1;
				char flag = 'N';
				while(i >= 0) {
					if(node_flag[residual_array_[i]] == 'S'	||
							node_flag[residual_array_[i]] == 'E') {
						if(node_flag[residual_array_[i]] == 'S') flag = 'S';
						else flag = 'E';
						break;
					} else
						-- i;
				} // while
				MPI_Allgather(&flag, 1, MPI_CHAR, recv_flag, 1, MPI_CHAR, mpi_data.comm());

				flag = 'N';
				i = mpi_data.rank() - 1;
				while(i >= 0) {
					if(recv_flag[i] == 'S' || recv_flag[i] == 'E') {
						if(recv_flag[i] == 'S') flag = 'N';
						else flag = 'E';
						break;
					} else
						-- i;
				} // while

				tree_node local_prefix_node, total_prefix_node;
				bool is_modified = false;
				zero_node(total_prefix_node);
				zero_node(local_prefix_node);
				local_prefix(upward_array, start_end, node_flag, total_prefix_node, flag,
						combine, mpi_data);
				parallel_prefix(node_flag, local_prefix_node, total_prefix_node, flag,
						is_modified, combine, mpi_data);
				if(is_modified) {
					update_prefix(upward_array, start_end, node_flag, local_prefix_node,
								is_modified, combine, mpi_data);
				} // if

				// mark the nodes in the chain as 'D'
				for(i = 0; i < residual_array_.size(); ++ i) {
					if(node_flag[residual_array_[i]] != 'D') {
						if(node_flag[residual_array_[i]] == 'S') {
							node_flag[residual_array_[i]] = 'D';
							start_end[residual_array_[i]] = 'S';
						} // if
						if(node_flag[residual_array_[i]] == 'U') {
							if(upward_array[residual_array_[i]] == 0) {
								node_flag[residual_array_[i]] = 'D';
								start_end[residual_array_[i]] = 'U';
							} // if
						} // if
					} // if
				} // for

				delete[] recv_flag;

				return true;
			} // segmented_parallel_prefix()


			template<typename CombineFunction>
			bool local_prefix(std::vector<int> &upward_array, std::vector<char> &start_end,
					std::vector<char> &node_flag,
					tree_node& total_prefix_node, char& flag,
					CombineFunction& combine, const MPI_data& mpi_data) {

				int i = 0, j = 0;
				while(i < residual_array_.size()) {
					if(node_flag[residual_array_[i]] == 'S') {
						flag = 'S';

						combine(total_prefix_node, node_list_[residual_array_[i]]);

						j = i + 1;
						while(j < residual_array_.size() && node_flag[residual_array_[j]] == 'D')
							++ j;
						if(j < residual_array_.size()) {
							combine(node_list_[residual_array_[j]],	node_list_[residual_array_[i]]);
							-- upward_array[residual_array_[j]];
							start_end[residual_array_[j]] = 'E';
						} // if
						i = j;
					} else if(node_flag[residual_array_[i]] == 'U') {
						j = i + 1;
						while(j < residual_array_.size() && node_flag[residual_array_[j]] == 'D')
							++ j;
						if(flag != 'E') {
							combine(total_prefix_node, node_list_[residual_array_[i]]);

							if(j < residual_array_.size()) {
								combine(node_list_[residual_array_[j]],	node_list_[residual_array_[i]]);
								if(flag == 'S') {
									-- upward_array[residual_array_[j]];
									start_end[residual_array_[j]] = 'E';
								} // if
							} // if
						} // if
						i = j;
					} else if(node_flag[residual_array_[i]] == 'E') {
						flag = 'E';
					 	++ i;
					} else
						++ i;
				} // while

				return true;
			} // local_prefix()


			template<typename CombineFunction>
			bool parallel_prefix(std::vector<char> &node_flag, tree_node& local_prefix_node,
					tree_node& total_prefix_node, char& flag, bool& is_modified,
					CombineFunction& combine, const MPI_data& mpi_data) {

				serial_tree_node* recv_buff;
				serial_tree_node* send_buff = new serial_tree_node;

				pack_node(send_buff[0], total_prefix_node);
				int recv_size = 0;
				serial_tree_node* temp_comm_node = new serial_tree_node;
				pack_node(temp_comm_node[0], total_prefix_node);
				twlib::all_gather(temp_comm_node, 1, recv_buff, recv_size, mpi_data.comm());

				char* recv_flags = new char[mpi_data.size()];
				MPI_Allgather(&flag, 1, MPI_CHAR, recv_flags, 1, MPI_CHAR, mpi_data.comm());

				int i = 0, j = 0, k = 0;
				while(i < residual_array_.size()) {
					if(node_flag[residual_array_[i]] == 'S') {
						zero_node(local_prefix_node);
						is_modified = false;
						return true;
					} else if(node_flag[residual_array_[i]] == 'U'	||
								node_flag[residual_array_[i]] == 'E') {
						j = mpi_data.rank() - 1;
						while(j >= 0) {
							if(recv_flags[j] == 'S') {
								tree_node temp_node;
								for(k = j; k < mpi_data.rank() ; ++ k) {
									unpack_node(recv_buff[k], temp_node);
									combine(local_prefix_node, temp_node);
								} // for
								is_modified = true;
								return true;
							} else if(recv_flags[j] == 'E') {
								zero_node(local_prefix_node);
								is_modified = false;
								return true;
							} else 
								-- j;
						} // while
						break;
					} else
						++ i;
				} // while

				delete[] recv_flags;
				delete[] recv_buff;
				delete temp_comm_node;
				delete send_buff;

				return true;
			} // parallel_prefix()


			template<typename CombineFunction>
			bool update_prefix(std::vector<int> &upward_array, std::vector<char> &start_end,
					std::vector<char> &node_flag,
					tree_node& local_prefix_node, bool& is_modified,
					CombineFunction& combine, const MPI_data& mpi_data) {

				int i = 0;
				while(i < residual_array_.size()) {
					if(node_flag[residual_array_[i]] != 'D') {
						if(node_flag[residual_array_[i]] == 'E') {
							combine(node_list_[residual_array_[i]], local_prefix_node);
							-- upward_array[residual_array_[i]];
							start_end[residual_array_[i]] = 'E';
							break;
						} else if(node_flag[residual_array_[i]] == 'U') {
							combine(node_list_[residual_array_[i]], local_prefix_node);
							-- upward_array[residual_array_[i]];
							start_end[residual_array_[i]] = 'U';
							++ i;
						} else
							break;
					} else
						++ i;
				} // while

				return true;
			} // update_prefix()


			template<typename CombineFunction>
			bool parallel_downward_accumulation(const int& iter, const tree_node& default_node,
					std::vector<int> &upward_array,
					CombineFunction& combine, const MPI_data& mpi_data) {

				// Leaf chains are expanded using the information of tree contraction
				// (iteration number and start_end flags)
				reverse_segmented_parallel_prefix(iter, default_node, combine, mpi_data);

				// Leaf nodes obtain parent information for accumulation

				int* send_count = new int[mpi_data.size()];
				int* ptr = new int[mpi_data.size()];
				for(int i = 0; i < mpi_data.size(); ++ i) {
					send_count[i] = 0;
					ptr[i] = 0;
				} // for

				for(int i = 0; i < residual_array_.size(); ++ i) {
					if(iteration_array_[residual_array_[i]] == iter &&
							start_end_[residual_array_[i]] != 'S' &&
							start_end_[residual_array_[i]] != 'U') {
						index_type parent = node_list_[residual_array_[i]].parent();
						if(parent.proc_id_ != mpi_data.rank())
							++ send_count[parent.proc_id_];
					} // if
				} // for

				int send_total = 0;
				for(int i = 0; i < mpi_data.size(); ++ i) {
				   ptr[i] = send_total;
				   send_total += send_count[i];
				}

				serial_tree_node *send_buff = new serial_tree_node[send_total];

				// accumulate and prepare the send_buff
				for(int i = 0; i < residual_array_.size(); ++ i) {
					if(iteration_array_[residual_array_[i]] == iter &&
							start_end_[residual_array_[i]] != 'S' &&
							start_end_[residual_array_[i]] != 'U') {
						index_type parent = node_list_[residual_array_[i]].parent();
						if(parent.proc_id_ == mpi_data.rank()) {
							combine(node_list_[residual_array_[i]],	node_list_[parent.index_]);
						} else {
							pack_node(send_buff[ptr[parent.proc_id_]], node_list_[residual_array_[i]]);
							++ ptr[parent.proc_id_];
						} // if-else
					} // if
				} // for

				// send and receive requests
				serial_tree_node *recv_buff;
				int *recv_count;
				int recv_size = 0;
				// perform the communication
				twlib::all_to_all(send_count, send_buff, recv_count, recv_buff,	recv_size,
						mpi_data.comm());

				// send parent data to remote children residual nodes
				serial_tree_node *temp_send = new serial_tree_node[recv_size];

				for(int i = 0; i < recv_size; ++ i) {
					-- upward_array[recv_buff[i].parent_.index_];
					pack_node(temp_send[i], node_list_[recv_buff[i].parent_.index_]);
				} // for

				serial_tree_node *temp_recv;
				int *temp_recv_count;
				int temp_recv_size = 0;
				// perform the communication
				twlib::all_to_all(recv_count, temp_send, temp_recv_count, temp_recv, temp_recv_size,
						mpi_data.comm());

				int temp_recv_total = 0;
				for(int i = 0; i < mpi_data.size(); ++ i) {
					ptr[i] = temp_recv_total;
					temp_recv_total += temp_recv_count[i];
				} // for

				for(int i = 0; i < residual_array_.size(); ++ i) {
					if(iteration_array_[residual_array_[i]] == iter &&
							start_end_[residual_array_[i]] != 'S' &&
							start_end_[residual_array_[i]] != 'U') {
						index_type parent = node_list_[residual_array_[i]].parent();
						if(parent.proc_id_ != mpi_data.rank()) {
							tree_node temp_node;
							unpack_node(temp_recv[ptr[parent.proc_id_]], temp_node);
							combine(node_list_[residual_array_[i]], temp_node);
							++ ptr[parent.proc_id_];
						} // if
					} // if
				} // for

				delete[] temp_recv_count;
				delete[] temp_recv;
				delete[] temp_send;
				delete[] recv_buff;
				delete[] recv_count;
				delete[] send_buff;
				delete[] ptr;
				delete[] send_count;

				return true;
			} // parallel_downward_accumulation()


			template<typename CombineFunction>
			bool reverse_segmented_parallel_prefix(const int& iter, const tree_node& default_node,
					CombineFunction& combine, const MPI_data& mpi_data) {

				tree_node total_prefix_node, local_prefix_node;
				total_prefix_node = default_node;
				zero_node(local_prefix_node);
				char flag = 'N';
				bool is_modified = false;

				reverse_local_prefix(iter, total_prefix_node, flag,	combine, mpi_data);
				reverse_parallel_prefix(iter, local_prefix_node, total_prefix_node,
						flag, is_modified, combine, mpi_data);
				if(is_modified)
					reverse_update_prefix(iter, local_prefix_node, combine, mpi_data);

				return true;
			} // reverse_segmented_parallel_prefix()


			template<typename CombineFunction>
			bool reverse_local_prefix(const int& iter, tree_node& total_prefix_node, char& flag,
					CombineFunction& combine, const MPI_data& mpi_data) {

				int i = residual_array_.size() - 1;
				while(i >= 0) {
					if(start_end_[residual_array_[i]] == 'E' &&
							iteration_array_[residual_array_[i]] == iter + 1) {
						flag = 'E';
						total_prefix_node = node_list_[residual_array_[i]];

						int j = i - 1;
						while(j >= 0 && start_end_[residual_array_[j]] == 'N') -- j;

						while(j >= 0) {
							if(iteration_array_[residual_array_[j]] == iter &&
									start_end_[residual_array_[j]] != 'E') {
								combine(node_list_[residual_array_[j]], node_list_[residual_array_[i]]);
								break;
							} else
								-- j;
						} // while
						i = j;
					} else if(start_end_[residual_array_[i]] == 'U' &&
								iteration_array_[residual_array_[i]] == iter) {
						flag = 'U';
						total_prefix_node = node_list_[residual_array_[i]];

						int j = i - 1;
						while(j >= 0 && start_end_[residual_array_[j]] == 'N') -- j;

						while(j >= 0) {
							if(iteration_array_[residual_array_[j]] == iter &&
									start_end_[residual_array_[j]] != 'E') {
								if(flag == 'E' || flag == 'U')
									combine(node_list_[residual_array_[j]],
											node_list_[residual_array_[i]]);
								break;
							} else
								-- j;
						} // while
						i = j;
					} else if(start_end_[residual_array_[i]] == 'S' &&
								iteration_array_[residual_array_[i]] == iter) {
						flag = 'S';
						-- i;
					} else
						-- i;
				} // while

				return true;
			} // reverse_local_prefix()


			template<typename CombineFunction>
			bool reverse_parallel_prefix(const int& iter, tree_node& local_prefix_node,
					tree_node& total_prefix_node, char& flag, bool& is_modified,
					CombineFunction& combine, const MPI_data& mpi_data) {

				serial_tree_node* recv_buff;
				int recv_size = 0;
				serial_tree_node* temp_comm_node = new serial_tree_node;
				pack_node(temp_comm_node[0], total_prefix_node);
				twlib::all_gather(temp_comm_node, 1, recv_buff, recv_size, mpi_data.comm());

				char* recv_flags = new char[mpi_data.size()];
				MPI_Allgather(&flag, 1, MPI_CHAR, recv_flags, 1, MPI_CHAR, mpi_data.comm());

				int i = residual_array_.size() - 1;
				is_modified = false;
				while(i >= 0) {
					if(iteration_array_[residual_array_[i]] == iter &&
							(start_end_[residual_array_[i]] == 'S' ||
							start_end_[residual_array_[i]] == 'U')) {
						int j = mpi_data.rank() + 1;
						while(j < mpi_data.size()) {
							if(recv_flags[j] == 'U') {
								tree_node temp_node;
								unpack_node(recv_buff[j], temp_node);
								combine(local_prefix_node, temp_node);
								is_modified = true;
								++ j;
							} else if(recv_flags[j] == 'E') {
								tree_node temp_node;
								unpack_node(recv_buff[j], temp_node);
								combine(local_prefix_node, temp_node);
								is_modified = true;
								break;
							} else
								++ j;
						} // while
						break;
					} else
						-- i;
				} // while

				delete[] recv_flags;
				delete[] recv_buff;
				delete temp_comm_node;

				return true;
			} // reverse_parallel_prefix()


			template<typename CombineFunction>
			bool reverse_update_prefix(const int& iter, tree_node& local_prefix_node,
					CombineFunction& combine, const MPI_data& mpi_data) {
				int i = residual_array_.size() - 1;
				while(i >= 0) {
					if(start_end_[residual_array_[i]] != 'N') {
						if(iteration_array_[residual_array_[i]] == iter &&
								(start_end_[residual_array_[i]] == 'S' ||
								start_end_[residual_array_[i]] == 'U')) {
							combine(node_list_[residual_array_[i]], local_prefix_node);
							if(start_end_[residual_array_[i]] == 'S') break;
							else -- i;
						} else
							break;
					} else
						-- i;
				} // while

				return true;
			} // reverse_update_prefix()


			bool zero_node(tree_node& node) const {
				// improve
				// ...
				value_type zero_value;
				memset((void*)&zero_value, 0, sizeof(value_type));	// value is POD !!!!!!
				node.set_value(zero_value);

				// the following not actually needed
				index_type zero_index;
				node.set_index(zero_index);
				node.set_parent(zero_index);
				node.clear_children();
				node.set_small_cell(0);
				node.set_large_cell(0);
				node.set_level(0);

				return true;
			} // zero_node()


			bool insert_data_points(const OctreeData& octree_data, const int& total_nodes) {

				for(int i = 0, j = 0; i < octree_data.n(); ++ j) {
					if(node_list_[j].is_leaf()) {
						while(octree_data[i].cell_key_ == node_list_[j].small_cell()) {
							node_list_[j].add_point(octree_data[i]);
							++ i;
							if(i >= octree_data.n()) break;
						} // while
					} // if
				} // for

				return true;
			} // insert_data_points()

			
			// construct postordering of the tree rooted at root
			// from temp_nodes into node_list_ from position count
			// root is the index of the root node
			bool construct_postorder(const std::vector<tree_node>& temp_nodes,
					const int& root, int& count) {

				int num_children = temp_nodes[root].num_children();
				if(num_children == 0) {
					tree_node temp1 = temp_nodes[root];
					temp1.clear_children();
					node_list_.push_back(temp1);
					++ count;
					return true;
				} // if

				typedef typename tree_node::const_children_iterator children_iter;
				std::pair<children_iter, children_iter> children = temp_nodes[root].children();

				std::vector<int> children_index;
				for(children_iter iter = children.first; iter != children.second; ++ iter) {
					construct_postorder(temp_nodes, (*iter).index_, count);
					children_index.push_back(count - 1);
				} // for

				// insert temp_nodes[root] and set its parent
				tree_node temp1; // = temp_nodes[root];	// check is this is working fine
				temp1 = temp_nodes[root];
				temp1.clear_children();
				node_list_.push_back(temp1);

				// set the parent child relations
				for(std::vector<int>::iterator iter = children_index.begin();
						iter != children_index.end(); ++ iter) {
					node_list_[*iter].set_parent(make_index(mpi_data_.rank(), count));
					node_list_[count].add_child(make_index(mpi_data_.rank(), *iter));
				} // for
				++ count;

				return true;
			} // construct_postorder()


			int get_root(const std::vector<tree_node>& temp_nodes) {
				int root_index = 0;

				while(temp_nodes[root_index].parent().index_ != -1)
					root_index = temp_nodes[root_index].parent().index_;

				return root_index;
			} // get_root()


			bool insert_leaf(tree_node& curr_leaf, std::vector<tree_node>& temp_nodes,
					const unsigned int i, unsigned int &j, const unsigned int m) {
				// i is the leaf to be inserted (curr_leaf == temp_nodes[i], when i < m)
				// j is the current number of nodes in the tree
				// m is the total number of leaf nodes

				// find the node k whose large cell contains the leaf i
				int k = i - 1;
				while( !is_contained(curr_leaf.small_cell(), temp_nodes[k].large_cell()) )
					k = temp_nodes[k].parent().index_; // tree is local

				index_type temp_index;

				// if leaf i is contained in small cell of k, insert i as child of k
				// else create a new internal node and make i and k the children of this node
				if( is_contained(curr_leaf.small_cell(), temp_nodes[k].small_cell()) ) {
					// node k will become the parent of node i

					// set the large cell of i
					uint64_t large_cell = compute_large_cell(
							temp_nodes[k].small_cell(),	curr_leaf.small_cell());
					curr_leaf.set_large_cell(large_cell);

					// set k as parent of i
					curr_leaf.set_parent(make_index(mpi_data_.rank(), k));

					// set i as child of k (if i is not borrwoed leaf)
					if(i != m)
						temp_nodes[k].add_child(make_index(mpi_data_.rank(), i));

				} else {
					// a new internal node j is created which is made the parent of nodes i and k

					uint64_t small_cell = compute_small_cell(
							curr_leaf.small_cell(), temp_nodes[k].small_cell());

					// make the new node j
					temp_nodes[j].set_small_cell(small_cell);
					temp_nodes[j].set_large_cell(temp_nodes[k].large_cell());

					// make parent of k as parent of j
					temp_index = temp_nodes[k].parent();
					temp_nodes[j].set_parent(temp_index);
					temp_nodes[j].set_parent_proc(mpi_data_.rank());

					// make i and k the children of j
					temp_nodes[j].clear_children();
					temp_nodes[j].add_child(make_index(mpi_data_.rank(), k));
					if(i != m)
						temp_nodes[j].add_child(make_index(mpi_data_.rank(), i));

					// make j the child of k's parent
					if(temp_index.index_ != -1) {
						// check why is k always the last child of its parent ???
						temp_nodes[temp_index.index_].update_last_child(
								make_index(mpi_data_.rank(), j));
					} // if

					// make j the parent of i and k
					uint64_t large_cell = compute_large_cell(
							temp_nodes[j].small_cell(), temp_nodes[i].small_cell());
					curr_leaf.set_large_cell(large_cell);
					curr_leaf.set_parent(make_index(mpi_data_.rank(), j));

					large_cell = compute_large_cell(
							temp_nodes[j].small_cell(), temp_nodes[k].small_cell());
					temp_nodes[k].set_large_cell(large_cell);
					temp_nodes[k].set_parent(make_index(mpi_data_.rank(), j));

					++j; // increase the total number nodes

				} // if-else

				return true;
			} // insert_leaf()


			// sorting for postordering
			// check how to do this using quick sort (currently similar to BHANU'S)
			void postorder_sort(int lower, int upper, std::vector<tree_node> &node_array) const {

				int par, temp_index;
				int  down = lower, up = upper;

				if(down >= up)
					return;

				tree_node* temp_node = new tree_node;

				if(up - down == 1) {	// only 2 elements
					if(node_array[up] < node_array[down]) {
						temp_node[0] = node_array[down];
						node_array[down] = node_array[up];
						node_array[up] = temp_node[0];
					} // if
					delete temp_node;
					return;
				} // if

				if(up - down == 2) {	// only 3 elements
					if(node_array[down + 1] < node_array[down]) {
						temp_node[0] = node_array[down];
						node_array[down] = node_array[down + 1];
						node_array[down + 1] = temp_node[0];
					} // if
					if(node_array[up] < node_array[down + 1]) {
						temp_node[0] = node_array[down + 1];
						node_array[down + 1] = node_array[up];
						node_array[up] = temp_node[0];

						if(node_array[down + 1] < node_array[down]) {
							temp_node[0] = node_array[down];
							node_array[down] = node_array[down + 1];
							node_array[down + 1] = temp_node[0];
						} // if
					} // if
					delete temp_node;
					return;
				} // if

				while(down < up) {
					while((node_array[down] < node_array[lower]) && (down < up)) ++ down;
					while((node_array[up] >= node_array[lower]) && (up > down)) -- up;

					if(down < up) {
						temp_node[0] = node_array[down];
						node_array[down] = node_array[up];
						node_array[up] = temp_node[0];
					} // if
				} // while

				temp_node[0] = node_array[lower];
				node_array[lower] = node_array[up];
				node_array[up] = temp_node[0];

				par = up;

				delete temp_node;

				if(lower < par - 1) postorder_sort(lower, par - 1, node_array);
				if((par + 1) < upper) postorder_sort(par + 1, upper, node_array);

				return;
			} // postorder_sort()


			bool pack_node(serial_tree_node& dest, const tree_node& src) const {

				typedef typename tree_node::const_children_iterator children_iter;
				typedef typename tree_node::const_point_iterator point_iter;

				dest.is_leaf_ = src.is_leaf();
				dest.is_root_ = src.is_root();
				dest.level_ = src.level();
				dest.value_ = src.value(); // value is POD !!!!!!!!
				dest.index_ = src.index();
				dest.parent_ = src.parent();
				std::pair<children_iter, children_iter> children = src.children(); 
				children_iter iter = children.first;
				int i = 0;
				for(; iter != children.second; ++iter, ++i)
					dest.children_[i] = *iter;
				dest.num_children_ = src.num_children();
				if(i > MAX_CHILDREN) {
					std::cout << mpi_data_.rank()
						<< ". ERROR: More children than what should be." << std::endl;
				} // if
				dest.ttd_.small_cell_ = src.small_cell();
				dest.ttd_.large_cell_ = src.large_cell();

				// pack the data points
				dest.num_points_ = src.num_points();
				if(dest.num_points_ != 0) {
					std::pair<point_iter, point_iter> points = src.points(); 
					point_iter p_iter = points.first;
					i = 0;
					for(; p_iter != points.second && i < MAX_POINTS; ++ p_iter, ++ i)
						dest.point_list_[i] = *p_iter;
					if(i > MAX_POINTS) {
						std::cout << "WARNING: More points than allowed. Points may be lost"
							<< std::endl;
						dest.num_points_ = MAX_POINTS;
					} // if
				} // if

				return true;
			} // pack_node()


			bool unpack_node(const serial_tree_node& src, tree_node& dest) const {

				if(src.is_leaf_) dest.make_leaf();
				if(src.is_root_) dest.make_root();
				dest.set_level(src.level_);
				dest.set_value(src.value_); // value is POD !!!!!
				dest.set_index(src.index_);
				dest.set_parent(src.parent_);
				dest.clear_children();
				for(int i = 0; i < src.num_children_; ++ i)
					dest.add_child(src.children_[i]);
				dest.set_small_cell(src.ttd_.small_cell_);
				dest.set_large_cell(src.ttd_.large_cell_);

				// unpack all the points, if leaf
				dest.clear_points();
				for(int i = 0; i < src.num_points_; ++ i) {
					octree_point temp = src.point_list_[i];
					dest.add_point(temp);
				} // for

				return true;
			} // unpack_node()


			// given proc_id and index, create the index_type
			index_type make_index(int p, int i) const {
				index_type new_index;
				new_index.index_ = i;
				new_index.proc_id_ = p;

				return new_index;
			} // make_index()


		private:
			// The nodes of the tree local to a processor
			std::vector<tree_node> node_list_;

			int levels_;		// number of levels of the tree (box size, not levels in tree)
			double domain_;		// domain size of the input data points
			MPI_data mpi_data_;	// put the mpi_data

			// the following are some other data/data structures required for
			// tree construction and computations: improve later !!!!
			std::vector<tree_node> boundaries_;	// vector to store the sorted order of leftmost
												// and rightmost nodes from all processors
			bool is_residual_tree_;			// flag to say if the residual tree is constructed
			int new_node_count_;			// for "new" in Bhaanu's code
			std::vector<int> new_node_index_;	// this is created during tree construction
			std::vector<int> residual_array_;	// contains the local indices of residual
												// nodes. this is created once and used in
												// basic accumulations
			bool is_accumulation_ready_;		// flag to say ifthe below information is constructed
			int accumulation_iter_;				// number of iterations for accumulations
			std::vector<int> iteration_array_;	// for iteration_arr in Bhaanu's code
			std::vector<char> start_end_;		// for start_end in Bhaanu's code
			std::vector<char> node_flag_;		// for node.type in Bhaanu's code

			// the following is to facilitate one-sided communication through RMA
			serial_tree_node *serialized_node_list_;
			MPI_Win node_list_win_;

	}; // BaseTree

} // namespace mpi
} // namespace tw

#endif // MPI_BASE_TREE_HPP
