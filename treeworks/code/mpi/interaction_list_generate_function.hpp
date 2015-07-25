/***
 *  Project: TreeWorks
 *
 *  File: interaction_list_generate_function.hpp
 *  Created: Aug 16, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#ifndef INTERACTION_LIST_GENERATE_FUNCTION_HPP
#define INTERACTION_LIST_GENERATE_FUNCTION_HPP

#include <deque>

#include "base_tree.hpp"

#include "user_value_type.hpp"
#include "user_input.hpp"

// This class is a "user-defined" generate function for the octree to generate the interaction set as
// the interaction-list (nodes that are in farfield of local_node, but their parents are in nearfield
// of local_node's parent.
class InteractionListGenerate {
	public:

		typedef InputData::Point point_type;
		typedef OctreeNodeValueType value_type;
		typedef tw::mpi::BaseTree<value_type, point_type> tree_type;
		typedef tree_type::index_type index_type;
		typedef tree_type::tree_node tree_node;
		typedef tree_node::const_children_iterator children_iter;

		// returns the dependency_flag
		bool operator()(const tree_type& t, const tree_node& local_node,
				std::back_insert_iterator<std::vector<tree_node> > &out) {

			if(local_node.is_root()) return false;

			tree_node curr = t[local_node.parent()];

			double min_cell = t.min_cell_size();
			double x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0, z1 = 0.0, z2 = 0.0;
			t.compute_node_coordinates(curr, x1, x2, y1, y2, z1, z2);
			point_type R1, R2;
			R1.x_ = x1 - (min_cell / 2);
			R1.y_ = y1 - (min_cell / 2);
			R1.z_ = z1 - (min_cell / 2);
			R2.x_ = x2 + (min_cell / 2);
			R2.y_ = y2 + (min_cell / 2);
			R2.z_ = z2 + (min_cell / 2);
			// rectangle formed by R1 and R2 is the rectangle slightly larger than curr
			
			//std::cout << "A" << std::endl;

			while((!curr.is_root()) && (!is_contained(R1, R2, curr, t))) {
				curr = t[curr.parent()];

				//std::cout << "AA" << std::endl;
			} // while

			//std::cout << "B" << std::endl;

			int level = curr.level();
			std::deque<tree_node> P;
			P.push_back(curr);
			std::deque<tree_node> A;
			while(level < local_node.level() - 1) {
				while(P.size() != 0) {
					tree_node v = P[0];
					P.pop_front();

					if(!v.is_leaf()) {
						std::pair<children_iter, children_iter> children = v.children();
						for(children_iter iter = children.first; iter != children.second; ++ iter) {
							tree_node child = t[(*iter)];
							if(intersect(R1, R2, child, t)) A.push_back(child);
							level = child.level();
						} // for
					} else {
						if(intersect(R1, R2, v, t)) A.push_back(v);
					} // if-else
				} // while

				P.clear();
				P = A;
				A.clear();

				int i = 0;
				for(; i < P.size(); ++ i) {
					if(!P[i].is_leaf()) break;
				} // for
				if(i == P.size()) break;
			} // while

			//std::cout << "C" << std::endl;

			while(P.size() != 0) {
				tree_node v = P[0];
				P.pop_front();

				if(!v.is_leaf()) {
					std::pair<children_iter, children_iter> children = v.children();
					for(children_iter iter = children.first; iter != children.second; ++ iter) {
						tree_node child = t[(*iter)];
						if(!t.adjacent(local_node, child)) A.push_back(child);
					} // for
				} else {
					if(!t.adjacent(local_node, v)) A.push_back(v);
				} // if-else
			} // while

			//std::cout << "D" << std::endl;

			while(A.size() != 0) {
				out ++ = A[0];
				A.pop_front();
			} // while

			return false;
		} // operator()()

	private:

		// check is R1-R2 is fully contained in node
		bool is_contained(const point_type& R1, const point_type& R2, const tree_node& node,
				const tree_type& t) const {

			double x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0, z1 = 0.0, z2 = 0.0;
			t.compute_node_coordinates(node, x1, x2, y1, y2, z1, z2);

			if(R1.x_ > x1 && R2.x_ < x2 && R1.y_ > y1 && R2.y_ < y2 && R1.z_ > z1 && R2.z_ < z2)
				return true;

			return false;
		} // is_contained()


		// check if R1-R2 intersect with node
		bool intersect(const point_type& R1, const point_type& R2, const tree_node& node,
				const tree_type& t) const {

			// if R1-R2 is contained in node then they obviously intersect
			if(is_contained(R1, R2, node, t)) return true;

			double x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0, z1 = 0.0, z2 = 0.0;
			t.compute_node_coordinates(node, x1, x2, y1, y2, z1, z2);

			// if node is contained in R1-R2 then they obviously intersect
			if(R1.x_ < x1 && R2.x_ > x2 && R1.y_ < y1 && R2.y_ > y2 && R1.z_ < z1 && R2.z_ > z2)
				return true;

			// in the general case, check if any one of the eight corners of node is inside R1-R2
			if(is_point_contained(x1, y1, z1, R1, R2)) return true;
			if(is_point_contained(x1, y1, z2, R1, R2)) return true;
			if(is_point_contained(x1, y2, z1, R1, R2)) return true;
			if(is_point_contained(x1, y2, z2, R1, R2)) return true;
			if(is_point_contained(x2, y1, z1, R1, R2)) return true;
			if(is_point_contained(x2, y1, z2, R1, R2)) return true;
			if(is_point_contained(x2, y2, z1, R1, R2)) return true;
			if(is_point_contained(x2, y2, z2, R1, R2)) return true;

			return false;
		} // intersect()


		// check if the given point lies with a rectangular region R1-R2
		bool is_point_contained(double x, double y, double z,
				const point_type& R1, const point_type& R2) const {

			if((R1.x_ < x) && (R2.x_ > x) && (R1.y_ < y && R2.y_ > y) && (R1.z_ < z && R2.z_ > z))
				return true;

			return false;
		} // is_point_contained()

}; // class InteractionListGenerate


#endif // INTERACTION_LIST_GENERATE_FUNCTION_HPP
