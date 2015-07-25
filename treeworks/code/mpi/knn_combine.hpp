/***
 *  Project: TreeWorks
 *
 *  File: user_combine_function.hpp
 *  Created: Mar 25, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#ifndef KNN_COMBINE_HPP
#define KNN_COMBINE_HPP

#include "base_tree.hpp"
#include "user_value_type.hpp"
#include "user_input.hpp"
#include "user_head.hpp"

class KNNCombine {
	private:

	public:
		typedef OctreeNodeValueType value_type;
		typedef InputData::Point point_type;
		typedef tw::mpi::BaseTree<value_type, point_type>::tree_node tree_node;

//		template<typename Tree>
//		bool operator()(typename Tree::tree_node& u, const typename Tree::tree_node& v) {
		bool operator()(tree_node& u, const tree_node& v) {

//			typedef OctreeNodeValueType value_type;
//			typedef InputData::Point point_type;
//			typedef typename Tree::tree_node tree_node;

			typedef tree_node::const_children_iterator children_iter;
			typedef tree_node::const_point_iterator point_iter;

			if(u.is_leaf()) {
				std::pair<point_iter, point_iter> u_points = u.points();
				std::pair<point_iter, point_iter> v_points = v.points();

				point_iter i = u_points.first;
				point_iter j = v_points.first;

				// for now assuming only one point per node

				double dist = distance((*i).data_point_.x_, (*i).data_point_.y_, (*i).data_point_.z_,
						(*j).data_point_.x_, (*j).data_point_.y_, (*j).data_point_.z_);

				value_type temp_value = u.value();
				if(dist < temp_value.d_k) {
					if(temp_value.count < K) {
						temp_value.nn_dist[temp_value.count] = dist;
						temp_value.nn_coord[temp_value.count].x_ = (*j).data_point_.x_;
						temp_value.nn_coord[temp_value.count].y_ = (*j).data_point_.y_;
						temp_value.nn_coord[temp_value.count].z_ = (*j).data_point_.z_;

						++ temp_value.count;
					} else { // remove the largest distance point and shift everything
						int i = 0;
						while(temp_value.nn_dist[i] < temp_value.d_k && i < K) ++ i;

						temp_value.nn_dist[i] = dist;
						temp_value.nn_coord[i].x_ = (*j).data_point_.x_;
						temp_value.nn_coord[i].y_ = (*j).data_point_.y_;
						temp_value.nn_coord[i].z_ = (*j).data_point_.z_;

						// find the Kth smallest distance (largest among the ones present)
						for(i = 0; i < K; ++ i) {
							if(temp_value.nn_dist[i] > dist)
								dist = temp_value.nn_dist[i];
						} // for

						temp_value.d_k = dist;
					} // if-else
				} // if

				u.set_value(temp_value);
			} // if

			return true;
		} // operator()()

}; // class KNNCombine

#endif // KNN_COMBINE_HPP
