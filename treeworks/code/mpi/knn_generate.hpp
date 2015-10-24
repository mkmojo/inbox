/***
 *  Project: TreeWorks
 * 
 *  File: knn_generate.hpp
 *  Created: Sep 6, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#ifndef KNN_GENERATE_HPP
#define KNN_GENERATE_HPP

#include "base_tree.hpp"

#include <climits>
#include <deque>
#include <set>

#include "user_head.hpp"

// generate function that returns the near nodes for KNN computations
class KNNGenerate {
	public:
		template<typename Tree>
		bool operator()(Tree& t, const typename Tree::tree_node& local_node,
				std::back_insert_iterator<std::vector<typename Tree::tree_node> > &out) {

			typedef typename Tree::tree_node tree_node;
			typedef typename Tree::index_type index_type;
			typedef typename Tree::tree_node::const_children_iterator children_iter;

			if(!local_node.is_leaf()) return false;

			tree_node curr;
			curr = local_node;

			// get the root
			while(!curr.is_root()) {
				curr = t[curr.parent()];
			} // while

			std::deque<tree_node> P;
			P.push_back(curr);
			std::deque<tree_node> A;

			// perform downward traversal
			while(1) {
				std::multiset<double> large_dist;
				for(int i = 0; i < P.size(); ++ i)
					large_dist.insert(t.largest_distance(curr, P[i]));
				unsigned int d_k = 0.0;
				if(P.size() < K) d_k = UINT_MAX;
				else {
					std::multiset<double>::iterator i = large_dist.begin();
					int j = 0;
					while(j != K - 1) { ++ i; ++ j; }
					d_k = (unsigned int) *i;
				} // if-else

//				for(std::multiset<double>::iterator i = large_dist.begin(); i != large_dist.end(); ++ i)
//					std::cout << *i << " ";
//				std::cout << std::endl;
//				std::cout << "*** C: d_k = " << d_k << std::endl;

				while(P.size() != 0) {
					tree_node w = P[0];

//					std::cout << "P = " << P.size() << ", A = " << A.size()
//						<< ", sm = " << t.smallest_distance(curr, w)
//						<< ", la = " << t.largest_distance(curr, w) << ", d_k = " << d_k
//						<< ", curr = " << t.cell_size(curr) << ", w = " << t.cell_size(w)
//						<< ", curr = " << curr.small_cell() << ", w = " << w.small_cell() << std::endl;

					P.pop_front();

					if(t.smallest_distance(curr, w) < d_k) {
						if(t.cell_size(curr) > t.cell_size(w) || w.is_leaf()) {
							A.push_back(w);
						} else { // add children of w to P
							std::pair<children_iter, children_iter> children = w.children();
							for(children_iter iter = children.first; iter != children.second; ++ iter) {
								tree_node child = t[(*iter)];
								P.push_back(child);
							} // for
						} // if-else
					} // if
				} // while

				if(curr == local_node || curr.is_leaf()) break;

				std::pair<children_iter, children_iter> children = curr.children();
				for(children_iter iter = children.first; iter != children.second; ++ iter) {
					tree_node child = t[(*iter)];
					if(child == local_node || child.is_ancestor(local_node)) {
						curr = child;
//						std::cout << "-----------------------------------------------------------\n";
						break;
					} // if
				} // for

				P.clear();
				P = A;
				A.clear();
			} // while

//			std::cout << "*** D: P = " << P.size() << ", " << A.size() << std::endl;

			while(A.size() != 0) {
				out ++ = A[0];
				A.pop_front();
			} // while

			// dependency flag is false
			return false;
		} // operator()()
}; // class KNNGenerate


#endif // KNN_GENERATE_HPP
