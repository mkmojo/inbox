/***
 *  Project: TreeWorks
 * 
 *  File: specialized_generate.hpp
 *  Created: Aug 23, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 *
 *  Note: Migrate adjacent_nodes_generate_function.hpp and
 *        interaction_list_generate_function.hpp also to this file.
 */

#ifndef SPECIALIZED_GENERATE_HPP
#define SPECIALIZED_GENERATE_HPP

#include "base_tree.hpp"


class UpwardAccumulationGenerate {
	public:
		template<typename Tree>
		bool operator()(Tree& t, const typename Tree::tree_node& local_node,
				std::back_insert_iterator<std::vector<typename Tree::tree_node> > &out) {

			// do nothing, and return the DEPENDENCY_FLAG as true
			return true;
		} // operator()()
}; // class UpwardAccumulationGenerate


class DownwardAccumulationGenerate {
	public:
		template<typename Tree>
		bool operator()(Tree& t, const typename Tree::tree_node& local_node,
				std::back_insert_iterator<std::vector<typename Tree::tree_node> > &out) {

			// do nothing, and return the DEPENDENCY_FLAG as true
			return true;
		} // operator()()
}; // class UpwardAccumulationGenerate


// generate function that returns the parent node
class DownwardAccumulationGenerate2 {
	public:
		template<typename Tree>
		bool operator()(Tree& t, const typename Tree::tree_node& local_node,
				std::back_insert_iterator<std::vector<typename Tree::tree_node> > &out) {

			typedef typename Tree::tree_node tree_node;
			typedef typename Tree::index_type index_type;

			index_type parent = local_node.parent();
			tree_node parent_node = t[parent];
			//tree_node parent_node = t[local_node.parent()];
			out ++ = parent_node;

			return true;
		} // operator()()
}; // class DownwardAccumulationGenerate


// generate function that returns the children nodes
class UpwardAccumulationGenerate2 {
	public:
		template<typename Tree>
		bool operator()(Tree& t, const typename Tree::tree_node& local_node,
				std::back_insert_iterator<std::vector<typename Tree::tree_node> > &out) {

			typedef typename Tree::tree_node tree_node;
			typedef typename Tree::index_type index_type;
			typedef typename Tree::tree_node::const_children_iterator children_iter;

			std::pair<children_iter, children_iter> children = local_node.children();
			for(children_iter iter = children.first; iter != children.second; ++ iter) {
				tree_node child = t[(*iter)];
				out ++ = child;
			} // for

			return true;
		} // operator()()
}; // class UpwardAccumulationGenerate2


// generate function that returns the node itself (for local computations)
class LocalComputationGenerate {
	public:
		template<typename Tree>
		bool operator()(Tree& t, const typename Tree::tree_node& local_node,
				std::back_insert_iterator<std::vector<typename Tree::tree_node> > &out) {

			out ++ = local_node;

			// the dependency_flag doesnt matter in this case, it can be either true or false
			return false;
		} // operator()()
}; // class LocalComputationGenerate


#endif // SPECIALIZED_GENERATE_HPP
