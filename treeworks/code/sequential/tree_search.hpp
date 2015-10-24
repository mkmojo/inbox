/***
 *  $Id: tree_search.hpp 63 2009-03-15 22:19:34Z asarje $
 **
 *  File: tree_search.hpp
 *  Created: Mar 11, 2009
 */

#ifndef SEQUENTIAL_TREE_SEARCH_HPP
#define SEQUENTIAL_TREE_SEARCH_HPP

#include <iterator>
#include <tuple>
#include <vector>

namespace tw
{
	namespace seq
	{

		template <typename Tree,
			typename QueryPredicate,
			typename SelectPredicate,
			typename OutputIterator>
		bool search_tree_base_(const Tree & t, const Tree::node_type u,
			QueryPredicate q, SelectPredicate pu,
			SelectPredicate pv, OutputIterator out) {

			std::vector<index_type> vlst;
			if(pu(t[u], q) == true) {
				*out++ = u;
				return true;
			} else {
				const_children_iterator i1, i2;
				std::tr1::tie(i1, i2) = u.children();
				for(; i1 != i2; ++i1) {
					index_type v = *i1;
					if(pv(t[v], q) == true)
						vlst.push_back(v);
				}
			}

			std::vector<index_type>::const_iterator i(vlst.begin ());
			for (; i != vlst.end(); ++i)
				search_tree_base_(t, *i, q, pu, pv, out);

			return true;
		} // search_tree_base_


		/**
		*  This function is a model of TreeSearchFunction.
		*  It performs sequential search on the tree @a t for the list of query
		*  predicates [@a first, @a last). [TODO: explain pu, pv, out]
		*  @param Tree must be a model of BaseTree.
		*  @param QueryIterator must be a model of ForwardIterator to a sequence
		*  of predicates modeling QueryPredicate.
		*  @param SelectPredicate must be a model of SelectPredicate.
		*  @param OutputIterator must be a model of OutputIterator to a sequence
		*  of objects modeling BackInsertionSequence.
		*  @return true on success, false otherwise.
		*  [TODO: do we need SelectPredicate, or is SelectPredicate equivalent
		*  to QueryPredicate? Is there any better way to return output than by
		*  "sequence of sequences"?]
		*/
		template <typename Tree,
				typename QueryIterator,
				typename SelectPredicate,
				typename OutputIterator>
		bool tree_search(const Tree & t, QueryIterator first,
				QueryIterator last, SelectPredicate pu,
				SelectPredicate pv, OutputIterator out) {

			typedef Tree::index_type index_type;
			typedef Tree::node_type node_type;

			typedef node_type::const_children_iterator const_children_iterator;

			for(; first != last; ++first, ++out) {
				index_type u = t.root();

				std::vector<index_type> vlst;
				std::back_insert_iterator<std::vector<index_type> >ivlst(vlst);

				if(search_tree_base_(t, u, *first, pu, pv, ivlst) == false)
					return false;
				for (std::vector<index_type>::const_iterator i(vlst.begin());
						i != vlst.end(); ++i)
					(*out).push_back(*i);
			}
			return true;
		} // search_tree

	} // namespace seq
} // namespace tw

#endif // SEQUENTIAL_TREE_SEARCH_HPP
