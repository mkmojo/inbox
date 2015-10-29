/***
 *  $Id: functional.hpp 59 2009-03-11 21:23:44Z zola $
 **
 *  File: functional.hpp
 *  Created: Mar 11, 2009
 */

#ifndef FUNCTIONAL_HPP
#define FUNCTIONAL_HPP

#include <functional>


namespace tw {

  template <typename T>
  struct query_predicate : public std::unary_function<T, bool> { };

  template <typename T, template <typename> class Pred>
  struct selection_predicate : public std::binary_function<T, Pred<T>, bool> { };

} // namespace tw

#endif // FUNCTIONAL_HPP
