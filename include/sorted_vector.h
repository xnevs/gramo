#ifndef SORTED_VECTOR_H_
#define SORTED_VECTOR_H_

#include <algorithm>
#include <vector>

template <
    typename T,
    typename Compare = std::less<T>>
struct sorted_vector {
  std::vector<T> V;
  Compare cmp;
  
  using size_type = typename std::vector<T>::size_type;
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;
  
  sorted_vector(const Compare& c = Compare())
      : cmp(c) {
  }
  
  template <class InputIterator>
  sorted_vector(InputIterator first, InputIterator last, Compare const & cmp = Compare())
      : V(first, last), cmp{cmp} {
    std::sort(begin(), end(), cmp);
  }
  
  size_type size() const {
    return V.size();
  }
  
  iterator begin() {
    return V.begin();
  }
  
  iterator end() {
    return V.end();
  }
  
  const_iterator begin() const {
    return V.begin();
  }
  
  const_iterator end() const {
    return V.end();
  }
  
  iterator insert(const T& t) {
    iterator i = lower_bound(begin(), end(), t, cmp);
    if (i == end() || cmp(t, *i))
      V.insert(i, t);
    return i;
  }
  
  iterator find(const T& t) {
    iterator i = lower_bound(begin(), end(), t, cmp);
    return i == end() || cmp(t, *i) ? end() : i;
  }
  
  void erase(iterator pos) {
    V.erase(pos);
  }
};

#endif  // SORTED_VECTOR_H_
