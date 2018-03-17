#ifndef REDUCED_COMPATIBILITY_LINKED_MATRIX_H_
#define REDUCED_COMPATIBILITY_LINKED_MATRIX_H_

#include <iterator>
#include <algorithm>
#include <vector>
#include <stack>

template <
    typename IndexG,
    typename IndexH>
class reduced_compatibility_linked_matrix {
 private:
  IndexG const m;
  IndexH const n;
  
  struct node {
    node * prev;
    node * next;
    IndexH idx;
    bool active;
    
    node()
        : prev{nullptr},
          next{nullptr},
          active{false} {
    }
  };
  
  std::vector<node> dummy;
  std::vector<node> data;
  
  std::vector<IndexH> count;
  
  std::stack<typename decltype(data)::size_type> history;
  std::stack<typename decltype(history)::size_type> shots;

 public:
  reduced_compatibility_linked_matrix(IndexG m, IndexH n)
      : m{m},
        n{n},
        dummy(m),
        data(m*n),
        count(m) {
    for (IndexG i=0; i<m; ++i) {
      for (IndexH j=0; j<n; ++j) {
        auto idx = i*n + j;
        data[idx].idx = j;
      }
    }
  }
  
  void init() {
    for (IndexG i=0; i<m; ++i) {
      node * prev = &dummy[i];
      for (IndexH j=0; j<n; ++j) {
        auto idx = i*n + j;
        if (get(i, j)) {
          prev->next = &data[idx];
          data[idx].prev = prev;
          prev = &data[idx];
        }
      }
    }
  }

  bool get(IndexG i, IndexH j) const {
    return data[i*n + j].active;
  }
  void set(IndexG i, IndexH j) {
    auto idx = i*n + j;
    data[idx].active = true;
    ++count[i];
  }
  void unset(IndexG i, IndexH j) {
    auto idx = i*n + j;
    if (data[idx].active) {
      history.push(idx);
      data[idx].active = false;
      if (data[idx].prev != nullptr) {
        data[idx].prev->next = data[idx].next;
      }
      if (data[idx].next != nullptr) {
        data[idx].next->prev = data[idx].prev;
      }
      --count[i];
    }
  }

  void advance() {
    shots.push(history.size());
  }
  void revert() {
    auto size = shots.top();
    shots.pop();
    while (history.size() > size) {
      auto idx = history.top();
      if (data[idx].prev != nullptr) {
        data[idx].prev->next = &data[idx];
      }
      if (data[idx].next != nullptr) {
        data[idx].next->prev = &data[idx];
      }
      data[idx].active = true;
			++count[idx/n];
      history.pop();
    }
  }
  
  struct iterator : public std::iterator<std::forward_iterator_tag, IndexG> {
    node * crnt;
    
    iterator(node * crnt)
        : crnt{crnt} {
    }
    
    IndexG & operator*() {
      return crnt->idx;
    }
    
    IndexG & operator++() {
      crnt = crnt->next;
      return crnt->idx;
    }
    
    bool operator!=(iterator const & other) {
      return crnt != other.crnt;
    }
  };
  
  iterator row_begin(IndexG i) {
    return iterator{dummy[i].next};
  }
  
  iterator row_end(IndexG i) {
    return iterator{nullptr};
  }
  
  IndexH num_candidates(IndexG i) const {
    return count[i];
  }
  
};

#endif  // REDUCED_COMPATIBILITY_LINKED_MATRIX_H_
