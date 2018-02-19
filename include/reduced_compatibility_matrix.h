#ifndef REDUCED_COMPATIBILITY_MATRIX_H_
#define REDUCED_COMPATIBILITY_MATRIX_H_

#include <iterator>
#include <algorithm>
#include <vector>
#include <stack>

template <
    typename IndexG,
    typename IndexH>
class reduced_compatibility_matrix {
 private:
  IndexG const m;
  IndexH const n;
  
  std::vector<bool> data;
  
  std::stack<typename decltype(data)::size_type> history;
  std::stack<typename decltype(history)::size_type> shots;

 public:
  reduced_compatibility_matrix(IndexG m, IndexH n)
      : m{m},
        n{n},
        data(m*n) {
  }

  bool get(IndexG i, IndexH j) const {
    return data[i*n + j];
  }
  void set(IndexG i, IndexH j) {
    auto idx = i*n + j;
    data[idx] = true;
  }
  void unset(IndexG i, IndexH j) {
    auto idx = i*n + j;
    if (data[idx]) {
      history.push(idx);
      data[idx] = false;
    }
  }
  
  bool possible(IndexG i) const {
    for (IndexH j=0; j<n; ++j) {
      if (get(i, j)) {
        return true;
      }
    }
    return false;
  }

  void advance() {
    shots.push(history.size());
  }
  void revert() {
    auto size = shots.top();
    shots.pop();
    while (history.size() > size) {
      auto idx = history.top();
      data[idx] = true;
      history.pop();
    }
  }
};

#endif  // REDUCED_COMPATIBILITY_MATRIX_H_
