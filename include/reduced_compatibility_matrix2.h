#ifndef REDUCED_COMPATIBILITY_MATRIX2_H_
#define REDUCED_COMPATIBILITY_MATRIX2_H_

#include <iterator>
#include <algorithm>
#include <vector>
#include <stack>

template <
    typename IndexG,
    typename IndexH>
class reduced_compatibility_matrix2 {
 private:
  IndexG const m;
  IndexH const n;
  
  std::vector<char> data;
  
  std::vector</*typename decltype(data)::size_type*/int> history;
  /*typename decltype(history)::size_type*/ int index;
  std::vector</*typename decltype(history)::size_type*/int> shots;
  /*typename decltype(shots)::size_type*/ int shotidx;

 public:
  reduced_compatibility_matrix2(IndexG m, IndexH n)
      : m{m},
        n{n},
        data(m*n),
        history(m*n),
        index{-1},
        shots(m*n),
        shotidx{-1} {
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
      data[idx] = false;
      history[++index] = idx;
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
    shots[++shotidx] = index;
  }
  void revert() {
    int stop = shots[shotidx--];
    while (index > stop) {
			data[history[index--]] = true;
    }
  }
};

#endif  // REDUCED_COMPATIBILITY_MATRIX2_H_
