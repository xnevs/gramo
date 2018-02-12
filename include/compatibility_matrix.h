#ifndef COMPATIBILITY_MATRIX_H_
#define COMPATIBILITY_MATRIX_H_

#include <iterator>
#include <algorithm>
#include <vector>

template <
    typename IndexG,
    typename IndexH>
class compatibility_matrix {
 private:
  IndexG const m;
  IndexH const n;

  IndexG l;
  std::vector<char> data;

 public:
  compatibility_matrix(IndexG m, IndexH n)
      : m{m},
        n{n},
        l{0},
        data((m+1)*m*n) {
  }

  bool get(IndexG i, IndexH j) const {
    return data[l*m*n + i*n + j];
  }
  void set(IndexG i, IndexH j) {
    data[l*m*n + i*n + j] = true;
  }
  void unset(IndexG i, IndexH j) {
    data[l*m*n + i*n + j] = false;
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
    auto s_it = std::next(std::begin(data),l*m*n);
    auto t_it = std::next(s_it,m*n);
    std::copy_n(s_it,m*n,t_it);
    ++l;
  }
  void revert() {
    --l;
  }

  /*
  void print() {
    std::cout << m << " " << n << std::endl;
   for (int i=0; i<m; ++i) {
      for (int j=0; j<n; ++j) {
        std::cout << (get(i,j) ? '1' : '0');
      }
      std::cout << std::endl;
    }
  }
  */
};

#endif  // COMPATIBILITY_MATRIX_H_
