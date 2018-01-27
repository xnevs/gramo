#ifndef COMPATIBILITY_MATRIX_H_
#define COMPATIBILITY_MATRIX_H_

#include <vector>
#include <algorithm>

template <
    typename IndexSmall,
    typename IndexLarge>
class compatibility_matrix {
 public:
  IndexSmall const m;
  IndexLarge const n;

  IndexSmall l;
  std::vector<char> data;

  compatibility_matrix(IndexSmall m, IndexLarge n)
      : m{m},
        n{n},
        l{0},
        data((m+1)*m*n) {
  }

  bool get(IndexSmall i, IndexLarge j) const {
    return data[l*m*n + i*n + j];
  }
  void set(IndexSmall i, IndexLarge j) {
    data[l*m*n + i*n + j] = true;
  }
  void unset(IndexSmall i, IndexLarge j) {
    data[l*m*n + i*n + j] = false;
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
