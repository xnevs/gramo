#ifndef PACKED_COMPATIBILITY_MATRIX_H_
#define PACKED_COMPATIBILITY_MATRIX_H_

#include <vector>

template <
    typename IndexG,
    typename IndexH>
class packed_compatibility_matrix {
 protected:
  using bin_type = unsigned char;
 public:
  IndexG const m;
  IndexH const n;
  
  decltype((m*n+sizeof(bin_type))/sizeof(bin_type)) const frame_size;

  IndexG l;
  std::vector<bin_type> data;

  packed_compatibility_matrix(IndexG m, IndexH n)
      : m{m},
        n{n},
        frame_size{(m*n+sizeof(bin_type))/sizeof(bin_type)},
        l{0},
        data((m+1)*frame_size) {
  }

  bool get(IndexG i, IndexH j) const {
    auto idx = i*n + j;
    return static_cast<bool>(data[l*frame_size + (idx / sizeof(bin_type))] & (static_cast<bin_type>(1) << (idx % sizeof(bin_type))));
  }
  void set(IndexG i, IndexH j) {
    auto idx = i*n + j;
    data[l*frame_size + (idx / sizeof(bin_type))] |= static_cast<bin_type>(1) << (idx % sizeof(bin_type));
  }
  void unset(IndexG i, IndexH j) {
    auto idx = i*n + j;
    data[l*frame_size + (idx / sizeof(bin_type))] &= ~(static_cast<bin_type>(1) << (idx % sizeof(bin_type)));
  }

  void advance() {
    auto s_it = std::next(std::begin(data),l*frame_size);
    auto t_it = std::next(s_it,frame_size);
    std::copy_n(s_it,frame_size,t_it);
    ++l;
  }
  void revert() {
    --l;
  }
};

#endif  // PACKED_COMPATIBILITY_MATRIX_H_
