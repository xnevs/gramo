#ifndef ORDERED_ADJACENCY_MATRIX_H_
#define ORDERED_ADJACENCY_MATRIX_H_

// TODO

class ordered_adjacency_matrix {
 public:
  using index_type = typename std::vector<bool>::size_type;
  
 private:
  index_type n;
  std::vector<bool> mat;
  std::vector<index_type> outdeg;
  std::vector<index_type> indeg;
  
  void set(index_type i, index_type j) {
    mat[i*n + j] = true;
  }
  
  bool get(index_type i, index_type j) const {
    return mat[i*n + j];
  }
  
 public:
  template <
      typename G_,
      typename IndexOrder>
  ordered_adjacency_matrix(
      G_ const & g_,
      IndexOrder const & index_order)
      : n{g_.size()},
        mat(n * n),
        outdeg(n),
        indeg(n) {
    std::vector<index_type> index_pos(n);
    for (index_type i=0; i<n; ++i) {
      index_pos[index_order[i]] = i;
    }
    for (int u=0; u<g_.size(); ++u) {
      auto uu = index_pos[u];
      for(auto v : g_[u]) {
        auto vv = index_pos[v];
        set(uu, vv);
        ++outdeg[uu];
        ++indeg[vv];
      }
    }
  }
  
  index_type num_vertices() const {
    return n;
  }
  
  bool edge(index_type u, index_type v) const {
    return get(u, v);
  }
  
  index_type out_degree(index_type u) const {
    return outdeg[u];
  }
  
  index_type in_degree(index_type u) const {
    return indeg[u];
  }
  
  index_type degree(index_type u) const {
    return out_degree(u) + in_degree(u);
  }
};

#endif  // ORDERED_ADJACENCY_MATRIX_H_
