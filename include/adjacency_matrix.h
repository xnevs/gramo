#ifndef ADJACENCY_MATRIX_H_
#define ADJACENCY_MATRIX_H_

class adjacency_matrix {
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
  template <typename SimpleAdjacencyList>
  adjacency_matrix(SimpleAdjacencyList const & g)
      : n{g.size()},
        mat(n * n),
        outdeg(n),
        indeg(n) {
    for (int u=0; u<g.size(); ++u) {
      for(auto v : g[u]) {
        set(u, v);
        ++outdeg[u];
        ++indeg[v];
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

#endif  // ADJACENCY_MATRIX_H_
