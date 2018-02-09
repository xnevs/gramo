#ifndef ORDRED_ADJACENCY_LISTMAT_H_
#define ORDRED_ADJACENCY_LISTMAT_H_

#include <vector>

template <typename Index>
class ordered_adjacency_listmat {
 public:
  using index_type = Index;
  
 private:
  index_type n;
 
  struct node {
    std::vector<index_type> out_before;
    std::vector<index_type> in_before;
  };
  std::vector<node> nodes;
  
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
      typename G,
      typename IndexOrder>
  ordered_adjacency_listmat(
      G const & g,
      IndexOrder const & index_order)
      : n{g.num_vertices()},
        nodes(n),
        mat(n * n),
        outdeg(n),
        indeg(n) {
    std::vector<index_type> index_pos(n);
    for (index_type i=0; i<n; ++i) {
      index_pos[index_order[i]] = i;
    }
    for (index_type u=0; u<n; ++u) {
      for (auto v : g.adjacent_vertices(u)) {
        if (index_pos[u] > index_pos[v]) {
          nodes[u].out_before.push_back(v);
        } else {
          nodes[v].in_before.push_back(u);
        }
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
  
  index_type out_degree_before(index_type u) const {
    return nodes[u].out_before.size();
  }
  
  index_type in_degree_before(index_type u) const {
    return nodes[u].in_before.size();
  }
  
  index_type degree_before(index_type u) const {
    return out_degree_before(u) + in_degree_before(u);
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
  
  std::vector<index_type> const & adjacent_vertices_before(index_type u) const {
    return nodes[u].out_before;
  }
  
  std::vector<index_type> const & inv_adjacent_vertices_before(index_type u) const {
    return nodes[u].in_before;
  }
};

#endif  // ORDERED_ADJACENCY_LISTMAT_H_
