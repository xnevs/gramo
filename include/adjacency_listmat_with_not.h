#ifndef ADJACENCY_LISTMAT_WITH_NOT_H_
#define ADJACENCY_LISTMAT_WITH_NOT_H_

#include <vector>

#include "graph_traits.h"

template <typename Index>
class adjacency_listmat_with_not {
 public:
  using directed_category = bidirectional_tag;
  
  using index_type = Index;
  using adjacent_vertices_container_type = std::vector<index_type>;
  
 private:
  index_type n;
 
  struct node {
    std::vector<index_type> out;
    std::vector<index_type> in;
    std::vector<index_type> not_out;
    std::vector<index_type> not_in;
  };
  std::vector<node> nodes;
  
  std::vector<bool> mat;
  
  void set(index_type i, index_type j) {
    mat[i*n + j] = true;
  }
  
  bool get(index_type i, index_type j) const {
    return mat[i*n + j];
  }
  
 public:
  template <typename G>
  explicit adjacency_listmat_with_not(G const & g)
      : n{g.num_vertices()},
        nodes(n),
        mat(n * n) {
    for (index_type u=0; u<n; ++u) {
      for (auto v : g.adjacent_vertices(u)) {
        nodes[u].out.push_back(v);
        nodes[v].in.push_back(u);
        set(u, v);
      }
    }
    for (index_type u=0; u<n; ++u) {
      for (index_type v=u+1; v<n; ++v) {
        if (!get(u, v)) {
          nodes[u].not_out.push_back(v);
          nodes[v].not_in.push_back(u);
        }
        if (!get(v, u)) {
          nodes[u].not_in.push_back(v);
          nodes[v].not_out.push_back(u);
        }
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
    return nodes[u].out.size();
  }
  
  index_type in_degree(index_type u) const {
    return nodes[u].in.size();
  }
  
  index_type degree(index_type u) const {
    return out_degree(u) + in_degree(u);
  }
  
  std::vector<index_type> const & adjacent_vertices(index_type u) const {
    return nodes[u].out;
  }
  
  std::vector<index_type> const & inv_adjacent_vertices(index_type u) const {
    return nodes[u].in;
  }
  
  std::vector<index_type> const & not_adjacent_vertices(index_type u) const {
    return nodes[u].not_out;
  }
  
  std::vector<index_type> const & not_inv_adjacent_vertices(index_type u) const {
    return nodes[u].not_in;
  }
};

#endif  // ADJACENCY_LISTMAT_WITH_NOT_H_
