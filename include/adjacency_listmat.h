#ifndef ADJACENCY_LISTMAT_H_
#define ADJACENCY_LISTMAT_H_

#include <cstddef>

#include <algorithm>
#include <utility>
#include <vector>

//template <typename Index>
class adjacency_listmat {
 public:
  using index_type = std::size_t;
  
 private:
  index_type n;
 
  struct node {
    std::vector<index_type> out;
    std::vector<index_type> in;
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
  template <typename G_>
  explicit adjacency_listmat(G_ const & g_)
      : n{g_.size()},
        nodes(n),
        mat(n * n) {
    for(index_type u=0; u<g_.size(); ++u) {
      for(auto v : g_[u]) {
        nodes[u].out.push_back(v);
        nodes[v].in.push_back(u);
        set(u, v);
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
};

#endif  // ADJACENCY_LISTMAT_H_
