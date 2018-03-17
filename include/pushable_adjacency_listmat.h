#ifndef PUSHABLE_ADJACENCY_LISTMAT_H_
#define PUSHABLE_ADJACENCY_LISTMAT_H_

#include <iterator>
#include <algorithm>
#include <vector>

#include <boost/range/iterator_range.hpp>

#include "graph_traits.h"

template <typename Index>
class pushable_adjacency_listmat {
 public:
  using directed_category = bidirectional_tag;
  
  using index_type = Index;
  using adjacent_vertices_container_type = std::vector<index_type>;
  
 private:
  index_type n;
 
  struct node {
    using iterator = typename std::vector<index_type>::iterator;
    using const_iterator = typename std::vector<index_type>::const_iterator;
    
    std::vector<index_type> out;
    std::vector<index_type> in;
    
    index_type vis;
    index_type neigh;
    
    node() : vis{0}, neigh{0} {
    }
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
  explicit pushable_adjacency_listmat(G const & g)
      : n{g.num_vertices()},
        nodes(n),
        mat(n * n) {
    for (index_type u=0; u<n; ++u) {
      for (auto v : g.adjacent_vertices(u)) {
        set(u, v);
        nodes[u].out.push_back(v);
        nodes[v].in.push_back(u);
      }
    }
  }
  
  void push(index_type u) {
    bool was_neigh = nodes[u].vis > 0;
    for (auto v : adjacent_vertices(u)) {
      if (was_neigh) {
        --nodes[v].neigh;
      }
      ++nodes[v].vis;
      if (nodes[v].vis == 1) {
        for (auto w : adjacent_vertices(v)) {
          ++nodes[w].neigh;
        }
        for (auto w : inv_adjacent_vertices(v)) {
          ++nodes[w].neigh;
        }
      }
    }
    for (auto v : inv_adjacent_vertices(u)) {
      if (was_neigh) {
        --nodes[v].neigh;
      }
      ++nodes[v].vis;
      if (nodes[v].vis == 1) {
        for (auto w : adjacent_vertices(v)) {
          ++nodes[w].neigh;
        }
        for (auto w : inv_adjacent_vertices(v)) {
          ++nodes[w].neigh;
        }
      }
    }
  }
  
  void pop(index_type u) {
    bool was_neigh = nodes[u].vis > 0;
    for (auto v : adjacent_vertices(u)) {
      if (was_neigh) {
        ++nodes[v].neigh;
      }
      --nodes[v].vis;
      if (nodes[v].vis == 0) {
        for (auto w : adjacent_vertices(v)) {
          --nodes[w].neigh;
        }
        for (auto w : inv_adjacent_vertices(v)) {
          --nodes[w].neigh;
        }
      }
    }
    for (auto v : inv_adjacent_vertices(u)) {
      if (was_neigh) {
        ++nodes[v].neigh;
      }
      --nodes[v].vis;
      if (nodes[v].vis == 0) {
        for (auto w : adjacent_vertices(v)) {
          --nodes[w].neigh;
        }
        for (auto w : inv_adjacent_vertices(v)) {
          --nodes[w].neigh;
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
  
  index_type degree_vis(index_type u) const {
    return nodes[u].vis;
  }
  
  index_type degree_neigh(index_type u) const {
    return nodes[u].neigh;
  }
  
  index_type degree_unv(index_type u) const {
    return degree(u) - degree_neigh(u) - degree_vis(u);
  }
  
  auto adjacent_vertices(index_type u) const {
    return boost::make_iterator_range(
        std::begin(nodes[u].out),
        std::end(nodes[u].out));
  }
  
  auto inv_adjacent_vertices(index_type u) const {
    return boost::make_iterator_range(
        std::begin(nodes[u].in),
        std::end(nodes[u].in));
  }
};

#endif  // PUSHABLE_ADJACENCY_LISTMAT_H_
