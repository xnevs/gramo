#ifndef ORDERED_ADJACENCY_LIST_WITH_NOT_AFTER_H_
#define ORDERED_ADJACENCY_LIST_WITH_NOT_AFTER_H_

#include <cstddef>

#include <vector>

#include "adjacency_matrix.h"

class ordered_adjacency_list_with_not_after {
 public:
  using index_type = std::size_t;
 
 private:
  struct node {
    std::vector<index_type> out_before;
    std::vector<index_type> in_before;
    std::vector<index_type> out_after;
    std::vector<index_type> in_after;
    std::vector<index_type> not_out_after;
    std::vector<index_type> not_in_after;
  };
  std::vector<node> nodes;
  
 public:
  template <
      typename G,
      typename IndexOrder>
  ordered_adjacency_list_with_not_after(
      G const & g,
      IndexOrder const & index_order)
      : nodes(g.num_vertices()) {
    auto n = g.num_vertices();
    std::vector<index_type> index_pos(n);
    for (index_type i=0; i<n; ++i) {
      index_pos[index_order[i]] = i;
    }
    
    for (index_type u=0; u<n; ++u) {
      auto u_pos = index_pos[u];
      for (index_type v=0; v<n; ++v) {
        if (g.edge(u, v)) {
          if (u_pos < index_pos[v]) {
            nodes[u].out_after.push_back(v);
            nodes[v].in_before.push_back(u);
          } else {
            nodes[v].in_after.push_back(u);
            nodes[u].out_before.push_back(v);
          }
        } else {
          if (u_pos < index_pos[v]) {
            nodes[u].not_out_after.push_back(v);
          } else {
            nodes[v].not_in_after.push_back(u);
          }
        }
      }
    }
  }
  
  index_type num_vertices() const {
    return nodes.size();
  }
  
  index_type out_degree_before(index_type u) const {
    return nodes[u].out_before.size();
  }
  
  index_type out_degree_after(index_type u) const {
    return nodes[u].out_after.size();
  }
  
  index_type out_degree(index_type u) const {
    return out_degree_before(u) + out_degree_after(u);
  }
  
  index_type in_degree_before(index_type u) const {
    return nodes[u].in_before.size();
  }
  
  index_type in_degree_after(index_type u) const {
    return nodes[u].in_after.size();
  }
  
  index_type in_degree(index_type u) const {
    return in_degree_before(u) + in_degree_after(u);
  }
  
  index_type degree(index_type u) const {
    return out_degree(u) + in_degree(u);
  }
  
  index_type not_out_degree_after(index_type u) const {
    return nodes[u].not_out_after.size();
  }
  
  index_type not_in_degree_after(index_type u) const {
    return nodes[u].not_in_after.size();
  }
  
  std::vector<index_type> const & adjacent_vertices_before(index_type u) const {
    return nodes[u].out_before;
  }
  
  std::vector<index_type> const & adjacent_vertices_after(index_type u) const {
    return nodes[u].out_after;
  }
  
  std::vector<index_type> const & inv_adjacent_vertices_before(index_type u) const {
    return nodes[u].in_before;
  }
  
  std::vector<index_type> const & inv_adjacent_vertices_after(index_type u) const {
    return nodes[u].in_after;
  }
  
  std::vector<index_type> const & not_adjacent_vertices_after(index_type u) const {
    return nodes[u].not_out_after;
  }
  
  std::vector<index_type> const & not_inv_adjacent_vertices_after(index_type u) const {
    return nodes[u].not_in_after;
  }
};

#endif  // ORDERED_ADJACENCY_LIST_WITH_NOT_AFTER_H_
