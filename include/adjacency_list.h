#ifndef ADJACENCY_LIST_H_
#define ADJACENCY_LIST_H_

#include <cstddef>

#include <algorithm>
#include <utility>
#include <vector>

//template <typename Index>
class adjacency_list {
 public:
  using index_type = std::size_t;
  
 private:
  struct node {
    std::vector<index_type> out;
    std::vector<index_type> in;
  };
  std::vector<node> nodes;
  
 public:
  template <typename G_>
  explicit adjacency_list(G_ const & g_)
      : nodes(g_.size()) {
    for(index_type i0=0; i0<g_.size(); ++i0) {
      for(auto i1 : g_[i0]) {
        nodes[i0].out.push_back(i1);
        nodes[i1].in.push_back(i0);
      }
    }
  }
  
  index_type num_vertices() const {
    return nodes.size();
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

#endif  // ADJACENCY_LIST_H_
