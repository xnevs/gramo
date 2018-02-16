#ifndef ADJACENCY_LIST_IN_ORDER_H_
#define ADJACENCY_LIST_IN_ORDER_H_

#include <vector>
#include <algorithm>

template <typename Index>
class adjacency_list_in_order {
 public:
  using index_type = Index;
  
 private:
  struct node {
    std::vector<index_type> out;
    std::vector<index_type> in;
  };
  std::vector<node> nodes;
  
 public:
  template <typename G>
  explicit adjacency_list_in_order(G const & g)
      : nodes(g.num_vertices()) {
    auto n = g.num_vertices();
    for (index_type u=0; u<n; ++u) {
      for (auto v : g.adjacent_vertices(u)) {
        nodes[u].out.push_back(v);
        nodes[v].in.push_back(u);
      }
    }
    for (auto & node : nodes) {
      std::sort(std::begin(node.out), std::end(node.out));
      std::sort(std::begin(node.in), std::end(node.in));
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
  
  std::vector<index_type> const & adjacent_vertices_in_order(index_type u) const {
    return nodes[u].out;
  }
  
  std::vector<index_type> const & inv_adjacent_vertices_in_order(index_type u) const {
    return nodes[u].in;
  }
};

#endif  // ADJACENCY_LIST_IN_ORDER_H_
