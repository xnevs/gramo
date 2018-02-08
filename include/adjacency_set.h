#ifndef ADJACENCY_SET_H_
#define ADJACENCY_SET_H_

#include <vector>
#include <set>

template <typename Index>
class adjacency_set {
 public:
  using index_type = Index;
  
 private:
  struct node {
    std::set<index_type> out;
    std::set<index_type> in;
  };
  std::vector<node> nodes;
  
 public:
  template <typename G>
  explicit adjacency_set(G const & g)
      : nodes(g.num_vertices()) {
    auto n = g.num_vertices();
    for(index_type u=0; u<n; ++u) {
      for(auto v : g.adjacent_vertices(u)) {
        nodes[u].out.insert(v);
        nodes[v].in.insert(u);
      }
    }
  }
  
  index_type num_vertices() const {
    return nodes.size();
  }
  
  bool edge(index_type u, index_type v) const {
    return nodes[u].out.find(v) != nodes[u].out.end();
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
  
  std::set<index_type> const & adjacent_vertices(index_type u) const {
    return nodes[u].out;
  }
  
  std::set<index_type> const & inv_adjacent_vertices(index_type u) const {
    return nodes[u].in;
  }
};

#endif  // ADJACENCY_SET_H_
