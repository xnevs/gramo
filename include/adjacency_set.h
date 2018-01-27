#ifndef ADJACENCY_SET_H_
#define ADJACENCY_SET_H_

#include <cstddef>

#include <algorithm>
#include <utility>
#include <vector>
#include <unordered_set>

//template <typename Index>
class adjacency_set {
 public:
  using index_type = std::size_t;
  
 private:
  struct node {
    std::unordered_set<index_type> out;
    std::unordered_set<index_type> in;
  };
  std::vector<node> nodes;
  
 public:
  template <typename SimpleAdjacencyList>
  adjacency_set(SimpleAdjacencyList const & g)
      : nodes(g.size()) {
    for(index_type i0=0; i0<g.size(); ++i0) {
      for(auto i1 : g[i0]) {
        nodes[i0].out.insert(i1);
        nodes[i1].in.insert(i0);
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
  
  std::unordered_set<index_type> const & adjacent_vertices(index_type u) const {
    return nodes[u].out;
  }
  
  std::unordered_set<index_type> const & inv_adjacent_vertices(index_type u) const {
    return nodes[u].in;
  }
};

#endif  // ADJACENCY_SET_H_
