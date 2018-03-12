#ifndef ORDERABLE_ADJACENCY_LISTMAT_WITH_RI_DEGREE_H_
#define ORDERABLE_ADJACENCY_LISTMAT_WITH_RI_DEGREE_H_

#include <iterator>
#include <algorithm>
#include <vector>

#include <boost/range/iterator_range.hpp>

#include "graph_traits.h"

template <typename Index>
class orderable_adjacency_listmat_with_ri_degree {
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
    iterator out_mid;
    iterator in_mid;
    
    const_iterator out_cmid() const {
      return out_mid;
    }
    const_iterator in_cmid() const {
      return in_mid;
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
  explicit orderable_adjacency_listmat_with_ri_degree(G const & g)
      : n{g.num_vertices()},
        nodes(n),
        mat(n * n) {
    for (index_type u=0; u<n; ++u) {
      for (auto v : g.adjacent_vertices(u)) {
        set(u, v);
      }
    }
    for (index_type u=0; u<n; ++u) {
      for (index_type v=0; v<n; ++v) {
        if (get(u, v)) {
          nodes[u].out.push_back(v);
          nodes[v].in.push_back(u);
        }
      }
    }
    for (index_type u=0; u<n; ++u) {
      std::sort(std::begin(nodes[u].out), std::end(nodes[u].out));
      std::sort(std::begin(nodes[u].in), std::end(nodes[u].in));
      nodes[u].out_mid = std::begin(nodes[u].out);
      nodes[u].in_mid = std::begin(nodes[u].in);
    }
  }
  
  void push(index_type u) {
    // TODO watch out for loops
    for (auto v : adjacent_vertices_after(u)) {
      auto u_it = std::lower_bound(
          nodes[v].in_mid,
          std::end(nodes[v].in),
          u);
      std::rotate(nodes[v].in_mid, u_it, std::next(u_it));
      ++nodes[v].in_mid;
    }
    for (auto v : inv_adjacent_vertices_after(u)) {
      auto u_it = std::lower_bound(
          nodes[v].out_mid,
          std::end(nodes[v].out),
          u);
      std::rotate(nodes[v].out_mid, u_it, std::next(u_it));
      ++nodes[v].out_mid;
    }
  }
  
  void pop(index_type u) {
    // TODO watch out for loops
    for (auto v : adjacent_vertices_after(u)) {
      auto it = std::lower_bound(
          nodes[v].in_mid,
          std::end(nodes[v].in),
          u);
      std::rotate(std::prev(nodes[v].in_mid), nodes[v].in_mid, it);
      --nodes[v].in_mid;
    }
    for (auto v : inv_adjacent_vertices_after(u)) {
      auto it = std::lower_bound(
          nodes[v].out_mid,
          std::end(nodes[v].out),
          u);
      std::rotate(std::prev(nodes[v].out_mid), nodes[v].out_mid, it);
      --nodes[v].out_mid;
    }
  }
  
  index_type num_vertices() const {
    return n;
  }
  
  bool edge(index_type u, index_type v) const {
    return get(u, v);
  }
  
  index_type out_degree_before(index_type u) const {
    return adjacent_vertices_before(u).size();
  }
  index_type in_degree_before(index_type u) const {
    return inv_adjacent_vertices_before(u).size();
  }
  index_type degree_before(index_type u) const {
    return out_degree_before(u) + in_degree_before(u);
  }
  
  index_type out_degree_after(index_type u) const {
    return adjacent_vertices_after(u).size();
  }
  index_type in_degree_after(index_type u) const {
    return inv_adjacent_vertices_after(u).size();
  }
  index_type degree_after(index_type u) const {
    return out_degree_after(u) + in_degree_after(u);
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
  
  index_type out_degree_after_neigh(index_type u) const {
    index_type deg = 0;
    for (auto v : adjacent_vertices_after(u)) {
      if (degree_before(v) > 0) {
        ++deg;
      }
    }
    return deg;
  }
  index_type in_degree_after_neigh(index_type u) const {
    index_type deg = 0;
    for (auto v : inv_adjacent_vertices_after(u)) {
      if (degree_before(v) > 0) {
        ++deg;
      }
    }
    return deg;
  }
  index_type degree_after_neigh(index_type u) const {
    return out_degree_after_neigh(u) + in_degree_after_neigh(u);
  }
  
  auto adjacent_vertices(index_type u) const {
    return boost::make_iterator_range(
        std::begin(nodes[u].out),
        std::end(nodes[u].out));
  }
  
  auto adjacent_vertices_before(index_type u) const {
    return boost::make_iterator_range(
        std::begin(nodes[u].out),
        nodes[u].out_cmid());
  }
  
  auto adjacent_vertices_after(index_type u) const {
    return boost::make_iterator_range(
        nodes[u].out_cmid(),
        std::end(nodes[u].out));
  }
  
  auto inv_adjacent_vertices(index_type u) const {
    return boost::make_iterator_range(
        std::begin(nodes[u].in),
        std::end(nodes[u].in));
  }
  
  auto inv_adjacent_vertices_before(index_type u) const {
    return boost::make_iterator_range(
        std::begin(nodes[u].in),
        nodes[u].in_cmid());
  }
  
  auto inv_adjacent_vertices_after(index_type u) const {
    return boost::make_iterator_range(
        nodes[u].in_cmid(),
        std::end(nodes[u].in));
  }
  
  auto sorted_adjacent_vertices_after(index_type u) const {
    return boost::make_iterator_range(
        nodes[u].out_cmid(),
        std::end(nodes[u].out));
  }
  
  auto sorted_inv_adjacent_vertices_after(index_type u) const {
    return boost::make_iterator_range(
        nodes[u].in_cmid(),
        std::end(nodes[u].in));
  }
  
  void print() const {
    for (index_type u=0; u<n; ++u) {
      std::cout << u << ":" << std::endl;
      std::cout << "  out:";
      for (auto v : adjacent_vertices_before(u)) {
        std::cout << " " << v;
      }
      std::cout << " |";
      for (auto v : adjacent_vertices_after(u)) {
        std::cout << " " << v;
      }
      std::cout << std::endl;
      
      std::cout << "  in:";
      for (auto v : inv_adjacent_vertices_before(u)) {
        std::cout << " " << v;
      }
      std::cout << " |";
      for (auto v : inv_adjacent_vertices_after(u)) {
        std::cout << " " << v;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
};

#endif  // ORDERABLE_ADJACENCY_LISTMAT_WITH_RI_DEGREE_H_
