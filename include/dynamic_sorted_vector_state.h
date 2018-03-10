#ifndef DYNAMIC_SORTED_VECTOR_STATE_H_
#define DYNAMIC_SORTED_VECTOR_STATE_H_

#include <iterator>
#include <vector>
#include <stack>

#include "sorted_vector.h"

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
class dynamic_sorted_vector_state_base {
 protected:
  using IndexG = typename G::index_type;
  using IndexH = typename H::index_type;
   
  IndexG m;
  IndexH n;
  
  G & g;
  H const & h;

  VertexEquivalencePredicate vertex_comp;
  EdgeEquivalencePredicate edge_comp;
  
  std::vector<IndexG> index_order_g;
  using x_it_type = typename decltype(index_order_g)::iterator;
  x_it_type x_it;
  
  std::vector<IndexH> map;
  std::vector<IndexG> inv;

  std::vector<sorted_vector<IndexH>> M;
  std::vector<std::stack<std::pair<IndexG,IndexH>>> changes;

 public:
  dynamic_sorted_vector_state_base(
      G & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp)
      : m{g.num_vertices()},
        n{h.num_vertices()},
        g{g},
        h{h},
        vertex_comp{vertex_comp},
        edge_comp{edge_comp},
        index_order_g(m),
        x_it{std::begin(index_order_g)},
        map(m, n),
        inv(n, m),
        M(m),
        changes(m) {
    std::iota(std::begin(index_order_g), std::end(index_order_g), 0);
    
    for (IndexG i=0; i<m; ++i) {
      for (IndexH j=0; j<n; ++j) {
        if (vertex_comp(i, j) &&
            g.out_degree(i) <= h.out_degree(j) &&
            g.in_degree(i) <= h.in_degree(j)) {
          M[i].insert(j);
        }
      }
    }
  }
  
  dynamic_sorted_vector_state_base(dynamic_sorted_vector_state_base const &) = delete;

  bool empty() const {
    return x_it == std::begin(index_order_g);
  }
  
  bool full() const {
    return x_it == std::end(index_order_g);
  }

  void prepare() {
    auto x_best_it = std::min_element(
        x_it,
        std::end(index_order_g),
        [this](auto const & a, auto const & b) {
          return M[a].size() < M[b].size();
        });
    std::rotate(x_it, x_best_it, std::next(x_best_it));
  }
  
  void forget() {
    auto it = std::lower_bound(
        std::next(x_it),
        std::end(index_order_g),
        *x_it);
    std::rotate(x_it, std::next(x_it), it);
  }
  
  auto const & candidates() const {
    auto x = *x_it;
    return M[x];
  }

  void advance() {
  }
  
  void revert() {
  }

  void push(IndexH y) {
    auto x = *x_it;
    map[x] = y;
    inv[y] = x;
    ++x_it;
  }
  
  void pop() {
    --x_it;
    auto x = *x_it;
    auto y = map[x];
    map[x] = n;
    inv[y] = m;
  }
};

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
class dynamic_sorted_vector_state_ind
  : public dynamic_sorted_vector_state_base<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate> {
 private:
  using base = dynamic_sorted_vector_state_base<
      G,
      H,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate>;

 protected:
  using typename base::IndexG;
  using typename base::IndexH;
  using typename base::x_it_type;

  using base::m;
  using base::n;
  using base::g;
  using base::h;
  using base::map;
  using base::inv;
  using base::index_order_g;
  using base::x_it;
  using base::M;
  using base::changes;

  void filter_after(x_it_type u_it, IndexH v) {
    auto u = *u_it;
    for (auto i_it=std::next(u_it); i_it!=std::end(index_order_g); ++i_it) {
      auto i = *i_it;
      auto v_pos = M[i].find(v);
      if (v_pos != M[i].end()) {
        M[i].erase(v_pos);
        changes[u].emplace(i, v);
      }
    }
  }

  void neighborhood_filter_after(IndexG u, IndexH v) {
    // TODO premakni v mono in extend
    for (auto i : g.adjacent_vertices_after(u)) {
      for (auto j : h.not_adjacent_vertices(v)) {
        if (inv[j] == m) {
          auto j_pos = M[i].find(j);
          if (j_pos != M[i].end()) {
            M[i].erase(j_pos);
            changes[u].emplace(i, j);
          }
        }
      }
    }
    for (auto i : g.inv_adjacent_vertices_after(u)) {
      for (auto j : h.not_inv_adjacent_vertices(v)) {
        if (inv[j] == m) {
          auto j_pos = M[i].find(j);
          if (j_pos != M[i].end()) {
            M[i].erase(j_pos);
            changes[u].emplace(i, j);
          }
        }
      }
    }
    
    {
      auto u_s_adj_a = g.sorted_adjacent_vertices_after(u);
      for (auto j : h.adjacent_vertices(v)) {
        if (inv[j] == m) {
          auto ni_it = std::begin(u_s_adj_a);
          auto ni_end = std::end(u_s_adj_a);
          for (auto i_it=std::next(x_it); i_it!=std::end(index_order_g); ++i_it) {
            auto i = *i_it;
            if (ni_it == ni_end || i < *ni_it ) {
              auto j_pos = M[i].find(j);
              if (j_pos != M[i].end()) {
                M[i].erase(j_pos);
                changes[u].emplace(i, j);
              }
            } else {
              ++ni_it;
            }
          }
        }
      }
    }
    
    {
      auto u_s_inv_adj_a = g.sorted_inv_adjacent_vertices_after(u);
      for (auto j : h.inv_adjacent_vertices(v)) {
        if (inv[j] == m) {
          auto ni_it = std::begin(u_s_inv_adj_a);
          auto ni_end = std::end(u_s_inv_adj_a);
          for (auto i_it=std::next(x_it); i_it!=std::end(index_order_g); ++i_it) {
            auto i = *i_it;
            if (ni_it == ni_end || i < *ni_it ) {
              auto j_pos = M[i].find(j);
              if (j_pos != M[i].end()) {
                M[i].erase(j_pos);
                changes[u].emplace(i, j);
              }
            } else {
              ++ni_it;
            }
          }
        }
      }
    } 
  }

 public:
  dynamic_sorted_vector_state_ind(
      G & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp)
      : dynamic_sorted_vector_state_base<
            G,
            H,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate>(g, h, vertex_comp, edge_comp) {
  }
  
  bool assign(IndexH y) {
    return true;
  }
  
  void push(IndexH y) {
    auto x = *x_it;
    
    g.push(x);
    
    filter_after(x_it, y);
    neighborhood_filter_after(x, y);
    base::push(y);
  }
  
  void pop() {
    base::pop();
    
    auto x = *x_it;
    while(!changes[x].empty()) {
      auto p = changes[x].top();
      changes[x].pop();
      // invariant: p.first comes after x in index_order_g (push() makes sure of that),
      // thus this doesn't invalidate any iterators held at previous levels of explore
      M[p.first].insert(p.second);
    }
    
    g.pop(x);
  }
};

#endif  // DYNAMIC_SORTED_VECTOR_STATE_H
