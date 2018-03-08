#ifndef REFINED_RI_STATE_H_
#define REFINED_RI_STATE_H_

#include <iterator>
#include <utility>
#include <vector>
#include <algorithm>
#include <numeric>

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderG>
class refined_ri_state_mono {
 protected:
  using IndexG = typename G::index_type;
  using IndexH = typename H::index_type;
  
  IndexG m;
  IndexH n;
  
  G const & g;
  H const & h;

  VertexEquivalencePredicate vertex_comp;
  EdgeEquivalencePredicate edge_comp;

  IndexOrderG const & index_order_g;
  typename IndexOrderG::const_iterator x_it;

  std::vector<std::pair<IndexG,bool>> g_parents;

  std::vector<IndexH> map;
  std::vector<IndexG> inv;
  
  using H_adjacent_vertices_container_type = typename H::adjacent_vertices_container_type;
  
  H_adjacent_vertices_container_type h_vertices;
  
  bool topology_condition(IndexG u, IndexH v) {
    for (auto i : g.adjacent_vertices(u)) {
      auto j = map[i];
      if (j != n) {
        if (!h.edge(v, j) || !edge_comp(u, i, v, j)) {
          return false;
        }
      } else {
        bool exists = false;
        for (auto j : h.adjacent_vertices(v)) {
          if (inv[j] == m &&
              g.out_degree(i) <= h.out_degree(j) &&
              g.in_degree(i) <= h.in_degree(j)) {
            exists = true;
            break;
          }
        }
        if (!exists) {
          return false;
        }
      }
    }
    for (auto i : g.inv_adjacent_vertices(u)) {
      auto j = map[i];
      if (j != n) {
        if (!h.edge(j, v) || !edge_comp(i, u, j, v)) {
          return false;
        }
      } else {
        bool exists = false;
        for (auto j : h.inv_adjacent_vertices(v)) {
          if (inv[j] == m &&
              g.out_degree(i) <= h.out_degree(j) &&
              g.in_degree(i) <= h.in_degree(j)) {
            exists = true;
            break;
          }
        }
        if (!exists) {
          return false;
        }
      }
    }
    return true;
  }

 public:
  refined_ri_state_mono(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : m{g.num_vertices()},
        n{h.num_vertices()},
        g{g},
        h{h},
        vertex_comp{vertex_comp},
        edge_comp{edge_comp},
        index_order_g{index_order_g},
        x_it{std::begin(index_order_g)},
        g_parents(m,{m,true}),
        map(m, n),
        inv(n, m),
        h_vertices(n) {
    for (auto i : index_order_g) {
      auto const & i_adj = g.adjacent_vertices(i);
      auto const & i_inv_adj = g.inv_adjacent_vertices(i);
      auto parent_it = std::find_if(std::begin(i_adj), std::end(i_adj), [this](auto ii) {
        return g_parents[ii].first != m;
      });
      if (parent_it != std::end(i_adj)) {
        g_parents[i] = {*parent_it, false};
      } else {
        parent_it = std::find_if(std::begin(i_inv_adj), std::end(i_inv_adj), [this](auto ii) {
          return g_parents[ii].first != m;
        });
        if (parent_it != std::end(i_inv_adj)) {
          g_parents[i] = {*parent_it, true};
        } else {
          g_parents[i] = {i, false};
        }
      }
    }
    
    std::iota(std::begin(h_vertices), std::end(h_vertices), 0);
  }
  
  refined_ri_state_mono(refined_ri_state_mono const &) = delete;

  bool empty() {
    return x_it == std::begin(index_order_g);
  }
  
  bool full() {
    return x_it == std::end(index_order_g);
  }
  
  H_adjacent_vertices_container_type const & candidates() {
    auto x = *x_it;
    auto parent = g_parents[x].first;
    auto out = g_parents[x].second;
    if (parent != x) {
      return out ? h.adjacent_vertices(map[parent]) : h.inv_adjacent_vertices(map[parent]);
    } else {
      return h_vertices;
    }
  }

  void advance() {
  }
  
  void revert() {
  }

  bool assign(IndexH y) {
    auto x = *x_it;
    return
        inv[y] == m &&
        vertex_comp(x, y) &&
        g.out_degree(x) <= h.out_degree(y) &&
        g.in_degree(x) <= h.in_degree(y) &&
        topology_condition(x, y);
  }

  void push(IndexH y) {
    auto x = *x_it;
    
    map[x] = y;
    inv[y] = x;
    
    ++x_it;
  }
  
  IndexH pop() {
    --x_it;
    
    auto x = *x_it;
    auto y = map[x];
    map[x] = n;
    inv[y] = m;
    return y;
  }
};

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderG>
class refined_ri_state_ind
  : public refined_ri_state_mono<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        IndexOrderG> {
 private:
  using base = refined_ri_state_mono<
      G,
      H,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      IndexOrderG>;
      
 protected:
  using IndexG = typename base::IndexG;
  using IndexH = typename base::IndexH;

  using base::m;
  using base::n;
  using base::g;
  using base::h;
  using base::x_it;
  using base::map;
  using base::inv;
  
  using base::vertex_comp;
  using base::edge_comp;
  
  std::vector<IndexG> g_out_count;
  std::vector<IndexG> g_in_count;
  
  std::vector<IndexH> h_out_count;
  std::vector<IndexH> h_in_count;
 
 
  bool topology_condition(IndexG u, IndexH v) {
    for (auto i : g.adjacent_vertices(u)) {
      auto j = map[i];
      if (j != n) {
        if (!h.edge(v, j) || !edge_comp(u, i, v, j)) {
          return false;
        }
      } else {
        bool exists = false;
        for (auto j : h.adjacent_vertices(v)) {
          if (inv[j] == m &&
              g_out_count[i] == h_out_count[j] &&
              g_in_count[i] == h_in_count[j] &&
              g.out_degree(i) <= h.out_degree(j) &&
              g.in_degree(i) <= h.in_degree(j)) {
            exists = true;
            break;
          }
        }
        if (!exists) {
          return false;
        }
      }
    }
    for (auto i : g.inv_adjacent_vertices(u)) {
      auto j = map[i];
      if (j != n) {
        if (!h.edge(j, v) || !edge_comp(i, u, j, v)) {
          return false;
        }
      } else {
        bool exists = false;
        for (auto j : h.inv_adjacent_vertices(v)) {
          if (inv[j] == m &&
              g_out_count[i] == h_out_count[j] &&
              g_in_count[i] == h_in_count[j] &&
              g.out_degree(i) <= h.out_degree(j) &&
              g.in_degree(i) <= h.in_degree(j)) {
            exists = true;
            break;
          }
        }
        if (!exists) {
          return false;
        }
      }
    }
    return true;
  }
 
 public:
  refined_ri_state_ind(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : refined_ri_state_mono<
            G,
            H,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            IndexOrderG>(g, h, vertex_comp, edge_comp, index_order_g),
        g_out_count(m),
        g_in_count(m),
        h_out_count(n),
        h_in_count(n) {
  }
 
  bool assign(IndexH y) {
    auto x = *x_it;
    return
        inv[y] == m &&
        vertex_comp(x, y) &&
        g_out_count[x] == h_out_count[y] &&
        g_in_count[x] == h_in_count[y] &&
        g.out_degree(x) <= h.out_degree(y) &&
        g.in_degree(x) <= h.in_degree(y) &&
        topology_condition(x, y);
  }
  
  void push(IndexH y) {
    auto x = *x_it;
    for (auto i : g.adjacent_vertices(x)) {
      ++g_in_count[i];
    }
    for (auto i : g.inv_adjacent_vertices(x)) {
      ++g_out_count[i];
    }
    for (auto j : h.adjacent_vertices(y)) {
      ++h_in_count[j];
    }
    for (auto j : h.inv_adjacent_vertices(y)) {
      ++h_out_count[j];
    }
    base::push(y);
  }
  
  void pop() {
    auto y = base::pop();
    auto x = *x_it;
    for (auto i : g.adjacent_vertices(x)) {
      --g_in_count[i];
    }
    for (auto i : g.inv_adjacent_vertices(x)) {
      --g_out_count[i];
    }
    for (auto j : h.adjacent_vertices(y)) {
      --h_in_count[j];
    }
    for (auto j : h.inv_adjacent_vertices(y)) {
      --h_out_count[j];
    }
  }
};

#endif  // REFINED_RI_STATE_H
