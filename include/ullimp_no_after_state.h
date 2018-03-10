#ifndef ULLIMP_NO_AFTER_STATE_H_
#define ULLIMP_NO_AFTER_STATE_H_

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
    typename CompatibilityMatrix,
    typename IndexOrderG>
class ullimp_no_after_state_mono {
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

  CompatibilityMatrix M;
  
  std::vector<IndexG> index_pos_g;
  
  void neighborhood_filter_after(IndexG u, IndexH v) {
    for (auto i : g.adjacent_vertices(u)) {
      if (map[i] == n) {
        for (auto j : h.not_adjacent_vertices(v)) {
          if (inv[j] == m) {
            M.unset(i, j);
          }
        }
      }
    }
    for (auto i : g.inv_adjacent_vertices(u)) {
      if (map[i] ==  n) {
        for (auto j : h.not_inv_adjacent_vertices(v)) {
          if (inv[j] == m) {
            M.unset(i, j);
          }
        }
      }
    }
  }
  
  bool ullmann_condition(IndexG u, IndexH v) {
    for (auto i : g.adjacent_vertices(u)) {
      bool all_false = true;
      for (auto j : h.adjacent_vertices(v)) {
        if (M.get(i, j)) {
          all_false = false;
          break;
        }
      }
      if (all_false) {
        return false;
      }
    }
    for (auto i : g.inv_adjacent_vertices(u)) {
      bool all_false = true;
      for (auto j : h.inv_adjacent_vertices(v)) {
        if (M.get(i, j)) {
          all_false = false;
          break;
        }
      }
      if (all_false) {
        return false;
      }
    }
    return true;
  }
  
  void refine() {
    bool changed;
    do {
      changed = false;
      for (IndexG u=0; u<m; ++u) {
        for (IndexH v=0; v<n; ++v) {
          if (M.get(u, v) && !ullmann_condition(u, v)) {
            M.unset(u, v);
            changed = true;
          }
        }
      }
    } while(changed);
  }
  
  bool partial_ullmann_condition(IndexG u, IndexH v) {
    for (auto i : g.adjacent_vertices(u)) {
      if (map[i] == n) {
        bool all_false = true;
        for (auto j : h.adjacent_vertices(v)) {
          if (M.get(i, j) && inv[j] == m) {
            all_false = false;
            break;
          }
        }
        if (all_false) {
          return false;
        }
      }
    }
    for (auto i : g.inv_adjacent_vertices(u)) {
      if (map[i] == n) {
        bool all_false = true;
        for (auto j : h.inv_adjacent_vertices(v)) {
          if (M.get(i, j) && inv[j] == m) {
            all_false = false;
            break;
          }
        }
        if (all_false) {
          return false;
        }
      }
    }
    return true;
  }
  
  void partial_refine(IndexG u, IndexH v) {
    for (auto i : g.adjacent_vertices(u)) {
      if (map[i] == n) {
        for (auto j : h.adjacent_vertices(v)) {
          if (M.get(i, j) && inv[j] == m && !partial_ullmann_condition(i, j)) {
            M.unset(i, j);
          }
        }
      }
    }
    for (auto i : g.inv_adjacent_vertices(u)) {
      if (map[i] == n) {
        for (auto j : h.inv_adjacent_vertices(v)) {
          if (M.get(i, j) && inv[j] == m && !partial_ullmann_condition(i, j)) {
            M.unset(i, j);
          }
        }
      }
    }
  }

 public:
  ullimp_no_after_state_mono(
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
        h_vertices(n),
        M(m, n),
        index_pos_g(m) {
    for (IndexG i=0; i<m; ++i) {
      index_pos_g[index_order_g[i]] = i;
    }
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
    
    for (IndexG i=0; i<m; ++i) {
      for (IndexH j=0; j<n; ++j) {
        if (vertex_comp(i, j) &&
            g.out_degree(i) <= h.out_degree(j) &&
            g.in_degree(i) <= h.in_degree(j)) {
          M.set(i, j);
        }
      }
    }
    
    refine();
  }
  
  ullimp_no_after_state_mono(ullimp_no_after_state_mono const &) = delete;

  bool empty() {
    return x_it == std::begin(index_order_g);
  }
  
  bool full() {
    return x_it == std::end(index_order_g);
  }
  
  void prepare() {
  }
  
  void forget() {
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
        M.get(x, y) &&
        topology_condition(x, y);
  }

  void push(IndexH y) {
    auto x = *x_it;
    
    map[x] = y;
    inv[y] = x;
    
    M.advance();
    neighborhood_filter_after(x, y);
    if (std::distance(std::begin(index_order_g), x_it) < m/2) {
      partial_refine(x, y);
    }
    
    ++x_it;
  }
  
  IndexH pop() {
    --x_it;
    
    auto x = *x_it;
    auto y = map[x];
    
    M.revert();
    
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
    typename CompatibilityMatrix,
    typename IndexOrderG>
class ullimp_no_after_state_ind
  : public ullimp_no_after_state_mono<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        CompatibilityMatrix,
        IndexOrderG> {
 private:
  using base = ullimp_no_after_state_mono<
      G,
      H,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      CompatibilityMatrix,
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
  using base::M;
  using base::index_pos_g;
  using base::vertex_comp;
  
  using base::partial_refine;
  
  std::vector<IndexG> g_out_count;
  std::vector<IndexG> g_in_count;
  
  std::vector<IndexH> h_out_count;
  std::vector<IndexH> h_in_count;
 
  void neighborhood_filter_after(IndexG u, IndexH v) {
    base::neighborhood_filter_after(u, v);
    
    for (auto j : h.adjacent_vertices(v)) {
      if (inv[j] == m) {
        for (auto i : g.not_adjacent_vertices(u)) {
          if (map[i] == n) {
            M.unset(i, j);
          }
        }
      }
    }
    for (auto j : h.inv_adjacent_vertices(v)) {
      if (inv[j] == m) {
        for (auto i : g.not_inv_adjacent_vertices(u)) {
          if (map[i] == n) {
            M.unset(i, j);
          }
        }
      }
    }
  }
 
 public:
  ullimp_no_after_state_ind(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : ullimp_no_after_state_mono<
            G,
            H,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            CompatibilityMatrix,
            IndexOrderG>(g, h, vertex_comp, edge_comp, index_order_g),
        g_out_count(m),
        g_in_count(m),
        h_out_count(n),
        h_in_count(n) {
    for (IndexG i=0; i<m; ++i) {
      auto i_pos = index_pos_g[i];
      auto const & i_adj = g.adjacent_vertices(i);
      g_out_count[i] = std::count_if(std::begin(i_adj), std::end(i_adj), [this, i_pos](auto ii) {
        return index_pos_g[ii] < i_pos;
      });
      auto const & i_inv_adj = g.inv_adjacent_vertices(i);
      g_in_count[i] = std::count_if(std::begin(i_inv_adj), std::end(i_inv_adj), [this, i_pos](auto ii) {
        return index_pos_g[ii] < i_pos;
      });
    }
  }
 
  bool assign(IndexH y) {
    auto x = *x_it;
    return
        M.get(x, y) &&
        inv[y] == m &&
        vertex_comp(x, y);
  }
  
  void push(IndexH y) {
    for (auto j : h.adjacent_vertices(y)) {
      ++h_in_count[j];
    }
    for (auto j : h.inv_adjacent_vertices(y)) {
      ++h_out_count[j];
    }
    
    auto x = *x_it;
    
    map[x] = y;
    inv[y] = x;
    
    M.advance();
    neighborhood_filter_after(x, y);
    partial_refine(x, y);
    
    ++x_it;
  }
  
  void pop() {
    --x_it;
    
    auto x = *x_it;
    auto y = map[x];
    
    M.revert();
    
    map[x] = n;
    inv[y] = m;
    
    for (auto j : h.adjacent_vertices(y)) {
      --h_in_count[j];
    }
    for (auto j : h.inv_adjacent_vertices(y)) {
      --h_out_count[j];
    }
  }
};

#endif  // ULLIMP_NO_AFTER_STATE_H
