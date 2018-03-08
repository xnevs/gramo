#ifndef RI_LOOKAHEAD_STATE_H_
#define RI_LOOKAHEAD_STATE_H_

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
class ri_lookahead_state_mono {
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
  
  template<typename Index>
  struct rank {
    Index vis_out_degree;
    Index vis_in_degree;
    Index neigh_out_degree;
    Index neigh_in_degree;
    Index unv_out_degree;
    Index unv_in_degree;
    
    rank(Index vod, Index vid, Index nod, Index nid, Index uod, Index uid)
        : vis_out_degree{vod},
          vis_in_degree{vid},
          neigh_out_degree{nod},
          neigh_in_degree{nid},
          unv_out_degree{uod},
          unv_in_degree{uid} {
    }
    
    rank()
        : vis_out_degree{0},
          vis_in_degree{0},
          neigh_out_degree{0},
          neigh_in_degree{0},
          unv_out_degree{0},
          unv_in_degree{0} {
    }
  };
  
  std::vector<rank<IndexG>> g_ranks;
  std::vector<rank<IndexH>> h_ranks;
  
  enum struct Flag {
    vis,
    neigh,
    unv
  };
  
  bool topology_condition(IndexG u, IndexH v) {
    for (auto i : g.adjacent_vertices(u)) {
      auto j = map[i];
      if (j != n) {
        if (!h.edge(v, j) || !edge_comp(u, i, v, j)) {
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
      }
    }
    return true;
  }

 public:
  ri_lookahead_state_mono(
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
        g_ranks(m),
        h_ranks(n) {
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
    
    for (auto u : index_order_g) {
      g_ranks[u].unv_out_degree = g.out_degree(u);
      g_ranks[u].unv_in_degree = g.in_degree(u);
    }
    
    std::vector<Flag> g_flags(m, Flag::unv);
    for(auto u : index_order_g) {
      if (g_flags[u] == Flag::unv) {
        for (auto i : g.inv_adjacent_vertices(u)) {
          if (g_flags[i] != Flag::vis) {
            --g_ranks[i].unv_out_degree;
          }
        }
        for (auto i : g.adjacent_vertices(u)) {
          if (g_flags[i] != Flag::vis) {
            --g_ranks[i].unv_in_degree;
          }
        }
      } else if (g_flags[u] == Flag::neigh) {
        for (auto i : g.inv_adjacent_vertices(u)) {
          if (g_flags[i] != Flag::vis) {
            --g_ranks[i].neigh_out_degree;
          }
        }
        for (auto i : g.adjacent_vertices(u)) {
          if (g_flags[i] != Flag::vis) {
            --g_ranks[i].neigh_in_degree;
          }
        }
      }
      g_flags[u] = Flag::vis;
      for (auto i : g.inv_adjacent_vertices(u)) {
        if (g_flags[i] != Flag::vis) {
          ++g_ranks[i].vis_out_degree;
          if (g_flags[i] == Flag::unv) {
            g_flags[i] = Flag::neigh;
            for (auto k : g.inv_adjacent_vertices(i)) {
              if (g_flags[k] != Flag::vis) {
                ++g_ranks[k].neigh_out_degree;
              }
            }
            for (auto k : g.adjacent_vertices(i)) {
              if (g_flags[k] != Flag::vis) {
                ++g_ranks[k].neigh_in_degree;
              }
            }
          }
        }
      }
      for (auto i : g.adjacent_vertices(u)) {
        if (g_flags[i] != Flag::vis) {
          ++g_ranks[i].vis_in_degree;
          if (g_flags[i] == Flag::unv) {
            g_flags[i] = Flag::neigh;
            for (auto k : g.inv_adjacent_vertices(i)) {
              if (g_flags[k] != Flag::vis) {
                ++g_ranks[k].neigh_out_degree;
              }
            }
            for (auto k : g.adjacent_vertices(i)) {
              if (g_flags[k] != Flag::vis) {
                ++g_ranks[k].neigh_in_degree;
              }
            }
          }
        }
      }
    }
    
    for (IndexH v=0; v<n; ++v) {
      h_ranks[v].unv_out_degree = h.out_degree(v);
      h_ranks[v].unv_in_degree = h.in_degree(v);
    }
  }
  
  ri_lookahead_state_mono(ri_lookahead_state_mono const &) = delete;

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
        g_ranks[x].vis_out_degree <= h_ranks[y].vis_out_degree &&
        g_ranks[x].vis_in_degree <= h_ranks[y].vis_in_degree &&
        g_ranks[x].neigh_out_degree <= h_ranks[y].neigh_out_degree &&
        g_ranks[x].neigh_in_degree <= h_ranks[y].neigh_in_degree &&
        g_ranks[x].unv_out_degree <= h_ranks[y].unv_out_degree &&
        g_ranks[x].unv_in_degree <= h_ranks[y].unv_in_degree &&
        topology_condition(x, y);
  }

  void push(IndexH y) {
    auto x = *x_it;
    
    map[x] = y;
    inv[y] = x;
    
    auto y_flag = h_ranks[y].vis_out_degree > 0 || h_ranks[y].vis_in_degree > 0 ? Flag::neigh : Flag::unv;
    
    if (y_flag == Flag::unv) {
      for (auto v : h.inv_adjacent_vertices(y)) {
        if (inv[v] == m) {
          --h_ranks[v].unv_out_degree;
        }
      }
      for (auto v : h.adjacent_vertices(y)) {
        if (inv[v] == m) {
          --h_ranks[v].unv_in_degree;
        }
      }
    } else if (y_flag == Flag::neigh) {
      for (auto v : h.inv_adjacent_vertices(y)) {
        if (inv[v] == m) {
          --h_ranks[v].neigh_out_degree;
        }
      }
      for (auto v : h.adjacent_vertices(y)) {
        if (inv[v] == m) {
          --h_ranks[v].neigh_in_degree;
        }
      }
    }
    for (auto v : h.inv_adjacent_vertices(y)) {
      if (inv[v] == m) {
        auto v_flag = h_ranks[v].vis_out_degree > 0 || h_ranks[v].vis_in_degree > 0 ? Flag::neigh : Flag::unv;
        if (v_flag == Flag::unv) {
          for (auto w : h.inv_adjacent_vertices(v)) {
            if (inv[w] == m) {
              ++h_ranks[w].neigh_out_degree;
            }
          }
          for (auto w : h.adjacent_vertices(v)) {
            if (inv[w] == m) {
              ++h_ranks[w].neigh_in_degree;
            }
          }
        }
        ++h_ranks[v].vis_out_degree;
      }
    }
    for (auto v : h.adjacent_vertices(y)) {
      if (inv[v] == m) {
        auto v_flag = h_ranks[v].vis_out_degree > 0 || h_ranks[v].vis_in_degree > 0 ? Flag::neigh : Flag::unv;
        if (v_flag == Flag::unv) {
          for (auto w : h.inv_adjacent_vertices(v)) {
            if (inv[w] == m) {
              ++h_ranks[w].neigh_out_degree;
            }
          }
          for (auto w : h.adjacent_vertices(v)) {
            if (inv[w] == m) {
              ++h_ranks[w].neigh_in_degree;
            }
          }
        }
        ++h_ranks[v].vis_in_degree;
      }
    }
    
    ++x_it;
  }
  
  IndexH pop() {
    --x_it;
    
    auto x = *x_it;
    auto y = map[x];
    
    auto y_flag = h_ranks[y].vis_out_degree > 0 || h_ranks[y].vis_in_degree > 0 ? Flag::neigh : Flag::unv;
    if (y_flag == Flag::unv) {
      for (auto v : h.inv_adjacent_vertices(y)) {
        if (inv[v] == m) {
          ++h_ranks[v].unv_out_degree;
        }
      }
      for (auto v : h.adjacent_vertices(y)) {
        if (inv[v] == m) {
          ++h_ranks[v].unv_in_degree;
        }
      }
    } else if (y_flag == Flag::neigh) {
      for (auto v : h.inv_adjacent_vertices(y)) {
        if (inv[v] == m) {
          ++h_ranks[v].neigh_out_degree;
        }
      }
      for (auto v : h.adjacent_vertices(y)) {
        if (inv[v] == m) {
          ++h_ranks[v].neigh_in_degree;
        }
      }
    }
    for (auto v : h.inv_adjacent_vertices(y)) {
      if (inv[v] == m) {
        --h_ranks[v].vis_out_degree;
        auto v_flag = h_ranks[v].vis_out_degree > 0 || h_ranks[v].vis_in_degree > 0 ? Flag::neigh : Flag::unv;
        if (v_flag == Flag::unv) {
          for (auto w : h.inv_adjacent_vertices(v)) {
            if (inv[w] == m) {
              --h_ranks[w].neigh_out_degree;
            }
          }
          for (auto w : h.adjacent_vertices(v)) {
            if (inv[w] == m) {
              --h_ranks[w].neigh_in_degree;
            }
          }
        }
      }
    }
    for (auto v : h.adjacent_vertices(y)) {
      if (inv[v] == m) {
        --h_ranks[v].vis_in_degree;
        auto v_flag = h_ranks[v].vis_out_degree > 0 || h_ranks[v].vis_in_degree > 0 ? Flag::neigh : Flag::unv;
        if (v_flag == Flag::unv) {
          for (auto w : h.inv_adjacent_vertices(v)) {
            if (inv[w] == m) {
              --h_ranks[w].neigh_out_degree;
            }
          }
          for (auto w : h.adjacent_vertices(v)) {
            if (inv[w] == m) {
              --h_ranks[w].neigh_in_degree;
            }
          }
        }
      }
    }
    
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
class ri_lookahead_state_ind
  : public ri_lookahead_state_mono<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        IndexOrderG> {
 private:
  using base = ri_lookahead_state_mono<
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
  using base::g_ranks;
  using base::h_ranks;
  using base::vertex_comp;
  
  using base::topology_condition;
 
 public:
  using base::ri_lookahead_state_mono;
 
  bool assign(IndexH y) {
    auto x = *x_it;
    return
        inv[y] == m &&
        vertex_comp(x, y) &&
        g_ranks[x].vis_out_degree == h_ranks[y].vis_out_degree &&
        g_ranks[x].vis_in_degree == h_ranks[y].vis_in_degree &&
        g_ranks[x].neigh_out_degree <= h_ranks[y].neigh_out_degree &&
        g_ranks[x].neigh_in_degree <= h_ranks[y].neigh_in_degree &&
        g_ranks[x].unv_out_degree <= h_ranks[y].unv_out_degree &&
        g_ranks[x].unv_in_degree <= h_ranks[y].unv_in_degree &&
        topology_condition(x, y);
  }
};

#endif  // RI_LOOKAHEAD_STATE_H
