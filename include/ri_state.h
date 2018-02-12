#ifndef RI_STATE_H_
#define RI_STATE_H_

#include <utility>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <vector>
#include <stack>

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderG>
class ri_state_mono {
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
  
  using H_adjacent_vertices_type = typename std::decay_t<decltype(h.adjacent_vertices(IndexH()))>;
  
  H_adjacent_vertices_type h_vertices;
  
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
    /*
    for (auto i : g.adjacent_vertices_before(u)) {
      auto j = map[i];
      if (!h.edge(v, j)) {
        return false;
      }
    }
    for (auto i : g.inv_adjacent_vertices_before(u)) {
      auto j = map[i];
      if (!h.edge(j, v)) {
        return false;
      }
    }
    return true;
    */
  }

 public:
  ri_state_mono(
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
      
    /*for (auto i : index_order_g) {
      for (auto pi : index_order_g) {
        if (g_parents[pi].first != m && (g.edge(i, pi) || g.edge(pi, i))) {
          g_parents[i] = {pi, g.edge(pi, i)};
          break;
        }
      }
      if (g_parents[i].first == m) {
        g_parents[i] = {i, false};
      }
    }*/
  }
  
  ri_state_mono(ri_state_mono const &) = delete;

  bool empty() {
    return x_it == std::begin(index_order_g);
  }
  bool full() {
    return x_it == std::end(index_order_g);
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
  IndexH pop() {
    --x_it;
    
    auto x = *x_it;
    auto y = map[x];
    map[x] = n;
    inv[y] = m;
    return y;
  }
  
  H_adjacent_vertices_type const & candidates() {
    auto x = *x_it;
    auto parent = g_parents[x].first;
    auto out = g_parents[x].second;
    if (parent != x) {
      return out ? h.adjacent_vertices(map[parent]) : h.inv_adjacent_vertices(map[parent]);
    } else {
      return h_vertices;
    }
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
};

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderG>
class ri_state_ind
  : public ri_state_mono<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        IndexOrderG> {
 private:
  using base = ri_state_mono<
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
  
  std::vector<IndexG> g_out_count;
  std::vector<IndexG> g_in_count;
  
  std::vector<IndexH> h_out_count;
  std::vector<IndexH> h_in_count;
 
 public:
  ri_state_ind(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : ri_state_mono<
            G,
            H,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            IndexOrderG>(g, h, vertex_comp, edge_comp, index_order_g),
        g_out_count(m),
        g_in_count(m),
        h_out_count(n),
        h_in_count(n) {
    std::vector<IndexG> index_pos_g(m);
    for (IndexG i=0; i<m; ++i) {
      index_pos_g[index_order_g[i]] = i;
    }
    for (IndexG i=0; i<m; ++i) {
      auto i_pos = index_pos_g[i];
      auto const & i_adj = g.adjacent_vertices(i);
      g_out_count[i] = std::count_if(std::begin(i_adj), std::end(i_adj), [i_pos, &index_pos_g](auto ii) {
        return index_pos_g[ii] < i_pos;
      });
      auto const & i_inv_adj = g.inv_adjacent_vertices(i);
      g_in_count[i] = std::count_if(std::begin(i_inv_adj), std::end(i_inv_adj), [i_pos, &index_pos_g](auto ii) {
        return index_pos_g[ii] < i_pos;
      });
      //g_out_count[i] = g.out_degree_before(i);
      //g_in_count[i] = g.in_degree_before(i);
    }
  }
  
  void push(IndexH y) {
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
    for (auto j : h.adjacent_vertices(y)) {
      --h_in_count[j];
    }
    for (auto j : h.inv_adjacent_vertices(y)) {
      --h_out_count[j];
    }
  }
 
  bool assign(IndexH y) {
    auto x = *x_it;
    return
        g_out_count[x] == h_out_count[y] &&
        g_in_count[x] == h_in_count[y] &&
        base::assign(y);
  }
};

#endif  // RI_STATE_H
