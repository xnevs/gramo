#ifndef RI_DYNAMIC_PARENT_STATE_H_
#define RI_DYNAMIC_PARENT_STATE_H_

#include <iterator>
#include <utility>
#include <vector>
#include <stack>
#include <algorithm>
#include <numeric>

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderG>
class ri_dynamic_parent_state_mono {
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

  std::stack<std::pair<IndexH,bool>> h_parents;

  std::vector<IndexH> map;
  std::vector<IndexG> inv;
  
  using H_adjacent_vertices_container_type = typename H::adjacent_vertices_container_type;
  
  H_adjacent_vertices_container_type h_vertices;
  
  bool topology_condition(IndexG u, IndexH v) {
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
  }

 public:
  ri_dynamic_parent_state_mono(
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
        map(m, n),
        inv(n, m),
        h_vertices(n) {
    h_parents.emplace(n, false);
    std::iota(std::begin(h_vertices), std::end(h_vertices), 0);
  }
  
  ri_dynamic_parent_state_mono(ri_dynamic_parent_state_mono const &) = delete;

  bool empty() {
    return x_it == std::begin(index_order_g);
  }
  
  bool full() {
    return x_it == std::end(index_order_g);
  }
  
  H_adjacent_vertices_container_type const & candidates() {
    auto parent = h_parents.top().first;
    auto out = h_parents.top().second;
    if (parent != n) {
      return out ? h.adjacent_vertices(parent) : h.inv_adjacent_vertices(parent);
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
    
    if (!full()) {
      x = *x_it;
      IndexH parent = n;
      bool out = false;
      IndexH size = n;
      for (auto u : g.adjacent_vertices_before(x)) {
        auto v = map[u];
        auto v_in_degree = h.in_degree(v);
        if (v_in_degree < size) {
          parent = v;
          out = false;
          size = v_in_degree;
        }
      }
      for (auto u : g.inv_adjacent_vertices_before(x)) {
        auto v = map[u];
        auto v_out_degree = h.out_degree(v);
        if (v_out_degree < size) {
          parent = v;
          out = true;
          size = v_out_degree;
        }
      }
      h_parents.emplace(parent, out);
    }
  }
  
  IndexH pop() {
    if (!full()) {
      h_parents.pop();
    }
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
class ri_dynamic_parent_state_ind
  : public ri_dynamic_parent_state_mono<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        IndexOrderG> {
 private:
  using base = ri_dynamic_parent_state_mono<
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
  ri_dynamic_parent_state_ind(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : ri_dynamic_parent_state_mono<
            G,
            H,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            IndexOrderG>(g, h, vertex_comp, edge_comp, index_order_g),
        g_out_count(m),
        g_in_count(m),
        h_out_count(n),
        h_in_count(n) {
    for (IndexG i=0; i<m; ++i) {
      g_out_count[i] = g.out_degree_before(i);
      g_in_count[i] = g.in_degree_before(i);
    }
  }
 
  bool assign(IndexH y) {
    auto x = *x_it;
    return
        g_out_count[x] == h_out_count[y] &&
        g_in_count[x] == h_in_count[y] &&
        base::assign(y);
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
};

#endif  // RI_DYNAMIC_PARENT_STATE_H
