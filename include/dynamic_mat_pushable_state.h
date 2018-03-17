#ifndef DYNAMIC_MAT_PUSHABLE_STATE_H_
#define DYNAMIC_MAT_PUSHABLE_STATE_H_

#include <iterator>
#include <vector>
#include <stack>

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename CompatibilityMatrix>
class dynamic_mat_pushable_state_base {
 protected:
  using IndexG = typename G::index_type;
  using IndexH = typename H::index_type;
   
  IndexG m;
  IndexH n;
  
  G & g;
  H const & h;

  VertexEquivalencePredicate vertex_comp;
  EdgeEquivalencePredicate edge_comp;
  
  std::vector<IndexH> map;
  std::vector<IndexG> inv;
  
  std::vector<IndexG> index_order_g;
  using x_it_type = typename decltype(index_order_g)::iterator;
  x_it_type x_it;

  CompatibilityMatrix M;
  
  using H_adjacent_vertices_container_type = typename H::adjacent_vertices_container_type;
  H_adjacent_vertices_container_type h_vertices;

  std::vector<std::pair<IndexH,bool>> h_parents;

 public:
  dynamic_mat_pushable_state_base(
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
        map(m, n),
        inv(n, m),
        index_order_g(m),
        x_it{std::begin(index_order_g)},
        M(m, n),
        h_vertices(n),
        h_parents(m) {
        
    std::iota(std::begin(index_order_g), std::end(index_order_g), 0);
        
    for (IndexG i=0; i<m; ++i) {
      for (IndexH j=0; j<n; ++j) {
        if (vertex_comp(i, j) &&
            g.out_degree(i) <= h.out_degree(j) &&
            g.in_degree(i) <= h.in_degree(j)) {
          M.set(i, j);
        }
      }
    }
    std::iota(h_vertices.begin(), h_vertices.end(), 0);
  }
  
  dynamic_mat_pushable_state_base(dynamic_mat_pushable_state_base const &) = delete;

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
          auto a_num = M.num_candidates(a);
          auto b_num = M.num_candidates(b);
          return
              a_num < b_num || (
              a_num == b_num && g.degree_neigh(a) > g.degree_neigh(b)) || (
              a_num == b_num && g.degree_neigh(a) == g.degree_neigh(b) && g.degree_unv(a) > g.degree_unv(b));
        });
        
    std::rotate(x_it, x_best_it, std::next(x_best_it));;
    
    auto x = *x_it;
    
    auto x_adj = g.adjacent_vertices(x);
    auto parent_it = std::find_if(std::begin(x_adj), std::end(x_adj), [this](auto u) {
      return map[u] != n;
    });
    if (parent_it != std::end(x_adj)) {
      h_parents[x] = {map[*parent_it], false};
    } else {
      auto x_inv_adj = g.inv_adjacent_vertices(x);
      parent_it = std::find_if(std::begin(x_inv_adj), std::end(x_inv_adj), [this](auto u) {
        return map[u] != n;
      });
      if (parent_it != std::end(x_inv_adj)) {
        h_parents[x] = {map[*parent_it], true};
      } else {
        h_parents[x] = {n, false};
      }
    }
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
    auto h_parent = h_parents[x].first;
    auto out = h_parents[x].second;
    if (h_parent != n) {
      return out ? h.adjacent_vertices(h_parent) : h.inv_adjacent_vertices(h_parent);
    } else {
      return h_vertices;
    }
  }

  void advance() {
    M.advance();
  }
  
  void revert() {
    M.revert();
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
    typename EdgeEquivalencePredicate,
    typename CompatibilityMatrix>
class dynamic_mat_pushable_state_ind
  : public dynamic_mat_pushable_state_base<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        CompatibilityMatrix> {
 private:
  using base = dynamic_mat_pushable_state_base<
      G,
      H,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      CompatibilityMatrix>;

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

  void filter_after(x_it_type u_it, IndexH v) {
    for (auto i_it=std::next(u_it); i_it!=std::end(index_order_g); ++i_it) {
      auto i = *i_it;
      M.unset(i, v);
    }
  }

  void neighborhood_filter_after(IndexG u, IndexH v) {
    for (auto i_it=std::next(x_it); i_it!=std::end(index_order_g); ++i_it) {
      auto i = *i_it;
      auto g_out = g.edge(u, i);
      if (g_out) {
        for (auto j : h.not_adjacent_vertices(v)) {
          if (inv[j] == m) {
            M.unset(i, j);
          }
        }
      } else {
        for (auto j : h.adjacent_vertices(v)) {
          if (inv[j] == m) {
            M.unset(i, j);
          }
        }
      }
      auto g_in = g.edge(i, u);
      if (g_in) {
        for (auto j : h.not_inv_adjacent_vertices(v)) {
          if (inv[j] == m) {
            M.unset(i, j);
          }
        }
      } else {
        for (auto j : h.inv_adjacent_vertices(v)) {
          if (inv[j] == m) {
            M.unset(i, j);
          }
        }
      }
    }
  }

 public:
  dynamic_mat_pushable_state_ind(
      G & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp)
      : dynamic_mat_pushable_state_base<
            G,
            H,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            CompatibilityMatrix>(g, h, vertex_comp, edge_comp) {
  }
  
  bool assign(IndexH y) {
    auto x = *x_it;
    return M.get(x, y);
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
    
    g.pop(x);
  }
};

#endif  // DYNAMIC_MAT_PUSHABLE_STATE_H
