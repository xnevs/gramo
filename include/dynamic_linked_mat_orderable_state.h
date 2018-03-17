#ifndef DYNAMIC_LINKED_MAT_ORDERABLE_STATE_H_
#define DYNAMIC_LINKED_MAT_ORDERABLE_STATE_H_

#include <iterator>
#include <vector>
#include <stack>

#include <boost/range/iterator_range.hpp>

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename CompatibilityMatrix>
class dynamic_linked_mat_orderable_state_base {
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

 public:
  dynamic_linked_mat_orderable_state_base(
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
        M(m, n) {
        
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
    M.init();
  }
  
  dynamic_linked_mat_orderable_state_base(dynamic_linked_mat_orderable_state_base const &) = delete;

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
          return M.num_candidates(a) < M.num_candidates(b);
        });
        
    std::rotate(x_it, x_best_it, std::next(x_best_it));;
  }
  
  void forget() {
    auto it = std::lower_bound(
        std::next(x_it),
        std::end(index_order_g),
        *x_it);
    std::rotate(x_it, std::next(x_it), it);
  }
  
  auto candidates() {
    auto x = *x_it;
    
    return boost::make_iterator_range(
        M.row_begin(x),
        M.row_end(x));
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
class dynamic_linked_mat_orderable_state_ind
  : public dynamic_linked_mat_orderable_state_base<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        CompatibilityMatrix> {
 private:
  using base = dynamic_linked_mat_orderable_state_base<
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
    // TODO premakni v mono in extend
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
      /*for (auto j_it=M.row_begin(i); j_it!=M.row_end(i); ++j_it) {
        auto j = *j_it;
        auto h_out = h.edge(v, j);
        auto h_in = h.edge(j, v);
        if (g_out != h_out || g_in != h_in) {
          M.unset(i, j);
        }
      }*/
    }
  }

 public:
  dynamic_linked_mat_orderable_state_ind(
      G & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp)
      : dynamic_linked_mat_orderable_state_base<
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
    
    filter_after(x_it, y);
    neighborhood_filter_after(x, y);
    //if (std::distance(std::begin(index_order_g), x_it) < m/2) {
    //partial_refine(x, y);
    //}
    base::push(y);
  }
};

#endif  // DYNAMIC_LINKED_MAT_ORDERABLE_STATE_H
