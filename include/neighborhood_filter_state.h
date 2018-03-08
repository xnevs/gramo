#ifndef NEIGHBORHOOD_FILTER_STATE_H_
#define NEIGHBORHOOD_FILTER_STATE_H_

#include <iterator>
#include <vector>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/adaptor/filtered.hpp>

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename CompatibilityMatrix,
    typename IndexOrderG>
class neighborhood_filter_state_base {
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

  CompatibilityMatrix M;
  
  std::vector<IndexH> map;
  std::vector<IndexG> inv;

 public:
  neighborhood_filter_state_base(
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
        M(m, n),
        map(m, n),
        inv(n, m) {
    for (IndexG i=0; i<m; ++i) {
      for (IndexH j=0; j<n; ++j) {
        if (vertex_comp(i, j) &&
            g.out_degree(i) <= h.out_degree(j) &&
            g.in_degree(i) <= h.in_degree(j)) {
          M.set(i, j);
        }
      }
    }
  }
  
  neighborhood_filter_state_base(neighborhood_filter_state_base const &) = delete;

  bool empty() const {
    return x_it == std::begin(index_order_g);
  }
  
  bool full() const {
    return x_it == std::end(index_order_g);
  }

  auto candidates() const {
    auto x = *x_it;
    boost::counting_iterator<IndexH> begin{0}, end{n};
    return boost::adaptors::filter(
        boost::make_iterator_range(begin, end),
        [this, x](auto y){return M.get(x, y) && inv[y] == m/* && M.get(x, y)*/;});
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
    typename CompatibilityMatrix,
    typename IndexOrderG>
class neighborhood_filter_state_ind
  : public neighborhood_filter_state_base<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        CompatibilityMatrix,
        IndexOrderG> {
 private:
  using base = neighborhood_filter_state_base<
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
  using base::index_order_g;
  using base::x_it;
  using base::M;
  using base::inv;
  
  void neighborhood_filter_after(IndexG u, IndexH v) {
    for (auto i : g.adjacent_vertices_after(u)) {
      for (auto j : h.not_adjacent_vertices(v)) {
        if (inv[j] == m) {
          M.unset(i, j);
        }
      }
    }
    for (auto i : g.inv_adjacent_vertices_after(u)) {
      for (auto j : h.not_inv_adjacent_vertices(v)) {
        if (inv[j] == m) {
          M.unset(i, j);
        }
      }
    }
    for (auto j : h.adjacent_vertices(v)) {
      if (inv[j] == m) {
        for (auto i : g.not_adjacent_vertices_after(u)) {
          M.unset(i, j);
        }
      }
    }
    for (auto j : h.inv_adjacent_vertices(v)) {
      if (inv[j] == m) {
        for (auto i : g.not_inv_adjacent_vertices_after(u)) {
          M.unset(i, j);
        }
      }
    }
  }

 public:
  neighborhood_filter_state_ind(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : neighborhood_filter_state_base<
            G,
            H,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            CompatibilityMatrix,
            IndexOrderG>(g, h, vertex_comp, edge_comp, index_order_g) {
  }
  
  bool assign(IndexH y) {
    return true;
  }
  
  void push(IndexH y) {
    auto x = *x_it;
    neighborhood_filter_after(x, y);
    base::push(y);
  }
};

#endif  // NEIGHBORHOOD_FILTER_STATE_H
