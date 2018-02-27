#ifndef ULLIMP_RI_STATE_H_
#define ULLIMP_RI_STATE_H_

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
class ullimp_ri_state_base {
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
  
  IndexG level;
  IndexG const cutoff;

 public:
  ullimp_ri_state_base(
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
        inv(n, m),
        level{0},
        cutoff{static_cast<IndexG>(m)} {
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
  
  ullimp_ri_state_base(ullimp_ri_state_base const &) = delete;

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
        [this, x](auto y){return/* M.get(x, y) &&*/ inv[y] == m && M.get(x, y);});
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
    ++level;
  }
  
  void pop() {
    --level;
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
class ullimp_ri_state_ind
  : public ullimp_ri_state_base<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        CompatibilityMatrix,
        IndexOrderG> {
 private:
  using base = ullimp_ri_state_base<
      G,
      H,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      CompatibilityMatrix,
      IndexOrderG>;

 protected:
  using IndexG = typename base::IndexG;
  using IndexH = typename base::IndexH;
  
  using base::vertex_comp;
  using base::edge_comp;

  using base::m;
  using base::n;
  using base::g;
  using base::h;
  using base::index_order_g;
  using base::x_it;
  using base::M;
  using base::map;
  using base::inv;
  using base::level;
  using base::cutoff;
  
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
  ullimp_ri_state_ind(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : ullimp_ri_state_base<
            G,
            H,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            CompatibilityMatrix,
            IndexOrderG>(g, h, vertex_comp, edge_comp, index_order_g) {
  }
  
  bool assign(IndexH y) {
    if (level <= cutoff) {
      return true;
    } else {
      auto x = *x_it;
      if (topology_condition(x, y)) {
        auto const & y_adj = h.adjacent_vertices(y);
        auto h_out_degree_before_y = std::count_if(std::begin(y_adj), std::end(y_adj), [this](auto j) {
          return inv[j] != m;
        });
        auto const & y_inv_adj = h.inv_adjacent_vertices(y);
        auto h_in_degree_before_y = std::count_if(std::begin(y_inv_adj), std::end(y_inv_adj), [this](auto j) {
          return inv[j] != m;
        });
        return
            g.out_degree_before(x) == h_out_degree_before_y &&
            g.in_degree_before(x) == h_in_degree_before_y;
      } else {
        return false;
      }
    }
  }
  
  void push(IndexH y) {
    auto x = *x_it;
    if (level < cutoff) {
      neighborhood_filter_after(x, y);
    }
    base::push(y);
  }
};

#endif  // ULLIMP_RI_STATE_H
