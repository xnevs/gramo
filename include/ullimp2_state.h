#ifndef ULLIMP2_STATE_H_
#define ULLIMP2_STATE_H_

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
class ullimp2_state_base {
 protected:
  using IndexG = typename G::index_type;
  using IndexH = typename H::index_type;
  
  using x_it_type = typename IndexOrderG::const_iterator;
   
  IndexG m;
  IndexH n;
  
  G const & g;
  H const & h;

  VertexEquivalencePredicate vertex_comp;
  EdgeEquivalencePredicate edge_comp;

  IndexOrderG const & index_order_g;
  x_it_type x_it;

  CompatibilityMatrix M;

 public:
  ullimp2_state_base(
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
        M(m, n) {
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
  
  ullimp2_state_base(ullimp2_state_base const &) = delete;

  bool empty() const {
    return x_it == std::begin(index_order_g);
  }
  
  bool full() const {
    return x_it == std::end(index_order_g);
  }
  
  void prepare() {
  }
  
  void forget() {
  }

  auto candidates() const {
    auto x = *x_it;
    boost::counting_iterator<IndexH> begin{0}, end{n};
    return boost::adaptors::filter(
        boost::make_iterator_range(begin, end),
        [this,x](auto y){return M.get(x, y);});
  }

  void advance() {
    M.advance();
  }
  
  void revert() {
    M.revert();
  }

  void push(IndexH y) {
    ++x_it;
  }
  
  void pop() {
    --x_it;
  }
};

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename CompatibilityMatrix,
    typename IndexOrderG>
class ullimp2_state_ind
  : public ullimp2_state_base<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        CompatibilityMatrix,
        IndexOrderG> {
 private:
  using base = ullimp2_state_base<
      G,
      H,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      CompatibilityMatrix,
      IndexOrderG>;

 protected:
  using IndexG = typename base::IndexG;
  using IndexH = typename base::IndexH;
  using x_it_type = typename base::x_it_type;

  using base::m;
  using base::n;
  using base::g;
  using base::h;
  using base::index_order_g;
  using base::x_it;
  using base::M;

  void neighborhood_filter_after(x_it_type u_it, IndexH v) {
    auto u = *u_it;
    for (auto i_it=std::next(u_it); i_it!=std::end(index_order_g); ++i_it) {
      auto i = *i_it;
      M.unset(i, v);
      bool out_g = g.edge(u, i);
      bool in_g = g.edge(i, u);
      for (IndexH j=0; j<n; ++j) {
        if (M.get(i, j)) {
          bool out_h = h.edge(v, j);
          bool in_h = h.edge(j, v);
          if (out_g != out_h || in_g != in_h) {
            M.unset(i, j);
          }
        }
      }
    }
  }

 public:
  ullimp2_state_ind(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : ullimp2_state_base<
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
    neighborhood_filter_after(x_it, y);
    base::push(y);
  }
};

#endif  // ULLIMP2_STATE_H
