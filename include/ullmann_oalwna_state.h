#ifndef ULLMANN_OALWNA_STATE_H_
#define ULLMANN_OALWNA_STATE_H_

#include <iterator>

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
class ullmann_oalwna_state_base {
 protected:
  using IndexG = typename G::index_type;
  using IndexH = typename H::index_type;
   
  IndexG m;
  IndexH n;

  G const & g;
  H const & h;

  VertexEquivalencePredicate vertex_comp;
  EdgeEquivalencePredicate edge_comp;

  CompatibilityMatrix M;

  IndexOrderG const & index_order_g;
  typename IndexOrderG::const_iterator x_it;

 public:
  ullmann_oalwna_state_base(
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
        M{m, n},
        index_order_g{index_order_g},
        x_it{std::begin(index_order_g)} {
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
  
  ullmann_oalwna_state_base(ullmann_oalwna_state_base const &) = delete;

  bool empty() {
    return x_it == std::begin(index_order_g);
  }
  
  bool full() {
    return x_it == std::end(index_order_g);
  }

  auto candidates() {
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
class ullmann_oalwna_state_mono
  : public ullmann_oalwna_state_base<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        CompatibilityMatrix,
        IndexOrderG> {
 private:
  using base = ullmann_oalwna_state_base<
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
  using base::M;
  using base::x_it;
  
  using base::vertex_comp;
  using base::edge_comp;

  void filter(IndexG i, IndexH j) {
    for (IndexG ii=0; ii<m; ++ii) {
      M.unset(ii, j);
    }
    for (IndexH jj=0; jj<n; ++jj) {
      M.unset(i, jj);
    }
    M.set(i, j);
  }
  
  bool ullmann_condition(IndexG i, IndexH j) {
    auto const & i_adj_a = g.adjacent_vertices_after(i);
    auto const & j_adj = h.adjacent_vertices(j);
    auto const & i_inv_adj_a = g.inv_adjacent_vertices_after(i);
    auto const & j_inv_adj = h.inv_adjacent_vertices(j);
    return
        std::all_of(std::begin(i_adj_a), std::end(i_adj_a), [this,&j_adj,i,j](auto ii) {
          return std::any_of(std::begin(j_adj), std::end(j_adj), [this,i,j,ii](auto jj) {
            return M.get(ii, jj) && edge_comp(i, ii, j, jj);
          });
        })
        &&
        std::all_of(std::begin(i_inv_adj_a), std::end(i_inv_adj_a), [this,&j_inv_adj,i,j](auto ii) {
          return std::any_of(std::begin(j_inv_adj), std::end(j_inv_adj), [this,i,j,ii](auto jj) {
            return M.get(ii, jj) && edge_comp(ii, i, jj, j);
          });
        });
  }
  
  bool refine() {
    bool change;
    do {
      change = false;
      for (IndexG i=0; i<m; ++i) {
        for (IndexH j=0; j<n; ++j) {
          if (M.get(i, j) && !ullmann_condition(i, j)) {
            M.unset(i, j);
            if (!M.possible(i)) {
              return false;
            }
            change = true;
          }
        }
      }
    } while (change);
    return true;
  }

 public:
  ullmann_oalwna_state_mono(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : ullmann_oalwna_state_base<
            G,
            H,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            CompatibilityMatrix,
            IndexOrderG>(g, h, vertex_comp, edge_comp, index_order_g) {
    refine();
  }

  bool assign(IndexH y) {
    auto x = *x_it;
    filter(x, y);
    bool success = refine();
    return success;
  }
};

#endif  // ULLMANN_OALWNA_STATE_H_
