#ifndef ULLMANN_STATE_H_
#define ULLMANN_STATE_H_

#include <iterator>
#include <numeric>
#include <stack>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/adaptor/filtered.hpp>

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    template <typename, typename> typename CompatibilityMatrix,
    typename IndexOrderG>
class ullmann_state_base {
 protected:
  using IndexG = typename G::index_type;
  using IndexH = typename H::index_type;
   
  IndexG m;
  IndexH n;
  
  G const & g;
  H const & h;

  VertexEquivalencePredicate vertex_comp;
  EdgeEquivalencePredicate edge_comp;

  CompatibilityMatrix<IndexG, IndexH> M;

  IndexOrderG const & index_order_g;
  typename IndexOrderG::const_iterator x_it;

 public:
  // TODO delete copy constructor
  ullmann_state_base(
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
  
  ullmann_state_base(ullmann_state_base const &) = delete;

  bool empty() {
    return x_it == std::begin(index_order_g);
  }
  bool full() {
    return x_it == std::end(index_order_g);
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


  auto candidates() {
    auto x = *x_it;
    boost::counting_iterator<IndexH> begin{0}, end{n};
    return boost::adaptors::filter(
        boost::make_iterator_range(begin, end),
        [this,x](auto y){return M.get(x, y);});
  }
};

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    template <typename, typename> typename CompatibilityMatrix,
    typename IndexOrderG>
class ullmann_state_mono
  : public ullmann_state_base<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        CompatibilityMatrix,
        IndexOrderG> {
 private:
  using base = ullmann_state_base<
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
  
  bool possible(IndexG i) {
    for (IndexH j=0; j<n; ++j) {
      if (M.get(i, j)) {
        return true;
      }
    }
    return false;
  }
  
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
    for (IndexG ii=0; ii<m; ++ii) {
      if (g.edge(i, ii)) {
        bool exists = false;
        for (IndexH jj=0; jj<n; ++jj) {
          if (h.edge(j, jj) && M.get(ii, jj)) {
            exists = true;
            break;
          }
        }
        if (!exists) {
          return false;
        }
      }
      if (g.edge(ii, i)) {
        bool exists = false;
        for (IndexH jj=0; jj<n; ++jj) {
          if (h.edge(jj, j) && M.get(ii, jj)) {
            exists = true;
            break;
          }
        }
        if (!exists) {
          return false;
        }
      }
    }
    return true;
  }
  
  bool refine() {
    bool change;
    do {
      change = false;
      for (IndexG i=0; i<m; ++i) {
        for (IndexH j=0; j<n; ++j) {
          if (M.get(i,j) && !ullmann_condition(i, j)) {
            M.unset(i, j);
            if (!possible(i)) {
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
  ullmann_state_mono(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : ullmann_state_base<
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

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    template <typename, typename> typename CompatibilityMatrix,
    typename IndexOrderG>
class ullmann_state_ind
  : public ullmann_state_base<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        CompatibilityMatrix,
        IndexOrderG> {
 private:
  using base = ullmann_state_base<
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

  bool possible(IndexG i) {
    for (IndexH j=0; j<n; ++j) {
      if (M.get(i, j)) {
        return true;
      }
    }
    return false;
  }

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
    for (IndexG ii=0; ii<m; ++ii) {
      bool out_g = g.edge(i, ii);
      bool in_g = g.edge(ii, i);
      bool exists_out = false;
      bool exists_in = false;
      for (IndexH jj=0; jj<n; ++jj) {
        if (M.get(ii, jj)) {
          bool out_h = h.edge(j, jj);
          bool in_h = h.edge(jj, j);
          if (out_g == out_h) {
            exists_out = true;
          }
          if (in_g == in_h) {
            exists_in = true;
          }
          if (exists_out && exists_in) {
            break;
          }
        }
      }
      if (!(exists_out && exists_in)) {
        return false;
      }
    }
    return true;
  }
  
  bool refine() {
    bool change;
    do {
      change = false;
      for (IndexG i=0; i<m; ++i) {
        for (IndexH j=0; j<n; ++j) {
          if (M.get(i, j) && !ullmann_condition(i, j)) {
            M.unset(i, j);
            if (!possible(i)) {
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
  ullmann_state_ind(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : ullmann_state_base<
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


#endif  // ULLMANN_STATE_H
