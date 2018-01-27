#ifndef ULLMANN_STATE_H_
#define ULLMANN_STATE_H_

#include <iterator>
#include <stack>

template <
    typename AdjacencyMatrixSmall,
    typename AdjacencyMatrixLarge,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    template <typename, typename> typename CompatibilityMatrix,
    typename IndexOrderSmall>
class ullmann_state_base {
 protected:
  using IndexSmall = typename AdjacencyMatrixSmall::index_type;
  using IndexLarge = typename AdjacencyMatrixLarge::index_type;
   
  IndexSmall m;
  IndexLarge n;

  AdjacencyMatrixSmall const & g;
  AdjacencyMatrixLarge const & h;

  VertexEquivalencePredicate vertex_comp;
  EdgeEquivalencePredicate edge_comp;

  CompatibilityMatrix<IndexSmall, IndexLarge> M;

  IndexOrderSmall const & index_order_small;
  typename IndexOrderSmall::const_iterator x_it;
  IndexSmall x;

  std::stack<IndexLarge> y_st;
  IndexLarge y;

 public:
  ullmann_state_base(
      AdjacencyMatrixSmall const & g,
      AdjacencyMatrixLarge const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderSmall const & index_order_small)
      : m{g.num_vertices()},
        n{h.num_vertices()},
        g{g},
        h{h},
        vertex_comp{vertex_comp},
        edge_comp{edge_comp},
        M{m, n},
        index_order_small{index_order_small} {
    for (IndexSmall i=0; i<m; ++i) {
      for (IndexLarge j=0; j<n; ++j) {
        if (vertex_comp(i, j) &&
            g.out_degree(i) <= h.out_degree(j) &&
            g.in_degree(i) <= h.in_degree(j)) {
          M.set(i, j);
        }
      }
    }
    
    x_it = std::begin(index_order_small);
    x = *x_it;
    y = 0;
    while (y < n && !M.get(x, y)) {
      ++y;
    }
  }

  bool empty() {
    return x_it == std::begin(index_order_small);
  }
  bool full() {
    return x_it == std::end(index_order_small);
  }

  void advance() {
    M.advance();
  }
  void revert() {
    M.revert();
  }

  void push() {
    ++x_it;
    y_st.push(y);
    if (x_it != std::end(index_order_small)) {
      x = *x_it;
      y = 0;
      while (y < n && !M.get(x, y)) {
        ++y;
      }
    }
  }
  void pop() {
    --x_it;
    x = *x_it;
    y = y_st.top();
    y_st.pop();
  }

  bool available() {
    return y != n;
  }
  void next() {
    do {
      ++y;
    } while (y < n && !M.get(x, y));
  }
};

template <
    typename AdjacencyMatrixSmall,
    typename AdjacencyMatrixLarge,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    template <typename, typename> typename CompatibilityMatrix,
    typename IndexOrderSmall>
class ullmann_state_mono
  : public ullmann_state_base<
        AdjacencyMatrixSmall,
        AdjacencyMatrixLarge,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        CompatibilityMatrix,
        IndexOrderSmall> {
 private:
  using base = ullmann_state_base<
      AdjacencyMatrixSmall,
      AdjacencyMatrixLarge,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      CompatibilityMatrix,
      IndexOrderSmall>;

 protected:
  using IndexSmall = typename base::IndexSmall;
  using IndexLarge = typename base::IndexLarge;

  using base::m;
  using base::n;
  using base::g;
  using base::h;
  using base::M;
  using base::x;
  using base::y;
  
  bool possible(IndexSmall i) {
    for (IndexLarge j=0; j<n; ++j) {
      if (M.get(i, j)) {
        return true;
      }
    }
    return false;
  }
  
  void filter(IndexSmall i, IndexLarge j) {
    for (IndexSmall ii=0; ii<m; ++ii) {
      M.unset(ii, j);
    }
    for (IndexLarge jj=0; jj<n; ++jj) {
      M.unset(i, jj);
    }
    M.set(i, j);
  }

  bool ullmann_condition(IndexSmall i, IndexLarge j) {
    for (IndexSmall ii=0; ii<m; ++ii) {
      if (g.edge(i, ii)) {
        bool exists = false;
        for (IndexLarge jj=0; jj<n; ++jj) {
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
        for (IndexLarge jj=0; jj<n; ++jj) {
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
      for (IndexSmall i=0; i<m; ++i) {
        for (IndexLarge j=0; j<n; ++j) {
          if (M.get(i,j) && !ullmann_condition(i, j)) {
            M.unset(i, j);
            if(!possible(i)) {
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
      AdjacencyMatrixSmall const & g,
      AdjacencyMatrixLarge const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderSmall const & index_order_small)
      : ullmann_state_base<
            AdjacencyMatrixSmall,
            AdjacencyMatrixLarge,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            CompatibilityMatrix,
            IndexOrderSmall>(g, h, vertex_comp, edge_comp, index_order_small) {
    refine();
    while (y < n && !M.get(x, y)) {
      ++y;
    }
  }

  bool assign() {
    filter(x, y);
    bool success = refine();
    return success;
  }
};

template <
    typename AdjacencyMatrixSmall,
    typename AdjacencyMatrixLarge,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    template <typename, typename> typename CompatibilityMatrix,
    typename IndexOrderSmall>
class ullmann_state_ind
  : public ullmann_state_base<
        AdjacencyMatrixSmall,
        AdjacencyMatrixLarge,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        CompatibilityMatrix,
        IndexOrderSmall> {
 private:
  using base = ullmann_state_base<
      AdjacencyMatrixSmall,
      AdjacencyMatrixLarge,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      CompatibilityMatrix,
      IndexOrderSmall>;

 protected:
  using IndexSmall = typename base::IndexSmall;
  using IndexLarge = typename base::IndexLarge;

  using base::m;
  using base::n;
  using base::g;
  using base::h;
  using base::M;
  using base::x;
  using base::y;

  bool possible(IndexSmall i) {
    for (IndexLarge j=0; j<n; ++j) {
      if (M.get(i, j)) {
        return true;
      }
    }
    return false;
  }

  void filter(IndexSmall i, IndexLarge j) {
    for (IndexSmall ii=0; ii<m; ++ii) {
      M.unset(ii, j);
    }
    for (IndexLarge jj=0; jj<n; ++jj) {
      M.unset(i, jj);
    }
    M.set(i, j);
  }

  bool ullmann_condition(IndexSmall i, IndexLarge j) {
    for (IndexSmall ii=0; ii<m; ++ii) {
      bool out_g = g.edge(i, ii);
      bool in_g = g.edge(ii, i);
      bool exists_out = false;
      bool exists_in = false;
      for (IndexLarge jj=0; jj<n; ++jj) {
        if (M.get(ii, jj)) {
          bool out_h = h.edge(j, jj);
          bool in_h = h.edge(jj, j);
          if(out_g == out_h) {
            exists_out = true;
          }
          if(in_g == in_h) {
            exists_in = true;
          }
          if(exists_out && exists_in) {
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
      for (IndexSmall i=0; i<m; ++i) {
        for (IndexLarge j=0; j<n; ++j) {
          if (M.get(i, j) && !ullmann_condition(i, j)) {
            M.unset(i, j);
            if(!possible(i)) {
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
      AdjacencyMatrixSmall const & g,
      AdjacencyMatrixLarge const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderSmall const & index_order_small)
      : ullmann_state_base<
            AdjacencyMatrixSmall,
            AdjacencyMatrixLarge,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            CompatibilityMatrix,
            IndexOrderSmall>(g, h, vertex_comp, edge_comp, index_order_small) {
    refine();
    while (y < n && !M.get(x, y)) {
      ++y;
    }
  }

  bool assign() {
    filter(x, y);
    bool success = refine();
    return success;
  }
};


#endif  // ULLMANN_STATE_H
