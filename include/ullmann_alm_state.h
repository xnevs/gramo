#ifndef ULLMANN_ALM_STATE_H_
#define ULLMANN_ALM_STATE_H_

#include <iterator>

#include "ullmann_state.h"

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    template <typename, typename> typename CompatibilityMatrix,
    typename IndexOrderG>
class ullmann_alm_state_mono
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
  
  void filter(IndexG i, IndexH j) {
    for (IndexG ii=0; ii<m; ++ii) {
      M.unset(ii, j);
    }
    for (IndexH jj=0; jj<n; ++jj) {
      M.unset(i, jj);
    }
    M.set(i, j);
  }
  
  bool ullmann_condition(IndexG u, IndexH v) {
    for (auto i : g.adjacent_vertices(u)) {
      bool all_false = true;
      for (auto j : h.adjacent_vertices(v)) {
        if (M.get(i, j)) {
          all_false = false;
          break;
        }
      }
      if (all_false) {
        return false;
      }
    }
    for (auto i : g.inv_adjacent_vertices(u)) {
      bool all_false = true;
      for (auto j : h.inv_adjacent_vertices(v)) {
        if (M.get(i, j)) {
          all_false = false;
          break;
        }
      }
      if (all_false) {
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
          if (M.get(i,j) && !ullmann_condition(i, j)) {
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
  ullmann_alm_state_mono(
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
class ullmann_alm_state_ind
  : public ullmann_alm_state_mono<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        CompatibilityMatrix,
        IndexOrderG> {
 private:
  using base = ullmann_alm_state_mono<
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

 public:
  ullmann_alm_state_ind(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : ullmann_alm_state_mono<
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

#endif  // ULLMANN_ALM_STATE_H
