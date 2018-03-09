#ifndef DYNAMIC_MAT_STATE_H_
#define DYNAMIC_MAT_STATE_H_

#include <iterator>
#include <vector>
#include <stack>
#include <set>
#include <unordered_set>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/adaptor/filtered.hpp>

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename CompatibilityMatrix>
class dynamic_mat_state_base {
 protected:
  using IndexG = typename G::index_type;
  using IndexH = typename H::index_type;
   
  IndexG m;
  IndexH n;
  
  G const & g;
  H const & h;

  VertexEquivalencePredicate vertex_comp;
  EdgeEquivalencePredicate edge_comp;
  
  std::vector<IndexH> map;
  std::vector<IndexG> inv;
  
  std::stack<IndexG> x_st;
  std::unordered_set<IndexG> available;

  CompatibilityMatrix M;
  
  using H_adjacent_vertices_container_type = typename H::adjacent_vertices_container_type;
  H_adjacent_vertices_container_type h_vertices;

  std::vector<std::pair<IndexH,bool>> h_parents;

 public:
  dynamic_mat_state_base(
      G const & g,
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
        M(m, n),
        h_vertices(n),
        h_parents(m) {
    for (IndexG i=0; i<m; ++i) {
      available.insert(i);
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
  
  dynamic_mat_state_base(dynamic_mat_state_base const &) = delete;

  bool empty() const {
    return available.size() == m;
  }
  
  bool full() const {
    return available.size() == 0;
  }

  void prepare() {
    auto x_it = std::min_element(
        available.begin(),
        available.end(),
        [this](auto const & a, auto const & b) {
          return M.num_candidates(a) < M.num_candidates(b);
        });
    x_st.push(*x_it);
    available.erase(*x_it);
    
    auto x = x_st.top();
    
    auto const & x_adj = g.adjacent_vertices(x);
    auto g_parent_it = std::find_if(
        std::begin(x_adj),
        std::end(x_adj),
        [this](auto u) {
          return map[u] != n;
        });
    if (g_parent_it != std::end(x_adj)) {
      h_parents[x] = {map[*g_parent_it], false};
    } else {
      auto const & x_inv_adj = g.inv_adjacent_vertices(x);
      g_parent_it = std::find_if(
          std::begin(x_inv_adj),
          std::end(x_inv_adj),
          [this](auto u) {
            return map[u] != n;
          });
      if (g_parent_it != std::end(x_inv_adj)) {
        h_parents[x] = {map[*g_parent_it], true};
      } else {
        h_parents[x] = {n, false};
      }
    }
  }
  
  void forget() {
    auto x = x_st.top();
    x_st.pop();
    available.insert(x);
  }
  
  auto const & candidates() const {
    auto x = x_st.top();
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
    auto x = x_st.top();
    map[x] = y;
    inv[y] = x;
  }
  
  void pop() {
    auto x = x_st.top();
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
class dynamic_mat_state_ind
  : public dynamic_mat_state_base<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        CompatibilityMatrix> {
 private:
  using base = dynamic_mat_state_base<
      G,
      H,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      CompatibilityMatrix>;

 protected:
  using IndexG = typename base::IndexG;
  using IndexH = typename base::IndexH;

  using base::m;
  using base::n;
  using base::g;
  using base::h;
  using base::map;
  using base::inv;
  using base::x_st;
  using base::available;
  using base::M;

  void filter_after(IndexG u, IndexH v) {
    for (auto i : available) {
      M.unset(i, v);
    }
  }

  void neighborhood_filter_after(IndexG u, IndexH v) {
    // TODO premakni v mono in extend
    for (auto i : g.adjacent_vertices(u)) {
      if (map[i] == n) {
        for (auto j : h.not_adjacent_vertices(v)) {
          if (inv[j] == m) {
            M.unset(i, j);
          }
        }
      }
    }
    for (auto i : g.inv_adjacent_vertices(u)) {
      if (map[i] ==  n) {
        for (auto j : h.not_inv_adjacent_vertices(v)) {
          if (inv[j] == m) {
            M.unset(i, j);
          }
        }
      }
    }
    
    for (auto j : h.adjacent_vertices(v)) {
      if (inv[j] == m) {
        for (auto i : g.not_adjacent_vertices(u)) {
          if (map[i] == n) {
            M.unset(i, j);
          }
        }
      }
    }
    for (auto j : h.inv_adjacent_vertices(v)) {
      if (inv[j] == m) {
        for (auto i : g.not_inv_adjacent_vertices(u)) {
          if (map[i] == n) {
            M.unset(i, j);
          }
        }
      }
    } 
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
  
  void refine() {
    bool changed;
    do {
      changed = false;
      for (IndexG u=0; u<m; ++u) {
        for (IndexH v=0; v<n; ++v) {
          if (M.get(u, v) && !ullmann_condition(u, v)) {
            M.unset(u, v);
            changed = true;
          }
        }
      }
    } while(changed);
  }
  
  bool partial_ullmann_condition(IndexG u, IndexH v) {
    for (auto i : g.adjacent_vertices(u)) {
      if (map[i] == n) {
        bool all_false = true;
        for (auto j : h.adjacent_vertices(v)) {
          if (M.get(i, j) && inv[j] == m) {
            all_false = false;
            break;
          }
        }
        if (all_false) {
          return false;
        }
      }
    }
    for (auto i : g.inv_adjacent_vertices(u)) {
      if (map[i] == n) {
        bool all_false = true;
        for (auto j : h.inv_adjacent_vertices(v)) {
          if (M.get(i, j) && inv[j] == m) {
            all_false = false;
            break;
          }
        }
        if (all_false) {
          return false;
        }
      }
    }
    return true;
  }
  
  void partial_refine(IndexG u, IndexH v) {
    for (auto i : g.adjacent_vertices(u)) {
      if (map[i] == n) {
        for (auto j : h.adjacent_vertices(v)) {
          if (M.get(i, j) && inv[j] == m && !partial_ullmann_condition(i, j)) {
            M.unset(i, j);
          }
        }
      }
    }
    for (auto i : g.inv_adjacent_vertices(u)) {
      if (map[i] == n) {
        for (auto j : h.inv_adjacent_vertices(v)) {
          if (M.get(i, j) && inv[j] == m && !partial_ullmann_condition(i, j)) {
            M.unset(i, j);
          }
        }
      }
    }
  }

 public:
  dynamic_mat_state_ind(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp)
      : dynamic_mat_state_base<
            G,
            H,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            CompatibilityMatrix>(g, h, vertex_comp, edge_comp) {
    refine();
  }
  
  bool assign(IndexH y) {
    auto x = x_st.top();
    return M.get(x, y);
  }
  
  void push(IndexH y) {
    auto x = x_st.top();
    filter_after(x, y);
    neighborhood_filter_after(x, y);
    if (available.size() >= m/2) {
      partial_refine(x, y);
    }
    base::push(y);
  }
};

#endif  // DYNAMIC_MAT_STATE_H
