#ifndef DYNAMIC_STATE_H_
#define DYNAMIC_STATE_H_

#include <iterator>
#include <vector>
#include <stack>
#include <set>

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
class dynamic_state_base {
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
  std::set<IndexG> available;

  std::vector<std::set<IndexH>> M;
  std::vector<std::stack<std::pair<IndexG,IndexH>>> changes;

 public:
  dynamic_state_base(
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
        M(m),
        changes(m) {
    for (IndexG i=0; i<m; ++i) {
      available.insert(i);
      for (IndexH j=0; j<n; ++j) {
        if (vertex_comp(i, j) &&
            g.out_degree(i) <= h.out_degree(j) &&
            g.in_degree(i) <= h.in_degree(j)) {
          M[i].insert(j);
        }
      }
    }
  }
  
  dynamic_state_base(dynamic_state_base const &) = delete;

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
          return M[a].size() < M[b].size();
        });
    x_st.push(*x_it);
    available.erase(*x_it);
  }
  
  void forget() {
    auto x = x_st.top();
    x_st.pop();
    available.insert(x);
  }
  
  auto const & candidates() const {
    auto x = x_st.top();
    return M[x];
  }

  void advance() {
  }
  
  void revert() {
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
    typename EdgeEquivalencePredicate>
class dynamic_state_ind
  : public dynamic_state_base<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate> {
 private:
  using base = dynamic_state_base<
      G,
      H,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate>;

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
  using base::changes;

  void filter_after(IndexG u, IndexH v) {
    for (auto i_it=available.begin(); i_it!=available.end(); ++i_it) {
      auto i = *i_it;
      auto v_pos = M[i].find(v);
      if (v_pos != M[i].end()) {
        M[i].erase(v_pos);
        changes[u].emplace(i, v);
      }
    }
  }

  void neighborhood_filter_after(IndexG u, IndexH v) {
    // TODO premakni v mono in extend
    for (auto i : g.adjacent_vertices(u)) {
      if (map[i] == n) {
        for (auto j : h.not_adjacent_vertices(v)) {
          if (inv[j] == m) {
            auto j_pos = M[i].find(j);
            if (j_pos != M[i].end()) {
              M[i].erase(j);
              changes[u].emplace(i, j);
            }
          }
        }
      }
    }
    for (auto i : g.inv_adjacent_vertices(u)) {
      if (map[i] ==  n) {
        for (auto j : h.not_inv_adjacent_vertices(v)) {
          if (inv[j] == m) {
            auto j_pos = M[i].find(j);
            if (j_pos != M[i].end()) {
              M[i].erase(j);
              changes[u].emplace(i, j);
            }
          }
        }
      }
    }
    
    for (auto j : h.adjacent_vertices(v)) {
      if (inv[j] == m) {
        for (auto i : g.not_adjacent_vertices(u)) {
          if (map[i] == n) {
            auto j_pos = M[i].find(j);
            if (j_pos != M[i].end()) {
              M[i].erase(j);
              changes[u].emplace(i, j);
            }
          }
        }
      }
    }
    for (auto j : h.inv_adjacent_vertices(v)) {
      if (inv[j] == m) {
        for (auto i : g.not_inv_adjacent_vertices(u)) {
          if (map[i] == n) {
            auto j_pos = M[i].find(j);
            if (j_pos != M[i].end()) {
              M[i].erase(j);
              changes[u].emplace(i, j);
            }
          }
        }
      }
    } 
  }

 public:
  dynamic_state_ind(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp)
      : dynamic_state_base<
            G,
            H,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate>(g, h, vertex_comp, edge_comp) {
  }
  
  bool assign(IndexH y) {
    return true;
  }
  
  void push(IndexH y) {
    auto x = x_st.top();
    filter_after(x, y);
    neighborhood_filter_after(x, y);
    base::push(y);
  }
  
  void pop() {
    base::pop();
    
    auto x = x_st.top();
    while(!changes[x].empty()) {
      auto p = changes[x].top();
      changes[x].pop();
      // invariant: p.first comes after x in index_order_g (push() makes sure of that),
      // thus this doesn't invalidate any iterators held at previous levels of explore
      M[p.first].insert(p.second);
    }
  }
};

#endif  // DYNAMIC_STATE_H
