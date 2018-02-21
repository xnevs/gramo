#ifndef ULLIMP4_STATE_H_
#define ULLIMP4_STATE_H_

#include <iterator>
#include <vector>
#include <stack>
#include <set>

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderG>
class ullimp4_state_base {
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

  std::vector<std::set<IndexH>> M;
  std::vector<std::stack<std::pair<IndexG,IndexH>>> changes;
  
  std::vector<IndexG> map;
  std::vector<IndexG> inv;

 public:
  ullimp4_state_base(
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
        M(m),
        changes(m),
        map(m, n),
        inv(n, m) {
    for (IndexG i=0; i<m; ++i) {
      for (IndexH j=0; j<n; ++j) {
        if (vertex_comp(i, j) &&
            g.out_degree(i) <= h.out_degree(j) &&
            g.in_degree(i) <= h.in_degree(j)) {
          M[i].insert(j);
        }
      }
    }
  }
  
  ullimp4_state_base(ullimp4_state_base const &) = delete;

  bool empty() const {
    return x_it == std::begin(index_order_g);
  }
  
  bool full() const {
    return x_it == std::end(index_order_g);
  }
  
  auto const & candidates() const {
    return M[*x_it];
  }

  void advance() {
  }
  
  void revert() {
  }

  void push(IndexH y) {
    auto x = *x_it;
    map[x] = y;
    inv[y] = x;
    ++x_it;
  }
  
  IndexH pop() {
    --x_it;
    auto x = *x_it;
    auto y = map[x];
    map[x] = n;
    inv[y] = m;
    return y;
  }
};

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderG>
class ullimp4_state_mono
  : public ullimp4_state_base<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        IndexOrderG> {
 private:
  using base = ullimp4_state_base<
      G,
      H,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
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
  using base::changes;
  using base::map;
  using base::inv;
  
  void filter_after(x_it_type u_it, IndexH v) {
    auto u = *u_it;
    for (auto i_it=std::next(u_it); i_it!=std::end(index_order_g); ++i_it) {
      auto i = *i_it;
      auto v_pos = M[i].find(v);
      if (v_pos != M[i].end()) {
        M[i].erase(v_pos);
        changes[u].emplace(i, v);
      }
    }
  }
  
  void neighborhood_filter_after(IndexG u, IndexH v) {
    for (auto i : g.adjacent_vertices_after(u)) {
      for (auto j_it=M[i].begin(); j_it!=M[i].end(); ) {
        auto j = *j_it;
        if (!h.edge(v, j)) {
          j_it = M[i].erase(j_it);
          changes[u].emplace(i, j);
        } else {
          ++j_it;
        }
      }
    }
    for (auto i : g.inv_adjacent_vertices_after(u)) {
      for (auto j_it=M[i].begin(); j_it!=M[i].end(); ) {
        auto j = *j_it;
        if (!h.edge(j, v)) {
          j_it = M[i].erase(j_it);
          changes[u].emplace(i, j);
        } else {
          ++j_it;
        }
      }
    }
  }

 public:
  ullimp4_state_mono(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : ullimp4_state_base<
            G,
            H,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            IndexOrderG>(g, h, vertex_comp, edge_comp, index_order_g) {
  }
  
  bool assign(IndexH y) {
    return true;
  }
  
  void push(IndexH y) {
    auto x = *x_it;
    filter_after(x_it, y);
    neighborhood_filter_after(x, y);
    base::push(y);
  }
  
  IndexH pop() {
    auto y = base::pop();
    auto x = *x_it;
    while(!changes[x].empty()) {
      auto p = changes[x].top();
      changes[x].pop();
      // invariant: p.first comes after x in index_order_g (push() makes sure of that),
      // thus this doesn't invalidate any iterators held at previous levels of explore
      M[p.first].insert(p.second);
    }
    return y;
  }
};

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderG>
class ullimp4_state_ind
  : public ullimp4_state_mono<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        IndexOrderG> {
 private:
  using base = ullimp4_state_mono<
      G,
      H,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
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
  using base::changes;
  
  std::vector<IndexH> h_out_degree_mapped;
  std::vector<IndexH> h_in_degree_mapped;

 public:
  ullimp4_state_ind(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : ullimp4_state_mono<
            G,
            H,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            IndexOrderG>(g, h, vertex_comp, edge_comp, index_order_g),
        h_out_degree_mapped(n),
        h_in_degree_mapped(n) {
  }
  
  bool assign(IndexH y) {
    auto x = *x_it;
    return
        g.out_degree_before(x) == h_out_degree_mapped[y] &&
        g.in_degree_before(x) == h_in_degree_mapped[y] &&
        base::assign(y);
  }
  
  void push(IndexH y) {
    for (auto v : h.adjacent_vertices(y)) {
      ++h_in_degree_mapped[v];
    }
    for (auto v : h.inv_adjacent_vertices(y)) {
      ++h_out_degree_mapped[v];
    }
    base::push(y);
  }
  
  void pop() {
    auto y = base::pop();
    for (auto v : h.adjacent_vertices(y)) {
      --h_in_degree_mapped[v];
    }
    for (auto v : h.inv_adjacent_vertices(y)) {
      --h_out_degree_mapped[v];
    }
  }
};

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderG>
class ullimp4_state_ind2
  : public ullimp4_state_mono<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        IndexOrderG> {
 private:
  using base = ullimp4_state_mono<
      G,
      H,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
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
  using base::changes;
  using base::map;
  using base::inv;

 public:
  ullimp4_state_ind2(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : ullimp4_state_mono<
            G,
            H,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            IndexOrderG>(g, h, vertex_comp, edge_comp, index_order_g) {
  }
  
  bool assign(IndexH y) {
    auto x = *x_it;
    IndexH h_out_degree_before_y = 0;
    IndexH h_in_degree_before_y = 0;
    for (auto v : h.adjacent_vertices(y)) {
      if (inv[v] != m) {
        ++h_out_degree_before_y;
      }
    }
    for (auto v : h.inv_adjacent_vertices(y)) {
      if (inv[v] != m) {
        ++h_in_degree_before_y;
      }
    }
    return
        g.out_degree_before(x) == h_out_degree_before_y &&
        g.in_degree_before(x) == h_in_degree_before_y;
        
  }
};

#endif  // ULLIMP4_STATE_H
