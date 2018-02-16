#ifndef ULLIMP_STATE_H_
#define ULLIMP_STATE_H_

#include <iterator>
#include <vector>
#include <unordered_set>

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderG>
class ullimp_state_base {
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

  std::vector<std::unordered_set<IndexH>> M;

 public:
  ullimp_state_base(
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
        M(m) {
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
  
  ullimp_state_base(ullimp_state_base const &) = delete;

  bool empty() const {
    return x_it == std::begin(index_order_g);
  }
  bool full() const {
    return x_it == std::end(index_order_g);
  }

  void advance() {
  }
  void revert() {
  }

  void push(IndexH y) {
    ++x_it;
  }
  void pop() {
    --x_it;
  }

  auto const & candidates() const {
    return M[*x_it];
  }
  
  bool assign(IndexH y) {
    return true;
  }
};

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderG>
class ullimp_state_ind
  : public ullimp_state_base<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        IndexOrderG> {
 private:
  using base = ullimp_state_base<
      G,
      H,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
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
  
  std::vector<std::vector<std::pair<IndexG,IndexH>>> changes;

 public:
  ullimp_state_ind(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : ullimp_state_base<
            G,
            H,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            IndexOrderG>(g, h, vertex_comp, edge_comp, index_order_g),
        changes(m) {
  }
  
  void push(IndexH y) {
    auto x = *x_it;
    for (auto u_it=std::next(x_it); u_it!=std::end(index_order_g); ++u_it) {
      auto u = *u_it;
      {
        auto y_pos = M[u].find(y);
        if (y_pos != M[u].end()) {
          M[u].erase(y_pos);
          changes[x].emplace_back(u, y);
        }
      }
      bool out_g = g.edge(x, u);
      bool in_g = g.edge(u, x);
      for (auto v_it=M[u].begin(); v_it!=M[u].end(); ) {
        auto v = *v_it;
        bool out_h = h.edge(y, v);
        bool in_h = h.edge(v, y);
        if (out_g != out_h || in_g != in_h) {
          v_it = M[u].erase(v_it);
          changes[x].emplace_back(u, v);
        } else {
          ++v_it;
        }
      }
    }
    
    base::push(y);
  }
  
  void pop() {
    base::pop();
    
    auto x = *x_it;
    for (auto p : changes[x]) {
      // invariant: p.first comes after x in index_order_g (push() makes sure of that),
      // thus this doesn't invalidate any iterators held at previous levels of explore
      M[p.first].insert(p.second);
    }
    changes[x].clear();
  }
};

#endif  // ULLIMP_STATE_H
