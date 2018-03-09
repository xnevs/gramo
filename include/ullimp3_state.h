#ifndef ULLIMP3_STATE_H_
#define ULLIMP3_STATE_H_

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
class ullimp3_state_base {
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

 public:
  ullimp3_state_base(
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
        changes(m) {
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
  
  ullimp3_state_base(ullimp3_state_base const &) = delete;

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
  
  auto const & candidates() const {
    return M[*x_it];
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
};

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderG>
class ullimp3_state_ind
  : public ullimp3_state_base<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        IndexOrderG> {
 private:
  using base = ullimp3_state_base<
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

  void neighborhood_filter_after(x_it_type u_it, IndexH v) {
    auto u = *u_it;
    for (auto i_it=std::next(u_it); i_it!=std::end(index_order_g); ++i_it) {
      auto i = *i_it;
      {
        auto v_pos = M[i].find(v);
        if (v_pos != M[i].end()) {
          M[i].erase(v_pos);
          changes[u].emplace(i, v);
        }
      }
      bool out_g = g.edge(u, i);
      bool in_g = g.edge(i, u);
      for (auto j_it=M[i].begin(); j_it!=M[i].end(); ) {
        auto j = *j_it;
        bool out_h = h.edge(v, j);
        bool in_h = h.edge(j, v);
        if (out_g != out_h || in_g != in_h) {
          j_it = M[i].erase(j_it);
          changes[u].emplace(i, j);
        } else {
          ++j_it;
        }
      }
    }
  }

 public:
  ullimp3_state_ind(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : ullimp3_state_base<
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
    neighborhood_filter_after(x_it, y);
    base::push(y);
  }
  
  void pop() {
    base::pop();
    
    auto x = *x_it;
    while(!changes[x].empty()) {
      auto p = changes[x].top();
      changes[x].pop();
      // invariant: p.first comes after x in index_order_g (push() makes sure of that),
      // thus this doesn't invalidate any iterators held at previous levels of explore
      M[p.first].insert(p.second);
    }
  }
};

#endif  // ULLIMP3_STATE_H
