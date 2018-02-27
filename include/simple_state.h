#ifndef SIMPLE_STATE_H_
#define SIMPLE_STATE_H_

#include <iterator>
#include <vector>
#include <algorithm>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/adaptor/filtered.hpp>

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderG>
class simple_state_mono {
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

  std::vector<IndexH> map;
  std::vector<IndexG> inv;
  
  bool topology_condition(IndexG u, IndexH v) {
    for (auto i : g.adjacent_vertices(u)) {
      auto j = map[i];
      if (j != n) {
        if (!h.edge(v, j) || !edge_comp(u, i, v, j)) {
          return false;
        }
      }
    }
    for (auto i : g.inv_adjacent_vertices(u)) {
      auto j = map[i];
      if (j != n) {
        if (!h.edge(j, v) || !edge_comp(i, u, j, v)) {
          return false;
        }
      }
    }
    return true;
  }

 public:
  simple_state_mono(
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
        map(m, n),
        inv(n, m) {
  }
  
  simple_state_mono(simple_state_mono const &) = delete;

  bool empty() {
    return x_it == std::begin(index_order_g);
  }
  
  bool full() {
    return x_it == std::end(index_order_g);
  }
  
  auto candidates() {
    boost::counting_iterator<IndexH> begin{0}, end{n};
    return boost::adaptors::filter(
        boost::make_iterator_range(begin, end),
        [this](auto y){return inv[y] == m;});
  }

  void advance() {
  }
  
  void revert() {
  }

  bool assign(IndexH y) {
    auto x = *x_it;
    return
        vertex_comp(x, y) &&
        g.out_degree(x) <= h.out_degree(y) &&
        g.in_degree(x) <= h.in_degree(y) &&
        topology_condition(x, y);
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
class simple_state_ind
  : public simple_state_mono<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        IndexOrderG> {
 private:
  using base = simple_state_mono<
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
  using base::x_it;
  using base::map;
  using base::inv;
  
  std::vector<IndexG> g_out_count;
  std::vector<IndexG> g_in_count;
  
  std::vector<IndexH> h_out_count;
  std::vector<IndexH> h_in_count;
 
 public:
  simple_state_ind(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : simple_state_mono<
            G,
            H,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            IndexOrderG>(g, h, vertex_comp, edge_comp, index_order_g),
        g_out_count(m),
        g_in_count(m),
        h_out_count(n),
        h_in_count(n) {
    std::vector<IndexG> index_pos_g(m);
    for (IndexG i=0; i<m; ++i) {
      index_pos_g[index_order_g[i]] = i;
    }
    for (IndexG i=0; i<m; ++i) {
      auto i_pos = index_pos_g[i];
      auto const & i_adj = g.adjacent_vertices(i);
      g_out_count[i] = std::count_if(std::begin(i_adj), std::end(i_adj), [i_pos, &index_pos_g](auto ii) {
        return index_pos_g[ii] < i_pos;
      });
      auto const & i_inv_adj = g.inv_adjacent_vertices(i);
      g_in_count[i] = std::count_if(std::begin(i_inv_adj), std::end(i_inv_adj), [i_pos, &index_pos_g](auto ii) {
        return index_pos_g[ii] < i_pos;
      });
    }
  }

  bool assign(IndexH y) {
    auto x = *x_it;
    return
        g_out_count[x] == h_out_count[y] &&
        g_in_count[x] == h_in_count[y] &&
        base::assign(y);
  }
  
  void push(IndexH y) {
    for (auto j : h.adjacent_vertices(y)) {
      ++h_in_count[j];
    }
    for (auto j : h.inv_adjacent_vertices(y)) {
      ++h_out_count[j];
    }
    base::push(y);
  }
  
  void pop() {
    auto y = base::pop();
    for (auto j : h.adjacent_vertices(y)) {
      --h_in_count[j];
    }
    for (auto j : h.inv_adjacent_vertices(y)) {
      --h_out_count[j];
    }
  }
};

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderG>
class simple_state_ind2
  : public simple_state_mono<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        IndexOrderG> {
 private:
  using base = simple_state_mono<
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
  using base::x_it;
  using base::map;
  using base::inv;
 
 public:
  using base::simple_state_mono;
 
  bool assign(IndexH y) {
    if (!base::assign(y)) {
      return false;
    }
    auto x = *x_it;
    auto const & x_adj = g.adjacent_vertices(x);
    auto count_g = std::count_if(std::begin(x_adj), std::end(x_adj), [this](auto i) {
      return map[i] != n;
    });
    auto const & y_adj = h.adjacent_vertices(y);
    auto count_h = std::count_if(std::begin(y_adj), std::end(y_adj), [this](auto j) {
      return inv[j] != m;
    });
    if (count_g != count_h) {
      return false;
    } else {
      auto const & x_inv_adj = g.inv_adjacent_vertices(x);
      count_g = std::count_if(std::begin(x_inv_adj), std::end(x_inv_adj), [this](auto i) {
        return map[i] != n;
      });
      auto const & y_inv_adj = h.inv_adjacent_vertices(y);
      count_h = std::count_if(std::begin(y_inv_adj), std::end(y_inv_adj), [this](auto j) {
        return inv[j] != m;
      });
      return count_g == count_h;
    }
  }
};

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderG>
class simple_state_ind3
  : public simple_state_mono<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        IndexOrderG> {
 private:
  using base = simple_state_mono<
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
  using base::x_it;
  using base::map;
  using base::inv;
  
  std::vector<IndexG> g_out_count;
  std::vector<IndexG> g_in_count;
 
  bool relative_degree_condition(IndexG u, IndexH v) {
    auto const & v_adj = h.adjacent_vertices(v);
    auto count_h = std::count_if(std::begin(v_adj), std::end(v_adj), [this](auto j) {
      return inv[j] != m;
    });
    if (g_out_count[u] != count_h) {
      return false;
    } else {
      auto const & v_inv_adj = h.inv_adjacent_vertices(v);
      count_h = std::count_if(std::begin(v_inv_adj), std::end(v_inv_adj), [this](auto j) {
        return inv[j] != m;
      });
      return g_in_count[u] == count_h;
    }
  }
 
 public:
  simple_state_ind3(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : simple_state_mono<
            G,
            H,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            IndexOrderG>(g, h, vertex_comp, edge_comp, index_order_g),
        g_out_count(m),
        g_in_count(m) {
    std::vector<IndexG> index_pos_g(m);
    for (IndexG i=0; i<m; ++i) {
      index_pos_g[index_order_g[i]] = i;
    }
    for (IndexG i=0; i<m; ++i) {
      auto i_pos = index_pos_g[i];
      auto const & i_adj = g.adjacent_vertices(i);
      g_out_count[i] = std::count_if(std::begin(i_adj), std::end(i_adj), [i_pos, &index_pos_g](auto ii) {
        return index_pos_g[ii] < i_pos;
      });
      auto const & i_inv_adj = g.inv_adjacent_vertices(i);
      g_in_count[i] = std::count_if(std::begin(i_inv_adj), std::end(i_inv_adj), [i_pos, &index_pos_g](auto ii) {
        return index_pos_g[ii] < i_pos;
      });
    }
  }
 
  bool assign(IndexH y) {
    auto x = *x_it;
    return base::assign(y) && relative_degree_condition(x, y);
  }
};

#endif  // SIMPLE_STATE_H
