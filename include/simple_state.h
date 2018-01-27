#ifndef SIMPLE_STATE_H_
#define SIMPLE_STATE_H_

#include <iterator>
#include <algorithm>
#include <vector>
#include <stack>
#include <unordered_set>

template <
    typename AdjacencySetSmall,
    typename AdjacencySetLarge,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderSmall>
class simple_state_mono {
 protected:
  using IndexSmall = typename AdjacencySetSmall::index_type;
  using IndexLarge = typename AdjacencySetLarge::index_type;
   
  IndexSmall m;
  IndexLarge n;
  
  AdjacencySetSmall const & g;
  AdjacencySetLarge const & h;

  VertexEquivalencePredicate vertex_comp;
  EdgeEquivalencePredicate edge_comp;

  IndexOrderSmall const & index_order_small;
  typename IndexOrderSmall::const_iterator x_it;
  IndexSmall x;

  std::stack<IndexLarge> y_st;
  IndexLarge y;

  std::vector<IndexLarge> map;
  std::vector<IndexSmall> inv;

 public:
  simple_state_mono(
      AdjacencySetSmall const & g,
      AdjacencySetLarge const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderSmall const & index_order_small)
      : m{g.num_vertices()},
        n{h.num_vertices()},
        g{g},
        h{h},
        vertex_comp{vertex_comp},
        edge_comp{edge_comp},
        index_order_small{index_order_small},
        x_it{std::begin(index_order_small)},
        x{*x_it},
        y{0},
        map(m, n),
        inv(n, m) {
  }

  bool empty() {
    return x_it == std::begin(index_order_small);
  }
  bool full() {
    return x_it == std::end(index_order_small);
  }

  void advance() {
  }
  void revert() {
  }

  void push() {
    map[x] = y;
    inv[y] = x;
    
    ++x_it;
    y_st.push(y);
    if (x_it != std::end(index_order_small)) {
      x = *x_it;
      y = 0;
      while (y < n && inv[y] != m) {
        ++y;
      }
    }
  }
  void pop() {
    --x_it;
    x = *x_it;
    y = y_st.top();
    y_st.pop();

    map[x] = n;
    inv[y] = m;
  }

  bool available() {
    return y != n;
  }
  void next() {
    do {
      ++y;
    } while (y < n && inv[y] != m);
  }

  bool assign() {
    if (!vertex_comp(x, y)) {
      return false;
    } else {
      for (auto i : g.adjacent_vertices(x)) {
        auto j = map[i];
        if (j != n) {
          if (!h.edge(y, j)) {
            return false;
          }
        }
      }
      for (auto i : g.inv_adjacent_vertices(x)) {
        auto j = map[i];
        if (j != n) {
          if (!h.edge(j, y)) {
            return false;
          }
        }
      }
      return true;
    }
  }
};

template <
    typename AdjacencySetSmall,
    typename AdjacencySetLarge,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderSmall>
class simple_state_ind
  : public simple_state_mono<
        AdjacencySetSmall,
        AdjacencySetLarge,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        IndexOrderSmall> {
 private:
  using base = simple_state_mono<
      AdjacencySetSmall,
      AdjacencySetLarge,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      IndexOrderSmall>;
      
 protected:
  using IndexSmall = typename base::IndexSmall;
  using IndexLarge = typename base::IndexLarge;

  using base::m;
  using base::n;
  using base::g;
  using base::h;
  using base::x;
  using base::y;
  using base::map;
  using base::inv;
  
  std::vector<IndexSmall> g_out_count;
  std::vector<IndexSmall> g_in_count;
  
  std::vector<IndexLarge> h_out_count;
  std::vector<IndexLarge> h_in_count;
 
 public:
  simple_state_ind(
      AdjacencySetSmall const & g,
      AdjacencySetLarge const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderSmall const & index_order_small)
      : simple_state_mono<
            AdjacencySetSmall,
            AdjacencySetLarge,
            VertexEquivalencePredicate,
            EdgeEquivalencePredicate,
            IndexOrderSmall>(g, h, vertex_comp, edge_comp, index_order_small),
        g_out_count(m),
        g_in_count(m),
        h_out_count(n),
        h_in_count(n) {
    std::vector<IndexSmall> index_pos_small(m);
    for (IndexSmall i=0; i<m; ++i) {
      index_pos_small[index_order_small[i]] = i;
    }
    for (IndexSmall i=0; i<m; ++i) {
      auto i_pos = index_pos_small[i];
      auto const & i_adj = g.adjacent_vertices(i);
      g_out_count[i] = std::count_if(std::begin(i_adj), std::end(i_adj), [i_pos, &index_pos_small](auto ii) {
        return index_pos_small[ii] < i_pos;
      });
      auto const & i_inv_adj = g.inv_adjacent_vertices(i);
      g_in_count[i] = std::count_if(std::begin(i_inv_adj), std::end(i_inv_adj), [i_pos, &index_pos_small](auto ii) {
        return index_pos_small[ii] < i_pos;
      });
    }
  }
        
 
  bool assign() {
    return base::assign() && g_out_count[x] == h_out_count[y] && g_in_count[x] == h_in_count[y];
  }
  
  void push() {
    for (auto j : h.adjacent_vertices(y)) {
      if (inv[j] == m) {
        ++h_in_count[j];
      }
    }
    for (auto j : h.inv_adjacent_vertices(y)) {
      if (inv[j] == m) {
        ++h_out_count[j];
      }
    }
    base::push();
  }
  void pop() {
    base::pop();
    for (auto j : h.adjacent_vertices(y)) {
      if (inv[j] == m) {
        --h_in_count[j];
      }
    }
    for (auto j : h.inv_adjacent_vertices(y)) {
      if (inv[j] == m) {
        --h_out_count[j];
      }
    }
  }
};

template <
    typename AdjacencySetSmall,
    typename AdjacencySetLarge,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderSmall>
class simple_state_ind2
  : public simple_state_mono<
        AdjacencySetSmall,
        AdjacencySetLarge,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        IndexOrderSmall> {
 private:
  using base = simple_state_mono<
      AdjacencySetSmall,
      AdjacencySetLarge,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      IndexOrderSmall>;
      
 protected:
  using IndexSmall = typename base::IndexSmall;
  using IndexLarge = typename base::IndexLarge;

  using base::m;
  using base::n;
  using base::g;
  using base::h;
  using base::x;
  using base::y;
  using base::map;
  using base::inv;
 
 public:
  using base::simple_state_mono;
 
  bool assign() {
    if (!base::assign()) {
      return false;
    }
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

#endif  // SIMPLE_STATE_H
