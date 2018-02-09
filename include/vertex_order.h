#ifndef VERTEX_ORDER_H_
#define VERTEX_ORDER_H_

#include <tuple>
#include <vector>
#include <algorithm>
#include <numeric>

template <typename G>
std::vector<typename G::index_type> vertex_order_DEG(G const & g) {
  std::vector<typename G::index_type> order(g.num_vertices());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&g](auto i, auto j) {
    return g.degree(i) > g.degree(j);
  });
  return order;
}

template <typename G>
typename G::index_type clustering1(typename G::index_type i, G const & g) {
  using Index = typename G::index_type;
  Index c = 0;
  for (auto w : g.adjacent_vertices(i)) {
    auto const & w_adj = g.adjacent_vertices(w);
    c += std::count_if(std::begin(w_adj), std::end(w_adj), [i, &g](auto r) {
      return g.edge(r, i) || g.edge(i, r);
    });
    auto const & w_inv_adj = g.inv_adjacent_vertices(w);
    c += std::count_if(std::begin(w_inv_adj), std::end(w_inv_adj), [i, &g](auto r) {
      return g.edge(r, i) || g.edge(i, r);
    });
  }
  for (auto w : g.inv_adjacent_vertices(i)) {
    auto const & w_adj = g.adjacent_vertices(w);
    c += std::count_if(std::begin(w_adj), std::end(w_adj), [i, &g](auto r) {
      return g.edge(r, i) || g.edge(i, r);
    });
    auto const & w_inv_adj = g.inv_adjacent_vertices(w);
    c += std::count_if(std::begin(w_inv_adj), std::end(w_inv_adj), [i, &g](auto r) {
      return g.edge(r, i) || g.edge(i, r);
    });
  }
  return c;
}

template <typename G>
std::vector<typename G::index_type> vertex_order_RDEG_CNC(G const & g) {
  using Index = typename G::index_type;
  
  auto n = g.num_vertices();
  
  std::vector<Index> vertex_order(n);
  
  std::vector<Index> clustdeg(n);
  for (Index i=0; i<n; ++i) {
    clustdeg[i] = clustering1(i, g) + g.degree(i);
  }
  
  std::vector<bool> avail(n, true);
  for (Index idx=0; idx<n; ++idx) {
    auto bestn = n;
    int bestv = -1;
    for (Index i=0; i<n; ++i) {
      if (avail[i]) {
        auto const & i_adj = g.adjacent_vertices(i);
        auto const & i_inv_adj = g.inv_adjacent_vertices(i);
        int rdeg =
            std::count_if(std::begin(i_adj), std::end(i_adj), [&avail](auto p) {
              return !avail[p];
            })
            +
            std::count_if(std::begin(i_inv_adj), std::end(i_inv_adj), [&avail](auto p) {
              return !avail[p];
            });
        if (bestn == n || (rdeg > bestv || (rdeg == bestv && clustdeg[i] > clustdeg[bestn]))) {
          bestn = i;
          bestv = rdeg;
        }
      }
    }
    avail[bestn] = false;
    vertex_order[idx] = bestn;
  }
  return vertex_order;
}

template <typename G>
std::vector<typename G::index_type> vertex_order_GreatestConstraintFirst(G const & g) {
  using Index = typename G::index_type;
  
  auto n = g.num_vertices();
  
  std::vector<Index> vertex_order(n);
  std::iota(std::begin(vertex_order), std::end(vertex_order), 0);
  
  enum struct Flag {
    vis,
    neigh,
    unv
  };
  
  std::vector<Flag> flags(n, Flag::unv);
  std::vector<std::tuple<Index,Index,Index>> ranks(n);
  
  for (Index i=0; i<n; ++i) {
    std::get<2>(ranks[i]) = g.degree(i);
  }
  
  for (Index m=0; m<n; ++m) {
    auto u_it = std::max_element(
        std::next(std::begin(vertex_order), m),
        std::end(vertex_order),
        [&ranks](auto u, auto v) {
          return ranks[u] < ranks[v];
        });
    Index u = *u_it;
    
    if (flags[u] == Flag::unv) {
      for (auto v : g.adjacent_vertices(u)) {
        --std::get<2>(ranks[v]);
      }
      for (auto v : g.inv_adjacent_vertices(u)) {
        --std::get<2>(ranks[v]);
      }
    } else if (flags[u] == Flag::neigh) {
      for (auto v : g.adjacent_vertices(u)) {
        --std::get<1>(ranks[v]);
      }
      for (auto v : g.inv_adjacent_vertices(u)) {
        --std::get<1>(ranks[v]);
      }
    }
    
    std::swap(vertex_order[m], *u_it);
    flags[u] = Flag::vis;
    
    for (auto v : g.adjacent_vertices(u)) {
      ++std::get<0>(ranks[v]);
      if (flags[v] == Flag::unv) {
        flags[v] = Flag::neigh;
        for (auto w : g.adjacent_vertices(v)) {
          ++std::get<1>(ranks[w]);
        }
        for (auto w : g.inv_adjacent_vertices(v)) {
          ++std::get<1>(ranks[w]);
        }
      }
    }
    for (auto v : g.inv_adjacent_vertices(u)) {
      ++std::get<0>(ranks[v]);
      if (flags[v] == Flag::unv) {
        flags[v] = Flag::neigh;
        for (auto w : g.adjacent_vertices(v)) {
          ++std::get<1>(ranks[w]);
        }
        for (auto w : g.inv_adjacent_vertices(v)) {
          ++std::get<1>(ranks[w]);
        }
      }
    }
  }
  return vertex_order;
}

#endif  // VERTEX_ORDER_H_
