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
        if(bestn == n || (rdeg > bestv || (rdeg == bestv && clustdeg[i] > clustdeg[bestn]))) {
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
  
  std::vector<bool> avail(n, true);
  
  Index u0 = 0;
  Index u0_deg = g.degree(u0);
  for (Index u=0; u<n; ++u) {
    auto u_deg = g.degree(u);
    if (u_deg > u0_deg) {
      u0 = u;
      u0_deg = u_deg;
    }
  }
  
  vertex_order[0] = u0;
  avail[u0] = false;
  
  std::vector<bool> avail_unv(n, true);
  
  for (Index m=1; m<n; ++m) {
    Index um = n; // null
    std::tuple<Index,Index,Index> um_rank{0,0,0};
    for (Index u=0; u<n; ++u) {
      if (avail[u]) {
        Index num_vis = 0;
        Index num_unv = 0;
        for (auto ui : g.adjacent_vertices(u)) {
          if (!avail[ui]) {
            ++num_vis;
          } else if(avail_unv[ui]) {
            ++num_unv;
          }
        }
        for (auto ui : g.inv_adjacent_vertices(u)) {
          if (!avail[ui]) {
            ++num_vis;
          } else if(avail_unv[ui]) {
            ++num_unv;
          }
        }
        
        std::vector<bool> V_neigh(n, false);
        for (auto v : g.adjacent_vertices(u)) {
          if (avail[v]) {
            for (auto ui : g.adjacent_vertices(v)) {
              if (!avail[ui]) {
                V_neigh[ui] = true;
              }
            }
            for (auto ui : g.inv_adjacent_vertices(v)) {
              if (!avail[ui]) {
                V_neigh[ui] = true;
              }
            }
          }
        }
        for (auto v : g.inv_adjacent_vertices(u)) {
          if (avail[v]) {
            for (auto ui : g.adjacent_vertices(v)) {
              if (!avail[ui]) {
                V_neigh[ui] = true;
              }
            }
            for (auto ui : g.inv_adjacent_vertices(v)) {
              if (!avail[ui]) {
                V_neigh[ui] = true;
              }
            }
          }
        }
        Index num_neigh = std::count(std::begin(V_neigh), std::end(V_neigh), true);
        
        auto u_rank = std::make_tuple(num_vis, num_neigh, num_unv);
        if (um == n || u_rank > um_rank) {
          um = u;
          um_rank = u_rank;  
        }
      }
    }
    vertex_order[m] = um;
    avail[um] = false;
    for (auto v : g.adjacent_vertices(um)) {
      avail_unv[v] = false;
    }
    for (auto v : g.inv_adjacent_vertices(um)) {
      avail_unv[v] = false;
    }
  }
  
  return vertex_order;
}

#endif  // VERTEX_ORDER_H_
