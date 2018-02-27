#ifndef RIIMP_STATE_H_
#define RIIMP_STATE_H_

#include <iterator>
#include <utility>
#include <vector>
#include <algorithm>
#include <numeric>

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderG>
class riimp_state_mono {
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

  std::vector<std::pair<IndexG,bool>> g_parents;

  std::vector<IndexH> map;
  std::vector<IndexG> inv;
  
  using H_adjacent_vertices_container_type = typename H::adjacent_vertices_container_type;
  
  H_adjacent_vertices_container_type h_vertices;

 public:
  riimp_state_mono(
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
        g_parents(m,{m,true}),
        map(m, n),
        inv(n, m),
        h_vertices(n) {
    std::iota(std::begin(h_vertices), std::end(h_vertices), 0);
  }
  
  riimp_state_mono(riimp_state_mono const &) = delete;

  bool empty() {
    return x_it == std::begin(index_order_g);
  }
  
  bool full() {
    return x_it == std::end(index_order_g);
  }
  
  H_adjacent_vertices_container_type candidates() {
    struct item {
      IndexH idx;
      H_adjacent_vertices_container_type const & vec;
      
      item(IndexH idx, H_adjacent_vertices_container_type const & vec)
          : idx{idx},
            vec{vec} {
      }
      
      IndexH operator*() const {
        return vec[idx];
      }
      
      item & operator++() {
        ++idx;
        return *this;
      }
      
      bool empty() const {
        return idx >= vec.size();
      }
    };
  
    struct compare {
      bool operator()(
          item const & a,
          item const & b) {
        if (!a.empty() && !b.empty()) {
          return *a < *b;
        } else {
          return a.empty() && !b.empty();
        }
      }
    };
    
    auto x = *x_it;
    std::vector<item> neighs;
    for (auto i : g.adjacent_vertices(x)) {
      auto j = map[i];
      if (j != n) {
        neighs.push_back(item(0, h.inv_adjacent_vertices(j)));
      }
    }
    for (auto i : g.inv_adjacent_vertices(x)) {
      auto j = map[i];
      if (j != n) {
        neighs.push_back(item(0, h.adjacent_vertices(j)));
      }
    }
    if (neighs.empty()) {
      return h_vertices;
    } else {
      H_adjacent_vertices_container_type result;
      while (true) {
        IndexH min = n;
        for (auto const & neigh : neighs) {
          if (neigh.empty()) {
            goto end;
          } else {
            auto i = *neigh;
            if (i < min) {
              min = i;
            }
          }
        }
        bool success = true;
        for (auto const & neigh : neighs) {
          if (*neigh != min) {
            success = false;
            break;
          }
        }
        if (success) {
          result.push_back(min);
        }
        for (auto & neigh : neighs) {
          if (*neigh == min) {
            ++neigh;
          }
        }
      }
 end:
      return result;
    }
  }

  void advance() {
  }
  
  void revert() {
  }

  bool assign(IndexH y) {
    auto x = *x_it;
    return
        inv[y] == m &&
        vertex_comp(x, y) &&
        g.out_degree(x) <= h.out_degree(y) &&
        g.in_degree(x) <= h.in_degree(y);
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
class riimp_state_ind
  : public riimp_state_mono<
        G,
        H,
        VertexEquivalencePredicate,
        EdgeEquivalencePredicate,
        IndexOrderG> {
 private:
  using base = riimp_state_mono<
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
  riimp_state_ind(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp,
      IndexOrderG const & index_order_g)
      : riimp_state_mono<
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

#endif  // RIIMP_STATE_H
