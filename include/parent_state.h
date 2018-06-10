#ifndef PARENT_STATE_H_
#define PARENT_STATE_H_

#include <iterator>
#include <vector>
#include <functional>
#include <algorithm>

#include <boost/range/iterator_range.hpp>

#include "index_heap.h"
#include "multi_stack.h"

template <
    typename G,
    typename H,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate,
    typename IndexOrderG>
class parent_state_mono {
 protected:
  using IndexG = typename G::index_type;
  using IndexH = typename H::index_type;
  
  IndexG const m;
  IndexH const n;
  
  G const & g;
  H const & h;

  VertexEquivalencePredicate vertex_comp;
  EdgeEquivalencePredicate edge_comp;
  
  std::vector<char> compatibility;
  multi_stack<std::pair<IndexG, IndexH>> compatibility_stack;
  
  std::vector<typename H::adjacent_vertices_container_type> initial_cands_vec;
  std::vector<boost::iterator_range<decltype(initial_cands_vec)::const_iterator>> initial_cands;
  std::vector<boost::iterator_range<decltype(initial_cands_vec)::const_iterator>> cands;
  multi_stack<IndexG> parent_stack;

  std::vector<IndexH> score;

  struct compare {
    std::vector<IndexH> const & score;
    
    compare(std::vector<IndexH> const & score)
        : score{score} {
    }
    
    bool operator()(IndexG a, IndexG b) {
      return score[a] > score[b];
    }
  };

  index_heap<IndexG, compare> g_heap;
  
 public:
  parent_state_mono(
      G const & g,
      H const & h,
      VertexEquivalencePredicate const & vertex_comp,
      EdgeEquivalencePredicate const & edge_comp)
      : m{g.num_vertices()},
        n{h.num_vertices()},
        g{g},
        h{h},
        compatibility(m*n),
        compatibility_stack(m*n),
        initial_cands_vec(m),
        initial_cands(m),
        cands(m),
        parent_stack(m),
        score(m),
        g_heap(m, compare(score)) {
    for (IndexG i=0; i<m; ++i) {
      for (IndexH j=0; j<m; ++j) {
        if (vertex_comp(i, j) &&
            g.out_degree(i) <= h.out_degree(j) &&
            g.in_degree(i) <= h.in_degree(j)) {
          auto idx = i*n + j;
          compatibility[idx] = true;
          initial_cands_vec[i].push_back(j);
        }
      }
      cands[i] = initial_cands[i] = boost::make_iterator_range(
          std::cbegin(initial_cands_vec[i]),
          std::cend(initial_cands_vec[i]));
      score[i] = n + cands[i].size();
    }
    g_heap.heapify();
  }
  
  parent_state_mono(parent_state_mono const &) = delete;

  bool empty() {
    return g_heap.full();
  }
  
  bool full() {
    return g_heap.empty();
  }
  
  void prepare() {
  }
  
  void forget() {
  }
  
  auto candidates() {
    auto x = g_heap.top();
    return cands[x];
  }

  void advance() {
    parent_stack.push_level();
    compatibility_stack.push_level();
  }
  
  void revert() {
    compatibility_stack.pop_level();
    parent_stack.pop_level();
  }

  bool assign(IndexH y) {
    auto x = g_heap.top();
    return compatibility[x*n+y];
  }

  void push(IndexH y) {
    
    // set parents
    // topology condition
    
    
    // neighborhood filter
  }
  
  IndexH pop() {
    while (!compatibility_stack.level_empty()) {
      IndexG i;
      IndexH j;
      std::tie(i, j) = compatibility_stack.top();
      compatibility_stack.pop();
      compatibility[i*n+j] = true;
      ++score[i];
    }
    while(!parent_stack.level_empty()) {
      auto i = parent_stack.top();
      parent_stack.pop();
      cands[i] = initial_cands[i];
      score[i] = n + cands[i].size();
    }
  }
};

#endif  // PARENT_STATE_H
