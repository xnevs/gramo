#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "ullmann_state.h"
#include "adjacency_matrix.h"
#include "simple_state.h"
#include "adjacency_set.h"
#include "compatibility_matrix.h"
#include "vertex_order.h"
#include "explore.h"

template <
    typename AdjacencySmall,
    typename AdjacencyLarge,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void ullmann_mono(
    AdjacencySmall const & g_,
    AdjacencyLarge const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_matrix g{g_};
  adjacency_matrix h{h_};
  
  auto index_order_small = vertex_order_DEG(g);
  
  ullmann_state_mono<
      adjacency_matrix,
      adjacency_matrix,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      compatibility_matrix,
      decltype(index_order_small)> S{g, h, vertex_comp, edge_comp, index_order_small};
  
  explore(S, callback);
}

template <
    typename AdjacencySmall,
    typename AdjacencyLarge,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void ullmann_ind(
    AdjacencySmall const & g_,
    AdjacencyLarge const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_matrix g{g_};
  adjacency_matrix h{h_};
  
  auto index_order_small = vertex_order_DEG(g);
  
  ullmann_state_ind<
      adjacency_matrix,
      adjacency_matrix,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      compatibility_matrix,
      decltype(index_order_small)> S{g, h, vertex_comp, edge_comp, index_order_small};
  
  explore(S, callback);
}

template <
    typename AdjacencySmall,
    typename AdjacencyLarge,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void simple_mono(
    AdjacencySmall const & g_,
    AdjacencyLarge const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_set g{g_};
  adjacency_set h{h_};

  auto index_order_small = vertex_order_DEG(g);
  
  simple_state_mono<
      adjacency_set,
      adjacency_set,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      decltype(index_order_small)> S{g, h, vertex_comp, edge_comp, index_order_small};
  
  explore(S, callback);
}

template <
    typename AdjacencySmall,
    typename AdjacencyLarge,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void simple_ind(
    AdjacencySmall const & g_,
    AdjacencyLarge const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_set g{g_};
  adjacency_set h{h_};

  auto index_order_small = vertex_order_RDEG_CNC(g);
  
  simple_state_ind<
      adjacency_set,
      adjacency_set,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      decltype(index_order_small)> S{g, h, vertex_comp, edge_comp, index_order_small};
  
  explore(S, callback);
}

template <
    typename AdjacencySmall,
    typename AdjacencyLarge,
    typename Callback,
    typename VertexEquivalencePredicate,
    typename EdgeEquivalencePredicate>
void simple_ind2(
    AdjacencySmall const & g_,
    AdjacencyLarge const & h_,
    Callback callback,
    VertexEquivalencePredicate vertex_comp,
    EdgeEquivalencePredicate edge_comp) {
  adjacency_set g{g_};
  adjacency_set h{h_};

  auto index_order_small = vertex_order_RDEG_CNC(g);
  
  simple_state_ind2<
      adjacency_set,
      adjacency_set,
      VertexEquivalencePredicate,
      EdgeEquivalencePredicate,
      decltype(index_order_small)> S{g, h, vertex_comp, edge_comp, index_order_small};
  
  explore(S, callback);
}

#endif  // FUNCTIONS_H_
